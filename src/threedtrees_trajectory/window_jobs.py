"""Lock a job's checkpoints and admit window workers within its CPU and memory budgets."""
from __future__ import annotations

import fcntl
import hashlib
import json
import multiprocessing
import os
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from contextlib import contextmanager

import numpy as np

from .resources import available_memory

SCHEMA = 3


class JobBusy(ValueError):
    """Another process holds this output's checkpoint lock."""


def identify(paths, parameters, time, xyz):
    inputs = []
    for path in paths:
        stat = path.stat()
        inputs.append((str(path.resolve()), stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns))
    key = hashlib.sha256(json.dumps((SCHEMA, inputs, parameters), sort_keys=True).encode())
    for array in (time, xyz):
        key.update(np.ascontiguousarray(array, dtype="<f8").tobytes())
    return key.hexdigest()[:24]


@contextmanager
def lock_job(work):
    work.mkdir(exist_ok=True)
    fd = os.open(work, os.O_RDONLY)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise JobBusy(f"Another occlusion run owns checkpoint {work.name}.") from error
        yield
    finally:
        os.close(fd)


def atomic_json(path, data):
    pending = path.with_suffix(".tmp")
    with pending.open("w", encoding="utf-8") as stream:
        json.dump(data, stream, allow_nan=False)
        stream.flush()
        os.fsync(stream.fileno())
    pending.replace(path)
    fd = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def window_files(work, index):
    return tuple(work / f"window-{index:05d}.{suffix}" for suffix in ("npy", "json"))


def finished(work, index, shape):
    counts, marker = window_files(work, index)
    if marker.is_file() and counts.is_file():
        try:
            array = np.load(counts, mmap_mode="r", allow_pickle=False)
            if array.shape == shape and array.dtype == np.uint32:
                return json.loads(marker.read_text(encoding="utf-8"))
        except (OSError, ValueError, EOFError):
            pass
    return None


def execute(work, specs, slots, worker, init, initargs, stage):
    """Admit jobs in the supplied priority order subject to CPU and memory reservations."""
    results, pending = {}, []
    for index, amount, shape, args in specs:
        previous = finished(work, index, shape)
        if previous is None:
            pending.append((index, amount, args))
        else:
            results[index] = previous
    memory = available_memory()
    stage.update(reused_windows=len(results), scheduled_windows=len(pending), worker_slots=slots,
                 worker_memory_bytes=memory, peak_concurrent_workers=0, peak_reserved_bytes=0)
    if pending and max(amount for _, amount, _ in pending) > memory:
        raise ValueError("A planned occlusion window exceeds available job memory. Lower --grid-memory-gib "
                         "or request more job memory.")
    if not pending:
        return results

    reserved, running = 0, {}
    try:
        with ProcessPoolExecutor(max_workers=min(len(pending), slots), mp_context=multiprocessing.get_context("spawn"),
                                 initializer=init, initargs=initargs) as pool:
            while pending or running:
                free = memory - reserved
                for job in pending[:]:
                    if len(running) >= slots:
                        break
                    index, amount, args = job
                    if amount <= free:
                        future = pool.submit(worker, index, *args)
                        running[future] = (index, amount)
                        pending.remove(job)
                        reserved += amount
                        free -= amount
                stage["peak_concurrent_workers"] = max(stage["peak_concurrent_workers"], len(running))
                stage["peak_reserved_bytes"] = max(stage["peak_reserved_bytes"], reserved)
                done, _ = wait(running, return_when=FIRST_COMPLETED)
                for future in done:
                    index, amount = running.pop(future)
                    reserved -= amount
                    results[index] = future.result()
    finally:
        stage["completed_windows"] = len(results)
    return results
