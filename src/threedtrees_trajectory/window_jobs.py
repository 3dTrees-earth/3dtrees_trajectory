"""Lock a job's checkpoints and admit window workers within its CPU and memory budgets."""
from __future__ import annotations

import errno
import fcntl
import hashlib
import json
import multiprocessing
import os
import signal
from collections import deque
from multiprocessing.connection import wait
from contextlib import contextmanager
from dataclasses import dataclass

import numpy as np

from .resources import available_memory

SCHEMA = 4


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


@contextmanager
def lock_output(output):
    """Serialize CLI publication for a single pair of occlusion destinations."""
    output.parent.mkdir(parents=True, exist_ok=True)
    work = output.parent / f".trajectory-output-{identify([], str(output), np.empty(0), np.empty(0))}"
    with lock_job(work):
        try:
            yield
        finally:
            work.rmdir()


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


@dataclass
class WorkerMemory:
    """Calibrate the planner's relative echo, voxel and column costs against peak RSS."""

    scale: float = 0.
    overhead = 256 * 2**20

    def observe(self, planned, result):
        peak = result.get("peak_rss_bytes", 0)
        if peak > 0:
            self.scale = max(self.scale, peak / planned)

    def reserve(self, planned):
        return self.overhead + int(1.25 * self.scale * planned) if self.scale else planned


def _run_window(sender, worker, init, initargs, index, args):
    with sender:
        try:
            if init is not None:
                init(*initargs)
            result = worker(index, *args)
        except Exception as error:
            sender.send((None, error))
        else:
            sender.send((result, None))


def execute(work, specs, slots, worker, init, initargs, stage):
    """Admit fresh processes by estimated RSS and retry suspected OOMs in isolation."""
    results, pending = {}, []
    model = WorkerMemory()
    for index, planned, shape, args in specs:
        previous = finished(work, index, shape)
        if previous is None:
            pending.append((index, planned, shape, args))
        else:
            results[index] = previous
            model.observe(planned, previous)
    memory = available_memory()
    stage.update(reused_windows=len(results), scheduled_windows=len(pending), worker_slots=slots,
                 worker_memory_bytes=memory, peak_concurrent_workers=0, peak_reserved_bytes=0,
                 window_reservations_bytes={}, memory_retries=0)
    context = multiprocessing.get_context("spawn")
    retry, running = deque(), {}
    reserved = 0

    def start(job, amount, isolated=False):
        nonlocal reserved
        index, _, _, args = job
        receiver, sender = context.Pipe(duplex=False)
        process = context.Process(target=_run_window, args=(sender, worker, init, initargs, index, args))
        try:
            process.start()
        except BaseException:
            receiver.close()
            raise
        finally:
            sender.close()
        running[receiver] = process, job, amount, isolated
        reserved += amount
        stage["window_reservations_bytes"].setdefault(index + 1, []).append(amount)
        stage["peak_concurrent_workers"] = max(stage["peak_concurrent_workers"], len(running))
        stage["peak_reserved_bytes"] = max(stage["peak_reserved_bytes"], reserved)

    def stopped(index):
        return ValueError(f"Occlusion window {index + 1} stopped. Repeat the command to resume from {work.name}.")

    try:
        while pending or running or retry:
            free = min(memory - reserved, available_memory())
            if retry:
                if not running:
                    job, previous = retry.popleft()
                    if free <= previous:
                        raise stopped(job[0])
                    start(job, free, isolated=True)
                    stage["memory_retries"] += 1
            elif not any(isolated for _, _, _, isolated in running.values()):
                for job in pending[:]:
                    if len(running) >= slots:
                        break
                    _, planned, _, _ = job
                    # A pessimistic estimate can still try the window with the entire job budget.
                    amount = min(model.reserve(planned), free) if not running and planned <= free else model.reserve(planned)
                    if amount <= free:
                        start(job, amount)
                        pending.remove(job)
                        free -= amount
            if not running:
                raise ValueError("A planned occlusion window exceeds available job memory. Lower --grid-memory-gib "
                                 "or request more job memory.")

            for receiver in wait(running):
                process, job, amount, isolated = running[receiver]
                index, planned, shape, _ = job
                try:
                    result, error = receiver.recv()
                except EOFError:
                    result, error = finished(work, index, shape), None
                process.join()
                exit_code = process.exitcode
                process.close()
                receiver.close()
                del running[receiver]
                reserved -= amount
                if result is not None:
                    results[index] = result
                    model.observe(planned, result)
                elif (exit_code == -signal.SIGKILL or isinstance(error, MemoryError)
                      or isinstance(error, OSError) and error.errno == errno.ENOMEM):
                    if isolated:
                        raise stopped(index) from error
                    retry.append((job, amount))
                elif error is not None:
                    raise error
                else:
                    raise stopped(index)
    finally:
        for receiver, (process, _, _, _) in running.items():
            if process.is_alive():
                process.terminate()
            process.join()
            process.close()
            receiver.close()
        stage["completed_windows"] = len(results)
    return results
