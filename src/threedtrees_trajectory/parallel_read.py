"""Decode and summarize independent LAS/LAZ record ranges in scanner order."""
from __future__ import annotations

import multiprocessing
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from dataclasses import replace
from itertools import islice

import numpy as np

from .cloud import RETURN_FIELDS, open_cloud, read_points_exact
from .coverage import UnusableTime
from .resources import available_memory

MAX_TAIL = 65_536
WORKER_BYTES = 1_500_000_000


class SerialReadRequired(Exception):
    """The acquisition's return grouping needs the sequential reader."""


def _tail_start(times):
    if len(times) < 2:
        return 0
    if times[-2] != times[-1]:
        return len(times) - 1
    # LAS preserves record order, so the equal-time suffix only needs a short probe.
    end = max(0, len(times) - MAX_TAIL)
    boundaries = np.flatnonzero(times[end:-1] != times[-1])
    if len(boundaries):
        return end + int(boundaries[-1]) + 1
    if end == 0:
        return 0
    raise SerialReadRequired


def _previous_tail(reader, start):
    if start == 0:
        return 0
    offset = min(start, 64)
    while offset <= min(start, MAX_TAIL):
        reader.seek(start - offset)
        tail = _tail_start(read_points_exact(reader, offset).gps_time)
        if tail or offset == start:
            return start - offset + tail
        offset = min(start, offset * 2)
    raise SerialReadRequired


def _collect_range(path, index, start, stop, total, shift, strides, grid, probe):
    from .pulses import _FilePulses
    from .window_input import ScanIndex

    with open_cloud(path, sequential=True) as reader:
        left = _previous_tail(reader, start)
        if reader.points_read != left:
            reader.seek(left)
        chunk = read_points_exact(reader, stop - left)
        if stop < total:
            chunk = chunk[:_tail_start(chunk.gps_time)]
    if not len(chunk) or (start == 0 and len(chunk) < probe):
        raise SerialReadRequired
    minimum, maximum = float(chunk.gps_time.min()), float(chunk.gps_time.max())
    if not np.isfinite([minimum, maximum]).all():
        raise UnusableTime("GPS times must be finite.")
    local = _FilePulses([path])
    local.collect(chunk, index, shift, strides, False)
    blocks, summaries = [], []
    if grid is not None:
        scan = ScanIndex([path], grid)
        scan.observe(0, chunk)
        blocks = [replace(block, start=block.start + left, stop=block.stop + left) for block in scan.blocks[0]]
        summaries = scan.summaries[0]
    return local, minimum, maximum, len(chunk), blocks, summaries


def probe_single_returns(path, size):
    from .pulses import _looks_single_return

    with open_cloud(path, sequential=True, fields=RETURN_FIELDS) as reader:
        return _looks_single_return(read_points_exact(reader, size))


def collect_ranges(path, index, total, batch, record_size, shift, strides, scan_index, slots, probe):
    from .pulses import _FilePulses

    # Reserve the decoded record plus the temporary endpoint, source and index arrays.
    reservation = max(WORKER_BYTES, batch * (record_size + 224))
    slots = min(slots, -(-total // batch), available_memory() // reservation)
    if slots < 2:
        raise SerialReadRequired
    local = _FilePulses([path])
    grid = scan_index.grid if scan_index is not None else None
    minimum, maximum, points = np.inf, -np.inf, 0
    blocks, summaries = [], []
    starts = iter(range(0, total, batch))

    def submit(pool, start):
        return pool.submit(_collect_range, path, index, start, min(start + batch, total),
                           total, shift, strides, grid, probe)

    try:
        with ProcessPoolExecutor(max_workers=slots, mp_context=multiprocessing.get_context("spawn")) as pool:
            pending = deque(submit(pool, start) for start in islice(starts, slots))
            while pending:
                part, low, high, count, part_blocks, part_summaries = pending.popleft().result()
                following = next(starts, None)
                if following is not None:
                    pending.append(submit(pool, following))
                local.append_batch(part)
                if local.stray_limit_exceeded:
                    raise SerialReadRequired
                minimum, maximum = min(minimum, low), max(maximum, high)
                points += count
                blocks.extend(part_blocks)
                summaries.extend(part_summaries)
    except BrokenProcessPool as error:
        raise ValueError("A scanner-order read worker stopped before completing its record range.") from error
    if points != total:
        raise SerialReadRequired
    if scan_index is not None:
        scan_index.blocks[index] = blocks
        scan_index.summaries[index] = summaries
    return local, minimum, maximum, points, slots
