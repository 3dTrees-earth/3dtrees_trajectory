"""Plan and trace disjoint GPS intervals, then sum their overlapping voxel counts."""
from __future__ import annotations

import math
import os
import resource
from bisect import bisect_right
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from time import perf_counter, process_time

import numpy as np

from .ray_walk import Columns, Walker, sum_block
from .resources import _proc_rss, _rss, available_memory
from .window_input import TRACE_BATCH, TRACE_LIMIT, group_multi, pack_echoes
from .window_jobs import atomic_json, window_files

GIB = 2**30
WORKING_BYTES = 2 * GIB
# Concatenation briefly holds two 24-byte copies per echo. The indirect
# sort needs an index and workspace, while grouping only expands one batch.
FRAGMENT_BYTES = 64
COUNT_BYTES = 12
COLUMN_BYTES = 24
MAX_WINDOWS = 10_000


@dataclass(frozen=True)
class TimeWindow:
    lower: float
    upper: float
    bounds: tuple[slice, slice, slice]
    blocks: tuple
    echoes: int
    count_cells: int

    @property
    def shape(self):
        return tuple(part.stop - part.start for part in self.bounds)

    @property
    def memory(self):
        return (WORKING_BYTES + COUNT_BYTES * self.count_cells + COLUMN_BYTES * math.prod(self.shape[:2])
                + FRAGMENT_BYTES * self.echoes)


def _path_bounds(block, lower, upper, time, xyz):
    first, last = max(lower, block.time_min), min(upper, block.time_max)
    left, right = np.searchsorted(time, [first, last])
    endpoints = np.array([np.interp([first, last], time, axis) for axis in xyz.T]).T
    points = np.concatenate((endpoints, xyz[left:right]))
    return points.min(axis=0), points.max(axis=0)


def _ray_bounds(block, lower, upper, grid, time, xyz):
    """Conservatively enclose rays up to the original grid exit."""
    sensor_lo, sensor_hi = _path_bounds(block, lower, upper, time, xyz)
    grid_lo, grid_hi = grid.origin, grid.origin + np.array(grid.shape) * grid.voxel
    echo_lo = block.lower
    echo_hi = block.upper + grid.voxel * 1e-7
    direction_lo, direction_hi = echo_lo - sensor_hi, echo_hi - sensor_lo
    limits = np.full(3, np.inf)
    forward, backward = direction_lo > 0, direction_hi < 0
    limits[forward] = np.maximum(0., grid_hi[forward] - echo_lo[forward]) / direction_lo[forward]
    limits[backward] = np.maximum(0., echo_hi[backward] - grid_lo[backward]) / -direction_hi[backward]
    extension = limits.min()
    if not np.isfinite(extension):
        return grid_lo, grid_hi
    low = np.minimum(sensor_lo, echo_lo + extension * np.minimum(direction_lo, 0))
    high = np.maximum(sensor_hi, echo_hi + extension * np.maximum(direction_hi, 0))
    return np.maximum(grid_lo, low), np.minimum(grid_hi, high)


def _plan_window(index, grid, time, xyz, lower, upper, layer_ranges):
    files = tuple(tuple(block for block in blocks if block.overlaps(lower, upper)) for blocks in index.blocks)
    blocks = [block for file in files for block in file]
    if not blocks:
        return None
    bounds = np.asarray([_ray_bounds(block, lower, upper, grid, time, xyz) for block in blocks])
    start = np.maximum(0, np.floor((bounds[:, 0].min(axis=0) - grid.origin) / grid.voxel).astype(int) - 1)
    stop = np.minimum(grid.shape, np.floor((bounds[:, 1].max(axis=0) - grid.origin) / grid.voxel).astype(int) + 2)
    extent = tuple(slice(int(lo), int(hi)) for lo, hi in zip(start, stop))
    cells = Columns.from_bounds(extent, *layer_ranges).cells
    return TimeWindow(lower, upper, extent, files, sum(block.stop - block.start for block in blocks), cells)


def plan(index, grid, time, xyz, time_bounds, slots, budget, layer_ranges):
    if budget <= WORKING_BYTES:
        raise ValueError("The occlusion job lacks memory for one GPS-time worker.")
    lower, upper = time_bounds
    samples = [block.sample for blocks in index.blocks for block in blocks
               if block.sample is not None and block.overlaps(lower, upper)]
    if samples and slots > 1:
        seconds = np.concatenate(samples)
        seconds = seconds[(lower <= seconds) & (seconds < upper)]
        cuts = np.quantile(seconds, np.arange(1, slots) / slots) if len(seconds) else []
        edges = np.unique(np.r_[lower, cuts, upper])
    else:
        edges = np.unique(np.linspace(lower, upper, slots + 1))
    pending = list(zip(edges[:-1], edges[1:]))[::-1]
    windows = []
    while pending:
        lower, upper = pending.pop()
        window = _plan_window(index, grid, time, xyz, lower, upper, layer_ranges)
        if window is None:
            continue
        if window.memory <= budget and window.echoes <= TRACE_LIMIT:
            windows.append(window)
            continue
        middle = lower / 2 + upper / 2
        if not lower < middle < upper or len(windows) + len(pending) >= MAX_WINDOWS - 1:
            raise ValueError("A GPS-time window exceeds the occlusion grid memory or native echo budget. "
                             "Increase --voxel or --grid-memory-gib, or restrict the acquisition.")
        pending.extend(((middle, upper), (lower, middle)))
    return windows


def init_worker(index, grid, time, xyz, work, single_files, layer_ranges, coverage):
    global _CONTEXT
    _CONTEXT = index, grid, time, xyz, work, single_files, layer_ranges, coverage


def trace_window(number, window):
    index, grid, time, xyz, work, single_files, layer_ranges, coverage = _CONTEXT
    start_wall, start_cpu = perf_counter(), process_time()
    try:
        Path("/proc/self/clear_refs").write_text("5", encoding="ascii")
    except OSError:
        pass
    columns = Columns.from_bounds(window.bounds, *layer_ranges)
    destination, marker = window_files(work, number)
    pending = destination.with_suffix(".part")
    counts = np.lib.format.open_memmap(pending, mode="w+", dtype=np.uint32, shape=(columns.cells, 3))
    tracer = Walker(grid, columns, counts)
    decoded = recovered = unpaired = traced = 0
    for file, blocks in enumerate(window.blocks):
        fragments = []
        scales, offsets = None, None
        for chunk in index.read(file, blocks):
            decoded += len(chunk)
            scales, offsets = chunk.scales, chunk.offsets
            t = chunk.gps_time
            keep = (t >= window.lower) & (t < window.upper)
            if not keep.any():
                continue
            force_single = file in single_files
            one = (np.ones(len(t), bool) if force_single else
                   (np.asarray(chunk.return_number) <= 1) & (np.asarray(chunk.number_of_returns) <= 1))
            single = keep & one
            if single.any():
                seconds = t[single]
                coords = [axis[single] * scale + offset + grid.voxel * 1e-7 for axis, scale, offset in
                          zip((chunk.X, chunk.Y, chunk.Z), scales, offsets)]
                tracer.add(coords, _positions(seconds, time, xyz))
                traced += len(seconds)
            if (keep & ~one).any():
                # Group in acquisition order before selecting this window.
                # Removing other GPS times first can turn separated echoes
                # into an apparently consecutive multi-return pulse.
                records = pack_echoes(chunk, ~one & (t >= coverage[0]) & (t < coverage[1]))
                sizes, complete = group_multi(records, acquisition_order=True)
                first = np.cumsum(sizes) - sizes
                clocks = records["time"][first]
                selected = (clocks >= window.lower) & (clocks < window.upper)
                paired = selected & complete
                if paired.any():
                    _trace_multi(tracer, records[np.repeat(paired, sizes)], sizes[paired],
                                 scales, offsets, grid.voxel, time, xyz)
                    traced += int(sizes[paired].sum())
                if (selected & ~complete).any():
                    fragments.append(records[np.repeat(selected & ~complete, sizes)])
        if fragments:
            records = np.concatenate(fragments)
            fragments.clear()
            joined = 0
            for batch in _sorted_fragments(records):
                sizes, complete = group_multi(batch, acquisition_order=False)
                if complete.any():
                    paired = np.repeat(complete, sizes)
                    _trace_multi(tracer, batch[paired], sizes[complete], scales, offsets, grid.voxel, time, xyz)
                    joined += int(paired.sum())
            recovered += joined
            unpaired += len(records) - joined
            traced += joined
            del records
    traversed, missing = map(int, tracer.stats)
    del tracer
    counts.flush()
    del counts
    with pending.open("rb") as stream:
        os.fsync(stream.fileno())
    pending.replace(destination)
    data = dict(window=number + 1, gps_from=window.lower, gps_until=window.upper, points_read=decoded,
                echoes_traced=traced, recovered_echoes=recovered, unpaired_echoes=unpaired,
                traversed_pulses=traversed, missing_returns=missing,
                grid_shape=window.shape, grid_voxels=math.prod(window.shape), count_voxels=columns.cells,
                wall_seconds=perf_counter() - start_wall, cpu_seconds=process_time() - start_cpu,
                peak_rss_bytes=_proc_rss('self', 'VmHWM') or _rss(), bytes_written=destination.stat().st_blocks * 512)
    atomic_json(marker, data)
    return data


def _sorted_fragments(records):
    order = np.lexsort((records["returns"] & 15, records["channel"], records["source"], records["time"]))

    def key(index):
        point = records[index]
        return point["time"], point["source"], point["channel"]

    start = 0
    while start < len(order):
        stop = min(start + TRACE_BATCH, len(order))
        batch = records[order[start:stop]]
        if stop < len(order) and (last := key(order[stop - 1])) == key(order[stop]):
            # Preserve a complete time/source/channel group across batch edges.
            different = np.flatnonzero((batch["time"] != last[0]) | (batch["source"] != last[1])
                                       | (batch["channel"] != last[2]))
            if not len(different):
                # A group of this size cannot describe a LAS pulse. Its echoes
                # still count toward the original fragment total.
                start = bisect_right(order, last, lo=stop, key=key)
                continue
            batch = batch[:different[-1] + 1]
            stop = start + len(batch)
        yield batch
        start = stop


def _positions(seconds, time, xyz):
    return (np.interp(seconds, time, axis) for axis in xyz.T)


def _trace_multi(tracer, records, sizes, scales, offsets, voxel, time, xyz):
    points = [records[axis] * scale + offset + voxel * 1e-7
              for axis, scale, offset in zip(("x", "y", "z"), scales, offsets)]
    pulses = np.r_[0, np.cumsum(sizes, dtype=np.int32)].astype(np.int32, copy=False)
    rn = records["returns"] & 15
    tracer.add(points, _positions(records["time"][pulses[:-1]], time, xyz), rn, pulses)


def count_blocks(work, windows, layer_ranges):
    """Build a reader for individual product blocks with shared, read-only count maps."""
    def open_window(number):
        values = np.load(window_files(work, number)[0], mmap_mode="r", allow_pickle=False)
        return Columns.from_bounds(windows[number].bounds, *layer_ranges), values

    # An XY tile can overlap every time window. A small LRU then rebuilds
    # whole-window column offsets on nearly every lookup. Retain all the
    # offsets and virtual maps when the job has room for them.
    sizes = [12 * math.prod(window.shape[:2]) + 8 for window in windows]
    cache_bytes = available_memory() // 2
    fd_limit, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
    map_limit = len(windows) if fd_limit == resource.RLIM_INFINITY else max(1, fd_limit - 64)
    if len(windows) <= map_limit and sum(sizes) <= cache_bytes:
        stored = tuple(open_window(number) for number in range(len(windows))).__getitem__
    else:
        stored = lru_cache(maxsize=max(1, min(32, map_limit, cache_bytes // max(sizes))))(open_window)

    def read(block):
        start = np.array([part.start for part in block], dtype=np.int32)
        stop = np.array([part.stop for part in block], dtype=np.int32)
        counts = np.zeros((*tuple(stop - start), 3), np.uint64)
        for number, window in enumerate(windows):
            if any(a.start >= b.stop or b.start >= a.stop for a, b in zip(window.bounds, block)):
                continue
            columns, values = stored(number)
            sum_block(columns, values, start, stop, counts)
        if counts.max(initial=0) > np.iinfo(np.uint32).max:
            raise ValueError("Too many pulses traversed a voxel for the 32-bit output count.")
        return np.moveaxis(counts.astype(np.uint32), -1, 0)

    return read


def summed_counts(work, windows, blocks, layer_ranges):
    read = count_blocks(work, windows, layer_ranges)
    for block in blocks:
        yield block[:2], block[2], read(block)
