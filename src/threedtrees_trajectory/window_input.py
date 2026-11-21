"""Seekable LAS/LAZ scan ranges, indexed in record order by GPS time."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cloud import TRACE_FIELDS, XYZ_TIME_FIELDS, is_copc, open_cloud, point_chunks, pulse_chunks, read_points_exact
from .pulses import TIME_TOLERANCE_S, _repair_times

INDEX_POINTS = 200_000
TRACE_BATCH = 500_000
TRACE_LIMIT = np.iinfo(np.int32).max


@dataclass(frozen=True)
class ScanBlock:
    start: int
    stop: int
    time_min: float
    time_max: float
    lower: np.ndarray
    upper: np.ndarray
    sample: np.ndarray | None = None

    def overlaps(self, lower, upper):
        return self.time_min < upper and self.time_max >= lower


class ScanIndex:
    """Keep metadata for source reads, including clock corrections accepted by the path reader."""

    def __init__(self, paths, grid):
        self.paths = paths
        self.grid = grid
        self.blocks = [[] for _ in paths]
        self.summaries = [[] for _ in paths]
        self.repaired = set()
        self.copc = False
        self.source_identity = None

    def observe(self, index, chunk):
        if self.copc:
            return
        blocks = self.blocks[index]
        start = blocks[-1].stop if blocks else 0
        size = len(chunk) if index in self.repaired else INDEX_POINTS
        for offset in range(0, len(chunk), size):
            block = chunk[offset:offset + size]
            time = block.gps_time
            raw = (block.X, block.Y, block.Z)
            lower = np.array([axis.min() for axis in raw]) * block.scales + block.offsets
            upper = np.array([axis.max() for axis in raw]) * block.scales + block.offsets
            blocks.append(ScanBlock(start + offset, start + offset + len(block),
                                    float(time.min()), float(time.max()), lower, upper,
                                    time[::max(1, len(block) // 16)].copy()))
            if self.grid is not None:
                self.summaries[index].append(_surface_summary(block, self.grid))

    def begin_repair(self, index):
        """Build the repaired index during the required fallback read."""
        original = self.blocks[index], self.summaries[index]
        self.blocks[index], self.summaries[index] = [], []
        self.repaired.add(index)
        return original

    def abandon_repair(self, index, original):
        self.blocks[index], self.summaries[index] = original
        self.repaired.remove(index)

    def populate(self):
        for index, path in enumerate(self.paths):
            with open_cloud(path, read_evlrs=False) as reader:
                if is_copc(reader.header):
                    raise ValueError("Occlusion requires scanner-order LAS/LAZ. Supply the standardized files instead of COPC.")
            for chunk in point_chunks(path, INDEX_POINTS, fields=XYZ_TIME_FIELDS):
                self.observe(index, chunk)

    def surfaces(self, lower, upper, stats):
        """Merge first-read column limits and decode only blocks crossing the path endpoints."""
        grid = self.grid
        nx, ny, nz = grid.shape
        low = np.full(nx * ny, nz, np.int32)
        high = np.full(nx * ny, -1, np.int32)
        points_read = points_kept = boundary_points = 0
        for file, blocks in enumerate(self.blocks):
            boundary = []
            for block, summary in zip(blocks, self.summaries[file], strict=True):
                count = block.stop - block.start
                points_read += count
                if not block.overlaps(lower, upper):
                    continue
                if block.time_min < lower or block.time_max >= upper:
                    boundary.append(block)
                    continue
                columns, z_min, z_max = summary
                low[columns] = np.minimum(low[columns], z_min)
                high[columns] = np.maximum(high[columns], z_max)
                points_kept += count
            for chunk in self.read(file, boundary, group_pulses=False):
                t = chunk.gps_time
                keep = (t >= lower) & (t < upper)
                boundary_points += len(t)
                points_kept += int(keep.sum())
                if keep.any():
                    columns, z_min, z_max = _surface_summary(chunk[keep], grid)
                    low[columns] = np.minimum(low[columns], z_min)
                    high[columns] = np.maximum(high[columns], z_max)
        if points_read != grid.points:
            raise ValueError("The occlusion index does not cover the point count declared by its headers.")
        stats.update(echoes_read=points_read, echoes_outside_path=points_read - points_kept,
                     boundary_points_decoded=boundary_points, bytes_read=0)
        return low.reshape(nx, ny), high.reshape(nx, ny)

    def release_surfaces(self):
        self.summaries.clear()

    def read(self, file, blocks, *, group_pulses=True):
        """Read selected blocks in scanner order, coalescing consecutive record ranges."""
        if not blocks:
            return
        fields = TRACE_FIELDS if group_pulses or file in self.repaired else XYZ_TIME_FIELDS
        with open_cloud(self.paths[file], sequential=True, fields=fields) as reader:
            if file in self.repaired:
                for block in blocks:
                    reader.seek(block.start)
                    yield _repair_times(read_points_exact(reader, block.stop - block.start))[0]
            elif group_pulses:
                for contiguous in _contiguous_ranges(blocks):
                    yield from pulse_chunks(_read_ranges(reader, contiguous), TRACE_BATCH)
            else:
                yield from _read_ranges(reader, blocks)


def _surface_summary(chunk, grid):
    """Lossless extrema for touched columns, with a bounded dense fast path."""
    nx, ny, nz = grid.shape
    x, y, z = [np.floor((axis * scale + offset - origin) / grid.voxel + 1e-8).astype(np.int64)
               for axis, scale, offset, origin in zip((chunk.X, chunk.Y, chunk.Z),
                                                      chunk.scales, chunk.offsets, grid.origin)]
    for axis, maximum in zip((x, y, z), (nx, ny, nz)):
        np.clip(axis, 0, maximum - 1, out=axis)
    x0, y0 = int(x.min()), int(y.min())
    width, height = int(x.max()) - x0 + 1, int(y.max()) - y0 + 1
    if width * height <= 2 * len(chunk):
        key = (x - x0) * height + y - y0
        low = np.full(width * height, nz, np.int32)
        high = np.full(width * height, -1, np.int32)
        np.minimum.at(low, key, z)
        np.maximum.at(high, key, z)
        occupied = np.flatnonzero(high >= 0)
        columns = (x0 + occupied // height) * ny + y0 + occupied % height
        low, high = low[occupied], high[occupied]
    else:
        key = x * ny + y
        order = np.argsort(key, kind="stable")
        ordered = key[order]
        starts = np.flatnonzero(np.r_[True, ordered[1:] != ordered[:-1]])
        columns = ordered[starts]
        low = np.minimum.reduceat(z[order], starts).astype(np.int32)
        high = np.maximum.reduceat(z[order], starts).astype(np.int32)
    return columns.astype(np.uint32) if nx * ny <= 2**32 else columns, low, high


def _contiguous_ranges(blocks):
    start = 0
    for number in range(1, len(blocks)):
        if blocks[number - 1].stop != blocks[number].start:
            yield blocks[start:number]
            start = number
    yield blocks[start:]


def _read_ranges(reader, blocks):
    for block in blocks:
        if reader.points_read != block.start:
            reader.seek(block.start)
        # Each indexed block is independently bounded. The pulse reader holds the tail.
        yield read_points_exact(reader, block.stop - block.start)


# Compact echo records for pulse grouping. Only fragments are retained between batches.
ECHO = np.dtype([("time", "<f8"), ("x", "<i4"), ("y", "<i4"), ("z", "<i4"),
                 ("source", "<u2"), ("returns", "u1"), ("channel", "u1")])


def pack_echoes(points, keep):
    records = np.empty(np.count_nonzero(keep), ECHO)
    for raw, axis in zip(("X", "Y", "Z"), ("x", "y", "z")):
        records[axis] = getattr(points, raw)[keep]
    records["time"] = points.gps_time[keep]
    records["source"] = points.point_source_id[keep]
    records["returns"] = np.asarray(points.return_number)[keep] | (np.asarray(points.number_of_returns)[keep] << 4)
    records["channel"] = np.asarray(points.scanner_channel)[keep] if points.point_format.id >= 6 else 0
    return records


def group_multi(records, *, acquisition_order):
    """Match the original complete runs, or join separated returns by timestamp and source."""
    rn, nr = records["returns"] & 15, records["returns"] >> 4
    time, source, channel = records["time"], records["source"], records["channel"]
    same = (time[1:] == time[:-1]) & (source[1:] == source[:-1]) & (channel[1:] == channel[:-1])
    if acquisition_order:
        index = np.arange(len(records))
        group_start = np.maximum.accumulate(np.where(np.r_[True, ~same], index, 0))
        previous_first = np.maximum.accumulate(np.where(rn == 1, index, -1))
        same &= ~((rn[1:] == 1) & (previous_first[:-1] >= group_start[1:]))
    starts = np.flatnonzero(np.r_[True, ~same])
    sizes = np.diff(np.r_[starts, len(records)])
    target = nr[starts]
    numbers = np.bitwise_or.reduceat(np.left_shift(np.uint32(1), rn), starts)
    complete = (sizes == target) & (numbers == np.left_shift(np.uint32(1), target + 1) - 2)
    return sizes, complete
