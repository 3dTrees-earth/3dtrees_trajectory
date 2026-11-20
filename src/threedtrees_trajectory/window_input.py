"""Seekable LAS/LAZ scan ranges, indexed in record order by GPS time."""
from __future__ import annotations

from dataclasses import dataclass

import laspy
import numpy as np

from .cloud import LAZ_FIELDS, is_copc, point_chunks, pulse_chunks
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

    def __init__(self, paths):
        self.paths = paths
        self.blocks = [[] for _ in paths]
        self.repaired = set()
        self.copc = False

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

    def reset(self, index):
        self.blocks[index].clear()
        self.repaired.add(index)

    def populate(self):
        for index, path in enumerate(self.paths):
            with laspy.open(path, read_evlrs=False) as reader:
                if is_copc(reader.header):
                    raise ValueError("Occlusion requires scanner-order LAS/LAZ. Supply the standardized files instead of COPC.")
            for chunk in point_chunks(path, INDEX_POINTS):
                self.observe(index, chunk)

    def read(self, file, blocks):
        """Read selected blocks in scanner order, coalescing consecutive record ranges."""
        if not blocks:
            return
        arguments = dict(laz_backend=laspy.LazBackend.Lazrs, decompression_selection=LAZ_FIELDS)
        with laspy.open(self.paths[file], **arguments) as reader:
            if file in self.repaired:
                for block in blocks:
                    reader.seek(block.start)
                    yield _repair_times(_read(reader, block.stop - block.start))[0]
            else:
                yield from pulse_chunks(_read_ranges(reader, blocks), TRACE_BATCH)


def _read(reader, count):
    chunk = reader.read_points(count)
    if len(chunk) != count:
        raise ValueError(f"Truncated point data: expected {count:,} records, read {len(chunk):,}.")
    return chunk


def _read_ranges(reader, blocks):
    for block in blocks:
        if reader.points_read != block.start:
            reader.seek(block.start)
        # Each indexed block is independently bounded. The pulse reader holds the tail.
        yield _read(reader, block.stop - block.start)


# Compact incomplete echoes only. Complete pulses go straight to traversal.
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
