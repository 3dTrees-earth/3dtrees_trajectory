"""Disk-backed echo routing. A tile receives every echo of a shot touching its buffered bounds."""
from __future__ import annotations

from itertools import groupby
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from .pulses import _time_hash
from .spool_codec import append_block, block_count, read_blocks

ECHO = np.dtype([("time", "<f8"), ("x", "<i4"), ("y", "<i4"), ("z", "<i4"),
                 ("source", "<u2"), ("returns", "u1"), ("flags", "u1")])
TILE_ECHO = np.dtype([(name if name != "source" else "file", dtype) for name, dtype in ECHO.descr])
START = np.uint8(128)
CHANNEL = np.uint8(3)
BATCH = 2_000_000
HASH_BITS = 6
HASH_FANOUT = 1 << HASH_BITS


def pack_echoes(points, keep=None):
    count = len(points) if keep is None else np.count_nonzero(keep)
    keep = slice(None) if keep is None else keep
    records = np.empty(count, ECHO)
    for raw, axis in zip(("X", "Y", "Z"), ("x", "y", "z")):
        records[axis] = getattr(points, raw)[keep]
    records["time"] = points.gps_time[keep]
    records["source"] = points.point_source_id[keep]
    records["returns"] = (np.asarray(points.return_number)[keep]
                          | (np.asarray(points.number_of_returns)[keep] << 4))
    records["flags"] = np.asarray(points.scanner_channel)[keep] if points.point_format.id >= 6 else 0
    return records


class EchoSource:
    """One compact recording of the source scan the path estimate already needs."""

    def __init__(self, directory):
        self.directory = directory
        self.bytes_written = 0

    def path(self, index):
        return self.directory / f"source-{index:05d}.bin"

    def observe(self, index, points):
        records = pack_echoes(points)
        _append(self.path(index), records)
        self.bytes_written += records.nbytes


def return_fields(records):
    packed = records["returns"]
    return packed & 15, packed >> 4


def singles(records):
    rn, nr = return_fields(records)
    return (rn <= 1) & (nr <= 1)


def raw_batches(path, *, contiguous=False, dtype=ECHO):
    """Read bounded chunks, optionally holding the trailing acquisition-time group."""
    tail = np.empty(0, dtype)
    with path.open("rb") as stream:
        while len(records := np.fromfile(stream, dtype, BATCH)):
            if not contiguous:
                yield records
                continue
            if len(tail):
                records = np.concatenate((tail, records))
            times = records["time"]
            last = times[-1]
            boundaries = np.flatnonzero(times[:-1] != last)
            stop = int(boundaries[-1]) + 1 if len(boundaries) else 0
            if len(records) - stop > 2 * BATCH:
                raise ValueError("Too many contiguous records share one GPS timestamp.")
            if stop:
                yield records[:stop]
            tail = records[stop:].copy()
    if len(tail):
        yield tail


def _time_key(time):
    # Positive and negative zero compare equal as timestamps.
    time = np.where(time == 0, 0.0, time)
    return _time_hash(time)


def _append(path, records):
    with path.open("ab") as stream:
        records.tofile(stream)


def _partitions(indices, count):
    order = np.argsort(indices, kind="stable")
    ends = np.cumsum(np.bincount(indices, minlength=count))
    starts = np.r_[0, ends[:-1]]
    for part in np.flatnonzero(ends > starts):
        yield int(part), order[starts[part]:ends[part]]


def _hash_partition(records, folder, shift):
    buckets = ((_time_key(records["time"]) >> shift) & (HASH_FANOUT - 1)).astype(np.intp)
    for part, index in _partitions(buckets, HASH_FANOUT):
        _append(folder / f"{part:02d}.bin", records[index])


def _hash_leaves(folder, shift):
    """Co-locate separated echoes by hash. Only tile spools are ordered by GPS time."""
    for path in sorted(folder.glob("*.bin")):
        if path.stat().st_size <= 2 * BATCH * ECHO.itemsize:
            yield np.fromfile(path, ECHO)
        else:
            if shift >= 64:
                raise ValueError("Too many separated echoes share one GPS timestamp to reconstruct pulses.")
            with TemporaryDirectory(dir=folder, prefix="partition-") as sub:
                for records in raw_batches(path):
                    _hash_partition(records, Path(sub), shift)
                path.unlink()
                yield from _hash_leaves(Path(sub), shift + HASH_BITS)
        path.unlink(missing_ok=True)


def group_multi(records, *, acquisition_order, hashes=None):
    """Return offsets, sizes and completeness for original runs or keyed groups."""
    rn, nr = return_fields(records)
    key = records["time"] if hashes is None else hashes
    source, channel = records["source"], records["flags"] & CHANNEL
    same = (key[1:] == key[:-1]) & (source[1:] == source[:-1]) & (channel[1:] == channel[:-1])
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
    return starts, sizes, complete


class TileSpools:
    def __init__(self, folder, grid, layout, time, points, *, compress=False):
        self.folder = folder
        self.grid = grid
        self.layout = layout
        starts = np.array([[s.start for s in tile.extent[:2]] for tile in layout])
        stops = np.array([[s.stop for s in tile.extent[:2]] for tile in layout])
        self.lower = grid.origin[:2] + starts * grid.voxel
        self.upper = grid.origin[:2] + stops * grid.voxel
        windows = min(1024, max(1, -(-points // (len(layout) * BATCH))))
        self.edges = np.linspace(time[0], time[-1], windows + 1)[1:-1]
        self.echoes = np.zeros(len(layout), np.int64)
        self.bytes_written = 0
        self.compress = compress

    def _write_window(self, directory, file, window, records):
        if self.compress:
            self.bytes_written += append_block(directory / f"{window:04d}-{file:05d}.spool", records)
        else:
            _append(directory / f"{window:04d}.bin", records)
            self.bytes_written += records.nbytes

    def _write(self, file, tile, records):
        if not len(records):
            return
        directory = self.folder / f"tile-{tile:05d}"
        directory.mkdir(exist_ok=True)
        records = records.view(TILE_ECHO)
        records["file"] = file
        if not len(self.edges):
            self._write_window(directory, file, 0, records)
        else:
            windows = np.searchsorted(self.edges, records["time"])
            for window, index in _partitions(windows, len(self.edges) + 1):
                self._write_window(directory, file, window, records[index])
        self.echoes[tile] += len(records)

    def route(self, file, records, scales, offsets, starts=None, sizes=None):
        if not len(records):
            return
        if starts is None:
            records["flags"] |= START
        else:
            records["flags"][starts] |= START
        if len(self.layout) == 1 and self.layout[0].shape[:2] == self.grid.shape[:2]:
            self._write(file, 0, records)
            return
        x, y = (records[axis] * scale + offset for axis, scale, offset in zip(("x", "y"), scales, offsets))
        lower, upper = np.array([x.min(), y.min()]), np.array([x.max(), y.max()])
        candidates = np.flatnonzero(((self.lower <= upper) & (self.upper >= lower)).all(axis=1))
        for tile in candidates:
            lo, hi = self.lower[tile], self.upper[tile]
            keep = (lo[0] <= x) & (x <= hi[0]) & (lo[1] <= y) & (y <= hi[1])
            if starts is not None:
                keep = np.repeat(np.logical_or.reduceat(keep, starts), sizes)
            self._write(file, tile, records[keep])

    def stage_file(self, file, path, scales, offsets, copc, time_bounds=None):
        recovered = unpaired = 0
        if not path.is_file():
            return recovered, unpaired
        with TemporaryDirectory(dir=self.folder, prefix="echo-hash-") as folder:
            folder = Path(folder)
            for records in raw_batches(path, contiguous=not copc):
                if time_bounds is not None:
                    low, high = time_bounds
                    keep = (records["time"] >= low) & (records["time"] <= high)
                    if not keep.all():
                        records = records[keep]
                one = singles(records)
                self.route(file, records[one], scales, offsets)
                multi = records[~one]
                if not len(multi):
                    continue
                if not copc:
                    starts, sizes, complete = group_multi(multi, acquisition_order=True)
                    keep = np.repeat(complete, sizes)
                    stable = multi[keep]
                    stable_size = sizes[complete]
                    self.route(file, stable, scales, offsets, np.cumsum(stable_size) - stable_size, stable_size)
                    multi = multi[~keep]
                if len(multi):
                    _hash_partition(multi, folder, 0)
            path.unlink()
            for records in _hash_leaves(folder, HASH_BITS):
                hashes = _time_key(records["time"])
                order = np.lexsort((return_fields(records)[0], records["flags"], records["source"], hashes))
                records, hashes = records[order], hashes[order]
                _, sizes, complete = group_multi(records, acquisition_order=False, hashes=hashes)
                good = sizes[complete]
                self.route(file, records[np.repeat(complete, sizes)], scales, offsets, np.cumsum(good) - good, good)
                unpaired += int(sizes[~complete].sum())
                if not copc:
                    recovered += int(good.sum())
        return recovered, unpaired


def _shot_batches(batches):
    tail = np.empty(0, TILE_ECHO)
    for records in batches:
        if len(tail):
            records = np.concatenate((tail, records))
        stop = int(np.flatnonzero(records["flags"] & START)[-1])
        if stop:
            yield records[:stop]
        tail = records[stop:].copy()
    if len(tail):
        yield tail


def _window_count(path):
    if path.suffix == ".spool":
        return block_count(path, TILE_ECHO)
    return path.stat().st_size // TILE_ECHO.itemsize


def _window_batches(paths):
    for path in paths:
        read = read_blocks if path.suffix == ".spool" else raw_batches
        yield from read(path, dtype=TILE_ECHO)


def _sorted_window(paths, scratch):
    """Sort one bounded time window, splitting dense windows on disk if needed."""
    count = sum(map(_window_count, paths))
    if not count:
        return
    if count <= 2 * BATCH:
        records = np.concatenate(list(_window_batches(paths)))
        yield records[np.argsort(records["time"], kind="stable")]
        return
    low, high = np.inf, -np.inf
    for records in _window_batches(paths):
        low, high = min(low, records["time"].min()), max(high, records["time"].max())
    if low == high:
        yield from _shot_batches(_window_batches(paths))
        return
    middle = low / 2 + high / 2
    if middle == high:
        middle = low
    with TemporaryDirectory(dir=scratch, prefix="dense-window-") as folder:
        before, after = (Path(folder) / name for name in ("before.bin", "after.bin"))
        for records in _window_batches(paths):
            keep = records["time"] <= middle
            if keep.any():
                _append(before, records[keep])
            if not keep.all():
                _append(after, records[~keep])
        yield from _sorted_window((before,), scratch)
        yield from _sorted_window((after,), scratch)


def tile_batches(directory, scratch):
    paths = sorted((*directory.glob("*.bin"), *directory.glob("*.spool")))
    for _, window_paths in groupby(paths, key=lambda path: path.stem.split("-", 1)[0]):
        yield from _sorted_window(tuple(window_paths), scratch)
