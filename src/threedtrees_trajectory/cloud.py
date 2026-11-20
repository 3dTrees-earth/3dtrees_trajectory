"""Bounded point streams in acquisition order for COPC and original order for LAS/LAZ."""
from __future__ import annotations

import mmap
from pathlib import Path
from shutil import disk_usage

import laspy
import numpy as np
from laspy.copc import Bounds, CopcReader, load_octree_for_query

CHUNK = 5_000_000
LAZ_FIELDS = (laspy.DecompressionSelection.base().decompress_z().decompress_gps_time()
              .decompress_scan_angle().decompress_point_source_id())
FORMAT = laspy.PointFormat(6)
TIME_BYTES = np.dtype(np.float64).itemsize
SORT_BYTES_PER_POINT = FORMAT.size + TIME_BYTES
SORT_LIMIT = 750_000_000_000


def is_copc(header) -> bool:
    return any(type(vlr).__name__ == "CopcInfoVlr" for vlr in header.vlrs)


def open_cloud(path: Path) -> laspy.LasReader:
    reader = laspy.open(path, decompression_selection=LAZ_FIELDS)
    if is_copc(reader.header):
        # Variable-sized COPC chunks cannot be decoded by laspy's parallel LAS reader.
        reader.laz_backend = laspy.LazBackend.Lazrs
    return reader


def _spatial_chunks(path: Path, bounds: Bounds | None, size: int):
    with CopcReader.open(path, decompression_selection=LAZ_FIELDS) as reader:
        if bounds is not None:
            bounds = bounds.ensure_3d(reader.header.mins, reader.header.maxs)
            lo = np.ceil((bounds.mins - reader.header.offsets) / reader.header.scales - 1e-7)
            hi = np.floor((bounds.maxs - reader.header.offsets) / reader.header.scales + 1e-7)
        nodes = load_octree_for_query(reader.source, reader.copc_info, reader.root_page, query_bounds=bounds)
        batch, count = [], 0
        for node in sorted(nodes, key=lambda n: n.offset):
            if batch and count + node.point_count > size:
                points = reader._fetch_and_decompress_points_of_nodes(batch)
                yield points if bounds is None else _in_bounds(points, ("X", "Y"), lo, hi)
                batch, count = [], 0
            if node.point_count > CHUNK:
                raise ValueError(f"{path.name}: COPC node exceeds the supported {CHUNK:,}-point stream budget.")
            batch.append(node)
            count += node.point_count
        if batch:
            points = reader._fetch_and_decompress_points_of_nodes(batch)
            yield points if bounds is None else _in_bounds(points, ("X", "Y"), lo, hi)


def _in_bounds(points, axes, lo, hi):
    keep = np.ones(len(points), bool)
    for k, axis in enumerate(axes):
        keep &= (getattr(points, axis) >= lo[k]) & (getattr(points, axis) <= hi[k])
    return points[keep]


def point_chunks(path: Path, size: int, bounds: Bounds | None = None, *, spatial: bool = False):
    if bounds is not None or spatial:
        yield from (p for p in _spatial_chunks(path, bounds, size) if len(p))
    else:
        with open_cloud(path) as reader:
            while reader.points_read < reader.header.point_count:
                expected = min(size, reader.header.point_count - reader.points_read)
                chunk = reader.read_points(expected)
                if len(chunk) != expected:
                    raise ValueError(f"Truncated point data: expected {expected:,} records, read {len(chunk):,}.")
                yield chunk


def pulse_chunks(chunks, size: int, *, time_tolerance: float = 0):
    """Keep adjacent returns together across reader batches, within the stream budget."""
    tail = None
    for chunk in chunks:
        if tail is not None:
            points = np.concatenate((tail.array, chunk.array))
            chunk = laspy.ScaleAwarePointRecord(points, chunk.point_format, chunk.scales, chunk.offsets)
        times = chunk.gps_time
        end = len(times) - 1
        lower = max(0, end - 64)
        while end > lower and abs(times[end] - times[end - 1]) <= time_tolerance:
            end -= 1
        if end == lower and lower:
            breaks = np.flatnonzero(np.abs(np.diff(times[:lower + 1])) > time_tolerance)
            end = int(breaks[-1]) + 1 if len(breaks) else 0
        if end:
            yield chunk[:end]
        elif len(chunk) >= size * 2:
            # An enormous same-time group cannot be joined as a pulse. Keep
            # the reader bounded so it can still diagnose an unusable clock.
            yield chunk
            tail = None
            continue
        tail = laspy.ScaleAwarePointRecord(chunk.array[end:].copy(), chunk.point_format, chunk.scales, chunk.offsets)
    if tail is not None and len(tail):
        yield tail


def _ordered(records) -> np.ndarray:
    if len(records) <= 1:
        return records
    flags = records["classification_flags"]
    rn = records["bit_fields"] & 15
    order = np.lexsort((rn, (flags >> 4) & 3, records["point_source_id"], records["gps_time"]))
    return records[order]


class TimeSortedCopc:
    """Sorted disk runs with a batch merge of disjoint, complete GPS-time groups."""

    def __init__(self, path: Path, folder: Path, size: int, bounds: Bounds | None = None, *, select=None,
                 observe=None):
        self.size = size
        with open_cloud(path) as reader:
            self.header = reader.header
            if not is_copc(self.header):
                raise ValueError(f"{path.name} is not a COPC point cloud.")
        if bounds is None and select is None:
            expected = self.header.point_count * SORT_BYTES_PER_POINT
            if expected > SORT_LIMIT or expected > disk_usage(folder).free:
                raise ValueError(f"{path.name}: COPC time sort needs {expected / 1e9:.1f} GB of temporary disk "
                                 "within the sort budget.")
        self.path = folder / "runs.bin"
        self.keys = folder / "times.bin"
        self.runs = []
        self.bytes_written = 0
        self.points_read = 0
        self.points_spooled = 0
        self.time_range = (np.inf, -np.inf)
        samples = []
        # Keep all sorted runs in one file. The number of runs can exceed the
        # process's file descriptor limit for billion-point acquisitions.
        with self.path.open("xb") as spool, self.keys.open("xb") as keys:
            for chunk in point_chunks(path, size, bounds, spatial=observe is not None):
                time = chunk.gps_time
                low, high = float(time.min()), float(time.max())
                if not np.isfinite([low, high]).all():
                    raise ValueError(f"{path.name}: GPS times must be finite.")
                self.time_range = (min(self.time_range[0], low), max(self.time_range[1], high))
                self.points_read += len(chunk)
                if observe is not None:
                    observe(0, chunk)
                if select is not None:
                    chunk = select(chunk)
                if not len(chunk):
                    continue
                if self.bytes_written + len(chunk) * SORT_BYTES_PER_POINT > SORT_LIMIT:
                    raise ValueError(f"{path.name}: COPC time sort exceeds the temporary disk budget.")
                data = _ordered(chunk.array if chunk.array.dtype == FORMAT.dtype() else _base_records(chunk))
                self.runs.append((self.points_spooled, len(data)))
                data.tofile(spool)
                # Searching an unaligned strided structured field makes NumPy
                # copy the entire field. Contiguous keys keep each search sparse.
                times = np.ascontiguousarray(data["gps_time"])
                times.tofile(keys)
                self.points_spooled += len(data)
                self.bytes_written += data.nbytes + times.nbytes
                samples.append(times[::max(1, size // 256)].copy())
        n_groups = max(1, -(-self.points_spooled // size))
        self.edges = (np.unique(np.quantile(np.concatenate(samples), np.arange(1, n_groups) / n_groups))
                      if self.runs else np.empty(0))

    def discard(self) -> None:
        self.path.unlink(missing_ok=True)
        self.keys.unlink(missing_ok=True)

    def _merge(self, arrays, times, starts, ends):
        count = int(np.sum(ends - starts))
        if count == 0:
            return
        if count > self.size:
            largest = int(np.argmax(ends - starts))
            values = times[largest]
            middle = values[(starts[largest] + ends[largest]) // 2]
            split = np.array([min(hi, max(lo, np.searchsorted(arr, middle, side="left")))
                              for arr, lo, hi in zip(times, starts, ends)], dtype=np.int64)
            if np.array_equal(split, starts):
                split = np.array([min(hi, max(lo, np.searchsorted(arr, middle, side="right")))
                                  for arr, lo, hi in zip(times, starts, ends)], dtype=np.int64)
            if np.array_equal(split, ends) or np.array_equal(split, starts):
                if count > max(4096, self.size):
                    raise ValueError("Too many COPC records share a GPS timestamp to reconstruct pulses safely.")
            else:
                yield from self._merge(arrays, times, starts, split)
                yield from self._merge(arrays, times, split, ends)
                return
        records = np.concatenate([arr[lo:hi] for arr, lo, hi in zip(arrays, starts, ends) if lo < hi])
        yield laspy.ScaleAwarePointRecord(_ordered(records), FORMAT, self.header.scales, self.header.offsets)

    def spooled_chunks(self):
        """Read each disk run without the merge when a scan only needs point values."""
        with self.path.open("rb") as stream:
            for _, size in self.runs:
                records = np.fromfile(stream, dtype=FORMAT.dtype(), count=size)
                if len(records) != size:
                    raise ValueError("A temporary COPC sort run was truncated.")
                yield laspy.ScaleAwarePointRecord(records, FORMAT, self.header.scales, self.header.offsets)

    def chunks(self):
        if not self.runs:
            return
        # Drop resident mapping pages after each batch while the file remains
        # available to later batches. This bounds resident sort memory.
        with self.path.open("rb") as spool, self.keys.open("rb") as key_file, \
             mmap.mmap(spool.fileno(), 0, access=mmap.ACCESS_READ) as records, \
             mmap.mmap(key_file.fileno(), 0, access=mmap.ACCESS_READ) as keys:
            arrays = [np.ndarray((size,), FORMAT.dtype(), buffer=records, offset=offset * FORMAT.size)
                      for offset, size in self.runs]
            times = [np.ndarray((size,), np.float64, buffer=keys, offset=offset * TIME_BYTES)
                     for offset, size in self.runs]
            starts = np.zeros(len(arrays), np.int64)
            try:
                for edge in (*self.edges, np.inf):
                    ends = np.array([np.searchsorted(t, edge) for t in times], dtype=np.int64)
                    for batch in self._merge(arrays, times, starts, ends):
                        yield batch
                        records.madvise(mmap.MADV_DONTNEED)
                        keys.madvise(mmap.MADV_DONTNEED)
                    starts = ends
            finally:
                times.clear()
                arrays.clear()


def _base_records(chunk):
    records = np.empty(len(chunk), FORMAT.dtype())
    for name in records.dtype.names:
        records[name] = chunk.array[name]
    return records
