"""Bounded point streams in acquisition order for COPC and original order for LAS/LAZ."""
from __future__ import annotations

import mmap
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory

import laspy
import numpy as np
from laspy.copc import Bounds, CopcReader, load_octree_for_query

CHUNK = 5_000_000
LAZ_FIELDS = (laspy.DecompressionSelection.base().decompress_z().decompress_gps_time()
              .decompress_scan_angle().decompress_point_source_id())
FORMAT = laspy.PointFormat(6)
SORT_BYTES_PER_POINT = FORMAT.size
SORT_LIMIT = 750_000_000_000


def is_copc(header) -> bool:
    return any(type(vlr).__name__ == "CopcInfoVlr" for vlr in header.vlrs)


def open_cloud(path: Path) -> laspy.LasReader:
    reader = laspy.open(path, decompression_selection=LAZ_FIELDS)
    if is_copc(reader.header):
        # Variable-sized COPC chunks cannot be decoded by laspy's parallel LAS reader.
        reader.laz_backend = laspy.LazBackend.Lazrs
    return reader


def _spatial_chunks(path: Path, bounds: Bounds, size: int):
    with CopcReader.open(path, decompression_selection=LAZ_FIELDS) as reader:
        bounds = bounds.ensure_3d(reader.header.mins, reader.header.maxs)
        nodes = load_octree_for_query(reader.source, reader.copc_info, reader.root_page, query_bounds=bounds)
        axes = ("X", "Y")
        lo = np.ceil((bounds.mins - reader.header.offsets) / reader.header.scales - 1e-7)
        hi = np.floor((bounds.maxs - reader.header.offsets) / reader.header.scales + 1e-7)
        batch, count = [], 0
        for node in sorted(nodes, key=lambda n: n.offset):
            if batch and count + node.point_count > size:
                yield _in_bounds(reader._fetch_and_decompress_points_of_nodes(batch), axes, lo, hi)
                batch, count = [], 0
            if node.point_count > CHUNK:
                raise ValueError(f"{path.name}: COPC node exceeds the supported {CHUNK:,}-point stream budget.")
            batch.append(node)
            count += node.point_count
        if batch:
            yield _in_bounds(reader._fetch_and_decompress_points_of_nodes(batch), axes, lo, hi)


def _in_bounds(points, axes, lo, hi):
    keep = np.ones(len(points), bool)
    for k, axis in enumerate(axes):
        keep &= (getattr(points, axis) >= lo[k]) & (getattr(points, axis) <= hi[k])
    return points[keep]


def point_chunks(path: Path, size: int, bounds: Bounds | None = None):
    if bounds is not None:
        yield from (p for p in _spatial_chunks(path, bounds, size) if len(p))
    else:
        with open_cloud(path) as reader:
            while reader.points_read < reader.header.point_count:
                expected = min(size, reader.header.point_count - reader.points_read)
                chunk = reader.read_points(expected)
                if len(chunk) != expected:
                    raise ValueError(f"Truncated point data: expected {expected:,} records, read {len(chunk):,}.")
                yield chunk


def pulse_chunks(chunks, size: int):
    """Keep every contiguous timestamp group in one batch, including on a chunk boundary."""
    tail = None
    for chunk in chunks:
        if tail is not None:
            points = np.concatenate((tail.array, chunk.array))
            chunk = laspy.ScaleAwarePointRecord(points, chunk.point_format, chunk.scales, chunk.offsets)
        times = chunk.gps_time
        end = len(times) - 1
        while end > 0 and times[end - 1] == times[-1]:
            end -= 1
        if end:
            yield chunk[:end]
        elif len(chunk) >= size * 2:
            raise ValueError("Too many records share one GPS timestamp to reconstruct pulses within the stream budget.")
        tail = chunk[end:].copy()
    if tail is not None:
        if len(tail):
            yield laspy.ScaleAwarePointRecord(tail.array, tail.point_format, chunk.scales, chunk.offsets)


def _sorted_records(points) -> np.ndarray:
    rn = np.asarray(points.return_number)
    order = np.lexsort((rn, np.asarray(points.scanner_channel), points.point_source_id, points.gps_time))
    out = np.empty(len(order), FORMAT.dtype())
    for name in out.dtype.names:
        if name in ("classification", "intensity", "user_data"):
            out[name] = 0
        else:
            out[name] = points.array[name][order]
    return out


def _ordered(records) -> np.ndarray:
    if len(records) <= 1:
        return records
    flags = records["classification_flags"]
    rn = records["bit_fields"] & 15
    order = np.lexsort((rn, (flags >> 4) & 3, records["point_source_id"], records["gps_time"]))
    return records[order]


class TimeSortedCopc:
    """Sorted disk runs with a batch merge of disjoint, complete GPS-time groups."""

    def __init__(self, path: Path, folder: Path, size: int, bounds: Bounds | None = None):
        self.size = size
        with open_cloud(path) as reader:
            self.header = reader.header
            if not is_copc(self.header):
                raise ValueError(f"{path.name} is not a COPC point cloud.")
        self.path = folder / "runs.bin"
        self.runs = []
        self.bytes_written = 0
        self.points_read = 0
        samples = []
        # Keep all sorted runs in one file. The number of runs can exceed the
        # process's file descriptor limit for billion-point acquisitions.
        with self.path.open("xb") as spool:
            for chunk in point_chunks(path, size, bounds):
                if not np.isfinite(chunk.gps_time).all():
                    raise ValueError(f"{path.name}: GPS times must be finite.")
                if self.bytes_written + len(chunk) * SORT_BYTES_PER_POINT > SORT_LIMIT:
                    raise ValueError(f"{path.name}: COPC time sort exceeds the temporary disk budget.")
                data = _sorted_records(chunk)
                self.runs.append((self.bytes_written, len(data)))
                data.tofile(spool)
                self.points_read += len(data)
                self.bytes_written += data.nbytes
                samples.append(data["gps_time"][::max(1, size // 256)].copy())
        n_groups = max(1, -(-self.points_read // size))
        self.edges = (np.unique(np.quantile(np.concatenate(samples), np.arange(1, n_groups) / n_groups))
                      if self.runs else np.empty(0))

    def discard(self) -> None:
        self.path.unlink(missing_ok=True)

    def _merge(self, arrays, starts, ends):
        count = int(np.sum(ends - starts))
        if count == 0:
            return
        if count > self.size:
            largest = int(np.argmax(ends - starts))
            values = arrays[largest]["gps_time"]
            middle = values[(starts[largest] + ends[largest]) // 2]
            split = np.array([min(hi, max(lo, np.searchsorted(arr["gps_time"], middle, side="left")))
                              for arr, lo, hi in zip(arrays, starts, ends)], dtype=np.int64)
            if np.array_equal(split, starts):
                split = np.array([min(hi, max(lo, np.searchsorted(arr["gps_time"], middle, side="right")))
                                  for arr, lo, hi in zip(arrays, starts, ends)], dtype=np.int64)
            if np.array_equal(split, ends) or np.array_equal(split, starts):
                if count > max(4096, self.size):
                    raise ValueError("Too many COPC records share a GPS timestamp to reconstruct pulses safely.")
            else:
                yield from self._merge(arrays, starts, split)
                yield from self._merge(arrays, split, ends)
                return
        records = np.concatenate([arr[lo:hi] for arr, lo, hi in zip(arrays, starts, ends) if lo < hi])
        yield laspy.ScaleAwarePointRecord(_ordered(records), FORMAT, self.header.scales, self.header.offsets)

    def chunks(self):
        if not self.runs:
            return
        # Drop resident mapping pages after each batch while the file remains
        # available to later batches. This bounds resident sort memory.
        with self.path.open("rb") as spool, mmap.mmap(spool.fileno(), 0, access=mmap.ACCESS_READ) as mapping:
            arrays = [np.ndarray((size,), FORMAT.dtype(), buffer=mapping, offset=offset)
                      for offset, size in self.runs]
            starts = np.zeros(len(arrays), np.int64)
            try:
                for edge in (*self.edges, np.inf):
                    ends = np.array([np.searchsorted(a["gps_time"], edge) for a in arrays], dtype=np.int64)
                    yield from self._merge(arrays, starts, ends)
                    mapping.madvise(mmap.MADV_DONTNEED)
                    starts = ends
            finally:
                arrays.clear()


@contextmanager
def sorted_copc(paths: list[Path], scratch: Path, size: int = CHUNK, resources=None):
    """Prepare a reusable stream only for a single merged COPC; LAZ does no disk sort."""
    if not paths:
        yield None
        return
    try:
        with open_cloud(paths[0]) as reader:
            first_is_copc = is_copc(reader.header)
    except (OSError, laspy.errors.LaspyException, ValueError, OverflowError, EOFError):
        # The pulse reader validates every header and formats the refusal.
        first_is_copc = False
    if not first_is_copc:
        yield None
        return
    if len(paths) != 1:
        raise ValueError("Use one merged COPC or a list of per-file LAS/LAZ inputs.")
    with TemporaryDirectory(dir=scratch, prefix=".copc-sort-") as folder:
        if resources is None:
            result = TimeSortedCopc(paths[0], Path(folder), size)
        else:
            with resources.stage("read/time-sort") as stats:
                result = TimeSortedCopc(paths[0], Path(folder), size)
                stats.update(bytes_read=paths[0].stat().st_size, bytes_written=result.bytes_written)
        yield result
