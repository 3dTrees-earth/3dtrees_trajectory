"""Read a point cloud into pulses: the first and last return of every laser shot.

A pulse is a run of consecutive records with the same GPS time and point source whose return numbers count up by one.
The time alone is not enough: multi-beam scanners give a dozen pulses one time stamp, and some writers fill
NumberOfReturns wrongly, so the file order is what ties a pulse together. Some software stores the points out of
acquisition order (2587, LiDAR360); records of a multi-return pulse that end up alone are then collected and joined
by GPS time and point source after reading.

Coordinates use float32 relative to a common kilometre-aligned origin. The default budgets keep approximately
MAX_MULTI multi-return and MAX_SINGLE single-return pulses with a deterministic timestamp hash. A cloud whose
first million points hold few multi-return pulses uses the MAX_SINGLE_ONLY single-return budget. ``full``
keeps every pulse.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import hypot
from pathlib import Path
from struct import error as StructError
from tempfile import TemporaryDirectory

import laspy
import numpy as np
from lazrs import LazrsError
from pyproj import CRS
from pyproj.exceptions import CRSError

from .cloud import (
    CHUNK, COPC_BATCH_BYTES, COPC_BATCH_NODES, FORMAT, TimeSortedCopc,
    _base_records, _ordered, is_copc, open_cloud, point_chunks, pulse_chunks,
)

from .coverage import UnusableTime

TIME_TOLERANCE_S = 1e-6
PROBE = 1_000_000  # points read first to tell a single-return cloud
MAX_MULTI = 40_000_000
MAX_SINGLE = 10_000_000  # single returns next to multi-return pulses
MAX_SINGLE_ONLY = 40_000_000  # clouds without multi-return pulses, where single returns are all there is
MAX_STRAY = 40_000_000  # separate limit for records that need a global sort to reconstruct pulses


@dataclass
class Pulses:
    """One row per pulse; coordinates relative to ``origin``. ``first`` and ``last`` are equal for single returns."""

    origin: np.ndarray  # (3,) float64
    time: np.ndarray  # float64, the cloud's own clock
    first: np.ndarray  # (n, 3) float32
    last: np.ndarray  # (n, 3) float32
    returns: np.ndarray  # uint8, records in the pulse
    source: np.ndarray  # uint16, PointSourceId
    file: np.ndarray  # uint16, index of the input file
    scan_angle: np.ndarray  # float32, degrees off nadir as written in the file (0 when absent)
    has_time: bool
    has_scan_angle: bool
    multi_share: float  # share of multi-return pulses in the whole cloud, before thinning

    def __len__(self) -> int:
        return len(self.time)



@dataclass
class CloudInfo:
    crs: str  # "EPSG:xxxx" or "unknown"
    time_system: str  # gps_adjusted | gps_sow | none
    points: int
    files: list[str]
    kept_every: tuple[int, int] = (1, 1)  # (multi-return, single-return) pulse strides
    points_read: int = 0
    read_workers: int = 1
    incomplete_multi_pulses: int = 0
    multi_pulses: int = 0
    joined_out_of_order_points: int = 0
    fragment_echoes: int = 0
    complete_multi_echoes: int = 0
    unpaired_multi_echoes: int = 0
    population_counts_estimated: bool = False
    return_modes: list[dict] = field(default_factory=list)
    has_gps_time: bool | None = None
    single_return_files: list[int] = field(default_factory=list)

    @classmethod
    def from_paths(cls, paths: list[Path]) -> CloudInfo:
        return cls("unknown", "none", 0, [str(path) for path in paths])


def _crs(header: laspy.LasHeader) -> CRS | None:
    try:
        return header.parse_crs()
    except (CRSError, ValueError):
        return None  # a malformed CRS record leaves the CRS unknown


def is_local_geocentric(crs: CRS, centre) -> bool:
    """Recognise Standardization's CRS placeholder using Raster's distance cutoff."""
    return crs.is_geocentric and hypot(*centre) < 1e6


def cloud_crs(crs: CRS | None, lower: np.ndarray, upper: np.ndarray) -> CRS | None:
    if crs is not None and crs.is_geocentric:
        if is_local_geocentric(crs, 0.5 * lower + 0.5 * upper):
            return None
        raise ValueError("Earth-centred coordinates: reproject to a projected CRS in metres")
    return crs


def _time_system(header: laspy.LasHeader, names: set[str]) -> str:
    if "gps_time" not in names:
        return "none"
    return "gps_adjusted" if header.global_encoding.gps_time_type == laspy.header.GpsTimeType.STANDARD else "gps_sow"


class _Store:
    """Growing column lists, concatenated once at the end."""

    def __init__(self, *names: str):
        self.columns = {name: [] for name in names}

    def add(self, **values: np.ndarray) -> None:
        for name, value in values.items():
            if len(value):
                self.columns[name].append(value)

    def extend(self, other: _Store) -> None:
        for name, parts in self.columns.items():
            parts.extend(other.columns[name])

    def array(self, name: str, dtype, width: int = 0) -> np.ndarray:
        parts = self.columns.pop(name)
        if len(parts) == 1:
            return parts[0].astype(dtype, copy=False)
        if parts:
            return np.concatenate(parts, dtype=dtype, casting="unsafe")
        return np.empty((0, width) if width else 0, dtype)


def _join_stray(stray: _Store, pulses: _Store) -> tuple[int, int]:
    """Recover separated returns only when their complete, unambiguous sequence is present."""
    t, xyz, number = stray.array("time", np.float64), stray.array("xyz", np.float32, 3), stray.array("number", np.uint8)
    source, file, angle = stray.array("source", np.uint16), stray.array("file", np.uint16), stray.array("angle", np.float32)
    total, channel = stray.array("total", np.uint8), stray.array("channel", np.uint8)
    order = np.lexsort((number, channel, t, source, file))
    t, source, file, channel = t[order], source[order], file[order], channel[order]
    same = (t[1:] == t[:-1]) & (source[1:] == source[:-1]) & (file[1:] == file[:-1]) & (channel[1:] == channel[:-1])
    # A repeated/skipped return number invalidates the complete group.
    sequential = np.diff(number[order].astype(np.int16)) == 1
    start = np.flatnonzero(np.r_[True, ~same])
    end = np.r_[start[1:], len(t)] - 1
    broken = np.logical_or.reduceat(np.r_[same & ~sequential, False], start)
    good = (~broken) & (number[order[start]] == 1) & (end - start + 1 == total[order[start]])
    start, end = start[good], end[good]
    pulses.add(time=t[start], first=xyz[order[start]], last=xyz[order[end]],
               returns=(end - start + 1).astype(np.uint8), source=source[start], file=file[start], angle=angle[order[start]])
    return len(start), int((end - start + 1).sum())


def _sample(time: np.ndarray, stride: int) -> np.ndarray:
    if stride == 1:
        return np.ones(len(time), bool)
    # A deterministic hash of the timestamp, independent of COPC order, file
    # boundaries, and whether a return had to be joined out of order.
    return _time_hash(time) % stride == 0


def _time_hash(time: np.ndarray) -> np.ndarray:
    key = time.view(np.uint64).copy()
    key ^= key >> 33
    key *= np.uint64(0xff51afd7ed558ccd)
    key ^= key >> 33
    return key


class _CopcSample:
    """Keep all echoes at sampled times and the earliest records for the stride probe."""

    def __init__(self, points, max_multi, max_single, full):
        self.strides = tuple(sorted({1 if full else max(1, -(-points // limit))
                                     for limit in (max_multi, max_single, MAX_SINGLE_ONLY)}))
        self.probe = np.empty(0, FORMAT.dtype())

    def __call__(self, points):
        time = points.gps_time
        if len(self.probe) >= PROBE:
            earliest = time <= self.probe["gps_time"].max()
        else:
            earliest = np.ones(len(time), bool)
        if earliest.any():
            possible = np.concatenate((self.probe, _base_records(points[earliest])))
            if len(possible) > PROBE:
                times = possible["gps_time"]
                threshold = np.partition(times, PROBE - 1)[PROBE - 1]
                possible = possible[times <= threshold]
                if len(possible) > PROBE:
                    possible = _ordered(possible)[:PROBE].copy()
            self.probe = possible
        if 1 in self.strides:
            return points
        hashes = _time_hash(time)
        keep = np.zeros(len(time), bool)
        for stride in self.strides:
            keep |= hashes % stride == 0
        return points[keep]

    def looks_single(self, prepared):
        ordered = _ordered(self.probe)[:PROBE]
        self.probe = None
        return _looks_single_return(laspy.ScaleAwarePointRecord(ordered, FORMAT, prepared.header.scales,
                                                               prepared.header.offsets))


def _looks_single_return(chunk: laspy.ScaleAwarePointRecord) -> bool:
    """Under 1 % multi-return records by declaration or consecutive return numbers."""
    probe = chunk[:PROBE]
    rn = np.asarray(probe.return_number)
    by_field = float((np.asarray(probe.number_of_returns) > 1).mean())
    by_order = float((rn[1:] == rn[:-1] + 1).mean()) if len(rn) > 1 else 0.0
    return max(by_field, by_order) < 0.01


def _coordinates(chunk: laspy.ScaleAwarePointRecord, index: np.ndarray, shift: np.ndarray) -> np.ndarray:
    xyz = np.empty((len(index), 3), np.float32)
    for axis, raw in enumerate((chunk.X, chunk.Y, chunk.Z)):
        xyz[:, axis] = raw[index] * chunk.scales[axis] + shift[axis]
    return xyz


def _angles(chunk: laspy.ScaleAwarePointRecord, index: np.ndarray) -> np.ndarray:
    if chunk.point_format.id >= 6:
        return chunk.scan_angle[index].astype(np.float32) * 0.006
    return chunk.scan_angle_rank[index].astype(np.float32)


def _collect_chunk(chunk: laspy.ScaleAwarePointRecord, time: np.ndarray, file: int, shift: np.ndarray,
                   strides: tuple[int, int], counts: np.ndarray, pulses: _Store, stray: _Store, info: CloudInfo,
                   copc: bool, sampled_stats: bool = False) -> int:
    multi_every, single_every = strides
    rn, nr = np.asarray(chunk.return_number), np.asarray(chunk.number_of_returns)
    source = chunk.point_source_id
    channel = np.asarray(chunk.scanner_channel) if chunk.point_format.id >= 6 else np.zeros(len(chunk), np.uint8)
    same = (time[1:] == time[:-1]) & (source[1:] == source[:-1]) & (channel[1:] == channel[:-1])
    joined = same & (rn[1:] == rn[:-1] + 1)
    start = np.flatnonzero(np.r_[True, ~joined])
    end = np.r_[start[1:], len(chunk)] - 1
    single = (start == end) & (rn[start] <= 1) & (nr[start] <= 1)
    complete = ((start != end) & (rn[start] == 1) & ((nr[start] <= 1) | (rn[end] == nr[start])))
    if copc:
        # Time/source/channel only identify one multi-return shot when each
        # declared return has a unique number in the group.
        repeated = same & (rn[1:] == rn[:-1]) & ((nr[1:] > 1) | (nr[:-1] > 1))
        if repeated.any():
            groups = np.cumsum(np.r_[0, ~same], dtype=np.int64)
            ambiguous = np.zeros(groups[-1] + 1, bool)
            ambiguous[groups[:-1][repeated]] = True
            complete &= ~ambiguous[groups[start]]
    incomplete = ~single & ~complete
    sizes = end - start + 1
    multi_sample = _sample(time[start], multi_every)
    single_sample = _sample(time[start], single_every)
    if sampled_stats:
        # The COPC sample contains whole timestamp groups. Estimate population
        # totals using each pulse category's own inclusion fraction.
        counts += [int((complete & multi_sample).sum()) * multi_every,
                   int((single & single_sample).sum()) * single_every]
        info.incomplete_multi_pulses += int((incomplete & multi_sample).sum()) * multi_every
        info.multi_pulses += int((complete & multi_sample).sum()) * multi_every
        info.fragment_echoes += int(sizes[incomplete & multi_sample].sum()) * multi_every
        info.complete_multi_echoes += int(sizes[complete & multi_sample].sum()) * multi_every
    else:
        counts += [int(complete.sum()), int(single.sum())]
        info.incomplete_multi_pulses += int(incomplete.sum())
        info.multi_pulses += int(complete.sum())
        info.fragment_echoes += int(sizes[incomplete].sum())
        info.complete_multi_echoes += int(sizes[complete].sum())
    keep = (complete & multi_sample) | (single & single_sample)
    lone_idx = None
    # Retain every record of any incomplete LAZ pulse so it can be repaired.
    if not copc and incomplete.any():
        lone_idx = np.flatnonzero(np.repeat(incomplete & multi_sample, sizes))
    start, end = start[keep], end[keep]
    first = _coordinates(chunk, start, shift)
    multi = start != end
    last = first
    if np.any(multi):
        last = first.copy()
        last[multi] = _coordinates(chunk, end[multi], shift)
    pulses.add(time=time[start], first=first, last=last, returns=(end - start + 1).astype(np.uint8),
               source=source[start], file=np.full(len(start), file, dtype=np.uint16), angle=_angles(chunk, start))
    if lone_idx is not None and len(lone_idx):
        stray.add(time=time[lone_idx], xyz=_coordinates(chunk, lone_idx, shift), number=rn[lone_idx], total=nr[lone_idx],
                  channel=channel[lone_idx], source=source[lone_idx], file=np.full(len(lone_idx), file, dtype=np.uint16),
                  angle=_angles(chunk, lone_idx))
        return len(lone_idx)
    return 0


def _repair_times(chunk):
    """Snap only a complete, unique sequence within 1 µs on the same source and scanner channel."""
    time = chunk.gps_time
    rn, nr = np.asarray(chunk.return_number), np.asarray(chunk.number_of_returns)
    source = chunk.point_source_id
    channel = np.asarray(chunk.scanner_channel) if chunk.point_format.id >= 6 else np.zeros(len(chunk), np.uint8)
    multi = np.flatnonzero(nr > 1)
    if len(multi) < 2:
        return chunk, 0
    order = multi[np.lexsort((time[multi], channel[multi], source[multi]))]
    t, src, channels = time[order], source[order], channel[order]
    adjacent = (src[1:] == src[:-1]) & (channels[1:] == channels[:-1]) & (np.diff(t) <= TIME_TOLERANCE_S)
    starts = np.flatnonzero(np.r_[True, ~adjacent])
    sizes = np.diff(np.r_[starts, len(t)])
    ends = starts + sizes - 1
    numbers = np.bitwise_or.reduceat(np.left_shift(np.uint32(1), rn[order]), starts)
    totals = nr[order]
    expected = totals[starts]
    good = ((sizes == expected) & (t[ends] > t[starts]) & (t[ends] - t[starts] <= TIME_TOLERANCE_S)
            & (np.minimum.reduceat(totals, starts) == np.maximum.reduceat(totals, starts))
            & (numbers == (np.left_shift(np.uint32(1), expected + 1) - 2)))
    if not good.any():
        return chunk, 0
    keep = np.repeat(good, sizes)
    chunk = laspy.ScaleAwarePointRecord(chunk.array.copy(), chunk.point_format, chunk.scales, chunk.offsets)
    chunk.gps_time[order[keep]] = np.repeat(t[starts[good]], sizes[good])
    return chunk, int(good.sum())


def _single_returns(chunk):
    chunk = laspy.ScaleAwarePointRecord(chunk.array.copy(), chunk.point_format, chunk.scales, chunk.offsets)
    chunk.return_number[:] = 1
    chunk.number_of_returns[:] = 1
    return chunk


_PARTIAL_COUNTS = ("incomplete_multi_pulses", "multi_pulses", "joined_out_of_order_points",
                   "fragment_echoes", "complete_multi_echoes", "unpaired_multi_echoes")


class _FilePulses:
    """Keep one file's original and recovered endpoints until its return mode is known."""

    def __init__(self, paths):
        self.info = CloudInfo.from_paths(paths)
        self.counts = np.zeros(2, np.int64)
        self.pulses = _Store("time", "first", "last", "returns", "source", "file", "angle")
        self.recovered = _Store("time", "first", "last", "returns", "source", "file", "angle")
        self.stray = _Store("time", "xyz", "number", "source", "file", "angle", "total", "channel")
        self.stray_count = 0
        self.stray_limit_exceeded = False

    def collect(self, chunk, index, shift, strides, grouped, sampled=False):
        self.stray_count += _collect_chunk(chunk, chunk.gps_time, index, shift, strides, self.counts, self.pulses,
                                           self.stray, self.info, grouped or self.stray_limit_exceeded, sampled)
        if self.stray_count > MAX_STRAY and not self.stray_limit_exceeded:
            self.stray.columns.clear()
            self.stray_limit_exceeded = True

    def append_batch(self, part):
        self.stray_count += part.stray_count
        self.stray_limit_exceeded = part.stray_limit_exceeded or self.stray_count > MAX_STRAY
        self.counts += part.counts
        self.pulses.extend(part.pulses)
        if not self.stray_limit_exceeded:
            self.stray.extend(part.stray)
        for key in _PARTIAL_COUNTS:
            setattr(self.info, key, getattr(self.info, key) + getattr(part.info, key))

    def finish(self, stride):
        info = self.info
        info.joined_out_of_order_points = self.stray_count
        recovered = 0
        if self.stray_count and not self.stray_limit_exceeded:
            count, echoes = _join_stray(self.stray, self.recovered)
            self.counts[0] += count * stride
            recovered = min(info.fragment_echoes, echoes * stride)
        info.unpaired_multi_echoes = info.fragment_echoes - recovered
        if self.stray_limit_exceeded:
            return False
        if not info.unpaired_multi_echoes:
            return True
        return self.counts[0] > 0 and info.unpaired_multi_echoes <= max(200, info.complete_multi_echoes + recovered)

    def has_angles(self):
        return any(np.any(np.isfinite(part) & (np.abs(part) > 0) & (np.abs(part) < 90))
                   for part in self.pulses.columns["angle"])


def _fallback_chunks(path, mode, copc):
    if mode == "single":
        for chunk in point_chunks(path, CHUNK, spatial=copc):
            yield _single_returns(chunk)
    else:
        # LAS/LAZ acquisition order lets consecutive chunks share a bounded
        # tail. For spatially ordered COPC only exact joins are trusted.
        for chunk in pulse_chunks(point_chunks(path, CHUNK), CHUNK, time_tolerance=TIME_TOLERANCE_S):
            yield _repair_times(chunk)[0]


def _strides(points, full, max_multi, max_single):
    return (1, 1) if full else (max(1, -(-points // max_multi)), max(1, -(-points // max_single)))


def read_pulses(paths: list[Path], max_multi: int = MAX_MULTI, max_single: int = MAX_SINGLE,
                full: bool = False, *, info: CloudInfo | None = None, prepared=None,
                scratch: Path | None = None, resources=None, scan_index=None,
                threads: int | None = None, _selection=None) -> tuple[Pulses, CloudInfo]:
    """Build pulse endpoints from strips or from a COPC whose sampled echoes are ordered by time."""
    info = CloudInfo.from_paths(paths) if info is None else info
    if len(paths) > np.iinfo(np.uint16).max + 1:
        raise ValueError("At most 65,536 input files can be combined.")
    origins = []
    maxima = []
    copc_paths = []
    for path in paths:
        try:
            if not path.is_file():
                raise ValueError("Input must be an existing regular file.")
            with open_cloud(path, read_evlrs=False) as reader:
                if reader.header.point_count:
                    origins.append(reader.header.mins)
                    maxima.append(reader.header.maxs)
                info.points += reader.header.point_count
                if is_copc(reader.header):
                    copc_paths.append(path)
        except (OSError, laspy.errors.LaspyException, StructError, ValueError, OverflowError, EOFError) as error:
            raise ValueError(f"{path}: unable to read the LAS/LAZ header ({error}).") from error
    if not info.points:
        raise ValueError("The cloud has no points.")

    if copc_paths and (len(paths) != 1 or len(copc_paths) != 1):
        raise ValueError("Use one merged COPC or a list of per-file LAS/LAZ inputs.")
    if scan_index is not None:
        scan_index.copc = bool(copc_paths)
    if copc_paths and prepared is None:
        from contextlib import nullcontext

        path = copc_paths[0]
        declared_points = info.points
        selection = _CopcSample(declared_points, max_multi, max_single, full)
        info.points = 0
        directory = scratch or Path.cwd()
        directory.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(dir=directory, prefix=".copc-sort-") as folder:
            with (resources.stage("COPC sample read/sort") if resources else nullcontext({})) as stage:
                sorted_input = TimeSortedCopc(path, Path(folder), CHUNK, select=selection,
                                                      spatial=threads is not None)
                if threads is not None:
                    stage.update(decoder="lazrs COPC hierarchy nodes", decoder_threads=threads,
                                 decode_batch_max_points=CHUNK, decode_batch_max_compressed_bytes=COPC_BATCH_BYTES,
                                 decode_batch_max_nodes=COPC_BATCH_NODES)
                if sorted_input.points_read != declared_points:
                    raise ValueError("The COPC scan did not decode the point count declared by the header.")
                stage.update(bytes_read=path.stat().st_size, bytes_written=sorted_input.bytes_written,
                             points_read=sorted_input.points_read, sampled_echoes=sorted_input.points_spooled)
            return read_pulses(paths, max_multi, max_single, full, info=info, prepared=sorted_input,
                               scratch=scratch, resources=resources, _selection=selection,
                               scan_index=scan_index, threads=threads)

    lower, upper = np.min(origins, axis=0), np.max(maxima, axis=0)
    origin = np.floor(lower / 1000) * 1000
    crs = None
    first_chunk = True
    counts = np.zeros(2, np.int64)
    pulses = _Store("time", "first", "last", "returns", "source", "file", "angle")
    recovered = _Store("time", "first", "last", "returns", "source", "file", "angle")
    for index, path in enumerate(paths):
        try:
            with open_cloud(path) as reader:
                header = reader.header
                if not header.point_count:
                    continue
                file_time = _time_system(header, set(header.point_format.dimension_names))
                if file_time == "none":
                    info.has_gps_time = False
                    raise UnusableTime("A moving sensor's path needs GPS time per point; this file has none.")
                info.has_gps_time = True
                if info.time_system not in ("none", file_time):
                    raise ValueError("Input clouds must use the same GPS time system.")
                info.time_system = file_time
                file_crs = _crs(header)
                if file_crs is not None:
                    if any(axis.unit_conversion_factor != 1.0 for axis in file_crs.axis_info):
                        raise ValueError("Known coordinate reference systems must use metres.")
                    if crs is not None and not crs.equals(file_crs, ignore_axis_order=True):
                        raise ValueError("Input clouds must use the same coordinate reference system.")
                    crs = file_crs
                    reference = cloud_crs(crs, lower, upper)
                    epsg = reference.to_epsg() if reference is not None else None
                    info.crs = f"EPSG:{epsg}" if epsg else "unknown"

                if not np.isfinite(header.offsets).all() or not np.isfinite(header.scales).all() or np.any(header.scales <= 0):
                    raise ValueError("LAS coordinates must have finite offsets and positive finite scales.")
                if not np.isfinite(header.maxs - header.mins).all():
                    raise ValueError("LAS coordinate bounds must be finite.")

                shift = header.offsets - origin
                local = _FilePulses([path])
                time_min, time_max = np.inf, -np.inf
                if prepared is not None:
                    info.points_read += prepared.points_read
                    time_min, time_max = prepared.time_range
                    if time_min == time_max:
                        raise UnusableTime("GPS time must vary across the input points. "
                                           "A moving delivered path cannot be matched to a constant point clock.")
                parallel = None
                if prepared is None and (threads or 1) > 1 and header.point_count >= 2 * CHUNK:
                    from .parallel_read import SerialReadRequired, collect_ranges, probe_single_returns

                    strides = info.kept_every
                    if first_chunk:
                        looks_single = probe_single_returns(path, min(PROBE, CHUNK))
                        strides = _strides(info.points, full, max_multi, MAX_SINGLE_ONLY if looks_single else max_single)
                    try:
                        parallel = collect_ranges(path, index, header.point_count, CHUNK, header.point_format.size,
                                                  shift, strides, scan_index, threads, min(PROBE, CHUNK))
                    except (SerialReadRequired, UnusableTime):
                        pass
                if parallel is not None:
                    local, time_min, time_max, points, slots = parallel
                    info.points_read += points
                    info.read_workers = max(info.read_workers, slots)
                    info.kept_every = strides
                    first_chunk = False
                else:
                    chunks = prepared.chunks() if prepared is not None else pulse_chunks(point_chunks(path, CHUNK), CHUNK)
                    for chunk in chunks:
                        if _selection is None:
                            info.points_read += len(chunk)
                        if scan_index is not None:
                            scan_index.observe(index, chunk)
                        if first_chunk:
                            looks_single = (_looks_single_return(chunk) if _selection is None
                                            else _selection.looks_single(prepared))
                            single_limit = MAX_SINGLE_ONLY if looks_single else max_single
                            info.kept_every = _strides(info.points, full, max_multi, single_limit)
                            info.population_counts_estimated = _selection is not None and info.kept_every != (1, 1)
                            first_chunk = False
                        low, high = float(chunk.gps_time.min()), float(chunk.gps_time.max())
                        if not np.isfinite([low, high]).all():
                            raise UnusableTime("GPS times must be finite.")
                        time_min, time_max = min(time_min, low), max(time_max, high)
                        local.collect(chunk, index, shift, info.kept_every, prepared is not None, _selection is not None)
                if time_min == time_max:
                    raise UnusableTime("GPS time must vary across the input points. "
                                       "A moving delivered path cannot be matched to a constant point clock.")

                mode = "paired"
                single_stride = 1 if full else max(1, -(-info.points // max_single))
                strides = (info.kept_every[0], single_stride)
                okay = local.finish(info.kept_every[0])
                if not okay:
                    exact_info = local.info
                    del local
                    if prepared is None:
                        mode = "time_tolerance"
                        local = _FilePulses([path])
                        original_index = scan_index.begin_repair(index) if scan_index is not None else None
                        for chunk in _fallback_chunks(path, mode, False):
                            local.collect(chunk, index, shift, strides, False)
                            if scan_index is not None:
                                scan_index.observe(index, chunk)
                        okay = local.finish(strides[0])
                        if not okay and scan_index is not None:
                            scan_index.abandon_repair(index, original_index)
                    if not okay:
                        mode = "single"
                        info.single_return_files.append(index)
                        local = _FilePulses([path])
                        for chunk in _fallback_chunks(path, mode, prepared is not None):
                            local.collect(chunk, index, shift, strides, True)
                        local.finish(strides[0])
                        local.info = exact_info
                if local.counts[0] == 0 and mode == "paired":
                    mode = "single"
                angles = local.has_angles()
                info.return_modes.append({"file": str(path), "index": index, "mode": mode, "has_scan_angle": angles,
                                          "unpaired_multi_echoes": local.info.unpaired_multi_echoes})
                for key in _PARTIAL_COUNTS:
                    setattr(info, key, getattr(info, key) + getattr(local.info, key))
                counts += local.counts
                pulses.extend(local.pulses)
                recovered.extend(local.recovered)
        except UnusableTime as error:
            raise UnusableTime(f"{path}: {error}") from error
        except (OSError, laspy.errors.LaspyException, LazrsError, StructError, EOFError, OverflowError, ValueError) as error:
            raise ValueError(f"{path}: {error}") from error
    pulses.extend(recovered)
    angle = pulses.array("angle", np.float32)
    first = pulses.array("first", np.float32, 3)
    returns = pulses.array("returns", np.uint8)
    last = pulses.array("last", np.float32, 3) if np.any(returns > 1) else first
    result = Pulses(origin, pulses.array("time", np.float64), first, last, returns,
                    pulses.array("source", np.uint16), pulses.array("file", np.uint16),
                    angle, True, bool(np.any(np.isfinite(angle) & (np.abs(angle) > 0) & (np.abs(angle) < 90))), float(counts[0] / max(counts.sum(), 1)))
    return result, info
