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

from dataclasses import dataclass
from pathlib import Path
from struct import error as StructError

import laspy
import numpy as np
from lazrs import LazrsError
from pyproj import CRS
from pyproj.exceptions import CRSError

from .cloud import CHUNK, is_copc, open_cloud, point_chunks, pulse_chunks, sorted_copc

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
    incomplete_multi_pulses: int = 0
    multi_pulses: int = 0
    joined_out_of_order_points: int = 0

    @classmethod
    def from_paths(cls, paths: list[Path]) -> CloudInfo:
        return cls("unknown", "none", 0, [str(path) for path in paths])


def _crs(header: laspy.LasHeader) -> CRS | None:
    try:
        return header.parse_crs()
    except (CRSError, ValueError):
        return None  # a malformed CRS record leaves the CRS unknown


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

    def array(self, name: str, dtype, width: int = 0) -> np.ndarray:
        parts = self.columns.pop(name)
        if len(parts) == 1:
            return parts[0].astype(dtype, copy=False)
        if parts:
            return np.concatenate(parts, dtype=dtype, casting="unsafe")
        return np.empty((0, width) if width else 0, dtype)


def _join_stray(stray: _Store, pulses: _Store) -> None:
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


def _sample(time: np.ndarray, stride: int) -> np.ndarray:
    if stride == 1:
        return np.ones(len(time), bool)
    # A deterministic hash of the timestamp, independent of COPC order, file
    # boundaries, and whether a return had to be joined out of order.
    key = time.view(np.uint64).copy()
    key ^= key >> 33
    key *= np.uint64(0xff51afd7ed558ccd)
    key ^= key >> 33
    return key % stride == 0


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
                   copc: bool) -> int:
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
    counts += [int(complete.sum()), int(single.sum())]
    info.incomplete_multi_pulses += int(incomplete.sum())
    info.multi_pulses += int(complete.sum())
    keep = (complete & _sample(time[start], multi_every)) | (single & _sample(time[start], single_every))
    lone_idx = None
    # Retain every record of any incomplete LAZ pulse so it can be repaired.
    if not copc and incomplete.any():
        sizes = end - start + 1
        lone_idx = np.flatnonzero(np.repeat(incomplete & _sample(time[start], multi_every), sizes))
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


def read_pulses(paths: list[Path], max_multi: int = MAX_MULTI, max_single: int = MAX_SINGLE,
                full: bool = False, *, info: CloudInfo | None = None, prepared=None) -> tuple[Pulses, CloudInfo]:
    """Build pulse endpoints from LAS/LAZ or a time-sorted COPC within the sampling budgets."""
    info = CloudInfo.from_paths(paths) if info is None else info
    if len(paths) > np.iinfo(np.uint16).max + 1:
        raise ValueError("At most 65,536 input files can be combined.")
    origins = []
    copc_paths = []
    for path in paths:
        try:
            if not path.is_file():
                raise ValueError("Input must be an existing regular file.")
            with laspy.open(path, read_evlrs=False) as reader:
                if reader.header.point_count:
                    origins.append(reader.header.mins)
                info.points += reader.header.point_count
                if is_copc(reader.header):
                    copc_paths.append(path)
        except (OSError, laspy.errors.LaspyException, StructError, ValueError, OverflowError, EOFError) as error:
            raise ValueError(f"{path}: unable to read the LAS/LAZ header ({error}).") from error
    if not info.points:
        raise ValueError("The cloud has no points.")

    if copc_paths and (len(paths) != 1 or len(copc_paths) != 1):
        raise ValueError("Use one merged COPC or a list of per-file LAS/LAZ inputs.")
    if copc_paths and prepared is None:
        # Library callers have the same semantics as the CLI.
        info.points = 0
        with sorted_copc(paths, Path.cwd(), CHUNK) as sorted_input:
            return read_pulses(paths, max_multi, max_single, full, info=info, prepared=sorted_input)

    origin = np.floor(np.min(origins, axis=0) / 1000) * 1000
    crs = None
    first_chunk = True
    counts = np.zeros(2, np.int64)
    pulses = _Store("time", "first", "last", "returns", "source", "file", "angle")
    stray = _Store("time", "xyz", "number", "source", "file", "angle", "total", "channel")
    stray_count = 0
    time_min, time_max = np.inf, -np.inf
    for index, path in enumerate(paths):
        try:
            with open_cloud(path) as reader:
                header = reader.header
                if not header.point_count:
                    continue
                file_time = _time_system(header, set(header.point_format.dimension_names))
                if file_time == "none":
                    raise ValueError("A moving sensor's path needs GPS time per point; this file has none.")
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
                    epsg = crs.to_epsg()
                    info.crs = f"EPSG:{epsg}" if epsg else "unknown"

                if not np.isfinite(header.offsets).all() or not np.isfinite(header.scales).all() or np.any(header.scales <= 0):
                    raise ValueError("LAS coordinates must have finite offsets and positive finite scales.")
                if not np.isfinite(header.maxs - header.mins).all():
                    raise ValueError("LAS coordinate bounds must be finite.")

                shift = header.offsets - origin
                chunks = prepared.chunks() if prepared is not None else pulse_chunks(point_chunks(path, CHUNK), CHUNK)
                for chunk in chunks:
                    info.points_read += len(chunk)
                    if first_chunk:
                        single_limit = MAX_SINGLE_ONLY if _looks_single_return(chunk) else max_single
                        multi_every = 1 if full else max(1, -(-info.points // max_multi))
                        single_every = 1 if full else max(1, -(-info.points // single_limit))
                        info.kept_every = (multi_every, single_every)
                        first_chunk = False
                    time = chunk.gps_time
                    low, high = float(time.min()), float(time.max())
                    if not np.isfinite([low, high]).all():
                        raise ValueError("GPS times must be finite.")
                    time_min, time_max = min(time_min, low), max(time_max, high)
                    stray_count += _collect_chunk(chunk, time, index, shift, info.kept_every, counts, pulses, stray,
                                                  info, prepared is not None)
                    if stray_count > MAX_STRAY:
                        raise ValueError("Too many separated multi-return records to join. Provide the original "
                                         "LAS/LAZ strips in acquisition order.")
        except (OSError, laspy.errors.LaspyException, LazrsError, StructError, EOFError, OverflowError, ValueError) as error:
            raise ValueError(f"{path}: {error}") from error
    if time_min == time_max:
        raise ValueError("GPS time must vary across the input points.")

    if stray_count:
        info.joined_out_of_order_points = stray_count
        before = sum(map(len, pulses.columns["time"]))
        _join_stray(stray, pulses)
        recovered = sum(map(len, pulses.columns["time"])) - before
        counts[0] += recovered * info.kept_every[0]
    if info.incomplete_multi_pulses > max(100, 2 * int(counts[0])):
        raise ValueError("Too many declared multi-return records cannot be paired. Check return numbering and GPS time "
                         "or provide the original LAS/LAZ strips in acquisition order.")
    angle = pulses.array("angle", np.float32)
    first = pulses.array("first", np.float32, 3)
    returns = pulses.array("returns", np.uint8)
    last = pulses.array("last", np.float32, 3) if np.any(returns > 1) else first
    result = Pulses(origin, pulses.array("time", np.float64), first, last, returns,
                    pulses.array("source", np.uint16), pulses.array("file", np.uint16),
                    angle, True, bool(np.any(angle != 0)), float(counts[0] / max(counts.sum(), 1)))
    return result, info
