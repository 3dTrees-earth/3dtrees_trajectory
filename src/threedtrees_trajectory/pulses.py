"""Read a point cloud into pulses: the first and last return of every laser shot.

A pulse is a run of consecutive records with the same GPS time and point source whose return numbers count up by one.
The time alone is not enough: multi-beam scanners give a dozen pulses one time stamp, and some writers fill
NumberOfReturns wrongly, so the file order is what ties a pulse together. Some software stores the points out of
acquisition order (2587, LiDAR360); records of a multi-return pulse that end up alone are then collected and joined
by GPS time and point source after reading.

Memory stays bounded: coordinates are float32 relative to ``origin`` (sub-millimetre within 10 km), and at most
MAX_MULTI multi-return and MAX_SINGLE single-return pulses are kept, as every n-th pulse so the whole acquisition
time stays covered; a cloud whose first million points hold no multi-return pulses keeps MAX_SINGLE_ONLY single
returns instead. ``full`` keeps every pulse.
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

CHUNK = 5_000_000
PROBE = 1_000_000  # points read first to tell a single-return cloud
MAX_MULTI = 40_000_000
MAX_SINGLE = 10_000_000  # single returns next to multi-return pulses
MAX_SINGLE_ONLY = 40_000_000  # clouds without multi-return pulses, where single returns are all there is
MAX_STRAY = 40_000_000  # separate limit for records that need a global sort to reconstruct pulses
LAZ_FIELDS = (laspy.DecompressionSelection.base().decompress_z().decompress_gps_time()
              .decompress_scan_angle().decompress_point_source_id())


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

    @classmethod
    def from_paths(cls, paths: list[Path]) -> CloudInfo:
        return cls("unknown", "none", 0, [str(path) for path in paths])


def _crs(header: laspy.LasHeader) -> CRS | None:
    try:
        return header.parse_crs()
    except (CRSError, ValueError):
        return None  # a malformed CRS record leaves the CRS unknown


def _open(path: Path) -> laspy.LasReader:
    reader = laspy.open(path, decompression_selection=LAZ_FIELDS)
    if any(type(vlr).__name__ == "CopcInfoVlr" for vlr in reader.header.vlrs):
        # COPC uses variable-sized chunks, which break the parallel decoder.
        reader.laz_backend = laspy.LazBackend.Lazrs
    return reader


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
    """Pulses from lone records: same file, source and GPS time, first and last by return number."""
    t, xyz, number = stray.array("time", np.float64), stray.array("xyz", np.float32, 3), stray.array("number", np.uint8)
    source, file, angle = stray.array("source", np.uint16), stray.array("file", np.uint16), stray.array("angle", np.float32)
    order = np.lexsort((number, t, source, file))
    t, source, file = t[order], source[order], file[order]
    same = (t[1:] == t[:-1]) & (source[1:] == source[:-1]) & (file[1:] == file[:-1])
    start = np.flatnonzero(np.r_[True, ~same])
    end = np.r_[start[1:], len(t)] - 1
    pulses.add(time=t[start], first=xyz[order[start]], last=xyz[order[end]], returns=np.minimum(end - start + 1, 255).astype(np.uint8),
               source=source[start], file=file[start], angle=angle[order[start]])


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
                   strides: tuple[int, int], counts: np.ndarray, pulses: _Store, stray: _Store) -> int:
    multi_every, single_every = strides
    rn = np.asarray(chunk.return_number)
    source = chunk.point_source_id
    joined = (time[1:] == time[:-1]) & (source[1:] == source[:-1]) & (rn[1:] == rn[:-1] + 1)
    start = np.flatnonzero(np.r_[True, ~joined])
    end = np.r_[start[1:], len(chunk)] - 1
    single = start == end
    # A lone record of a multi-return pulse has its partners elsewhere in the file.
    lone = single & (np.asarray(chunk.number_of_returns)[start] > 1)
    keep = np.zeros(len(start), bool)
    multi_idx, single_idx = np.flatnonzero(~single), np.flatnonzero(single & ~lone)
    counts += [len(multi_idx) + int(lone.sum()) // 2, len(single_idx)]
    keep[multi_idx[::multi_every]] = True
    keep[single_idx[::single_every]] = True
    lone_idx = start[lone]
    if multi_every > 1:
        # Thin by time so a pulse's records are kept or dropped together.
        lone_idx = lone_idx[time[lone_idx].view(np.int64) % multi_every == 0]
    start, end = start[keep], end[keep]
    first = _coordinates(chunk, start, shift)
    multi = start != end
    last = first
    if np.any(multi):
        last = first.copy()
        last[multi] = _coordinates(chunk, end[multi], shift)
    pulses.add(time=time[start], first=first, last=last, returns=(end - start + 1).astype(np.uint8),
               source=source[start], file=np.full(len(start), file, dtype=np.uint16), angle=_angles(chunk, start))
    if len(lone_idx):
        stray.add(time=time[lone_idx], xyz=_coordinates(chunk, lone_idx, shift), number=rn[lone_idx],
                  source=source[lone_idx], file=np.full(len(lone_idx), file, dtype=np.uint16),
                  angle=_angles(chunk, lone_idx))
    return len(lone_idx)


def read_pulses(paths: list[Path], max_multi: int = MAX_MULTI, max_single: int = MAX_SINGLE,
                full: bool = False, *, info: CloudInfo | None = None) -> tuple[Pulses, CloudInfo]:
    """Read LAS/LAZ pulse endpoints and retain every n-th pulse past the sampling budgets."""
    info = CloudInfo.from_paths(paths) if info is None else info
    if len(paths) > np.iinfo(np.uint16).max + 1:
        raise ValueError("At most 65,536 input files can be combined.")
    for path in paths:
        try:
            if not path.is_file():
                raise ValueError("Input must be an existing regular file.")
            with laspy.open(path, read_evlrs=False) as reader:
                info.points += reader.header.point_count
        except (OSError, laspy.errors.LaspyException, StructError, ValueError, OverflowError, EOFError) as error:
            raise ValueError(f"{path}: unable to read the LAS/LAZ header ({error}).") from error
    if not info.points:
        raise ValueError("The cloud has no points.")

    origin = crs = None
    counts = np.zeros(2, np.int64)
    pulses = _Store("time", "first", "last", "returns", "source", "file", "angle")
    stray = _Store("time", "xyz", "number", "source", "file", "angle")
    stray_count = 0
    time_min, time_max = np.inf, -np.inf
    for index, path in enumerate(paths):
        try:
            with _open(path) as reader:
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

                first_file = origin is None
                if first_file:
                    origin = header.offsets.copy()
                shift = header.offsets - origin
                while reader.points_read < header.point_count:
                    expected = min(CHUNK, header.point_count - reader.points_read)
                    chunk = reader.read_points(expected)
                    info.points_read += len(chunk)
                    if len(chunk) != expected:
                        raise ValueError(f"Truncated point data: expected {expected:,} records, read {len(chunk):,}.")
                    if first_file:
                        single_limit = MAX_SINGLE_ONLY if _looks_single_return(chunk) else max_single
                        multi_every = 1 if full else max(1, -(-info.points // max_multi))
                        single_every = 1 if full else max(1, -(-info.points // single_limit))
                        info.kept_every = (multi_every, single_every)
                        first_file = False
                    time = chunk.gps_time
                    low, high = float(time.min()), float(time.max())
                    if not np.isfinite([low, high]).all():
                        raise ValueError("GPS times must be finite.")
                    time_min, time_max = min(time_min, low), max(time_max, high)
                    stray_count += _collect_chunk(chunk, time, index, shift, info.kept_every, counts, pulses, stray)
                    if stray_count > MAX_STRAY:
                        raise ValueError("Too many separated multi-return records to join. Provide the original "
                                         "LAS/LAZ strips in acquisition order.")
        except (OSError, laspy.errors.LaspyException, LazrsError, StructError, EOFError, OverflowError, ValueError) as error:
            raise ValueError(f"{path}: {error}") from error
    if time_min == time_max:
        raise ValueError("GPS time must vary across the input points.")

    if stray_count:
        _join_stray(stray, pulses)
    angle = pulses.array("angle", np.float32)
    first = pulses.array("first", np.float32, 3)
    returns = pulses.array("returns", np.uint8)
    last = pulses.array("last", np.float32, 3) if np.any(returns > 1) else first
    result = Pulses(origin, pulses.array("time", np.float64), first, last, returns,
                    pulses.array("source", np.uint16), pulses.array("file", np.uint16),
                    angle, True, bool(np.any(angle != 0)), float(counts[0] / max(counts.sum(), 1)))
    return result, info
