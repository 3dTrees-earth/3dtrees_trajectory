"""Read a point cloud into pulses: the first and last return of every laser shot.

A pulse is a run of consecutive records with the same GPS time and point source whose return numbers count up by one.
The time alone is not enough: multi-beam scanners give a dozen pulses one time stamp, and some writers fill
NumberOfReturns wrongly, so the file order is what ties a pulse together. Some software stores the points out of
acquisition order (2587, LiDAR360); records of a multi-return pulse that end up alone are then collected and joined
by GPS time and point source after reading.

Memory stays bounded: coordinates are float32 relative to ``origin`` (sub-millimetre within 10 km), and at most
MAX_MULTI multi-return and MAX_SINGLE single-return pulses are kept, as every n-th pulse so the whole acquisition
time stays covered. A cloud whose first million points hold no multi-return pulses is read faster: evenly spread
blocks of BLOCK points, MAX_SINGLE_ONLY in total, so most of it is never decoded. ``full`` reads every point.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import laspy
import numpy as np

CHUNK = 5_000_000
BLOCK = 200_000  # points per contiguous block when a single-return cloud is sampled; small enough to spread
MAX_MULTI = 40_000_000
MAX_SINGLE = 10_000_000  # single returns next to multi-return pulses
MAX_SINGLE_ONLY = 40_000_000  # clouds without multi-return pulses, where single returns are all there is


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

    def __len__(self) -> int:
        return len(self.time)

    def subset(self, index: np.ndarray) -> "Pulses":
        return Pulses(self.origin, self.time[index], self.first[index], self.last[index], self.returns[index],
                      self.source[index], self.file[index], self.scan_angle[index], self.has_time, self.has_scan_angle)

    @property
    def multi(self) -> np.ndarray:
        return self.returns > 1


@dataclass
class CloudInfo:
    crs: str  # "EPSG:xxxx" or "unknown"
    time_system: str  # gps_adjusted | gps_sow | none
    points: int
    files: list[str]
    kept_every: tuple[int, int, int]  # (multi-return pulse, single-return pulse, block) strides


def _crs(header: laspy.LasHeader) -> str:
    try:
        crs = header.parse_crs()
    except Exception:  # noqa: BLE001 - broken VLRs are common; the CRS is then unknown, not fatal
        return "unknown"
    if crs is None:
        return "unknown"
    epsg = crs.to_epsg()
    return f"EPSG:{epsg}" if epsg else "unknown"


def _backend() -> laspy.LazBackend | None:
    """Parallel LAZ decoding when available (threads: RAYON_NUM_THREADS, set by the CLI)."""
    available = laspy.LazBackend.detect_available()
    return laspy.LazBackend.LazrsParallel if laspy.LazBackend.LazrsParallel in available else None


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
            self.columns[name].append(value)

    def array(self, name: str, dtype, width: int = 0) -> np.ndarray:
        parts = self.columns[name]
        if parts:
            return np.concatenate(parts).astype(dtype, copy=False)
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
    pulses.add(time=t[start], first=xyz[order[start]], last=xyz[order[end]], returns=(end - start + 1),
               source=source[start], file=file[start], angle=angle[order[start]])


def _multi_share(path: Path) -> float:
    """Share of records in multi-return pulses among the first BLOCK points (tied by time, source, return number)."""
    with laspy.open(path, laz_backend=_backend()) as reader:
        names = set(reader.header.point_format.dimension_names)
        if "gps_time" not in names or reader.header.point_count == 0:
            return 0.0
        chunk = reader.read_points(min(BLOCK, reader.header.point_count))
    t, rn, src = np.asarray(chunk.gps_time), np.asarray(chunk.return_number), np.asarray(chunk.point_source_id)
    joined = (t[1:] == t[:-1]) & (src[1:] == src[:-1]) & (rn[1:] == rn[:-1] + 1)
    return float(joined.mean()) if len(joined) else 0.0


def _chunks(reader: laspy.LasReader, every_block: int):
    """All points in CHUNK pieces, or every ``every_block``-th block of BLOCK points."""
    if every_block <= 1:
        yield from reader.chunk_iterator(CHUNK)
        return
    for start in range(0, reader.header.point_count, BLOCK * every_block):
        reader.seek(start)
        yield reader.read_points(BLOCK)


def read_pulses(paths: list[Path], max_multi: int = MAX_MULTI, max_single: int = MAX_SINGLE,
                full: bool = False) -> tuple[Pulses, CloudInfo]:
    """Pulses of all files; every n-th pulse (or block, for single-return clouds) beyond the budgets."""
    total = 0
    origin = None
    for path in paths:
        with laspy.open(path) as reader:
            total += reader.header.point_count
            origin = np.array(reader.header.offsets, float) if origin is None else origin
    single_only = bool(paths) and _multi_share(paths[0]) < 0.01
    every_block = 1 if full or not single_only else max(1, -(-total // MAX_SINGLE_ONLY))
    multi_every = 1 if full else max(1, -(-total // max_multi))  # counts every record as a pulse: the worst case
    single_every = 1 if full or single_only else max(1, -(-total // max_single))
    pulses = _Store("time", "first", "last", "returns", "source", "file", "angle")
    stray = _Store("time", "xyz", "number", "source", "file", "angle")
    crs, time_system, has_time, has_angle = None, None, True, True
    for index, path in enumerate(paths):
        with laspy.open(path, laz_backend=_backend()) as reader:
            header = reader.header
            names = set(header.point_format.dimension_names)
            has_time &= "gps_time" in names
            has_angle &= bool({"scan_angle", "scan_angle_rank"} & names)
            file_crs = _crs(header)
            crs = file_crs if crs in (None, "unknown") else crs
            file_time = _time_system(header, names)
            time_system = file_time if time_system in (None, file_time) else "mixed"
            scale = np.asarray(header.scales)
            shift = np.asarray(header.offsets) - origin
            for chunk in _chunks(reader, every_block):
                n = len(chunk)
                t = np.asarray(chunk.gps_time) if "gps_time" in names else np.zeros(n)
                rn = np.asarray(chunk.return_number)
                src = np.asarray(chunk.point_source_id)
                if "gps_time" in names:
                    joined = (t[1:] == t[:-1]) & (src[1:] == src[:-1]) & (rn[1:] == rn[:-1] + 1)
                else:
                    joined = np.zeros(n - 1, bool)  # without time, returns can't be tied together
                start = np.flatnonzero(np.r_[True, ~joined])
                end = np.r_[start[1:], n] - 1
                single = start == end
                # a lone record of a multi-return pulse: its partners are elsewhere in the file
                lone = single & (np.asarray(chunk.number_of_returns)[start] > 1) & ("gps_time" in names)
                keep = np.zeros(len(start), bool)
                multi_idx, single_idx = np.flatnonzero(~single), np.flatnonzero(single & ~lone)
                keep[multi_idx[::multi_every]] = True
                keep[single_idx[::single_every]] = True
                lone_idx = start[lone]
                if multi_every > 1:  # thin by time, so a pulse's records are kept or dropped together
                    lone_idx = lone_idx[t[lone_idx].view(np.int64) % multi_every == 0]
                start, end = start[keep], end[keep]
                used = np.unique(np.r_[start, end, lone_idx])
                xyz = np.zeros((n, 3), np.float32)
                for axis, raw in enumerate((chunk.X, chunk.Y, chunk.Z)):
                    xyz[used, axis] = np.asarray(raw)[used] * scale[axis] + shift[axis]
                if "scan_angle" in names:
                    angle = np.asarray(chunk["scan_angle"]).astype(np.float32) * 0.006
                elif "scan_angle_rank" in names:
                    angle = np.asarray(chunk["scan_angle_rank"]).astype(np.float32)
                else:
                    angle = np.zeros(n, np.float32)
                pulses.add(time=t[start], first=xyz[start], last=xyz[end], returns=end - start + 1, source=src[start],
                           file=np.full(len(start), index), angle=angle[start])
                stray.add(time=t[lone_idx], xyz=xyz[lone_idx], number=rn[lone_idx], source=src[lone_idx],
                          file=np.full(len(lone_idx), index), angle=angle[lone_idx])
    if stray.columns["time"] and sum(len(a) for a in stray.columns["time"]):
        _join_stray(stray, pulses)
    angle = pulses.array("angle", np.float32)
    result = Pulses(origin if origin is not None else np.zeros(3), pulses.array("time", np.float64),
                    pulses.array("first", np.float32, 3), pulses.array("last", np.float32, 3),
                    pulses.array("returns", np.uint8), pulses.array("source", np.uint16), pulses.array("file", np.uint16),
                    angle, has_time, has_angle and bool(np.any(angle != 0)))
    info = CloudInfo(crs or "unknown", time_system or "none", total, [str(p) for p in paths],
                     (multi_every, single_every, every_block))
    return result, info
