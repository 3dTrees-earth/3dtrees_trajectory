"""Read a point cloud into pulses: the first and last return of every laser shot, in file order.

A pulse is a run of consecutive records with the same GPS time and point source whose return numbers count up by one.
The time alone is not enough: multi-beam scanners give a dozen pulses one time stamp, and some writers fill
NumberOfReturns wrongly, so the file order is what ties a pulse together.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import laspy
import numpy as np

CHUNK = 5_000_000


@dataclass
class Pulses:
    """One row per pulse. ``first`` and ``last`` are equal for single-return pulses."""

    time: np.ndarray  # float64, the cloud's own clock
    first: np.ndarray  # (n, 3) float64
    last: np.ndarray  # (n, 3) float64
    returns: np.ndarray  # int16, records in the pulse
    source: np.ndarray  # int32, PointSourceId
    file: np.ndarray  # int16, index of the input file
    scan_angle: np.ndarray  # float32, degrees off nadir as written in the file (0 when absent)
    has_time: bool
    has_scan_angle: bool

    def __len__(self) -> int:
        return len(self.time)

    def subset(self, mask: np.ndarray) -> "Pulses":
        return Pulses(self.time[mask], self.first[mask], self.last[mask], self.returns[mask], self.source[mask],
                      self.file[mask], self.scan_angle[mask], self.has_time, self.has_scan_angle)

    @property
    def multi(self) -> np.ndarray:
        return self.returns > 1


@dataclass
class CloudInfo:
    crs: str  # "EPSG:xxxx" or "unknown"
    time_system: str  # gps_adjusted | gps_sow | none
    points: int
    files: list[str]


def _crs(header: laspy.LasHeader) -> str:
    try:
        crs = header.parse_crs()
    except Exception:  # noqa: BLE001 - broken VLRs are common; the CRS is then unknown, not fatal
        return "unknown"
    if crs is None:
        return "unknown"
    epsg = crs.to_epsg()
    return f"EPSG:{epsg}" if epsg else "unknown"


MAX_PULSES = 100_000_000  # about 7 GB of pulse arrays; larger clouds keep every n-th pulse


def read_pulses(paths: list[Path], keep_every: int = 0) -> tuple[Pulses, CloudInfo]:
    """Pulses of all files, in file order. ``keep_every`` > 1 keeps every n-th pulse; 0 picks the smallest n that
    keeps at most MAX_PULSES (counting every record as a pulse, the worst case)."""
    if keep_every <= 0:
        total = 0
        for path in paths:
            with laspy.open(path) as reader:
                total += reader.header.point_count
        keep_every = max(1, -(-total // MAX_PULSES))
    parts: dict[str, list[np.ndarray]] = {k: [] for k in ("time", "first", "last", "returns", "source", "file", "angle")}
    crs, time_system, points, has_time, has_angle = None, None, 0, True, True
    for index, path in enumerate(paths):
        with laspy.open(path) as reader:
            header = reader.header
            names = set(header.point_format.dimension_names)
            has_time &= "gps_time" in names
            has_angle &= bool({"scan_angle", "scan_angle_rank"} & names)
            file_crs = _crs(header)
            crs = file_crs if crs in (None, "unknown") else crs
            if "gps_time" not in names:
                file_time = "none"
            elif header.global_encoding.gps_time_type == laspy.header.GpsTimeType.STANDARD:
                file_time = "gps_adjusted"
            else:
                file_time = "gps_sow"
            time_system = file_time if time_system in (None, file_time) else "mixed"
            points += header.point_count
            for chunk in reader.chunk_iterator(CHUNK):
                n = len(chunk)
                t = np.asarray(chunk.gps_time) if "gps_time" in names else np.zeros(n)
                rn = np.asarray(chunk.return_number).astype(np.int16)
                src = np.asarray(chunk.point_source_id).astype(np.int32)
                joined = (t[1:] == t[:-1]) & (src[1:] == src[:-1]) & (rn[1:] == rn[:-1] + 1)
                if "gps_time" not in names:
                    joined[:] = False  # without time, returns can't be tied together
                start = np.flatnonzero(np.r_[True, ~joined])
                end = np.r_[start[1:], n] - 1
                if keep_every > 1:
                    start, end = start[::keep_every], end[::keep_every]
                xyz = np.c_[chunk.x, chunk.y, chunk.z]
                if "scan_angle" in names:
                    angle = np.asarray(chunk["scan_angle"]).astype(np.float32) * 0.006
                elif "scan_angle_rank" in names:
                    angle = np.asarray(chunk["scan_angle_rank"]).astype(np.float32)
                else:
                    angle = np.zeros(n, np.float32)
                parts["time"].append(t[start])
                parts["first"].append(xyz[start])
                parts["last"].append(xyz[end])
                parts["returns"].append((end - start + 1).astype(np.int16))
                parts["source"].append(src[start])
                parts["file"].append(np.full(len(start), index, np.int16))
                parts["angle"].append(angle[start])
    cat = {k: np.concatenate(v) if v else np.empty(0) for k, v in parts.items()}
    pulses = Pulses(cat["time"], cat["first"].reshape(-1, 3), cat["last"].reshape(-1, 3), cat["returns"], cat["source"],
                    cat["file"], cat["angle"], has_time, has_angle and bool(np.any(cat["angle"] != 0)))
    return pulses, CloudInfo(crs or "unknown", time_system or "none", points, [str(p) for p in paths])
