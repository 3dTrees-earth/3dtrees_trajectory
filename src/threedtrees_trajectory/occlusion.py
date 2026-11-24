"""Occlusion of a moving-sensor point cloud: which voxels the lidar saw as filled or empty, and which stayed hidden.

Every pulse is traced from the sensor (the trajectory interpolated at its GPS time) through its returns and on to the
grid edge with OccPy's voxel traversal (Kuekenbrink et al. 2017, the compiled ``raytr`` core of ``occpy_ls``).
Each voxel collects hits, passes (miss) and passes behind the last return (occluded); OccPy's rule turns the counts
into a class: filled (any hit), empty (passed, no hit), occluded (only behind returns), unobserved (no pulse).

Heights above ground use the supplied terrain raster where it has data. The fallback is the lowest filled voxel
of each column, with columns more than 2 m above their 5 m neighbourhood median replaced by that median.

Outputs: a COPC with one point per voxel between the ground and the canopy top that is not empty (class, counts,
occluded fraction, height above ground), and a COG with the voxel columns as pixels and one band per voxel layer
above ground holding the occluded share of the observed voxels, plus a last band with the column's unobserved share.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import shutil
import signal
import subprocess
from contextlib import contextmanager, nullcontext
from dataclasses import dataclass
from pathlib import Path
from struct import pack
from tempfile import TemporaryDirectory

import laspy
import numpy as np
from laspy.vlrs.vlrlist import VLRList

from .cloud import is_copc, open_cloud
from .output import delivered_miss_flag
from .pulses import Pulses, _crs, _time_system, cloud_crs, is_local_geocentric
from .resources import OutOfMemory, RunResources, available_memory, decoding_threads, raise_if_oom
from . import ray_walk, time_windows
from .window_input import ScanIndex
from .window_jobs import execute, identify, lock_job

FILLED, EMPTY, OCCLUDED, UNOBSERVED = 1, 2, 3, 4  # OccPy's classes
# Native grids, list conversion, NumPy counts and output windows all overlap
# with the active chunk. Bound the estimated working memory before allocating.
GIB = 2**30
DEFAULT_BUDGET = 16 * GIB
HARD_BUDGET = 64 * GIB
WINDOW = 128  # XY block size for deriving and writing the products
SPIKE_M = 2.0  # a lowest-hit ground this far above its neighbourhood is vegetation, not ground
SPIKE_WINDOW_M = 5.0
TIME_MARGIN = 1.0  # seconds a pulse may lie outside the trajectory's span
UNTWINE = os.environ.get("TRAJECTORY_UNTWINE", "/opt/untwine/bin/untwine")
LAS14_EVLR_COUNT_OFFSET = 243


@dataclass
class Grid:
    origin: np.ndarray  # min x, y, z of the grid
    shape: tuple[int, int, int]
    voxel: float
    points: int = 0


class PathCoverageError(ValueError):
    """The path and cloud clocks do not overlap."""


class TerrainError(ValueError):
    """The supplied terrain raster cannot be used."""


def read_trajectory(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Time and xyz of a path in the 3Dtrees trajectory CSV (``#`` header lines, then time,x,y,z first)."""
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = csv.reader(line for line in stream if line.strip() and not line.lstrip().startswith("#"))
        if next(rows, [])[:4] != ["time", "x", "y", "z"]:
            raise ValueError(f"{path.name}: expected a 3Dtrees trajectory CSV with time,x,y,z columns")
        values = []
        for row in rows:
            if len(row) < 4 or not row[0].strip():
                raise ValueError(f"{path.name}: a path with a time for every position is needed for occlusion")
            values.append(row[:4])
    if not values:
        raise ValueError(f"{path.name}: a path with a time for every position is needed for occlusion")
    data = np.array(values, np.float64)
    order = np.argsort(data[:, 0], kind="stable")
    return data[order, 0], data[order, 1:4]


GPS_EPOCH_UNIX = 315_964_800  # 1980-01-06 in Unix seconds
WEEK_S = 604_800
# GPS - UTC in seconds, from the Unix time it applies (1999 onwards; no leap second since 2017).
LEAP_SECONDS = ((915_148_800, 13), (1_136_073_600, 14), (1_230_768_000, 15), (1_341_100_800, 16),
                (1_435_708_800, 17), (1_483_228_800, 18))
PATH_CLOCK_MIN_SHARE = 0.5  # share of the cloud's time span the path must cover before its clock is trusted


def read_trajectory_meta(path: Path) -> dict[str, str]:
    """The ``# key: value`` header lines of a 3Dtrees trajectory CSV."""
    meta = {}
    with path.open(encoding="utf-8-sig") as stream:
        for line in stream:
            if line.strip() and not line.lstrip().startswith("#"):
                break
            key, _, value = line.lstrip("# ").partition(":")
            meta[key.strip()] = value.strip()
    return meta


def read_trajectory_crs(path: Path) -> tuple[str | None, str]:
    """The path's ``# crs`` (``EPSG:n`` or ``local``; legacy ``unknown`` reads as local; None without a crs line)
    and its ``# vertical_crs``."""
    meta = read_trajectory_meta(path)
    crs = meta.get("crs")
    if crs is not None and crs.lower() in ("", "unknown", "local"):
        crs = "local"
    return crs, meta.get("vertical_crs") or "unknown"


def _is_epsg(label: str) -> bool:
    return label.upper().startswith("EPSG:")


def align_path_crs(xyz: np.ndarray, path_crs: str | None, path_vertical: str, cloud_crs: str,
                   cloud_vertical: str) -> tuple[np.ndarray, str | None, dict | None]:
    """The path's positions in the cloud's CRS, a warning when its CRS was converted or assumed, and what was
    converted. Raises ValueError when local and EPSG coordinates meet or the conversion fails."""
    if path_crs is None:
        return xyz, f"The delivered path declares no CRS; the cloud's CRS ({cloud_crs}) was assumed.", None
    path_local, cloud_local = not _is_epsg(path_crs), not _is_epsg(cloud_crs)
    if path_local and cloud_local:
        return xyz, None, None
    if path_local != cloud_local:
        raise ValueError(f"The delivered path is in {path_crs} coordinates and the point cloud in {cloud_crs}; "
                         "they cannot be aligned.")
    if path_crs.upper() == cloud_crs.upper():
        return xyz, None, None
    from pyproj import Transformer
    from pyproj.exceptions import ProjError

    vertical = _is_epsg(path_vertical) and _is_epsg(cloud_vertical)
    source, target = (f"{path_crs}+{path_vertical}", f"{cloud_crs}+{cloud_vertical}") if vertical else (path_crs, cloud_crs)
    try:
        transformer = Transformer.from_crs(source, target, always_xy=True)
        if vertical:
            x, y, z = transformer.transform(xyz[:, 0], xyz[:, 1], xyz[:, 2], errcheck=True)
        else:
            (x, y), z = transformer.transform(xyz[:, 0], xyz[:, 1], errcheck=True), xyz[:, 2]
    except ProjError as error:
        raise ValueError(f"The delivered path's CRS ({path_crs}) could not be converted to the cloud's CRS "
                         f"({cloud_crs}): {error}") from error
    out = np.c_[x, y, z].astype(np.float64)
    if not np.isfinite(out).all():
        raise ValueError(f"The delivered path's CRS ({path_crs}) could not be converted to the cloud's CRS ({cloud_crs}).")
    return (out, f"The delivered path's CRS ({path_crs}) was converted to the cloud's CRS ({cloud_crs}).",
            {"from": source, "to": target, "vertical": vertical})


def read_trajectory_clock(path: Path) -> tuple[str, int | None]:
    """The ``# time_system`` and ``# gps_week`` lines of a 3Dtrees trajectory CSV."""
    meta = read_trajectory_meta(path)
    week = meta.get("gps_week", "")
    return meta.get("time_system", "none") or "none", int(week) if week.isdigit() else None


def path_to_cloud_clock(time: np.ndarray, system: str, week: int | None, cloud_system: str) -> np.ndarray | None:
    """The path's times on the cloud's clock (gps_adjusted or gps_sow), or None when they can't be converted."""
    if system == cloud_system or cloud_system not in ("gps_adjusted", "gps_sow"):
        return None
    if system == "unix_utc":
        starts, offsets = np.array(LEAP_SECONDS).T
        index = np.searchsorted(starts, time, side="right") - 1
        if (index < 0).any():
            return None
        gps = time - GPS_EPOCH_UNIX + offsets[index]
    elif system == "gps_unix":
        gps = time - GPS_EPOCH_UNIX
    elif system == "gps_adjusted":
        gps = time + 1e9
    elif system == "gps_sow" and week is not None:
        gps = week * WEEK_S + time
    else:
        return None
    if cloud_system == "gps_adjusted":
        return gps - 1e9
    return gps - math.floor(gps[0] / WEEK_S) * WEEK_S


def _clock_share(time: np.ndarray, cloud_span: tuple[float, float]) -> float:
    low, high = cloud_span
    overlap = min(high, time[-1] + TIME_MARGIN) - max(low, time[0] - TIME_MARGIN)
    return max(0.0, overlap) / max(high - low, TIME_MARGIN)


def align_path_clock(time: np.ndarray, system: str, week: int | None, cloud_system: str,
                     cloud_span: tuple[float, float]) -> tuple[np.ndarray, dict | None]:
    """Use the path's times as they are; if they cover less than half the cloud's time span, try the path's
    declared clock converted to the cloud's. Returns the times and, when converted, what was converted."""
    raw = _clock_share(time, cloud_span)
    if raw >= PATH_CLOCK_MIN_SHARE:
        return time, None
    converted = path_to_cloud_clock(time, system, week, cloud_system)
    if converted is None or _clock_share(converted, cloud_span) < PATH_CLOCK_MIN_SHARE:
        return time, None
    return converted, {"from": system, "to": cloud_system, "raw_share": round(raw, 4),
                       "converted_share": round(_clock_share(converted, cloud_span), 4)}


def _grid(paths: list[Path], voxel: float) -> tuple[Grid, str | None]:
    lo, hi = np.full(3, np.inf), np.full(3, -np.inf)
    system, crs = None, None
    point_count = 0
    for path in paths:
        with open_cloud(path) as reader:
            header = reader.header
            if not header.point_count:
                continue
            point_count += header.point_count
            file_system = _time_system(header, set(header.point_format.dimension_names))
            if file_system == "none":
                raise ValueError(f"{path.name} has no GPS time, so its pulses can't be placed on the path")
            if system is not None and system != file_system:
                raise ValueError("Input clouds must use the same GPS time system.")
            system = file_system
            file_crs = _crs(header)
            if crs is not None and file_crs is not None and crs != file_crs:
                raise ValueError("Input clouds must use the same coordinate reference system.")
            crs = crs or file_crs
            lo = np.minimum(lo, header.mins)
            hi = np.maximum(hi, header.maxs)
    if not point_count:
        raise ValueError("The cloud has no points.")
    crs = cloud_crs(crs, lo, hi)
    lo = np.floor(lo / voxel) * voxel
    # Include returns at the maximum, even when it falls on a voxel face.
    shape = tuple(math.floor(float(n)) + 1 for n in np.round((hi - lo) / voxel, 9))
    return Grid(lo, shape, voxel, point_count), crs.to_wkt() if crs else None


def validate_trajectory(time: np.ndarray, xyz: np.ndarray) -> None:
    if len(time) < 2 or xyz.shape != (len(time), 3):
        raise ValueError("Occlusion needs a path with at least two positions.")
    if not np.isfinite(time).all() or not np.isfinite(xyz).all() or not np.all(time[1:] > time[:-1]):
        raise ValueError("Occlusion needs finite positions and strictly increasing GPS times.")
    if not np.any(np.ptp(xyz, axis=0) > 0):
        raise ValueError("Occlusion needs a moving sensor path; static (TLS) scans are not supported.")


def _positions(t: np.ndarray, time: np.ndarray, xyz: np.ndarray) -> tuple:
    return tuple(np.interp(t, time, xyz[:, k]) for k in range(3))


def delivered_miss(pulses: Pulses, time: np.ndarray, xyz: np.ndarray) -> dict:
    """Median and p95 perpendicular miss to the sampled usable return lines."""
    indices = np.flatnonzero((pulses.returns > 1) & (pulses.time >= time[0] - TIME_MARGIN)
                             & (pulses.time <= time[-1] + TIME_MARGIN))
    indices = indices[::max(1, -(-len(indices) // 1_000_000))]
    first, last = pulses.first[indices].astype(float), pulses.last[indices].astype(float)
    delta = last - first
    length = np.linalg.norm(delta, axis=1)
    usable = length >= 0.5
    sensor = np.column_stack(_positions(pulses.time[indices][usable], time, xyz)) - pulses.origin
    miss = np.linalg.norm(np.cross(sensor - first[usable], delta[usable]), axis=1) / length[usable]
    if not len(miss):
        return dict(usable_pulses=0, threshold_m=0.1, median_m=None, p95_m=None, threshold_exceeded=None)
    median, p95 = np.percentile(miss, [50, 95])
    return dict(usable_pulses=len(miss), threshold_m=0.1, median_m=float(median), p95_m=float(p95),
                threshold_exceeded=bool(median > 0.1))


def _time_mask(t, time):
    return (t >= time[0] - TIME_MARGIN) & (t <= time[-1] + TIME_MARGIN)


def classify(nhit: np.ndarray, nmiss: np.ndarray, nocc: np.ndarray) -> np.ndarray:
    """OccPy's rule: filled where hit, else empty where passed, else occluded where hidden, else unobserved."""
    out = np.full(nhit.shape, UNOBSERVED, np.uint8)
    out[nocc > 0] = OCCLUDED
    out[nmiss > 0] = EMPTY
    out[nhit > 0] = FILLED
    return out


def _fill_nearest(surface: np.ndarray) -> np.ndarray:
    from scipy import ndimage

    missing = np.isnan(surface)
    if missing.all():
        return surface
    if missing.any():
        index = ndimage.distance_transform_edt(missing, return_distances=False, return_indices=True)
        surface = surface[tuple(index)]
    return surface


def _xy_windows(grid: Grid):
    nx, ny, _ = grid.shape
    for x in range(0, nx, WINDOW):
        for stop in range(ny, 0, -WINDOW):
            yield slice(x, min(x + WINDOW, nx)), slice(max(0, stop - WINDOW), stop)


def _surface_indices(nhit: np.ndarray, grid: Grid) -> tuple[np.ndarray, np.ndarray]:
    low = np.empty(grid.shape[:2], np.int32)
    high = np.empty_like(low)
    for xy in _xy_windows(grid):
        hits = nhit[xy] > 0
        has_hits = hits.any(axis=2)
        low[xy] = np.where(has_hits, hits.argmax(axis=2), -1)
        high[xy] = np.where(has_hits, grid.shape[2] - 1 - hits[:, :, ::-1].argmax(axis=2), -1)
    return low, high


def _ground_from_lowest(lowest: np.ndarray, grid: Grid, raster_columns=None) -> tuple[np.ndarray, int]:
    """Ground from the lowest filled voxel per column with neighbourhood spike correction."""
    from scipy import ndimage

    has = lowest >= 0
    ground = np.where(has, grid.origin[2] + (lowest + 0.5) * grid.voxel, np.nan)
    ground = _fill_nearest(ground)
    size = max(3, int(round(SPIKE_WINDOW_M / grid.voxel)) | 1)
    median = ndimage.median_filter(ground, size=size, mode="nearest")
    spikes = has & (ground > median + SPIKE_M)
    if raster_columns is not None:
        spikes &= ~raster_columns
    ground[spikes] = median[spikes]
    return ground, int(spikes.sum())


@contextmanager
def open_dtm(path: Path, crs_wkt: str | None):
    """Open a terrain raster in the cloud coordinates, with a warning if assumed."""
    import rasterio
    from pyproj import CRS
    from rasterio.errors import RasterioError

    try:
        with rasterio.open(path) as raster:
            raster_crs = raster.crs
            if raster_crs is not None:
                left, bottom, right, top = raster.bounds
                reference = CRS.from_wkt(raster_crs.to_wkt())
                if is_local_geocentric(reference, ((left + right) / 2, (bottom + top) / 2)):
                    raster_crs = None
            warning = None
            if crs_wkt and raster_crs:
                if raster_crs != rasterio.crs.CRS.from_wkt(crs_wkt):
                    raise TerrainError(f"{path.name} and the point cloud must use the same coordinate reference system.")
            elif crs_wkt:
                warning = f"DTM {path.name} has no usable CRS; assuming point-cloud coordinates."
            elif raster_crs:
                warning = f"The point cloud has no usable CRS; assuming DTM {path.name} uses the same coordinates."
            yield raster, warning
    except (OSError, RasterioError) as error:
        raise TerrainError(f"Could not sample terrain raster {path.name}: {error}") from error


def _sample_dtm(path: Path, grid: Grid, crs_wkt: str | None) -> tuple[np.ndarray, list[str]]:
    """Sample the raster cell containing each voxel column centre."""
    from rasterio.transform import rowcol
    from rasterio.windows import Window

    nx, ny = grid.shape[:2]
    cx = grid.origin[0] + (np.arange(nx) + 0.5) * grid.voxel
    cy = grid.origin[1] + (np.arange(ny) + 0.5) * grid.voxel
    ground = np.full((nx, ny), np.nan)
    covered = False
    with open_dtm(path, crs_wkt) as (raster, warning):
        for xy in _xy_windows(grid):
            x, y = np.meshgrid(cx[xy[0]], cy[xy[1]], indexing="ij")
            rows, cols = rowcol(raster.transform, x, y)
            rows, cols = rows.reshape(x.shape), cols.reshape(x.shape)
            inside = (rows >= 0) & (rows < raster.height) & (cols >= 0) & (cols < raster.width)
            if not inside.any():
                continue
            covered = True
            rows, cols = rows[inside], cols[inside]
            row, col = int(rows.min()), int(cols.min())
            window = Window(col, row, int(cols.max()) - col + 1, int(rows.max()) - row + 1)
            values = raster.read(1, window=window, masked=True)[rows - row, cols - col]
            elevations = values.astype(np.float64).filled(np.nan) * raster.scales[0] + raster.offsets[0]
            ground[xy][inside] = np.where(np.isfinite(elevations), elevations, np.nan)
    if not covered:
        return ground, [f"DTM {path.name} does not cover any point-cloud voxel columns; using lowest hits."]
    return ground, [warning] if warning else []


def _stage(resources, name, **kwargs):
    return resources.stage(name, **kwargs) if resources is not None else nullcontext({})


def prepare_index(paths, voxel):
    for path in paths:
        with open_cloud(path, read_evlrs=False) as reader:
            if is_copc(reader.header):
                index = ScanIndex(paths, None)
                index.copc = True
                return index, None
    source_identity = _source_identity(paths)
    grid, crs_wkt = _grid(paths, voxel)
    nx, ny, nz = grid.shape
    if nx * ny * 96 > HARD_BUDGET - time_windows.WORKING_BYTES or nz > np.iinfo(np.int32).max:
        raise ValueError(f"The occlusion extent {grid.shape} exceeds the 64 GiB surface budget; "
                         "increase --voxel or restrict the acquisition.")
    index = ScanIndex(paths, grid)
    index.source_identity = source_identity
    return index, crs_wkt


def _source_identity(paths):
    return identify(paths, "source", np.empty(0), np.empty(0))


def _validate_source(index):
    if _source_identity(index.paths) != index.source_identity:
        raise ValueError("The LAS/LAZ collection changed during the run. Retry with stable source files.")


def _crop(grid, low, high, terrain=None, raster_columns=None):
    nx, ny, nz = grid.shape
    if (high < 0).all():
        raise PathCoverageError("The point cloud and path have no overlapping GPS times. "
                         "Check GPS week seconds versus adjusted GPS time.")
    low[low == nz] = -1
    if raster_columns is not None and raster_columns.all():
        ground, spikes = terrain, 0
    else:
        ground, spikes = _ground_from_lowest(low, grid, raster_columns)
        if raster_columns is not None:
            np.copyto(ground, terrain, where=raster_columns)
    top_z = grid.origin[2] + (int(high.max()) + 0.5) * grid.voxel
    offset = max(0, math.floor((float(ground.min()) - 1 - grid.origin[2]) / grid.voxel))
    offset = min(offset, nz - 1)
    stop = min(nz, max(offset + 1, math.ceil((top_z + 1 - grid.origin[2]) / grid.voxel)))
    origin = grid.origin.copy()
    origin[2] += offset * grid.voxel
    return Grid(origin, (nx, ny, stop - offset), grid.voxel), high - offset, ground, spikes


def _product_blocks(grid, top, ground):
    for xy in _xy_windows(grid):
        high = int(top[xy].max(initial=-1))
        if high < 0:
            continue
        # Include the voxel at the first ground crossing. The product mask
        # uses float32 heights, so keep one extra layer for boundary rounding.
        low = max(0, math.floor((float(ground[xy].min()) - grid.origin[2]) / grid.voxel) - 1)
        if low <= high:
            yield (*xy, slice(low, min(grid.shape[2], high + 2)))


def _copc_header(grid: Grid, crs_wkt: str | None) -> laspy.LasHeader:
    header = laspy.LasHeader(point_format=6, version="1.4")
    header.offsets = grid.origin
    header.scales = np.full(3, 0.001)
    header.add_extra_dims([laspy.ExtraBytesParams("occlusion", np.uint8, "1 filled 3 occluded 4 unobserved"),
                           laspy.ExtraBytesParams("nhit", np.uint32, "returns registered here"),
                           laspy.ExtraBytesParams("nmiss", np.uint32, "pulses passing through"),
                           laspy.ExtraBytesParams("nocc", np.uint32, "pulses passing behind returns"),
                           laspy.ExtraBytesParams("occluded_fraction", np.float32, "nocc / (nhit + nmiss + nocc)"),
                           laspy.ExtraBytesParams("height_above_ground", np.float32, "voxel centre above ground (m)")])
    if crs_wkt:
        from pyproj import CRS

        header.add_crs(CRS.from_wkt(crs_wkt))
    return header


def _index_copc(flat: Path, copc: Path, size: int, work: Path, *, tags: dict[str, str] | None = None) -> None:
    if not size:
        copc.write_bytes(b"")
        return
    env = dict(os.environ, PROJ_DATA=str(Path(UNTWINE).parents[1] / "share" / "proj"), TMPDIR=str(work))
    env.pop("PROJ_LIB", None)
    done = subprocess.run([UNTWINE, "--files", str(flat), "--output_dir", str(copc), "--temp_dir", str(work)],
                          env=env, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if done.returncode != 0 or not copc.is_file():
        if done.returncode in (-signal.SIGKILL, 128 + signal.SIGKILL):
            raise OutOfMemory("The Untwine COPC writer was killed (SIGKILL).")
        stderr = done.stderr.strip()
        raise_if_oom(RuntimeError(stderr))
        raise ValueError(f"could not write the occlusion COPC: {stderr[-500:]}")
    with copc.open("r+b" if tags else "rb") as stream:
        header = laspy.LasHeader.read_from(stream, read_evlrs=False)
        if header.point_count != size:
            raise ValueError("the occlusion COPC lost voxels")
        if tags:
            if header.version.minor != 4 or not header.number_of_evlrs:
                raise ValueError("the occlusion COPC has no LAS 1.4 extended records")
            # Untwine ends the file with the hierarchy EVLR. Append the tags
            # there so all existing COPC chunk and hierarchy offsets stay valid.
            record = laspy.VLR("3dtrees", 1, "3Dtrees occlusion metadata",
                               json.dumps(tags, separators=(",", ":"), sort_keys=True).encode("utf-8"))
            stream.seek(0, os.SEEK_END)
            VLRList([record]).write_to(stream, as_extended=True)
            stream.seek(LAS14_EVLR_COUNT_OFFSET)
            stream.write(pack("<I", header.number_of_evlrs + 1))


def _voxel_records(header, grid: Grid, xy, centres, heights, classes, inside, counts):
    keep = inside & (classes != EMPTY)
    ix, iy, iz = np.nonzero(keep)
    if not len(ix):
        return
    xy0 = grid.origin[:2] + np.array([part.start for part in xy]) * grid.voxel
    points = laspy.ScaleAwarePointRecord.zeros(len(ix), header=header)
    points.x = xy0[0] + (ix + 0.5) * grid.voxel
    points.y = xy0[1] + (iy + 0.5) * grid.voxel
    points.z = centres[iz]
    points["occlusion"] = classes[keep]
    points["height_above_ground"] = heights[keep]
    hit, miss, occ = (arr[keep] for arr in counts)
    points["nhit"], points["nmiss"], points["nocc"] = hit, miss, occ
    total = hit.astype(np.float64) + miss + occ
    points["occluded_fraction"] = np.divide(occ, total, out=np.full(len(ix), np.nan, np.float32), where=total > 0)
    return points


def _band_shares(classes: np.ndarray, first_layer: int, ground_layers: np.ndarray, inside: np.ndarray,
                 n_bands: int) -> tuple:
    nx, ny, _ = classes.shape
    x, y, z = np.nonzero(inside & (classes != UNOBSERVED))
    band = first_layer + z - ground_layers[x, y]
    occluded = classes[x, y, z] == OCCLUDED
    shares = np.full((n_bands + 1, nx, ny), np.nan, np.float32)
    # Each band contains exactly one voxel layer per column.
    shares[band, x, y] = occluded
    column = inside.sum(axis=2)
    unseen = (inside & (classes == UNOBSERVED)).sum(axis=2)
    np.divide(unseen, column, out=shares[-1], where=column > 0)
    seen = np.bincount(band, minlength=n_bands)
    hidden = np.bincount(band[occluded], minlength=n_bands)
    return shares.transpose(0, 2, 1)[:, ::-1, :], seen, hidden


def _ordered_blocks(blocks, classify_block, workers):
    from collections import deque
    from concurrent.futures import ThreadPoolExecutor
    from itertools import islice

    if workers == 1:
        yield from map(classify_block, blocks)
        return
    with ThreadPoolExecutor(workers) as pool:
        pending = deque(pool.submit(classify_block, block) for block in islice(blocks, workers))
        while pending:
            product = pending.popleft().result()
            if (block := next(blocks, None)) is not None:
                pending.append(pool.submit(classify_block, block))
            yield product


def _write_products(grid, layout, layer_ranges, top, ground, copc, cog, crs_wkt, tags, work, resources, slots,
                    raster_columns):
    import rasterio
    from rasterio.shutil import copy as copy_raster
    from rasterio.transform import from_origin
    from rasterio.windows import Window

    nx, ny, nz = grid.shape
    layers = np.arange(nz)
    centres = grid.origin[2] + (layers + 0.5) * grid.voxel
    # Lowest-hit estimates lie on voxel centres. Rounding preserves their band boundaries.
    ground_layers = np.rint((ground - centres[0]) / grid.voxel)
    if raster_columns is not None:
        # For an off-centre DTM, the first centre at or above ground belongs to band zero.
        nearest_centres = grid.origin[2] + (ground_layers + 0.5) * grid.voxel
        ground_layers += raster_columns & (nearest_centres < ground)
        ground_layers = np.minimum(ground_layers, nz)
        if np.any(ground_layers < -65_533):
            raise TerrainError("Terrain elevations are too far below the voxel grid for the COG height bands.")
    ground_layers = ground_layers.astype(np.int32)
    covered = top >= 0
    n_bands = int((top[covered] - ground_layers[covered]).max(initial=0)) + 1
    if raster_columns is not None and n_bands >= 65_535:
        raise TerrainError("Terrain elevations require too many COG height bands; increase --voxel.")
    names = [f"occluded_{i * grid.voxel:.12g}-{(i + 1) * grid.voxel:.12g}m" for i in range(n_bands)]
    names.append("unobserved_column")
    north = grid.origin[1] + ny * grid.voxel
    profile = dict(driver="GTiff", width=nx, height=ny, count=len(names), dtype="float32", nodata=np.nan,
                   transform=from_origin(grid.origin[0], north, grid.voxel, grid.voxel), crs=crs_wkt,
                   tiled=True, blockxsize=128, blockysize=128, compress="DEFLATE", predictor=3,
                   BIGTIFF="IF_SAFER", SPARSE_OK=True)
    wide_raster = len(names) > 32
    if wide_raster:
        # GDAL builds overviews band by band. Keep each read from decoding
        # the other bands of the same intermediate tile.
        profile["interleave"] = "band"
        # The COG writer compresses the final tiles. Skip the intermediate
        # compression pass when the uncompressed raster fits on scratch.
        if nx * ny * len(names) * 4 < shutil.disk_usage(work).free // 2:
            profile.update(compress="NONE", predictor=1)
    cache_mib = 512 if wide_raster else 128
    cog_threads = min(slots, 4) if wide_raster else 1
    flat = work / "voxels.las"
    raster = work / "voxel_shares.tif"
    classes_total = np.zeros(5, np.int64)
    seen_total, hidden_total = np.zeros(n_bands, np.int64), np.zeros(n_bands, np.int64)

    def classify_block(block):
        xs, ys, z = block
        counts = read_counts(block)
        hits = counts[0] > 0
        local_top = np.where(hits.any(axis=2), z.stop - 1 - hits[:, :, ::-1].argmax(axis=2), -1)
        heights = (centres[z] - ground[xs, ys, None]).astype(np.float32)
        threshold = np.minimum(top[xs, ys], local_top)
        below_canopy = layers[z] <= threshold[:, :, None]
        # COPC uses float32 heights while the COG uses integer height bands.
        copc_inside = (heights >= 0) & below_canopy
        ground_block = ground_layers[xs, ys]
        inside = (layers[z] >= ground_block[:, :, None]) & below_canopy
        if not (inside.any() or copc_inside.any()):
            return None
        classes = classify(*counts)
        totals = np.bincount(classes[inside], minlength=5)
        points = _voxel_records(header, grid, (xs, ys), centres[z], heights, classes, copc_inside, counts)
        shares, seen, hidden = _band_shares(classes, z.start, ground_block, inside, n_bands)
        window = Window(xs.start, ny - ys.stop, xs.stop - xs.start, ys.stop - ys.start)
        return window, points, shares, totals, seen, hidden

    # Count sums, classification, and record construction share read-only buffers.
    # Keep memory bounded and commit windows to both output streams in grid order.
    blocks = _product_blocks(grid, top, ground)
    header = _copc_header(grid, crs_wkt)
    bytes_per_worker = WINDOW**2 * (72 * nz + 4 * len(names))
    with rasterio.Env(GDAL_CACHEMAX=cache_mib * 1024**2, CPL_TMPDIR=str(work)):
        with laspy.open(flat, "w", header=header) as writer, rasterio.open(raster, "w", **profile) as out:
            out.descriptions = tuple(names)
            out.update_tags(OCCLUSION_BAND_M=str(grid.voxel), **tags,
                            OCCLUSION_VALUE="occluded / (filled + empty + occluded)")
            with _stage(resources, "classify voxels + write products", scratch=work.parent):
                read_counts = time_windows.count_blocks(work.parent, layout, layer_ranges)
                workers = min(max(1, slots - 1), max(1, available_memory() // bytes_per_worker))
                for product in _ordered_blocks(blocks, classify_block, workers):
                    if product is None:
                        continue
                    window, points, shares, totals, seen, hidden = product
                    if points is not None:
                        writer.write_points(points)
                    out.write(shares, window=window)
                    classes_total += totals
                    seen_total += seen
                    hidden_total += hidden
                del read_counts
            size = writer.header.point_count
        with _stage(resources, "COPC write", children=True, scratch=work.parent) as stage:
            stage["bytes_read"] = flat.stat().st_size
            _index_copc(flat, copc, size, work, tags=tags)
            stage["bytes_written"] = copc.stat().st_size
        with _stage(resources, "COG write", scratch=work.parent) as stage:
            stage["bytes_read"] = raster.stat().st_size
            stage.update(threads=cog_threads, bands=len(names))
            copy_raster(raster, cog, driver="COG", compress="DEFLATE", predictor=3, blocksize=512,
                        BIGTIFF="IF_SAFER", NUM_THREADS=str(cog_threads))
            stage["bytes_written"] = cog.stat().st_size
    stats = {name: int(classes_total[code]) for name, code in
             (("filled", FILLED), ("empty", EMPTY), ("occluded", OCCLUDED), ("unobserved", UNOBSERVED))}
    observed = sum(seen_total)
    return {
        "voxels_ground_to_canopy": stats,
        "occluded_share": round(stats["occluded"] / int(observed), 4) if observed else None,
        "cog_band_step_m": grid.voxel,
        "cog_bands": names,
        "occluded_share_by_height_band": [int(hidden) / int(seen) if seen else None
                                          for seen, hidden in zip(seen_total, hidden_total)],
        "copc_points": int(size),
    }


def run(paths: list[Path], time: np.ndarray, xyz: np.ndarray, copc: Path, cog: Path, voxel: float = 0.5,
        scratch: Path | None = None, source: str = "estimated", *,
        dtm: Path | None = None, quality: str | None = None, delivered_check: dict | None = None,
        path_method: str | None = None, path_quality: str | None = None, path_warnings=(),
        budget: int = DEFAULT_BUDGET, resources: RunResources | None = None,
        threads: int | None = None, prepared=None, single_files=(), pathset=None) -> dict:
    """Trace disjoint GPS-time windows and sum the voxel counts before classification. ``pathset`` (several
    delivered paths, each pulse segment on its own path) replaces ``time`` and ``xyz``."""
    if pathset is None:
        validate_trajectory(time, xyz)
        span = time[0], time[-1]
    else:
        span = pathset.time_range
        time, xyz = np.empty(0), np.empty((0, 3))
    scratch = scratch or copc.parent
    copc.parent.mkdir(parents=True, exist_ok=True)
    cog.parent.mkdir(parents=True, exist_ok=True)
    slots = decoding_threads(threads)
    budget = min(budget, HARD_BUDGET)
    if prepared is None:
        index, crs_wkt = prepare_index(paths, voxel)
        if not index.copc:
            index.populate()
    else:
        index, crs_wkt = prepared
    if index.copc:
        raise ValueError("Occlusion requires scanner-order LAS/LAZ. Supply the standardized files instead of COPC.")
    _validate_source(index)
    full_grid = index.grid
    time_bounds = span[0] - TIME_MARGIN, np.nextafter(span[1] + TIME_MARGIN, np.inf)
    with _stage(resources, "finalize extent/crop") as extent:
        low, high = index.surfaces(*time_bounds, extent)
        index.release_surfaces()
        terrain, terrain_warnings = _sample_dtm(dtm, full_grid, crs_wkt) if dtm is not None else (None, [])
        raster_columns = np.isfinite(terrain) if terrain is not None else None
        terrain_digest = hashlib.sha256(terrain).hexdigest() if terrain is not None else None
        grid, top, ground, spikes = _crop(full_grid, low, high, terrain, raster_columns)
        del terrain
        del low, high
        extent.update(grid_shape=list(grid.shape), grid_voxels=math.prod(grid.shape))
    ground_source = "lowest_hit"
    ground_shares = None
    if raster_columns is not None:
        raster_count = int(np.count_nonzero(raster_columns))
        total_columns = raster_columns.size
        ground_shares = {"dtm": raster_count / total_columns, "lowest_hit": (total_columns - raster_count) / total_columns}
        if raster_count:
            ground_source = "dtm" if raster_count == total_columns else "dtm+lowest_hit"
    layer_ranges = ray_walk.layers(grid, top, ground)
    with _stage(resources, "plan GPS-time windows") as planning:
        layout = time_windows.plan(index, grid, time, xyz, time_bounds, slots, min(budget, available_memory()), layer_ranges,
                                   pathset=pathset)
        planning.update(windows=len(layout), estimated_checkpoint_bytes=sum(12 * w.count_cells for w in layout))
    identity = identify(paths,
                        [str(copc.resolve()), str(cog.resolve()), voxel, source, quality, list(single_files),
                         [(w.lower, w.upper, [(s.start, s.stop) for s in w.bounds]) for w in layout],
                         *([terrain_digest] if terrain_digest is not None else []),
                         *([pathset.digest()] if pathset is not None else [])], time, xyz)
    work = scratch / f".occlusion-{identity}"
    with lock_job(work):
        for orphan in work.glob("merge-*"):
            if orphan.is_dir():
                shutil.rmtree(orphan)
        tags = dict(OCCLUSION_GROUND=ground_source, OCCLUSION_TRAJECTORY=source,
                    OCCLUSION_PATH_METHOD=path_method or source, OCCLUSION_QUALITY=quality or path_quality or "unavailable")
        if ground_source == "dtm+lowest_hit":
            tags["OCCLUSION_GROUND_DTM_SHARE"] = str(ground_shares["dtm"])
            tags["OCCLUSION_GROUND_LOWEST_HIT_SHARE"] = str(ground_shares["lowest_hit"])
        warnings = [*path_warnings, *terrain_warnings]
        if warnings:
            tags["OCCLUSION_WARNINGS"] = json.dumps(warnings, ensure_ascii=False, separators=(",", ":"))
        if delivered_check is not None:
            if delivered_check.get("crs_warning"):
                warnings.append(delivered_check["crs_warning"])
                tags["OCCLUSION_WARNINGS"] = json.dumps(warnings, ensure_ascii=False, separators=(",", ":"))
            for warning in delivered_check.get("warnings", ()):
                warnings.append(warning)
                tags["OCCLUSION_WARNINGS"] = json.dumps(warnings, ensure_ascii=False, separators=(",", ":"))
            if "coverage" in delivered_check:
                tags["OCCLUSION_DELIVERED_PATH_FILES"] = str(delivered_check["path_files"])
                tags["OCCLUSION_DELIVERED_PATH_COVERAGE"] = str(delivered_check["coverage"])
            if delivered_check.get("clock_converted"):
                clock = delivered_check["clock_converted"]
                warnings.append(f"The delivered path's times ({clock['from']}) did not match the cloud's GPS time; "
                                f"they were converted to {clock['to']}.")
                tags["OCCLUSION_WARNINGS"] = json.dumps(warnings, ensure_ascii=False, separators=(",", ":"))
            tags["OCCLUSION_DELIVERED_MISS_FLAG"] = delivered_miss_flag(delivered_check)
            flag = delivered_check["threshold_exceeded"]
            if flag is not None:
                tags["OCCLUSION_DELIVERED_MISS_MEDIAN_M"] = str(delivered_check["median_m"])
        specs = [(number, window.memory, (window.count_cells, 3), (window,)) for number, window in enumerate(layout)]
        specs.sort(key=lambda job: (-layout[job[0]].echoes, -job[1]))
        with _stage(resources, "traversal", children=True, scratch=work) as stage:
            windows = execute(work, specs, slots, time_windows.trace_window, time_windows.init_worker,
                              (index, grid, time, xyz, work, single_files, layer_ranges, time_bounds, pathset), stage)
            stage["windows"] = [windows[i] | {"estimated_worker_bytes": w.memory} for i, w in enumerate(layout)]
            stage["bytes_written"] = sum(w["bytes_written"] for w in windows.values())
        with TemporaryDirectory(dir=work, prefix="merge-") as folder, \
             TemporaryDirectory(prefix=".occlusion-", dir=cog.parent) as cog_folder:
            merge, staged_cog = Path(folder), Path(cog_folder) / "occlusion.tif"
            staged_copc = merge / "voxels.copc.laz"
            products = _write_products(grid, layout, layer_ranges, top, ground, staged_copc, staged_cog,
                                       crs_wkt, tags, merge, resources, slots, raster_columns)
            _validate_source(index)
            os.replace(staged_copc, copc)
            os.replace(staged_cog, cog)
        shutil.rmtree(work)
    return {
        "method": "OccPy-compatible whole-grid traversal with compact window counts",
        "trajectory": source,
        "single_return_fallback_files": [str(paths[index]) for index in single_files],
        **({"quality": quality} if quality else {}),
        **({"delivered_path_check": delivered_check} if delivered_check is not None else {}),
        "voxel_m": voxel,
        "grid": {"origin": [round(float(v), 3) for v in grid.origin], "shape": list(grid.shape),
                 "uncropped_shape": list(full_grid.shape), "voxels": math.prod(grid.shape),
                 "uncropped_voxels": math.prod(full_grid.shape)},
        "time_windows": {"windows": len(layout), "budget_bytes": budget,
                         "peak_concurrent_workers": stage["peak_concurrent_workers"],
                         "reused_windows": stage["reused_windows"]},
        "ground": ground_source,
        **({"warnings": terrain_warnings} if terrain_warnings else {}),
        **({"ground_column_shares": ground_shares} if ground_shares is not None else {}),
        "ground_spikes_replaced": spikes,
        **{key: extent[key] for key in ("echoes_read", "echoes_outside_path")},
        **{key: sum(w[key] for w in windows.values()) for key in
           ("traversed_pulses", "missing_returns", "recovered_echoes", "unpaired_echoes")},
        **products,
    }
