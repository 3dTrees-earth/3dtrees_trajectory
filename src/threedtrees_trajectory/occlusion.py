"""Occlusion of a moving-sensor point cloud: which voxels the lidar saw as filled or empty, and which stayed hidden.

Every pulse is traced from the sensor (the trajectory interpolated at its GPS time) through its returns and on to the
grid edge with OccPy's voxel traversal (Kuekenbrink et al. 2017, the compiled ``raytr`` core of ``occpy_ls``).
Each voxel collects hits, passes (miss) and passes behind the last return (occluded); OccPy's rule turns the counts
into a class: filled (any hit), empty (passed, no hit), occluded (only behind returns), unobserved (no pulse).

Heights above ground come from a DTM when one is given, else from the lowest filled voxel of each column, with
columns more than 2 m above their 5 m neighbourhood replaced by it (no pulse reached the ground there).

Outputs: a COPC with one point per voxel between the ground and the canopy top that is not empty (class, counts,
occluded fraction, height above ground), and a COG with the voxel columns as pixels and one band per metre above
ground holding the occluded share of the observed voxels, plus a last band with the column's unobserved share.
"""
from __future__ import annotations

import csv
import math
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

import laspy
import numpy as np

from .pulses import _crs, _open, _time_system

FILLED, EMPTY, OCCLUDED, UNOBSERVED = 1, 2, 3, 4  # OccPy's classes
CHUNK = 5_000_000
# The wheel converts each native count grid to nested Python lists before
# NumPy can copy it. Bound both the cells and the number of column objects.
MAX_VOXELS = 1_000_000_000
MAX_COLUMNS = 3_000_000
WINDOW = 128  # XY block size for deriving and writing the products
BAND_M = 1.0  # height band of the COG
SPIKE_M = 2.0  # a lowest-hit ground this far above its neighbourhood is vegetation, not ground
SPIKE_WINDOW_M = 5.0
TIME_MARGIN = 1.0  # seconds a pulse may lie outside the trajectory's span
UNTWINE = os.environ.get("TRAJECTORY_UNTWINE", "/opt/untwine/bin/untwine")


@dataclass
class Grid:
    origin: np.ndarray  # min x, y, z of the grid
    shape: tuple[int, int, int]
    voxel: float


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


def _grid(paths: list[Path], voxel: float) -> tuple[Grid, str | None]:
    lo, hi = np.full(3, np.inf), np.full(3, -np.inf)
    system, crs = None, None
    point_count = 0
    for path in paths:
        with laspy.open(path) as reader:
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
    if point_count > np.iinfo(np.int32).max:
        raise ValueError("The cloud exceeds OccPy's 32-bit echo count limit; process smaller acquisitions.")
    lo = np.floor(lo / voxel) * voxel
    # Include returns at the maximum, even when it falls on a voxel face.
    shape = tuple(math.floor(float(n)) + 1 for n in np.round((hi - lo) / voxel, 9))
    if math.prod(shape) > MAX_VOXELS or shape[0] * shape[1] > MAX_COLUMNS:
        raise ValueError(f"The occlusion grid {shape} is too large for OccPy. Increase --voxel or use a smaller area.")
    return Grid(lo, shape, voxel), crs.to_wkt() if crs else None


def _single_return(reader: laspy.LasReader) -> bool:
    by_return = np.asarray(reader.header.number_of_points_by_return)
    if by_return.sum() == reader.header.point_count:
        return not by_return[1:].any()
    # Some LAS writers omit counts or leave them stale after filtering.
    try:
        return not any((np.asarray(chunk.return_number) > 1).any() or
                       (np.asarray(chunk.number_of_returns) > 1).any() for chunk in reader.chunk_iterator(CHUNK))
    finally:
        reader.seek(0)


def _positions(t: np.ndarray, time: np.ndarray, xyz: np.ndarray) -> tuple:
    return tuple(np.interp(t, time, xyz[:, k]) for k in range(3))


def _trace_file(path: Path, grid: Grid, time: np.ndarray, xyz: np.ndarray, counters: dict) -> tuple:
    from raytr import PyRaytracer

    tracer = PyRaytracer()
    hi = grid.origin + np.array(grid.shape) * grid.voxel
    tracer.defineGrid(grid.origin, hi, *grid.shape, grid.voxel)
    span = (time[0] - TIME_MARGIN, time[-1] + TIME_MARGIN)
    with _open(path) as reader:
        single = _single_return(reader)
        deferred = False
        previous = -np.inf
        for chunk in reader.chunk_iterator(CHUNK):
            t = np.asarray(chunk.gps_time, np.float64)
            if not np.isfinite(t).all():
                raise ValueError(f"{path.name}: GPS times must be finite.")
            keep = (t >= span[0]) & (t <= span[1])
            counters["echoes_read"] += len(t)
            counters["echoes_outside_path"] += int((~keep).sum())
            if not keep.any():
                continue
            t = np.ascontiguousarray(t[keep])
            x, y, z = (np.ascontiguousarray(np.asarray(axis)[keep], np.float64) for axis in (chunk.x, chunk.y, chunk.z))
            sx, sy, sz = _positions(t, time, xyz)
            if single:
                tracer.doRaytracing_singleReturnPulses(x, y, z, sx, sy, sz, t)
                continue
            if not deferred and (np.any(t[1:] < t[:-1]) or t[0] < previous):
                deferred = True
            previous = t[-1]
            tracer.addPointData(x, y, z, sx, sy, sz, t,
                                np.ascontiguousarray(np.asarray(chunk.return_number)[keep]),
                                np.ascontiguousarray(np.asarray(chunk.number_of_returns)[keep]))
            if not deferred:
                tracer.doRaytracing()
                tracer.clearPulseDataset()
        if deferred:
            tracer.doRaytracing()
            tracer.clearPulseDataset()
    counters["traversed_pulses"] += int(tracer.get_num_traversed_pulses())
    counters["missing_returns"] += int(tracer.get_num_missing_returns())
    return (np.asarray(tracer.getNhit(), np.int32), np.asarray(tracer.getNmiss(), np.int32),
            np.asarray(tracer.getNocc(), np.int32))


def _trace_grid(paths: list[Path], grid: Grid, time: np.ndarray, xyz: np.ndarray, counters: dict) -> tuple:
    # Use the acquisition bounds for traversal. Clipping per XY tile can
    # change how OccPy handles a return before entry or off the fitted ray.
    # Each file has independent pulse state, including incomplete echoes.
    counts = _trace_file(paths[0], grid, time, xyz, counters)
    for path in paths[1:]:
        for total, values in zip(counts, _trace_file(path, grid, time, xyz, counters)):
            total += values
    return counts


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


def ground_from_dtm(path: Path, grid: Grid, crs_wkt: str | None = None) -> np.ndarray:
    """DTM height at every column centre (nearest cell); columns outside the DTM take their nearest neighbour."""
    import rasterio

    nx, ny, _ = grid.shape
    cx = grid.origin[0] + (np.arange(nx) + 0.5) * grid.voxel
    cy = grid.origin[1] + (np.arange(ny) + 0.5) * grid.voxel
    with rasterio.open(path) as dtm:
        if crs_wkt and dtm.crs and dtm.crs != rasterio.crs.CRS.from_wkt(crs_wkt):
            raise ValueError(f"{path.name} and the point cloud must use the same coordinate reference system.")
        rows, cols = rasterio.transform.rowcol(dtm.transform, *np.meshgrid(cx, cy, indexing="ij"))
        rows, cols = rows.reshape(nx, ny), cols.reshape(nx, ny)
        inside = (rows >= 0) & (rows < dtm.height) & (cols >= 0) & (cols < dtm.width)
        ground = np.full((nx, ny), np.nan)
        if inside.any():
            r, c = rows[inside], cols[inside]
            row, col = int(r.min()), int(c.min())
            window = rasterio.windows.Window(col, row, int(c.max()) - col + 1, int(r.max()) - row + 1)
            band = dtm.read(1, window=window, masked=True).astype(np.float64).filled(np.nan)
            band[~np.isfinite(band)] = np.nan
            ground[inside] = band[r - row, c - col]
    if np.isnan(ground).all():
        raise ValueError(f"{path.name} does not cover the point cloud with valid elevations.")
    return _fill_nearest(ground)


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


def _ground_from_lowest(lowest: np.ndarray, grid: Grid) -> tuple[np.ndarray, int]:
    """Ground from the lowest filled voxel per column with neighbourhood spike correction."""
    from scipy import ndimage

    has = lowest >= 0
    ground = np.where(has, grid.origin[2] + (lowest + 0.5) * grid.voxel, np.nan)
    ground = _fill_nearest(ground)
    size = max(3, int(round(SPIKE_WINDOW_M / grid.voxel)) | 1)
    median = ndimage.median_filter(ground, size=size, mode="nearest")
    spikes = has & (ground > median + SPIKE_M)
    ground[spikes] = median[spikes]
    return ground, int(spikes.sum())


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


def _index_copc(flat: Path, copc: Path, size: int, work: Path) -> None:
    if not size:
        copc.write_bytes(b"")
        return
    env = dict(os.environ, PROJ_DATA=str(Path(UNTWINE).parents[1] / "share" / "proj"), TMPDIR=str(work))
    env.pop("PROJ_LIB", None)
    done = subprocess.run([UNTWINE, "--files", str(flat), "--output_dir", str(copc), "--temp_dir", str(work)],
                          env=env, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if done.returncode != 0 or not copc.is_file():
        raise ValueError(f"could not write the occlusion COPC: {done.stderr.strip()[-500:]}")
    with laspy.open(copc) as check:
        if check.header.point_count != size:
            raise ValueError("the occlusion COPC lost voxels")


def _write_voxels(writer, grid: Grid, xy, centres, heights, classes, inside, counts) -> None:
    keep = inside & (classes != EMPTY)
    ix, iy, iz = np.nonzero(keep)
    if not len(ix):
        return
    xy0 = grid.origin[:2] + np.array([part.start for part in xy]) * grid.voxel
    points = laspy.ScaleAwarePointRecord.zeros(len(ix), header=writer.header)
    points.x = xy0[0] + (ix + 0.5) * grid.voxel
    points.y = xy0[1] + (iy + 0.5) * grid.voxel
    points.z = centres[iz]
    points["occlusion"] = classes[keep]
    points["height_above_ground"] = heights[keep]
    hit, miss, occ = (arr[keep] for arr in counts)
    points["nhit"], points["nmiss"], points["nocc"] = hit, miss, occ
    total = hit.astype(np.float64) + miss + occ
    points["occluded_fraction"] = np.divide(occ, total, out=np.full(len(ix), np.nan, np.float32), where=total > 0)
    writer.write_points(points)


def _band_shares(classes: np.ndarray, heights: np.ndarray, inside: np.ndarray, n_bands: int) -> tuple:
    nx, ny, _ = classes.shape
    x, y, z = np.nonzero(inside & (classes != UNOBSERVED))
    band = (heights[x, y, z] / BAND_M).astype(np.intp)
    key = (band * nx + x) * ny + y
    shape = (n_bands, nx, ny)
    seen = np.bincount(key, minlength=math.prod(shape)).reshape(shape)
    hidden = np.bincount(key[classes[x, y, z] == OCCLUDED], minlength=math.prod(shape)).reshape(shape)
    shares = np.full((n_bands + 1, nx, ny), np.nan, np.float32)
    np.divide(hidden, seen, out=shares[:-1], where=seen > 0)
    column = inside.sum(axis=2)
    unseen = (inside & (classes == UNOBSERVED)).sum(axis=2)
    np.divide(unseen, column, out=shares[-1], where=column > 0)
    return shares.transpose(0, 2, 1)[:, ::-1, :], seen.sum(axis=(1, 2)), hidden.sum(axis=(1, 2))


def _write_products(grid, counts, top, ground, copc, cog, crs_wkt, source, ground_source, work):
    import rasterio
    from rasterio.shutil import copy as copy_raster
    from rasterio.transform import from_origin
    from rasterio.windows import Window

    nx, ny, nz = grid.shape
    layers = np.arange(nz)
    centres = grid.origin[2] + (layers + 0.5) * grid.voxel
    max_height = (centres[np.maximum(top, 0)] - ground)[top >= 0].astype(np.float32).max(initial=0)
    n_bands = math.floor(float(max_height) / BAND_M) + 1
    names = [f"occluded_{b:g}-{b + BAND_M:g}m" for b in np.arange(n_bands) * BAND_M] + ["unobserved_column"]
    north = grid.origin[1] + ny * grid.voxel
    profile = dict(driver="GTiff", width=nx, height=ny, count=len(names), dtype="float32", nodata=np.nan,
                   transform=from_origin(grid.origin[0], north, grid.voxel, grid.voxel), crs=crs_wkt,
                   tiled=True, blockxsize=128, blockysize=128, compress="DEFLATE", predictor=3,
                   BIGTIFF="IF_SAFER", SPARSE_OK=True)
    flat = work / "voxels.las"
    raster = work / "voxel_shares.tif"
    classes_total = np.zeros(5, np.int64)
    seen_total, hidden_total = np.zeros(n_bands, np.int64), np.zeros(n_bands, np.int64)
    with rasterio.Env(GDAL_CACHEMAX=128 * 1024**2, CPL_TMPDIR=str(work)):
        with laspy.open(flat, "w", header=_copc_header(grid, crs_wkt)) as writer, rasterio.open(raster, "w", **profile) as out:
            out.descriptions = tuple(names)
            out.update_tags(OCCLUSION_BAND_M=str(BAND_M), OCCLUSION_GROUND=ground_source, OCCLUSION_TRAJECTORY=source,
                            OCCLUSION_VALUE="occluded / (filled + empty + occluded)")
            for xy in _xy_windows(grid):
                local_counts = tuple(count[xy] for count in counts)
                heights = (centres - ground[xy][:, :, None]).astype(np.float32)
                inside = (heights >= 0) & (layers <= top[xy][:, :, None])
                if not inside.any():
                    continue
                classes = classify(*local_counts)
                classes_total += np.bincount(classes[inside], minlength=5)
                _write_voxels(writer, grid, xy, centres, heights, classes, inside, local_counts)
                shares, seen, hidden = _band_shares(classes, heights, inside, n_bands)
                xs, ys = xy
                out.write(shares, window=Window(xs.start, ny - ys.stop, xs.stop - xs.start, ys.stop - ys.start))
                seen_total += seen
                hidden_total += hidden
            size = writer.header.point_count
        copy_raster(raster, cog, driver="COG", compress="DEFLATE", predictor=3, blocksize=512, BIGTIFF="IF_SAFER")
    _index_copc(flat, copc, size, work)
    stats = {name: int(classes_total[code]) for name, code in
             (("filled", FILLED), ("empty", EMPTY), ("occluded", OCCLUDED), ("unobserved", UNOBSERVED))}
    observed = sum(seen_total)
    return {
        "voxels_ground_to_canopy": stats,
        "occluded_share": round(stats["occluded"] / int(observed), 4) if observed else None,
        "cog_bands": names,
        "occluded_share_by_height_band": [int(hidden) / int(seen) if seen else None
                                          for seen, hidden in zip(seen_total, hidden_total)],
        "copc_points": int(size),
    }


def run(paths: list[Path], time: np.ndarray, xyz: np.ndarray, copc: Path, cog: Path, voxel: float = 0.5,
        dtm: Path | None = None, scratch: Path | None = None, source: str = "estimated") -> dict:
    """Trace, classify and write both products; returns the report section."""
    if len(time) < 2 or xyz.shape != (len(time), 3):
        raise ValueError("Occlusion needs a path with at least two positions.")
    if not np.isfinite(time).all() or not np.isfinite(xyz).all() or not np.all(time[1:] > time[:-1]):
        raise ValueError("Occlusion needs finite positions and strictly increasing GPS times.")
    if not np.any(np.ptp(xyz, axis=0) > 0):
        raise ValueError("Occlusion needs a moving sensor path; static (TLS) scans are not supported.")
    grid, crs_wkt = _grid(paths, voxel)
    counters = {"echoes_read": 0, "echoes_outside_path": 0, "traversed_pulses": 0, "missing_returns": 0}
    counts = _trace_grid(paths, grid, time, xyz, counters)
    if counters["echoes_read"] == counters["echoes_outside_path"]:
        raise ValueError("The point cloud and path have no overlapping GPS times. Check GPS week seconds versus adjusted GPS time.")

    low, top = _surface_indices(counts[0], grid)
    if dtm is not None:
        ground, spikes, ground_source = ground_from_dtm(dtm, grid, crs_wkt), 0, "dtm"
    else:
        (ground, spikes), ground_source = _ground_from_lowest(low, grid), "lowest_hit"
    copc.parent.mkdir(parents=True, exist_ok=True)
    cog.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".occlusion-", dir=scratch or copc.parent) as folder, \
         TemporaryDirectory(prefix=".occlusion-", dir=cog.parent) as cog_folder:
        work = Path(folder)
        staged_copc = work / "voxels.copc.laz"
        staged_cog = Path(cog_folder) / "occlusion.tif"
        products = _write_products(grid, counts, top, ground, staged_copc, staged_cog, crs_wkt,
                                   source, ground_source, work)
        os.replace(staged_copc, copc)
        os.replace(staged_cog, cog)
    return {
        "method": "OccPy voxel traversal (occpy_ls 0.3.0 raytr core)",
        "trajectory": source,
        "voxel_m": voxel,
        "grid": {"origin": [round(float(v), 3) for v in grid.origin], "shape": list(grid.shape)},
        "ground": ground_source,
        "ground_spikes_replaced": spikes,
        **counters,
        **products,
    }
