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

import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import laspy
import numpy as np

from .pulses import _open

FILLED, EMPTY, OCCLUDED, UNOBSERVED = 1, 2, 3, 4  # OccPy's classes
CHUNK = 5_000_000
MAX_VOXELS = 300_000_000  # per tile: OccPy keeps three int32 grids, about 3.6 GB at this size
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
    rows = [line for line in path.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#")]
    if not rows or rows[0].split(",")[:4] != ["time", "x", "y", "z"]:
        raise ValueError(f"{path.name}: expected a 3Dtrees trajectory CSV with time,x,y,z columns")
    values = [row.split(",")[:4] for row in rows[1:]]
    if not values or any(not v[0] for v in values):
        raise ValueError(f"{path.name}: a path with a time for every position is needed for occlusion")
    data = np.array(values, float)
    order = np.argsort(data[:, 0], kind="stable")
    return data[order, 0], data[order, 1:4]


def _grid(paths: list[Path], voxel: float) -> Grid:
    lo, hi = np.full(3, np.inf), np.full(3, -np.inf)
    for path in paths:
        with laspy.open(path) as reader:
            lo = np.minimum(lo, reader.header.mins)
            hi = np.maximum(hi, reader.header.maxs)
    lo = np.floor(lo / voxel) * voxel
    shape = tuple(int(n) for n in np.maximum(np.ceil((hi - lo) / voxel - 1e-9), 1))
    return Grid(lo, shape, voxel)


def _tiles(grid: Grid) -> list[tuple[int, int, int, int]]:
    """Index ranges (x0, x1, y0, y1) of XY tiles that keep each traversal grid under MAX_VOXELS."""
    nx, ny, nz = grid.shape
    side = max(1, int(np.sqrt(MAX_VOXELS / nz)))
    return [(x0, min(x0 + side, nx), y0, min(y0 + side, ny))
            for x0 in range(0, nx, side) for y0 in range(0, ny, side)]


def _single_return(reader: laspy.LasReader) -> bool:
    by_return = np.asarray(reader.header.number_of_points_by_return)
    return by_return.sum() > 0 and by_return[1:].sum() == 0


def _trace_tile(paths, grid: Grid, tile, time, xyz, counters) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    from raytr import PyRaytracer

    x0, x1, y0, y1 = tile
    lo = grid.origin + np.array([x0, y0, 0]) * grid.voxel
    shape = (x1 - x0, y1 - y0, grid.shape[2])
    hi = lo + np.array(shape) * grid.voxel
    tracer = PyRaytracer()
    tracer.defineGrid(lo.astype(np.float64), hi.astype(np.float64), *shape, float(grid.voxel))
    span = (time[0] - TIME_MARGIN, time[-1] + TIME_MARGIN)
    for path in paths:
        with _open(path) as reader:
            single = _single_return(reader)
            names = set(reader.header.point_format.dimension_names)
            if "gps_time" not in names:
                raise ValueError(f"{Path(path).name} has no GPS time, so its pulses can't be placed on the path")
            deferred = False
            previous = -np.inf
            for chunk in reader.chunk_iterator(CHUNK):
                t = np.asarray(chunk.gps_time, np.float64)
                keep = (t >= span[0]) & (t <= span[1])
                counters["outside_path"] += int((~keep).sum())
                if not keep.any():
                    continue
                t = np.ascontiguousarray(t[keep])
                x = np.ascontiguousarray(np.asarray(chunk.x)[keep], np.float64)
                y = np.ascontiguousarray(np.asarray(chunk.y)[keep], np.float64)
                z = np.ascontiguousarray(np.asarray(chunk.z)[keep], np.float64)
                sx, sy, sz = (np.ascontiguousarray(np.interp(t, time, xyz[:, k])) for k in range(3))
                if single:
                    tracer.doRaytracing_singleReturnPulses(x, y, z, sx, sy, sz, t)
                    continue
                if not deferred and (np.any(np.diff(t) < 0) or t[0] < previous):
                    deferred = True  # unsorted: hold the file's pulses and trace them once at its end
                previous = t[-1]
                tracer.addPointData(x, y, z, sx, sy, sz, t,
                                    np.ascontiguousarray(np.asarray(chunk.return_number)[keep]),
                                    np.ascontiguousarray(np.asarray(chunk.number_of_returns)[keep]))
                if not deferred:
                    tracer.doRaytracing()
                    tracer.clearPulseDataset()
            if not single:
                tracer.doRaytracing()
                tracer.clearPulseDataset()
    counters["traversed_pulses"] += int(tracer.get_num_traversed_pulses())
    counters["missing_returns"] += int(tracer.get_num_missing_returns())
    return (np.asarray(tracer.getNhit(), np.int32), np.asarray(tracer.getNmiss(), np.int32),
            np.asarray(tracer.getNocc(), np.int32))


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


def ground_from_hits(filled: np.ndarray, grid: Grid) -> tuple[np.ndarray, int]:
    """Lowest filled voxel centre per column; spikes over the neighbourhood median are replaced by it."""
    from scipy import ndimage

    has = filled.any(axis=2)
    lowest = np.where(has, filled.argmax(axis=2), 0)
    ground = np.where(has, grid.origin[2] + (lowest + 0.5) * grid.voxel, np.nan)
    ground = _fill_nearest(ground)
    size = max(3, int(round(SPIKE_WINDOW_M / grid.voxel)) | 1)
    median = ndimage.median_filter(ground, size=size, mode="nearest")
    spikes = has & (ground > median + SPIKE_M)
    ground[spikes] = median[spikes]
    return ground, int(spikes.sum())


def ground_from_dtm(path: Path, grid: Grid) -> np.ndarray:
    """DTM height at every column centre (nearest cell); columns outside the DTM take their nearest neighbour."""
    import rasterio

    nx, ny, _ = grid.shape
    cx = grid.origin[0] + (np.arange(nx) + 0.5) * grid.voxel
    cy = grid.origin[1] + (np.arange(ny) + 0.5) * grid.voxel
    with rasterio.open(path) as dtm:
        band = dtm.read(1, masked=True).astype(np.float64).filled(np.nan)
        rows, cols = rasterio.transform.rowcol(dtm.transform, *np.meshgrid(cx, cy, indexing="ij"))
        rows, cols = np.asarray(rows), np.asarray(cols)
        inside = (rows >= 0) & (rows < band.shape[0]) & (cols >= 0) & (cols < band.shape[1])
        ground = np.full((nx, ny), np.nan)
        ground[inside] = band[rows[inside], cols[inside]]
    if np.isnan(ground).all():
        raise ValueError(f"{path.name} does not cover the point cloud")
    return _fill_nearest(ground)


def _crs_wkt(paths: list[Path]) -> str | None:
    with laspy.open(paths[0]) as reader:
        try:
            crs = reader.header.parse_crs()
        except Exception:  # noqa: BLE001 - an unreadable CRS leaves the outputs without one
            return None
    return crs.to_wkt() if crs is not None else None


def write_copc(path: Path, points: dict, grid: Grid, crs_wkt: str | None, scratch: Path) -> None:
    """The kept voxels as a LAS 1.4 cloud, indexed into a COPC by untwine (which keeps the extra dimensions)."""
    header = laspy.LasHeader(point_format=6, version="1.4")
    header.offsets = grid.origin
    header.scales = np.full(3, 0.001)
    header.add_extra_dims([laspy.ExtraBytesParams("occlusion", np.uint8, "1 filled 3 occluded 4 unobserved"),
                           laspy.ExtraBytesParams("nhit", np.uint32, "pulses with a return here"),
                           laspy.ExtraBytesParams("nmiss", np.uint32, "pulses passing through"),
                           laspy.ExtraBytesParams("nocc", np.uint32, "pulses passing behind returns"),
                           laspy.ExtraBytesParams("occluded_fraction", np.float32, "nocc / (nhit + nmiss + nocc)"),
                           laspy.ExtraBytesParams("height_above_ground", np.float32, "voxel centre above ground (m)")])
    if crs_wkt:
        from pyproj import CRS

        header.add_crs(CRS.from_wkt(crs_wkt))
    las = laspy.LasData(header)
    las.x, las.y, las.z = points["x"], points["y"], points["z"]
    for name in ("occlusion", "nhit", "nmiss", "nocc", "occluded_fraction", "height_above_ground"):
        las[name] = points[name]
    flat = scratch / "voxels.las"
    las.write(flat)
    partial = Path(f"{path}.partial.copc.laz")
    partial.unlink(missing_ok=True)
    env = dict(os.environ, PROJ_DATA=str(Path(UNTWINE).parents[1] / "share" / "proj"))
    env.pop("PROJ_LIB", None)
    done = subprocess.run([UNTWINE, "--files", str(flat), "--output_dir", str(partial), "--temp_dir", str(scratch)],
                          env=env, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if done.returncode != 0 or not partial.is_file():
        raise ValueError(f"could not write the occlusion COPC: {done.stderr.strip()[-500:]}")
    with laspy.open(partial) as check:
        if check.header.point_count != len(points["x"]):
            raise ValueError("the occlusion COPC lost voxels")
    os.replace(partial, path)


def write_cog(path: Path, bands: np.ndarray, names: list[str], grid: Grid, crs_wkt: str | None, tags: dict) -> None:
    """Bands (count, rows, cols) north-up, float32, NaN as nodata, written as a Cloud Optimized GeoTIFF."""
    import rasterio
    from rasterio.transform import from_origin

    _, rows, cols = bands.shape
    top = grid.origin[1] + grid.shape[1] * grid.voxel
    profile = dict(driver="COG", width=cols, height=rows, count=len(bands), dtype="float32", nodata=np.nan,
                   transform=from_origin(grid.origin[0], top, grid.voxel, grid.voxel),
                   crs=crs_wkt, compress="DEFLATE", predictor=3, blocksize=512, BIGTIFF="IF_SAFER")
    partial = Path(f"{path}.partial.tif")
    with rasterio.open(partial, "w", **profile) as out:
        for index, (band, name) in enumerate(zip(bands, names), start=1):
            out.write(band.astype(np.float32), index)
            out.set_band_description(index, name)
        out.update_tags(**{key: json.dumps(value) if not isinstance(value, str) else value
                           for key, value in tags.items()})
    os.replace(partial, path)


def run(paths: list[Path], time: np.ndarray, xyz: np.ndarray, copc: Path, cog: Path, voxel: float = 0.5,
        dtm: Path | None = None, scratch: Path | None = None, source: str = "estimated") -> dict:
    """Trace, classify and write both products; returns the report section."""
    if len(time) < 2:
        raise ValueError("occlusion needs a path with at least two positions")
    grid = _grid(paths, voxel)
    counters = {"outside_path": 0, "traversed_pulses": 0, "missing_returns": 0}
    nx, ny, nz = grid.shape
    nhit = np.zeros(grid.shape, np.int32)
    nmiss = np.zeros(grid.shape, np.int32)
    nocc = np.zeros(grid.shape, np.int32)
    tiles = _tiles(grid)
    for tile in tiles:
        x0, x1, y0, y1 = tile
        h, m, o = _trace_tile(paths, grid, tile, time, xyz, counters)
        nhit[x0:x1, y0:y1], nmiss[x0:x1, y0:y1], nocc[x0:x1, y0:y1] = h, m, o
    counters["outside_path"] //= len(tiles)  # every tile reads every pulse
    classes = classify(nhit, nmiss, nocc)
    filled = classes == FILLED
    if dtm is not None:
        ground, spikes, ground_source = ground_from_dtm(dtm, grid), 0, "dtm"
    else:
        (ground, spikes), ground_source = ground_from_hits(filled, grid), "lowest_hit"
    has = filled.any(axis=2)
    top = np.where(has, nz - 1 - filled[:, :, ::-1].argmax(axis=2), -1)  # highest filled voxel index
    centres = grid.origin[2] + (np.arange(nz) + 0.5) * grid.voxel
    hag = (centres[None, None, :] - ground[:, :, None]).astype(np.float32)
    inside = (hag >= 0) & (np.arange(nz)[None, None, :] <= top[:, :, None])  # ground to canopy top

    keep = inside & (classes != EMPTY)
    ix, iy, iz = np.nonzero(keep)
    total = nhit[keep].astype(np.float64) + nmiss[keep] + nocc[keep]
    points = {"x": grid.origin[0] + (ix + 0.5) * voxel, "y": grid.origin[1] + (iy + 0.5) * voxel,
              "z": centres[iz], "occlusion": classes[keep], "nhit": nhit[keep].astype(np.uint32),
              "nmiss": nmiss[keep].astype(np.uint32), "nocc": nocc[keep].astype(np.uint32),
              "occluded_fraction": np.where(total > 0, nocc[keep] / np.maximum(total, 1), np.nan).astype(np.float32),
              "height_above_ground": hag[keep].astype(np.float32)}

    n_bands = int(np.floor(float(hag[inside].max()) / BAND_M)) + 1 if inside.any() else 1
    shares = np.full((n_bands + 1, ny, nx), np.nan, np.float32)
    ox, oy, oz = np.nonzero(inside & (classes != UNOBSERVED))
    key = (np.floor(hag[ox, oy, oz] / BAND_M).astype(np.int64) * nx + ox) * ny + oy
    seen = np.bincount(key, minlength=n_bands * nx * ny).reshape(n_bands, nx, ny)
    occ = np.bincount(key, weights=classes[ox, oy, oz] == OCCLUDED, minlength=n_bands * nx * ny)
    occ = occ.reshape(n_bands, nx, ny)
    del ox, oy, oz, key
    for b in range(n_bands):
        shares[b] = np.flipud(np.where(seen[b] > 0, occ[b] / np.maximum(seen[b], 1), np.nan).T)
    column = inside.sum(axis=2)
    unseen = (inside & (classes == UNOBSERVED)).sum(axis=2)
    shares[n_bands] = np.flipud(np.where(column > 0, unseen / np.maximum(column, 1), np.nan).T)
    names = [f"occluded_{b:g}-{b + BAND_M:g}m" for b in np.arange(n_bands) * BAND_M] + ["unobserved_column"]

    crs_wkt = _crs_wkt(paths)
    scratch = scratch or copc.parent
    work = scratch / "occlusion_work"
    work.mkdir(parents=True, exist_ok=True)
    try:
        counts = {name: int((classes[inside] == code).sum())
                  for name, code in (("filled", FILLED), ("empty", EMPTY), ("occluded", OCCLUDED),
                                     ("unobserved", UNOBSERVED))}
        seen_total = counts["filled"] + counts["empty"] + counts["occluded"]
        section = {
            "method": "OccPy voxel traversal (occpy_ls 0.3.0 raytr core)",
            "trajectory": source,
            "voxel_m": voxel,
            "grid": {"origin": [round(float(v), 3) for v in grid.origin], "shape": list(grid.shape),
                     "tiles": len(tiles)},
            "ground": ground_source,
            "ground_spikes_replaced": spikes,
            "pulses_outside_path": counters["outside_path"],
            "traversed_pulses": counters["traversed_pulses"],
            "missing_returns": counters["missing_returns"],
            "voxels_ground_to_canopy": counts,
            "occluded_share": round(counts["occluded"] / seen_total, 4) if seen_total else None,
            "cog_bands": names,
            "copc_points": int(len(points["x"])),
        }
        if len(points["x"]):
            write_copc(copc, points, grid, crs_wkt, work)
        else:
            copc.write_bytes(b"")
        write_cog(cog, shares, names, grid, crs_wkt,
                  {"OCCLUSION_BAND_M": str(BAND_M), "OCCLUSION_GROUND": ground_source,
                   "OCCLUSION_TRAJECTORY": source, "OCCLUSION_VALUE": "occluded / (filled + empty + occluded)"})
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return section
