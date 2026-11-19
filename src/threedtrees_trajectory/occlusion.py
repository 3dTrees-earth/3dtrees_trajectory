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
import json
import math
import multiprocessing
import os
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

import laspy
import numpy as np

from .cloud import CHUNK, is_copc, point_chunks
from .pulses import MAX_STRAY, Pulses, _crs, _sample, _time_system
from .resources import RunResources, available_memory, decoding_threads
from .tile_jobs import atomic_json, discard_finished_spools, execute, identify, lock_job, save_tile, tile_files
from .tile_spool import ECHO, START, TileSpools, pack_echoes, raw_batches, return_fields, tile_batches

FILLED, EMPTY, OCCLUDED, UNOBSERVED = 1, 2, 3, 4  # OccPy's classes
# Native grids, list conversion, NumPy counts and output windows all overlap
# with the active chunk. Bound the estimated working memory before allocating.
GIB = 2**30
DEFAULT_BUDGET = 16 * GIB
HARD_BUDGET = 64 * GIB
GRID_BYTES_PER_VOXEL = 72
GRID_BYTES_PER_COLUMN = 320
WORKING_BYTES = 2 * GIB
MAX_FRAGMENT_ECHOES = MAX_STRAY
FRAGMENT_BYTES_PER_ECHO = 112
NATIVE_ECHO_LIMIT = int(np.iinfo(np.int32).max)
ECHO_BIN = 16
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
    points: int = 0


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


def _stage(resources, name, **kwargs):
    return resources.stage(name, **kwargs) if resources is not None else nullcontext({})


def _extent_points(paths, time, stats, spool, spatial, transforms, source):
    for file, path in enumerate(paths):
        if source is not None:
            raw = source.path(file)
            if not raw.is_file():
                continue
            scales, offsets = transforms[file]
            for records in raw_batches(raw):
                t = records["time"]
                keep = _time_mask(t, time)
                stats["echoes_read"] = stats.get("echoes_read", 0) + len(t)
                stats["echoes_outside_path"] = stats.get("echoes_outside_path", 0) + len(t) - int(keep.sum())
                if keep.any():
                    coord = np.column_stack([records[axis][keep] * scale + offset
                                             for axis, scale, offset in zip(("x", "y", "z"), scales, offsets)])
                    yield coord, t[keep], (records["returns"][keep] & 15)
            stats["bytes_read"] = stats.get("bytes_read", 0) + raw.stat().st_size
            continue
        with (spool(file) if spool is not None else nullcontext()) as stream:
            for chunk in point_chunks(path, CHUNK, spatial=spatial):
                t = chunk.gps_time
                keep = _time_mask(t, time)
                if stream is not None:
                    stats["echoes_read"] = stats.get("echoes_read", 0) + len(t)
                    stats["echoes_outside_path"] = stats.get("echoes_outside_path", 0) + len(t) - int(keep.sum())
                    if not np.isfinite([t.min(), t.max()]).all():
                        raise ValueError(f"{path.name}: GPS times must be finite.")
                    records = pack_echoes(chunk, keep)
                    records.tofile(stream)
                    stats["bytes_written"] = stats.get("bytes_written", 0) + records.nbytes
                if keep.any():
                    coord = np.column_stack([np.asarray(axis)[keep] for axis in (chunk.x, chunk.y, chunk.z)])
                    yield coord, t[keep], np.asarray(chunk.return_number)[keep]
        stats["bytes_read"] = stats.get("bytes_read", 0) + path.stat().st_size


def _extent_pass(paths, grid, time, xyz, dtm, crs_wkt, include_sample, stats, echo_bins=None, spool=None, spatial=False,
                 transforms=None, source=None):
    nx, ny, nz = grid.shape
    if nx * ny * 96 > HARD_BUDGET - WORKING_BYTES or nz > np.iinfo(np.int32).max:
        raise ValueError(f"The occlusion extent {grid.shape} exceeds the 64 GiB surface budget; "
                         "increase --voxel or restrict the acquisition.")
    low = np.full(nx * ny, nz, np.int32)
    high = np.full(nx * ny, -1, np.int32)
    sample_parts = []
    for coord, pulse_t, rn in _extent_points(paths, time, stats, spool, spatial, transforms, source):
        cells = np.floor((coord - grid.origin) / grid.voxel + 1e-8).astype(np.int64)
        for axis, maximum in enumerate(grid.shape):
            np.clip(cells[:, axis], 0, maximum - 1, out=cells[:, axis])
        flat = cells[:, 0] * ny + cells[:, 1]
        np.minimum.at(low, flat, cells[:, 2])
        np.maximum.at(high, flat, cells[:, 2])
        if echo_bins is not None:
            bins = (cells[:, 0] // ECHO_BIN) * echo_bins.shape[1] + cells[:, 1] // ECHO_BIN
            echo_bins += np.bincount(bins, minlength=echo_bins.size).reshape(echo_bins.shape)
        if include_sample:
            samples = (rn <= 1) & _sample(pulse_t, 100)
            if samples.any():
                d = np.column_stack(_positions(pulse_t[samples], time, xyz)) - coord[samples]
                sample_parts.append((flat[samples], (np.hypot(d[:, 0], d[:, 1]) /
                                                    np.maximum(np.abs(d[:, 2]), 1e-9)).astype(np.float32)))
    if (high < 0).all():
        raise ValueError("The point cloud and path have no overlapping GPS times. "
                         "Check GPS week seconds versus adjusted GPS time.")
    low[low == nz] = -1
    low, high = low.reshape(nx, ny), high.reshape(nx, ny)
    if dtm is None:
        ground, spikes = _ground_from_lowest(low, grid)
    else:
        ground, spikes = ground_from_dtm(dtm, grid, crs_wkt), 0
        stats["bytes_read"] += dtm.stat().st_size
    top_z = grid.origin[2] + (high + 0.5) * grid.voxel
    offset = max(0, math.floor((float(ground.min()) - 1 - grid.origin[2]) / grid.voxel))
    offset = min(offset, nz - 1)
    stop = min(nz, max(offset + 1, math.ceil((float(top_z.max()) + 1 - grid.origin[2]) / grid.voxel)))
    origin = grid.origin.copy()
    origin[2] += offset * grid.voxel
    cropped = Grid(origin, (nx, ny, stop - offset), grid.voxel)
    if sample_parts:
        height = np.maximum(top_z - ground, 0).ravel()
        samples = np.concatenate([height[index] * slope for index, slope in sample_parts])
        percentile = float(np.percentile(samples, 99.9))
    else:
        percentile = 0.0
    return cropped, high - offset, ground, spikes, percentile


def _grid_bytes(shape, whole_columns=0, fragment_bytes=0):
    nx, ny, _ = shape
    return WORKING_BYTES + fragment_bytes + GRID_BYTES_PER_VOXEL * math.prod(shape) + GRID_BYTES_PER_COLUMN * nx * ny + 32 * whole_columns


@dataclass(frozen=True)
class Tile:
    core: tuple[slice, slice]
    extent: tuple[slice, slice, slice]
    echo_upper_bound: int | None = None

    @property
    def shape(self) -> tuple[int, int, int]:
        return tuple(part.stop - part.start for part in self.extent)


def _tile_layout(grid, top, ground, budget, buffer_m, fragment_bytes=0, echo_bins=None):
    nx, ny, nz = grid.shape
    core = (slice(0, nx), slice(0, ny))
    total_echoes = int(echo_bins.sum()) if echo_bins is not None else None
    if (_grid_bytes(grid.shape, nx * ny, fragment_bytes) <= budget
            and (total_echoes is None or total_echoes <= NATIVE_ECHO_LIMIT)):
        return [Tile(core, (*core, slice(0, nz)), total_echoes)], 0
    prefix = (np.pad(echo_bins.cumsum(axis=0).cumsum(axis=1), ((1, 0), (1, 0)))
              if echo_bins is not None else None)
    halo = min(math.ceil(buffer_m / grid.voxel), math.floor(100.0 / grid.voxel + 1e-12))
    minimum = _grid_bytes((min(nx, 1 + 2 * halo), min(ny, 1 + 2 * halo), nz), nx * ny, fragment_bytes)
    budget = max(budget, min(minimum, HARD_BUDGET))
    layout, pending = [], [core]
    while pending:
        core = pending.pop()
        if (top[core] < 0).all():
            continue
        buffered = tuple(slice(max(0, part.start - halo), min(size, part.stop + halo))
                         for part, size in zip(core, (nx, ny)))
        low = math.floor((float(ground[buffered].min()) - 1 - grid.origin[2]) / grid.voxel)
        bottom = min(nz - 1, max(0, low))
        high = math.ceil(float(top[buffered].max()) + 0.5 + 1 / grid.voxel)
        upper_bound = None
        if prefix is not None:
            # The spatial query includes points on the upper faces. Cover
            # every intersecting bin to bound the number of input echoes.
            x, y = buffered
            x0, y0 = x.start // ECHO_BIN, y.start // ECHO_BIN
            x1 = min(prefix.shape[0] - 1, x.stop // ECHO_BIN + 1)
            y1 = min(prefix.shape[1] - 1, y.stop // ECHO_BIN + 1)
            upper_bound = int(prefix[x1, y1] - prefix[x0, y1] - prefix[x1, y0] + prefix[x0, y0])
        tile = Tile(core, (*buffered, slice(bottom, min(nz, max(bottom + 1, high)))), upper_bound)
        amount = _grid_bytes(tile.shape, nx * ny, fragment_bytes)
        if amount <= budget and (upper_bound is None or upper_bound <= NATIVE_ECHO_LIMIT):
            layout.append(tile)
        else:
            axis = int(core[1].stop - core[1].start > core[0].stop - core[0].start)
            part = core[axis]
            middle = (part.start + part.stop) // 2
            if middle == part.start:
                if upper_bound is not None and upper_bound > NATIVE_ECHO_LIMIT:
                    raise ValueError("A buffered tile exceeds the native 2.1 billion-echo limit. Split the acquisition.")
                raise ValueError(f"Even a one-column occlusion tile with the {buffer_m:.1f} m buffer "
                                 "would exceed 64 GiB; increase --voxel or restrict the acquisition.")
            for split in (slice(middle, part.stop), slice(part.start, middle)):
                pending.append((split, core[1]) if axis == 0 else (core[0], split))
        if len(layout) + len(pending) > 10_000:
            raise ValueError("The occlusion grid requires more than 10,000 tiles; increase --voxel or the grid budget.")
    return layout, halo


def _read_storage_headers(paths):
    transforms, spatial = [], False
    for path in paths:
        with laspy.open(path, read_evlrs=False) as reader:
            transforms.append((reader.header.scales.tolist(), reader.header.offsets.tolist()))
            spatial |= is_copc(reader.header)
    if spatial and len(paths) != 1:
        raise ValueError("Use one merged COPC or a list of per-file LAS/LAZ inputs.")
    return transforms, spatial


def _route_strip(arguments):
    spool, grid, layout, time, points, file, path, transform, time_bounds = arguments
    routing = TileSpools(spool, grid, layout, time, points, compress=True)
    recovered, unpaired = routing.stage_file(file, path, *transform, False, time_bounds)
    return recovered, unpaired, routing.bytes_written, routing.echoes


def _route_file_spools(spool, grid, layout, time, points, paths, transforms, time_bounds, spatial, slots, stage):
    compressed = points > NATIVE_ECHO_LIMIT
    route_slots = min(slots, max(1, available_memory() // (2 * GIB))) if compressed else 1
    workers = min(route_slots, len(paths)) if not spatial else 1
    recovered = unpaired = written = 0
    if workers > 1:
        calls = [(spool, grid, layout, time, points, index, path, transform, time_bounds)
                 for index, (path, transform) in enumerate(zip(paths, transforms))]
        tile_echoes = np.zeros(len(layout), np.int64)
        with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as pool:
            for joined, missing, size, echoes in pool.map(_route_strip, calls):
                recovered += joined
                unpaired += missing
                written += size
                tile_echoes += echoes
    else:
        routing = TileSpools(spool, grid, layout, time, points, compress=compressed)
        for index, (path, transform) in enumerate(zip(paths, transforms)):
            joined, missing = routing.stage_file(index, path, *transform, spatial, time_bounds,
                                                 hash_slots=route_slots if spatial else 1)
            recovered += joined
            unpaired += missing
        written, tile_echoes = routing.bytes_written, routing.echoes
        workers = max(workers, routing.hash_workers)
    echoes_routed = int(tile_echoes.sum())
    stage.update(bytes_written=written, echoes_routed=echoes_routed, uncompressed_bytes=echoes_routed * ECHO.itemsize,
                 compressed=compressed, routing_workers=workers)
    return recovered, unpaired, tile_echoes


def _prepare_job(work, paths, time, xyz, voxel, dtm, budget, resources, fragment_echoes, source, slots):
    manifest = work / "prepared.json"
    spool = work / "spools"
    surfaces = work / "surfaces.npz"
    if manifest.is_file():
        data = json.loads(manifest.read_text(encoding="utf-8"))
    else:
        if spool.is_dir():
            shutil.rmtree(spool)
        spool.mkdir()
        transforms, spatial = _read_storage_headers(paths)
        full_grid, crs_wkt = _grid(paths, voxel)
        # Preserve the established layout for previously supported strips.
        # Fragment buffers no longer contribute to worker reservations.
        legacy_allowance = 0 if spatial else FRAGMENT_BYTES_PER_ECHO * min(fragment_echoes or 0, MAX_FRAGMENT_ECHOES)
        planning_budget = min(budget, available_memory())
        with _stage(resources, "grid extent/crop pass") as stage:
            if planning_budget <= WORKING_BYTES:
                raise ValueError("The occlusion job lacks memory for one tile's point buffers.")
            echo_bins = None
            if full_grid.points > NATIVE_ECHO_LIMIT:
                dimensions = tuple(-(-length // ECHO_BIN) for length in full_grid.shape[:2])
                echo_bins = np.zeros(dimensions, np.int64)
            sample = (_grid_bytes(full_grid.shape, math.prod(full_grid.shape[:2]), legacy_allowance)
                      > planning_budget or echo_bins is not None)
            raw_paths = [(source.path(index) if source is not None else spool / f"source-{index:05d}.bin")
                         for index in range(len(paths))]
            grid, top, ground, spikes, percentile = _extent_pass(
                paths, full_grid, time, xyz, dtm, crs_wkt, sample, stage,
                echo_bins=echo_bins, spool=lambda file: raw_paths[file].open("xb"), spatial=spatial,
                transforms=transforms, source=source)
            if stage["echoes_read"] != full_grid.points:
                raise ValueError("The occlusion scan did not decode the point count declared by the headers.")
            layout, halo = _tile_layout(grid, top, ground, planning_budget, min(100.0, percentile),
                                        legacy_allowance, echo_bins)
            stage.update(grid_shape=list(grid.shape), grid_voxels=math.prod(grid.shape),
                         planned_tiles=len(layout), buffer_m=halo * voxel)
            counters = {key: stage.get(key, 0) for key in ("echoes_read", "echoes_outside_path")}
        with _stage(resources, "route tile spools", children=True) as stage:
            time_bounds = (time[0] - TIME_MARGIN, time[-1] + TIME_MARGIN) if source is not None else None
            try:
                recovered, unpaired, tile_echoes = _route_file_spools(
                    spool, grid, layout, time, full_grid.points, raw_paths, transforms, time_bounds, spatial, slots, stage)
            except BrokenProcessPool as error:
                raise ValueError("An occlusion spool worker stopped during source routing.") from error
            counters.update(recovered_echoes=recovered, unpaired_echoes=unpaired)
        if any(tile_echoes > NATIVE_ECHO_LIMIT):
            raise ValueError("A buffered tile including complete pulses exceeds the native 2.1 billion-echo limit. "
                             "Lower --grid-memory-gib or split the acquisition.")
        data = dict(origin=grid.origin.tolist(), shape=grid.shape, full_shape=full_grid.shape, crs=crs_wkt,
                    transforms=transforms, spikes=spikes, percentile=percentile if sample else None, halo=halo,
                    counters=counters, tile_echoes=tile_echoes.tolist(),
                    tiles=[dict(core=[(s.start, s.stop) for s in tile.core],
                                extent=[(s.start, s.stop) for s in tile.extent],
                                echo_upper_bound=tile.echo_upper_bound) for tile in layout])
        with surfaces.open("wb") as stream:
            np.savez(stream, top=top, ground=ground)
            stream.flush()
            os.fsync(stream.fileno())
        atomic_json(manifest, data)
    with np.load(surfaces, allow_pickle=False) as arrays:
        top, ground = arrays["top"], arrays["ground"]
    grid = Grid(np.array(data["origin"]), tuple(data["shape"]), voxel)
    layout = [Tile(tuple(slice(*s) for s in tile["core"]), tuple(slice(*s) for s in tile["extent"]),
                   tile["echo_upper_bound"]) for tile in data["tiles"]]
    return data, grid, layout, spool, top, ground


def _init_traversal(work, spool, time, xyz, transforms):
    global _TILE_CONTEXT
    _TILE_CONTEXT = work, spool, time, xyz, np.asarray(transforms)


def _trace_tile(index, tile, core, expected_echoes):
    from time import perf_counter, process_time

    from raytr import PyRaytracer

    from .resources import _rss

    work, spool, time, xyz, transforms = _TILE_CONTEXT
    start_wall, start_cpu = perf_counter(), process_time()
    tracer = PyRaytracer()
    tracer.defineGrid(tile.origin, tile.origin + np.array(tile.shape) * tile.voxel, *tile.shape, tile.voxel)
    points = 0
    with TemporaryDirectory(dir=work, prefix=f"worker-{index:05d}-") as folder:
        for records in tile_batches(spool / f"tile-{index:05d}", Path(folder)):
            t = np.ascontiguousarray(records["time"])
            files = records["file"]
            if np.all(files == files[0]):
                scales, offsets = transforms[files[0]]
                coords = tuple(records[axis] * scales[k] + offsets[k] + tile.voxel * 1e-7
                               for k, axis in enumerate(("x", "y", "z")))
            else:
                coords = tuple(records[axis] * transforms[files, 0, k] + transforms[files, 1, k] + tile.voxel * 1e-7
                               for k, axis in enumerate(("x", "y", "z")))
            sensor = _positions(t, time, xyz)
            rn, nr = return_fields(records)
            single = (rn <= 1) & (nr <= 1)
            if single.any():
                tracer.doRaytracing_singleReturnPulses(*(axis[single] for axis in (*coords, *sensor)), t[single])
            multi = ~single
            if multi.any():
                shot = np.cumsum((records["flags"] & START) != 0, dtype=np.float64)
                tracer.addPointData(*(axis[multi] for axis in (*coords, *sensor)), shot[multi], rn[multi], nr[multi])
                tracer.doRaytracing()
                tracer.clearPulseDataset()
            points += len(t)
    if points != expected_echoes:
        raise ValueError(f"Occlusion tile {index + 1} spool has {points:,} of {expected_echoes:,} expected echoes.")
    traversed = int(tracer.get_num_traversed_pulses())
    missing = int(tracer.get_num_missing_returns())
    counts = tuple(np.asarray(values(), np.int32)[core] for values in (tracer.getNhit, tracer.getNmiss, tracer.getNocc))
    del tracer
    if any(values.min(initial=0) < 0 for values in counts):
        raise ValueError("Too many pulses traverse the same voxel for the native 32-bit count.")
    data = dict(tile=index + 1, traversed_pulses=traversed, missing_returns=missing, points_read=points,
                grid_shape=tile.shape, grid_voxels=math.prod(tile.shape), wall_seconds=perf_counter() - start_wall,
                cpu_seconds=process_time() - start_cpu, peak_rss_bytes=_rss(), bytes_written=sum(a.nbytes for a in counts))
    save_tile(work, index, counts, data)
    return data


def _tile_spec(grid, index, plan, expected_echoes):
    start = np.array([part.start for part in plan.extent])
    tile = Grid(grid.origin + start * grid.voxel, plan.shape, grid.voxel)
    core = tuple(slice(part.start - offset, part.stop - offset)
                 for part, offset in zip(plan.core, start[:2]))
    shape = tuple(part.stop - part.start for part in (*core, plan.extent[2]))
    return index, _grid_bytes(tile.shape), shape, (tile, core, expected_echoes)


def _completed_tiles(work, layout):
    for index, plan in enumerate(layout):
        counts, _ = tile_files(work, index)
        yield plan.core, plan.extent[2], tuple(np.load(counts, mmap_mode="r", allow_pickle=False))


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


def _write_products(grid, tiles, top, ground, copc, cog, crs_wkt, tags, work, resources, slots):
    import rasterio
    from rasterio.shutil import copy as copy_raster
    from rasterio.transform import from_origin
    from rasterio.windows import Window

    nx, ny, nz = grid.shape
    layers = np.arange(nz)
    centres = grid.origin[2] + (layers + 0.5) * grid.voxel
    max_height = (centres[np.clip(top, 0, nz - 1)] - ground)[top >= 0].astype(np.float32).max(initial=0)
    n_bands = math.floor(float(max_height) / BAND_M) + 1
    names = [f"occluded_{b:g}-{b + BAND_M:g}m" for b in np.arange(n_bands) * BAND_M] + ["unobserved_column"]
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
    cache_mib = 512 if wide_raster else 128
    cog_threads = min(slots, 4) if wide_raster else 1
    flat = work / "voxels.las"
    raster = work / "voxel_shares.tif"
    classes_total = np.zeros(5, np.int64)
    seen_total, hidden_total = np.zeros(n_bands, np.int64), np.zeros(n_bands, np.int64)
    with rasterio.Env(GDAL_CACHEMAX=cache_mib * 1024**2, CPL_TMPDIR=str(work)):
        with laspy.open(flat, "w", header=_copc_header(grid, crs_wkt)) as writer, rasterio.open(raster, "w", **profile) as out:
            out.descriptions = tuple(names)
            out.update_tags(OCCLUSION_BAND_M=str(BAND_M), **tags, OCCLUSION_VALUE="occluded / (filled + empty + occluded)")
            for number, (core, z, counts) in enumerate(tiles):
                with _stage(resources, "classify + ground") as stage:
                    stage["tile"] = number + 1
                    offset = np.array([part.start for part in core])
                    local_grid = Grid(grid.origin, counts[0].shape, grid.voxel)
                    for xy in _xy_windows(local_grid):
                        local_counts = tuple(count[xy] for count in counts)
                        xs, ys = tuple(slice(part.start + start, part.stop + start) for part, start in zip(xy, offset))
                        hits = local_counts[0] > 0
                        local_top = np.where(hits.any(axis=2), z.stop - 1 - hits[:, :, ::-1].argmax(axis=2), -1)
                        heights = (centres[z] - ground[xs, ys, None]).astype(np.float32)
                        threshold = np.minimum(top[xs, ys], local_top)
                        inside = (heights >= 0) & (layers[z] <= threshold[:, :, None])
                        if not inside.any():
                            continue
                        classes = classify(*local_counts)
                        classes_total += np.bincount(classes[inside], minlength=5)
                        _write_voxels(writer, grid, (xs, ys), centres[z], heights, classes, inside, local_counts)
                        shares, seen, hidden = _band_shares(classes, heights, inside, n_bands)
                        out.write(shares, window=Window(xs.start, ny - ys.stop, xs.stop - xs.start, ys.stop - ys.start))
                        seen_total += seen
                        hidden_total += hidden
                    del counts, local_counts
            size = writer.header.point_count
        with _stage(resources, "COPC write", children=True) as stage:
            stage["bytes_read"] = flat.stat().st_size
            _index_copc(flat, copc, size, work)
            stage["bytes_written"] = copc.stat().st_size
        with _stage(resources, "COG write") as stage:
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
        "cog_bands": names,
        "occluded_share_by_height_band": [int(hidden) / int(seen) if seen else None
                                          for seen, hidden in zip(seen_total, hidden_total)],
        "copc_points": int(size),
    }


def run(paths: list[Path], time: np.ndarray, xyz: np.ndarray, copc: Path, cog: Path, voxel: float = 0.5,
        dtm: Path | None = None, scratch: Path | None = None, source: str = "estimated", *,
        quality: str | None = None, delivered_check: dict | None = None,
        budget: int = DEFAULT_BUDGET, resources: RunResources | None = None, fragment_echoes: int | None = None,
        threads: int | None = None, source_spool=None) -> dict:
    """Prepare whole-pulse spools, trace resumable tiles, then merge both products."""
    validate_trajectory(time, xyz)
    scratch = scratch or copc.parent
    copc.parent.mkdir(parents=True, exist_ok=True)
    cog.parent.mkdir(parents=True, exist_ok=True)
    budget = min(budget, HARD_BUDGET)
    slots = decoding_threads(threads)
    identity = identify([*paths, *([dtm] if dtm is not None else [])],
                        [str(copc.resolve()), str(cog.resolve()), voxel, budget, source, quality, fragment_echoes,
                         dtm is not None], time, xyz)
    work = scratch / f".occlusion-{identity}"
    with lock_job(work):
        for prefix in ("worker-*", "merge-*"):
            for orphan in work.glob(prefix):
                if orphan.is_dir():
                    shutil.rmtree(orphan)
        data, grid, layout, spool, top, ground = _prepare_job(
            work, paths, time, xyz, voxel, dtm, budget, resources, fragment_echoes, source_spool, slots)
        ground_source = "dtm" if dtm is not None else "lowest_hit"
        tags = dict(OCCLUSION_GROUND=ground_source, OCCLUSION_TRAJECTORY=source)
        if quality is not None:
            tags["OCCLUSION_QUALITY"] = quality
        if delivered_check is not None:
            flag = delivered_check["threshold_exceeded"]
            tags["OCCLUSION_DELIVERED_MISS_FLAG"] = "unavailable" if flag is None else str(flag).lower()
            if flag is not None:
                tags["OCCLUSION_DELIVERED_MISS_MEDIAN_M"] = str(delivered_check["median_m"])
        specs = [_tile_spec(grid, index, plan, echoes)
                 for index, (plan, echoes) in enumerate(zip(layout, data["tile_echoes"]))]
        priority = sorted(specs, key=lambda job: (-data["tile_echoes"][job[0]], -job[1]))
        try:
            with _stage(resources, "traversal", children=True) as stage:
                tiles = execute(work, priority, slots, _trace_tile, _init_traversal,
                                (work, spool, time, xyz, data["transforms"]), stage)
                stage["tiles"] = [tiles[i] | {"estimated_grid_bytes": specs[i][1]} for i in range(len(layout))]
        except BrokenProcessPool as error:
            raise ValueError(f"An occlusion worker stopped. Repeat the command to resume from {work.name}.") from error
        finally:
            discard_finished_spools(spool, work, layout)
        with TemporaryDirectory(dir=work, prefix="merge-") as folder, \
             TemporaryDirectory(prefix=".occlusion-", dir=cog.parent) as cog_folder:
            merge, staged_cog = Path(folder), Path(cog_folder) / "occlusion.tif"
            staged_copc = merge / "voxels.copc.laz"
            products = _write_products(grid, _completed_tiles(work, layout), top, ground, staged_copc, staged_cog,
                                       data["crs"], tags, merge, resources, slots)
            os.replace(staged_copc, copc)
            os.replace(staged_cog, cog)
        shutil.rmtree(work)
    return {
        "method": "OccPy voxel traversal (occpy_ls 0.3.0 raytr core)",
        "trajectory": source,
        **({"quality": quality} if quality else {}),
        **({"delivered_path_check": delivered_check} if delivered_check is not None else {}),
        "voxel_m": voxel,
        "grid": {"origin": [round(float(v), 3) for v in grid.origin], "shape": list(grid.shape),
                 "uncropped_shape": list(data["full_shape"]), "voxels": math.prod(grid.shape),
                 "uncropped_voxels": math.prod(data["full_shape"])},
        "tiling": {"tiles": len(layout), "buffer_m": data["halo"] * grid.voxel,
                   "sample_p99_9_m": data["percentile"], "budget_bytes": budget,
                   "peak_concurrent_workers": stage["peak_concurrent_workers"], "reused_tiles": stage["reused_tiles"]},
        "ground": ground_source,
        "ground_spikes_replaced": data["spikes"],
        **data["counters"],
        **{key: sum(tile[key] for tile in tiles.values()) for key in ("traversed_pulses", "missing_returns")},
        **products,
    }
