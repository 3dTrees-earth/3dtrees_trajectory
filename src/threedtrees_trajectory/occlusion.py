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
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

import laspy
import numpy as np

from laspy.copc import Bounds

from .cloud import CHUNK, TimeSortedCopc, is_copc, point_chunks, pulse_chunks, sorted_copc
from .pulses import Pulses, _crs, _sample, _time_system
from .resources import RunResources

FILLED, EMPTY, OCCLUDED, UNOBSERVED = 1, 2, 3, 4  # OccPy's classes
# Native grids, list conversion, NumPy counts and output windows all overlap
# with the active chunk. Bound the estimated working memory before allocating.
GIB = 2**30
DEFAULT_BUDGET = 16 * GIB
HARD_BUDGET = 64 * GIB
GRID_BYTES_PER_VOXEL = 72
GRID_BYTES_PER_COLUMN = 320
WORKING_BYTES = 2 * GIB
MAX_FRAGMENT_ECHOES = 5_000_000
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
    lo = np.floor(lo / voxel) * voxel
    # Include returns at the maximum, even when it falls on a voxel face.
    shape = tuple(math.floor(float(n)) + 1 for n in np.round((hi - lo) / voxel, 9))
    return Grid(lo, shape, voxel), crs.to_wkt() if crs else None


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


def _xy_mask(chunk, bounds, time, keep):
    """Include all returns of a nearby pulse, including those just outside the buffer."""
    x, y = np.asarray(chunk.x), np.asarray(chunk.y)
    lo, hi = bounds.mins, bounds.maxs
    nearby = ((lo[0] <= x) & (x <= hi[0]) & (lo[1] <= y) & (y <= hi[1]))
    if not nearby.any():
        return nearby
    starts = np.flatnonzero(np.r_[True, time[1:] != time[:-1]])
    return keep & np.repeat(np.logical_or.reduceat(nearby, starts), np.diff(np.r_[starts, len(time)]))


def _trace_multi(tracer, coords, sensor, time, rn, nr, source, channel, ordered=False) -> np.ndarray:
    """Trace complete shots. OccPy retains incomplete shots after clearPulseDataset()."""
    same = (time[1:] == time[:-1]) & (source[1:] == source[:-1]) & (channel[1:] == channel[:-1])
    if not ordered:
        # On original strips, a second first return can start another shot
        # at the same GPS time. Spatially ordered data cannot disambiguate it.
        index = np.arange(len(time))
        group_start = np.maximum.accumulate(np.where(np.r_[True, ~same], index, 0))
        previous_first = np.maximum.accumulate(np.where(rn == 1, index, -1))
        same &= ~((rn[1:] == 1) & (previous_first[:-1] >= group_start[1:]))
    starts = np.flatnonzero(np.r_[True, ~same])
    counts = np.diff(np.r_[starts, len(time)])
    target = nr[starts]
    numbers = np.bitwise_or.reduceat(np.left_shift(np.uint32(1), rn), starts)
    complete = (counts == target) & (numbers == np.left_shift(np.uint32(1), target + 1) - 2)
    keep = np.repeat(complete, counts)
    if keep.any():
        shot = np.repeat(np.arange(len(starts), dtype=np.float64), counts)[keep]
        tracer.addPointData(*(axis[keep] for axis in (*coords, *sensor)), shot, rn[keep], nr[keep])
        tracer.doRaytracing()
        tracer.clearPulseDataset()
    return ~keep


def _recover_fragments(tracer, fragments, time, xyz) -> int:
    values = [np.concatenate(parts) for parts in zip(*fragments)]
    x, y, z, t, rn, nr, source, channel = values
    order = np.lexsort((rn, channel, source, t))
    coords = [axis[order] for axis in (x, y, z)]
    t = t[order]
    unpaired = _trace_multi(tracer, coords, _positions(t, time, xyz), t, rn[order], nr[order],
                            source[order], channel[order], ordered=True)
    return int(unpaired.sum())


def _trace_file(path: Path, tracer, grid: Grid, time: np.ndarray, xyz: np.ndarray, counters: dict,
                traced_records: int, prepared=None, bounds=None) -> int:
    chunks = prepared.chunks() if prepared is not None else pulse_chunks(point_chunks(path, CHUNK), CHUNK)
    fragments = []
    num_fragments = unpaired = 0
    for chunk in chunks:
        t = chunk.gps_time
        if not np.isfinite(t).all():
            raise ValueError(f"{path.name}: GPS times must be finite.")
        keep = _time_mask(t, time)
        counters["echoes_read"] += len(t)
        counters["echoes_outside_path"] += int((~keep).sum())
        if bounds is not None and prepared is None:
            keep = _xy_mask(chunk, bounds, t, keep)
        if not keep.any():
            continue
        t = np.ascontiguousarray(t[keep])
        traced_records += len(t)
        if traced_records > np.iinfo(np.int32).max:
            raise ValueError("More than 2.1 billion echoes in one traversal. Reduce --grid-memory-gib "
                             "to use smaller tiles or split the acquisition.")
        # OccPy normalises by the box size before flooring. Nudge exact faces
        # inward by a negligible fraction of a voxel so tiling and z cropping
        # use a consistent half-open lattice at representable coordinates.
        coords = tuple(np.asarray(axis)[keep] + grid.voxel * 1e-7 for axis in (chunk.x, chunk.y, chunk.z))
        sensor = _positions(t, time, xyz)
        rn, nr = np.asarray(chunk.return_number)[keep], np.asarray(chunk.number_of_returns)[keep]
        single = (rn <= 1) & (nr <= 1)
        if single.any():
            tracer.doRaytracing_singleReturnPulses(*(axis[single] for axis in (*coords, *sensor)), t[single])
        if single.all():
            continue
        multi = ~single
        source = chunk.point_source_id[keep][multi]
        channel = np.asarray(chunk.scanner_channel)[keep][multi] if chunk.point_format.id >= 6 else np.zeros(len(source), np.uint8)
        multi_coords = [axis[multi] for axis in coords]
        shot = (t[multi], rn[multi], nr[multi], source, channel)
        pending = _trace_multi(tracer, multi_coords, [axis[multi] for axis in sensor], *shot,
                               ordered=prepared is not None)
        count = int(pending.sum())
        if prepared is not None:
            unpaired += count
        elif count:
            num_fragments += count
            if num_fragments > MAX_FRAGMENT_ECHOES:
                raise ValueError(f"{path.name}: too many separated multi-return echoes for the LAZ stream. "
                                 "Provide acquisition-order strips or a merged COPC.")
            fragments.append(tuple(axis[pending] for axis in (*multi_coords, *shot)))
    if fragments:
        unpaired = _recover_fragments(tracer, fragments, time, xyz)
        counters["recovered_echoes"] = counters.get("recovered_echoes", 0) + num_fragments - unpaired
    counters["unpaired_echoes"] = counters.get("unpaired_echoes", 0) + unpaired
    return traced_records


def _trace_grid(paths: list[Path], grid: Grid, time: np.ndarray, xyz: np.ndarray, counters: dict,
                prepared=None, bounds=None) -> tuple:
    from raytr import PyRaytracer

    tracer = PyRaytracer()
    hi = grid.origin + np.array(grid.shape) * grid.voxel
    tracer.defineGrid(grid.origin, hi, *grid.shape, grid.voxel)
    traced_records = 0
    for path in paths:
        traced_records = _trace_file(path, tracer, grid, time, xyz, counters, traced_records, prepared, bounds)
    counters["traversed_pulses"] += int(tracer.get_num_traversed_pulses())
    counters["missing_returns"] += int(tracer.get_num_missing_returns())
    counts = tuple(np.asarray(values(), np.int32) for values in (tracer.getNhit, tracer.getNmiss, tracer.getNocc))
    if any(values.min(initial=0) < 0 for values in counts):
        raise ValueError("Too many pulses traverse the same voxel for the native 32-bit count.")
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


def _stage(resources, name, **kwargs):
    return resources.stage(name, **kwargs) if resources is not None else nullcontext({})


def _extent_pass(paths, grid, time, xyz, dtm, crs_wkt, include_sample, stats):
    nx, ny, nz = grid.shape
    if nx * ny * 96 > HARD_BUDGET - WORKING_BYTES or nz > np.iinfo(np.int32).max:
        raise ValueError(f"The occlusion extent {grid.shape} exceeds the 64 GiB surface budget; "
                         "increase --voxel or restrict the acquisition.")
    low = np.full(nx * ny, nz, np.int32)
    high = np.full(nx * ny, -1, np.int32)
    sample_parts = []
    for path in paths:
        for chunk in point_chunks(path, CHUNK):
            t = chunk.gps_time
            keep = _time_mask(t, time)
            if not keep.any():
                continue
            coord = np.column_stack([np.asarray(axis)[keep] for axis in (chunk.x, chunk.y, chunk.z)])
            cells = np.floor((coord - grid.origin) / grid.voxel + 1e-8).astype(np.int64)
            for axis, maximum in enumerate(grid.shape):
                np.clip(cells[:, axis], 0, maximum - 1, out=cells[:, axis])
            flat = cells[:, 0] * ny + cells[:, 1]
            np.minimum.at(low, flat, cells[:, 2])
            np.maximum.at(high, flat, cells[:, 2])
            if include_sample:
                pulse_t = t[keep]
                samples = (np.asarray(chunk.return_number)[keep] <= 1) & _sample(pulse_t, 100)
                if samples.any():
                    d = np.column_stack(_positions(pulse_t[samples], time, xyz)) - coord[samples]
                    sample_parts.append((flat[samples], (np.hypot(d[:, 0], d[:, 1]) /
                                                        np.maximum(np.abs(d[:, 2]), 1e-9)).astype(np.float32)))
        stats["bytes_read"] = stats.get("bytes_read", 0) + path.stat().st_size
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


def _grid_bytes(shape, whole_columns=0):
    nx, ny, _ = shape
    return WORKING_BYTES + GRID_BYTES_PER_VOXEL * math.prod(shape) + GRID_BYTES_PER_COLUMN * nx * ny + 32 * whole_columns


def _tile_layout(grid, budget, buffer_m):
    nx, ny, nz = grid.shape
    if _grid_bytes(grid.shape, nx * ny) <= budget:
        return [(slice(0, nx), slice(0, ny))], 0
    halo = math.ceil(buffer_m / grid.voxel)
    # Use the full length if one XY axis already fits across the buffered tile.
    def footprint(edge):
        return min(nx, edge + 2 * halo), min(ny, edge + 2 * halo), nz

    # Use the hard allowance for the smallest tile when the requested buffer
    # and the global surfaces alone exceed the soft target.
    minimum = _grid_bytes(footprint(1), nx * ny)
    if minimum > HARD_BUDGET:
        raise ValueError(f"Even a one-column occlusion tile with the {buffer_m:.1f} m buffer "
                         "would exceed 64 GiB; increase --voxel or restrict the acquisition.")
    budget = max(budget, minimum)
    low, high = 1, max(nx, ny)
    while low < high:
        middle = (low + high + 1) // 2
        if _grid_bytes(footprint(middle), nx * ny) <= budget:
            low = middle
        else:
            high = middle - 1
    if -(-nx // low) * -(-ny // low) > 10_000:
        raise ValueError("The occlusion grid requires more than 10,000 tiles; increase --voxel or the grid budget.")
    return [(slice(x, min(nx, x + low)), slice(y, min(ny, y + low)))
            for x in range(0, nx, low) for y in range(0, ny, low)], halo


def _traversal_tiles(paths, grid, time, xyz, layout, halo, prepared, work, resources, counters):
    tiled = len(layout) > 1
    for number, core in enumerate(layout):
        x, y = core
        start = np.maximum([x.start - halo, y.start - halo], 0)
        end = np.minimum([x.stop + halo, y.stop + halo], grid.shape[:2])
        origin = grid.origin.copy()
        origin[:2] += start * grid.voxel
        tile = Grid(origin, (*map(int, end - start), grid.shape[2]), grid.voxel)
        bounds = Bounds(origin[:2], grid.origin[:2] + end * grid.voxel) if tiled else None
        with _stage(resources, "traversal") as stage:
            stage.update(tile=number + 1, grid_voxels=math.prod(tile.shape),
                         estimated_grid_bytes=_grid_bytes(tile.shape, math.prod(grid.shape[:2])))
            if bounds is not None and prepared is not None:
                with TemporaryDirectory(dir=work, prefix="tile-sort-") as folder:
                    selected = TimeSortedCopc(paths[0], Path(folder), CHUNK, bounds)
                    stage["bytes_written"] = selected.bytes_written
                    stage["points_read"] = selected.points_read
                    counts = _trace_grid(paths, tile, time, xyz, counters, prepared=selected)
            else:
                stage["bytes_read"] = (prepared.bytes_written if prepared is not None
                                       else sum(p.stat().st_size for p in paths))
                counts = _trace_grid(paths, tile, time, xyz, counters, prepared=prepared, bounds=bounds)
        keep = (slice(x.start - start[0], x.stop - start[0]), slice(y.start - start[1], y.stop - start[1]))
        yield core, tuple(a[keep] for a in counts)
        del counts


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


def _write_products(grid, tiles, top, ground, copc, cog, crs_wkt, tags, work, resources):
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
    flat = work / "voxels.las"
    raster = work / "voxel_shares.tif"
    classes_total = np.zeros(5, np.int64)
    seen_total, hidden_total = np.zeros(n_bands, np.int64), np.zeros(n_bands, np.int64)
    with rasterio.Env(GDAL_CACHEMAX=128 * 1024**2, CPL_TMPDIR=str(work)):
        with laspy.open(flat, "w", header=_copc_header(grid, crs_wkt)) as writer, rasterio.open(raster, "w", **profile) as out:
            out.descriptions = tuple(names)
            out.update_tags(OCCLUSION_BAND_M=str(BAND_M), **tags, OCCLUSION_VALUE="occluded / (filled + empty + occluded)")
            for number, (core, counts) in enumerate(tiles):
                with _stage(resources, "classify + ground") as stage:
                    stage["tile"] = number + 1
                    offset = np.array([part.start for part in core])
                    local_grid = Grid(grid.origin, counts[0].shape, grid.voxel)
                    for xy in _xy_windows(local_grid):
                        local_counts = tuple(count[xy] for count in counts)
                        xs, ys = tuple(slice(part.start + start, part.stop + start) for part, start in zip(xy, offset))
                        hits = local_counts[0] > 0
                        local_top = np.where(hits.any(axis=2), nz - 1 - hits[:, :, ::-1].argmax(axis=2), -1)
                        heights = (centres - ground[xs, ys, None]).astype(np.float32)
                        threshold = np.minimum(top[xs, ys], local_top)
                        inside = (heights >= 0) & (layers <= threshold[:, :, None])
                        if not inside.any():
                            continue
                        classes = classify(*local_counts)
                        classes_total += np.bincount(classes[inside], minlength=5)
                        _write_voxels(writer, grid, (xs, ys), centres, heights, classes, inside, local_counts)
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
            copy_raster(raster, cog, driver="COG", compress="DEFLATE", predictor=3, blocksize=512, BIGTIFF="IF_SAFER")
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
        quality: str | None = None, delivered_check: dict | None = None, prepared: TimeSortedCopc | None = None,
        budget: int = DEFAULT_BUDGET, resources: RunResources | None = None) -> dict:
    """Trace, classify and write both products; returns the report section."""
    validate_trajectory(time, xyz)
    scratch = scratch or copc.parent
    copc.parent.mkdir(parents=True, exist_ok=True)
    cog.parent.mkdir(parents=True, exist_ok=True)
    counters = {"echoes_read": 0, "echoes_outside_path": 0, "traversed_pulses": 0, "missing_returns": 0,
                "unpaired_echoes": 0, "recovered_echoes": 0}
    budget = min(budget, HARD_BUDGET)
    with TemporaryDirectory(prefix=".occlusion-", dir=scratch) as folder, \
         TemporaryDirectory(prefix=".occlusion-", dir=cog.parent) as cog_folder, \
         (sorted_copc(paths, Path.cwd(), CHUNK, resources) if prepared is None else nullcontext(prepared)) as stream:
        work = Path(folder)
        with _stage(resources, "grid extent/crop pass") as stage:
            full_grid, crs_wkt = _grid(paths, voxel)
            sample = _grid_bytes(full_grid.shape, math.prod(full_grid.shape[:2])) > budget
            grid, top, ground, spikes, percentile = _extent_pass(paths, full_grid, time, xyz, dtm, crs_wkt, sample, stage)
        ground_source = "dtm" if dtm is not None else "lowest_hit"
        buffer = min(100.0, percentile)
        layout, halo = _tile_layout(grid, budget, buffer)
        if stream is not None and len(layout) > 1:
            # Tiles reread and sort their COPC windows. The acquisition sort
            # used to estimate/check the path is no longer needed.
            stream.discard()
        tags = dict(OCCLUSION_GROUND=ground_source, OCCLUSION_TRAJECTORY=source)
        if quality is not None:
            tags["OCCLUSION_QUALITY"] = quality
        if delivered_check is not None:
            flag = delivered_check["threshold_exceeded"]
            tags["OCCLUSION_DELIVERED_MISS_FLAG"] = "unavailable" if flag is None else str(flag).lower()
            if flag is not None:
                tags["OCCLUSION_DELIVERED_MISS_MEDIAN_M"] = str(delivered_check["median_m"])
        staged_copc, staged_cog = work / "voxels.copc.laz", Path(cog_folder) / "occlusion.tif"
        tiles = _traversal_tiles(paths, grid, time, xyz, layout, halo, stream, work, resources, counters)
        products = _write_products(grid, tiles, top, ground, staged_copc, staged_cog, crs_wkt, tags, work, resources)
        os.replace(staged_copc, copc)
        os.replace(staged_cog, cog)
    return {
        "method": "OccPy voxel traversal (occpy_ls 0.3.0 raytr core)",
        "trajectory": source,
        **({"quality": quality} if quality else {}),
        **({"delivered_path_check": delivered_check} if delivered_check is not None else {}),
        "voxel_m": voxel,
        "grid": {"origin": [round(float(v), 3) for v in grid.origin], "shape": list(grid.shape),
                 "uncropped_shape": list(full_grid.shape), "voxels": math.prod(grid.shape),
                 "uncropped_voxels": math.prod(full_grid.shape)},
        "tiling": {"tiles": len(layout), "buffer_m": halo * grid.voxel, "sample_p99_9_m": percentile if sample else None,
                   "budget_bytes": budget},
        "ground": ground_source,
        "ground_spikes_replaced": spikes,
        **counters,
        **products,
    }
