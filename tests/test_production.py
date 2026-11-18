import json
import math

import laspy
import numpy as np
import pytest
from laspy.copc import Bounds

from synthetic import flight
from test_occlusion import needs_untwine
from threedtrees_trajectory import occlusion, pulses
from threedtrees_trajectory.cli import main
from threedtrees_trajectory.cloud import TimeSortedCopc


def _delivered(path, time, xyz):
    np.savetxt(path, np.column_stack((time, xyz)), delimiter=",", header="time,x,y,z", comments="")


@needs_untwine
def test_spatial_copc_rebuilds_pulses_across_small_sort_runs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    original = tmp_path / "flight.las"
    copc = tmp_path / "flight.copc.laz"
    flight(original, seconds=6, rate=500)
    with laspy.open(original) as reader:
        occlusion._index_copc(original, copc, reader.header.point_count, tmp_path)

    work = tmp_path / "sort"
    work.mkdir()
    source = TimeSortedCopc(copc, work, 15)
    assert len(source.runs) > 100
    first = [chunk.gps_time.copy() for chunk in source.chunks()]
    assert max(map(len, first)) <= 15
    assert all(p[-1] < q[0] for p, q in zip(first, first[1:]))
    gps = np.concatenate(first)
    assert np.all(gps[1:] >= gps[:-1])
    assert source.points_read == len(gps)

    monkeypatch.setattr(pulses, "CHUNK", 15)
    a, _ = pulses.read_pulses([original], max_multi=1000, max_single=1000)
    b, _ = pulses.read_pulses([copc], max_multi=1000, max_single=1000)
    for name in ("origin", "time", "first", "last", "returns"):
        np.testing.assert_array_equal(getattr(a, name), getattr(b, name))
    for index, cloud in enumerate((original, copc)):
        assert main([str(cloud), "-o", str(tmp_path / f"path-{index}.csv")]) == 0
    def positions(path):
        return [row for row in path.read_bytes().splitlines() if not row.startswith(b"#")]

    assert positions(tmp_path / "path-0.csv") == positions(tmp_path / "path-1.csv")
    assert not list(tmp_path.glob(".copc-sort-*"))


@needs_untwine
def test_copc_tile_stream_reads_only_buffered_bounds(tmp_path):
    original, copc = tmp_path / "flight.las", tmp_path / "flight.copc.laz"
    flight(original, seconds=2, rate=500)
    data = laspy.read(original)
    occlusion._index_copc(original, copc, len(data.points), tmp_path)
    midpoint = np.median(np.asarray(data.x))
    bounds = Bounds(np.array([midpoint, data.header.mins[1]]), data.header.maxs[:2])
    selected = TimeSortedCopc(copc, tmp_path, 101, bounds)
    gps = np.concatenate([p.gps_time for p in selected.chunks()])
    np.testing.assert_array_equal(gps, np.sort(data.gps_time[np.asarray(data.x) >= midpoint]))


def test_laz_traversal_recovers_bounded_separated_echoes(tmp_path, monkeypatch):
    original, separated = tmp_path / "flight.las", tmp_path / "separated.laz"
    time, xyz = flight(original, seconds=1, rate=100)
    data = laspy.read(original)
    order = np.argsort(np.asarray(data.return_number), kind="stable")
    data.points = data.points[order]
    data.write(separated)
    grid, _ = occlusion._grid([original], 1)
    monkeypatch.setattr(occlusion, "CHUNK", 47)

    def trace(path):
        stats = dict(echoes_read=0, echoes_outside_path=0, traversed_pulses=0, missing_returns=0)
        return occlusion._trace_grid([path], grid, time, xyz, stats), stats

    expected, _ = trace(original)
    actual, recovered = trace(separated)
    assert recovered["recovered_echoes"] > 100 and recovered["unpaired_echoes"] == 0
    endpoints, info = pulses.read_pulses([separated], full=True)
    assert len(endpoints) == len(time) and info.unpaired_multi_echoes == 0
    for a, b in zip(expected, actual):
        np.testing.assert_array_equal(a, b)


@needs_untwine
def test_crop_keeps_products_inside_vegetation_band(tmp_path, monkeypatch):
    import rasterio

    cloud = tmp_path / "flight.las"
    time, xyz = flight(cloud, seconds=2, rate=100)
    data = laspy.read(cloud)
    points = data.points.copy()
    data.points = laspy.ScaleAwarePointRecord(np.concatenate([points.array, points.array[:2]]),
                                               points.point_format, data.header.scales, data.header.offsets)
    data.gps_time[-2:] = time[-1] + 10
    data.z[-2:] = [-20, 80]
    data.write(cloud)
    crop = occlusion._extent_pass

    def without_crop(paths, grid, *args, **kwargs):
        clipped, top, ground, spikes, percentile = crop(paths, grid, *args, **kwargs)
        offset = round((clipped.origin[2] - grid.origin[2]) / grid.voxel)
        return grid, np.where(top >= 0, top + offset, top), ground, spikes, percentile

    results = []
    for index in range(2):
        if index:
            monkeypatch.setattr(occlusion, "_extent_pass", without_crop)
        copc, cog = (tmp_path / f"{index}.{suffix}" for suffix in ("copc.laz", "tif"))
        report = occlusion.run([cloud], time, xyz, copc, cog, voxel=1)
        with rasterio.open(cog) as tif:
            points = laspy.read(copc)
            xyz_key = np.round(points.xyz / 0.001).astype(np.int64)
            order = np.lexsort(xyz_key.T[::-1])
            results.append((report, tif.read(), xyz_key[order], points.points.array[order]))
    (cropped, crop_tif, crop_xyz, crop_pts), (full, full_tif, full_xyz, full_pts) = results
    assert cropped["grid"]["voxels"] < full["grid"]["voxels"]
    assert cropped["voxels_ground_to_canopy"] == full["voxels_ground_to_canopy"]
    np.testing.assert_array_equal(crop_tif, full_tif)
    np.testing.assert_array_equal(crop_xyz, full_xyz)
    for name in ("occlusion", "nhit", "nmiss", "nocc"):
        np.testing.assert_array_equal(crop_pts[name], full_pts[name])


@needs_untwine
def test_buffered_tiles_match_single_grid_on_short_rays(tmp_path, monkeypatch):
    import rasterio

    cloud = tmp_path / "flight.las"
    time, xyz = flight(cloud, seconds=2, rate=500)
    reference = occlusion.run([cloud], time, xyz, tmp_path / "full.copc.laz", tmp_path / "full.tif", voxel=1)
    nx, ny, nz = reference["grid"]["shape"]
    monkeypatch.setattr(occlusion, "WORKING_BYTES", 0)
    budget = occlusion._grid_bytes((math.ceil(nx * 0.8), math.ceil(ny * 0.8), nz), nx * ny)
    tiled = occlusion.run([cloud], time, xyz, tmp_path / "tile.copc.laz", tmp_path / "tile.tif", voxel=1, budget=budget)
    assert tiled["tiling"]["tiles"] > 1
    with rasterio.open(tmp_path / "full.tif") as first, rasterio.open(tmp_path / "tile.tif") as second:
        a, b = first.read(), second.read()
    valid = np.isfinite(a) | np.isfinite(b)
    assert np.isclose(a[valid], b[valid], atol=0.01).mean() >= 0.99
    assert abs(reference["occluded_share"] - tiled["occluded_share"]) <= 0.01


def test_tiled_plan_crops_the_local_terrain_envelope(monkeypatch):
    monkeypatch.setattr(occlusion, "WORKING_BYTES", 0)
    grid = occlusion.Grid(np.zeros(3), (80, 60, 150), 1)
    ground = np.broadcast_to(np.where(np.arange(80) < 40, 1, 120)[:, None], grid.shape[:2])
    top = ground + 10
    budget = occlusion._grid_bytes((40, 60, 16), 80 * 60)
    layout, halo = occlusion._tile_layout(grid, top, ground, budget, 0)
    assert len(layout) == 2 and halo == 0
    for tile in layout:
        xy = tile.extent[:2]
        assert tile.extent[2].start == ground[xy].min() - 1
        assert tile.extent[2].stop > top[xy].max()
        assert tile.shape[2] < 16


@needs_untwine
def test_delivered_miss_flag_is_written_without_replacing_path(tmp_path):
    import rasterio

    cloud = tmp_path / "flight.las"
    time, xyz = flight(cloud, seconds=3, rate=500)
    delivered, cog, copc, report = [tmp_path / name for name in
                                   ("path.csv", "out.tif", "out.copc.laz", "report.json")]
    _delivered(delivered, time[::50], xyz[::50] + [0.5, 0, 0])
    assert main([str(cloud), "--trajectory", str(delivered), "--occlusion-cog", str(cog), "--occlusion-copc", str(copc),
                 "--report", str(report)]) == 0
    data = json.loads(report.read_text())
    assert data["occlusion"]["trajectory"] == "delivered"
    check = data["occlusion"]["delivered_path_check"]
    assert check["threshold_exceeded"] and check["median_m"] > 0.1
    assert {s["name"] for s in data["stages"]} >= {"read/time-sort", "platform check", "delivered-path check",
                                                  "grid extent/crop pass", "traversal", "classify + ground",
                                                  "COPC write", "COG write"}
    with rasterio.open(cog) as raster:
        assert raster.tags()["OCCLUSION_DELIVERED_MISS_FLAG"] == "true"
        assert raster.tags()["OCCLUSION_TRAJECTORY"] == "delivered"


@needs_untwine
def test_single_return_scan_angle_occlusion_is_rough(tmp_path):
    import rasterio

    cloud = tmp_path / "single.laz"
    flight(cloud, seconds=5, rate=1400, single_return=True)
    trajectory, cog, copc = [tmp_path / name for name in ("path.csv", "out.tif", "out.copc.laz")]
    assert main([str(cloud), "-o", str(trajectory), "--occlusion-cog", str(cog), "--occlusion-copc", str(copc)]) == 0
    data = json.loads(trajectory.with_suffix(".report.json").read_text())
    assert data["method"] == "scan_angle" and data["occlusion"]["quality"] == "rough"
    assert data["occlusion"]["traversed_pulses"] > 0
    with rasterio.open(cog) as raster:
        assert raster.tags()["OCCLUSION_QUALITY"] == "rough"


@needs_untwine
def test_single_returns_without_an_angle_need_a_delivered_path(tmp_path):
    cloud = tmp_path / "single.laz"
    time, xyz = flight(cloud, seconds=5, rate=500, single_return=True)
    data = laspy.read(cloud)
    data.scan_angle[:] = 0
    data.write(cloud)
    path = tmp_path / "estimated.csv"
    report = tmp_path / "report.json"
    products = ["--occlusion-copc", str(tmp_path / "o.copc.laz"), "--occlusion-cog", str(tmp_path / "o.tif")]
    args = [str(cloud), "--report", str(report), *products]
    assert main([*args, "-o", str(path)]) == 2
    assert "Single returns without scan angle" in json.loads(report.read_text())["error"]

    delivered = tmp_path / "delivered.csv"
    _delivered(delivered, time[::50], xyz[::50])
    assert main([*args, "--trajectory", str(delivered)]) == 0
    data = json.loads(report.read_text())["occlusion"]
    assert data["trajectory"] == "delivered"
    assert data["delivered_path_check"]["threshold_exceeded"] is None
