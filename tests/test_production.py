import json

import laspy
import numpy as np
import pytest

from legacy_occlusion import trace_reference
from synthetic import flight
from test_occlusion import needs_untwine
from threedtrees_trajectory import occlusion, pulses
from threedtrees_trajectory.cli import main
from threedtrees_trajectory.cloud import TimeSortedCopc


def _delivered(path, time, xyz):
    np.savetxt(path, np.column_stack((time, xyz)), delimiter=",", header="time,x,y,z", comments="")


def _copc_tags(path):
    with laspy.open(path) as reader:
        record, = (item for item in reader.evlrs if (item.user_id, item.record_id) == ("3dtrees", 1))
        return json.loads(record.record_data_bytes())


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
    assert source.keys.stat().st_size == source.points_read * np.dtype(np.float64).itemsize
    first = [chunk.gps_time.copy() for chunk in source.chunks()]
    assert max(map(len, first)) <= 15
    assert all(p[-1] < q[0] for p, q in zip(first, first[1:]))
    gps = np.concatenate(first)
    assert np.all(gps[1:] >= gps[:-1])
    assert source.points_read == len(gps)
    from_spool = np.concatenate([chunk.gps_time for chunk in source.spooled_chunks()])
    np.testing.assert_array_equal(np.sort(from_spool), gps)

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


def test_laz_traversal_recovers_bounded_separated_echoes(tmp_path, monkeypatch):
    original, separated = tmp_path / "flight.las", tmp_path / "separated.laz"
    time, xyz = flight(original, seconds=1, rate=100)
    data = laspy.read(original)
    order = np.argsort(np.asarray(data.return_number), kind="stable")
    data.points = data.points[order]
    data.write(separated)
    grid, _ = occlusion._grid([original], 1)
    monkeypatch.setattr("legacy_occlusion.CHUNK", 47)

    def trace(path):
        stats = dict(echoes_read=0, echoes_outside_path=0, traversed_pulses=0, missing_returns=0)
        return trace_reference([path], grid, time, xyz, stats), stats

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
    crop = occlusion._crop

    def without_crop(grid, *args, **kwargs):
        clipped, top, ground, spikes = crop(grid, *args, **kwargs)
        offset = round((clipped.origin[2] - grid.origin[2]) / grid.voxel)
        return grid, np.where(top >= 0, top + offset, top), ground, spikes

    results = []
    for index in range(2):
        if index:
            monkeypatch.setattr(occlusion, "_crop", without_crop)
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


def test_window_plan_splits_on_memory():
    from pathlib import Path

    from threedtrees_trajectory.time_windows import plan
    from threedtrees_trajectory.window_input import ScanBlock, ScanIndex

    grid = occlusion.Grid(np.zeros(3), (400, 10, 60), 1)
    index = ScanIndex([Path("acquisition.laz")], None)
    index.blocks[0] = [
        ScanBlock(0, 100, 0., 1., np.array([1., 1., 0.]), np.array([199., 2., 2.])),
        ScanBlock(100, 200, 2., 3., np.array([201., 1., 0.]), np.array([399., 2., 2.])),
    ]
    time = np.array([0., 1., 2., 3.])
    xyz = np.array([[1., 1., 3.], [199., 1., 3.], [201., 1., 3.], [399., 1., 3.]])
    layers = np.zeros(grid.shape[:2], np.int32), np.full(grid.shape[:2], grid.shape[2], np.int32)
    actual = plan(index, grid, time, xyz, (0., 4.), 1, occlusion.DEFAULT_BUDGET, layers)
    assert len(actual) == 1
    split = plan(index, grid, time, xyz, (0., 4.), 1, actual[0].memory - 1, layers)
    assert len(split) == 2
    assert sum(window.echoes for window in split) == 200


@needs_untwine
@pytest.mark.parametrize(("debug", "combined"), [(True, False), (False, False), (False, True)])
def test_delivered_miss_flag_is_written_without_replacing_path(tmp_path, debug, combined):
    import rasterio

    cloud = tmp_path / "flight.las"
    time, xyz = flight(cloud, seconds=3, rate=500)
    delivered, cog, copc, report = [tmp_path / name for name in
                                   ("path.csv", "out.tif", "out.copc.laz", "report.json")]
    _delivered(delivered, time[::50], xyz[::50] + [0.5, 0, 0])
    original = delivered.read_bytes()
    output = tmp_path / "estimated.csv"
    args = [str(cloud), "--trajectory", str(delivered), "--occlusion-cog", str(cog), "--occlusion-copc", str(copc)]
    if combined:
        args += ["-o", str(output)]
    if debug:
        args += ["--debug", "--report", str(report)]

    assert main(args) == 0
    assert delivered.read_bytes() == original
    if combined:
        assert "# delivered_path_miss_flag: true" in output.read_text().splitlines()
    if debug:
        data = json.loads(report.read_text())
        assert data["occlusion"]["trajectory"] == "delivered"
        check = data["occlusion"]["delivered_path_check"]
        assert check["threshold_exceeded"] and check["median_m"] > 0.1
        assert {s["name"] for s in data["stages"]} >= {"read/time-sort", "platform check", "delivered-path check",
                                                      "finalize extent/crop", "traversal", "classify + ground",
                                                      "COPC write", "COG write"}
    else:
        assert not list(tmp_path.rglob("*.json"))
    with rasterio.open(cog) as raster:
        tags = raster.tags()
    for metadata in (tags, _copc_tags(copc)):
        assert metadata["OCCLUSION_DELIVERED_MISS_FLAG"] == "true"
        assert metadata["OCCLUSION_TRAJECTORY"] == "delivered"
        assert metadata["OCCLUSION_PATH_METHOD"] == "delivered"
    with laspy.CopcReader.open(copc) as reader:
        assert len(reader.query()) == reader.header.point_count


@needs_untwine
def test_single_return_scan_angle_occlusion_is_rough(tmp_path):
    import rasterio

    cloud = tmp_path / "single.laz"
    flight(cloud, seconds=5, rate=1400, single_return=True)
    trajectory, cog, copc = [tmp_path / name for name in ("path.csv", "out.tif", "out.copc.laz")]
    assert main([str(cloud), "--debug", "-o", str(trajectory), "--occlusion-cog", str(cog), "--occlusion-copc", str(copc)]) == 0
    data = json.loads(trajectory.with_suffix(".report.json").read_text())
    assert data["method"] == "scan_angle" and data["occlusion"]["quality"] == "rough"
    assert data["occlusion"]["traversed_pulses"] > 0
    csv_lines = trajectory.read_text().splitlines()
    assert data["warnings"]
    assert [line.removeprefix("# warning: ") for line in csv_lines if line.startswith("# warning: ")] == data["warnings"]
    with rasterio.open(cog) as raster:
        tags = raster.tags()
    for metadata in (tags, _copc_tags(copc)):
        assert metadata["OCCLUSION_QUALITY"] == "rough"
        assert metadata["OCCLUSION_PATH_METHOD"] == "scan_angle"
        assert metadata["OCCLUSION_TRAJECTORY"] == "estimated"
        assert json.loads(metadata["OCCLUSION_WARNINGS"]) == data["warnings"]


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
    args = [str(cloud), "--debug", "--report", str(report), *products]
    assert main([*args, "-o", str(path)]) == 4
    assert "Single returns without scan angle" in json.loads(report.read_text())["error"]

    delivered = tmp_path / "delivered.csv"
    _delivered(delivered, time[::50], xyz[::50])
    assert main([*args, "--trajectory", str(delivered)]) == 0
    data = json.loads(report.read_text())["occlusion"]
    assert data["trajectory"] == "delivered"
    assert data["delivered_path_check"]["threshold_exceeded"] is None
