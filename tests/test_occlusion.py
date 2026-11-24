import json
import os
from pathlib import Path

import laspy
import numpy as np
import pytest

from legacy_occlusion import trace_reference
from synthetic import _write, flight, scan
from threedtrees_trajectory import occlusion
from threedtrees_trajectory.cli import main

UNTWINE = Path(occlusion.UNTWINE)
needs_untwine = pytest.mark.skipif(not os.access(UNTWINE, os.X_OK), reason="untwine not installed")


def test_classify_follows_occpy():
    nhit = np.array([1, 0, 0, 0])
    nmiss = np.array([5, 3, 0, 0])
    nocc = np.array([2, 1, 4, 0])
    assert occlusion.classify(nhit, nmiss, nocc).tolist() == [1, 2, 3, 4]


def test_ground_from_hits_replaces_vegetation_spikes():
    grid = occlusion.Grid(np.zeros(3), (30, 30, 40), 0.5)
    filled = np.zeros(grid.shape, bool)
    filled[:, :, 0] = True  # ground in the lowest layer
    filled[15, 15, 0] = False
    filled[15, 15, 20] = True  # no pulse reached the ground here: lowest hit is a branch 10 m up
    lowest, _ = occlusion._surface_indices(filled, grid)
    ground, spikes = occlusion._ground_from_lowest(lowest, grid)
    assert spikes == 1
    assert ground[15, 15] == pytest.approx(0.25)
    assert np.allclose(ground, 0.25)


def test_read_trajectory_needs_times(tmp_path):
    path = tmp_path / "stations.csv"
    path.write_text("# kind: stations\ntime,x,y,z,roll,pitch,yaw,name\n,1,2,3,,,,scan0\n")
    with pytest.raises(ValueError, match="time for every position"):
        occlusion.read_trajectory(path)


def test_occlusion_needs_both_products(tmp_path, capsys):
    flight(tmp_path / "f.laz")
    out = tmp_path / "p.csv"
    assert main([str(tmp_path / "f.laz"), "-o", str(out), "--occlusion-copc", str(tmp_path / "o.copc.laz")]) == 2
    assert "both products" in capsys.readouterr().err


@needs_untwine
def test_estimated_path_gives_trajectory_2_report_copc_and_cog(tmp_path):
    import csv
    import rasterio

    flight(tmp_path / "f.laz", seconds=6, rate=500)
    out, copc, cog = tmp_path / "p.csv", tmp_path / "o.copc.laz", tmp_path / "o.tif"
    assert main([str(tmp_path / "f.laz"), "--debug", "-o", str(out), "--occlusion-copc", str(copc), "--occlusion-cog", str(cog)]) == 0
    lines = out.read_text().splitlines()
    assert lines[0] == "# format: 3dtrees-trajectory/2"
    header, *rows = csv.reader(line for line in lines if not line.startswith("#"))
    assert header == ["time", "x", "y", "z", "roll", "pitch", "yaw", "name"]
    assert rows and all(len(row) == 8 and row[4:] == ["", "", "", ""] for row in rows)

    report = json.loads(out.with_suffix(".report.json").read_text())
    assert report["exit_code"] == 0 and report["positions"] == len(rows)
    assert report["method"] == "lines" and report["quality"] == "good"
    assert report["resources"]["points_read"] == report["points"]
    stages = {stage["name"] for stage in report["stages"]}
    assert stages >= {"read/time-sort", "platform check", "path estimation", "finalize extent/crop",
                      "plan GPS-time windows", "traversal", "classify voxels + write products", "COPC write", "COG write"}
    for stage in report["stages"]:
        assert stage["start_seconds"] >= 0
        assert stage["wall_seconds"] >= 0 and stage["cpu_seconds"] >= 0
        assert stage["peak_rss_bytes"] > 0

    section = report["occlusion"]
    assert section["trajectory"] == "estimated" and section["ground"] == "lowest_hit"
    assert section["cog_band_step_m"] == 0.5
    with laspy.open(copc) as reader:
        assert reader.header.point_count == section["copc_points"] > 0
        assert any(type(vlr).__name__ == "CopcInfoVlr" for vlr in reader.header.vlrs)
        points = reader.read()
    hit, miss, hidden = (np.asarray(points[name]) for name in ("nhit", "nmiss", "nocc"))
    expected_class = np.select([hit > 0, miss > 0, hidden > 0], [1, 2, 3], default=4)
    np.testing.assert_array_equal(points["occlusion"], expected_class)
    assert set(np.unique(points["occlusion"])) <= {1, 3, 4}
    assert set(section["voxels_ground_to_canopy"]) == {"filled", "empty", "occluded", "unobserved"}
    assert np.all(points["height_above_ground"] >= 0)
    with rasterio.open(cog) as raster:
        assert raster.tags(ns="IMAGE_STRUCTURE")["LAYOUT"] == "COG"
        assert float(raster.tags()["OCCLUSION_BAND_M"]) == section["cog_band_step_m"]
        assert raster.res == (0.5, 0.5)
        height_bands = int(np.max(points["height_above_ground"]) / 0.5) + 1
        assert raster.count == height_bands + 1
        names = [f"occluded_{i / 2:g}-{(i + 1) / 2:g}m" for i in range(height_bands)] + ["unobserved_column"]
        assert list(raster.descriptions) == section["cog_bands"] == names
        assert len(section["occluded_share_by_height_band"]) == height_bands
        values = raster.read()
    finite = values[np.isfinite(values)]
    assert finite.size and finite.min() >= 0 and finite.max() <= 1


@needs_untwine
@pytest.mark.parametrize("voxel", [1, 0.3])
def test_delivered_path_is_used_without_estimating(tmp_path, voxel):
    import rasterio

    time, sensor = flight(tmp_path / "f.laz")
    path = tmp_path / "delivered.csv"
    rows = "\n".join(f"{t:.6f},{x:.4f},{y:.4f},{z:.4f}" for t, (x, y, z) in zip(time[::400], sensor[::400]))
    path.write_text(f"# format: 3dtrees-trajectory/1\ntime,x,y,z\n{rows}\n")
    copc, cog, report = tmp_path / "o.copc.laz", tmp_path / "o.tif", tmp_path / "r.json"
    assert main([str(tmp_path / "f.laz"), "--trajectory", str(path), "--debug", "--report", str(report),
                 "--occlusion-copc", str(copc), "--occlusion-cog", str(cog), "--voxel", str(voxel)]) == 0
    data = json.loads(report.read_text())
    assert data["trajectory"] == "delivered" and data["occlusion"]["trajectory"] == "delivered"
    assert copc.stat().st_size > 0 and cog.stat().st_size > 0
    assert data["occlusion"]["cog_band_step_m"] == voxel
    with rasterio.open(cog) as raster:
        assert float(raster.tags()["OCCLUSION_BAND_M"]) == voxel
        assert raster.res == (voxel, voxel)
        assert raster.descriptions[:2] == (f"occluded_0-{voxel:g}m", f"occluded_{voxel:g}-{2 * voxel:g}m")


def test_band_shares_use_column_relative_layers():
    classes = np.array([[[1, 3, 2, 4]], [[4, 1, 4, 3]]])
    ground_layers = np.array([[2], [3]])
    inside = np.arange(2, 6) >= ground_layers[:, :, None]
    values, seen, hidden = occlusion._band_shares(classes, 2, ground_layers, inside, 4)
    np.testing.assert_allclose(values[:, 0], [[0, 0], [1, np.nan], [0, 1], [np.nan, np.nan], [0.25, 1 / 3]])
    np.testing.assert_array_equal(seen, [2, 1, 2, 0])
    np.testing.assert_array_equal(hidden, [0, 1, 1, 0])


def _terrain_fixture(tmp_path, coverage):
    import rasterio
    from pyproj import CRS
    from rasterio.transform import from_origin

    cloud = tmp_path / "flight.laz"
    flight(cloud, seconds=6, rate=500)
    points = laspy.read(cloud)
    points.header.add_crs(CRS.from_epsg(32632))
    points.write(cloud)
    grid, _ = occlusion._grid([cloud], 0.5)
    x, y = grid.origin[:2]
    nx, ny = grid.shape[:2]
    width, height = (nx + 1) // 2, (ny + 1) // 2
    if coverage == "partial":
        width -= 1
    transform = from_origin(x, y + height, 1, 1)
    terrain = np.full((height, width), 0.375, np.float32)
    if coverage == "partial":
        terrain[:, width // 2] = -9999
    elif coverage == "missing":
        terrain[:] = -9999
    path = tmp_path / "ground.tif"
    crs = "EPSG:32633" if coverage == "mismatch" else "EPSG:32632"
    with rasterio.open(path, "w", driver="GTiff", width=width, height=height, count=1, dtype=terrain.dtype,
                       crs=crs, transform=transform, nodata=-9999) as raster:
        raster.write(terrain, 1)
    return cloud, path, grid, width


@needs_untwine
@pytest.mark.parametrize(("coverage", "source"), [
    ("full", "dtm"), ("partial", "dtm+lowest_hit"), ("missing", "lowest_hit"),
])
def test_terrain_ground_in_both_products(tmp_path, coverage, source):
    import csv

    import rasterio
    from rasterio.transform import rowcol

    cloud, terrain, grid, width = _terrain_fixture(tmp_path, coverage)
    out, copc, cog = tmp_path / "p.csv", tmp_path / "occlusion.copc.laz", tmp_path / "occlusion.tif"
    assert main([str(cloud), "-o", str(out), "--debug", "--dtm", str(terrain),
                 "--occlusion-copc", str(copc), "--occlusion-cog", str(cog)]) == 0
    plot = laspy.read(copc)
    record = next(vlr for vlr in plot.evlrs if vlr.user_id == "3dtrees" and vlr.record_id == 1)
    tags = json.loads(record.record_data)
    report = json.loads(out.with_suffix(".report.json").read_text())["occlusion"]
    column = np.floor(np.asarray(plot.x) - grid.origin[0]).astype(int)
    covered = np.ones(len(plot.points), bool)
    if coverage == "partial":
        covered = (column < width) & (column != width // 2)
        assert np.any(column >= width) and np.any(column == width // 2)
    elif coverage == "missing":
        covered[:] = False
    assert covered.any() == (coverage != "missing")
    assert (~covered).any() == (coverage != "full")
    ground = plot.z - plot["height_above_ground"]
    np.testing.assert_allclose(ground[covered], 0.375)
    if (~covered).any():
        control_copc = tmp_path / "control.copc.laz"
        assert main([str(cloud), "-o", str(tmp_path / "control.csv"), "--occlusion-copc", str(control_copc),
                     "--occlusion-cog", str(tmp_path / "control.tif")]) == 0
        control = laspy.read(control_copc)
        lowest_ground = np.full(grid.shape[:2], np.nan)
        x, y = np.floor((control.xyz[:, :2] - grid.origin[:2]) / grid.voxel).astype(int).T
        lowest_ground[x, y] = control.z - control["height_above_ground"]
        x, y = np.floor((plot.xyz[~covered, :2] - grid.origin[:2]) / grid.voxel).astype(int).T
        np.testing.assert_array_equal(ground[~covered], lowest_ground[x, y])
    assert tags["OCCLUSION_GROUND"] == report["ground"] == source
    csv_lines = (line for line in out.read_text().splitlines() if not line.startswith("#"))
    assert next(csv.reader(csv_lines)) == ["time", "x", "y", "z", "roll", "pitch", "yaw", "name"]
    with rasterio.open(cog) as raster:
        assert raster.tags()["OCCLUSION_GROUND"] == source
        if coverage == "partial":
            share = (min(grid.shape[0], 2 * width) - 2) / grid.shape[0]
            assert report["ground_column_shares"] == pytest.approx({"dtm": share, "lowest_hit": 1 - share})
            for tag_name, key in (("OCCLUSION_GROUND_DTM_SHARE", "dtm"),
                                  ("OCCLUSION_GROUND_LOWEST_HIT_SHARE", "lowest_hit")):
                assert raster.tags()[tag_name] == tags[tag_name]
                assert float(tags[tag_name]) == pytest.approx(report["ground_column_shares"][key])
        observed = np.isin(plot["occlusion"], [occlusion.FILLED, occlusion.OCCLUDED])
        rows, cols = rowcol(raster.transform, plot.x[observed], plot.y[observed])
        bands = (plot["height_above_ground"][observed] / 0.5).astype(int)
        np.testing.assert_array_equal(raster.read()[bands, rows, cols],
                                      plot["occlusion"][observed] == occlusion.OCCLUDED)


def test_terrain_crs_mismatch_exits_two_in_combined_run(tmp_path, capsys):
    cloud, terrain, _, _ = _terrain_fixture(tmp_path, "mismatch")
    out, copc, cog = tmp_path / "p.csv", tmp_path / "occlusion.copc.laz", tmp_path / "occlusion.tif"
    assert main([str(cloud), "-o", str(out), "--debug", "--dtm", str(terrain),
                 "--occlusion-copc", str(copc), "--occlusion-cog", str(cog)]) == 2
    reason = "must use the same coordinate reference system"
    assert reason in capsys.readouterr().err
    report = json.loads(out.with_suffix(".report.json").read_text())
    assert report["exit_code"] == 2 and reason in report["occlusion"]["error"]
    assert out.stat().st_size == copc.stat().st_size == cog.stat().st_size == 0


@pytest.mark.parametrize("offset", [100, 400_000_000])
def test_rejects_path_without_clock_overlap(tmp_path, offset):
    cloud = tmp_path / "flight.las"
    time, sensor = flight(cloud, seconds=2, rate=100)
    path, copc, cog, report = [tmp_path / name for name in ("delivered.csv", "out.laz", "out.tif", "report.json")]
    np.savetxt(path, np.column_stack((time[::10] + offset, sensor[::10])), delimiter=",",
               header="time,x,y,z", comments="")
    assert main([str(cloud), "--trajectory", str(path), "--debug", "--report", str(report),
                 "--occlusion-copc", str(copc), "--occlusion-cog", str(cog)]) == 2
    reason = json.loads(report.read_text())["occlusion"]["error"]
    assert "no overlapping GPS times" in reason
    assert reason.endswith("Keep the delivered path as supplementary data on the dataset; it will be added once supported.")
    assert copc.stat().st_size == cog.stat().st_size == 0


def test_delivered_path_does_not_bypass_static_cloud_refusal(tmp_path):
    cloud = tmp_path / "scan.laz"
    scan(cloud, [1000, 2000, 5], rays=2000)
    path = tmp_path / "path.csv"
    path.write_text("time,x,y,z\n100000,1000,2000,5\n100001,1001,2000,5\n")
    report = tmp_path / "report.json"
    assert main([str(cloud), "--trajectory", str(path), "--debug", "--report", str(report),
                 "--occlusion-copc", str(tmp_path / "out.laz"), "--occlusion-cog", str(tmp_path / "out.tif")]) == 2
    assert "static (TLS)" in report.read_text()


def test_streamed_and_unsorted_multi_return_counts_agree(tmp_path, monkeypatch):
    sorted_cloud, unsorted_cloud = tmp_path / "sorted.las", tmp_path / "unsorted.las"
    time, sensor = flight(sorted_cloud, seconds=2, rate=101)
    las = laspy.read(sorted_cloud)
    order = np.argsort(las.gps_time, kind="stable")[::-1]
    las.points = las.points[order]
    las.write(unsorted_cloud)
    monkeypatch.setattr("legacy_occlusion.CHUNK", 77)
    grid, _ = occlusion._grid([sorted_cloud], 0.5)

    def trace(paths):
        stats = dict(echoes_read=0, echoes_outside_path=0, traversed_pulses=0, missing_returns=0)
        return trace_reference(paths, grid, time, sensor, stats)

    for expected, actual, combined in zip(trace([sorted_cloud]), trace([unsorted_cloud]),
                                          trace([sorted_cloud, unsorted_cloud])):
        np.testing.assert_array_equal(expected, actual)
        np.testing.assert_array_equal(expected * 2, combined)


@pytest.mark.parametrize("return_number", [0, 1])
def test_single_returns_are_independent_even_with_a_shared_time(tmp_path, return_number):
    cloud = tmp_path / "single.las"
    hits = np.array([[0.25, 0.25, 0.25], [1.25, 0.25, 0.25], [2.25, 0.25, 0.25]])
    returns = np.full(len(hits), return_number, np.uint8)
    _write(cloud, hits, np.array([11., 11., 12.]), returns, returns, np.zeros(len(hits)))
    grid, _ = occlusion._grid([cloud], 0.5)
    stats = dict(echoes_read=0, echoes_outside_path=0, traversed_pulses=0, missing_returns=0)
    hit, _, _ = trace_reference([cloud], grid, np.array([10., 13.]),
                                      np.array([[0.28, 0.31, 2.25], [2.31, 0.30, 2.25]]), stats)
    assert hit[:, 0, 0].tolist() == [1, 0, 1, 0, 1]
    assert stats["traversed_pulses"] == 3


@needs_untwine
def test_product_windows_preserve_cog_and_copc(tmp_path, monkeypatch):
    import rasterio

    cloud = tmp_path / "flight.las"
    time, sensor = flight(cloud, seconds=2, rate=100)
    copcs = [tmp_path / f"{i}.copc.laz" for i in range(2)]
    cogs = [tmp_path / f"{i}.tif" for i in range(2)]
    results = []
    for i, block in enumerate((256, 16)):
        monkeypatch.setattr(occlusion, "WINDOW", block)
        results.append(occlusion.run([cloud], time, sensor, copcs[i], cogs[i]))
    assert results[0] == results[1]
    with rasterio.open(cogs[0]) as single, rasterio.open(cogs[1]) as windows:
        np.testing.assert_array_equal(single.read(), windows.read())
    points = [laspy.read(copc).points.array for copc in copcs]
    points = [np.sort(array, order=("X", "Y", "Z")) for array in points]
    for name in points[0].dtype.names:
        np.testing.assert_array_equal(points[0][name], points[1][name])


def test_slam_unix_path_converts_to_gps_week_seconds():
    # Dataset 2018 (GeoSLAM via CloudCompare): path in Unix UTC, cloud in GPS week seconds.
    unix = np.array([1761057992.41066, 1761058886.62656])
    sow = occlusion.path_to_cloud_clock(unix, "unix_utc", None, "gps_sow")
    np.testing.assert_allclose(sow, unix - 315_964_800 + 18 - 2389 * 604_800)
    adjusted = occlusion.path_to_cloud_clock(unix, "unix_utc", None, "gps_adjusted")
    np.testing.assert_allclose(adjusted, unix - 315_964_800 + 18 - 1e9)
    np.testing.assert_allclose(occlusion.path_to_cloud_clock(sow, "gps_sow", 2389, "gps_adjusted"), adjusted)
    assert occlusion.path_to_cloud_clock(sow, "gps_sow", None, "gps_adjusted") is None
    assert occlusion.path_to_cloud_clock(sow, "gps_sow", 2389, "gps_sow") is None
    assert occlusion.path_to_cloud_clock(unix, "none", None, "gps_sow") is None


def test_path_clock_is_converted_only_when_raw_times_miss_the_cloud():
    time = np.arange(300000.0, 300020.0, 0.5)
    same, clock = occlusion.align_path_clock(time, "unix_utc", None, "gps_sow", (300000.0, 300020.0))
    assert clock is None and same is time
    unix = time + 2389 * 604_800 + 315_964_800 - 18
    converted, clock = occlusion.align_path_clock(unix, "unix_utc", None, "gps_sow", (300000.0, 300020.0))
    np.testing.assert_allclose(converted, time)
    assert clock["from"] == "unix_utc" and clock["to"] == "gps_sow" and clock["raw_share"] == 0
    kept, clock = occlusion.align_path_clock(unix, "none", None, "gps_sow", (300000.0, 300020.0))
    assert clock is None and kept is unix


def test_delivered_slam_path_in_unix_time_is_aligned(tmp_path):
    cloud = tmp_path / "flight.laz"
    time, sensor = flight(cloud, seconds=4, rate=500)  # GPS week seconds from 300000
    path, report = tmp_path / "slam.csv", tmp_path / "report.json"
    unix = time[::25] + 2389 * 604_800 + 315_964_800 - 18
    with path.open("w") as stream:
        stream.write("# format: 3dtrees-trajectory/1\n# time_system: unix_utc\ntime,x,y,z\n")
        np.savetxt(stream, np.column_stack((unix, sensor[::25])), delimiter=",", fmt="%.6f")
    main([str(cloud), "--trajectory", str(path), "--debug", "--report", str(report),
          "--occlusion-copc", str(tmp_path / "out.copc.laz"), "--occlusion-cog", str(tmp_path / "out.tif")])
    check = json.loads(report.read_text())["occlusion"]["delivered_path_check"]
    assert check["clock_converted"]["from"] == "unix_utc" and check["clock_converted"]["to"] == "gps_sow"
    assert check["median_m"] < 0.01  # the converted path sits on the pulse lines


@needs_untwine
def test_delivered_slam_path_in_unix_time_runs_with_warning(tmp_path):
    import rasterio

    cloud = tmp_path / "flight.laz"
    time, sensor = flight(cloud, seconds=4, rate=500)
    path = tmp_path / "slam.csv"
    with path.open("w") as stream:
        stream.write("# time_system: unix_utc\ntime,x,y,z\n")
        np.savetxt(stream, np.column_stack((time[::25] + 2389 * 604_800 + 315_964_800 - 18, sensor[::25])),
                   delimiter=",", fmt="%.6f")
    cog = tmp_path / "out.tif"
    assert main([str(cloud), "--trajectory", str(path), "--occlusion-copc", str(tmp_path / "out.copc.laz"),
                 "--occlusion-cog", str(cog)]) == 0
    with rasterio.open(cog) as raster:
        assert "converted to gps_sow" in raster.tags()["OCCLUSION_WARNINGS"]
