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
    assert main([str(tmp_path / "f.laz"), "-o", str(out), "--occlusion-copc", str(copc), "--occlusion-cog", str(cog)]) == 0
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
    assert stages >= {"read/time-sort", "platform check", "path estimation", "grid extent/crop pass",
                      "plan GPS-time windows", "traversal", "classify + ground", "COPC write", "COG write"}
    for stage in report["stages"]:
        assert stage["start_seconds"] >= 0
        assert stage["wall_seconds"] >= 0 and stage["cpu_seconds"] >= 0
        assert stage["peak_rss_bytes"] > 0

    section = report["occlusion"]
    assert section["trajectory"] == "estimated" and section["ground"] == "lowest_hit"
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
        assert float(raster.tags()["OCCLUSION_BAND_M"]) == 1
        assert raster.res == (0.5, 0.5)
        names = [f"occluded_{band}-{band + 1}m" for band in range(raster.count - 1)] + ["unobserved_column"]
        assert list(raster.descriptions) == section["cog_bands"] == names
        values = raster.read()
    finite = values[np.isfinite(values)]
    assert finite.size and finite.min() >= 0 and finite.max() <= 1


@needs_untwine
def test_delivered_path_is_used_without_estimating(tmp_path):
    time, sensor = flight(tmp_path / "f.laz")
    path = tmp_path / "delivered.csv"
    rows = "\n".join(f"{t:.6f},{x:.4f},{y:.4f},{z:.4f}" for t, (x, y, z) in zip(time[::400], sensor[::400]))
    path.write_text(f"# format: 3dtrees-trajectory/1\ntime,x,y,z\n{rows}\n")
    copc, cog, report = tmp_path / "o.copc.laz", tmp_path / "o.tif", tmp_path / "r.json"
    assert main([str(tmp_path / "f.laz"), "--trajectory", str(path), "--report", str(report),
                 "--occlusion-copc", str(copc), "--occlusion-cog", str(cog), "--voxel", "1"]) == 0
    data = json.loads(report.read_text())
    assert data["trajectory"] == "delivered" and data["occlusion"]["trajectory"] == "delivered"
    assert copc.stat().st_size > 0 and cog.stat().st_size > 0


@pytest.mark.parametrize("offset", [100, 400_000_000])
def test_rejects_path_without_clock_overlap(tmp_path, offset):
    cloud = tmp_path / "flight.las"
    time, sensor = flight(cloud, seconds=2, rate=100)
    path, copc, cog, report = [tmp_path / name for name in ("delivered.csv", "out.laz", "out.tif", "report.json")]
    np.savetxt(path, np.column_stack((time[::10] + offset, sensor[::10])), delimiter=",",
               header="time,x,y,z", comments="")
    assert main([str(cloud), "--trajectory", str(path), "--report", str(report),
                 "--occlusion-copc", str(copc), "--occlusion-cog", str(cog)]) == 2
    assert "no overlapping GPS times" in json.loads(report.read_text())["occlusion"]["error"]
    assert copc.stat().st_size == cog.stat().st_size == 0


def test_delivered_path_does_not_bypass_static_cloud_refusal(tmp_path):
    cloud = tmp_path / "scan.laz"
    scan(cloud, [1000, 2000, 5], rays=2000)
    path = tmp_path / "path.csv"
    path.write_text("time,x,y,z\n100000,1000,2000,5\n100001,1001,2000,5\n")
    report = tmp_path / "report.json"
    assert main([str(cloud), "--trajectory", str(path), "--report", str(report),
                 "--occlusion-copc", str(tmp_path / "out.laz"), "--occlusion-cog", str(tmp_path / "out.tif")]) == 2
    assert "static (TLS)" in report.read_text()


def test_streamed_and_unsorted_multi_return_counts_agree(tmp_path, monkeypatch):
    sorted_cloud, unsorted_cloud = tmp_path / "sorted.las", tmp_path / "unsorted.las"
    time, sensor = flight(sorted_cloud, seconds=2, rate=101)
    las = laspy.read(sorted_cloud)
    order = np.argsort(las.gps_time, kind="stable")[::-1]
    las.points = las.points[order]
    las.write(unsorted_cloud)
    monkeypatch.setattr(occlusion, "CHUNK", 77)
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


def test_dtm_without_coverage_is_rejected(tmp_path):
    import rasterio
    from rasterio.transform import from_origin

    dtm = tmp_path / "out-of-extent.tif"
    with rasterio.open(dtm, "w", driver="GTiff", width=3, height=3, count=1, dtype="float32",
                       transform=from_origin(50, 60, 1, 1)) as out:
        out.write(np.ones((3, 3), np.float32), 1)
    with pytest.raises(ValueError, match="does not cover"):
        occlusion.ground_from_dtm(dtm, occlusion.Grid(np.zeros(3), (5, 5, 40), 0.5))


def test_ground_above_canopy_gives_empty_copc_and_readable_cog(tmp_path):
    import rasterio
    from rasterio.transform import from_origin

    cloud, dtm, copc, cog = [tmp_path / name for name in ("flight.las", "dtm.tif", "out.laz", "out.tif")]
    time, sensor = flight(cloud, seconds=1, rate=101)
    with rasterio.open(dtm, "w", driver="GTiff", width=20, height=10, count=1, dtype="float32",
                       transform=from_origin(900, 2050, 10, 10)) as out:
        out.write(np.full((10, 20), 1000, np.float32), 1)
    stats = occlusion.run([cloud], time, sensor, copc, cog, dtm=dtm)
    assert stats["copc_points"] == copc.stat().st_size == 0
    with rasterio.open(cog) as raster:
        assert raster.tags(ns="IMAGE_STRUCTURE")["LAYOUT"] == "COG"
        assert raster.descriptions == ("occluded_0-1m", "unobserved_column")
        assert np.isnan(raster.read()).all()


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
