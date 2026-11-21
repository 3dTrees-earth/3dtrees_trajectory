import json

import laspy
import numpy as np
import pytest

from synthetic import flight
from threedtrees_trajectory import pulses
from threedtrees_trajectory.cli import main
from test_occlusion import needs_untwine


def test_time_tolerance_uses_channel_and_crosses_read_boundaries(tmp_path, monkeypatch):
    cloud = tmp_path / "flight.laz"
    flight(cloud, seconds=6, rate=400)
    las = laspy.read(cloud)
    las.gps_time[np.asarray(las.return_number) == 2] += 0.5e-6
    las.write(cloud)
    monkeypatch.setattr(pulses, "CHUNK", 313)
    sample, info = pulses.read_pulses([cloud])
    assert info.return_modes[0]["mode"] == "time_tolerance"
    assert info.unpaired_multi_echoes == 0
    assert np.count_nonzero(sample.returns == 2) > 1000

    las = laspy.read(cloud)
    las.scanner_channel[np.asarray(las.return_number) == 2] = 1
    las.write(cloud)
    sample, info = pulses.read_pulses([cloud])
    assert info.return_modes[0]["mode"] == "single"
    assert len(sample) == len(las.points)
    assert np.all(sample.returns == 1)


def test_unpaired_file_falls_back_without_changing_other_file(tmp_path):
    good, broken = tmp_path / "good.laz", tmp_path / "broken.laz"
    flight(good, seconds=5, rate=400)
    flight(broken, seconds=5, rate=400, seed=1)
    las = laspy.read(broken)
    las.gps_time[:] += 5
    las.number_of_returns[:] = 15
    las.write(broken)
    sample, info = pulses.read_pulses([good, broken])
    assert [item["mode"] for item in info.return_modes] == ["paired", "single"]
    assert info.single_return_files == [1]
    assert np.any(sample.returns[sample.file == 0] > 1)
    assert np.all(sample.returns[sample.file == 1] == 1)


@needs_untwine
@pytest.mark.parametrize("delivered", [False, True])
def test_unpaired_points_reach_single_return_traversal(tmp_path, delivered):
    cloud, output = tmp_path / "flight.laz", tmp_path / "estimated.csv"
    time, xyz = flight(cloud, seconds=5, rate=500)
    las = laspy.read(cloud)
    las.number_of_returns[:] = 15
    if delivered:
        las.scan_angle[:] = 0
    las.write(cloud)
    args = [str(cloud), "-o", str(output), "--voxel", "2",
            "--occlusion-cog", str(tmp_path / "o.tif"), "--occlusion-copc", str(tmp_path / "o.copc.laz")]
    if delivered:
        path = tmp_path / "delivered.csv"
        np.savetxt(path, np.column_stack((time[::25], xyz[::25])), delimiter=",", header="time,x,y,z", comments="")
        args += ["--trajectory", str(path)]
    assert main(["--debug", *args]) == 0
    report = json.loads(output.with_suffix(".report.json").read_text())
    assert report["pulse_check"]["return_modes"][0]["mode"] == "single"
    assert report["occlusion"]["quality"] == "rough"
    assert report["occlusion"]["single_return_fallback_files"] == [str(cloud)]
    assert report["occlusion"]["trajectory"] == ("delivered" if delivered else "estimated")
    assert report["occlusion"]["traversed_pulses"] > 0
    assert report["occlusion"]["missing_returns"] == 0


def test_constant_clock_refuses_a_moving_delivered_path(tmp_path, monkeypatch):
    cloud, report = tmp_path / "flight.laz", tmp_path / "report.json"
    flight(cloud, seconds=2, rate=400)
    las = laspy.read(cloud)
    las.gps_time[:] = 300000
    las.write(cloud)
    delivered = tmp_path / "delivered.csv"
    delivered.write_text("time,x,y,z\n299999,0,0,10\n300001,1,0,10\n")
    monkeypatch.setattr(pulses, "CHUNK", 100)
    assert main([str(cloud), "--debug", "--report", str(report), "--trajectory", str(delivered),
                 "--occlusion-cog", str(tmp_path / "o.tif"), "--occlusion-copc", str(tmp_path / "o.copc.laz")]) == 3
    assert "cannot be matched to a constant point clock" in json.loads(report.read_text())["error"]


@needs_untwine
@pytest.mark.parametrize("declared_single", [False, True])
def test_single_return_file_fills_missing_time_windows(tmp_path, declared_single):
    good, bad = tmp_path / "good.laz", tmp_path / "bad.laz"
    flight(good, seconds=6, rate=500)
    flight(bad, seconds=6, rate=500, single_return=declared_single)
    cloud = laspy.read(bad)
    cloud.gps_time[:] += 6
    if not declared_single:
        cloud.number_of_returns[:] = 15
    cloud.write(bad)
    out = tmp_path / "trajectory.csv"
    assert main([str(good), str(bad), "--debug", "-o", str(out), "--voxel", "2",
                 "--occlusion-cog", str(tmp_path / "o.tif"), "--occlusion-copc", str(tmp_path / "o.copc.laz")]) == 0
    report = json.loads(out.with_suffix(".report.json").read_text())
    assert report["method"] == "lines_scan_angle" and report["quality"] == "rough"
    assert report["occlusion"]["quality"] == "rough"
    assert report["time_span"][0] < 300001 < 300011 < report["time_span"][1]


def test_small_all_unpaired_file_still_takes_single_return_fallback(tmp_path):
    cloud = tmp_path / "few.las"
    flight(cloud, seconds=1, rate=30)
    data = laspy.read(cloud)
    data.number_of_returns[:] = 10
    data.write(cloud)
    sample, info = pulses.read_pulses([cloud])
    assert info.single_return_files == [0]
    assert len(sample) == len(data.points)
