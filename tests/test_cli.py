import json
import os
import subprocess
import sys
from pathlib import Path

import laspy
import numpy as np
import pytest
from pyproj import CRS

from synthetic import flight
from threedtrees_trajectory import pulses
from threedtrees_trajectory.cli import main
from threedtrees_trajectory.estimate import Result, estimate
from threedtrees_trajectory.output import write_csv
from threedtrees_trajectory.pulses import CloudInfo
from threedtrees_trajectory.solve import Fit

POINT_FORMATS = [(pf, "1.2" if pf < 4 else "1.3" if pf < 6 else "1.4") for pf in range(11)]


def _report(output):
    return json.loads(output.with_suffix(".report.json").read_text())


def _assert_refusal(output, reason):
    assert output.read_bytes() == b""
    assert output.with_suffix(".geojson").read_bytes() == b""
    report = _report(output)
    assert reason in report["error"]
    assert report["positions"] == 0
    assert report["quality"] == "none"
    assert report["resources"]["peak_rss_bytes"] > 0
    assert report["resources"]["wall_seconds"] >= 0
    assert report["resources"]["cpu_seconds"] >= 0
    return report


@pytest.mark.parametrize(("point_format", "version"), POINT_FORMATS)
def test_las_point_formats(tmp_path, point_format, version):
    # Cover both readers, all LAS point formats and all supported file versions.
    path = tmp_path / ("input.las" if point_format % 2 == 0 else "input.laz")
    time, sensor = flight(path, seconds=6, rate=400, point_format=point_format, version=version)
    output = tmp_path / "output.csv"
    assert path.stat().st_size < 1_000_000
    if point_format in (0, 2):
        assert main([str(path), "-o", str(output)]) == 2
        _assert_refusal(output, "needs GPS time")
        return
    assert main([str(path), "-o", str(output)]) == 0
    rows = np.loadtxt(output, delimiter=",", skiprows=14, usecols=(0, 1, 2, 3))
    assert len(rows) == 12
    assert np.max(np.abs(rows[:, 1:] - np.column_stack([
        np.interp(rows[:, 0], time, coord) for coord in sensor.T
    ]))) < 0.05
    report = _report(output)
    assert report["crs"] == "unknown"
    assert report["resources"]["points_read"] == report["points"]


def test_multifile_with_matched_crs_and_different_las_versions(tmp_path, monkeypatch):
    first, second = (tmp_path / name for name in ("first.las", "second.laz"))
    time, sensor = flight(first, seconds=10, rate=400, point_format=1, version="1.2")
    las = laspy.read(first)
    cut = np.searchsorted(las.gps_time, time[2000])
    tail = laspy.convert(las, point_format_id=8, file_version="1.4")
    tail.points = tail.points[cut:]
    tail.header.add_crs(CRS.from_epsg(32632))
    tail.write(second)
    las.points = las.points[:cut]
    las.write(first)
    monkeypatch.setattr(pulses, "CHUNK", 313)
    monkeypatch.setattr(pulses, "PROBE", 100)
    output = tmp_path / "output.csv"
    assert main(["-o", str(output), str(first), str(second)]) == 0
    rows = np.loadtxt(output, delimiter=",", skiprows=14, usecols=(0, 1, 2, 3))
    assert len(rows) == 20
    assert np.max(np.abs(rows[:, 1:] - np.column_stack([
        np.interp(rows[:, 0], time, coord) for coord in sensor.T
    ]))) < 0.05
    assert (_report(output)["resources"]["points_read"] ==
            len(las.points) + len(tail.points))
    layer = json.loads(output.with_suffix(".geojson").read_text())
    assert layer["features"][0]["properties"]["crs"] == "EPSG:4326"


@pytest.mark.parametrize(("scenario", "reason"), [
    ("empty", "has no points"),
    ("missing", "existing regular file"),
    ("invalid_header", "unable to read the LAS/LAZ header"),
    ("truncated", "Truncated point data"),
    ("constant_time", "GPS time must vary"),
    ("nonfinite_time", "GPS times must be finite"),
    ("mixed_clocks", "same GPS time system"),
    ("mixed_crs", "same coordinate reference system"),
    ("nonmetric_crs", "must use metres"),
    ("no_fit", "No reliable position could be fitted"),
])
def test_bad_input_has_refusal_and_clears_stale_outputs(tmp_path, scenario, reason):
    input_path, other = tmp_path / "input.las", tmp_path / "other.las"
    paths = [str(input_path)]
    match scenario:
        case "empty":
            laspy.create(point_format=6).write(input_path)
        case "missing":
            pass
        case "invalid_header":
            input_path.write_bytes(b"LASF")
        case "truncated":
            flight(input_path, seconds=2, rate=400)
            data = input_path.read_bytes()
            input_path.write_bytes(data[:-30])
        case "constant_time" | "nonfinite_time":
            flight(input_path, seconds=2, rate=400)
            cloud = laspy.read(input_path)
            cloud.gps_time[:] = 0 if scenario == "constant_time" else np.nan
            cloud.write(input_path)
        case "mixed_clocks" | "mixed_crs":
            flight(input_path, seconds=2, rate=400)
            cloud = laspy.read(input_path)
            if scenario == "mixed_clocks":
                cloud.header.global_encoding.gps_time_type = laspy.header.GpsTimeType.STANDARD
            else:
                cloud.header.add_crs(CRS.from_epsg(32632))
                cloud.write(input_path)
                cloud.header.add_crs(CRS.from_epsg(32633))
            cloud.write(other)
            paths.append(str(other))
        case "nonmetric_crs":
            flight(input_path, seconds=2, rate=400)
            cloud = laspy.read(input_path)
            cloud.header.add_crs(CRS.from_epsg(4326))
            cloud.write(input_path)
        case "no_fit":
            flight(input_path, seconds=1, rate=100, single_return=True)
    output = tmp_path / "output.csv"
    output.write_text("stale")
    output.with_suffix(".geojson").write_text("stale")
    assert main([*paths, "-o", str(output)]) == 2
    _assert_refusal(output, reason)


@pytest.mark.parametrize(("option", "reason"), [
    ("--interval=inf", "positive finite"),
    ("--interval=0", "positive finite"),
    ("--interval=abc", "positive finite"),
    ("--threads=0", "positive integer"),
    ("--threads=2.5", "positive integer"),
])
def test_invalid_arguments_have_a_refusal(tmp_path, option, reason):
    output = tmp_path / "output.csv"
    assert main(["input.las", "-o", str(output), option]) == 2
    assert _assert_refusal(output, reason)["resources"]["points_read"] == 0


def test_outputs_cannot_alias_the_input(tmp_path):
    input_path = tmp_path / "input.las"
    input_path.write_bytes(b"unchanged")
    output = tmp_path / "alias.las"
    output.symlink_to(input_path)
    assert main([str(input_path), "-o", str(output)]) == 2
    assert input_path.read_bytes() == b"unchanged"
    assert "differ from input" in _report(output)["error"]


def test_zero_gps_time_stays_in_csv(tmp_path):
    info = CloudInfo("unknown", "gps_sow", 1, ["input.laz"])
    result = Result("lines", [Fit(0.0, np.zeros(3), 50, 0.0)])
    output = tmp_path / "output.csv"
    write_csv(output, info, result)
    assert "# time_system: gps_sow" in output.read_text()
    assert output.read_text().splitlines()[-1].startswith("0.000000,0.0000,")


def test_decoder_threads_do_not_change_trajectory_bytes(tmp_path):
    path = tmp_path / "compressed.laz"
    flight(path, seconds=6, rate=1200)
    results = []
    for name, slots in (("one", "1"), ("four", "4")):
        output = tmp_path / name / "trajectory.csv"
        env = os.environ | {"GALAXY_SLOTS": slots}
        entrypoint = Path(__file__).resolve().parents[1] / "src" / "run.py"
        subprocess.run([sys.executable, str(entrypoint), str(path), "-o", str(output), "--threads", "16"],
                       env=env, check=True)
        results.append((output.read_bytes(), output.with_suffix(".geojson").read_bytes()))
        assert 1 <= _report(output)["resources"]["decoding_threads"] <= int(slots)
    assert results[0] == results[1]


def test_one_position_skips_the_map_layer(tmp_path):
    path = tmp_path / "input.laz"
    flight(path, seconds=6, rate=1600, single_return=True)
    output = tmp_path / "output.csv"
    assert main([str(path), "-o", str(output), "--interval", "10"]) == 0
    assert _report(output)["positions"] == 1
    assert output.with_suffix(".geojson").read_bytes() == b""


def test_position_fit_handles_sampled_chunks(tmp_path, monkeypatch):
    path = tmp_path / "input.laz"
    time, sensor = flight(path, seconds=6, rate=1200)
    monkeypatch.setattr(pulses, "CHUNK", 801)
    cloud, info = pulses.read_pulses([path], max_multi=3000, max_single=3000)
    assert info.kept_every[0] > 1 and info.points_read == info.points
    assert len(cloud) < info.points // 2
    result = estimate(cloud)
    truth = np.column_stack([np.interp([f.time for f in result.fits], time, coord) for coord in sensor.T])
    assert len(result.fits) == 12
    assert np.max(np.abs(np.array([f.xyz for f in result.fits]) - truth)) < 0.05
