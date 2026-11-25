import json
import os
import subprocess
import sys
from concurrent.futures.process import BrokenProcessPool
from errno import ENOMEM
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


def test_debug_report_is_opt_in(tmp_path):
    cloud = tmp_path / "flight.laz"
    flight(cloud, seconds=6, rate=400)
    output = tmp_path / "path.csv"
    layer = output.with_suffix(".geojson")
    default_report = output.with_suffix(".report.json")
    explicit_report = tmp_path / "diagnostics" / "report.json"
    args = [str(cloud), "-o", str(output)]

    assert main(args) == 0
    assert set(tmp_path.iterdir()) == {cloud, output, layer}
    products = (output.read_bytes(), layer.read_bytes())

    reports = []
    for options, path in [([], default_report), (["--report", str(explicit_report)], explicit_report)]:
        assert main([*args, "--debug", *options]) == 0
        assert (output.read_bytes(), layer.read_bytes()) == products
        report = json.loads(path.read_text())
        assert report["exit_code"] == 0
        assert report["method"] == "lines" and report["quality"] == "good"
        assert report["resources"]["points_read"] == report["points"]
        assert {stage["name"] for stage in report["stages"]} >= {"read/time-sort", "platform check", "path estimation"}
        for stage in report["stages"]:
            assert {"wall_seconds", "cpu_seconds", "peak_rss_bytes"} <= stage.keys()
        reports.append({key: value for key, value in report.items() if key not in ("resources", "stages")})
        path.unlink()
    assert reports[0] == reports[1]


@pytest.mark.parametrize("ecef", [False, True])
def test_geocentric_crs_near_origin_is_a_local_placeholder(tmp_path, ecef):
    cloud = tmp_path / "flight.las"
    flight(cloud, seconds=6, rate=400)
    points = laspy.read(cloud)
    points.header.offsets += [(6_370_000 if ecef else 0) - 500, -1500, 10]
    points.points.offsets = points.header.offsets.copy()
    points.header.add_crs(CRS.from_epsg(4978))
    points.write(cloud)
    output = tmp_path / "path.csv"

    assert main([str(cloud), "-o", str(output), "--debug"]) == (2 if ecef else 0)
    if ecef:
        _assert_refusal(output, "Earth-centred coordinates: reproject to a projected CRS in metres")
    else:
        assert "# crs: local\n" in output.read_text()
        layer = json.loads(output.with_suffix(".geojson").read_text())["features"][0]
        assert layer["properties"]["source_crs"] == "local"
        x, y, _ = layer["geometry"]["coordinates"][0][0]
        assert x == pytest.approx(500, abs=0.02)
        assert 500 <= y < 503


def test_report_path_requires_debug(tmp_path, capsys):
    output = tmp_path / "path.csv"
    output.write_text("existing output")
    report = tmp_path / "standalone.json"
    with pytest.raises(SystemExit) as exit:
        main(["unused.las", "-o", str(output), "--report", str(report)])
    assert exit.value.code == 2
    assert "--report requires --debug" in capsys.readouterr().err
    assert output.read_text() == "existing output"
    assert not list(tmp_path.rglob("*.json"))


@pytest.mark.parametrize(("constant_clock", "code", "reason"), [
    (True, 3, "GPS time must vary"),
    (False, 4, "No reliable position could be fitted"),
])
def test_refusal_on_stderr_and_empty_products_without_debug(tmp_path, capsys, constant_clock, code, reason):
    cloud = tmp_path / "flight.las"
    flight(cloud, seconds=1, rate=100, single_return=not constant_clock)
    if constant_clock:
        points = laspy.read(cloud)
        points.gps_time[:] = 0
        points.write(cloud)
    output, geojson, copc, cog = [tmp_path / name for name in ("path.csv", "path.geojson", "out.copc.laz", "out.tif")]
    for path in (output, geojson, copc, cog):
        path.write_text("stale product")

    assert main([str(cloud), "-o", str(output), "--occlusion-copc", str(copc), "--occlusion-cog", str(cog)]) == code
    stderr = capsys.readouterr().err
    assert reason in stderr
    assert "Keep the delivered path as supplementary data on the dataset" in stderr
    assert all(path.read_bytes() == b"" for path in (output, geojson, copc, cog))
    assert not list(tmp_path.rglob("*.json"))


@pytest.mark.parametrize(("module_name", "function_name", "failure"), [
    ("pulses", "read_pulses", MemoryError()),
    ("estimate", "estimate", BrokenProcessPool("Path worker disappeared")),
    ("occlusion", "read_trajectory", OSError(ENOMEM, "Cannot allocate memory")),
    ("occlusion", "run", RuntimeError("std::bad_alloc")),
])
def test_memory_failure_is_a_job_refusal(tmp_path, monkeypatch, capsys, module_name, function_name, failure):
    from importlib import import_module

    def exhaust_memory(*args, **kwargs):
        raise failure

    module = import_module(f"threedtrees_trajectory.{module_name}")
    monkeypatch.setattr(module, function_name, exhaust_memory)
    cloud = tmp_path / "flight.las"
    flight(cloud, seconds=6, rate=400)
    out, copc, cog = tmp_path / "path.csv", tmp_path / "out.copc.laz", tmp_path / "out.tif"
    for path in (out, out.with_suffix(".geojson"), copc, cog):
        path.write_text("stale")
    args = [str(cloud), "-o", str(out), "--debug", "--occlusion-copc", str(copc), "--occlusion-cog", str(cog)]
    if function_name == "read_trajectory":
        args += ["--trajectory", str(tmp_path / "delivered.csv")]

    assert main(args) == 127
    report = _assert_refusal(out, "out of memory:")
    assert report["exit_code"] == 127
    assert "rerun with more memory" in capsys.readouterr().err
    assert copc.read_bytes() == cog.read_bytes() == b""


@pytest.mark.parametrize(("point_format", "version"), POINT_FORMATS)
def test_las_point_formats(tmp_path, point_format, version):
    # Cover both readers, all LAS point formats and all supported file versions.
    path = tmp_path / ("input.las" if point_format % 2 == 0 else "input.laz")
    time, sensor = flight(path, seconds=6, rate=400, point_format=point_format, version=version)
    output = tmp_path / "output.csv"
    assert path.stat().st_size < 1_000_000
    if point_format in (0, 2):
        assert main([str(path), "--debug", "-o", str(output)]) == 3
        _assert_refusal(output, "needs GPS time")
        return
    assert main([str(path), "--debug", "-o", str(output)]) == 0
    rows = np.loadtxt(output, delimiter=",", comments=("#", "time,"), usecols=(0, 1, 2, 3))
    assert len(rows) == 12
    assert np.max(np.abs(rows[:, 1:] - np.column_stack([
        np.interp(rows[:, 0], time, coord) for coord in sensor.T
    ]))) < 0.05
    report = _report(output)
    assert report["crs"] == "local"
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
    assert main(["--debug", "-o", str(output), str(first), str(second)]) == 0
    rows = np.loadtxt(output, delimiter=",", comments=("#", "time,"), usecols=(0, 1, 2, 3))
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
    ("truncated_laz", "failed to fill whole buffer"),
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
        case "truncated_laz":
            input_path = input_path.with_suffix(".laz")
            paths = [str(input_path)]
            flight(input_path, seconds=2, rate=400)
            data = input_path.read_bytes()
            input_path.write_bytes(data[:len(data) // 2])
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
    expected = {"constant_time": 3, "nonfinite_time": 3, "no_fit": 4}.get(scenario, 2)
    assert main([*paths, "--debug", "-o", str(output)]) == expected
    report = _assert_refusal(output, reason)
    assert report["exit_code"] == expected
    if expected in (3, 4):
        assert report["error"].endswith(
            "Keep the delivered path as supplementary data on the dataset; it will be added once supported."
        )


@pytest.mark.parametrize(("option", "reason"), [
    ("--interval=inf", "positive finite"),
    ("--interval=0", "positive finite"),
    ("--interval=abc", "positive finite"),
    ("--threads=0", "positive integer"),
    ("--threads=2.5", "positive integer"),
    ("--dtm=terrain.tif", "--dtm requires --occlusion-copc and --occlusion-cog"),
])
def test_invalid_arguments_have_a_refusal(tmp_path, option, reason):
    output = tmp_path / "output.csv"
    assert main(["input.las", "--debug", "-o", str(output), option]) == 2
    assert _assert_refusal(output, reason)["resources"]["points_read"] == 0


def test_outputs_cannot_alias_the_input(tmp_path):
    input_path = tmp_path / "input.las"
    input_path.write_bytes(b"unchanged")
    output = tmp_path / "alias.las"
    output.symlink_to(input_path)
    assert main([str(input_path), "--debug", "-o", str(output)]) == 2
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
        subprocess.run([sys.executable, str(entrypoint), str(path), "--debug", "-o", str(output), "--threads", "16"],
                       env=env, check=True)
        results.append((output.read_bytes(), output.with_suffix(".geojson").read_bytes()))
        assert 1 <= _report(output)["resources"]["decoding_threads"] <= int(slots)
    assert results[0] == results[1]


def test_one_position_skips_the_map_layer(tmp_path):
    path = tmp_path / "input.laz"
    flight(path, seconds=6, rate=1600, single_return=True)
    output = tmp_path / "output.csv"
    assert main([str(path), "--debug", "-o", str(output), "--interval", "10"]) == 0
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


def test_returns_stored_apart_with_missing_middle_returns_still_give_lines(tmp_path):
    """LiDAR360-style export: first and last return of a shot stored apart, return 2 of 3 missing, many lone
    last returns. The shots with first + last still give the lines method."""
    from synthetic import _records, _hits, _write
    rng = np.random.default_rng(3)
    time = 300000.0 + np.arange(0, 8, 1 / 400)
    sensor = np.c_[np.full(len(time), 1000.0), 2000.0 + 5.0 * (time - time[0]), np.full(len(time), 60.0)]
    angle = rng.uniform(-40, 40, len(time))
    a = np.radians(angle)
    direction = np.c_[np.sin(a), np.zeros(len(a)), -np.cos(a)]
    first, last, two = _hits(sensor, direction, rng, 1.0)
    keep = two & (rng.random(len(time)) < 0.3)  # 30 % of shots keep their first return
    xyz = np.r_[last, first[keep]]  # all last returns first, then the first returns: stored apart
    t = np.r_[time, time[keep]]
    rn = np.r_[np.full(len(time), 3), np.ones(keep.sum())].astype(np.uint8)
    nr = np.full(len(t), 3, np.uint8)
    cloud = tmp_path / "apart.las"
    _write(cloud, xyz, t, rn, nr, np.zeros(len(t)))
    output = tmp_path / "path.csv"
    assert main([str(cloud), "-o", str(output), "--debug"]) == 0
    report = _report(output)
    assert report["method"] == "lines"
    assert report["pulse_check"]["return_modes"][0]["line_pulses"] > 0
    fits = np.array([[float(v) for v in line.split(",")[:4]] for line in output.read_text().splitlines()
                     if line and line[0].isdigit()])
    truth = np.c_[[np.interp(fits[:, 0], time, sensor[:, k]) for k in range(3)]].T
    assert np.median(np.linalg.norm(fits[:, 1:] - truth, axis=1)) < 0.05
