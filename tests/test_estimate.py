import json

import numpy as np

from threedtrees_trajectory.cli import main

from synthetic import flight, scan


def _read(path):
    """Column names and rows: time, x, y, z, roll, pitch, yaw, name."""
    lines = [l for l in path.read_text().splitlines() if not l.startswith("#")]
    return lines[0].split(","), [l.split(",") for l in lines[1:]]


def test_flight_multi_return_recovers_path_to_centimetres(tmp_path):
    time, sensor = flight(tmp_path / "f.laz")
    assert main([str(tmp_path / "f.laz"), "-o", str(tmp_path / "p.csv")]) == 0
    columns, rows = _read(tmp_path / "p.csv")
    assert columns == ["time", "x", "y", "z", "roll", "pitch", "yaw", "name"]
    est = np.array([r[:4] for r in rows], float)
    assert all(not any(r[4:]) for r in rows)  # a path: no attitude, no names
    truth = np.c_[np.interp(est[:, 0], time, sensor[:, 0]), np.interp(est[:, 0], time, sensor[:, 1]), np.interp(est[:, 0], time, sensor[:, 2])]
    assert len(est) >= 35
    assert np.abs(est[:, 1:] - truth).max() < 0.05
    report = json.loads((tmp_path / "p.report.json").read_text())
    assert report["method"] == "lines" and "kind" not in report
    text = (tmp_path / "p.csv").read_text()
    assert "# quality: good" in text and "# kind: path" in text and "# origin: estimated" in text and "miss_m" not in text
    assert report["median_miss_m"] < 0.05 and report["typical_error_m"] == "0.01"


def test_flight_single_return_uses_scan_angle_and_warns(tmp_path):
    time, sensor = flight(tmp_path / "f.laz", single_return=True)
    assert main([str(tmp_path / "f.laz"), "-o", str(tmp_path / "p.csv")]) == 0
    est = np.array([r[:4] for r in _read(tmp_path / "p.csv")[1]], float)
    h = np.hypot(est[:, 1] - np.interp(est[:, 0], time, sensor[:, 0]), est[:, 2] - np.interp(est[:, 0], time, sensor[:, 1]))
    assert np.median(h) < 1.0
    assert np.median(np.abs(est[:, 3] - 60.0)) < 0.5
    report = json.loads((tmp_path / "p.report.json").read_text())
    assert report["method"] == "scan_angle" and report["warnings"]


def test_static_scans_are_refused(tmp_path, capsys):
    """1.0.0 estimates paths only: a TLS scan (multi-return lines meeting in one point) is refused, with a report."""
    scan(tmp_path / "a.laz", (3.0, 4.0, 1.5), seed=1)
    scan(tmp_path / "b.laz", (25.0, -7.0, 1.6), seed=2)
    assert main([str(tmp_path / "a.laz"), str(tmp_path / "b.laz"), "-o", str(tmp_path / "s.csv")]) == 2
    assert "static (TLS)" in capsys.readouterr().err
    assert (tmp_path / "s.csv").stat().st_size == 0
    assert (tmp_path / "s.geojson").stat().st_size == 0
    assert json.loads((tmp_path / "s.report.json").read_text())["quality"] == "none"


def test_single_return_without_scan_angle_is_refused(tmp_path, capsys):
    import laspy

    flight(tmp_path / "f.laz", single_return=True)
    las = laspy.read(tmp_path / "f.laz")
    las.scan_angle = np.zeros(len(las.points), np.int16)
    las.write(tmp_path / "g.laz")
    assert main([str(tmp_path / "g.laz"), "-o", str(tmp_path / "p.csv"), "--no-platform-check"]) == 4
    assert "trajectory file" in capsys.readouterr().err


def test_platform_is_detected_for_single_returns(tmp_path, capsys):
    scan(tmp_path / "s.laz", (3.0, 4.0, 1.5), single_return=True, seed=1, rays=200000)
    flight(tmp_path / "f.laz", single_return=True)
    assert main([str(tmp_path / "s.laz"), "-o", str(tmp_path / "s.csv")]) == 2
    assert "ground density" in capsys.readouterr().err
    assert main([str(tmp_path / "f.laz"), "-o", str(tmp_path / "f.csv")]) == 0
    assert json.loads((tmp_path / "f.report.json").read_text())["method"] == "scan_angle"


def test_geojson_map_layer_splits_at_gaps(tmp_path):
    flight(tmp_path / "f.laz")
    assert main([str(tmp_path / "f.laz"), "-o", str(tmp_path / "p.csv")]) == 0
    layer = json.loads((tmp_path / "p.geojson").read_text())
    feature = layer["features"][0]
    assert feature["geometry"]["type"] == "MultiLineString"
    assert feature["properties"]["estimated"] is True and feature["properties"]["quality"] == "good"
    assert sum(len(line) for line in feature["geometry"]["coordinates"]) <= 110  # 100 m flown, about 1 vertex per metre
