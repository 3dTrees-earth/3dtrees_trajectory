import json

import numpy as np

from threedtrees_trajectory.cli import main

from synthetic import flight, scan


def _read(path):
    """Column names and rows: time, x, y, z, speed_m_s, scan."""
    lines = [l for l in path.read_text().splitlines() if not l.startswith("#")]
    return lines[0].split(","), [l.split(",") for l in lines[1:]]


def test_flight_multi_return_recovers_path_to_centimetres(tmp_path):
    time, sensor = flight(tmp_path / "f.laz")
    assert main([str(tmp_path / "f.laz"), "-o", str(tmp_path / "p.csv")]) == 0
    columns, rows = _read(tmp_path / "p.csv")
    assert columns == ["time", "x", "y", "z", "speed_m_s", "scan"]
    est = np.array([r[:4] for r in rows], float)
    assert all(r[4] and not r[5] for r in rows)  # a path has speeds and no scan names
    truth = np.c_[np.interp(est[:, 0], time, sensor[:, 0]), np.interp(est[:, 0], time, sensor[:, 1]), np.interp(est[:, 0], time, sensor[:, 2])]
    assert len(est) >= 35
    assert np.abs(est[:, 1:] - truth).max() < 0.05
    report = json.loads((tmp_path / "p.report.json").read_text())
    assert report["platform"] == "mobile" and report["method"] == "lines"
    text = (tmp_path / "p.csv").read_text()
    assert "# quality: good" in text and "# platform: mobile" in text and "miss_m:" in text


def test_flight_single_return_uses_scan_angle_and_warns(tmp_path):
    time, sensor = flight(tmp_path / "f.laz", single_return=True)
    assert main([str(tmp_path / "f.laz"), "-o", str(tmp_path / "p.csv"), "--platform", "mobile"]) == 0
    est = np.array([r[:4] for r in _read(tmp_path / "p.csv")[1]], float)
    h = np.hypot(est[:, 1] - np.interp(est[:, 0], time, sensor[:, 0]), est[:, 2] - np.interp(est[:, 0], time, sensor[:, 1]))
    assert np.median(h) < 1.0
    assert np.median(np.abs(est[:, 3] - 60.0)) < 0.5
    report = json.loads((tmp_path / "p.report.json").read_text())
    assert report["method"] == "scan_angle" and report["warnings"]


def test_tls_one_file_per_scan_gives_each_station(tmp_path):
    scan(tmp_path / "a.laz", (3.0, 4.0, 1.5), seed=1)
    scan(tmp_path / "b.laz", (25.0, -7.0, 1.6), seed=2)
    assert main([str(tmp_path / "a.laz"), str(tmp_path / "b.laz"), "-o", str(tmp_path / "s.csv")]) == 0
    columns, rows = _read(tmp_path / "s.csv")
    assert columns == ["time", "x", "y", "z", "speed_m_s", "scan"]
    xyz = np.array([r[1:4] for r in rows], float)
    assert [r[5] for r in rows] == ["scan0", "scan1"] and not any(r[4] for r in rows)  # stations: names, no speed
    assert np.abs(xyz - [[3.0, 4.0, 1.5], [25.0, -7.0, 1.6]]).max() < 0.02


def test_tls_merged_cloud_is_split_at_time_jumps(tmp_path):
    import laspy

    scan(tmp_path / "a.laz", (0.0, 0.0, 1.5), seed=1, t0=100000.0)
    scan(tmp_path / "b.laz", (20.0, 0.0, 1.5), seed=2, t0=50000.0)  # earlier time: merged files jump backwards
    a, b = laspy.read(tmp_path / "a.laz"), laspy.read(tmp_path / "b.laz")
    merged = laspy.LasData(laspy.LasHeader(point_format=6, version="1.4"))
    merged.header.offsets, merged.header.scales = [-1000.0, -1000.0, -10.0], [0.001, 0.001, 0.001]
    for name in ("x", "y", "z", "gps_time", "return_number", "number_of_returns"):
        setattr(merged, name, np.concatenate([np.asarray(getattr(a, name)), np.asarray(getattr(b, name))]))
    merged.write(tmp_path / "m.laz")
    assert main([str(tmp_path / "m.laz"), "-o", str(tmp_path / "s.csv"), "--platform", "static"]) == 0
    xyz = np.array([r[1:4] for r in _read(tmp_path / "s.csv")[1]], float)
    assert np.abs(xyz - [[0.0, 0.0, 1.5], [20.0, 0.0, 1.5]]).max() < 0.02


def test_tls_single_return_falls_back_to_density(tmp_path):
    scan(tmp_path / "a.laz", (3.0, 4.0, 1.5), single_return=True, seed=1, rays=200000)
    assert main([str(tmp_path / "a.laz"), "-o", str(tmp_path / "s.csv"), "--platform", "static"]) == 0
    xyz = np.array([r[1:4] for r in _read(tmp_path / "s.csv")[1]], float)
    assert np.hypot(xyz[0, 0] - 3.0, xyz[0, 1] - 4.0) < 1.0
    report = json.loads((tmp_path / "s.report.json").read_text())
    assert report["method"] == "density"


def test_single_return_without_scan_angle_is_refused(tmp_path, capsys):
    import laspy

    flight(tmp_path / "f.laz", single_return=True)
    las = laspy.read(tmp_path / "f.laz")
    las.scan_angle = np.zeros(len(las.points), np.int16)
    las.write(tmp_path / "g.laz")
    assert main([str(tmp_path / "g.laz"), "-o", str(tmp_path / "p.csv"), "--platform", "mobile"]) == 2
    assert "trajectory file" in capsys.readouterr().err


def test_platform_is_detected_for_single_returns(tmp_path):
    scan(tmp_path / "s.laz", (3.0, 4.0, 1.5), single_return=True, seed=1, rays=200000)
    flight(tmp_path / "f.laz", single_return=True)
    assert main([str(tmp_path / "s.laz"), "-o", str(tmp_path / "s.csv")]) == 0
    assert json.loads((tmp_path / "s.report.json").read_text())["platform"] == "static"
    assert main([str(tmp_path / "f.laz"), "-o", str(tmp_path / "f.csv")]) == 0
    assert json.loads((tmp_path / "f.report.json").read_text())["platform"] == "mobile"


def test_geojson_map_layer_splits_at_gaps(tmp_path):
    flight(tmp_path / "f.laz")
    assert main([str(tmp_path / "f.laz"), "-o", str(tmp_path / "p.csv")]) == 0
    layer = json.loads((tmp_path / "p.geojson").read_text())
    feature = layer["features"][0]
    assert feature["geometry"]["type"] == "MultiLineString"
    assert feature["properties"]["estimated"] is True and feature["properties"]["quality"] == "good"
    assert sum(len(line) for line in feature["geometry"]["coordinates"]) <= 110  # 100 m flown, about 1 vertex per metre
