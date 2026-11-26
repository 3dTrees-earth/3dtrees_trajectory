import json

import laspy
import numpy as np
import pytest

from synthetic import flight, scan
from threedtrees_trajectory import cli, full_path
from threedtrees_trajectory.cli import main


def _run_path(clouds, out):
    assert main([*map(str, clouds), "--full", "--threads", "1", "--debug", "-o", str(out)]) == 0
    return json.loads(out.with_suffix(".report.json").read_text())


@pytest.mark.parametrize("separated", [False, True])
def test_full_windows_match_whole_cloud(tmp_path, monkeypatch, separated):
    cloud = tmp_path / "flight.las"
    flight(cloud, seconds=12, rate=500)
    if separated:
        las = laspy.read(cloud)
        order = np.argsort(np.asarray(las.return_number), kind="stable")
        las.points = las.points[order]
        # All second returns occur after the first pass through the GPS clock.
        las.write(cloud)

    monkeypatch.setattr(full_path, "CHUNK", 3601)
    expected, actual = tmp_path / "whole.csv", tmp_path / "streamed.csv"
    with monkeypatch.context() as patch:
        patch.setattr(full_path, "try_estimate", lambda *args: None)
        patch.setattr(cli, "FULL_IN_MEMORY", True)
        reference = _run_path([cloud], expected)
    report = _run_path([cloud], actual)
    assert report["stages"][0]["full_window_fit"] == "streamed"
    assert report["pulse_check"] == reference["pulse_check"]
    assert actual.read_bytes() == expected.read_bytes()
    assert not list(tmp_path.glob(".path-windows-*"))


def test_full_merges_overlapping_strips_before_fitting(tmp_path, monkeypatch):
    clouds = [tmp_path / name for name in ("a.las", "b.las")]
    for i, cloud in enumerate(clouds):
        flight(cloud, seconds=8, rate=500, seed=i)
    second = laspy.read(clouds[1])
    second.points = second.points[np.argsort(np.asarray(second.return_number), kind="stable")]
    second.write(clouds[1])

    monkeypatch.setattr(full_path, "CHUNK", 3001)
    expected, actual = tmp_path / "whole.csv", tmp_path / "streamed.csv"
    with monkeypatch.context() as patch:
        patch.setattr(full_path, "try_estimate", lambda *args: None)
        patch.setattr(cli, "FULL_IN_MEMORY", True)
        reference = _run_path(clouds, expected)
    report = _run_path(clouds, actual)
    assert report["stages"][0]["full_window_fit"] == "streamed"
    assert report["pulse_check"] == reference["pulse_check"]
    assert actual.read_bytes() == expected.read_bytes()


def test_full_still_refuses_static_cloud(tmp_path):
    cloud = tmp_path / "static.las"
    scan(cloud, (1000, 2000, 30), rays=35_000)
    out = tmp_path / "path.csv"
    assert main([str(cloud), "--full", "--debug", "-o", str(out)]) == 2
    assert not out.stat().st_size
    assert "static (TLS) scan" in out.with_suffix(".report.json").read_text()


def test_all_pulses_is_the_default_and_sample_opts_out(tmp_path):
    cloud = tmp_path / "flight.las"
    flight(cloud, seconds=6, rate=400)
    full, sample = tmp_path / "full.csv", tmp_path / "sample.csv"
    assert main([str(cloud), "--threads", "1", "--debug", "-o", str(full)]) == 0
    assert main([str(cloud), "--threads", "1", "--debug", "-o", str(sample), "--sample"]) == 0
    full_report = json.loads(full.with_suffix(".report.json").read_text())
    sample_report = json.loads(sample.with_suffix(".report.json").read_text())
    assert full_report["stages"][0]["full_window_fit"] == "streamed"
    assert "full_window_fit" not in sample_report["stages"][0]
    assert main([str(cloud), "--full", "--threads", "1", "-o", str(tmp_path / "alias.csv")]) == 0  # old flag still accepted
    assert (tmp_path / "alias.csv").read_bytes() == full.read_bytes()
