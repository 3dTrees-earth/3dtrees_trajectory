"""Several delivered paths: each pulse segment (file, point source, contiguous GPS time) uses the path that covers it
and best explains its return lines."""
import json

import laspy
import numpy as np
import pytest

from synthetic import flight
from threedtrees_trajectory.cli import main

needs_untwine = pytest.mark.skipif(
    not __import__("shutil").which("untwine") and not __import__("os").environ.get("TRAJECTORY_UNTWINE"),
    reason="untwine not installed")


def _path(path, time, xyz, step=10):
    with path.open("w") as stream:
        stream.write("# format: 3dtrees-trajectory/1\ntime,x,y,z\n")
        np.savetxt(stream, np.column_stack((time[::step], xyz[::step])), delimiter=",", fmt="%.6f")
    return path


def _split(cloud, at, first, second):
    """Write the points of ``cloud`` before / from GPS time ``at`` into two files."""
    las = laspy.read(cloud)
    for target, keep in ((first, las.gps_time < at), (second, las.gps_time >= at)):
        part = laspy.LasData(las.header)
        part.points = las.points[keep]
        part.write(target)


def _run(tmp_path, clouds, paths):
    report = tmp_path / "report.json"
    csv = tmp_path / "path.csv"
    code = main([*map(str, clouds), "--trajectory", *map(str, paths), "--debug", "--report", str(report),
                 "-o", str(csv), "--occlusion-copc", str(tmp_path / "o.copc.laz"),
                 "--occlusion-cog", str(tmp_path / "o.tif")])
    data = json.loads(report.read_text())
    return code, data["occlusion"], csv.read_text()


@needs_untwine
def test_one_file_two_paths_split_in_time(tmp_path):
    cloud = tmp_path / "f.laz"
    time, sensor = flight(cloud, seconds=8, rate=500)
    mid = len(time) // 2
    a = _path(tmp_path / "a.csv", time[:mid], sensor[:mid])
    b = _path(tmp_path / "b.csv", time[mid:], sensor[mid:])
    code, section, text = _run(tmp_path, [cloud], [a, b])
    check = section["delivered_path_check"]
    assert code == 0 and "error" not in section
    assert check["coverage"] == 1.0 and check["path_files"] == 2
    assert {s["path"] for s in check["segments"]} == {"a.csv", "b.csv"}
    assert check["median_m"] < 0.01
    assert "# delivered_path_files: 2" in text and "# delivered_path_coverage: 1.0" in text


@needs_untwine
def test_two_files_one_path(tmp_path):
    whole = tmp_path / "whole.laz"
    time, sensor = flight(whole, seconds=8, rate=500)
    _split(whole, time[len(time) // 2], tmp_path / "f1.laz", tmp_path / "f2.laz")
    path = _path(tmp_path / "p.csv", time, sensor)
    code, section, text = _run(tmp_path, [tmp_path / "f1.laz", tmp_path / "f2.laz"], [path])
    check = section["delivered_path_check"]
    assert code == 0 and check["coverage"] == 1.0 and check["path_files"] == 1
    assert "segments" not in check  # one path over everything: traced as before
    assert check["median_m"] < 0.01


@needs_untwine
def test_two_files_two_paths(tmp_path):
    whole = tmp_path / "whole.laz"
    time, sensor = flight(whole, seconds=8, rate=500)
    mid = len(time) // 2
    _split(whole, time[mid], tmp_path / "f1.laz", tmp_path / "f2.laz")
    a = _path(tmp_path / "a.csv", time[:mid], sensor[:mid])
    b = _path(tmp_path / "b.csv", time[mid:], sensor[mid:])
    code, section, _ = _run(tmp_path, [tmp_path / "f1.laz", tmp_path / "f2.laz"], [b, a])
    check = section["delivered_path_check"]
    assert code == 0 and check["coverage"] == 1.0
    assert {(s["file"], s["path"]) for s in check["segments"]} == {(0, "a.csv"), (1, "b.csv")}


@needs_untwine
def test_overlapping_paths_the_one_that_fits_wins(tmp_path):
    cloud = tmp_path / "f.laz"
    time, sensor = flight(cloud, seconds=6, rate=500)
    right = _path(tmp_path / "right.csv", time, sensor)
    off = _path(tmp_path / "off.csv", time, sensor + [5.0, 0, 0])
    code, section, _ = _run(tmp_path, [cloud], [off, right])  # the wrong one first: order must not decide
    check = section["delivered_path_check"]
    assert code == 0 and {s["path"] for s in check["segments"]} == {"right.csv"}
    miss = check["segments"][0]["miss_m"]
    assert miss["right.csv"] < 0.01 < 1 < miss["off.csv"]
    assert check["median_m"] < 0.01


@needs_untwine
def test_week_seconds_collision_separated_by_point_source(tmp_path):
    """Two days in one file with the same GPS week seconds, told apart by point source and position."""
    one, two = tmp_path / "day1.laz", tmp_path / "day2.laz"
    time, sensor = flight(one, seconds=6, rate=500, seed=0)
    flight(two, seconds=6, rate=500, seed=1)
    first, second = laspy.read(one), laspy.read(two)
    first.point_source_id[:] = 1
    second.point_source_id[:] = 2
    merged = laspy.LasData(first.header)
    merged.header.offsets = first.header.offsets
    shifted = np.c_[second.x + 400.0, second.y, second.z]
    xyz = np.r_[np.c_[first.x, first.y, first.z], shifted]
    merged.points = laspy.ScaleAwarePointRecord.zeros(len(xyz), header=merged.header)
    merged.x, merged.y, merged.z = xyz.T
    for name in ("gps_time", "return_number", "number_of_returns", "scan_angle", "point_source_id"):
        setattr(merged, name, np.r_[np.asarray(first[name]), np.asarray(second[name])])
    cloud = tmp_path / "both.laz"
    merged.write(cloud)
    p1 = _path(tmp_path / "day1.csv", time, sensor)
    p2 = _path(tmp_path / "day2.csv", time, sensor + [400.0, 0, 0])
    code, section, _ = _run(tmp_path, [cloud], [p2, p1])
    check = section["delivered_path_check"]
    assert code == 0 and check["coverage"] == 1.0
    assert {(s["point_source"], s["path"]) for s in check["segments"]} == {(1, "day1.csv"), (2, "day2.csv")}
    assert check["median_m"] < 0.01


@needs_untwine
def test_duplicate_path_is_used_once(tmp_path):
    cloud = tmp_path / "f.laz"
    time, sensor = flight(cloud, seconds=6, rate=500)
    a = _path(tmp_path / "a.csv", time, sensor)
    b = _path(tmp_path / "b.csv", time, sensor)
    code, section, text = _run(tmp_path, [cloud], [a, b])
    check = section["delivered_path_check"]
    assert code == 0 and check["files"][1]["duplicate_of"] == "a.csv"
    assert {s["path"] for s in check["segments"]} == {"a.csv"}
    assert "b.csv repeats a.csv" in text


@needs_untwine
@pytest.mark.parametrize(("share", "refused"), [(0.7, False), (0.3, True)])
def test_uncovered_pulses_are_skipped_or_refused(tmp_path, share, refused):
    cloud = tmp_path / "f.laz"
    time, sensor = flight(cloud, seconds=10, rate=400)
    keep = time <= time[0] + share * (time[-1] - time[0])
    path = _path(tmp_path / "p.csv", time[keep], sensor[keep])
    code, section, text = _run(tmp_path, [cloud], [path])
    assert code == 0
    if refused:
        assert "had a delivered path covering their time" in section["error"]
        assert "supplementary" in section["error"]
        assert (tmp_path / "o.tif").read_bytes() == b""
    else:
        check = section["delivered_path_check"]
        assert 0.6 < check["coverage"] < 0.8
        assert "of pulses had no delivered path covering their time" in text
        assert (tmp_path / "o.tif").stat().st_size > 0
