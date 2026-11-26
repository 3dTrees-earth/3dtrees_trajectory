"""CRS labels (``local`` for cloud coordinates) and delivered paths converted to the cloud's CRS."""
import json

import laspy
import numpy as np
import pytest
from pyproj import CRS, Transformer

from synthetic import flight
from threedtrees_trajectory import occlusion
from threedtrees_trajectory.cli import main
from threedtrees_trajectory.pulses import crs_labels

SHIFT = np.array([500_000.0, 5_300_000.0, 0.0])  # the synthetic flight placed in UTM 32N


def _utm_flight(tmp_path, epsg=32632):
    cloud = tmp_path / "flight.las"
    time, sensor = flight(cloud, seconds=4, rate=500)
    points = laspy.read(cloud)
    points.header.offsets += SHIFT
    points.points.offsets = points.header.offsets.copy()
    if epsg:
        points.header.add_crs(CRS.from_epsg(epsg))
    points.write(cloud)
    return cloud, time, sensor + (SHIFT if epsg else 0)


def _path(tmp_path, time, xyz, crs_line):
    path = tmp_path / "delivered.csv"
    with path.open("w") as stream:
        stream.write("# format: 3dtrees-trajectory/1\n" + (f"# crs: {crs_line}\n" if crs_line else "") + "time,x,y,z\n")
        np.savetxt(stream, np.column_stack((time[::25], xyz[::25])), delimiter=",", fmt="%.9f")
    return path


def _run(tmp_path, cloud, path):
    report = tmp_path / "report.json"
    code = main([str(cloud), "--trajectory", str(path), "--debug", "--report", str(report),
                 "--occlusion-copc", str(tmp_path / "out.copc.laz"), "--occlusion-cog", str(tmp_path / "out.tif")])
    return code, json.loads(report.read_text())["occlusion"]


def test_crs_labels():
    assert crs_labels(None) == ("local", "unknown")
    assert crs_labels(CRS.from_epsg(32632)) == ("EPSG:32632", "unknown")
    assert crs_labels(CRS("EPSG:25832+EPSG:7837")) == ("EPSG:25832", "EPSG:7837")


@pytest.mark.parametrize("target", ["EPSG:4326", "EPSG:32633"])
def test_delivered_path_in_another_crs_is_converted(tmp_path, target):
    cloud, time, sensor = _utm_flight(tmp_path)
    x, y = Transformer.from_crs("EPSG:32632", target, always_xy=True).transform(sensor[:, 0], sensor[:, 1])
    code, section = _run(tmp_path, cloud, _path(tmp_path, time, np.c_[x, y, sensor[:, 2]], target))
    check = section["delivered_path_check"]
    assert check["crs_converted"]["from"] == target and check["crs_converted"]["to"] == "EPSG:32632"
    assert check["median_m"] < 0.01  # back on the pulse lines after the conversion
    assert f"CRS ({target}) was converted to the cloud's CRS (EPSG:32632)" in check["crs_warning"]


def test_local_path_and_local_cloud_are_used_as_they_are(tmp_path):
    cloud, time, sensor = _utm_flight(tmp_path, epsg=None)
    code, section = _run(tmp_path, cloud, _path(tmp_path, time, sensor + SHIFT, "local"))
    check = section["delivered_path_check"]
    assert "crs_converted" not in check and "crs_warning" not in check and check["median_m"] < 0.01


@pytest.mark.parametrize("legacy", ["local", "unknown"])
def test_local_path_on_an_epsg_cloud_is_refused(tmp_path, legacy):
    cloud, time, sensor = _utm_flight(tmp_path)
    code, section = _run(tmp_path, cloud, _path(tmp_path, time, sensor, legacy))
    assert "cannot be aligned" in section["error"] and "supplementary" in section["error"]
    assert (tmp_path / "out.tif").read_bytes() == b""


def test_path_without_crs_line_assumes_the_cloud_crs(tmp_path):
    import rasterio

    cloud, time, sensor = _utm_flight(tmp_path)
    code, section = _run(tmp_path, cloud, _path(tmp_path, time, sensor, None))
    assert code == 0 and section["delivered_path_check"]["median_m"] < 0.01
    with rasterio.open(tmp_path / "out.tif") as raster:
        assert "declares no CRS; the cloud's CRS (EPSG:32632) was assumed" in raster.tags()["OCCLUSION_WARNINGS"]


def test_epsg_path_on_a_local_cloud_is_refused():
    with pytest.raises(ValueError, match="cannot be aligned"):
        occlusion.align_path_crs(np.zeros((2, 3)), "EPSG:32632", "unknown", "local", "unknown")


def test_local_cloud_writes_local(tmp_path):
    cloud = tmp_path / "flight.las"
    flight(cloud, seconds=6, rate=400)
    output = tmp_path / "path.csv"
    assert main([str(cloud), "-o", str(output)]) == 0
    assert "# crs: local\n" in output.read_text()
    assert not output.with_suffix(".geojson").exists()
