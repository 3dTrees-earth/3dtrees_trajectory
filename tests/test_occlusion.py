import json
import os
from pathlib import Path

import laspy
import numpy as np
import pytest

from synthetic import flight
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
    ground, spikes = occlusion.ground_from_hits(filled, grid)
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
def test_estimated_path_gives_copc_and_cog(tmp_path):
    import rasterio

    flight(tmp_path / "f.laz")
    out, copc, cog = tmp_path / "p.csv", tmp_path / "o.copc.laz", tmp_path / "o.tif"
    assert main([str(tmp_path / "f.laz"), "-o", str(out), "--occlusion-copc", str(copc), "--occlusion-cog", str(cog),
                 "--voxel", "1"]) == 0
    section = json.loads(out.with_suffix(".report.json").read_text())["occlusion"]
    assert section["trajectory"] == "estimated" and section["ground"] == "lowest_hit"
    with laspy.open(copc) as reader:
        assert reader.header.point_count == section["copc_points"] > 0
        assert any(type(vlr).__name__ == "CopcInfoVlr" for vlr in reader.header.vlrs)
        points = reader.read()
    assert set(np.unique(points["occlusion"])) <= {1, 3, 4}
    assert np.all(points["height_above_ground"] >= 0)
    with rasterio.open(cog) as raster:
        assert raster.count == len(section["cog_bands"])
        assert raster.descriptions[0] == "occluded_0-1m" and raster.descriptions[-1] == "unobserved_column"
        assert raster.res == (1.0, 1.0)
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
