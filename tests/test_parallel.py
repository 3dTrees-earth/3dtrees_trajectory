import os
import signal
from concurrent.futures.process import BrokenProcessPool

import laspy
import numpy as np
import pytest

from legacy_occlusion import trace_reference
from synthetic import flight
from test_occlusion import needs_untwine
from threedtrees_trajectory import cloud as cloud_io, occlusion, pulses, tile_jobs
from threedtrees_trajectory.cli import main
from threedtrees_trajectory.cloud import TimeSortedCopc
from threedtrees_trajectory.pulses import read_pulses
from threedtrees_trajectory.resources import RunResources
from threedtrees_trajectory.tile_spool import ECHO, START, TileSpools, pack_echoes, tile_batches


def _worker_checkpoint_test(index, work, kill):
    if kill:
        os.kill(os.getpid(), signal.SIGKILL)
    values = np.full((2, 2, 2), index, np.int32)
    data = {"tile": index + 1}
    tile_jobs.save_tile(work, index, (values,) * 3, data)
    return data


def _specs(work, kill):
    return [(index, 1, (2, 2, 2), (work, kill and index == 1)) for index in range(2)]


def test_reuse_finished_tile_after_a_worker_is_killed(tmp_path):
    with pytest.raises(BrokenProcessPool):
        tile_jobs.execute(tmp_path, _specs(tmp_path, True), 1, _worker_checkpoint_test, None, (), {})
    first, _ = tile_jobs.tile_files(tmp_path, 0)
    committed = first.stat().st_mtime_ns
    stage = {}
    tiles = tile_jobs.execute(tmp_path, _specs(tmp_path, False), 2, _worker_checkpoint_test, None, (), stage)
    assert stage["reused_tiles"] == 1 and stage["scheduled_tiles"] == 1
    assert tiles == {0: {"tile": 1}, 1: {"tile": 2}}
    assert first.stat().st_mtime_ns == committed


@needs_untwine
def test_parallel_and_sequential_tiles_write_identical_products(tmp_path, monkeypatch):
    cloud = tmp_path / "flight.las"
    time, xyz = flight(cloud, seconds=2, rate=300)
    grid, _ = occlusion._grid([cloud], 1)
    nx, ny, nz = grid.shape
    monkeypatch.setattr(occlusion, "WORKING_BYTES", 0)
    budget = occlusion._grid_bytes((nx * 3 // 4, ny, nz), nx * ny)
    results = []
    for count in (1, 2):
        copc, cog = (tmp_path / f"{count}.{ext}" for ext in ("copc.laz", "tif"))
        report = occlusion.run([cloud], time, xyz, copc, cog, voxel=1, budget=budget, threads=count)
        assert report["tiling"]["tiles"] >= 2 and report["tiling"]["peak_concurrent_workers"] == count
        results.append((copc.read_bytes(), cog.read_bytes(), report["voxels_ground_to_canopy"]))
    assert results[0] == results[1]


@needs_untwine
def test_separated_strip_tile_spool_matches_original_order(tmp_path):
    cloud, separated = tmp_path / "flight.las", tmp_path / "separated.laz"
    time, xyz = flight(cloud, seconds=2, rate=300)
    data = laspy.read(cloud)
    data.points = data.points[np.argsort(np.asarray(data.return_number), kind="stable")]
    data.write(separated)
    results = []
    for path in (cloud, separated):
        copc, cog = (tmp_path / (path.stem + ext) for ext in (".copc.laz", ".tif"))
        report = occlusion.run([path], time, xyz, copc, cog, voxel=1)
        results.append((copc.read_bytes(), cog.read_bytes(), report))
    assert results[0][:2] == results[1][:2]
    assert results[1][2]["recovered_echoes"] > 100
    assert results[1][2]["unpaired_echoes"] == 0


@needs_untwine
def test_copc_sample_sort_retains_reference_path_inputs(tmp_path, monkeypatch):
    cloud, copc = tmp_path / "flight.las", tmp_path / "flight.copc.laz"
    flight(cloud, seconds=6, rate=500)
    with laspy.open(cloud) as src:
        occlusion._index_copc(cloud, copc, src.header.point_count, tmp_path)
    full = tmp_path / "full"
    full.mkdir()
    old_sort = TimeSortedCopc(copc, full, 250)
    monkeypatch.setattr(pulses, "MAX_SINGLE_ONLY", 1500)
    reference, _ = read_pulses([copc], max_multi=1500, max_single=250, prepared=old_sort)
    metrics = RunResources()
    actual, info = read_pulses([copc], max_multi=1500, max_single=250, scratch=tmp_path, resources=metrics)
    for field in ("time", "first", "last", "returns", "scan_angle", "file", "source"):
        np.testing.assert_array_equal(getattr(actual, field), getattr(reference, field))
    assert info.points_read == old_sort.points_read
    sample, = metrics.stages
    assert sample["points_read"] > sample["sampled_echoes"] > 0


@needs_untwine
@pytest.mark.parametrize("spatial", [False, True])
def test_cli_decodes_source_once_for_path_and_tiles(tmp_path, monkeypatch, spatial):
    original = tmp_path / "flight.laz"
    flight(original, seconds=4, rate=500)
    source = original
    if spatial:
        source = tmp_path / "flight.copc.laz"
        with laspy.open(original) as src:
            occlusion._index_copc(original, source, src.header.point_count, tmp_path)
    reads = []
    read = cloud_io.point_chunks

    def counted(*args, **kwargs):
        reads.append(kwargs.get("spatial", False))
        yield from read(*args, **kwargs)

    for module in (cloud_io, pulses, occlusion):
        monkeypatch.setattr(module, "point_chunks", counted)
    assert main([str(source), "-o", str(tmp_path / "path.csv"),
                 "--occlusion-copc", str(tmp_path / "out.copc.laz"),
                 "--occlusion-cog", str(tmp_path / "out.tif")]) == 0
    assert reads == [spatial]
    assert (tmp_path / "out.tif").stat().st_size


@pytest.mark.parametrize("compress", [False, True])
def test_spool_routes_whole_separated_pulse_to_both_tiles(tmp_path, compress):
    grid = occlusion.Grid(np.zeros(3), (10, 1, 4), 1)
    left, right = (occlusion.Tile((slice(start, stop), slice(0, 1)),
                                 (slice(start, stop), slice(0, 1), slice(0, 4)))
                   for start, stop in ((0, 5), (5, 10)))
    raw = tmp_path / "raw.bin"
    records = np.zeros(3, ECHO)
    records["time"], records["x"], records["returns"] = [11, 10, 11], [2, 8, 7], [33, 17, 34]
    records.tofile(raw)
    routing = TileSpools(tmp_path, grid, [left, right], np.array([10., 11.]), 3, compress=compress)
    assert routing.stage_file(0, raw, np.ones(3), np.zeros(3), False) == (2, 0)
    for index, expected in enumerate(([2, 7], [8, 2, 7])):
        batch, = tile_batches(tmp_path / f"tile-{index:05d}", tmp_path)
        assert batch["x"].tolist() == expected
        assert batch["flags"][-2:].tolist() == [START, 0]


@pytest.mark.parametrize("compress", [False, True])
def test_dense_tile_window_is_time_sorted_in_bounded_batches(tmp_path, monkeypatch, compress):
    from threedtrees_trajectory import tile_spool

    monkeypatch.setattr(tile_spool, "BATCH", 3)
    grid = occlusion.Grid(np.zeros(3), (2, 1, 1), 1)
    core = (slice(0, 2), slice(0, 1))
    tile = occlusion.Tile(core, (*core, slice(0, 1)))
    raw = tmp_path / "raw.bin"
    records = np.zeros(16, ECHO)
    records["time"] = np.r_[np.arange(8, 16), np.arange(8)]
    records["returns"] = 17
    records.tofile(raw)
    routing = TileSpools(tmp_path, grid, [tile], np.array([0., 15.]), 1, compress=compress)
    routing.stage_file(0, raw, np.ones(3), np.zeros(3), False)
    batches = list(tile_batches(tmp_path / "tile-00000", tmp_path))
    assert max(map(len, batches)) <= 6
    np.testing.assert_array_equal(np.concatenate(batches)["time"], np.arange(16))


def test_parallel_strip_spools_merge_each_window_in_time_order(tmp_path):
    grid = occlusion.Grid(np.zeros(3), (2, 1, 1), 1)
    core = (slice(0, 2), slice(0, 1))
    tile = occlusion.Tile(core, (*core, slice(0, 1)))
    paths = [tmp_path / f"strip-{i}.bin" for i in range(2)]
    for index, path in enumerate(paths):
        records = np.zeros(2, ECHO)
        records["time"] = [5.0 + index, 2.0 + index]
        records["returns"] = 17
        records.tofile(path)
    stage = {}
    recovered, unpaired, echoes = occlusion._route_file_spools(
        tmp_path, grid, [tile], np.array([0., 20_480.]), occlusion.NATIVE_ECHO_LIMIT + 1, paths,
        [(np.ones(3), np.zeros(3))] * 2, None, False, 2, stage)
    assert (recovered, unpaired, echoes.tolist()) == (0, 0, [4])
    assert stage["routing_workers"] == 2
    batch, = tile_batches(tmp_path / "tile-00000", tmp_path)
    np.testing.assert_array_equal(batch["time"], [2., 3., 5., 6.])
    np.testing.assert_array_equal(batch["file"], [0, 1, 0, 1])


def test_parallel_copc_partitions_join_echoes_before_tile_sort(tmp_path):
    grid = occlusion.Grid(np.zeros(3), (20, 1, 1), 1)
    core = (slice(0, 20), slice(0, 1))
    tile = occlusion.Tile(core, (*core, slice(0, 1)))
    path = tmp_path / "spatial.bin"
    records = np.zeros(6, ECHO)
    records["time"] = [12, 11, 10, 10, 11, 12]
    records["x"] = [3, 2, 1, 11, 12, 13]
    records["returns"] = [33, 33, 33, 34, 34, 34]
    records.tofile(path)
    stage = {}
    _, unpaired, echoes = occlusion._route_file_spools(
        tmp_path, grid, [tile], np.array([0., 20_480.]), occlusion.NATIVE_ECHO_LIMIT + 1, [path],
        [(np.ones(3), np.zeros(3))], None, True, 2, stage)
    assert (unpaired, echoes.tolist()) == (0, [6])
    assert stage["routing_workers"] == 2
    batch, = tile_batches(tmp_path / "tile-00000", tmp_path)
    np.testing.assert_array_equal(batch["time"], [10, 10, 11, 11, 12, 12])
    np.testing.assert_array_equal(batch["x"], [1, 11, 2, 12, 3, 13])


def test_tile_spool_counts_match_acquisition_order_traversal(tmp_path):
    cloud = tmp_path / "flight.las"
    time, xyz = flight(cloud, seconds=2, rate=300)
    grid, _ = occlusion._grid([cloud], 0.5)
    counters = dict(echoes_read=0, echoes_outside_path=0, traversed_pulses=0, missing_returns=0)
    expected = trace_reference([cloud], grid, time, xyz, counters)
    core, z = tuple(slice(0, size) for size in grid.shape[:2]), slice(0, grid.shape[2])
    tile = occlusion.Tile(core, (*core, z))
    points = laspy.read(cloud).points
    raw = tmp_path / "echoes.bin"
    pack_echoes(points).tofile(raw)
    routing = TileSpools(tmp_path, grid, [tile], time, len(points))
    routing.stage_file(0, raw, points.scales, points.offsets, False)
    occlusion._init_traversal(tmp_path, tmp_path, time, xyz, [(points.scales, points.offsets)])
    occlusion._trace_tile(0, grid, core, int(routing.echoes[0]))
    counts, _ = tile_jobs.tile_files(tmp_path, 0)
    for old, new in zip(expected, np.load(counts, mmap_mode="r", allow_pickle=False)):
        np.testing.assert_array_equal(old, new)
