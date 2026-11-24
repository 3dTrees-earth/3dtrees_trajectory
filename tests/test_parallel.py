import os
import signal
from errno import ENOMEM, ENOSPC
from types import SimpleNamespace

import laspy
import numpy as np
import pytest

from legacy_occlusion import trace_reference
from synthetic import flight
from test_occlusion import needs_untwine
from threedtrees_trajectory import cloud as cloud_io, occlusion, parallel_read, pulses, time_windows, window_jobs
from threedtrees_trajectory.cli import main
from threedtrees_trajectory.cloud import TimeSortedCopc
from threedtrees_trajectory.pulses import read_pulses
from threedtrees_trajectory.resources import OutOfMemory, RunResources
from threedtrees_trajectory.window_input import ECHO, ScanIndex, group_multi


def _worker_checkpoint_test(index, work, failure, peak):
    if failure == "kill":
        os.kill(os.getpid(), signal.SIGKILL)
    if failure == "invalid":
        raise ValueError("invalid input")
    if failure == "full":
        raise OSError(ENOSPC, "No space left on device")
    if failure == "sigbus":
        os.kill(os.getpid(), signal.SIGBUS)
    if failure == "numpy":
        np.empty(2**60, dtype=np.uint8)
    if failure == "mmap":
        raise OSError(ENOMEM, "Cannot allocate memory")
    counts, marker = window_jobs.window_files(work, index)
    np.save(counts, np.full((8, 3), index, np.uint32))
    data = {"window": index + 1, "peak_rss_bytes": peak}
    window_jobs.atomic_json(marker, data)
    return data


def _spec(index, work, failure=None, *, planned=2**30, peak=2**28):
    return index, planned, (8, 3), (work, failure, peak)


def test_reuse_finished_window_after_a_worker_is_killed(tmp_path, monkeypatch):
    memory = 5 * 2**30
    monkeypatch.setattr(window_jobs, "available_memory", lambda: memory)
    stage = {}
    with pytest.raises(OutOfMemory, match=r"Occlusion window 2 was killed \(SIGKILL\)"):
        window_jobs.execute(tmp_path, [_spec(0, tmp_path), _spec(1, tmp_path, "kill")], 1,
                            _worker_checkpoint_test, None, (), stage)
    assert stage["completed_windows"] == 1
    first, _ = window_jobs.window_files(tmp_path, 0)
    committed = first.stat().st_mtime_ns
    stage = {}
    windows = window_jobs.execute(tmp_path, [_spec(0, tmp_path), _spec(1, tmp_path)], 2,
                                  _worker_checkpoint_test, None, (), stage)
    assert stage["reused_windows"] == 1 and stage["scheduled_windows"] == 1
    assert [result["window"] for result in windows.values()] == [1, 2]
    assert first.stat().st_mtime_ns == committed


def test_peak_rss_admits_more_windows(tmp_path, monkeypatch):
    gib = 2**30
    memory = int(2.75 * gib)
    monkeypatch.setattr(window_jobs, "available_memory", lambda: memory)
    specs = [_spec(0, tmp_path, planned=2 * gib), _spec(1, tmp_path, planned=2 * gib),
             _spec(2, tmp_path, planned=2 * gib), _spec(3, tmp_path, planned=2 * gib)]
    stage = {}
    windows = window_jobs.execute(tmp_path, specs, 3, _worker_checkpoint_test, None, (), stage)
    assert sorted(windows) == list(range(4))
    reservations = stage["window_reservations_bytes"]
    assert reservations[1] == [2 * gib]
    assert reservations[2][0] < gib
    assert stage["peak_concurrent_workers"] == 3
    assert stage["peak_reserved_bytes"] == 2 * gib


def test_invalid_input_from_a_worker_remains_a_value_error(tmp_path):
    stage = {}
    with pytest.raises(ValueError, match="invalid input"):
        window_jobs.execute(tmp_path, [_spec(0, tmp_path, "invalid")], 1, _worker_checkpoint_test, None, (), stage)
    assert stage["completed_windows"] == 0


@pytest.mark.parametrize("failure", ["numpy", "mmap"])
def test_worker_allocation_failure_is_oom(tmp_path, failure, monkeypatch):
    monkeypatch.setattr(window_jobs, "available_memory", lambda: 2**30)
    with pytest.raises(OutOfMemory, match="Occlusion window 1"):
        window_jobs.execute(tmp_path, [_spec(0, tmp_path, failure)], 1, _worker_checkpoint_test, None, (), {})


def test_memory_pressure_admits_one_worker_at_a_time(tmp_path, monkeypatch):
    gib = 2**30
    memory = iter((8 * gib, gib, gib))
    monkeypatch.setattr(window_jobs, "available_memory", lambda: next(memory))
    stage = {}
    windows = window_jobs.execute(tmp_path, [_spec(0, tmp_path), _spec(1, tmp_path)], 2,
                                  _worker_checkpoint_test, None, (), stage)
    assert list(windows) == [0, 1]
    assert stage["peak_concurrent_workers"] == 1


@pytest.mark.parametrize(("failure", "message"), [("full", "No space left"), ("sigbus", "Check free scratch space")])
def test_checkpoint_write_failure_is_clear(tmp_path, failure, message):
    stage = {}
    with pytest.raises(OSError, match=message):
        window_jobs.execute(tmp_path, [_spec(0, tmp_path, failure)], 1, _worker_checkpoint_test, None, (), stage)
    assert window_jobs.finished(tmp_path, 0, (8, 3)) is None


def test_fragment_batches_keep_group_boundaries(monkeypatch):
    monkeypatch.setattr(time_windows, "TRACE_BATCH", 17)
    sizes = [2] * 10 + [45] + [3] * 10
    original = np.zeros(sum(sizes), dtype=ECHO)
    original["time"] = np.repeat(np.arange(len(sizes)), sizes)
    original["returns"] = np.concatenate([
        np.arange(1, size + 1) | (size << 4) if size <= 15 else np.full(size, 0x21)
        for size in sizes
    ])
    shuffled = original[np.random.default_rng(0).permutation(len(original))]
    recovered = []
    for batch in time_windows._sorted_fragments(shuffled):
        counts, valid = group_multi(batch, acquisition_order=False)
        recovered.append(batch[np.repeat(valid, counts)])
    np.testing.assert_array_equal(np.concatenate(recovered), original[original["time"] != 10])


def test_product_blocks_reuse_column_offsets_for_many_overlapping_windows(tmp_path, monkeypatch):
    block = (slice(0, 1),) * 3
    windows = [SimpleNamespace(bounds=block, shape=(1, 1, 1))] * 40
    for index in range(len(windows)):
        counts, _ = window_jobs.window_files(tmp_path, index)
        np.save(counts, np.ones((1, 3), np.uint32))
    ranges = np.zeros((1, 1), np.int32), np.ones((1, 1), np.int32)
    from_bounds = time_windows.Columns.from_bounds
    calls = 0

    def count_calls(*args):
        nonlocal calls
        calls += 1
        return from_bounds(*args)

    monkeypatch.setattr(time_windows.Columns, "from_bounds", count_calls)
    monkeypatch.setattr(time_windows, "available_memory", lambda: 2**30)
    read = time_windows.count_blocks(tmp_path, windows, ranges)
    np.testing.assert_array_equal(read(block), np.full((3, 1, 1, 1), 40, np.uint32))
    np.testing.assert_array_equal(read(block), np.full((3, 1, 1, 1), 40, np.uint32))
    assert calls == len(windows)


def _kill_second_window(index, window):
    if index == 1:
        os.kill(os.getpid(), signal.SIGKILL)
    return time_windows.trace_window(index, window)


@needs_untwine
def test_interrupted_occlusion_resumes_with_identical_products(tmp_path, monkeypatch):
    cloud = tmp_path / "flight.laz"
    time, xyz = flight(cloud, seconds=2, rate=300)
    copc, cog = tmp_path / "resumed.copc.laz", tmp_path / "resumed.tif"

    def interrupt(work, specs, slots, worker, init, initargs, stage):
        return window_jobs.execute(work, specs, slots, _kill_second_window, init, initargs, stage)

    with monkeypatch.context() as patch:
        patch.setattr(occlusion, "execute", interrupt)
        with pytest.raises(OutOfMemory, match="Occlusion window 2 was killed"):
            occlusion.run([cloud], time, xyz, copc, cog, voxel=1, threads=2)

    resumed = occlusion.run([cloud], time, xyz, copc, cog, voxel=1, threads=2)
    clean_copc, clean_cog = tmp_path / "clean.copc.laz", tmp_path / "clean.tif"
    occlusion.run([cloud], time, xyz, clean_copc, clean_cog, voxel=1, threads=2)
    assert resumed["time_windows"]["reused_windows"] == 1
    assert copc.read_bytes() == clean_copc.read_bytes()
    assert cog.read_bytes() == clean_cog.read_bytes()


def test_concurrent_cli_run_preserves_shared_outputs(tmp_path):
    cloud = tmp_path / "cloud.las"
    flight(cloud, seconds=2, rate=300)
    copc, cog = tmp_path / "voxels.copc.laz", tmp_path / "occlusion.tif"
    path = tmp_path / "path.csv"
    for product in (copc, cog, path):
        product.write_bytes(b"owner output")
    with window_jobs.lock_output(copc):
        assert main([str(cloud), "-o", str(path), "--occlusion-copc", str(copc), "--occlusion-cog", str(cog)]) == 2
    assert all(product.read_bytes() == b"owner output" for product in (copc, cog, path))


@needs_untwine
def test_parallel_and_sequential_windows_write_identical_products(tmp_path, monkeypatch):
    # Admit both tiny fixture windows despite their fixed production reservations.
    monkeypatch.setattr(window_jobs, "available_memory", lambda: 6 * 2**30)
    cloud = tmp_path / "flight.las"
    time, xyz = flight(cloud, seconds=2, rate=300)
    results = []
    for count in (1, 2):
        copc, cog = (tmp_path / f"{count}.{ext}" for ext in ("copc.laz", "tif"))
        report = occlusion.run([cloud], time, xyz, copc, cog, voxel=1, threads=count)
        assert report["time_windows"]["windows"] == count and report["time_windows"]["peak_concurrent_workers"] == count
        results.append((copc.read_bytes(), cog.read_bytes(), report["voxels_ground_to_canopy"]))
    assert results[0] == results[1]


@needs_untwine
def test_separated_strip_windows_match_original_order(tmp_path):
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
def test_copc_occlusion_refuses_and_preserves_estimated_path(tmp_path):
    original, source = tmp_path / "flight.laz", tmp_path / "flight.copc.laz"
    flight(original, seconds=4, rate=500)
    with laspy.open(original) as src:
        occlusion._index_copc(original, source, src.header.point_count, tmp_path)
    assert main([str(source), "--debug", "-o", str(tmp_path / "path.csv"),
                 "--occlusion-copc", str(tmp_path / "out.copc.laz"),
                 "--occlusion-cog", str(tmp_path / "out.tif")]) == 0
    assert (tmp_path / "path.csv").stat().st_size
    assert "Supply the standardized files instead of COPC" in (tmp_path / "path.report.json").read_text()
    assert "# warning: Occlusion: Occlusion requires scanner-order LAS/LAZ." in (tmp_path / "path.csv").read_text()
    assert not (tmp_path / "out.tif").stat().st_size


@pytest.mark.parametrize("voxel", [0.3, 0.5, 1.0])
def test_time_window_counts_match_whole_grid(tmp_path, monkeypatch, voxel):
    from threedtrees_trajectory import window_input

    monkeypatch.setattr(window_input, "INDEX_POINTS", 51)
    cloud = tmp_path / "flight.laz"
    time, xyz = flight(cloud, seconds=2, rate=300)
    grid, _ = occlusion._grid([cloud], voxel)
    counters = dict(echoes_read=0, echoes_outside_path=0, traversed_pulses=0, missing_returns=0)
    expected = trace_reference([cloud], grid, time, xyz, counters)
    index = ScanIndex([cloud], None)
    index.populate()
    bounds = time[0] - 1, time[-1] + 1
    layers = np.zeros(grid.shape[:2], np.int32), np.full(grid.shape[:2], grid.shape[2], np.int32)
    windows = time_windows.plan(index, grid, time, xyz, bounds, 6, occlusion.DEFAULT_BUDGET, layers)
    assert len(windows) >= 2
    time_windows.init_worker(index, grid, time, xyz, tmp_path, (), layers, bounds)
    stats = [time_windows.trace_window(number, window) for number, window in enumerate(windows)]
    block = tuple(slice(0, size) for size in grid.shape)
    (_, _, summed), = time_windows.summed_counts(tmp_path, windows, [block], layers)
    np.testing.assert_array_equal(summed, expected)
    assert sum(window["traversed_pulses"] for window in stats) == len(time)


@pytest.fixture
def four_read_workers(monkeypatch):
    # Give tiny fixtures the production planner's reservation for four decoders.
    monkeypatch.setattr(parallel_read, "available_memory", lambda: 4 * parallel_read.WORKER_BYTES)


def test_parallel_record_reads_preserve_order_and_extent(tmp_path, monkeypatch, four_read_workers):
    cloud = tmp_path / "flight.laz"
    time, _ = flight(cloud, seconds=3, rate=400)
    monkeypatch.setattr(pulses, "CHUNK", 313)
    monkeypatch.setattr(pulses, "PROBE", 100)
    expected_index, _ = occlusion.prepare_index([cloud], 0.5)
    actual_index, _ = occlusion.prepare_index([cloud], 0.5)
    expected, serial = read_pulses([cloud], max_multi=700, scan_index=expected_index, threads=1)
    actual, parallel = read_pulses([cloud], max_multi=700, scan_index=actual_index, threads=4)
    assert parallel.read_workers == 4 and parallel.points_read == serial.points_read
    assert parallel.return_modes == serial.return_modes
    for field in ("time", "first", "last", "returns", "source", "file", "scan_angle"):
        np.testing.assert_array_equal(getattr(actual, field), getattr(expected, field))
    clock = time[100], time[-100]
    for first, second in zip(expected_index.surfaces(*clock, {}), actual_index.surfaces(*clock, {})):
        np.testing.assert_array_equal(first, second)


@pytest.fixture
def pf7_cloud(tmp_path):
    original, compressed = (tmp_path / name for name in ("pf7.las", "pf7.laz"))
    flight(original, seconds=2, rate=300, point_format=7)
    data = laspy.read(original)
    _, pulse = np.unique(data.gps_time, return_inverse=True)
    data.scanner_channel = (pulse % 4).astype(np.uint8)
    data.point_source_id = (1 + pulse % 3).astype(np.uint16)
    data.red = (1 + pulse).astype(np.uint16)
    data.intensity = (2 + pulse).astype(np.uint16)
    data.add_extra_dim(laspy.ExtraBytesParams(name="unused", type=np.uint16))
    data.unused = (3 + pulse).astype(np.uint16)
    data.write(original)
    data.write(compressed)
    return original, compressed, data


def test_pf7_selective_read_preserves_required_values(pf7_cloud, monkeypatch, four_read_workers):
    original, compressed, data = pf7_cloud
    monkeypatch.setattr(pulses, "CHUNK", 313)
    monkeypatch.setattr(pulses, "PROBE", 100)
    expected, _ = read_pulses([original], threads=1)
    for threads in (1, 4):
        actual, info = read_pulses([compressed], threads=threads)
        if threads == 4:
            assert info.read_workers == 4
        for field in ("time", "first", "last", "returns", "source", "file", "scan_angle"):
            np.testing.assert_array_equal(getattr(actual, field), getattr(expected, field))

    with cloud_io.open_cloud(compressed) as reader:
        chunk = reader.read()
    for field in ("red", "intensity", "unused"):
        assert not np.array_equal(chunk[field], data[field])

    index = ScanIndex([compressed], None)
    index.populate()
    chunk, = index.read(0, index.blocks[0], group_pulses=False)
    for field in ("X", "Y", "Z", "gps_time", "scanner_channel", "return_number", "number_of_returns"):
        np.testing.assert_array_equal(chunk[field], data[field])
    assert not np.array_equal(chunk.scan_angle, data.scan_angle)
    traced = list(index.read(0, index.blocks[0]))
    np.testing.assert_array_equal(np.concatenate([batch.point_source_id for batch in traced]), data.point_source_id)


def test_first_read_surfaces_decode_only_coverage_boundary(tmp_path, monkeypatch):
    from threedtrees_trajectory import window_input

    monkeypatch.setattr(window_input, "INDEX_POINTS", 51)
    cloud = tmp_path / "flight.laz"
    time, _ = flight(cloud, seconds=2, rate=300)
    index, _ = occlusion.prepare_index([cloud], 0.5)
    read_pulses([cloud], scan_index=index)
    low, high = time[100], time[-100]
    stage = {}
    minimum, maximum = index.surfaces(low, high, stage)
    assert 0 < stage["boundary_points_decoded"] <= 2 * window_input.INDEX_POINTS
    assert stage["echoes_read"] == index.grid.points
    las = laspy.read(cloud)
    points = las.xyz[(las.gps_time >= low) & (las.gps_time < high)]
    cells = np.floor((points - index.grid.origin) / index.grid.voxel + 1e-8).astype(int)
    np.testing.assert_array_equal(np.unique(cells[:, :2], axis=0), np.argwhere(maximum >= 0))
    assert (minimum[maximum >= 0] <= maximum[maximum >= 0]).all()


def test_copc_node_window_keeps_file_order_and_respects_all_budgets(monkeypatch, tmp_path):
    from types import SimpleNamespace

    monkeypatch.setattr(cloud_io, "COPC_BATCH_BYTES", 20)
    monkeypatch.setattr(cloud_io, "COPC_BATCH_NODES", 2)
    nodes = [
        SimpleNamespace(offset=offset, point_count=points, byte_size=length)
        for offset, points, length in [(40, 7, 10), (10, 2, 8), (70, 1, 6), (0, 0, 0), (30, 7, 11), (20, 2, 8), (60, 1, 15)]
    ]
    batches = cloud_io._node_batches(nodes, 10, tmp_path / "flight.copc.laz")
    assert [[node.offset for node in batch] for batch in batches] == [[10, 20], [30], [40], [60], [70]]
