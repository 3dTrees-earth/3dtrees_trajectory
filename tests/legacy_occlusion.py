"""The earlier acquisition-order traversal for small reference tests."""
import numpy as np

from threedtrees_trajectory import occlusion
from threedtrees_trajectory.cloud import CHUNK, point_chunks, pulse_chunks
from threedtrees_trajectory.occlusion import _positions, _time_mask

FRAGMENT_BATCH = 1_000_000


def _trace_multi(tracer, coords, sensor, time, rn, nr, source, channel, ordered=False) -> np.ndarray:
    """Trace complete shots. OccPy retains incomplete shots after clearPulseDataset()."""
    same = (time[1:] == time[:-1]) & (source[1:] == source[:-1]) & (channel[1:] == channel[:-1])
    if not ordered:
        # On original strips, a second first return can start another shot
        # at the same GPS time. Spatially ordered data cannot disambiguate it.
        index = np.arange(len(time))
        group_start = np.maximum.accumulate(np.where(np.r_[True, ~same], index, 0))
        previous_first = np.maximum.accumulate(np.where(rn == 1, index, -1))
        same &= ~((rn[1:] == 1) & (previous_first[:-1] >= group_start[1:]))
    starts = np.flatnonzero(np.r_[True, ~same])
    counts = np.diff(np.r_[starts, len(time)])
    target = nr[starts]
    numbers = np.bitwise_or.reduceat(np.left_shift(np.uint32(1), rn), starts)
    complete = (counts == target) & (numbers == np.left_shift(np.uint32(1), target + 1) - 2)
    keep = np.repeat(complete, counts)
    if keep.any():
        shot = np.repeat(np.arange(len(starts), dtype=np.float64), counts)[keep]
        tracer.addPointData(*(axis[keep] for axis in (*coords, *sensor)), shot, rn[keep], nr[keep])
        tracer.doRaytracing()
        tracer.clearPulseDataset()
    return ~keep


def _recover_fragments(tracer, fragments, time, xyz) -> int:
    x, y, z, t, rn, nr, source, channel = (np.concatenate(parts) for parts in zip(*fragments))
    fragments.clear()
    order = np.lexsort((rn, channel, source, t))
    ordered_time = t[order]
    start = unpaired = 0
    while start < len(order):
        nominal_stop = min(len(order), start + FRAGMENT_BATCH)
        stop = int(np.searchsorted(ordered_time, ordered_time[nominal_stop - 1], side="right"))
        if stop - start > 2 * FRAGMENT_BATCH:
            raise ValueError("Too many separated echoes share a GPS time to reconstruct the LAZ pulses safely.")
        idx, times = order[start:stop], ordered_time[start:stop]
        coords = [axis[idx] for axis in (x, y, z)]
        missing = _trace_multi(tracer, coords, _positions(times, time, xyz), times, rn[idx], nr[idx],
                                source[idx], channel[idx], ordered=True)
        unpaired += int(missing.sum())
        start = stop
    return unpaired


def _trace_file(path, tracer, grid, time, xyz, counters, traced_records):
    chunks = pulse_chunks(point_chunks(path, CHUNK), CHUNK)
    fragments = []
    num_fragments = unpaired = 0
    for chunk in chunks:
        t = chunk.gps_time
        if not np.isfinite(t).all():
            raise ValueError(f"{path.name}: GPS times must be finite.")
        keep = _time_mask(t, time)
        counters["echoes_read"] += len(t)
        counters["echoes_outside_path"] += int((~keep).sum())
        if not keep.any():
            continue
        t = np.ascontiguousarray(t[keep])
        traced_records += len(t)
        # OccPy normalises by the box size before flooring. Nudge exact faces
        # inward by a negligible fraction of a voxel so tiling and z cropping
        # use a consistent half-open lattice at representable coordinates.
        coords = tuple(np.asarray(axis)[keep] + grid.voxel * 1e-7 for axis in (chunk.x, chunk.y, chunk.z))
        sensor = _positions(t, time, xyz)
        rn, nr = np.asarray(chunk.return_number)[keep], np.asarray(chunk.number_of_returns)[keep]
        single = (rn <= 1) & (nr <= 1)
        if single.any():
            tracer.doRaytracing_singleReturnPulses(*(axis[single] for axis in (*coords, *sensor)), t[single])
        if single.all():
            continue
        multi = ~single
        source = chunk.point_source_id[keep][multi]
        channel = np.asarray(chunk.scanner_channel)[keep][multi] if chunk.point_format.id >= 6 else np.zeros(len(source), np.uint8)
        multi_coords = [axis[multi] for axis in coords]
        shot = (t[multi], rn[multi], nr[multi], source, channel)
        pending = _trace_multi(tracer, multi_coords, [axis[multi] for axis in sensor], *shot)
        count = int(pending.sum())
        if count:
            num_fragments += count
            fragments.append(tuple(axis[pending] for axis in (*multi_coords, *shot)))
    if fragments:
        unpaired = _recover_fragments(tracer, fragments, time, xyz)
        counters["recovered_echoes"] = counters.get("recovered_echoes", 0) + num_fragments - unpaired
    counters["unpaired_echoes"] = counters.get("unpaired_echoes", 0) + unpaired
    return traced_records


def trace_reference(paths, grid, time, xyz, counters):
    from raytr import PyRaytracer

    tracer = PyRaytracer()
    hi = grid.origin + np.array(grid.shape) * grid.voxel
    tracer.defineGrid(grid.origin, hi, *grid.shape, grid.voxel)
    traced_records = 0
    for path in paths:
        traced_records = _trace_file(path, tracer, grid, time, xyz, counters, traced_records)
    counters["traversed_pulses"] += int(tracer.get_num_traversed_pulses())
    counters["missing_returns"] += int(tracer.get_num_missing_returns())
    counts = tuple(np.asarray(values(), np.int32) for values in (tracer.getNhit, tracer.getNmiss, tracer.getNocc))
    if any(values.min(initial=0) < 0 for values in counts):
        raise ValueError("Too many pulses traverse the same voxel for the native 32-bit count.")
    return counts
