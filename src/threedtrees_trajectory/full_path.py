"""Fit all paired LAS/LAZ pulses without collecting whole-cloud endpoints."""
from __future__ import annotations

import multiprocessing
from collections import defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from functools import partial
from itertools import groupby, islice
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from . import solve
from .cloud import is_copc, open_cloud, read_points_exact
from .estimate import (
    MIN_LINE_POSITIONS, MULTI_SHARE_MIN, UNRELIABLE_LINE_MISS,
    detect_platform, line_result, world_result,
)
from .parallel_read import SerialReadRequired, WORKER_BYTES, index_range, read_range
from .pulses import (
    CHUNK, CloudInfo, Pulses, _Store,
    _angles, _coordinates, _crs, _join_stray, _pulse_groups, _time_system, cloud_crs, crs_labels,
)
from .resources import available_memory

PLATFORM_POINTS = 1_000_000
MIN_LINES = 50
PULSE_COLUMNS = ("time", "first", "last", "returns", "source", "file", "angle")
FRAGMENT_COLUMNS = ("time", "xyz", "number", "source", "file", "angle", "total", "channel")
LINE_RECORD = np.dtype([("time", "f8"), ("first", "f4", 3), ("last", "f4", 3), ("source", "u2")])
STRAY_RECORD = np.dtype([("time", "f8"), ("xyz", "f4", 3), ("number", "u1"), ("total", "u1"),
                        ("source", "u2"), ("channel", "u1")])


class WindowReplayRequired(Exception):
    """The windows overlap too heavily for efficient range seeks."""


@dataclass
class PulseSummary:
    kept: int
    multi_share: float
    has_scan_angle: bool
    has_time: bool = True

    def __len__(self):
        return self.kept


@dataclass
class Window:
    start: int
    stop: int
    fits: list[solve.Fit]
    recovered: int
    echoes: int


@dataclass(frozen=True)
class Input:
    path: Path
    index: int
    points: int
    shift: np.ndarray
    record_size: int


@dataclass
class DiskWindow:
    path: Path
    offset: int
    lines: int
    fragments: int
    file: int


@dataclass
class Batch:
    low: float
    high: float
    points: int
    singles: int
    complete: int
    incomplete: int
    fragments: int
    multi_echoes: int
    angles: bool
    sample: _Store
    windows: dict[int, Window | DiskWindow]
    blocks: list
    surfaces: list


@dataclass
class FileScan:
    input: Input
    windows: dict[int, list[Window | DiskWindow]]
    sample: _Store
    counts: tuple[int, int, int, int, int]
    angles: bool
    blocks: list
    surfaces: list


def _endpoints(store, chunk, start, end, shift, *, angles=False):
    first = _coordinates(chunk, start, shift)
    last = first
    multi = start != end
    if multi.all():
        last = _coordinates(chunk, end, shift)
    elif multi.any():
        last = first.copy()
        last[multi] = _coordinates(chunk, end[multi], shift)
    store.add(time=chunk.gps_time[start], first=first, last=last, returns=(end - start + 1).astype(np.uint8),
              source=chunk.point_source_id[start], angle=_angles(chunk, start) if angles else np.empty(0, np.float32))


def _platform_sample(store, origin, summary):
    time = store.array("time", np.float64)
    return Pulses(origin, time, store.array("first", np.float32, 3),
                  store.array("last", np.float32, 3), store.array("returns", np.uint8),
                  store.array("source", np.uint16), store.array("file", np.uint16),
                  store.array("angle", np.float32), True, summary.has_scan_angle, summary.multi_share)


def _fit(lines, stray):
    joined = echoes = 0
    if stray.columns["time"]:
        joined, echoes = _join_stray(stray, lines)
    return _line_fits(lines), joined, echoes


def _line_fits(lines):
    time = lines.array("time", np.float64)
    first = lines.array("first", np.float32, 3).astype(np.float64)
    line = first - lines.array("last", np.float32, 3)
    length = np.sqrt(np.einsum("ij,ij->i", line, line))
    index = np.flatnonzero(length >= 1.0)
    fits = []
    if len(index) >= MIN_LINES:
        sources = lines.array("source", np.uint16)[index]
        groups = [index]
        if sources.min() != sources.max():
            order = np.argsort(sources, kind="stable")
            sources = sources[order]
            groups = np.split(index[order], np.flatnonzero(sources[1:] != sources[:-1]) + 1)
        for group in groups:
            if len(group) >= MIN_LINES:
                fit = solve.fit_lines(first[group], line[group] / length[group, None], time[group],
                                      True, MIN_LINES, compact=True)
                if fit is not None:
                    fits.append(fit)
    return fits


def _append_groups(lines, stray, chunk, shift, groups, indices):
    start, end, _, complete, rn, nr, source, channel = groups
    paired = indices[complete[indices]]
    _endpoints(lines, chunk, start[paired], end[paired], shift)
    fragments = indices[~complete[indices]]
    lo, hi = start[fragments], end[fragments]
    if len(lo):
        if np.any(lo != hi):
            lengths = hi - lo + 1
            offsets = np.cumsum(lengths) - lengths
            lo = np.repeat(lo - offsets, lengths) + np.arange(lengths.sum())
        stray.add(time=chunk.gps_time[lo], xyz=_coordinates(chunk, lo, shift), number=rn[lo], total=nr[lo],
                  source=source[lo], channel=channel[lo], file=np.zeros(len(lo), np.uint16), angle=_angles(chunk, lo))


def _window_ids(times, interval):
    return np.floor(times / interval).astype(np.int64)


def _spool(stream, lines, stray, file):
    offset = stream.tell()
    counts = []
    for store, dtype in ((lines, LINE_RECORD), (stray, STRAY_RECORD)):
        records = np.empty(sum(map(len, store.columns["time"])), dtype)
        for name in dtype.names:
            records[name] = store.array(name, records[name].dtype, width=3 if records[name].ndim == 2 else 0)
        records.tofile(stream)
        counts.append(len(records))
    return DiskWindow(Path(stream.name), offset, *counts, file)


def _read_window(item):
    key, parts = item
    lines, recovered = _Store(*PULSE_COLUMNS), _Store(*PULSE_COLUMNS)
    counts = {}
    for file, local in groupby(parts, key=lambda part: part.file):
        stray = _Store(*FRAGMENT_COLUMNS)
        for part in local:
            with part.path.open("rb") as stream:
                stream.seek(part.offset)
                paired = np.fromfile(stream, dtype=LINE_RECORD, count=part.lines)
                fragments = np.fromfile(stream, dtype=STRAY_RECORD, count=part.fragments)
            if len(paired) != part.lines or len(fragments) != part.fragments:
                raise ValueError("A temporary path window was truncated.")
            lines.add(**{name: paired[name] for name in LINE_RECORD.names})
            stray.add(**{name: fragments[name] for name in STRAY_RECORD.names},
                      file=np.full(len(fragments), file, np.uint16), angle=np.zeros(len(fragments), np.float32))
        if stray.columns["time"]:
            counts[file] = _join_stray(stray, recovered)
    # The whole-cloud reader appends recovered lines after regular endpoints.
    lines.extend(recovered)
    return key, _line_fits(lines), counts


def _scan(source, interval, stride, spool, grid, record):
    path, total, shift = source.path, source.points, source.shift
    left, chunk = read_range(path, record, min(record + CHUNK, total), total)
    if not len(chunk):
        raise SerialReadRequired
    time = chunk.gps_time
    low, high = float(time.min()), float(time.max())
    if not np.isfinite([low, high]).all() or max(abs(low), abs(high)) / np.iinfo(np.int64).max >= interval:
        raise SerialReadRequired
    groups = _pulse_groups(chunk, time, False)
    start, end, single, complete, *_ = groups
    incomplete = ~single & ~complete
    fragments = int((end[incomplete] - start[incomplete] + 1).sum())
    sample = _Store(*PULSE_COLUMNS)
    kept = np.flatnonzero(~incomplete)
    subset = kept[::stride]
    _endpoints(sample, chunk, start[subset], end[subset], shift, angles=True)
    sample.add(file=np.full(len(subset), source.index, np.uint16))
    raw_angles = chunk.scan_angle if chunk.point_format.id >= 6 else chunk.scan_angle_rank
    limit = 15_000 if chunk.point_format.id >= 6 else 90
    angles = raw_angles[start[kept]]
    has_angles = bool(np.any((angles != 0) & (angles > -limit) & (angles < limit)))

    indices = np.flatnonzero(~single)
    keys = _window_ids(time[start[indices]], interval)
    if np.any(keys[1:] < keys[:-1]):
        order = np.argsort(keys, kind="stable")
        keys, indices = keys[order], indices[order]
    if spool is None:
        windows = _collect_windows(chunk, groups, keys, indices, left, shift)
    else:
        with (spool / f"{source.index}-{record}.bin").open("xb") as stream:
            windows = _collect_windows(chunk, groups, keys, indices, left, shift, stream, source.index)
    blocks, surfaces = index_range(path, left, chunk, grid) if grid is not None else ([], [])
    return Batch(low, high, len(chunk), int(single.sum()), int(complete.sum()), int(incomplete.sum()),
                 fragments, int((end[complete] - start[complete] + 1).sum()), has_angles, sample, windows, blocks, surfaces)


def _collect_windows(chunk, groups, keys, indices, left, shift, stream=None, file=0):
    if not len(indices):
        return {}
    start, end = groups[:2]
    windows = {}
    edges = np.r_[0, 1 + np.flatnonzero(keys[1:] != keys[:-1]), len(keys)]
    for lo, hi in zip(edges[:-1], edges[1:]):
        pick = indices[lo:hi]
        lines, stray = _Store(*PULSE_COLUMNS), _Store(*FRAGMENT_COLUMNS)
        _append_groups(lines, stray, chunk, shift, groups, pick)
        if stream is None:
            window = Window(left + int(start[pick[0]]), left + int(end[pick[-1]]) + 1, *_fit(lines, stray))
        else:
            window = _spool(stream, lines, stray, file)
        windows[int(keys[lo])] = window
    return windows


def _refit(path, shift, interval, item):
    key, parts = item
    lines, stray = _Store(*PULSE_COLUMNS), _Store(*FRAGMENT_COLUMNS)
    with open_cloud(path, sequential=True) as reader:
        for part in parts:
            reader.seek(part.start)
            chunk = read_points_exact(reader, part.stop - part.start)
            groups = _pulse_groups(chunk, chunk.gps_time, False)
            start, _, single, *_ = groups
            pick = np.flatnonzero(~single & (_window_ids(chunk.gps_time[start], interval) == key))
            _append_groups(lines, stray, chunk, shift, groups, pick)
    return key, Window(parts[0].start, parts[-1].stop, *_fit(lines, stray))


def _bounded_map(pool, workers, fn, items):
    if pool is None:
        yield from map(fn, items)
        return
    items = iter(items)
    pending = deque(pool.submit(fn, value) for value in islice(items, workers))
    while pending:
        yield pending.popleft().result()
        if (value := next(items, None)) is not None:
            pending.append(pool.submit(fn, value))


def _scan_file(pool, workers, source, interval, stride, grid, spool):
    scan = partial(_scan, source, interval, stride, spool, grid)
    sample = _Store(*PULSE_COLUMNS)
    windows = defaultdict(list)
    low, high = np.inf, -np.inf
    points = singles = complete = incomplete = fragments = multi_echoes = 0
    angles = False
    blocks, surfaces = [], []
    for batch in _bounded_map(pool, workers, scan, range(0, source.points, CHUNK)):
        low, high = min(low, batch.low), max(high, batch.high)
        points += batch.points
        singles += batch.singles
        complete += batch.complete
        incomplete += batch.incomplete
        fragments += batch.fragments
        multi_echoes += batch.multi_echoes
        angles |= batch.angles
        sample.extend(batch.sample)
        blocks.extend(batch.blocks)
        surfaces.extend(batch.surfaces)
        for key, window in batch.windows.items():
            windows[key].append(window)
    if points != source.points or low == high:
        raise SerialReadRequired
    return FileScan(source, windows, sample, (singles, complete, incomplete, fragments, multi_echoes),
                    angles, blocks, surfaces)


def _fit_in_place(pool, workers, scan, interval):
    tasks = [(key, parts) for key, parts in scan.windows.items() if len(parts) > 1]
    spans = [sum(p.stop - p.start for p in parts) for _, parts in tasks]
    if spans and (max(spans) > CHUNK or sum(spans) > scan.input.points // 4):
        raise WindowReplayRequired
    fit = partial(_refit, scan.input.path, scan.input.shift, interval)
    for key, window in _bounded_map(pool, workers, fit, tasks):
        scan.windows[key] = [window]
    recovered, echoes = 0, 0
    fits = []
    for key in sorted(scan.windows):
        window, = scan.windows[key]
        fits.extend(window.fits)
        recovered += window.recovered
        echoes += window.echoes
    return fits, {scan.input.index: (recovered, echoes)}


def _fit_spooled(pool, workers, scans):
    windows = defaultdict(list)
    for scan in scans:
        for key, parts in scan.windows.items():
            windows[key].extend(parts)
    tasks = [(key, windows[key]) for key in sorted(windows)]
    if any(sum(p.lines + p.fragments for p in parts) > CHUNK for _, parts in tasks):
        raise SerialReadRequired
    counts = defaultdict(lambda: (0, 0))
    fits = []
    for _, fitted, by_file in _bounded_map(pool, workers, _read_window, tasks):
        fits.extend(fitted)
        for file, (joined, echoes) in by_file.items():
            old_joined, old_echoes = counts[file]
            counts[file] = old_joined + joined, old_echoes + echoes
    return fits, counts


def _needs_spool(source):
    if source.points <= CHUNK:
        return False
    with open_cloud(source.path, sequential=True) as reader:
        chunk = read_points_exact(reader, min(source.points, 65_536))
    _, _, single, complete, *_ = _pulse_groups(chunk, chunk.gps_time, False)
    return np.count_nonzero(~single & ~complete) > len(single) // 20


def _fit_ranges(inputs, interval, check_platform, workers, grid, scratch):
    total = sum(source.points for source in inputs)
    stride = max(1, -(-total // PLATFORM_POINTS)) if check_platform else CHUNK

    def scan(pool):
        if len(inputs) == 1 and not _needs_spool(inputs[0]):
            try:
                first = _scan_file(pool, workers, inputs[0], interval, stride, grid, None)
                fits, counts = _fit_in_place(pool, workers, first, interval)
                return [first], fits, counts
            except WindowReplayRequired:
                # Re-read each input range once to group widely separated windows on disk.
                pass
        scratch.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(dir=scratch, prefix=".path-windows-") as folder:
            scans = [_scan_file(pool, workers, source, interval, stride, grid, Path(folder)) for source in inputs]
            fits, counts = _fit_spooled(pool, workers, scans)
            return scans, fits, counts

    if workers == 1:
        return scan(None)
    with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        return scan(pool)


def _inspect_inputs(paths):
    headers = []
    for index, path in enumerate(paths):
        with open_cloud(path) as reader:
            header = reader.header
            if is_copc(header):
                raise SerialReadRequired
            if header.point_count:
                headers.append((index, path, header))
    if not headers:
        raise SerialReadRequired
    lower = np.min([header.mins for _, _, header in headers], axis=0)
    upper = np.max([header.maxs for _, _, header in headers], axis=0)
    origin = np.floor(lower / 1000) * 1000
    info = CloudInfo.from_paths(paths)
    inputs = []
    crs = None
    for index, path, header in headers:
        time_system = _time_system(header, set(header.point_format.dimension_names))
        if time_system == "none" or info.time_system not in ("none", time_system):
            raise SerialReadRequired
        info.time_system = time_system
        file_crs = _crs(header)
        if file_crs is not None:
            if (any(axis.unit_conversion_factor != 1.0 for axis in file_crs.axis_info)
                    or crs is not None and not crs.equals(file_crs, ignore_axis_order=True)):
                raise SerialReadRequired
            crs = file_crs
        if not np.isfinite(header.offsets).all() or not np.isfinite(header.scales).all() or np.any(header.scales <= 0):
            raise SerialReadRequired
        if not np.isfinite(header.maxs - header.mins).all():
            raise SerialReadRequired
        inputs.append(Input(path, index, header.point_count, header.offsets - origin, header.point_format.size))
        info.points += header.point_count
    info.crs, info.vertical_crs = crs_labels(cloud_crs(crs, lower, upper))
    info.has_gps_time = True
    return info, inputs, origin


def _summarize(scans, recovered, info):
    kept = multis = 0
    angles = False
    for scan in scans:
        singles, complete, incomplete, fragments, multi_echoes = scan.counts
        file = scan.input.index
        joined, echoes = recovered.get(file, (0, 0))
        unpaired, multi = fragments - echoes, complete + joined
        # The legacy reader retries these files with timestamp repair or scan-angle fitting.
        if not multi or unpaired > max(200, multi_echoes + echoes):
            raise SerialReadRequired
        multis += multi
        kept += singles + multi
        angles |= scan.angles
        info.multi_pulses += complete
        info.incomplete_multi_pulses += incomplete
        info.joined_out_of_order_points += fragments
        info.fragment_echoes += fragments
        info.complete_multi_echoes += multi_echoes
        info.unpaired_multi_echoes += unpaired
        info.return_modes.append({"file": str(scan.input.path), "index": file, "mode": "paired",
                                  "has_scan_angle": scan.angles, "unpaired_multi_echoes": unpaired,
                                  "line_pulses": multi})
    if not kept or multis / kept < MULTI_SHARE_MIN:
        raise SerialReadRequired
    return PulseSummary(kept, multis / kept, angles)


def try_estimate(paths, interval, workers, check_platform, scratch, scan_index=None):
    """Fit ordinary paired LAS/LAZ by window, with the existing reader as a compatibility fallback."""
    if len(paths) > np.iinfo(np.uint16).max + 1:
        return None
    try:
        info, inputs, origin = _inspect_inputs(paths)
        reservation = max(WORKER_BYTES, CHUNK * (max(source.record_size for source in inputs) + 224))
        workers = max(1, min(workers, max(-(-source.points // CHUNK) for source in inputs),
                             available_memory() // reservation))
        grid = scan_index.grid if scan_index is not None else None
        scans, fits, recovered = _fit_ranges(inputs, interval, check_platform, workers, grid, scratch)
        summary = _summarize(scans, recovered, info)
    except SerialReadRequired:
        return None
    fits = solve.filter_line_fits(fits)
    if len(fits) < MIN_LINE_POSITIONS:
        return None
    result = line_result(fits, "given")
    if not result.fits or np.median([fit.miss for fit in result.fits]) > UNRELIABLE_LINE_MISS:
        return None
    if check_platform:
        sample = _Store(*PULSE_COLUMNS)
        for scan in scans:
            sample.extend(scan.sample)
        evidence_pulses = _platform_sample(sample, origin, summary)
        if len(evidence_pulses) < 500:
            return None
        platform, evidence = detect_platform(evidence_pulses)
        if platform == "static":
            return None
        result.evidence = evidence

    result = world_result(result, origin)
    info.points_read = info.points
    info.read_workers = workers
    if scan_index is not None:
        for scan in scans:
            scan_index.blocks[scan.input.index] = scan.blocks
            scan_index.summaries[scan.input.index] = scan.surfaces
    return result, summary, info
