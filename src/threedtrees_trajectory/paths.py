"""Several delivered paths for one cloud: which path each part of the cloud was scanned from.

The matching unit is a pulse segment: pulses of one input file and point source in one contiguous GPS-time stretch
(split at gaps over ``SEGMENT_GAP_S``). That covers one merged file with several flights, several files with one path,
and several files with several paths alike. Each segment uses the path that covers at least half of its time span;
where several do (a duplicate export, a second platform at the same time, the same week seconds on another day), the
path whose positions best explain the segment's pulse lines wins.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field

import numpy as np

SEGMENT_GAP_S = 5.0  # a pause in GPS time that starts a new segment, and a gap inside a path that it does not cover
SEGMENT_EDGE_S = 1.0  # widening of a segment's span when tracing, for pulses thinned away when it was found
MIN_SEGMENT_COVERAGE = 0.5  # share of a segment's time span a path must cover to be a candidate
MIN_PULSE_COVERAGE = 0.5  # share of pulses that must get a path, else occlusion is refused
MISS_SAMPLE = 20_000  # multi-return pulses per segment used to choose between overlapping paths
DUPLICATE_M = 0.01  # paths with the same times and positions within this distance are one path


@dataclass
class DeliveredPath:
    name: str
    time: np.ndarray
    xyz: np.ndarray

    def pieces(self):
        """(start, end) of the stretches without a gap over ``SEGMENT_GAP_S``, each widened by ``SEGMENT_EDGE_S``
        (the margin a pulse may lie beyond a path's ends, as for one path)."""
        cut = np.flatnonzero(np.diff(self.time) > SEGMENT_GAP_S) + 1
        starts = np.r_[0, cut]
        ends = np.r_[cut - 1, len(self.time) - 1]
        return self.time[starts] - SEGMENT_EDGE_S, self.time[ends] + SEGMENT_EDGE_S

    def covered(self, lo: float, hi: float) -> float:
        """Share of [lo, hi] inside this path's stretches."""
        if hi <= lo:
            starts, ends = self.pieces()
            return float(np.any((starts <= lo) & (ends >= lo)))
        starts, ends = self.pieces()
        overlap = np.clip(np.minimum(ends, hi) - np.maximum(starts, lo), 0, None).sum()
        return float(overlap / (hi - lo))


@dataclass
class Segment:
    file: int
    source: int
    start: float
    end: float
    pulses: int
    lo: float = 0.0  # span used when tracing: widened by SEGMENT_EDGE_S where the pulses stop, not at path cuts
    hi: float = 0.0
    path: int = -1  # index into PathSet.paths, -1 when no path covers it
    candidates: list[int] = field(default_factory=list)
    misses: dict[int, float | None] = field(default_factory=dict)


def same_path(a: DeliveredPath, b: DeliveredPath) -> bool:
    return (len(a.time) == len(b.time) and np.array_equal(a.time, b.time)
            and float(np.max(np.linalg.norm(a.xyz - b.xyz, axis=1))) < DUPLICATE_M)


def segments(pulses, paths=()) -> tuple[list[Segment], np.ndarray, list[int]]:
    """Pulse segments: one input file and point source, split where GPS time pauses for more than 5 s and where any
    path starts or ends (so each part lies wholly inside or outside every path). Returns the segments, the pulse order
    they refer to (file, source, time) and each segment's first position in that order."""
    if not len(pulses):
        return [], np.empty(0, np.int64), []
    order = np.lexsort((pulses.time, pulses.source, pulses.file))
    t, source, file = pulses.time[order], pulses.source[order], pulses.file[order]
    cut = np.flatnonzero((file[1:] != file[:-1]) | (source[1:] != source[:-1]) | (np.diff(t) > SEGMENT_GAP_S)) + 1
    edges = np.unique(np.concatenate([np.concatenate(p.pieces()) for p in paths])) if len(paths) else np.empty(0)
    found, first = [], []
    for a, b in zip(np.r_[0, cut], np.r_[cut, len(t)]):
        inner = edges[(edges > t[a]) & (edges < t[b - 1])]
        bounds = np.r_[a, a + np.searchsorted(t[a:b], inner, side="left"), b]
        bounds = np.unique(bounds)
        for u, v in zip(bounds[:-1], bounds[1:]):
            start, end = float(t[u]), float(t[v - 1])
            lo = start - SEGMENT_EDGE_S if u == a else float(inner[inner <= t[u]].max())
            hi = end + SEGMENT_EDGE_S if v == b else float(inner[inner > t[v - 1]].min())
            if u != a:
                found[-1].hi = lo  # contiguous with the previous part: no unassigned sliver between them
            found.append(Segment(int(file[u]), int(source[u]), start, end, int(v - u), lo, hi))
            first.append(int(u))
    return found, order, first


def _miss(first, last, sensor) -> np.ndarray:
    delta = last - first
    length = np.linalg.norm(delta, axis=1)
    usable = length >= 0.5
    return np.linalg.norm(np.cross(sensor[usable] - first[usable], delta[usable]), axis=1) / length[usable]


def _sample(indices: np.ndarray, size: int) -> np.ndarray:
    """A fixed, evenly spaced subset (deterministic)."""
    if len(indices) <= size:
        return indices
    return indices[np.linspace(0, len(indices) - 1, size).astype(np.int64)]


class PathSet:
    """Delivered paths plus the path each pulse segment uses."""

    def __init__(self, paths: list[DeliveredPath], chosen: list[Segment]):
        self.paths = paths
        self.segments = chosen
        self.table: dict[int, dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]]] = {}
        by_key: dict[tuple[int, int], list[Segment]] = {}
        for segment in chosen:
            if segment.path >= 0:
                by_key.setdefault((segment.file, segment.source), []).append(segment)
        for (file, source), items in by_key.items():
            items.sort(key=lambda s: s.start)
            self.table.setdefault(file, {})[source] = (
                np.array([s.lo for s in items]),
                np.array([s.hi for s in items]),
                np.array([s.path for s in items], np.int64))
        used = sorted({s.path for s in chosen if s.path >= 0})
        self.used = used

    @property
    def time_range(self) -> tuple[float, float]:
        spans = [(s.lo, s.hi) for s in self.segments if s.path >= 0]
        return min(a for a, _ in spans), max(b for _, b in spans)

    def digest(self) -> str:
        key = hashlib.sha256()
        for path in self.paths:
            key.update(np.ascontiguousarray(path.time, "<f8").tobytes())
            key.update(np.ascontiguousarray(path.xyz, "<f8").tobytes())
        for s in self.segments:
            key.update(np.array([s.file, s.source, s.lo, s.hi, s.path], "<f8").tobytes())
        return key.hexdigest()[:24]

    def assign(self, seconds: np.ndarray, file: int, source: np.ndarray) -> np.ndarray:
        """Path index per pulse, -1 where no chosen segment holds it."""
        out = np.full(len(seconds), -1, np.int64)
        table = self.table.get(file)
        if not table or not len(seconds):
            return out
        source = np.asarray(source)
        for value in np.unique(source):
            entry = table.get(int(value))
            if entry is None:
                continue
            lo, hi, path = entry
            mask = np.flatnonzero(source == value)
            t = seconds[mask]
            k = np.searchsorted(lo, t, side="right") - 1
            ok = (k >= 0) & (t <= hi[np.clip(k, 0, None)])
            out[mask[ok]] = path[k[ok]]
        return out

    def positions(self, seconds: np.ndarray, which: np.ndarray):
        """Sensor x, y, z per pulse from its assigned path (``which`` >= 0)."""
        x, y, z = (np.empty(len(seconds)) for _ in range(3))
        for k in np.unique(which):
            mask = which == k
            path = self.paths[int(k)]
            for axis, values in zip((x, y, z), path.xyz.T):
                axis[mask] = np.interp(seconds[mask], path.time, values)
        return x, y, z

    def bounds(self, first: float, last: float):
        """Box around every used path's positions between ``first`` and ``last`` (conservative)."""
        boxes = []
        for k in self.used:
            path = self.paths[k]
            if path.time[-1] < first or path.time[0] > last:
                continue
            left, right = np.searchsorted(path.time, [first, last])
            ends = np.array([np.interp([first, last], path.time, axis) for axis in path.xyz.T]).T
            points = np.concatenate((ends, path.xyz[left:right]))
            boxes.append((points.min(axis=0), points.max(axis=0)))
        if not boxes:
            points = np.concatenate([self.paths[k].xyz for k in self.used])
            return points.min(axis=0), points.max(axis=0)
        return np.min([b[0] for b in boxes], axis=0), np.max([b[1] for b in boxes], axis=0)


def choose(pulses, paths: list[DeliveredPath]) -> tuple[list[Segment], float, list[str]]:
    """Each segment's path, the share of pulses with a path, and warnings."""
    found, order, starts = segments(pulses, paths)
    warnings = []
    ambiguous = 0
    for segment, start in zip(found, starts):
        segment.candidates = [k for k, p in enumerate(paths) if p.covered(segment.start, segment.end) >= MIN_SEGMENT_COVERAGE]
        if len(segment.candidates) == 1:
            segment.path = segment.candidates[0]
        elif len(segment.candidates) > 1:
            members = order[start:start + segment.pulses]
            lines = _sample(members[pulses.returns[members] > 1], MISS_SAMPLE)
            first, last = pulses.first[lines].astype(float), pulses.last[lines].astype(float)
            best = None
            for k in segment.candidates:
                path = paths[k]
                sensor = np.column_stack([np.interp(pulses.time[lines], path.time, axis) for axis in path.xyz.T]) - pulses.origin
                miss = _miss(first, last, sensor)
                segment.misses[k] = float(np.median(miss)) if len(miss) else None
                if segment.misses[k] is not None and (best is None or segment.misses[k] < segment.misses[best]):
                    best = k
            if best is None:
                best = segment.candidates[0]
                ambiguous += 1
            segment.path = best
    total = sum(s.pulses for s in found)
    covered = sum(s.pulses for s in found if s.path >= 0)
    share = covered / total if total else 0.0
    if total and covered < total:
        warnings.append(f"{1 - share:.0%} of pulses had no delivered path covering their time; they are not traced.")
    if ambiguous:
        warnings.append(f"{ambiguous} segment(s) had several delivered paths and no return lines to choose between "
                        "them; the first path given was used.")
    return found, share, warnings


def assigned_miss(pulses, pathset: PathSet) -> dict:
    """Median and p95 miss of the return lines to their own segment's path (all assigned segments)."""
    which = np.full(len(pulses), -1, np.int64)
    for file in np.unique(pulses.file):
        mask = np.flatnonzero(pulses.file == file)
        which[mask] = pathset.assign(pulses.time[mask], int(file), pulses.source[mask])
    lines = np.flatnonzero((pulses.returns > 1) & (which >= 0))
    lines = lines[::max(1, -(-len(lines) // 1_000_000))]
    if not len(lines):
        return dict(usable_pulses=0, threshold_m=0.1, median_m=None, p95_m=None, threshold_exceeded=None)
    sensor = np.column_stack(pathset.positions(pulses.time[lines], which[lines])) - pulses.origin
    miss = _miss(pulses.first[lines].astype(float), pulses.last[lines].astype(float), sensor)
    if not len(miss):
        return dict(usable_pulses=0, threshold_m=0.1, median_m=None, p95_m=None, threshold_exceeded=None)
    median, p95 = np.percentile(miss, [50, 95])
    return dict(usable_pulses=len(miss), threshold_m=0.1, median_m=float(median), p95_m=float(p95),
                threshold_exceeded=bool(median > 0.1))
