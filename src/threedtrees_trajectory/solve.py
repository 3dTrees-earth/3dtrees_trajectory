"""Sensor positions from pulses.

lines       A pulse's first and last return lie on a line through the sensor. The point (and, for a moving sensor, the
            velocity) closest to all lines of a time window is the sensor position. Lines that miss it are dropped and
            the fit repeated. Needs multi-return pulses; centimetre accuracy on the datasets tested.
scan_angle  Single returns with a scan angle: the sensor sits where every point is seen at its scan angle off nadir
            (a cone around each point). Level flight and a line scanner are assumed; about 1 m horizontally and
            0.1 m vertically on the RIEGL and DJI flights tested, much worse on YellowScan (8 m).

A fixed position from all lines of a scan (``station_from_lines``) and the ground-density concentration only serve to
recognise a static (TLS) cloud, which 1.0.0 refuses.
"""
from __future__ import annotations

import multiprocessing
from collections import deque
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from dataclasses import dataclass
from itertools import islice

import numpy as np
from scipy import ndimage
from scipy.optimize import least_squares

from .pulses import Pulses
from .resources import OutOfMemory


@dataclass
class Fit:
    time: float
    xyz: np.ndarray
    pulses: int
    miss: float  # median distance of the used pulses from the fit, metres (method-specific)


def _line_vectors(p: Pulses, index: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    first = p.first[index].astype(np.float64)
    line = first - p.last[index]
    length = np.linalg.norm(line, axis=1)
    return first, line, length


def _line_indices(p: Pulses, min_length: float, group: np.ndarray | None = None) -> np.ndarray:
    """Find usable lines in bounded batches while preserving the original pulse order."""
    parts = []
    count = len(p) if group is None else len(group)
    block = 500_000
    for lo in range(0, count, block):
        hi = min(lo + block, count)
        if group is None:
            index = lo + np.flatnonzero(p.returns[lo:hi] > 1)
        else:
            index = group[lo:hi]
            index = index[p.returns[index] > 1]
        _, _, length = _line_vectors(p, index)
        parts.append(index[length >= min_length])
    return np.concatenate(parts) if parts else np.empty(0, dtype=np.intp)


def _closest_point(origin: np.ndarray, direction: np.ndarray, dt: np.ndarray | None, cond_max: float):
    """Point (and velocity when dt is given) closest to the lines; None when the lines don't fix it."""
    projector = np.eye(3)[None] - direction[:, :, None] * direction[:, None, :]
    if dt is None:
        matrix = projector.sum(axis=0)
        rhs = np.einsum("nij,nj->i", projector, origin)
    else:
        jacobian = np.concatenate([projector, projector * dt[:, None, None]], axis=2)
        matrix = np.einsum("nij,nik->jk", jacobian, jacobian)
        rhs = np.einsum("nij,ni->j", jacobian, np.einsum("nij,nj->ni", projector, origin))
    if not np.all(np.isfinite(matrix)) or np.linalg.cond(matrix) > cond_max:
        return None
    solution = np.linalg.solve(matrix, rhs)
    return solution[:3], (solution[3:] if dt is not None else np.zeros(3))


def fit_lines(origin: np.ndarray, direction: np.ndarray, time: np.ndarray, moving: bool, min_pulses: int) -> Fit | None:
    """One position from a set of pulse lines; the lines that miss it by more than 3 median misses are dropped."""
    centre = origin.mean(axis=0)
    p = origin - centre
    t0 = float(time.mean())
    dt = time - t0 if moving else None
    keep = np.ones(len(p), bool)
    fit, miss = None, None
    for _ in range(5):
        fit = _closest_point(p[keep], direction[keep], dt[keep] if moving else None, 1e10 if moving else 1e6)
        if fit is None:
            return None
        offset = p - (fit[0] + (dt[:, None] * fit[1] if moving else 0))
        miss = np.linalg.norm(offset - (offset * direction).sum(axis=1)[:, None] * direction, axis=1)
        limit = 3 * max(float(np.median(miss[keep])), 0.01)
        new_keep = miss <= limit
        if new_keep.sum() < min_pulses:
            return None
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep
    return Fit(t0, fit[0] + centre, int(keep.sum()), float(np.median(miss[keep])))


def _fit_windows(fit, starts, length, minimum, workers):
    windows = [(lo, hi) for lo, hi in zip(starts, np.r_[starts[1:], length]) if hi - lo >= minimum]
    if workers == 1 or len(windows) <= 1:
        return [value for window in windows if (value := fit(window)) is not None]
    with ThreadPoolExecutor(max_workers=min(workers, len(windows))) as pool:
        return [value for value in pool.map(fit, windows) if value is not None]


def path_from_lines(p: Pulses, interval: float = 0.5, min_pulses: int = 50, min_length: float = 1.0, *,
                    workers: int = 1) -> list[Fit]:
    """Moving sensor: one position (with velocity) per time window and point source."""
    index = _line_indices(p, min_length)
    if len(index) < min_pulses:
        return []
    source = p.source[index]
    window = np.floor(p.time[index] / interval).astype(np.int64)
    order = np.lexsort((source, window))
    index, source, window = index[order], source[order], window[order]
    starts = np.flatnonzero(np.r_[True, (window[1:] != window[:-1]) | (source[1:] != source[:-1])])
    def fit_window(window):
        lo, hi = window
        selected = index[lo:hi]
        origin, line, length = _line_vectors(p, selected)
        return fit_lines(origin, line / length[:, None], p.time[selected], True, min_pulses)

    fits = _fit_windows(fit_window, starts, len(index), min_pulses, workers)
    # A window whose lines don't meet (median miss 0.7-30 m against 0.02 m) fitted noise or mixed pulses: drop it.
    if fits:
        limit = max(0.25, 10 * float(np.median([f.miss for f in fits])))
        fits = [f for f in fits if f.miss <= limit]
    return sorted(fits, key=lambda f: f.time)


def station_from_lines(p: Pulses, min_pulses: int = 2000, min_length: float = 0.5,
                       max_pulses: int = 200_000, *, group: np.ndarray | None = None) -> Fit | None:
    """One fixed position from an even sample of lines to recognise a static scanner."""
    index = _line_indices(p, min_length, group)
    if len(index) < min_pulses:
        return None
    if len(index) > max_pulses:
        index = index[np.random.default_rng(0).choice(len(index), max_pulses, replace=False)]
    origin, line, length = _line_vectors(p, index)
    return fit_lines(origin, line / length[:, None], p.time[index], False, min_pulses)


def path_from_scan_angle(p: Pulses, interval: float = 0.5, min_pulses: int = 200, *,
                         group: np.ndarray | None = None, workers: int = 1) -> list[Fit]:
    """Moving sensor, single returns: per window, the position from which each point is seen at its scan angle."""
    # Preserve the mean reduction order, but expand coordinates and angles only for the current window.
    points = p.first if group is None else p.first[group]
    angles = p.scan_angle if group is None else p.scan_angle[group]
    if not len(points):
        return []
    centre = points.astype(np.float64).mean(axis=0)
    usable = np.flatnonzero(np.isfinite(angles) & (np.abs(angles) < 90))
    if group is not None:
        usable = group[usable]
    if len(usable) < min_pulses:
        return []
    order = usable[np.argsort(p.time[usable], kind="stable")]
    time = p.time[order]
    window = np.floor(time / interval).astype(np.int64)
    starts = np.flatnonzero(np.r_[True, window[1:] != window[:-1]])
    def arguments(lo, hi):
        selected = order[lo:hi]
        return p.first[selected], p.scan_angle[selected], time[lo:hi], centre

    windows = [(lo, hi) for lo, hi in zip(starts, np.r_[starts[1:], len(time)]) if hi - lo >= min_pulses]
    tasks = (arguments(lo, hi) for lo, hi in windows)
    if workers == 1 or len(windows) <= 1:
        return [fit for args in tasks if (fit := _fit_scan_angle(*args)) is not None]
    workers = min(workers, len(windows))
    fits = []
    # Spawn clean workers so native libraries and the caller's threads aren't inherited.
    # Keep only a bounded number of windows in flight, including the serialized inputs.
    try:
        with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("spawn")) as pool:
            pending = deque(pool.submit(_fit_scan_angle, *args) for args in islice(tasks, workers))
            while pending:
                fit = pending.popleft().result()
                if fit is not None:
                    fits.append(fit)
                if (args := next(tasks, None)) is not None:
                    pending.append(pool.submit(_fit_scan_angle, *args))
    except BrokenProcessPool as error:
        raise OutOfMemory("A scan-angle path worker stopped before completing its fit.") from error
    return fits


def _fit_scan_angle(points, angles, time, centre):
    x, y, z = (points.astype(np.float64) - centre).T
    tn = np.tan(np.radians(np.abs(angles.astype(float))))
    t0 = time.mean()
    dt = time - t0

    def residual(q):
        return np.hypot(x - q[0] - q[3] * dt, y - q[1] - q[4] * dt) - (q[2] - z) * tn

    high = np.percentile(z, 95)
    start = [np.median(x), np.median(y), high + 50, 0.0, 0.0]
    result = least_squares(residual, start, loss="soft_l1", f_scale=0.5)
    if not result.success or result.x[2] <= high:
        return None
    return Fit(float(t0), result.x[:3] + centre, len(time), float(np.median(np.abs(result.fun))))


def _peak(points: np.ndarray, cell: float, smooth: float, extent_max: float) -> tuple[float, float]:
    """Densest spot of 2D points on a grid of ``cell``, smoothed by ``smooth``; the grid covers the central 98 % of the
    points and at most ``extent_max`` metres a side, so stray far points can't blow up its size."""
    lo = np.maximum(np.percentile(points, 1, axis=0), np.median(points, axis=0) - extent_max / 2)
    hi = np.minimum(np.percentile(points, 99, axis=0), lo + extent_max)
    inside = np.all((points >= lo) & (points <= hi), axis=1)
    idx = np.floor((points[inside] - lo) / cell).astype(np.int64)
    shape = np.floor((hi - lo) / cell).astype(np.int64) + 1
    image = np.bincount(idx[:, 0] * shape[1] + idx[:, 1], minlength=int(shape[0] * shape[1])).reshape(shape).astype(float)
    image = ndimage.gaussian_filter(image, smooth / cell)
    i, j = np.unravel_index(np.argmax(image), image.shape)
    return lo[0] + (i + 0.5) * cell, lo[1] + (j + 0.5) * cell


def _near_ground(xyz: np.ndarray, ground_band: float) -> np.ndarray:
    """Points within the ground band above the 5th percentile height per 1 m column."""
    xyz = xyz.astype(np.float64)
    gx, gy = np.floor(xyz[:, 0]).astype(np.int64), np.floor(xyz[:, 1]).astype(np.int64)
    key = (gx - gx.min()) * (gy.max() - gy.min() + 1) + (gy - gy.min())
    order = np.lexsort((xyz[:, 2], key))
    k_sorted = key[order]
    starts = np.flatnonzero(np.r_[True, k_sorted[1:] != k_sorted[:-1]])
    counts = np.diff(np.r_[starts, len(order)])
    ground = np.empty(len(xyz))
    ground[order] = np.repeat(xyz[order, 2][starts + (counts * 0.05).astype(int)], counts)
    return xyz[xyz[:, 2] - ground < ground_band]


def density_concentration(p: Pulses, radius: float = 5.0, max_points: int = 2_000_000,
                          *, group: np.ndarray | None = None) -> float:
    """How strongly near-ground points pile up around one spot: their share within ``radius`` of the density peak,
    over the share an even spread over the covered ground would give. A static scanner's density falls with range,
    so its scans score far higher than a flight or a walk, whose density is even along the path."""
    size = len(p) if group is None else len(group)
    if size > max_points:
        pick = np.random.default_rng(0).choice(size, max_points, replace=False)
        group = pick if group is None else group[pick]
    near = _near_ground(p.first if group is None else p.first[group], 0.5)
    if len(near) < 1000:
        return 0.0
    cx, cy = _peak(near[:, :2], 0.5, 2.0, 500.0)
    share = float((np.hypot(near[:, 0] - cx, near[:, 1] - cy) < radius).mean())
    covered = len(np.unique(np.floor(near[:, :2]).astype(np.int64), axis=0))  # m2 of ground with points
    return share / min(1.0, np.pi * radius**2 / max(covered, 1))


def drop_outliers(fits: list[Fit], window: float = 1.5, limit: float = 1.0) -> list[Fit]:
    """Positions more than ``limit`` metres off the path their neighbours describe: the median position of the
    positions up to ``window`` seconds before and after, interpolated to this position's time. A window fit can lock
    onto a few stray lines (a bird, noise, another pulse's echo) and land metres off; its neighbours don't."""
    if len(fits) < 5:
        return fits
    time = np.array([f.time for f in fits])
    xyz = np.array([f.xyz for f in fits])
    keep = np.ones(len(fits), bool)
    before_start = np.searchsorted(time, time - window, side="left")
    before_end = np.searchsorted(time, time, side="left")
    after_start = np.searchsorted(time, time, side="right")
    after_end = np.searchsorted(time, time + window, side="right")
    for i in np.flatnonzero((before_end - before_start >= 2) & (after_end - after_start >= 2)):
        t = time[i]
        before = slice(before_start[i], before_end[i])
        after = slice(after_start[i], after_end[i])
        tb, ta = np.median(time[before]), np.median(time[after])
        pb, pa = np.median(xyz[before], axis=0), np.median(xyz[after], axis=0)
        expected = pb + (pa - pb) * (t - tb) / (ta - tb)
        keep[i] = np.linalg.norm(xyz[i] - expected) <= limit
    return [f for f, k in zip(fits, keep) if k]
