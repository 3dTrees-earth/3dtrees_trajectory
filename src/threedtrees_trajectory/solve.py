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

from dataclasses import dataclass

import numpy as np
from scipy import ndimage
from scipy.optimize import least_squares

from .pulses import Pulses


@dataclass
class Fit:
    time: float
    xyz: np.ndarray
    pulses: int
    miss: float  # median distance of the used pulses from the fit, metres (method-specific)


def _lines(p: Pulses, min_length: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    first, last = p.first.astype(np.float64), p.last.astype(np.float64)
    line = first - last
    length = np.linalg.norm(line, axis=1)
    keep = p.multi & (length >= min_length)
    return first[keep], line[keep] / length[keep, None], p.time[keep]


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


def path_from_lines(p: Pulses, interval: float = 0.5, min_pulses: int = 50, min_length: float = 1.0) -> list[Fit]:
    """Moving sensor: one position (with velocity) per time window and point source."""
    origin, direction, time = _lines(p, min_length)
    source = p.source[p.multi & (np.linalg.norm(p.first.astype(np.float64) - p.last, axis=1) >= min_length)]
    window = np.floor(time / interval).astype(np.int64)
    order = np.lexsort((source, window))
    origin, direction, time, source, window = origin[order], direction[order], time[order], source[order], window[order]
    starts = np.flatnonzero(np.r_[True, (window[1:] != window[:-1]) | (source[1:] != source[:-1])])
    fits = []
    for lo, hi in zip(starts, np.r_[starts[1:], len(time)]):
        if hi - lo < min_pulses:
            continue
        fit = fit_lines(origin[lo:hi], direction[lo:hi], time[lo:hi], True, min_pulses)
        if fit is not None:
            fits.append(fit)
    # A window whose lines don't meet (median miss 0.7-30 m against 0.02 m) fitted noise or mixed pulses: drop it.
    if fits:
        limit = max(0.25, 10 * float(np.median([f.miss for f in fits])))
        fits = [f for f in fits if f.miss <= limit]
    return sorted(fits, key=lambda f: f.time)


def station_from_lines(p: Pulses, min_pulses: int = 2000, min_length: float = 0.5, max_pulses: int = 200_000) -> Fit | None:
    """One fixed position from the lines of one scan, to recognise a static scanner (an even sample of at most ``max_pulses``; more add
    nothing to a fit that is already good to about 1 cm)."""
    origin, direction, time = _lines(p, min_length)
    if len(origin) < min_pulses:
        return None
    if len(origin) > max_pulses:
        pick = np.random.default_rng(0).choice(len(origin), max_pulses, replace=False)
        origin, direction, time = origin[pick], direction[pick], time[pick]
    return fit_lines(origin, direction, time, False, min_pulses)


def path_from_scan_angle(p: Pulses, interval: float = 0.5, min_pulses: int = 200) -> list[Fit]:
    """Moving sensor, single returns: per window, the position from which each point is seen at its scan angle."""
    xyz, time = p.first.astype(np.float64), p.time
    tangent = np.tan(np.radians(np.abs(p.scan_angle.astype(float))))
    centre = xyz.mean(axis=0)
    rel = xyz - centre
    order = np.argsort(time, kind="stable")
    rel, time, tangent = rel[order], time[order], tangent[order]
    window = np.floor(time / interval).astype(np.int64)
    starts = np.flatnonzero(np.r_[True, window[1:] != window[:-1]])
    fits = []
    for lo, hi in zip(starts, np.r_[starts[1:], len(time)]):
        if hi - lo < min_pulses:
            continue
        x, y, z, tn = rel[lo:hi, 0], rel[lo:hi, 1], rel[lo:hi, 2], tangent[lo:hi]
        dt = time[lo:hi] - time[lo:hi].mean()

        def residual(q, x=x, y=y, z=z, tn=tn, dt=dt):
            return np.hypot(x - q[0] - q[3] * dt, y - q[1] - q[4] * dt) - (q[2] - z) * tn

        start = [np.median(x), np.median(y), np.percentile(z, 95) + 50, 0.0, 0.0]
        result = least_squares(residual, start, loss="soft_l1", f_scale=0.5)
        if not result.success or result.x[2] <= np.percentile(z, 95):
            continue
        fits.append(Fit(float(time[lo:hi].mean()), result.x[:3] + centre, int(hi - lo), float(np.median(np.abs(result.fun)))))
    return fits


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


def _near_ground(p: Pulses, ground_band: float, max_points: int) -> np.ndarray:
    """Points within ``ground_band`` above the local ground (5th percentile height per 1 m column), from an even
    sample of at most ``max_points`` points."""
    xyz = p.first
    if len(xyz) > max_points:
        xyz = xyz[np.random.default_rng(0).choice(len(xyz), max_points, replace=False)]
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


def density_concentration(p: Pulses, radius: float = 5.0, max_points: int = 2_000_000) -> float:
    """How strongly near-ground points pile up around one spot: their share within ``radius`` of the density peak,
    over the share an even spread over the covered ground would give. A static scanner's density falls with range,
    so its scans score far higher than a flight or a walk, whose density is even along the path."""
    near = _near_ground(p, 0.5, max_points)
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
    for i, t in enumerate(time):
        before = (time >= t - window) & (time < t)
        after = (time > t) & (time <= t + window)
        if before.sum() < 2 or after.sum() < 2:
            continue
        tb, ta = np.median(time[before]), np.median(time[after])
        pb, pa = np.median(xyz[before], axis=0), np.median(xyz[after], axis=0)
        expected = pb + (pa - pb) * (t - tb) / (ta - tb)
        keep[i] = np.linalg.norm(xyz[i] - expected) <= limit
    return [f for f, k in zip(fits, keep) if k]
