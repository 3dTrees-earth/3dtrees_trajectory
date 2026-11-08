"""Sensor positions from pulses.

lines       A pulse's first and last return lie on a line through the sensor. The point (and, for a moving sensor, the
            velocity) closest to all lines of a time window is the sensor position. Lines that miss it are dropped and
            the fit repeated. Needs multi-return pulses; centimetre accuracy on the datasets tested.
scan_angle  Single returns with a scan angle: the sensor sits where every point is seen at its scan angle off nadir
            (a cone around each point). Level flight and a line scanner are assumed; about 1 m horizontally and
            0.1 m vertically on the RIEGL and DJI flights tested, much worse on YellowScan (8 m).
density     Single returns of one static scan: the scanner stands where ground points are densest (they thin out
            with range). Rough, 0.5-0.7 m on the FARO and BLK360 scans tested, and only inside the scanned area.
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
    line = p.first - p.last
    length = np.linalg.norm(line, axis=1)
    keep = p.multi & (length >= min_length)
    return p.first[keep], line[keep] / length[keep, None], p.time[keep]


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
    source = p.source[p.multi & (np.linalg.norm(p.first - p.last, axis=1) >= min_length)]
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
    return sorted(fits, key=lambda f: f.time)


def station_from_lines(p: Pulses, min_pulses: int = 2000, min_length: float = 0.5) -> Fit | None:
    """Static scanner: one position from all lines of one scan."""
    origin, direction, time = _lines(p, min_length)
    if len(origin) < min_pulses:
        return None
    return fit_lines(origin, direction, time, False, min_pulses)


def path_from_scan_angle(p: Pulses, interval: float = 0.5, min_pulses: int = 200) -> list[Fit]:
    """Moving sensor, single returns: per window, the position from which each point is seen at its scan angle."""
    xyz, time = p.first, p.time
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


def station_from_density(p: Pulses, cell: float = 0.1, ground_band: float = 0.5, smooth: float = 1.0,
                         max_points: int = 5_000_000) -> Fit | None:
    """Static scanner, single returns: the peak of near-ground point density (density falls with range).

    Coarse peak on a 0.5 m grid first, then the fine peak on a ``cell`` grid within 10 m of it; at most ``max_points``
    points (an even sample) are used, so memory stays bounded for any scan size."""
    xyz = p.first
    if len(xyz) < 1000:
        return None
    if len(xyz) > max_points:
        xyz = xyz[np.random.default_rng(0).choice(len(xyz), max_points, replace=False)]
    # local ground: 5th percentile height per 1 m column; near-ground = within ground_band above it
    gx, gy = np.floor(xyz[:, 0]).astype(np.int64), np.floor(xyz[:, 1]).astype(np.int64)
    key = (gx - gx.min()) * (gy.max() - gy.min() + 1) + (gy - gy.min())
    order = np.lexsort((xyz[:, 2], key))
    k_sorted = key[order]
    starts = np.flatnonzero(np.r_[True, k_sorted[1:] != k_sorted[:-1]])
    counts = np.diff(np.r_[starts, len(order)])
    ground = np.empty(len(xyz))
    ground[order] = np.repeat(xyz[order, 2][starts + (counts * 0.05).astype(int)], counts)
    near = xyz[xyz[:, 2] - ground < ground_band]
    if len(near) < 1000:
        return None
    cx, cy = _peak(near[:, :2], 0.5, 2.0, 500.0)
    close = near[np.hypot(near[:, 0] - cx, near[:, 1] - cy) < 10]
    x, y = _peak(close[:, :2], cell, smooth, 20.0)  # smoothing bridges the blind spot under the tripod
    z = float(np.median(near[np.hypot(near[:, 0] - x, near[:, 1] - y) < 2, 2]))
    return Fit(float(np.median(p.time)) if p.has_time else 0.0, np.array([x, y, z]), len(near), float("nan"))
