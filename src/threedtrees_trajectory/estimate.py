"""Check that a cloud comes from a moving sensor and run the path solvers."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import solve
from .coverage import NoTrajectory, UnusableTime
from .pulses import Pulses

MULTI_SHARE_MIN = 0.005  # below this share of multi-return pulses the line method has too little to work with
STATIC_MISS_MAX = 0.05  # one position explaining all lines within 5 cm = a static scan
SCAN_GAP = 30.0  # seconds; a jump in file-order time this large separates scans of a merged TLS cloud
UNRELIABLE_LINE_MISS = 0.5
MIN_LINE_POSITIONS = 10  # fewer positions from pulse lines than this: fall back to scan angles
CONCENTRATION_MIN = 15.0  # ground density piling up this much around one spot = a static scanner
STATIC_REFUSAL = ("This looks like a static (TLS) scan ({evidence}). This tool estimates paths of moving "
                  "sensors (ULS and MLS).")


@dataclass
class Result:
    method: str  # lines | scan_angle | lines_scan_angle
    fits: list[solve.Fit]
    warnings: list[str] = field(default_factory=list)
    evidence: str = "given"  # why the cloud counts as a moving sensor


def _largest_block(p: Pulses) -> np.ndarray | None:
    """Indices of the largest file, point source and time block, or None for all pulses."""
    jumps = np.abs(np.diff(p.time)) > SCAN_GAP
    files_change = np.any(p.file != p.file[0])
    sources_change = np.any(p.source != p.source[0])
    time_jumps = np.any(jumps)
    if not (files_change or sources_change or time_jumps):
        return None
    keys = []
    if time_jumps:
        keys.append(np.r_[0, np.cumsum(jumps)])
    if sources_change:
        keys.append(p.source)
    if files_change:
        keys.append(p.file)
    order = np.lexsort(keys)
    boundary = np.zeros(len(p) - 1, dtype=bool)
    for key in keys:
        ordered = key[order]
        boundary |= ordered[1:] != ordered[:-1]
    starts = np.r_[0, np.flatnonzero(boundary) + 1, len(p)]
    largest = np.argmax(np.diff(starts))
    return order[starts[largest]:starts[largest + 1]].copy()


def detect_platform(p: Pulses) -> tuple[str, str]:
    """mobile or static, and the evidence, from the largest scan block.

    Multi-return pulses with time: static when one position explains their lines within 5 cm. Otherwise: static when
    near-ground points pile up around one spot (CONCENTRATION_MIN times an even spread), as a fixed scanner's do.
    Checked on 3Dtrees data: TLS 22 to 128, ULS 1.8 to 7.7, MLS 4.3."""
    block = _largest_block(p)
    if p.multi_share >= MULTI_SHARE_MIN:
        fit = solve.station_from_lines(p, min_pulses=500, group=block)
        if fit is not None:
            platform = "static" if fit.miss <= STATIC_MISS_MAX else "mobile"
            return platform, f"one position explains the pulse lines within {fit.miss:.3f} m"
    concentration = solve.density_concentration(p, group=block)
    platform = "static" if concentration >= CONCENTRATION_MIN else "mobile"
    return platform, f"ground density around its peak is {concentration:.1f} times an even spread"


def estimate(p: Pulses, check_platform: bool = True, interval: float = 0.5, *, use_lines: bool = True,
             angle_files=(), workers: int = 1) -> Result:
    """Path positions in the cloud's coordinates (the solvers work relative to ``p.origin``)."""
    result = _estimate(p, check_platform, interval, use_lines, workers)
    if (result.method == "lines" and angle_files
            and np.median([fit.miss for fit in result.fits]) <= UNRELIABLE_LINE_MISS):
        result = _fill_angle_windows(p, result, interval, angle_files, workers)
    return world_result(result, p.origin)


def world_result(result, origin):
    if not result.fits:
        raise NoTrajectory("No reliable position could be fitted. Please provide the trajectory file.")
    result.fits = [solve.Fit(f.time, f.xyz + origin, f.pulses, f.miss) for f in result.fits]
    if not all(np.isfinite(f.time) and np.isfinite(f.xyz).all() for f in result.fits):
        raise ValueError("The estimated path contains non-finite coordinates or times.")
    return result


def _fill_angle_windows(p, result, interval, angle_files, workers):
    """Use paired pulse lines where available and fit missing windows from the affected file."""
    by_window = {int(np.floor(fit.time / interval)): fit for fit in result.fits}
    line_windows = set(by_window)
    for file in angle_files:
        group = np.flatnonzero(p.file == file)
        if not len(group):
            continue
        windows = np.unique(np.floor(p.time[group] / interval).astype(np.int64))
        if all(window in by_window for window in windows):
            continue
        for fit in solve.path_from_scan_angle(p, interval=interval, group=group, workers=workers):
            by_window.setdefault(int(np.floor(fit.time / interval)), fit)
        if not any(window in by_window for window in windows):
            raise NoTrajectory("No reliable scan-angle position could be fitted for a file with unpaired returns. "
                               "Please provide the delivered trajectory file.")
    extra = len(by_window) - len(line_windows)
    if extra:
        result.method = "lines_scan_angle"
        result.fits = [by_window[window] for window in sorted(by_window)]
        result.warnings.append(f"Scan angles supplied {extra} positions in windows without a pulse-line fit.")
    return result


def _estimate(p: Pulses, check_platform: bool, interval: float, use_lines: bool, workers: int) -> Result:
    if len(p) == 0:
        raise ValueError("The cloud has no points.")
    if not p.has_time:
        raise UnusableTime("A moving sensor's path needs GPS time per point; this cloud has none.")
    if not np.isfinite(interval) or interval <= 0:
        raise ValueError("The position interval must be a positive finite number of seconds.")
    if max(abs(p.time.min()), abs(p.time.max())) / np.iinfo(np.int64).max >= interval:
        raise ValueError("The position interval is too small for this cloud's GPS times.")
    multi_share = p.multi_share
    warnings = []
    evidence = "given"
    if check_platform:
        platform, evidence = detect_platform(p)
        if platform == "static":
            raise ValueError(STATIC_REFUSAL.format(evidence=evidence))
    lines = solve.path_from_lines(p, interval=interval, workers=workers) if use_lines and multi_share >= MULTI_SHARE_MIN else []
    if len(lines) >= MIN_LINE_POSITIONS:  # a gappy path from lines beats a full one from scan angles (metres off)
        return line_result(lines, evidence)
    if p.has_scan_angle:
        warnings.append("Single returns only: path from scan angles, about 1 m horizontally and 0.1 m vertically on "
                        "RIEGL and DJI flights, much worse where the scan angle isn't measured off nadir (YellowScan 8 m).")
        return Result("scan_angle", solve.path_from_scan_angle(p, interval=interval, workers=workers), warnings, evidence)
    if multi_share >= MULTI_SHARE_MIN:
        raise NoTrajectory("Too few reliable positions from multi-return pulses and no usable scan angle. "
                         "Please provide the trajectory file.")
    raise NoTrajectory("Single returns without scan angle: the path can't be estimated from this cloud. "
                     "Please provide the trajectory file.")


def line_result(lines, evidence):
    kept = solve.drop_outliers(lines)
    warnings = []
    if len(kept) < len(lines):
        warnings.append(f"{len(lines) - len(kept)} of {len(lines)} positions were more than 1 m off the path of "
                        "their neighbours and were dropped.")
    return Result("lines", kept, warnings, evidence)
