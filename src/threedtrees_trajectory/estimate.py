"""Choose the platform and method for a cloud, split static scans, and run the solvers."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import solve
from .pulses import Pulses

MULTI_SHARE_MIN = 0.02  # below this share of multi-return pulses the line method has too little to work with
STATIC_MISS_MAX = 0.05  # one position explaining all lines within 5 cm = a static scan
SCAN_GAP = 30.0  # seconds; a jump in file-order time this large starts a new scan in a merged TLS cloud


@dataclass
class Result:
    platform: str  # mobile | static
    method: str  # lines | scan_angle | density
    fits: list[solve.Fit]
    names: list[str] = field(default_factory=list)  # static: one name per station
    warnings: list[str] = field(default_factory=list)


def split_scans(p: Pulses, by: str) -> list[tuple[str, np.ndarray]]:
    """Index masks, one per scan. by = file | source | time | auto.

    auto, with GPS time: every file and point source cut at time jumps over SCAN_GAP (merged clouds concatenate scans,
    and a source id, where set, may hold several scans). Without time: point sources if there are several (a
    registered project exported as spatial tiles keeps one source per scan), else files (one file per scan), else one
    scan."""
    jump = np.r_[0, np.cumsum(np.abs(np.diff(p.time)) > SCAN_GAP)] if p.has_time else np.zeros(len(p), np.int64)
    if by == "auto":
        if p.has_time:
            key = np.unique(np.c_[p.file, p.source, jump], axis=0, return_inverse=True)[1].ravel()
            return [(f"scan{k}", key == k) for k in np.unique(key)]
        by = "source" if len(np.unique(p.source)) > 1 else "file" if len(np.unique(p.file)) > 1 else "one"
    if by == "file":
        return [(f"file{k}", p.file == k) for k in np.unique(p.file)]
    if by == "source":
        return [(f"source{k}", p.source == k) for k in np.unique(p.source)]
    if by == "time":
        return [(f"scan{k}", jump == k) for k in np.unique(jump)]
    return [("scan0", np.ones(len(p), bool))]


def merge_stations(fits: list[solve.Fit], names: list[str], radius: float = 0.05):
    """One station per position: scans split into several time blocks give the same position more than once."""
    kept_fits, kept_names = [], []
    for fit, name in sorted(zip(fits, names), key=lambda pair: -pair[0].pulses):
        if all(np.linalg.norm(fit.xyz - k.xyz) > radius for k in kept_fits):
            kept_fits.append(fit)
            kept_names.append(name)
    order = np.argsort([names.index(n) for n in kept_names])
    return [kept_fits[i] for i in order], [kept_names[i] for i in order]


def detect_platform(p: Pulses) -> str:
    """Static when one position explains the multi-return lines of the first scan within 5 cm, else mobile."""
    if not p.has_time:
        return "static"
    _, mask = split_scans(p, "auto")[0]
    fit = solve.station_from_lines(p.subset(mask), min_pulses=500)
    if fit is not None and fit.miss <= STATIC_MISS_MAX:
        return "static"
    return "mobile"


def estimate(p: Pulses, platform: str = "auto", split: str = "auto", interval: float = 0.5) -> Result:
    if len(p) == 0:
        raise ValueError("The cloud has no points.")
    multi_share = float(p.multi.mean())
    platform = detect_platform(p) if platform == "auto" else platform
    warnings = []
    if platform == "mobile":
        if not p.has_time:
            raise ValueError("A moving sensor's path needs GPS time per point; this cloud has none.")
        if multi_share >= MULTI_SHARE_MIN:
            fits = solve.path_from_lines(p, interval=interval)
            return Result("mobile", "lines", fits, warnings=warnings)
        if p.has_scan_angle:
            warnings.append("Single returns only: path from scan angles, about 1 m horizontally and 0.1 m vertically on "
                            "RIEGL and DJI flights, much worse where the scan angle isn't measured off nadir (YellowScan 8 m).")
            return Result("mobile", "scan_angle", solve.path_from_scan_angle(p, interval=interval), warnings=warnings)
        raise ValueError("Single returns without scan angle: the path can't be estimated from this cloud. "
                         "Please provide the trajectory file.")
    method = "lines" if p.has_time and multi_share >= MULTI_SHARE_MIN else "density"
    scans = split_scans(p, split)
    if len(scans) == 1 and split == "auto":
        warnings.append("No way to tell scans apart (one file, one point source, no time jumps): treated as one scan.")
    if method == "density":
        warnings.append("Single returns: each scanner position is the peak of near-ground point density, about 0.5-0.7 m "
                        "off on the scans tested; z is the ground under the scanner, not the scanner itself.")
    fits, names = [], []
    for name, mask in scans:
        part = p.subset(mask)
        fit = solve.station_from_lines(part) if method == "lines" else solve.station_from_density(part)
        if fit is None:
            warnings.append(f"{name}: no position (too few usable pulses).")
            continue
        fits.append(fit)
        names.append(name)
    if method == "lines":
        fits, names = merge_stations(fits, names)
    return Result("static", method, fits, names, warnings)
