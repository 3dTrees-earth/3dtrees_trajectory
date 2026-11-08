"""Choose the platform and method for a cloud, split static scans, and run the solvers."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import solve
from .pulses import Pulses

MULTI_SHARE_MIN = 0.02  # below this share of multi-return pulses the line method has too little to work with
STATIC_MISS_MAX = 0.05  # one position explaining all lines within 5 cm = a static scan
SCAN_GAP = 30.0  # seconds; a jump in file-order time this large starts a new scan in a merged TLS cloud
CONCENTRATION_MIN = 15.0  # ground density piling up this much around one spot = a static scanner


@dataclass
class Result:
    platform: str  # mobile | static
    method: str  # lines | scan_angle | density
    fits: list[solve.Fit]
    names: list[str] = field(default_factory=list)  # static: one name per station
    warnings: list[str] = field(default_factory=list)
    evidence: str = "given"  # why this platform


MIN_SCAN_PULSES = 2000  # smaller groups are fragments, not scans


def _time_blocks(time: np.ndarray) -> np.ndarray:
    """Scan block per pulse. Merged clouds concatenate scans, so a jump of SCAN_GAP in file order starts a new one.
    When that leaves over 10 % of pulses in crumbs (points stored out of acquisition order), blocks are gaps in sorted
    time instead."""
    block = np.r_[0, np.cumsum(np.abs(np.diff(time)) > SCAN_GAP)]
    sizes = np.bincount(block)
    if len(time) and sizes[sizes >= MIN_SCAN_PULSES].sum() >= 0.9 * len(time):
        return block
    order = np.argsort(time, kind="stable")
    sorted_block = np.r_[0, np.cumsum(np.diff(time[order]) > SCAN_GAP)]
    block = np.empty(len(time), np.int64)
    block[order] = sorted_block
    return block


def split_scans(p: Pulses, by: str) -> tuple[list[str], np.ndarray]:
    """Scan names and one label per pulse (index into names). by = file | source | time | auto.

    auto, with GPS time: every file and point source cut into time blocks (a source id, where set, may hold several
    scans). Without time: point sources if there are several (a registered project exported as spatial tiles keeps one
    source per scan), else files (one file per scan), else one scan."""
    if by == "auto":
        if p.has_time:
            keys = np.c_[p.file, p.source, _time_blocks(p.time)]
            label = np.unique(keys, axis=0, return_inverse=True)[1].ravel()
            return [f"scan{k}" for k in range(label.max() + 1)], label
        by = "source" if len(np.unique(p.source)) > 1 else "file" if len(np.unique(p.file)) > 1 else "one"
    if by in ("file", "source"):
        values, label = np.unique(p.file if by == "file" else p.source, return_inverse=True)
        return [f"{by}{v}" for v in values], label.ravel()
    if by == "time" and p.has_time:
        label = _time_blocks(p.time)
        return [f"scan{k}" for k in range(label.max() + 1)], label
    return ["scan0"], np.zeros(len(p), np.int64)


def scan_groups(p: Pulses, by: str):
    """(name, pulse indices) per scan with at least MIN_SCAN_PULSES pulses, and the number of smaller fragments."""
    names, label = split_scans(p, by)
    order = np.argsort(label, kind="stable")
    bounds = np.r_[0, np.cumsum(np.bincount(label, minlength=len(names)))]
    groups = [(names[k], order[bounds[k]:bounds[k + 1]]) for k in range(len(names))]
    big = [g for g in groups if len(g[1]) >= MIN_SCAN_PULSES]
    return big, len(groups) - len(big)


def merge_stations(fits: list[solve.Fit], names: list[str], radius: float = 0.05):
    """One station per position: scans split into several time blocks give the same position more than once."""
    kept_fits, kept_names = [], []
    for fit, name in sorted(zip(fits, names), key=lambda pair: -pair[0].pulses):
        if all(np.linalg.norm(fit.xyz - k.xyz) > radius for k in kept_fits):
            kept_fits.append(fit)
            kept_names.append(name)
    order = np.argsort([names.index(n) for n in kept_names])
    return [kept_fits[i] for i in order], [kept_names[i] for i in order]


def detect_platform(p: Pulses) -> tuple[str, str]:
    """mobile or static, and the evidence, from the largest scan group.

    Multi-return pulses with time: static when one position explains their lines within 5 cm. Otherwise: static when
    near-ground points pile up around one spot (CONCENTRATION_MIN times an even spread), as a fixed scanner's do.
    Checked on 3Dtrees data: TLS 22 to 128, ULS 1.8 to 7.7, MLS 4.3."""
    names, label = split_scans(p, "auto")
    part = p.subset(label == np.argmax(np.bincount(label)))
    if part.has_time and float(part.multi.mean()) >= MULTI_SHARE_MIN:
        fit = solve.station_from_lines(part, min_pulses=500)
        if fit is not None:
            platform = "static" if fit.miss <= STATIC_MISS_MAX else "mobile"
            return platform, f"one position explains the pulse lines within {fit.miss:.3f} m"
    concentration = solve.density_concentration(part)
    platform = "static" if concentration >= CONCENTRATION_MIN else "mobile"
    return platform, f"ground density around its peak is {concentration:.1f} times an even spread"


def estimate(p: Pulses, platform: str = "auto", split: str = "auto", interval: float = 0.5) -> Result:
    if len(p) == 0:
        raise ValueError("The cloud has no points.")
    multi_share = float(p.multi.mean())
    warnings = []
    evidence = "given"
    if platform == "auto":
        platform, evidence = detect_platform(p)
    if platform == "mobile":
        if not p.has_time:
            raise ValueError("A moving sensor's path needs GPS time per point; this cloud has none.")
        if multi_share >= MULTI_SHARE_MIN:
            fits = solve.path_from_lines(p, interval=interval)
            return Result("mobile", "lines", fits, warnings=warnings, evidence=evidence)
        if p.has_scan_angle:
            warnings.append("Single returns only: path from scan angles, about 1 m horizontally and 0.1 m vertically on "
                            "RIEGL and DJI flights, much worse where the scan angle isn't measured off nadir (YellowScan 8 m).")
            return Result("mobile", "scan_angle", solve.path_from_scan_angle(p, interval=interval), warnings=warnings,
                          evidence=evidence)
        raise ValueError("Single returns without scan angle: the path can't be estimated from this cloud. "
                         "Please provide the trajectory file.")
    method = "lines" if p.has_time and multi_share >= MULTI_SHARE_MIN else "density"
    scans, fragments = scan_groups(p, split)
    if len(scans) == 1 and split == "auto":
        warnings.append("No way to tell scans apart (one file, one point source, no time jumps): treated as one scan.")
    if fragments:
        warnings.append(f"{fragments} fragments under {MIN_SCAN_PULSES} pulses were left out.")
    if method == "density":
        warnings.append("Single returns: each scanner position is the peak of near-ground point density, about 0.5-0.7 m "
                        "off on the scans tested; z is the ground under the scanner, not the scanner itself.")
    fits, names = [], []
    for name, index in scans:
        part = p.subset(index)
        fit = solve.station_from_lines(part) if method == "lines" else solve.station_from_density(part)
        if fit is None:
            warnings.append(f"{name}: no position (too few usable pulses).")
            continue
        fits.append(fit)
        names.append(name)
    if method == "lines":
        fits, names = merge_stations(fits, names)
    return Result("static", method, fits, names, warnings, evidence)
