"""Check that a cloud comes from a moving sensor and run the path solvers. Static (TLS) clouds are refused in 1.0.0:
the output format can hold scan positions, but this version estimates paths only."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import solve
from .pulses import Pulses

MULTI_SHARE_MIN = 0.005  # below this share of multi-return pulses the line method has too little to work with
STATIC_MISS_MAX = 0.05  # one position explaining all lines within 5 cm = a static scan
SCAN_GAP = 30.0  # seconds; a jump in file-order time this large separates scans of a merged TLS cloud
MIN_LINE_POSITIONS = 10  # fewer positions from pulse lines than this: fall back to scan angles
CONCENTRATION_MIN = 15.0  # ground density piling up this much around one spot = a static scanner
STATIC_REFUSAL = ("This looks like a static (TLS) scan ({evidence}). Trajectory 1.0.0 estimates paths of moving "
                  "sensors (ULS, MLS) only.")


@dataclass
class Result:
    method: str  # lines | scan_angle
    fits: list[solve.Fit]
    warnings: list[str] = field(default_factory=list)
    evidence: str = "given"  # why the cloud counts as a moving sensor


def _largest_block(p: Pulses) -> np.ndarray:
    """Pulses of the largest file x point source x time block: one scan, if the cloud is a merged TLS project."""
    block = np.r_[0, np.cumsum(np.abs(np.diff(p.time)) > SCAN_GAP)] if p.has_time else np.zeros(len(p), np.int64)
    label = np.unique(np.c_[p.file, p.source, block], axis=0, return_inverse=True)[1].ravel()
    return label == np.argmax(np.bincount(label))


def detect_platform(p: Pulses) -> tuple[str, str]:
    """mobile or static, and the evidence, from the largest scan block.

    Multi-return pulses with time: static when one position explains their lines within 5 cm. Otherwise: static when
    near-ground points pile up around one spot (CONCENTRATION_MIN times an even spread), as a fixed scanner's do.
    Checked on 3Dtrees data: TLS 22 to 128, ULS 1.8 to 7.7, MLS 4.3."""
    if not p.has_time:
        return "static", "no GPS time, so it can't be a path"
    part = p.subset(_largest_block(p))
    if part.has_time and p.multi_share >= MULTI_SHARE_MIN:
        fit = solve.station_from_lines(part, min_pulses=500)
        if fit is not None:
            platform = "static" if fit.miss <= STATIC_MISS_MAX else "mobile"
            return platform, f"one position explains the pulse lines within {fit.miss:.3f} m"
    concentration = solve.density_concentration(part)
    platform = "static" if concentration >= CONCENTRATION_MIN else "mobile"
    return platform, f"ground density around its peak is {concentration:.1f} times an even spread"


def estimate(p: Pulses, check_platform: bool = True, interval: float = 0.5) -> Result:
    """Path positions in the cloud's coordinates (the solvers work relative to ``p.origin``)."""
    result = _estimate(p, check_platform, interval)
    result.fits = [solve.Fit(f.time, f.xyz + p.origin, f.pulses, f.miss) for f in result.fits]
    return result


def _estimate(p: Pulses, check_platform: bool, interval: float) -> Result:
    if len(p) == 0:
        raise ValueError("The cloud has no points.")
    multi_share = p.multi_share
    warnings = []
    evidence = "given"
    if check_platform:
        platform, evidence = detect_platform(p)
        if platform == "static":
            raise ValueError(STATIC_REFUSAL.format(evidence=evidence))
    if not p.has_time:
        raise ValueError("A moving sensor's path needs GPS time per point; this cloud has none.")
    lines = solve.path_from_lines(p, interval=interval) if multi_share >= MULTI_SHARE_MIN else []
    if len(lines) >= MIN_LINE_POSITIONS:  # a gappy path from lines beats a full one from scan angles (metres off)
        kept = solve.drop_outliers(lines)
        if len(kept) < len(lines):
            warnings.append(f"{len(lines) - len(kept)} of {len(lines)} positions were more than 1 m off the path of "
                            "their neighbours and were dropped.")
        return Result("lines", kept, warnings, evidence)
    if p.has_scan_angle:
        warnings.append("Single returns only: path from scan angles, about 1 m horizontally and 0.1 m vertically on "
                        "RIEGL and DJI flights, much worse where the scan angle isn't measured off nadir (YellowScan 8 m).")
        return Result("scan_angle", solve.path_from_scan_angle(p, interval=interval), warnings, evidence)
    raise ValueError("Single returns without scan angle: the path can't be estimated from this cloud. "
                     "Please provide the trajectory file.")
