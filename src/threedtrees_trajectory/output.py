"""Write results as the 3Dtrees trajectory CSV (``3dtrees-trajectory/1``, as the platform's converter writes it) and a
JSON quality report."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import __version__
from .estimate import Result
from .pulses import CloudInfo

FORMAT = "3dtrees-trajectory/1"


def _header(info: CloudInfo, result: Result) -> list[str]:
    time_system = info.time_system if result.platform == "mobile" else "none"
    lines = [f"# format: {FORMAT}", f"# crs: {info.crs}", "# vertical_crs: unknown", f"# time_system: {time_system}"]
    if time_system == "gps_sow":
        lines.append("# gps_week: unknown")
    lines += [f"# source_layout: estimated_{result.method}", f"# source_file: {Path(info.files[0]).name}",
              f"# converter: 3dtrees_trajectory {__version__}"]
    return lines


def write_csv(path: Path, info: CloudInfo, result: Result) -> None:
    rows = _header(info, result)
    if result.platform == "mobile":
        rows.append("time,x,y,z")
        rows += [f"{f.time:.6f},{f.xyz[0]:.4f},{f.xyz[1]:.4f},{f.xyz[2]:.4f}" for f in result.fits]
    else:
        rows.append("name,x,y,z")
        rows += [f"{n},{f.xyz[0]:.4f},{f.xyz[1]:.4f},{f.xyz[2]:.4f}" for n, f in zip(result.names, result.fits)]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def quality(result: Result, miss: np.ndarray) -> str:
    """good: pulse lines meet within 10 cm (cm-accurate on every checked dataset); rough: scan-angle or density
    positions (0.5-8 m); unreliable: pulse lines miss by more than 0.5 m (1815, whose return numbers are broken,
    sits at 6.7 m and 25 m off); none: no positions."""
    if not result.fits:
        return "none"
    if result.method != "lines":
        return "rough"
    median = float(np.nanmedian(miss))
    return "good" if median <= 0.1 else "unreliable" if median > 0.5 else "rough"


def report(info: CloudInfo, result: Result, multi_share: float) -> dict:
    miss = np.array([f.miss for f in result.fits], float)
    out = {
        "tool": f"3dtrees_trajectory {__version__}",
        "files": info.files,
        "points": info.points,
        "kept_every": dict(zip(("multi_return_pulse", "single_return_pulse"), info.kept_every)),
        "crs": info.crs,
        "time_system": info.time_system,
        "platform": result.platform,
        "platform_evidence": result.evidence,
        "method": result.method,
        "multi_return_share": round(multi_share, 4),
        "positions": len(result.fits),
        "median_miss_m": None if not np.isfinite(miss).any() else round(float(np.nanmedian(miss)), 4),
        "quality": quality(result, miss),
        "warnings": result.warnings,
    }
    if result.platform == "mobile" and result.fits:
        t = np.array([f.time for f in result.fits])
        out["time_span"] = [round(float(t.min()), 6), round(float(t.max()), 6)]
        gaps = np.diff(t)
        out["gaps_over_5s"] = int((gaps > 5).sum())
        out["note"] = ("Positions exist only while the cloud has points. Turns and transit legs cut from the cloud "
                       "have no positions.")
    return out


def refusal(info: CloudInfo, pulses, reason: str) -> dict:
    """Report for a cloud the tool can't estimate from: what it found, and why it stopped."""
    return {
        "tool": f"3dtrees_trajectory {__version__}",
        "files": info.files,
        "points": info.points,
        "crs": info.crs,
        "time_system": info.time_system,
        "multi_return_share": round(pulses.multi_share, 4),
        "has_gps_time": pulses.has_time,
        "has_scan_angle": pulses.has_scan_angle,
        "positions": 0,
        "quality": "none",
        "error": reason,
    }


def write_report(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
