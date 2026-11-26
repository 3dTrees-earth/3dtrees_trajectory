"""Write results as the 3Dtrees trajectory CSV (``3dtrees-trajectory/estimated``: the columns delivered paths and scan positions
share; this tool writes paths) and a JSON quality report."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from . import __version__
from .estimate import UNRELIABLE_LINE_MISS, Result
from .publish import write_text
from .pulses import CloudInfo, Pulses

FORMAT = "3dtrees-trajectory/estimated"  # uploaded paths are /1, or /2 when rows have names
COLUMNS = "time,x,y,z,roll,pitch,yaw,name"


TYPICAL_ERROR_M = {"lines": "0.01", "scan_angle": "1 horizontal, 0.1 vertical",
                   "lines_scan_angle": "0.01 for lines and 1 horizontal, 0.1 vertical for scan angles"}


def _header(info: CloudInfo, result: Result, quality_label: str) -> list[str]:
    """The standard's keys, then the estimate's own (``method``, ``quality``)."""
    time_system = info.time_system
    lines = [f"# format: {FORMAT}", "# kind: path", "# origin: estimated", f"# crs: {info.crs}",
             f"# vertical_crs: {info.vertical_crs}", f"# time_system: {time_system}"]
    if time_system == "gps_sow":
        lines.append("# gps_week: unknown")
    lines += ["# attitude: none", f"# source_layout: estimated_{result.method}",
              f"# source_file: {Path(info.files[0]).name}", f"# converter: 3dtrees_trajectory {__version__}",
              f"# method: {result.method}", f"# quality: {quality_label}"]
    return lines


def delivered_miss_flag(check: dict) -> str:
    flag = check["threshold_exceeded"]
    return "unavailable" if flag is None else str(flag).lower()


def write_csv(path: Path, info: CloudInfo, result: Result, *, delivered_check: dict | None = None) -> None:
    """``time,x,y,z,roll,pitch,yaw,name``, the layout delivered paths and scan positions share; the name comes last so
    readers that take time, x, y, z by position are unaffected. A path leaves the name empty. Roll, pitch and yaw stay
    empty: estimation gives positions only. Speed is not stored (viewers derive it), nor how well the pulse lines meet
    (in the debug report)."""
    rows = _header(info, result, estimate_quality(result))
    for warning in result.warnings:
        rows.append("# warning: " + " ".join(warning.splitlines()))
    if delivered_check is not None:
        rows.append(f"# delivered_path_miss_flag: {delivered_miss_flag(delivered_check)}")
        if "path_files" in delivered_check:
            rows.append(f"# delivered_path_files: {delivered_check['path_files']}")
            rows.append(f"# delivered_path_coverage: {delivered_check['coverage']}")
        for warning in delivered_check.get("warnings", ()):
            rows.append("# warning: " + " ".join(warning.splitlines()))
    rows.append(COLUMNS)
    timed = info.time_system != "none"
    for f in result.fits:
        x, y, z = f.xyz
        rows.append(f"{f'{f.time:.6f}' if timed else ''},{x:.4f},{y:.4f},{z:.4f},,,,")
    write_text(path, "\n".join(rows) + "\n")


def quality(result: Result, miss: np.ndarray) -> str:
    """good: pulse lines meet within 10 cm (cm-accurate on every checked dataset); rough: 0.1-0.5 m, or scan-angle
    positions (1-8 m); unreliable: pulse lines miss by more than 0.5 m (1815, whose return numbers are broken,
    sits at 6.7 m and 25 m off); none: no positions."""
    if not result.fits:
        return "none"
    if result.method != "lines":
        return "rough"
    median = float(np.nanmedian(miss))
    return "good" if median <= 0.1 else "unreliable" if median > UNRELIABLE_LINE_MISS else "rough"


def estimate_quality(result: Result) -> str:
    return quality(result, np.array([f.miss for f in result.fits], float))


def report(info: CloudInfo, result: Result, multi_share: float) -> dict:
    miss = np.array([f.miss for f in result.fits], float)
    out = {
        "tool": f"3dtrees_trajectory {__version__}",
        "files": info.files,
        "points": info.points,
        "kept_every": dict(zip(("multi_return_pulse", "single_return_pulse"), info.kept_every)),
        "crs": info.crs,
        "time_system": info.time_system,
        "platform_evidence": result.evidence,
        "method": result.method,
        "multi_return_share": round(multi_share, 4),
        "positions": len(result.fits),
        "median_miss_m": None if not np.isfinite(miss).any() else round(float(np.nanmedian(miss)), 4),
        "typical_error_m": TYPICAL_ERROR_M[result.method],
        "quality": quality(result, miss),
        "warnings": result.warnings,
    }
    if result.fits:
        t = np.array([f.time for f in result.fits])
        out["time_span"] = [round(float(t.min()), 6), round(float(t.max()), 6)]
        gaps = np.diff(t)
        out["gaps_over_5s"] = int((gaps > 5).sum())
        out["note"] = ("Positions exist only while the cloud has points. Turns and transit legs cut from the cloud "
                       "have no positions.")
    return out


def refusal(info: CloudInfo, pulses: Pulses | None, reason: str) -> dict:
    """Report for a cloud the tool can't estimate from: what it found, and why it stopped."""
    return {
        "tool": f"3dtrees_trajectory {__version__}",
        "files": info.files,
        "points": info.points,
        "crs": info.crs,
        "time_system": info.time_system,
        "multi_return_share": round(pulses.multi_share, 4) if pulses is not None else None,
        "has_gps_time": pulses.has_time if pulses is not None else info.has_gps_time,
        "has_scan_angle": pulses.has_scan_angle if pulses is not None else None,
        "positions": 0,
        "pairing": {
            "complete_multi_pulses": info.multi_pulses,
            "incomplete_multi_pulses": info.incomplete_multi_pulses,
            "out_of_order_records_sampled": info.joined_out_of_order_points,
            "unpaired_multi_echoes": info.unpaired_multi_echoes,
            "return_modes": info.return_modes,
        },
        "quality": "none",
        "error": reason,
    }
