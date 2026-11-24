"""Write results as the 3Dtrees trajectory CSV (``3dtrees-trajectory/estimated``: the columns delivered paths and scan positions
share; this tool writes paths) and a JSON quality report."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import __version__
from .estimate import UNRELIABLE_LINE_MISS, Result
from .publish import write_text
from .pulses import LOCAL, CloudInfo, Pulses

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


def geojson(info: CloudInfo, result: Result, gap: float = 5.0, spacing: float = 1.0) -> dict:
    """A light map layer: the path as lines split where positions are more than ``gap`` seconds apart (no line drawn
    across unscanned turns), thinned to about one vertex per ``spacing`` metres. WGS84 when
    the cloud's CRS is known, else the cloud's own coordinates with ``crs`` set to ``local`` in the properties."""
    xyz = np.array([f.xyz for f in result.fits], float).reshape(-1, 3)
    lonlat = xyz[:, :2]
    known = info.crs not in (LOCAL, "unknown") and len(xyz)
    if known:
        from pyproj import Transformer
        lon, lat = Transformer.from_crs(info.crs, "EPSG:4326", always_xy=True).transform(xyz[:, 0], xyz[:, 1], errcheck=True)
        lonlat = np.c_[lon, lat]
    digits = 7 if known else 3
    properties = {"estimated": True, "method": result.method,
                  "quality": quality(result, np.array([f.miss for f in result.fits], float)),
                  "crs": "EPSG:4326" if known else LOCAL, "source_crs": info.crs}
    features = []
    if len(xyz):
        time = np.array([f.time for f in result.fits])
        lines = []
        for part in np.split(np.arange(len(time)), np.flatnonzero(np.diff(time) > gap) + 1):
            keep, last = [part[0]], xyz[part[0]]
            for i in part[1:]:
                if np.linalg.norm(xyz[i] - last) >= spacing or i == part[-1]:
                    keep.append(i)
                    last = xyz[i]
            if len(keep) > 1:
                lines.append([[round(float(a), digits), round(float(b), digits), round(float(z), 2)]
                              for (a, b), z in zip(lonlat[keep], xyz[keep, 2])])
        if lines:
            features.append({"type": "Feature", "geometry": {"type": "MultiLineString", "coordinates": lines},
                             "properties": properties | {"time_start": float(time.min()), "time_end": float(time.max()),
                                                         "time_system": info.time_system, "segments": len(lines)}})
    return {"type": "FeatureCollection", "features": features}


def write_geojson(path: Path, data: dict) -> None:
    write_text(path, json.dumps(data, separators=(",", ":"), allow_nan=False) + "\n")


def geojson_from_csv(path: Path) -> dict:
    """The map layer for a path CSV this tool wrote (header keys crs, kind, method, time_system)."""
    from .solve import Fit

    meta, rows = {}, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            key, _, value = line[1:].partition(":")
            meta[key.strip()] = value.strip()
        elif line and not line.startswith("time,"):
            rows.append(line.split(","))
    if meta.get("kind", "path") != "path":
        raise ValueError(f"{path}: kind {meta['kind']!r}; this version draws paths only")
    fits = [Fit(float(r[0]), np.array(r[1:4], float), 0, float("nan")) for r in rows]
    result = Result(meta.get("method", "lines"), fits)
    crs = meta.get("crs", LOCAL)
    info = CloudInfo(LOCAL if crs in ("", "unknown") else crs, meta.get("time_system", "none"), 0, [meta.get("source_file", "")], (1, 1))
    return geojson(info, result)


def main_geojson(argv: list[str] | None = None) -> int:
    """3dtrees-trajectory-geojson: write <name>.geojson next to each given trajectory CSV."""
    import sys
    for name in (argv if argv is not None else sys.argv[1:]):
        path = Path(name)
        write_geojson(path.with_suffix(".geojson"), geojson_from_csv(path))
        print(path.with_suffix(".geojson"))
    return 0
