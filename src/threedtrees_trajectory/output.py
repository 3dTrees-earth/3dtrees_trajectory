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


TYPICAL_ERROR_M = {"lines": "0.01", "scan_angle": "1 horizontal, 0.1 vertical", "density": "0.5 horizontal"}


def _header(info: CloudInfo, result: Result, quality_label: str, median_miss: float | None) -> list[str]:
    """The converter's keys, then the estimate's own (extra ``#`` lines; the platform's reader skips unknown keys)."""
    time_system = info.time_system if result.platform == "mobile" else "none"
    lines = [f"# format: {FORMAT}", f"# crs: {info.crs}", "# vertical_crs: unknown", f"# time_system: {time_system}"]
    if time_system == "gps_sow":
        lines.append("# gps_week: unknown")
    lines += [f"# source_layout: estimated_{result.method}", f"# source_file: {Path(info.files[0]).name}",
              f"# converter: 3dtrees_trajectory {__version__}",
              "# estimated: true", f"# method: {result.method}", f"# quality: {quality_label}",
              f"# typical_error_m: {TYPICAL_ERROR_M[result.method]}"]
    if median_miss is not None:
        lines.append(f"# median_miss_m: {median_miss:.4f}")
    if result.method == "density":
        lines.append("# z: ground under the scanner")
    return lines


def write_csv(path: Path, info: CloudInfo, result: Result) -> None:
    """Columns exactly as the standard (time or name, x, y, z), so the file reads like a delivered one. How well the
    pulse lines meet is only summarised in the header (``median_miss_m``)."""
    miss = np.array([f.miss for f in result.fits], float)
    median = float(np.nanmedian(miss)) if np.isfinite(miss).any() else None
    rows = _header(info, result, quality(result, miss), median)
    rows.append(("time" if result.platform == "mobile" else "name") + ",x,y,z")
    labels = [f"{f.time:.6f}" for f in result.fits] if result.platform == "mobile" else result.names
    for label, f in zip(labels, result.fits):
        rows.append(f"{label},{f.xyz[0]:.4f},{f.xyz[1]:.4f},{f.xyz[2]:.4f}")
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


def geojson(info: CloudInfo, result: Result, gap: float = 5.0, spacing: float = 1.0) -> dict:
    """A light map layer: the path as lines split where positions are more than ``gap`` seconds apart (no line drawn
    across unscanned turns), thinned to about one vertex per ``spacing`` metres; TLS stations as points. WGS84 when
    the cloud's CRS is known, else the cloud's own coordinates with ``crs`` set to unknown in the properties."""
    xyz = np.array([f.xyz for f in result.fits], float).reshape(-1, 3)
    lonlat = xyz[:, :2]
    known = info.crs != "unknown" and len(xyz)
    if known:
        from pyproj import Transformer
        lon, lat = Transformer.from_crs(info.crs, "EPSG:4326", always_xy=True).transform(xyz[:, 0], xyz[:, 1])
        lonlat = np.c_[lon, lat]
    digits = 7 if known else 3
    properties = {"estimated": True, "method": result.method,
                  "quality": quality(result, np.array([f.miss for f in result.fits], float)),
                  "crs": "EPSG:4326" if known else "unknown (cloud coordinates)", "source_crs": info.crs}
    features = []
    if result.platform == "mobile" and len(xyz):
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
        features.append({"type": "Feature", "geometry": {"type": "MultiLineString", "coordinates": lines},
                         "properties": properties | {"time_start": float(time.min()), "time_end": float(time.max()),
                                                     "time_system": info.time_system, "segments": len(lines)}})
    else:
        for name, (a, b), z in zip(result.names, lonlat, xyz[:, 2]):
            features.append({"type": "Feature", "geometry": {"type": "Point", "coordinates": [
                round(float(a), digits), round(float(b), digits), round(float(z), 2)]}, "properties": properties | {"name": name}})
    return {"type": "FeatureCollection", "features": features}


def write_geojson(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, separators=(",", ":")) + "\n", encoding="utf-8")


def geojson_from_csv(path: Path) -> dict:
    """The map layer for a trajectory CSV this tool wrote (header keys crs, method, quality, time_system)."""
    from .solve import Fit

    meta, rows = {}, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            key, _, value = line[1:].partition(":")
            meta[key.strip()] = value.strip()
        elif line and not line.startswith(("time,", "name,")):
            rows.append(line.split(","))
    mobile = path.read_text(encoding="utf-8").count("\ntime,") > 0
    fits = [Fit(float(r[0]) if mobile else 0.0, np.array(r[1:4], float), 0,
                float(r[4]) if len(r) > 4 else float("nan")) for r in rows]
    result = Result("mobile" if mobile else "static", meta.get("method", "lines"), fits,
                    [] if mobile else [r[0] for r in rows])
    info = CloudInfo(meta.get("crs", "unknown"), meta.get("time_system", "none"), 0, [meta.get("source_file", "")], (1, 1))
    return geojson(info, result)


def main_geojson(argv: list[str] | None = None) -> int:
    """3dtrees-trajectory-geojson: write <name>.geojson next to each given trajectory CSV."""
    import sys
    for name in (argv if argv is not None else sys.argv[1:]):
        path = Path(name)
        write_geojson(path.with_suffix(".geojson"), geojson_from_csv(path))
        print(path.with_suffix(".geojson"))
    return 0
