"""3dtrees-trajectory: the sensor path of a moving lidar (ULS, MLS, ALS) estimated from the point cloud alone."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from . import __version__


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="3dtrees-trajectory", description=__doc__)
    parser.add_argument("clouds", nargs="+", type=Path, help="LAS/LAZ files of one flight or walk (several files are read together)")
    parser.add_argument("-o", "--output", type=Path, required=True, help="3Dtrees trajectory CSV")
    parser.add_argument("--report", type=Path, help="JSON quality report (default: next to the CSV)")
    parser.add_argument("--geojson", type=Path, help="map layer of the path (default: next to the CSV)")
    parser.add_argument("--no-platform-check", action="store_true",
                        help="skip the check that refuses static (TLS) clouds")
    parser.add_argument("--interval", type=float, default=0.5, help="seconds per position (default 0.5)")
    parser.add_argument("--threads", type=int, default=8, help="threads for LAZ decoding (default 8)")
    parser.add_argument("--full", action="store_true", help="read and keep every point (slower, more memory)")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args(argv)
    os.environ["RAYON_NUM_THREADS"] = str(args.threads)  # read by the LAZ decoder when its thread pool starts

    from .estimate import estimate
    from .output import geojson, refusal, report, write_csv, write_geojson, write_report
    from .pulses import read_pulses

    pulses, info = read_pulses(args.clouds, full=args.full)
    try:
        result = estimate(pulses, check_platform=not args.no_platform_check, interval=args.interval)
    except ValueError as error:
        print(f"error: {error}", file=sys.stderr)
        write_report(args.report or args.output.with_suffix(".report.json"), refusal(info, pulses, str(error)))
        return 2
    write_csv(args.output, info, result)
    write_geojson(args.geojson or args.output.with_suffix(".geojson"), geojson(info, result))
    data = report(info, result, pulses.multi_share)
    write_report(args.report or args.output.with_suffix(".report.json"), data)
    print(f"path ({result.evidence}) via {result.method}: {len(result.fits)} positions -> {args.output}")
    for warning in result.warnings:
        print(f"warning: {warning}", file=sys.stderr)
    return 0 if result.fits else 1


if __name__ == "__main__":
    sys.exit(main())
