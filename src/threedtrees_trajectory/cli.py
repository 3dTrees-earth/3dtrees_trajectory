"""3dtrees-trajectory: sensor path (ULS, MLS, ALS) or scanner positions (TLS) estimated from a point cloud alone."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from . import __version__


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="3dtrees-trajectory", description=__doc__)
    parser.add_argument("clouds", nargs="+", type=Path, help="LAS/LAZ files; for TLS one file per scan works best")
    parser.add_argument("-o", "--output", type=Path, required=True, help="3Dtrees trajectory CSV")
    parser.add_argument("--report", type=Path, help="JSON quality report (default: next to the CSV)")
    parser.add_argument("--platform", choices=["auto", "mobile", "static"], default="auto",
                        help="mobile = ULS/MLS/ALS path, static = TLS scan positions (default: detected)")
    parser.add_argument("--split", choices=["auto", "file", "source", "time"], default="auto",
                        help="static only: how scans are told apart (default: point sources, then files, then time jumps)")
    parser.add_argument("--interval", type=float, default=0.5, help="mobile only: seconds per position (default 0.5)")
    parser.add_argument("--threads", type=int, default=8, help="threads for LAZ decoding (default 8)")
    parser.add_argument("--full", action="store_true", help="read and keep every point (slower, more memory)")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args(argv)
    os.environ["RAYON_NUM_THREADS"] = str(args.threads)  # read by the LAZ decoder when its thread pool starts

    from .estimate import estimate
    from .output import refusal, report, write_csv, write_report
    from .pulses import read_pulses

    pulses, info = read_pulses(args.clouds, full=args.full)
    try:
        result = estimate(pulses, platform=args.platform, split=args.split, interval=args.interval)
    except ValueError as error:
        print(f"error: {error}", file=sys.stderr)
        write_report(args.report or args.output.with_suffix(".report.json"), refusal(info, pulses, str(error)))
        return 2
    write_csv(args.output, info, result)
    data = report(info, result, pulses.multi_share)
    write_report(args.report or args.output.with_suffix(".report.json"), data)
    print(f"{result.platform} ({result.evidence}) via {result.method}: {len(result.fits)} positions -> {args.output}")
    for warning in result.warnings:
        print(f"warning: {warning}", file=sys.stderr)
    return 0 if result.fits else 1


if __name__ == "__main__":
    sys.exit(main())
