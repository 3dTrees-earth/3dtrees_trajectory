"""Estimate the sensor path of moving lidar (ULS, MLS, ALS) from its point cloud."""
from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

from . import __version__
from .resources import RunResources, decoding_threads


def _validate_paths(outputs: tuple[Path, ...], clouds: list[Path]) -> None:
    """Keep output replacement away from inputs and other destinations."""
    for index, path in enumerate(outputs):
        for other in outputs[:index]:
            if _paths_overlap(path, other):
                raise ValueError("CSV, GeoJSON and report need separate output paths.")
    if not set(outputs).isdisjoint(clouds):
        raise ValueError("Output paths must differ from input paths.")
    # Publishing replaces the directory entry, so a hard link at an output
    # destination leaves the existing input inode untouched.


def _paths_overlap(first: Path, second: Path) -> bool:
    return first == second or first in second.parents or second in first.parents


def main(argv: list[str] | None = None) -> int:
    resources = RunResources()
    parser = argparse.ArgumentParser(prog="3dtrees-trajectory", description=__doc__, allow_abbrev=False)
    parser.add_argument("clouds", nargs="*", type=Path, help="LAS/LAZ files of one flight or walk")
    parser.add_argument("-o", "--output", type=Path, help="3Dtrees trajectory CSV")
    parser.add_argument("--report", type=Path, help="JSON quality report (default next to the CSV)")
    parser.add_argument("--geojson", type=Path, help="map layer of the path (default next to the CSV)")
    parser.add_argument("--no-platform-check", action="store_true", help="skip static (TLS) cloud detection")
    parser.add_argument("--interval", default="0.5", help="seconds per position (default 0.5)")
    parser.add_argument("--threads", help="maximum LAZ decoding threads (defaults to GALAXY_SLOTS or 8)")
    parser.add_argument("--full", action="store_true", help="keep all pulses (requires more memory)")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args, unknown = parser.parse_known_intermixed_args(argv)
    if args.output is None and args.report is None:
        parser.error("the following argument is required: -o/--output")

    report_path = args.report or args.output.with_suffix(".report.json")
    geojson_path = args.geojson or (args.output.with_suffix(".geojson") if args.output is not None else None)
    products = (args.output, geojson_path) if args.output is not None else ()
    threads = None
    error_message = None
    safe_products = False
    safe_report = False
    try:
        clouds = [path.resolve() for path in args.clouds]
        resolved_report = report_path.resolve()
        resolved_products = tuple(path.resolve() for path in products)
        _validate_paths((resolved_report,), clouds)
        safe_report = not any(_paths_overlap(resolved_report, product) for product in resolved_products)
        _validate_paths((*resolved_products, resolved_report), clouds)
        safe_products = bool(products)
        if unknown:
            raise ValueError(f"Unrecognized arguments: {' '.join(unknown)}")
        if args.output is None:
            raise ValueError("Specify the trajectory CSV with -o/--output.")
        if not args.clouds:
            raise ValueError("Provide at least one LAS/LAZ file.")
        try:
            interval = float(args.interval)
        except ValueError:
            interval = math.nan
        if not math.isfinite(interval) or interval <= 0:
            raise ValueError("--interval must be a positive finite number of seconds.")
        try:
            requested = int(args.threads) if args.threads is not None else None
        except ValueError:
            requested = 0
        if requested is not None and requested < 1:
            raise ValueError("--threads must be a positive integer.")
        threads = decoding_threads(requested)
    except (OSError, RuntimeError, ValueError) as error:
        error_message = str(error)

    # Initialise native pools before NumPy, SciPy and lazrs are imported. The
    # solver runs after decoding and uses one BLAS thread for stable reductions.
    os.environ["PROJ_NETWORK"] = "OFF"
    os.environ["RAYON_NUM_THREADS"] = str(threads or 1)
    for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS"):
        os.environ[variable] = "1"

    from pyproj.exceptions import ProjError

    from .estimate import estimate
    from .output import geojson, refusal, report, write_csv, write_empty, write_geojson, write_report
    from .pulses import CloudInfo, read_pulses

    info = CloudInfo.from_paths(args.clouds)
    pulses = None
    if error_message is None:
        try:
            pulses, _ = read_pulses(args.clouds, full=args.full, info=info)
            result = estimate(pulses, check_platform=not args.no_platform_check, interval=interval)
            try:
                layer = geojson(info, result)
            except ProjError:
                layer = None
                result.warnings.append("The map layer could not be projected to WGS84.")
            write_csv(args.output, info, result)
            if layer is None or not layer["features"]:
                write_empty(geojson_path)
            else:
                write_geojson(geojson_path, layer)
            data = report(info, result, pulses.multi_share)
        except (OSError, ValueError, MemoryError) as error:
            error_message = str(error) or "The job ran out of memory."
    if error_message is not None:
        print(f"error: {error_message}", file=sys.stderr)
        data = refusal(info, pulses, error_message)
        if safe_products:
            for path in products:
                try:
                    write_empty(path)
                except OSError as error:
                    print(f"error: cannot publish empty output {path}: {error}", file=sys.stderr)
    if safe_report:
        data["resources"] = resources.report(threads, info.points_read)
        try:
            write_report(report_path, data)
        except OSError as error:
            print(f"error: cannot write report {report_path}: {error}", file=sys.stderr)
            return 2
    if error_message is not None:
        return 2
    print(f"path ({result.evidence}) via {result.method}: {len(result.fits)} positions -> {args.output}")
    for warning in result.warnings:
        print(f"warning: {warning}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
