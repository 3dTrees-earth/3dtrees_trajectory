"""Estimate the sensor path of moving lidar (ULS, MLS, ALS) from its point cloud, and map what the sensor could not
see (occlusion) along that path or along a delivered one."""
from __future__ import annotations

import argparse
import math
import os
import sys
from contextlib import ExitStack
from pathlib import Path

from . import __version__
from .resources import RunResources, decoding_threads


def _validate_paths(outputs: tuple[Path, ...], inputs: list[Path]) -> None:
    """Keep output replacement away from inputs and other destinations."""
    for index, path in enumerate(outputs):
        for other in outputs[:index]:
            if _paths_overlap(path, other):
                raise ValueError("All products and the report need separate output paths.")
    if any(_paths_overlap(output, item) for output in outputs for item in inputs):
        raise ValueError("Output paths must differ from input paths.")
    # Publishing replaces the directory entry, so a hard link at an output
    # destination leaves the existing input inode untouched.


def _paths_overlap(first: Path, second: Path) -> bool:
    return first == second or first in second.parents or second in first.parents


def main(argv: list[str] | None = None) -> int:
    resources = RunResources()
    parser = argparse.ArgumentParser(prog="3dtrees-trajectory", description=__doc__, allow_abbrev=False)
    parser.add_argument("clouds", nargs="*", type=Path, help="per-file LAS/LAZ or one merged COPC from the same acquisition")
    parser.add_argument("-o", "--output", type=Path, help="3Dtrees trajectory CSV")
    parser.add_argument("--report", type=Path, help="JSON quality report (default next to the CSV)")
    parser.add_argument("--geojson", type=Path, help="map layer of the path (default next to the CSV)")
    parser.add_argument("--no-platform-check", action="store_true", help="skip static (TLS) cloud detection when estimating a path")
    parser.add_argument("--interval", default="0.5", help="seconds per position (default 0.5)")
    parser.add_argument("--threads", help="maximum LAZ decoding threads (defaults to GALAXY_SLOTS or 8)")
    parser.add_argument("--full", action="store_true", help="keep all pulses (requires more memory)")
    parser.add_argument("--trajectory", type=Path,
                        help="delivered path (3Dtrees trajectory CSV); used for occlusion instead of estimating one")
    parser.add_argument("--occlusion-copc", type=Path, help="occlusion voxels as COPC (with --occlusion-cog)")
    parser.add_argument("--occlusion-cog", type=Path, help="occluded share per metre above ground as a multiband COG")
    parser.add_argument("--voxel", default="0.5", help="occlusion voxel size in metres (default 0.5)")
    parser.add_argument("--grid-memory-gib", type=float, default=16,
                        help="soft occlusion grid budget in GiB before tiling (default 16, maximum 64)")
    parser.add_argument("--dtm", type=Path, help="terrain raster for heights above ground (default: lowest hits)")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args, unknown = parser.parse_known_intermixed_args(argv)
    if args.output is None and args.report is None:
        parser.error("the following argument is required: -o/--output")

    occlusion_paths = tuple(p for p in (args.occlusion_copc, args.occlusion_cog) if p is not None)
    report_path = args.report or args.output.with_suffix(".report.json")
    geojson_path = args.geojson or (args.output.with_suffix(".geojson") if args.output is not None else None)
    products = ((args.output, geojson_path) if args.output is not None else ()) + occlusion_paths
    threads = None
    error_message = None
    safe_products = False
    safe_report = False
    try:
        clouds = [path.resolve() for path in args.clouds]
        inputs = clouds + [item.resolve() for item in (args.trajectory, args.dtm) if item is not None]
        resolved_report = report_path.resolve()
        resolved_products = tuple(path.resolve() for path in products)
        _validate_paths((resolved_report,), inputs)
        safe_report = not any(_paths_overlap(resolved_report, product) for product in resolved_products)
        _validate_paths((*resolved_products, resolved_report), inputs)
        safe_products = bool(products)
        if unknown:
            raise ValueError(f"Unrecognized arguments: {' '.join(unknown)}")
        if args.output is None and args.trajectory is None:
            raise ValueError("Specify the trajectory CSV with -o/--output.")
        if len(occlusion_paths) == 1:
            raise ValueError("Occlusion writes both products: give --occlusion-copc and --occlusion-cog.")
        if args.trajectory is not None and not occlusion_paths:
            raise ValueError("A delivered --trajectory is only used for occlusion: add --occlusion-copc and --occlusion-cog.")
        try:
            voxel = float(args.voxel)
        except ValueError:
            voxel = math.nan
        if not math.isfinite(voxel) or not 0.05 <= voxel <= 10:
            raise ValueError("--voxel must be between 0.05 and 10 m.")
        if not math.isfinite(args.grid_memory_gib) or not 0 < args.grid_memory_gib <= 64:
            raise ValueError("--grid-memory-gib must be positive and at most 64.")
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

    from .cloud import sorted_copc
    from .estimate import STATIC_REFUSAL, detect_platform, estimate
    from .output import geojson, refusal, report, write_csv, write_empty, write_geojson, write_report
    from .pulses import CloudInfo, read_pulses

    info = CloudInfo.from_paths(args.clouds)
    pulses = None
    result = None
    data = None
    with ExitStack() as scratch:
        if error_message is None:
            try:
                prepared = scratch.enter_context(sorted_copc(args.clouds, Path.cwd(), resources=resources))
                read_stage = "read/time-sort merge" if prepared is not None else "read/time-sort"
                with resources.stage(read_stage) as stage:
                    pulses, _ = read_pulses(args.clouds, full=args.full and args.output is not None, info=info,
                                            prepared=prepared)
                    stage["bytes_read"] = (prepared.bytes_written if prepared is not None
                                           else sum(path.stat().st_size for path in args.clouds))
                if not len(pulses):
                    raise ValueError("No usable pulses could be reconstructed from the input.")
                evidence = "given"
                with resources.stage("platform check"):
                    if args.trajectory is not None or not args.no_platform_check:
                        platform, evidence = detect_platform(pulses)
                        if platform == "static":
                            raise ValueError(STATIC_REFUSAL.format(evidence=evidence))
                data = {"tool": f"3dtrees_trajectory {__version__}", "files": info.files, "trajectory": "delivered"}
                if args.output is not None:
                    with resources.stage("path estimation"):
                        try:
                            result = estimate(pulses, check_platform=False, interval=interval)
                            result.evidence = evidence
                            data_estimate = report(info, result, pulses.multi_share)
                            if data_estimate["quality"] == "unreliable":
                                result = None
                                raise ValueError("Pulse lines miss the fitted path by more than 0.5 m at the median. "
                                                 "Check return pairing or provide the delivered trajectory.")
                        except ValueError as error:
                            if args.trajectory is None:
                                raise
                            result = None
                            data["estimated_path_error"] = str(error)
                    if result is None:
                        write_empty(args.output)
                        write_empty(geojson_path)
                    else:
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
                data["pulse_check"] = {
                    "kept": len(pulses), "multi_return_share": round(pulses.multi_share, 6),
                    "has_scan_angle": pulses.has_scan_angle, "complete_multi_pulses": info.multi_pulses,
                    "incomplete_multi_pulses": info.incomplete_multi_pulses,
                    "out_of_order_records_sampled": info.joined_out_of_order_points,
                }
            except (OSError, ValueError, MemoryError) as error:
                error_message = str(error) or "The job ran out of memory."
        if error_message is None and occlusion_paths:
            path_inputs = _occlusion_path(args, result, pulses, resources)
            pulses = None
            data["occlusion"] = _occlusion(args, path_inputs, prepared, resources)
            if "error" in data["occlusion"]:
                if result is None:
                    error_message = data["occlusion"]["error"]
                else:
                    result.warnings.append(f"Occlusion: {data['occlusion']['error']}")
                    for path in occlusion_paths:
                        try:
                            write_empty(path)
                        except OSError as error:
                            error_message = f"cannot publish empty output {path}: {error}"
    if error_message is not None:
        print(f"error: {error_message}", file=sys.stderr)
        details = {field: data[field] for field in ("occlusion", "pulse_check") if data is not None and field in data}
        data = refusal(info, pulses, error_message) | details
        if safe_products:
            for path in products:
                try:
                    write_empty(path)
                except OSError as error:
                    print(f"error: cannot publish empty output {path}: {error}", file=sys.stderr)
    if safe_report:
        data["resources"] = resources.report(threads, info.points_read)
        data["stages"] = resources.stages
        try:
            write_report(report_path, data)
        except OSError as error:
            print(f"error: cannot write report {report_path}: {error}", file=sys.stderr)
            return 2
    if error_message is not None:
        return 2
    if result is not None:
        print(f"path ({result.evidence}) via {result.method}: {len(result.fits)} positions -> {args.output}")
        for warning in result.warnings:
            print(f"warning: {warning}", file=sys.stderr)
    if occlusion_paths and "error" not in data["occlusion"]:
        print(f"occlusion: {data['occlusion']['occluded_share']} of the observed voxels hidden -> "
              f"{args.occlusion_copc}, {args.occlusion_cog}")
    return 0


def _occlusion_path(args, result, pulses, resources):
    import numpy as np

    from . import occlusion

    if args.trajectory is not None:
        with resources.stage("delivered-path check") as stats:
            try:
                time, xyz = occlusion.read_trajectory(args.trajectory)
                occlusion.validate_trajectory(time, xyz)
                check = occlusion.delivered_miss(pulses, time, xyz)
                stats["bytes_read"] = args.trajectory.stat().st_size
                return time, xyz, "delivered", check
            except (OSError, ValueError, MemoryError) as error:
                return {"error": str(error) or "The delivered-path check ran out of memory."}
    time = np.array([f.time for f in result.fits], float)
    xyz = np.array([f.xyz for f in result.fits], float).reshape(-1, 3)
    return time, xyz, result.method, None


def _occlusion(args, path_inputs, prepared, resources) -> dict:
    """Run the occlusion part; a failure is reported in its section instead of failing the trajectory."""
    from laspy.errors import LaspyException
    from lazrs import LazrsError

    from . import occlusion

    if isinstance(path_inputs, dict):
        return path_inputs
    time, xyz, method, check = path_inputs
    try:
        return occlusion.run(
            args.clouds, time, xyz, args.occlusion_copc, args.occlusion_cog, voxel=float(args.voxel),
            dtm=args.dtm, scratch=args.occlusion_copc.resolve().parent, source="delivered" if check else "estimated",
            quality="rough" if method == "scan_angle" else None, delivered_check=check, prepared=prepared,
            budget=int(args.grid_memory_gib * 2**30), resources=resources,
        )
    except (OSError, ValueError, MemoryError, ImportError, LaspyException, LazrsError) as error:
        return {**({"delivered_path_check": check} if check else {}),
                "error": str(error) or "The occlusion step ran out of memory."}


if __name__ == "__main__":
    sys.exit(main())
