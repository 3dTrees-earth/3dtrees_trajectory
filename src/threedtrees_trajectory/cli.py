"""Estimate the sensor path of moving lidar (ULS, MLS) from its point cloud, and map what the sensor could not
see (occlusion) along that path or along a delivered one."""
from __future__ import annotations

import argparse
import math
import os
import sys
from contextlib import ExitStack
from pathlib import Path
from struct import error as StructError

from . import __version__
from .coverage import NoTrajectory, UnusableTime
from .resources import RunResources, decoding_threads


SUPPLEMENTARY_PATH_NOTE = "Keep the delivered path as supplementary data on the dataset; it will be added once supported."


def _supplementary(reason: str) -> str:
    return f"{reason} {SUPPLEMENTARY_PATH_NOTE}"


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
    with ExitStack() as output_locks:
        return _run(argv, output_locks)


def _run(argv, output_locks) -> int:
    resources = RunResources()
    parser = argparse.ArgumentParser(prog="3dtrees-trajectory", description=__doc__, allow_abbrev=False)
    parser.add_argument("clouds", nargs="*", type=Path, help="standardized LAS/LAZ files from one dataset, in collection order")
    parser.add_argument("-o", "--output", type=Path, help="3Dtrees trajectory CSV")
    parser.add_argument("--debug", action="store_true", help="write a JSON quality and resource report")
    parser.add_argument("--report", type=Path,
                        help="debug report path (requires --debug, defaults next to the CSV, required without a CSV)")
    parser.add_argument("--geojson", type=Path, help="map layer of the path (default next to the CSV)")
    parser.add_argument("--no-platform-check", action="store_true", help="skip static (TLS) cloud detection when estimating a path")
    parser.add_argument("--interval", default="0.5", help="seconds per position (default 0.5)")
    parser.add_argument("--threads", help="maximum decoding threads and GPS-time window workers (defaults to GALAXY_SLOTS or 10)")
    parser.add_argument("--full", action="store_true", help="keep all pulses (requires more memory)")
    parser.add_argument("--trajectory", type=Path,
                        help="delivered path (3Dtrees trajectory CSV); used for occlusion instead of estimating one")
    parser.add_argument("--occlusion-copc", type=Path, help="occlusion voxels as COPC (with --occlusion-cog)")
    parser.add_argument("--occlusion-cog", type=Path,
                        help="occluded share above ground as a COG, with one height band per voxel layer (default 0.5 m)")
    parser.add_argument("--dtm", type=Path,
                        help="terrain raster in the cloud CRS for occlusion (nearest cell, lowest hits in gaps)")
    parser.add_argument("--voxel", default="0.5", help="occlusion voxel size and COG height-band step in metres (default 0.5)")
    parser.add_argument("--grid-memory-gib", type=float, default=16,
                        help="soft occlusion worker budget in GiB before further GPS-time splitting (default 16, maximum 64)")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args, unknown = parser.parse_known_intermixed_args(argv)
    if args.report is not None and not args.debug:
        parser.error("--report requires --debug")
    if args.output is None and args.trajectory is None and args.report is None:
        parser.error("the following argument is required: -o/--output")
    if args.debug and args.output is None and args.report is None:
        parser.error("--debug without -o/--output requires --report PATH")

    occlusion_paths = tuple(p for p in (args.occlusion_copc, args.occlusion_cog) if p is not None)
    report_path = (args.report or args.output.with_suffix(".report.json")) if args.debug else None
    geojson_path = args.geojson or (args.output.with_suffix(".geojson") if args.output is not None else None)
    products = ((args.output, geojson_path) if args.output is not None else ()) + occlusion_paths
    threads = None
    error_message = None
    exit_code = 2
    safe_products = False
    safe_report = False
    try:
        clouds = [path.resolve() for path in args.clouds]
        inputs = clouds + [path.resolve() for path in (args.trajectory, args.dtm) if path is not None]
        resolved_products = tuple(path.resolve() for path in products)
        destinations = resolved_products
        if report_path is not None:
            resolved_report = report_path.resolve()
            _validate_paths((resolved_report,), inputs)
            safe_report = not any(_paths_overlap(resolved_report, product) for product in resolved_products)
            destinations += (resolved_report,)
        _validate_paths(destinations, inputs)
        safe_products = bool(products)
        if unknown:
            raise ValueError(f"Unrecognized arguments: {' '.join(unknown)}")
        if args.output is None and args.trajectory is None:
            raise ValueError("Specify the trajectory CSV with -o/--output.")
        if len(occlusion_paths) == 1:
            raise ValueError("Occlusion writes both products: give --occlusion-copc and --occlusion-cog.")
        if args.trajectory is not None and not occlusion_paths:
            raise ValueError("A delivered --trajectory is only used for occlusion: add --occlusion-copc and --occlusion-cog.")
        if args.dtm is not None and not occlusion_paths:
            raise ValueError("--dtm requires --occlusion-copc and --occlusion-cog.")
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

    from laspy.errors import LaspyException
    from lazrs import LazrsError
    from pyproj.exceptions import ProjError

    from .estimate import STATIC_REFUSAL, detect_platform, estimate
    from .occlusion import TerrainError, open_dtm, prepare_index
    from .output import estimate_quality, geojson, refusal, report, write_csv, write_empty, write_geojson, write_report
    from .pulses import CloudInfo, read_pulses
    from .window_jobs import JobBusy, lock_output

    info = CloudInfo.from_paths(args.clouds)
    pulses = None
    result = None
    data = None
    indexed_clouds = None
    layer = None
    path_inputs = None
    delivered_check = None
    path_warnings = []
    if error_message is None:
        try:
            with resources.stage("read/time-sort", children=True) as stage:
                work = (report_path or args.output or args.occlusion_copc).resolve().parent
                if occlusion_paths:
                    output_locks.enter_context(lock_output(args.occlusion_copc.resolve()))

                    try:
                        indexed_clouds = prepare_index(args.clouds, voxel)
                    except (OSError, ValueError, LaspyException, StructError, EOFError) as error:
                        # Let the trajectory reader retain its own refusal codes and messages.
                        indexed_clouds = {"error": str(error)}
                    if args.dtm is not None and isinstance(indexed_clouds, tuple) and not indexed_clouds[0].copc:
                        with open_dtm(args.dtm, indexed_clouds[1]):
                            pass
                pulses, _ = read_pulses(args.clouds, full=args.full and args.output is not None, info=info,
                                        scratch=work, resources=resources,
                                        scan_index=indexed_clouds[0] if isinstance(indexed_clouds, tuple) else None,
                                        threads=threads)
                stage.update(bytes_read=sum(path.stat().st_size for path in args.clouds),
                             parallel_record_workers=info.read_workers)
            if not len(pulses):
                raise NoTrajectory("No usable pulses could be reconstructed from the input.")
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
                        if any(mode["mode"] == "single" and not mode["has_scan_angle"]
                               for mode in info.return_modes):
                            raise NoTrajectory("Single returns without scan angle: no usable pulse lines or scan "
                                               "angles in one or more files. Please provide the delivered trajectory file.")
                        angle_files = [mode["index"] for mode in info.return_modes if mode["mode"] == "single"]
                        result = estimate(pulses, check_platform=False, interval=interval, angle_files=angle_files, workers=threads)
                        if estimate_quality(result) == "unreliable":
                            if not pulses.has_scan_angle:
                                raise NoTrajectory("Pulse lines miss the fitted path by more than 0.5 m at the "
                                                   "median and there is no usable scan angle.")
                            result = estimate(pulses, check_platform=False, interval=interval, use_lines=False, workers=threads)
                            info.single_return_files = list(range(len(args.clouds)))
                            result.warnings.append("Pulse lines missed the fitted path by more than 0.5 m at the median.")
                        result.evidence = evidence
                    except ValueError as error:
                        if args.trajectory is None:
                            raise
                        result = None
                        path_warnings.append(f"Estimated path: {error}")
                        data["estimated_path_error"] = str(error)
                if result is not None:
                    try:
                        layer = geojson(info, result)
                    except ProjError:
                        result.warnings.append("The map layer could not be projected to WGS84.")
                    path_warnings = result.warnings
                    if safe_report:
                        data = report(info, result, pulses.multi_share)
            if occlusion_paths:
                path_inputs = _occlusion_path(args, result, pulses, resources)
                if not isinstance(path_inputs, dict):
                    delivered_check = path_inputs[3]
            if args.output is not None:
                if result is None:
                    write_empty(args.output)
                else:
                    write_csv(args.output, info, result, delivered_check=delivered_check)
                if layer is None or not layer["features"]:
                    write_empty(geojson_path)
                else:
                    write_geojson(geojson_path, layer)
            if safe_report:
                data["pulse_check"] = {
                    "kept": len(pulses), "multi_return_share": round(pulses.multi_share, 6),
                    "has_scan_angle": pulses.has_scan_angle, "complete_multi_pulses": info.multi_pulses,
                    "incomplete_multi_pulses": info.incomplete_multi_pulses,
                    "out_of_order_records_sampled": info.joined_out_of_order_points,
                    "unpaired_multi_echoes": info.unpaired_multi_echoes,
                    "population_counts_estimated": info.population_counts_estimated,
                    "return_modes": info.return_modes,
                }
        except (NoTrajectory, UnusableTime) as error:
            error_message = _supplementary(str(error))
            exit_code = error.exit_code
        except JobBusy as error:
            safe_products = safe_report = False
            error_message = str(error)
        except TerrainError as error:
            error_message = str(error)
            data = {"occlusion": {"error": error_message}}
        except (OSError, ValueError, MemoryError, LaspyException, LazrsError, StructError, EOFError) as error:
            error_message = str(error) or "The job could not read its input or ran out of memory."
    if error_message is None and occlusion_paths:
        pulses = None
        path_quality = estimate_quality(result) if result is not None and args.trajectory is None else None
        try:
            data["occlusion"] = _occlusion(args, path_inputs, resources, info, requested, indexed_clouds,
                                           path_quality, path_warnings)
        except TerrainError as error:
            error_message = str(error)
            data["occlusion"] = {"error": error_message}
        if error_message is None and "error" in data["occlusion"]:
            if result is None:
                error_message = data["occlusion"]["error"]
            else:
                result.warnings.append(f"Occlusion: {data['occlusion']['error']}")
                try:
                    write_csv(args.output, info, result, delivered_check=delivered_check)
                except OSError as error:
                    error_message = f"cannot publish trajectory metadata {args.output}: {error}"
                for path in occlusion_paths:
                    try:
                        write_empty(path)
                    except OSError as error:
                        error_message = f"cannot publish empty output {path}: {error}"
    if error_message is not None:
        print(f"error: {error_message}", file=sys.stderr)
        if safe_report:
            details = {field: data[field] for field in ("occlusion", "pulse_check") if data is not None and field in data}
            data = refusal(info, pulses, error_message) | details
        if safe_products:
            for path in products:
                try:
                    write_empty(path)
                except OSError as error:
                    print(f"error: cannot publish empty output {path}: {error}", file=sys.stderr)
    if safe_report:
        data["exit_code"] = exit_code if error_message is not None else 0
        data["resources"] = resources.report(threads, info.points_read)
        data["stages"] = resources.stages
        try:
            write_report(report_path, data)
        except OSError as error:
            print(f"error: cannot write report {report_path}: {error}", file=sys.stderr)
            return 2
    if error_message is not None:
        return exit_code
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
                return {"error": _supplementary(str(error) or "The delivered-path check ran out of memory.")}
    time = np.array([f.time for f in result.fits], float)
    xyz = np.array([f.xyz for f in result.fits], float).reshape(-1, 3)
    return time, xyz, result.method, None


def _occlusion(args, path_inputs, resources, info, threads, indexed_clouds, path_quality, path_warnings) -> dict:
    """Run the occlusion part; a failure is reported in its section instead of failing the trajectory."""
    from laspy.errors import LaspyException
    from lazrs import LazrsError

    from . import occlusion

    if isinstance(path_inputs, dict):
        return path_inputs
    time, xyz, method, check = path_inputs
    if isinstance(indexed_clouds, dict):
        return indexed_clouds
    try:
        return occlusion.run(
            args.clouds, time, xyz, args.occlusion_copc, args.occlusion_cog, voxel=float(args.voxel),
            dtm=args.dtm, scratch=args.occlusion_copc.resolve().parent, source="delivered" if check else "estimated",
            quality="rough" if method in ("scan_angle", "lines_scan_angle") or info.single_return_files else None,
            path_method=method, path_quality=path_quality, path_warnings=path_warnings,
            delivered_check=check, threads=threads, prepared=indexed_clouds, single_files=info.single_return_files,
            budget=int(args.grid_memory_gib * 2**30), resources=resources,
        )
    except occlusion.TerrainError:
        raise
    except (OSError, ValueError, MemoryError, ImportError, LaspyException, LazrsError) as error:
        reason = str(error) or "The occlusion step ran out of memory."
        if check is not None and isinstance(error, occlusion.PathCoverageError):
            reason = _supplementary(reason)
        return {**({"delivered_path_check": check} if check else {}), "error": reason}


if __name__ == "__main__":
    sys.exit(main())
