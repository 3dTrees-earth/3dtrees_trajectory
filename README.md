# 3dtrees_trajectory

Estimate the sensor path of a moving lidar acquisition (ULS, MLS) from its standardized point cloud, and map what
the sensor could not see (occlusion) along that path or along one or more delivered paths. Version 1.0.0, MIT.

- **Path estimation** from multi-return pulse lines (`lines`, cm-accurate on checked datasets) or from scan angles
  of single returns (`scan_angle`, about 1 m). Output: a 3Dtrees trajectory CSV (`3dtrees-trajectory/estimated`).
- **Occlusion** by tracing every pulse from the sensor through a voxel grid (OccPy 0.3.0 traversal): a voxel COPC
  with classes filled / empty / occluded / unobserved, and a COG with the occluded share per height band above
  ground.
- **Delivered paths**: one or more vendor trajectories (3Dtrees trajectory CSVs) can drive the occlusion instead of
  the estimate. Each part of the cloud uses the path that covers its time and fits its pulse lines best; paths in
  another EPSG CRS or on another GPS clock are converted.

## Quick start

```bash
# path only
3dtrees-trajectory standardized/strip1.laz standardized/strip2.laz -o trajectory.csv

# path + occlusion, with the Raster DTM
3dtrees-trajectory standardized/*.laz -o trajectory.csv --dtm raster/dataset_dtm.tif \
  --occlusion-copc occlusion.copc.laz --occlusion-cog occlusion.tif

# occlusion along delivered paths (one or more files)
3dtrees-trajectory standardized/*.laz --trajectory flight1.trajectory.csv flight2.trajectory.csv \
  --occlusion-copc occlusion.copc.laz --occlusion-cog occlusion.tif
```

In Docker: `docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/data" -w /data
ghcr.io/3dtrees-earth/3dtrees_trajectory:1.0.0 python /src/run.py …`.

## Inputs in one paragraph

All standardized LAS/LAZ files of one acquisition, in collection order (1 to 65,536 files, LAS 1.2–1.4, point
formats with GPS time). GPS time must be finite and vary; files share one clock and metric coordinates. Keep the
scanner's record order, ReturnNumber, NumberOfReturns, ScanAngle, PointSourceID and scanner channel. A cloud without
a usable CRS (none, or Standardization's EPSG:4978 placeholder near the origin) is **local**: outputs say
`crs: local` and stay in cloud coordinates. Optional: one terrain raster (`--dtm`, may have no CRS) and one or more
delivered paths (`--trajectory`). Details: [docs/inputs.md](docs/inputs.md).

**Options** (also plain fields in the Galaxy form): `--voxel` voxel edge and COG band step (default 0.5 m; Galaxy 0.25–2 m),
`--interval` seconds between path positions (default 0.5). The path is fitted from **all pulses** by default, streamed
in GPS-time windows with bounded memory; `--sample` fits it from a bounded sample instead (`--full` is still accepted).

## Outputs

| File | Contents |
| --- | --- |
| `trajectory.csv` | `# format: 3dtrees-trajectory/estimated`, CRS, clock, method, quality, warnings, delivered-path coverage; columns `time,x,y,z,roll,pitch,yaw,name` (attitude and name empty) |
| `occlusion.copc.laz` | one point per filled, occluded or unobserved voxel: `occlusion` class, `nhit`, `nmiss`, `nocc`, `occluded_fraction`, `height_above_ground` |
| `occlusion.tif` | COG, one band per height layer (`occluded_0-0.5m`, …) plus `unobserved_column`; provenance tags |
| `*.report.json` | with `--debug`: quality, refusal reason, pairing, per-segment path choice, stage resources |

Details: [docs/outputs.md](docs/outputs.md).

## Exit codes

| Code | Meaning |
| --- | --- |
| `0` | Path estimated and/or occlusion written |
| `2` | Invalid input or arguments, static (TLS) acquisition, Earth-centred coordinates, requested occlusion that cannot run (e.g. COPC input, no delivered-path coverage), failed output |
| `3` | GPS time missing, non-finite or constant |
| `4` | No usable path (no pulse lines, no scan angle, no delivered path) |
| `127` | Out of memory or a lost worker; rerun with more memory |

## Documentation

| Page | Contents |
| --- | --- |
| [docs/inputs.md](docs/inputs.md) | point clouds, CRS and clock rules, DTM, delivered path format |
| [docs/path-estimation.md](docs/path-estimation.md) | `lines`, `scan_angle`, pulse pairing (incl. returns stored apart), quality, static check |
| [docs/delivered-paths.md](docs/delivered-paths.md) | CRS conversion, clock fallback, several paths per dataset, coverage and refusal |
| [docs/occlusion.md](docs/occlusion.md) | traversal, voxel classes, ground, products, GPS-time windows, checkpoints |
| [docs/outputs.md](docs/outputs.md) | every header key, COPC dimension, COG band and tag, debug report |
| [docs/cli.md](docs/cli.md) | all options with ranges, Galaxy wrapper options |
| [docs/resources.md](docs/resources.md) | threads, memory, voxel-size scaling, scratch, Galaxy TPV rule, benchmarks |
| [docs/validation.md](docs/validation.md) | accuracy against delivered paths, test matrix, product equality |
| [docs/platform.md](docs/platform.md) | 3Dtrees platform gate, workflows, Galaxy, CI and releases |
| [docs/limitations.md](docs/limitations.md) | what 1.0.0 does not do, and the 1.1 plans |
| [docs/REPORT.md](docs/REPORT.md) | development report: stages, streaming evaluation, round-by-round measurements |

## Development

```bash
make test    # Docker suite (2 CPUs, 4 GiB): synthetic LAS/LAZ with known positions
make image   # runtime image 3dtrees_trajectory:dev
```

CI tests every pull request and push; it publishes `ghcr.io/3dtrees-earth/3dtrees_trajectory:<version>` for a
`v*` tag and `:latest` for `main`. See [docs/platform.md](docs/platform.md).

## Credits and licences

This tool is MIT-licensed. Its code was written for 3Dtrees; no third-party source code was copied into it except as
noted below.

- **Occlusion ray walk** (`src/threedtrees_trajectory/ray_walk.cpp`): a C++ port of the voxel traversal and echo
  classification of [OccPy](https://github.com/dkueken/OccPy) 0.3.0 (MIT, © 2023 Daniel Kükenbrink; notice kept in the
  file header). Method: Kükenbrink et al. (2017), *Remote Sensing of Environment*,
  [doi:10.1016/j.rse.2016.10.023](https://doi.org/10.1016/j.rse.2016.10.023), based on Amanatides & Woo (1987).
- **Path from pulse lines** (`solve.py`, method `lines`): the sensor lies where the lines through the first and last
  return of multi-return pulses meet, per time window — Gatziolis & McGaughey (2019), *Remote Sensing*,
  [doi:10.3390/rs11192258](https://doi.org/10.3390/rs11192258); the same idea underlies lidR's `track_sensor()`
  (Roussel2020) and PDAL's `filters.trajectory` (Karney 2022,
  [arXiv:2208.12116](https://arxiv.org/abs/2208.12116)). Our implementation is new Python code; no lidR (GPL-3) or PDAL
  (BSD-3) code is included.
- **Path from scan angles** (`solve.py`, method `scan_angle`): a cone fit of single returns on their recorded scan angle,
  written for this tool; Karney (2022) uses scan angles in a related way.
