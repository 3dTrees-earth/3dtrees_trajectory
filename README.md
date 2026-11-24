# 3dtrees_trajectory

The tool supports ULS and MLS. It estimates the flight or walking path from a standardized point cloud
collection and writes a 3Dtrees trajectory CSV and a map layer. It can also map occlusion along the estimated or
delivered path. Use `--debug` to write a separate quality and resource report. The current package version is 1.0.0.

## Run

```bash
3dtrees-trajectory standardized/strip1.las standardized/strip2.laz -o trajectory/trajectory.csv
```

Pass all standardized LAS/LAZ files from one dataset to a single invocation, in consistent collection order.
The collection can contain one to 65,536 files. The parallel chunked LAZ reader preserves the scanner's
record order within each file, including time blocks. Keep pulse siblings and their metadata together during
standardization.

The reader supports LAS/LAZ 1.2–1.4 and point formats 0–10. Each file needs finite, varying GPS time.
All files must share a GPS time system and metric coordinates. If the CRS is missing, the tool assumes
metres and preserves the input coordinates. The tool refuses point formats 0 and 2 because they lack GPS time.
Like Raster, Trajectory treats a geocentric CRS as a local placeholder when the cloud bounds are centred less
than 1,000 km from the origin. These products use local coordinates with unknown CRS. For Earth-centred
coordinates, reproject the cloud to a projected CRS in metres.

### Outputs

| File | Contents |
| --- | --- |
| `trajectory.csv` | `3dtrees-trajectory/2` with columns `time,x,y,z,roll,pitch,yaw,name`. Time and XYZ use the cloud's clock and coordinates. Attitude and name are empty. The `#` header records method, quality, source, and any warnings or delivered-path miss flag. |
| `trajectory.geojson` | Path segments separated at gaps over 5 seconds, thinned to about one vertex per metre. Uses WGS84 when the CRS is known and cloud coordinates otherwise. |

`--geojson` overrides the default map layer path next to the CSV. The tool creates output directories
and publishes each file by rename. The GeoJSON file is zero bytes when there are no drawable segments or
reprojection fails. `3dtrees-trajectory-geojson trajectory/trajectory.csv` can also draw the layer from an
existing trajectory CSV.

### Debug report

`--debug` writes `trajectory.report.json` next to the CSV. It records method, quality, warnings, flags,
refusal reason, exit code, point counts, per-file return modes, and stage and job resources.
Use `--debug --report PATH` to choose the report location. `--report PATH` requires `--debug` and
is required for an occlusion-only debug run. The Galaxy wrapper should publish JSON only when debug is enabled.

### Exit codes

| Code | Meaning |
| --- | --- |
| `0` | The path was estimated or the delivered-path occlusion products were written. In a combined job, the estimated path may succeed with an occlusion warning on stderr and in its CSV header. |
| `2` | Invalid arguments or data, a detected static (TLS) acquisition, or a failed output operation. |
| `3` | GPS time is missing, nonfinite or constant in a moving acquisition. A constant point clock cannot be aligned with a moving delivered path. |
| `4` | The cloud has no usable estimate and no delivered path to support occlusion. |
| `127` | The tool ran out of memory or lost a processing worker. Rerun with more memory. |

Refusals print their reason to stderr and record it in the debug report when requested. Safe, writable product
destinations receive zero-byte files. Keep a delivered path as supplementary data when the tool cannot use it.
If the kernel kills the main process, the caller still sees signal 9 (shell exit status 137).

Supply an output CSV for an estimated path. The CSV, GeoJSON, report and input paths must be distinct.
If you set `--interval`, use a positive finite number of seconds (default 0.5). An interval must be long enough for
the cloud's GPS times to fit into a signed 64-bit window index. `--threads` takes a positive integer.

## Occlusion

```bash
3dtrees-trajectory standardized/strip1.las standardized/strip2.laz -o trajectory.csv \
  --occlusion-copc occlusion.copc.laz --occlusion-cog occlusion.tif

3dtrees-trajectory standardized/strip1.las standardized/strip2.laz --trajectory delivered.csv \
  --dtm raster/dataset_dtm.tif --occlusion-copc occlusion.copc.laz --occlusion-cog occlusion.tif
```

Supply both occlusion destinations. When `--trajectory` is supplied, occlusion uses that path. With no
`-o` the tool checks the cloud and writes the occlusion products without fitting a path.
A delivered moving-path CSV starts with `time,x,y,z`, optionally preceded by `#` metadata. It needs at least two
positions, finite coordinates, and unique timestamps in the cloud's GPS clock and coordinate system.
Resolve week seconds versus adjusted GPS time before running. Both modes require a moving sensor.

The compiled walker uses the OccPy 0.3.0 grid traversal and linearly interpolates the sensor position for each pulse. Echoes up
to one second outside the path take the nearest endpoint. More distant echoes are counted and skipped.
The tool refuses disjoint time ranges. A voxel is filled (1) if there is a registered return, empty (2)
if seen clear, occluded (3) if only traversed behind returns, or unobserved (4) in the absence of traversal.

Voxel size defaults to 0.5 m. `--voxel` accepts 0.05–10 m. With `--dtm PATH`, the tool uses the terrain
raster cell containing each voxel column centre (nearest-neighbour sampling). The raster tool's 1 m
`<prefix>_dtm.tif` is a supported input. If both inputs declare a real CRS, they must agree. A DTM without a
CRS or with a local placeholder is sampled in the cloud coordinates. When one input has a real CRS and the
other has no usable CRS, the products and debug report record the coordinate assumption. A DTM outside
the cloud triggers a warning and falls back to lowest hits. Nodata and uncovered columns also use lowest hits.
That estimate starts at the centre of the lowest filled voxel per column, fills empty
columns from the nearest estimate and replaces heights more than 2 m above the local 5 m median.
The highest filled voxel caps each column.

| File | Contents |
| --- | --- |
| `occlusion.copc.laz` | One point at each filled, occluded or unobserved voxel centre from ground through canopy top. Extra dimensions are `occlusion`, `nhit`, `nmiss`, `nocc`, `occluded_fraction` and `height_above_ground`. The fraction is `nocc / (nhit + nmiss + nocc)`, with NaN for unobserved voxels. |
| `occlusion.tif` | North-up COG with one pixel per voxel column and one height band per voxel layer. The default bands are `occluded_0-0.5m`, `occluded_0.5-1m` and so on. Each value is the occluded voxel count divided by the sum of filled, empty and occluded voxels, with NaN when the sum is zero. The final `unobserved_column` band gives the unobserved share of all voxels from ground through canopy top. |

With step `s`, zero-based band `k` covers `[k·s, (k+1)·s)` metres above ground. The pixel size and height-band
step both follow `--voxel` and default to 0.5 m. The COG records the step in `OCCLUSION_BAND_M` and the debug report
records it in `occlusion.cog_band_step_m`.

The COG and COPC carry `OCCLUSION_GROUND` with the value `dtm`, `lowest_hit`, or `dtm+lowest_hit`. Mixed
runs include `OCCLUSION_GROUND_DTM_SHARE` and `OCCLUSION_GROUND_LOWEST_HIT_SHARE`, measured over all
columns in the XY grid. The debug report records `occlusion.ground` and, when a DTM is supplied,
`occlusion.ground_column_shares`.

If there are no qualifying voxels, the COPC is zero bytes and the COG is valid with NaN values. On an
occlusion error both destinations are zero bytes and stderr gives the reason. A combined job can
still return its estimated trajectory after a traversal error. An invalid terrain input exits 2.
The debug report includes traversal counts, in-range class totals and the pooled occluded share per height band.

The first read builds the trajectory sample, indexes record ranges by GPS time and collects their column extrema.
The path selects which echoes overlap the occlusion run. Blocks that intersect a path boundary are decoded to
apply the exact clock limits. Separated returns are joined by GPS time, point source and scanner channel within
each file.

Occlusion uses independent workers over disjoint GPS-time windows. Workers seek the original LAS/LAZ and
trace each complete pulse once. They use the same voxel lattice and store compact counts between the ground
and canopy. The planner splits a window further in time when its working memory or echo count exceeds the
budget. Set the planning ceiling with `--grid-memory-gib` (16 GiB by default). At runtime the scheduler reserves
the estimated memory of each window within the available cgroup memory and `--threads` limit. Clean file
cache is reclaimable. Dirty and writeback pages count toward the working budget. The scheduler adjusts
subsequent reservations using each completed window's peak RSS. If a worker runs out of memory, the job exits
with code 127 so the scheduler or operator can retry it with more memory.

Workers merge counts and classify independent voxel blocks. The writer publishes COPC records and raster windows
in grid order. Repeat the same command to resume an interrupted run and reuse finished windows. Keep the
`.occlusion-*` checkpoint beside the COPC until both outputs are published. A checkpoint write failure
from exhausted storage reports exit code 2.
A merged COPC can supply the trajectory estimate. For occlusion, supply the standardized LAS/LAZ collection.

The delivered path is used for every occlusion ray when supplied. `OCCLUSION_DELIVERED_MISS_FLAG` records whether
the median perpendicular distance from that path to usable pulse lines exceeds 0.1 m. Its values are `true`,
`false`, and `unavailable` when there are no usable lines. A CSV produced in the same run records the result
as `# delivered_path_miss_flag:`. OccPy traces single-return shots independently.

Both occlusion products record the path source (`OCCLUSION_TRAJECTORY`), path method
(`OCCLUSION_PATH_METHOD`), quality (`OCCLUSION_QUALITY`), any warnings (`OCCLUSION_WARNINGS`), and the
delivered-path miss flag when applicable. A scan-angle estimate is marked `rough`. The COG stores these as
dataset tags. The COPC stores the same keys in the `3dtrees` extended variable length record (record ID 1).

COPC indexing and COG conversion use temporary files beside the product destinations. The job needs space for
the uncompressed voxel LAS, intermediate TIFF and final outputs. Wide rasters use a band-interleaved
intermediate and up to four job slots for the COG conversion. They skip intermediate compression when scratch
has room for the uncompressed raster. The image includes Untwine 1.5.1. For local
runs, install the package with `python -m pip install .` using a C++17 compiler and set
`TRAJECTORY_UNTWINE` to an executable Untwine 1.5.1 installation.

Untwine randomizes COPC overview sampling. For regression checks, compare decoded records in voxel-centre order.
[Untwine 1.5.1 sampling implementation](https://github.com/hobuinc/untwine/blob/1.5.1/bu/Processor.cpp#L99-L126)

## Platform pre-gate

Use the standardization metadata to check that GPS time is present and varying. Treat
`max(NumberOfReturns) > 1` as evidence of multi-returns. A single-return cloud needs usable scan angles or an
available delivered path for occlusion. Check for a scan-angle dimension and actual nonzero values, since a LAS
point format can contain an angle field filled with zero. Multi-return clouds can estimate a path without angles
when they have enough usable pulse lines. Preserve GPS time, ReturnNumber, NumberOfReturns, ScanAngle
(or ScanAngleRank for legacy formats), PointSourceID and scanner channel wherever the input format provides
them. The tool checks usable pulses, the multi-return share, return pairing,
and static-scan evidence and prints its refusal to stderr. The platform gate can route clearly unsuitable
datasets before scheduling a job.

## How it works

| Input evidence | Method |
| --- | --- |
| Enough multi-return pulses with GPS time | `lines` fits position and velocity in each window using pulse endpoints. Pairing uses time, point source and scanner channel within each file. |
| Too many unpaired LAS/LAZ returns | The reader tries a 1 µs time tolerance for complete, unique return sequences on the same source and channel. |
| Usable scan angles and GPS time | `scan_angle` fits the angle off nadir, assuming a level line scanner. If pairing failed, the affected file goes to the single-return traversal. |
| Missing path windows in a file with usable single returns | `lines_scan_angle` keeps reliable line fits and uses scan angles to fill missing time windows. The combined estimate is marked `rough`. |
| A delivered path in the cloud's clock and CRS | Occlusion uses the delivered positions. An estimated trajectory is optional. |
| Too little geometric information for either estimator | The refusal message asks for the delivered trajectory. |

Line fits that miss their pulse lines by more than ten times the median, with a threshold of at least 0.25 m, are
dropped. The tool also drops positions over 1 m from the path described by nearby fits. The line quality labels are
`good` through a 0.1 m median miss and `rough` through 0.5 m. A line estimate above 0.5 m falls back to scan angles if available. Scan-angle fits are labelled `rough`. On the archived 3Dtrees comparison, most YellowScan multi-return paths were accurate to a few
centimetres. Scan-angle estimates had about 1 m horizontal error on RIEGL and DJI flights and about 8 m on
YellowScan.

The tool checks for static scans using the largest block by file, point source and time. It looks for a single
position that explains the pulse lines within 5 cm. For clouds without enough pulse lines, it checks whether
near-ground points concentrate around one location at least 15 times more than an even spread would.
`--no-platform-check` is a CLI diagnostic for known moving acquisitions. Static scan position estimation is outside
the path tool's scope.

## Resources and limits

- The reader decodes up to five million points per batch. The job caps decoding, path fitting and window workers
  at the configured `GALAXY_SLOTS` or `--threads`, CPU affinity and cgroup quota. The default is ten slots.
  Large files decode independent record ranges and merge them in scanner order. The path solver uses threads
  for pulse-line fits and separate processes for scan-angle fits. It merges both in time-window order.
  Compressed LAS 1.4 point formats 6–10 decode only
  coordinates, GPS time, return numbers, scanner channel, point source ID and scan angle during path
  estimation. Window tracing skips scan angle, and ordinary extent rereads use only positions and GPS time
  plus the LAZ base layer. Point formats 0–5 use full decoding.
- The defaults retain approximately 40 million multi-return and 10 million single-return pulses, or 40 million
  single-return pulses for single-return acquisitions. `--full` keeps every pulse and uses memory proportional
  to input size.
- GPS-window scratch contains durable compact counts at 12 bytes per stored voxel plus output intermediates.
  Count arrays are sparse files until the walker touches those cells. The debug report records recovered and unpaired
  echoes, executed and reused windows, per-window RSS, sampled aggregate RSS, allocated scratch, peak worker
  concurrency and memory reservations.
- `resources` records wall time, CPU, peak RSS, points read, decoder threads and cgroup measurements.
  Cgroup memory includes file cache.

- Positions cover sampled acquisition times. GPS week rolls, overlapping independent flights with the same
  clock and identifiers, and widely separated acquisitions should be split before processing.

## Galaxy resource proposal

Allocate ten CPU slots for trajectory and occlusion together. Start with the total standardized point count, then
check count cells in the ground-to-canopy envelope, the extent of the merge blocks and the number of raster bands.

| Input points | RAM (GiB) | Scratch (GiB) | Wall limit |
| ---: | ---: | ---: | ---: |
| Up to 75 million | 32 | 16 | 1 h |
| >75 million to 750 million | 64 | 96 | 2 h |
| >750 million to 2.1 billion | 96 | 256 | 5 h |
| >2.1 billion to 7.1 billion | 96 | 256 | 6 h |

The standardized LAZ replay for 2027 (29.7 M points) took 78 s and peaked at 3.0 GiB of cgroup memory.
Historical scanner-order LAZ benchmarks provide headroom estimates for the larger classes:

| Dataset | Points | Read stage (s) | Full job (s) | Peak cgroup RAM (GiB) | Peak scratch (GiB) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1807 | 52.7 M | 9.0 | 101 | 2.7 | 2.2 |
| 2014, multi-return replay | 52.7 M | 10.8 | 152 | 21.8 | 6.5 |
| 2014, single-return replay | 33.8 M | 6.4 | 153 | 24.0 | 7.5 |
| 1844 | 678.1 M | 101.4 | 994 | 34.5 | 32.5 |
| 2585 | 7.04 B | 1,238.8 | 27,953 | 190.0 | 567.9 |

The 2014 single-return replay reached 24.01 GiB, so the smallest class allows 32 GiB. The 1844 replay
fits the 64 GiB class with room for the parent process and file cache. The 750 M–2.1 B class is interpolated.
The round-10 raw 2585 replay took 2 h 11 min with ten CPU slots, a 50 GiB cap and about 50 GiB of
allocated scratch. A replay from completed checkpoints measured 1 min 13 s for the final merge.
The resource class is an estimate for acquisitions with a comparable ground-to-canopy footprint.
A standardized 2585 replay is still required for Galaxy validation.

The standardized 1788 collection has 24.7 M points and a 2,000 × 2,000 × 1,278 grid at 0.5 m.
Use `--debug` for current window memory and scratch measurements. [Stage resources](docs/REPORT.md#3-stage-resources)

## Archived trajectory-only replay

A replay of 36 archived 3Dtrees ULS/MLS datasets matched the time and XYZ bytes of all 35 prior paths.
Dataset 2018 reproduced its refusal because it has single returns without a usable scan angle.
The archive and replay ran separately on geos. Replay jobs used eight decoder slots with a ten-CPU affinity mask.

| Dataset | Points read | Archived → replay wall time | Archived → replay peak RSS |
| --- | ---: | ---: | ---: |
| 1807 YellowScan | 52.7 M | 23.14 → 11.28 s | 2.21 → 1.06 GB |
| 2018 MLS | 1.95 B | 188.48 → 161.52 s | 6.86 → 2.70 GB |
| 2585 DJI | 7.04 B | 646.63 → 712.40 s | 3.57 → 2.34 GB |

Across all 36 jobs the median wall time was 17.06 s, compared with 26.97 s in the archive. Peak replay process
RSS was 3.56 GB. The largest input charged 77.30 GB to its cgroup, which also accounted for file cache.
For 2585, the delivered-path median horizontal error remains 37.8 m. Quality labels measure pulse fit residuals.

## Development and container

```bash
make test
make image
```

The Dockerfile pins both base image digests, the Untwine and system Expat versions, and the Linux Python runtime packages.

The Docker suite runs with two CPUs and a 4 GiB limit by default. Tests generate small LAS/LAZ clouds with known positions and verify the coordinates, refusals, resource reporting,
and byte-for-byte stability across decoder thread counts.

Call `python /src/run.py` from a writable job directory. The Docker image defaults to a non-root user. For a
local bind-mounted job directory, run under the job owner's UID and GID.

```bash
docker run --rm --name trajectory-local --user "$(id -u):$(id -g)" \
  --label "owner=$(id -un)" --label created-by=trajectory-local \
  --mount "type=bind,source=$PWD,target=/data" \
  3dtrees_trajectory:dev python /src/run.py standardized/strip1.las standardized/strip2.laz \
    -o trajectory/trajectory.csv
```

The Galaxy wrapper should accept one collection of standardized LAS/LAZ files per dataset and pass all its
members to one tool job. Pass `--threads "${GALAXY_SLOTS:-4}"` from the job allocation and collect the requested
CSV, GeoJSON and occlusion products.
Keep sampling and platform diagnostics in the CLI. File metadata and geometry refusals are checked at runtime.
CI tests the image on pull requests, pushes to `main`, and pushes of `v*` tags. After the tests pass, it
publishes `ghcr.io/3dtrees-earth/3dtrees_trajectory:sha-<short-sha>` and either `:latest` for `main` or
the package version (for example, `:1.0.0`) for a matching release tag. Galaxy converts that Docker image
to Singularity when preparing the tool container.

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
