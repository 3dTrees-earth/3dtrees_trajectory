# 3dtrees_trajectory

Estimate the flight or walking path of a moving lidar sensor (ULS, MLS, ALS) from LAS or LAZ point clouds. The tool
writes a 3Dtrees trajectory CSV, a quality and resource report, and a map layer. It can also map occlusion along
the estimated or delivered path. The current package version is
0.1.0.

## Run

```bash
3dtrees-trajectory flight.laz -o trajectory/trajectory.csv
3dtrees-trajectory strip1.las strip2.laz -o trajectory/trajectory.csv
```

Combine strips from the same acquisition and use a consistent file order. The reader supports LAS/LAZ 1.2–1.4,
point formats 0–10, and COPC. Inputs need finite, varying GPS times in the same GPS time system and a shared metric
coordinate system. When the CRS is missing, the tool assumes metre units and keeps the input coordinates. Inputs in
formats 0 and 2 have no GPS time and receive a refusal report.

### Outputs

| File | Contents |
| --- | --- |
| `trajectory.csv` | `3dtrees-trajectory/2` with columns `time,x,y,z,roll,pitch,yaw,name`. Time and XYZ use the cloud's clock and coordinates. Attitude and name are empty. |
| `trajectory.report.json` | Method, quality, reasons for refusing the input, point counts, and job resources. |
| `trajectory.geojson` | Path segments separated at gaps over 5 seconds, thinned to about one vertex per metre. Uses WGS84 when the CRS is known and cloud coordinates otherwise. |

`--report` and `--geojson` override the default paths next to the CSV. The tool creates their parent directories
and publishes each file by rename. The GeoJSON file is zero bytes when there are no drawable segments or
reprojection fails. `3dtrees-trajectory-geojson trajectory/trajectory.csv` can also draw the layer from an
existing trajectory CSV.

### Exit codes

| Code | Meaning |
| --- | --- |
| `0` | The path was estimated or the delivered-path occlusion products were written. In a combined job, the estimated path may succeed with an occlusion error in the report. |
| `2` | Invalid arguments, unreadable or unsuitable cloud, static scan, no reliable fit, or an output write failed. The JSON report records the reason when the report path is available and writable. Safe, writable CSV and GeoJSON destinations receive zero-byte files. |

Supply one to 65,536 per-file LAS/LAZ strips from one acquisition, or one merged COPC (`.copc.laz`).
Supply an output CSV for an estimated path. The CSV, GeoJSON, report and input paths must be distinct.
If you set `--interval`, use a positive finite number of seconds (default 0.5). An interval must be long enough for
the cloud's GPS times to fit into a signed 64-bit window index. `--threads` takes a positive integer.

## Occlusion

```bash
3dtrees-trajectory flight.laz -o trajectory.csv \
  --occlusion-copc occlusion.copc.laz --occlusion-cog occlusion.tif

3dtrees-trajectory flight.laz --trajectory delivered.csv --report occlusion.report.json \
  --occlusion-copc occlusion.copc.laz --occlusion-cog occlusion.tif
```

Supply both occlusion destinations. When `--trajectory` is supplied, occlusion uses that path. With no
`-o` the tool checks the cloud and writes the occlusion products and report without fitting a path.
The delivered CSV starts with `time,x,y,z`, optionally preceded by `#` metadata. It needs at least two
positions, finite coordinates, and unique timestamps in the cloud's GPS clock and coordinate system.
Resolve week seconds versus adjusted GPS time before running. Both modes require a moving sensor.

The compiled OccPy 0.3.0 traversal linearly interpolates the sensor position for each pulse. Echoes up
to one second outside the path take the nearest endpoint. More distant echoes are counted and skipped.
The tool refuses disjoint time ranges. A voxel is filled (1) if there is a registered return, empty (2)
if seen clear, occluded (3) if only traversed behind returns, or unobserved (4) in the absence of traversal.

Voxel size defaults to 0.5 m. `--voxel` accepts 0.05–10 m. Ground is the centre of the lowest filled voxel
per column. Gaps take the nearest available ground estimate. Ground more than 2 m above its 5 m
neighbourhood median is replaced by that median. With `--dtm terrain.tif`, ground comes from band 1 at
each column centre. The DTM must cover part of the cloud with valid elevations and share its CRS if both
declare one. Uncovered columns use the nearest sampled elevation. The highest filled voxel caps each
column.

| File | Contents |
| --- | --- |
| `occlusion.copc.laz` | One point at each filled, occluded or unobserved voxel centre from ground through canopy top. Extra dimensions are `occlusion`, `nhit`, `nmiss`, `nocc`, `occluded_fraction` and `height_above_ground`. The fraction is `nocc / (nhit + nmiss + nocc)`, with NaN for unobserved voxels. |
| `occlusion.tif` | North-up COG with one pixel per voxel column. The 1 m bands are named `occluded_0-1m`, `occluded_1-2m` and so on. Each value is the occluded voxel count divided by the sum of filled, empty and occluded voxels, with NaN when the sum is zero. The final `unobserved_column` band gives the unobserved share of all voxels from ground through canopy top. |

If there are no qualifying voxels, the COPC is zero bytes and the COG is valid with NaN values. On an
occlusion error both destinations are zero bytes and the report gives the reason. A combined job can
still return its estimated trajectory. The report includes traversal counts, in-range class totals,
and the pooled occluded share per height band.

A first pass derives the lowest and highest in-range point voxel in each column. Traversal keeps the vertical
band from estimated ground minus 1 m to canopy top plus 1 m, on the original voxel lattice. The products cap
each column at its highest registered hit. If the estimate for the cropped grid exceeds 16 GiB
(`--grid-memory-gib`), the tool traces buffered XY tiles and writes their nonoverlapping cores to the same products.
Each tile uses its local ground and canopy envelope to size the vertical grid.
It estimates the buffer from the 99.9th percentile of horizontal ray distance through the vegetation band on a
deterministic 1% pulse sample, rounds to a voxel boundary, and caps the applied halo at 100 m. An unworkable grid under the
64 GiB planning ceiling fails with a report. COPC tiles use octree spatial reads. LAS/LAZ tiles scan the strips in
their original order. Rays longer than the buffer and non-collinear returns can change results near tile seams.

The delivered path is used for every occlusion ray when supplied. The report and COG metadata flag a median
perpendicular distance above 0.1 m from the delivered path to usable pulse lines. With only single returns there
are no return lines for this check. OccPy traces single-return shots independently. An occlusion run from a
scan-angle estimate is marked `rough` in both outputs.

COPC indexing and COG conversion use temporary files beside the product destinations. The job needs space for
the uncompressed voxel LAS, intermediate TIFF and final outputs. The image includes Untwine 1.5.1. For local
runs, set `TRAJECTORY_UNTWINE` to an executable installation and install the compiled OccPy wheel with
`python -m pip install --no-deps occpy_ls==0.3.0`.

## Platform pre-gate

Use the standardization metadata to check that GPS time is present and varying. Treat
`max(NumberOfReturns) > 1` as evidence of multi-returns. A single-return cloud needs usable scan angles or an
available delivered path for occlusion. Check for a scan-angle dimension and actual nonzero values, since a LAS
point format can contain an angle field filled with zero. Multi-return clouds can estimate a path without angles
when they have enough usable pulse lines. The tool checks usable pulses, the multi-return share, return pairing,
and static-scan evidence and records its refusal in JSON. The platform gate can route clearly unsuitable
datasets before scheduling a job.

## How it works

| Input evidence | Method |
| --- | --- |
| Enough multi-return pulses with GPS time | `lines` fits the position and velocity in each time window from the line through each pulse's first and last returns. |
| Usable scan angles and GPS time | `scan_angle` uses each point's angle off nadir. It assumes a line scanner in level flight. |
| Too little geometric information for either method | The report asks for the delivered trajectory. |

Line fits that miss their pulse lines by more than ten times the median, with a threshold of at least 0.25 m, are
dropped. The tool also drops positions over 1 m from the path described by nearby fits. The line quality labels are
`good` through a 0.1 m median miss and `rough` through 0.5 m. A line estimate above 0.5 m is refused
as unreliable. Scan-angle fits are labelled `rough`. On the archived 3Dtrees comparison, most YellowScan multi-return paths were accurate to a few
centimetres. Scan-angle estimates had about 1 m horizontal error on RIEGL and DJI flights and about 8 m on
YellowScan.

The tool checks for static scans using the largest block by file, point source and time. It looks for a single
position that explains the pulse lines within 5 cm. For clouds without enough pulse lines, it checks whether
near-ground points concentrate around one location at least 15 times more than an even spread would.
`--no-platform-check` is a CLI diagnostic for known moving acquisitions. Static scan position estimation is outside
the path tool's scope.

## Resources and limits

- The reader decodes up to five million points per batch and selects the LAZ fields needed by the solver where
  the point format permits. Non-COPC LAZ decoding uses `GALAXY_SLOTS`, or a lower `--threads` limit, capped by
  CPU affinity and the cgroup v2 quota. The default is eight. COPC sequential decoding and numerical reductions
  are serial.
- The COPC input is externally sorted by GPS time, point source, scanner channel and return number in the job
  working directory. Temporary runs use about 30 bytes per input point. When tiling, the acquisition sort
  is released before sorting the active COPC tile. The tool removes the temporary runs on exit. GPS time groups stay together across disk runs
  and read chunks. A merged COPC must preserve return numbers, GPS time and source/channel identity. Duplicate
  multi-return numbers with the same time, source and channel are ambiguous.
- Per-file LAS/LAZ retains the original point order and avoids the time sort. The trajectory reader can recover
  complete pulses whose returns were stored apart, within a bounded sample. Traversal can rejoin up to 40
  million separated echoes per input strip in bounded batches and refuses larger fragment sets. It reports recovered and skipped
  echoes. Heavy thinning, return renumbering, and
  clock resets can remove the information needed to reconstruct pulses.
- A timestamp hash samples the acquisition independently of file order. Default budgets are approximately
  40 million multi-return and 10 million single-return pulses, or 40 million single-return pulses in a cloud
  that starts with very few multi-returns. `--full` uses memory proportional to the input. The solver uses
  float32 endpoints relative to a common kilometre-aligned origin.
- `resources` records wall time, CPU, peak RSS, points read, decoder threads and cgroup measurements when
  available. `stages` records per-stage wall, self and child CPU, peak process RSS, and bytes when available.
  Linux resets the RSS high-water mark between stages. Cgroup memory includes file cache, so use a dedicated
  cgroup for each job when measuring.
- Positions cover sampled acquisition times. GPS week rolls, overlapping independent flights with the same
  clock and identifiers, and widely separated acquisitions should be split before processing.

## Container replay

The runtime source was mounted read-only into the trajectory image with Untwine 1.5.1. Each replay had
eight decoder slots and a ten-CPU affinity mask. The wall values are from the instrumented tool. Cgroup peaks
include cache. For tiled jobs, decoded points, traversed pulses and other counters are summed across buffered
tiles, so one point can contribute to multiple tiles.

| Dataset | Input | Points | Grid voxels | Tiles | Wall (s) | Peak RSS (GiB) | Cgroup peak (GiB) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1807 YellowScan | LAZ | 52.7 M | 15.8 M | 1 | 96.9 | 2.39 | 2.33 |
| 1807 YellowScan | COPC | 52.7 M | 15.8 M | 1 | 144.0 | 3.66 | 3.64 |
| 1807 YellowScan | LAZ | 52.7 M | 15.8 M | 2 | 122.4 | 2.36 | 2.75 |
| 1807 YellowScan | COPC | 52.7 M | 15.8 M | 2 | 176.9 | 3.16 | 3.19 |
| 2014 DJI L3 | LAZ | 52.7 M | 773.2 M | 9 | 451.9 | 6.61 | 6.76 |
| 2014 DJI L3 | COPC | 52.7 M | 773.2 M | 9 | 452.8 | 9.30 | 10.18 |

## Galaxy resource proposal

Use the cropped buffered tile size for memory planning. The tool estimates `2 GiB + 72 × tile voxels
+ 320 × tile XY columns + 32 × acquisition XY columns` bytes, plus a reserve for separated LAZ echoes. Allow separate headroom for the path solver,
the job runtime and COPC file cache. COPC acquisition sorting also needs about 30 GB of working disk per
billion points. With tiled LAZ, budget time for a complete scan per tile. With tiled COPC, budget for one
acquisition sort and spatial reads and sorts for the buffered tiles.

A candidate TPV rule would start at four slots for small jobs and eight for large LAZ decoding or Untwine
indexing. Estimate wall time from total decoded points and total allocated buffered voxels. Estimate memory
from the largest tile and the measured solver RSS. Validate both against dedicated job cgroup peaks and add
site-specific headroom before setting job classes. The exact Galaxy memory classes and the platform's default
input representation are deployment decisions.

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
uv venv .venv
uv pip install -e '.[test]'
.venv/bin/python -m pytest -q
docker build -t 3dtrees_trajectory:dev .
```

Tests generate small LAS/LAZ clouds with known positions and verify the coordinates, refusals, resource reporting,
and byte-for-byte stability across decoder thread counts.

Call `python /src/run.py` from a writable job directory. The Docker image defaults to a non-root user. For a
local bind-mounted job directory, run under the job owner's UID and GID.

```bash
docker run --rm --name trajectory-local --user "$(id -u):$(id -g)" \
  --label "owner=$(id -un)" --label created-by=trajectory-local \
  --mount "type=bind,source=$PWD,target=/data" \
  3dtrees_trajectory:dev python /src/run.py flight.laz -o trajectory/trajectory.csv
```

A future Galaxy wrapper can pass the selected cloud files and output paths, with a positive finite interval if
exposed. Pass `--threads "${GALAXY_SLOTS:-4}"` from the job allocation and collect the three fixed outputs.
Keep sampling and platform diagnostics in the CLI. File metadata and geometry refusals are checked at runtime.
A published release tag that matches the CLI and package versions publishes
`ghcr.io/3dtrees-earth/3dtrees_trajectory:<version>` and a SIF built from the same Docker digest at
`oras://ghcr.io/3dtrees-earth/3dtrees_trajectory:<version>-sif`.
