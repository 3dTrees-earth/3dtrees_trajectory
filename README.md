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

The CLI decodes each source point once for the trajectory sample and occlusion. It reads COPC spatial
nodes and records compact echoes in scratch. From those records, it derives ground and canopy extents and routes
whole pulses to each buffered XY tile they touch. Separated returns are joined by GPS time, source and scanner
channel within each input file. Each tile uses its own local ground and canopy envelope on the original voxel
lattice. Traversal keeps the vegetation band from ground minus 1 m to canopy top plus 1 m. The products cap
each column at its highest registered hit.

The tool tiles when the cropped grid exceeds the per-tile planning budget (`--grid-memory-gib`, default 16 GiB)
or a job needs to keep a grid within OccPy's 32-bit echo limits. The buffer uses the 99.9th percentile of horizontal
ray distance through the vegetation band on a deterministic 1% pulse sample. It rounds to the voxel boundary
and caps the halo at 100 m. Any return inside the tile's buffered bounds brings the whole pulse into that tile.
The tool writes each tile's nonoverlapping core to the products. Rays longer than the buffer can change results
near tile seams.

Separate worker processes trace tiles, with one OccPy grid per worker. The tool schedules tiles within the
job's CPU and memory budgets, giving priority to tiles with more routed echoes. It commits each tile's core
counts to scratch before merging. To resume
after an interruption, give the command the same inputs, parameters and output paths. Keep the
`.occlusion-*` directory beside the COPC until the merge finishes. The tool cleans up the checkpoint
after publishing both products.

The delivered path is used for every occlusion ray when supplied. The report and COG metadata flag a median
perpendicular distance above 0.1 m from the delivered path to usable pulse lines. With only single returns there
are no return lines for this check. OccPy traces single-return shots independently. An occlusion run from a
scan-angle estimate is marked `rough` in both outputs.

COPC indexing and COG conversion use temporary files beside the product destinations. The job needs space for
the uncompressed voxel LAS, intermediate TIFF and final outputs. Wide rasters use a band-interleaved
intermediate and up to four job slots for the COG conversion. The image includes Untwine 1.5.1. For local
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

- The reader decodes up to five million points per batch. The job uses up to `GALAXY_SLOTS` or the
  lower `--threads` limit for decoding and tile workers, also capped by CPU affinity and the cgroup v2
  quota. The default is eight slots. The path solver uses single-threaded numerical reductions.
- COPC trajectory estimation orders the sampled GPS-time groups by time, source, scanner channel and
  return number. Its temporary sort uses 38 bytes per candidate echo, including contiguous time keys.
  Pulse population counts from thinned COPC inputs are estimates and carry a flag in the report.
  The defaults retain approximately 40 million multi-return and 10 million single-return pulses,
  or 40 million single-return pulses for single-return acquisitions. `--full` explicitly samples
  every pulse and uses memory proportional to input size.
- Occlusion routes 24-byte echoes from the source recording into buffered tile spools. A bounded hash
  partition joins separated echoes across input batches. Workers order each tile by time in bounded
  windows. For acquisitions above the native 32-bit echo limit, tile spools use lossless compressed
  blocks to limit scratch use as buffers duplicate echoes. Large LAS/LAZ inputs route separate strips
  concurrently, and large COPC inputs join independent hash partitions concurrently. Both use the job
  caps. Ambiguous or incomplete return sequences are reported as skipped echoes.
  Each durable tile uses 12 bytes per core voxel for the three count arrays until merge.
- `resources` records wall time, CPU, peak RSS, points read, decoder threads and cgroup measurements
  when available. `stages` records wall, self and child CPU, peak process RSS and I/O where available.
  Traversal also reports concurrent workers, reused tiles, sampled aggregate RSS and per-tile timings.
  Cgroup memory includes file cache.
- Positions cover sampled acquisition times. GPS week rolls, overlapping independent flights with the same
  clock and identifiers, and widely separated acquisitions should be split before processing.

## Sequential container replay baseline

The earlier runtime source was mounted read-only into the trajectory image with Untwine 1.5.1. Each replay had
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
| 1844 YellowScan | LAZ | 678.1 M | 203.6 M | 2 | 1440.1 | 7.58 | 14.96 |
| 1844 YellowScan | COPC | 678.1 M | 203.6 M | 2 | 1901.9 | 6.24 | 32.57 |


The 7.04-billion-point 2585 replay produced trajectories from per-file LAZ and merged COPC. The LAZ
occlusion step refused a strip with more than 40 million separated echoes in a tile. The COPC process
was stopped in its fifth tile after eleven hours at the operator’s request.

## Parallel container replay

The source scan and tile worker replays used at most ten slots and 190 GiB of Docker memory. The table
compares tool wall times against the sequential runs above. Peak process RSS is the largest individual
process measurement. Cgroup peaks include the resident caches for the compact source and tile spools.

| Acquisition | Wall (s) | Speedup | Tiles / peak workers | Peak process RSS (GiB) | Cgroup peak (GiB) | Scratch peak (GiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1807 LAZ | 101.3 | 0.96× | 1 / 1 | 1.27 | 2.73 | 2.17 |
| 1807 COPC | 120.7 | 1.19× | 1 / 1 | 1.27 | 3.23 | 2.31 |
| 1844 LAZ | 994.3 | 1.45× | 2 / 2 | 5.22 | 34.49 | 32.49 |
| 1844 COPC | 1294.6 | 1.47× | 2 / 2 | 5.22 | 39.12 | 30.82 |
| 2014 multi LAZ | 152.1 | 2.97× | 8 / 8 | 3.96 | 21.75 | 6.48 |
| 2014 single LAZ | 153.3 | 1.91× | 4 / 4 | 7.83 | 24.01 | 7.53 |
| 2585 LAZ | 27953.0 | — | 88 / 10 | 7.33 | 190.00 | 567.87 |
| 2585 COPC | 20169.5 | — | 54 / 9 | 9.89 | 190.00 | 460.79 |

Both 1807 replays and the 1844 LAZ replay reproduce the earlier products exactly. The new 1844 COPC
reproduces the earlier LAZ voxel classes, per-voxel counts and COG values. Its trajectory CSV is
byte-identical to the earlier COPC CSV. The earlier COPC clipped some siblings at the outer tile
buffer. Including complete pulses changes 19 values relative to that COG and 16 classes on shared
emitted voxel centres.

The 2585 LAZ replay decoded 7.04 billion points and rejoined 3.31 billion separated echoes. Its
compressed tile spools occupied 382 GB for 916 GB of routed records. The measured run spent 6,204 s in the
earlier 64-band COG conversion. In a separate Docker probe, the new band-interleaved intermediate
converted to a COG in 51 s and matched every base raster value.

The 2585 COPC replay also completed, with 307.7 GB of compressed tile spools and a 51 s COG
conversion. It started during the LAZ COG conversion and kept the nine traversal slots set at
occlusion startup. Both 2585 trajectory CSVs match their archived outputs byte for byte.
The speedup column uses completed historical replays.

## Galaxy resource proposal

Budget each running tile at `2 GiB + 72 × buffered tile voxels + 320 × buffered XY columns`.
Add the reservations for all concurrently running tiles, the parent's resident memory and headroom
for the runtime and file cache. Plan the path solver separately since its main allocations precede
traversal. The scheduler uses the smallest applicable cgroup memory limit and the available host memory.

For small or single-tile jobs, start with four slots for decoding and indexing. For tiles large enough
to dominate wall time, start with eight slots and allow ten for very large acquisitions when memory can
support that many active grids. Estimate traversal wall time from the slowest group of tiles running
concurrently, using buffered echo counts and measured per-tile CPU. Size decoding, routing and product
merge separately using source point count, routed echo count and exported voxel count.

Allow scratch for the 24-byte source recording, the echoes copied to tile spools, temporary
partition data, and any sampled COPC sort still active. During merge, allow for the durable count
arrays, uncompressed voxel LAS, the intermediate TIFF, and the outputs. Measure both peak cgroup
memory and peak scratch in a dedicated job before choosing Galaxy resource classes.

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
