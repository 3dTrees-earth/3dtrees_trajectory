# 3dtrees_trajectory

The 1.0.0 scope covers ULS and MLS. The tool estimates the flight or walking path from a standardized point cloud
collection and writes a 3Dtrees trajectory CSV, a quality and resource report, and a map layer. It can also map
occlusion along the estimated or delivered path. The current package version is 0.1.0.

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
metres and preserves the input coordinates. Point formats 0 and 2 receive a refusal report because
they lack GPS time.

### Outputs

| File | Contents |
| --- | --- |
| `trajectory.csv` | `3dtrees-trajectory/2` with columns `time,x,y,z,roll,pitch,yaw,name`. Time and XYZ use the cloud's clock and coordinates. Attitude and name are empty. |
| `trajectory.report.json` | Method, quality, per-file return modes, exit code, refusal reason, point counts, and stage and job resources. |
| `trajectory.geojson` | Path segments separated at gaps over 5 seconds, thinned to about one vertex per metre. Uses WGS84 when the CRS is known and cloud coordinates otherwise. |

`--report` and `--geojson` override the default paths next to the CSV. The tool creates their parent directories
and publishes each file by rename. The GeoJSON file is zero bytes when there are no drawable segments or
reprojection fails. `3dtrees-trajectory-geojson trajectory/trajectory.csv` can also draw the layer from an
existing trajectory CSV.

### Exit codes

| Code | Meaning |
| --- | --- |
| `0` | The path was estimated or the delivered-path occlusion products were written. In a combined job, the estimated path may succeed with an occlusion error in the report. |
| `2` | Invalid arguments or data, a detected static (TLS) acquisition, or a failed output operation. |
| `3` | GPS time is missing, nonfinite or constant in a moving acquisition. A constant point clock cannot be aligned with a moving delivered path. |
| `4` | The cloud has no usable estimate and no delivered path to support occlusion. |

For a refusal, the JSON records the reason when the report destination is writable. Safe, writable product
destinations receive zero-byte files.

Supply an output CSV for an estimated path. The CSV, GeoJSON, report and input paths must be distinct.
If you set `--interval`, use a positive finite number of seconds (default 0.5). An interval must be long enough for
the cloud's GPS times to fit into a signed 64-bit window index. `--threads` takes a positive integer.

## Occlusion

```bash
3dtrees-trajectory standardized/strip1.las standardized/strip2.laz -o trajectory.csv \
  --occlusion-copc occlusion.copc.laz --occlusion-cog occlusion.tif

3dtrees-trajectory standardized/strip1.las standardized/strip2.laz --trajectory delivered.csv --report occlusion.report.json \
  --occlusion-copc occlusion.copc.laz --occlusion-cog occlusion.tif
```

Supply both occlusion destinations. When `--trajectory` is supplied, occlusion uses that path. With no
`-o` the tool checks the cloud and writes the occlusion products and report without fitting a path.
A delivered moving-path CSV starts with `time,x,y,z`, optionally preceded by `#` metadata. It needs at least two
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

The CLI normally reads the standardized collection once for both the trajectory sample and occlusion. It records
compact echoes in scratch. Testing a fallback can read an affected file again. From those records, it derives
ground and canopy extents and routes
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
when they have enough usable pulse lines. Preserve GPS time, ReturnNumber, NumberOfReturns, ScanAngle
(or ScanAngleRank for legacy formats), PointSourceID and scanner channel wherever the input format provides
them. The tool checks usable pulses, the multi-return share, return pairing,
and static-scan evidence and records its refusal in JSON. The platform gate can route clearly unsuitable
datasets before scheduling a job.

## How it works

| Input evidence | Method |
| --- | --- |
| Enough multi-return pulses with GPS time | `lines` fits position and velocity in each window using pulse endpoints. Pairing uses time, point source and scanner channel within each file. |
| Too many unpaired LAS/LAZ returns | The reader tries a 1 µs time tolerance for complete, unique return sequences on the same source and channel. |
| Usable scan angles and GPS time | `scan_angle` fits the angle off nadir, assuming a level line scanner. If pairing failed, the affected file goes to the single-return traversal. |
| Missing path windows in a file with usable single returns | `lines_scan_angle` keeps reliable line fits and uses scan angles to fill missing time windows. The combined estimate is marked `rough`. |
| A delivered path in the cloud's clock and CRS | Occlusion uses the delivered positions. An estimated trajectory is optional. |
| Too little geometric information for either estimator | The report asks for the delivered trajectory. |

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

- The reader decodes up to five million points per batch. The job uses up to `GALAXY_SLOTS` or the
  lower `--threads` limit for decoding and tile workers, also capped by CPU affinity and the cgroup v2
  quota. The default is eight slots. The path solver uses single-threaded numerical reductions.
- LAZ uses parallel lazrs decoding with bounded chunk reads in original record order.
  The defaults retain approximately 40 million multi-return and 10 million single-return pulses,
  or 40 million single-return pulses for single-return acquisitions. `--full` explicitly samples
  every pulse and uses memory proportional to input size.
- Occlusion routes 24-byte echoes from the source recording into buffered tile spools. A bounded hash
  partition joins separated echoes across input batches. Workers order each tile by time in bounded
  windows. For acquisitions above the native 32-bit echo limit, tile spools use lossless compressed
  blocks to limit scratch use as buffers duplicate echoes. Large LAS/LAZ inputs route separate strips
  concurrently within the job caps. The report records unpaired echoes and the files traced in single-return mode.
  Each durable tile uses 12 bytes per core voxel for the three count arrays until merge.
- `resources` records wall time, CPU, peak RSS, points read, decoder threads and cgroup measurements
  when available. `stages` records wall, self and child CPU, peak process RSS and I/O where available.
  Traversal also reports concurrent workers, reused tiles, sampled aggregate RSS and per-tile timings.
  Cgroup memory includes file cache.
- Positions cover sampled acquisition times. GPS week rolls, overlapping independent flights with the same
  clock and identifiers, and widely separated acquisitions should be split before processing.

## Galaxy resource proposal

Allocate ten CPU slots for trajectory and occlusion together. Start with the total standardized point count, then
check the buffered grid, the number of times each echo is routed, and the number of raster bands.

| Input points | RAM (GiB) | Scratch (GiB) | Wall limit |
| ---: | ---: | ---: | ---: |
| Up to 75 million | 32 | 16 | 1 h |
| >75 million to 750 million | 64 | 96 | 2 h |
| >750 million to 2.1 billion | 96 | 256 | 5 h |
| >2.1 billion to 7.1 billion | 140 | 750 | 10 h |

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
The raw 2585 replay reached its 190 GiB cap including file cache and used 34.5 GiB of sampled aggregate RSS
in traversal. The proposed 140 GiB class needs a full standardized 2585 run at that cap before deployment.

The standardized 1788 collection had 24.7 M points and a 2,000 × 2,000 × 1,278 grid at 0.5 m. It finished
123 tiles and reached 76.5 GiB of cgroup memory including cache. Start comparable wide grids in the 96 GiB
class. Budget each active tile at `2 GiB + 72 × buffered tile voxels + 320 × buffered XY columns`. Sum
concurrent reservations and leave room for the parent and file cache. Scratch covers the 24-byte source
records, buffered echo copies, hash partitions, durable tile counts and product intermediates.

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

The Dockerfile pins both base image digests, the Untwine and system Expat versions, and the Linux Python runtime packages.

Tests generate small LAS/LAZ clouds with known positions and verify the coordinates, refusals, resource reporting,
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
members to one tool job. Pass `--threads "${GALAXY_SLOTS:-4}"` from the job allocation and collect the three fixed outputs.
Keep sampling and platform diagnostics in the CLI. File metadata and geometry refusals are checked at runtime.
A published release tag that matches the CLI and package versions publishes
`ghcr.io/3dtrees-earth/3dtrees_trajectory:<version>` and a SIF built from the same Docker digest at
`oras://ghcr.io/3dtrees-earth/3dtrees_trajectory:<version>-sif`.
