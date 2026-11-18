# 3dtrees_trajectory

Estimate the flight or walking path of a moving lidar sensor (ULS, MLS, ALS) from LAS or LAZ point clouds. The tool
writes a 3Dtrees trajectory CSV, a quality and resource report, and a map layer. The current package version is
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
| `0` | At least one valid position was estimated. |
| `2` | Invalid arguments, unreadable or unsuitable cloud, static scan, no reliable fit, or an output write failed. The JSON report records the reason when the report path is available and writable. Safe, writable CSV and GeoJSON destinations receive zero-byte files. |

Supply one to 65,536 LAS/LAZ inputs and an output CSV path. The CSV, GeoJSON, report and input paths must be distinct.
If you set `--interval`, use a positive finite number of seconds (default 0.5). An interval must be long enough for
the cloud's GPS times to fit into a signed 64-bit window index. `--threads` takes a positive integer.

## How it works

| Input evidence | Method |
| --- | --- |
| Enough multi-return pulses with GPS time | `lines` fits the position and velocity in each time window from the line through each pulse's first and last returns. |
| Usable scan angles and GPS time | `scan_angle` uses each point's angle off nadir. It assumes a line scanner in level flight. |
| Too little geometric information for either method | The report asks for the delivered trajectory. |

Line fits that miss their pulse lines by more than ten times the median, with a threshold of at least 0.25 m, are
dropped. The tool also drops positions over 1 m from the path described by nearby fits. The line quality labels are
`good` through a 0.1 m median miss, `rough` through 0.5 m, and `unreliable` above 0.5 m. Scan-angle fits are
labelled `rough`. On the archived 3Dtrees comparison, most YellowScan multi-return paths were accurate to a few
centimetres. Scan-angle estimates had about 1 m horizontal error on RIEGL and DJI flights and about 8 m on
YellowScan.

The tool checks for static scans using the largest block by file, point source and time. It looks for a single
position that explains the pulse lines within 5 cm. For clouds without enough pulse lines, it checks whether
near-ground points concentrate around one location at least 15 times more than an even spread would.
`--no-platform-check` is a CLI diagnostic for known moving acquisitions. Static scan position estimation is outside
the path tool's scope.

## Resources and limits

- The reader decodes five million points at a time and selects only the LAZ fields used by the solver when the
  point format permits it. Non-COPC LAZ decoding uses `GALAXY_SLOTS` (default 8 locally) or the lower
  `--threads` limit, capped by CPU affinity and available cgroup v2 CPU quota. COPC decoding and numerical
  library reductions are serial.
- The default retains roughly up to 40 million multi-return pulses and 10 million single returns, or 40 million
  single returns when the first block looks single-return. The tool thins by stride across the acquisition.
  `--full` retains all pulses and requires memory proportional to the cloud size. The path fit keeps the selected
  pulse columns and processes one fitting window at a time. Endpoint coordinates use float32 relative to the first
  nonempty file's LAS offset. Widely separated clouds reduce sub-metre precision.
- The report's `resources` object includes `points_read`, `wall_seconds`, `cpu_seconds`, `peak_rss_bytes`
  and `decoding_threads`. On cgroup v2 it also records current and peak memory, the local limit, and the change
  in CPU usage. Cgroup values describe the entire cgroup, including charged file cache and other processes in a
  shared cgroup. Use a dedicated job/container cgroup for per-run measurements.
- Original LAS/LAZ strips retain the clearest pulse order. The tool can join separated returns by file, point
  source and GPS time, up to 40 million retained stray records. Shared timestamps in multibeam or COPC exports
  can mix pulses, and some pulses spanning read-chunk boundaries can be lost. Thinning or resampling upstream can
  remove pulse structure.
- Paths cover the times represented by points. Turns or transit legs omitted from the cloud leave gaps. GPS week
  rolls and files from different acquisitions should be processed separately.

## Local replay

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
