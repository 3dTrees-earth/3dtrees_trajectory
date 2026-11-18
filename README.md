# 3dtrees_trajectory

Estimates the flight or walking path of a moving lidar sensor (ULS, MLS, ALS) from the point cloud alone. The output
is the 3Dtrees trajectory CSV (`3dtrees-trajectory/2`), the layout delivered trajectories use on the platform, plus a
JSON quality report and a GeoJSON map layer.

```bash
3dtrees-trajectory flight.laz -o flight_trajectory.csv
3dtrees-trajectory strip1.laz strip2.laz strip3.laz -o flight_trajectory.csv
```

Static (TLS) scans are refused in this version. The CSV layout can hold scan positions (`kind: stations`), so a later
version can add them without a format change.

## How it works

| Cloud | Method | Accuracy on 3Dtrees data |
| --- | --- | --- |
| Multi-return pulses, GPS time | `lines`: a pulse's first and last return lie on a line through the sensor; per 0.5 s window, the position and velocity closest to all lines | under 0.5 m on 31 of 34 delivered trajectories, most YellowScan flights about 0.01 m |
| Single returns, scan angle | `scan_angle`: the position from which every point is seen at its scan angle off nadir | about 1 m horizontally, 0.1 m vertically (RIEGL, DJI); 8 m on YellowScan |
| Single returns without scan angle, or no GPS time | refused | provide the trajectory |

The GeoJSON layer draws the path as lines, split where positions are more than 5 s apart (so no line crosses
unscanned turns), about one vertex per metre, in WGS84 when the cloud's CRS is known.
`3dtrees-trajectory-geojson file.csv` makes the same layer from an existing path CSV.

## Static scans are refused

On the largest block of the cloud (one file, point source and time block):

1. With multi-return pulses and GPS time: static when one position explains all pulse lines within 5 cm.
2. Otherwise: static when near-ground points pile up around one spot at least 15 times more than an even spread
   would (a fixed scanner's density falls with range). Measured on 3Dtrees data: TLS 22 to 128, ULS 1.8 to 7.7,
   MLS 4.3.

A static cloud ends with exit code 2 and a report that names the evidence. `--no-platform-check` skips the check.

Multi-file datasets are read together. Multi-return records that were stored apart from each other (some software
sorts points spatially) are joined by GPS time and point source after reading.

## Output and quality

- Columns `time,x,y,z,roll,pitch,yaw,name`, name last so readers of the first four columns are unaffected. Roll, pitch
  and yaw stay empty (positions only), the name is empty for a path, and speed is not stored (viewers derive it).
- The header carries the standard keys plus `kind: path`, `origin: estimated`, `method` and `quality`; the report adds
  `median_miss_m` (how far the pulse lines miss the positions, median) and `typical_error_m`.
- `quality` is `good` when the lines meet within 10 cm, `rough` for 0.1 to 0.5 m and for scan-angle paths, and
  `unreliable` when the lines miss by more than 0.5 m (broken return numbers).
- A time window whose pulse lines miss the fitted position by more than 10 times the median (at least 0.25 m) is
  dropped, and so is a position more than 1 m off the path of its neighbours (within 1.5 s).

## Speed and memory

LAZ is decoded in parallel (`--threads`, default 8; COPC is read serially). Coordinates are kept as float32 relative
to the first file's offset. At most 40 M multi-return and 10 M single-return pulses are kept (40 M single returns when
a cloud has no multi-return pulses), as every n-th pulse so the whole acquisition stays covered. `--full` keeps every
pulse.

| Dataset | Points | Before | Now |
| --- | --- | --- | --- |
| 1844 YellowScan, 3 flights | 678 M | 146 s, 13.5 GB | 119 s, 3.8 GB |
| 2019 RIEGL ULS | 93 M | 69 s, 6.9 GB | 45 s, 3.5 GB |
| 2018 MLS, single return | 2 B | 244 s, 19.7 GB | 213 s, 7.0 GB |

Ground-truth backfill (36 flights and walks with delivered paths, `/mnt/data/kg281/3dtrees_trajectories` on geos):
every run under 8 GB peak, most 20-60 s, the 70 GB DJI dataset 2585 12 min.

## Limits

- Positions exist only while the cloud has points. Many exports cut turns and transit legs, so those parts of a
  flight have no positions.
- Thinned or resampled clouds lose the pulse structure (re-sorted ones are joined again by time); the report's
  `multi_return_share` shows what was found.
- A pulse split across two read chunks (5 M points) is lost; a handful per chunk.
- COPC works but loses accuracy: it stores points by octree node, not in acquisition order. Returns are then joined by
  GPS time alone, which mixes pulses of multi-beam scanners (1807 YellowScan: 0.001 m median from LAZ, 0.075 m median
  and 5 m at the 95th percentile from COPC). Original per-flight LAS/LAZ files are the better input.
- Time stays in the cloud's own clock (GPS week seconds or adjusted GPS time, read from the LAS header). The CRS is
  the cloud's.

## Development

```bash
uv venv .venv && uv pip install -e ".[test]"
.venv/bin/python -m pytest
```

The tests build synthetic clouds with a known sensor (`tests/synthetic.py`); no real data is checked in.
