# 3dtrees_trajectory

Estimates where the lidar sensor was from the point cloud alone: the flight or walking path for ULS, MLS and ALS, or
the scanner positions for TLS. The output is the 3Dtrees trajectory CSV (`3dtrees-trajectory/1`), the same format
the platform's trajectory converter writes, plus a JSON quality report.

```bash
3dtrees-trajectory flight.laz -o flight_trajectory.csv
3dtrees-trajectory scan01.laz scan02.laz scan03.laz -o scan_positions.csv
```

## How it works

| Cloud | Method | Accuracy on 3Dtrees data |
| --- | --- | --- |
| Moving sensor, multi-return pulses, GPS time | `lines`: a pulse's first and last return lie on a line through the sensor; per 0.5 s window, the position and velocity closest to all lines | under 0.5 m on 31 of 34 delivered trajectories, most YellowScan flights about 0.01 m |
| Moving sensor, single returns, scan angle | `scan_angle`: the position from which every point is seen at its scan angle off nadir | about 1 m horizontally, 0.1 m vertically (RIEGL, DJI); 8 m on YellowScan |
| Static scanner, multi-return pulses, GPS time | `lines` per scan | about 1 cm; 21 of 22 TERN plots find every listed station |
| Static scanner, single returns | `density`: the peak of near-ground point density per scan | 0.5 to 0.7 m (FARO, BLK360); z is the ground under the scanner |
| Single returns without scan angle (moving), or no way to tell scans apart (static) | refused or flagged | provide the trajectory or one file per scan |

Next to the CSV and the report, the tool writes a GeoJSON map layer: the path as lines (split where positions are
more than 5 s apart, so no line crosses unscanned turns, about one vertex per metre) or the TLS stations as points,
in WGS84 when the cloud's CRS is known. `3dtrees-trajectory-geojson file.csv` makes the same layer from an existing
trajectory CSV.

## Platform and scans, detected

No platform has to be given. On the largest scan group of the cloud:

1. With multi-return pulses and GPS time: static when one position explains all pulse lines within 5 cm.
2. Otherwise: static when near-ground points pile up around one spot at least 15 times more than an even spread
   would (a fixed scanner's density falls with range). Measured on 3Dtrees data: TLS 22 to 128, ULS 1.8 to 7.7,
   MLS 4.3.

`--platform mobile|static` overrides this; the report names the evidence.

Multi-file datasets are read together. Scans are told apart automatically:

- With GPS time: by file, point source id and time blocks. A block ends where the file-order time jumps by more than
  30 s (scans concatenated into one cloud); if the points are stored out of acquisition order, blocks are gaps in
  sorted time instead. Stations found twice (one scan cut into several blocks or tiles) are merged within 5 cm.
- Without GPS time: by point source id when there are several (a project exported as spatial tiles keeps one source
  per scan), else by file, else one scan.

`--split file|source|time` overrides this. Groups under 2,000 pulses are left out as fragments.

Multi-return records that were stored apart from each other (some software sorts points spatially) are joined by GPS
time and point source after reading.

## Quality checks

- Moving sensor: a time window whose pulse lines miss the fitted position by more than 10 times the median (at least
  0.25 m) is dropped, and so is a position more than 1 m off the path of its neighbours (within 1.5 s).
- The CSV header carries `method`, `quality`, `typical_error_m` and `median_miss_m`, and each row a `miss_m` column
  after z (how far the pulse lines miss that position). The report's `quality` is `good` when the lines meet within 10 cm, `rough` for scan-angle and density positions,
  and `unreliable` when the lines miss by more than 0.5 m (broken return numbers).

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
| 3582 TERN TLS, merged | 66 M | 136 s, 8.7 GB | 37 s, 2.1 GB |
| 1887 FARO TLS, 30 tiles | 330 M | 160 s, 14.1 GB | 145 s, 4.3 GB |

Ground-truth backfill (65 datasets, `/mnt/data/kg281/3dtrees_trajectories` on geos): every run under 8 GB peak, most
20-60 s, the 70 GB DJI dataset 2585 12 min.

## Limits

- Positions exist only while the cloud has points. Many exports cut turns and transit legs, so those parts of a
  flight have no positions.
- Thinned or resampled clouds lose the pulse structure (re-sorted ones are joined again by time); the report's
  `multi_return_share` shows what was found.
- A pulse split across two read chunks (5 M points) is lost; a handful per chunk.
- COPC works but loses accuracy: it stores points by octree node, not in acquisition order. Returns are then joined by
  GPS time alone, which mixes pulses of multi-beam scanners (1807 YellowScan: 0.001 m median from LAZ, 0.075 m median
  and 5 m at the 95th percentile from COPC), and merged TLS scans whose clocks overlap can no longer be told apart
  (3582 is then read as a flight). Original per-flight or per-scan LAS/LAZ files are the better input.
- Time stays in the cloud's own clock (GPS week seconds or adjusted GPS time, read from the LAS header). The CRS is
  the cloud's.

## Development

```bash
uv venv .venv && uv pip install -e ".[test]"
.venv/bin/python -m pytest
```

The tests build synthetic clouds with a known sensor (`tests/synthetic.py`); no real data is checked in.
