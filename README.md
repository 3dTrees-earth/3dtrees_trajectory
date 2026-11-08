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

## Limits

- Positions exist only while the cloud has points. Many exports cut turns and transit legs, so those parts of a
  flight have no positions.
- Thinned or resampled clouds lose the pulse structure (re-sorted ones are joined again by time); the report's
  `multi_return_share` shows what was found.
- A pulse split across two read chunks (5 M points) is lost; a handful per chunk.
- Time stays in the cloud's own clock (GPS week seconds or adjusted GPS time, read from the LAS header). The CRS is
  the cloud's.

## Development

```bash
uv venv .venv && uv pip install -e ".[test]"
.venv/bin/python -m pytest
```

The tests build synthetic clouds with a known sensor (`tests/synthetic.py`); no real data is checked in.
