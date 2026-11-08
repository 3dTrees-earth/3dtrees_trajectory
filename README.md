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

Static scans are told apart by point source id, then by input file, then by jumps in GPS time (merged clouds that
were concatenated scan by scan). `--split` overrides this.

The platform is detected: when one position explains the multi-return lines of the first scan within 5 cm, the
cloud is static. `--platform mobile|static` overrides this.

## Limits

- Positions exist only while the cloud has points. Many exports cut turns and transit legs, so those parts of a
  flight have no positions.
- The cloud must keep its points in their original order. Re-sorted, thinned or merged-and-resampled clouds lose the
  pulse structure; the report's `multi_return_share` shows what was found.
- A pulse split across two read chunks (5 M points) is lost; a handful per chunk.
- Time stays in the cloud's own clock (GPS week seconds or adjusted GPS time, read from the LAS header). The CRS is
  the cloud's.

## Development

```bash
uv venv .venv && uv pip install -e ".[test]"
.venv/bin/python -m pytest
```

The tests build synthetic clouds with a known sensor (`tests/synthetic.py`); no real data is checked in.
