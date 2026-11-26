# Validation

[Back to README](../README.md)

## Estimated path against delivered paths

Matched by GPS time; the delivered path is interpolated at each estimated time (re-cut 1.0.0).

| Dataset | Method | Median h / v | p95 h | Within 0.1 m | Note |
| --- | --- | --- | --- | --- | --- |
| 1807 YellowScan (PPK) | lines | 1 mm / 1 mm | 6 mm | 99.7 % | max 0.13 m |
| 1844 | lines | 2 mm / 1 mm | 7 mm | 100 % | |
| 2014 | lines | 6.4 / 5.5 cm | 7.4 cm | 97 % | constant offset (lever arm/datum); 1.5 cm after removing it |
| 3773 | lines | 5.4 / 2.5 cm | 6.8 cm | 100 % | constant offset; 3.3 cm after removing it |
| 1925 RIEGL | lines | 24.8 / 1 cm | 25 cm | 0 % | +85 ms clock shift in the delivered path; 1.1 cm after correcting it |
| 2587 LiDAR360 | lines (returns stored apart) | 1.7 mm / 1.4 mm | 7 mm | | first + last returns only |
| 2014 single-return copy | scan_angle | 1.04 / 0.13 m | 1.57 m | 0 % | |
| 3773 single-return copy | scan_angle | 1.23 / 0.04 m | 1.77 m | 0 % | |
| 1925 single-return copy | scan_angle | 1.24 / 0.06 m | 15.4 m | 0.3 % | |
| 1807 single-return copy | scan_angle | 8.4 / 0.46 m | 13.1 m | 0 % | YellowScan does not record the angle off nadir |

## Multi-path assignment

Delivered-path check and segment choice on all paths of each dataset: 1791, 1818, 1844, 1926, 1930, 2018, 2019 and 1804 (two heads) at 100 % coverage; 2585 at 95 % with the miss flag set by the 0.001° rounding of its DJI POS. Details: [delivered-paths.md](delivered-paths.md#examples-from-production-data).

## Test matrix (local, 21 cases)

All cases run the CLI as the Galaxy wrapper does (estimate + occlusion, DTM when present, all delivered paths) and the
3Dtrees platform gate on the same files. All behave as expected:

| Cell | Dataset | Outcome |
| --- | --- | --- |
| ULS multi-return, no path / 1 path / 1 path + Raster DTM | 1921, 1807 | exit 0, lines good, ground `dtm+lowest_hit` |
| 1 file + 2 paths; 2 files + 3 paths; 2 files + 1 path | 1926, 2019, 1921 split | exit 0, coverage 1.0 |
| overlapping paths (two heads) | 1804 | exit 0, 12 paths, miss 6 mm |
| WGS84 paths (+ DTM), week seconds | 2014, 3773 | exit 0, converted, miss 5–7 cm |
| cloud without CRS + EPSG path | 1921 no-CRS | exit 0, `crs: local`, occlusion refused (gate drops the path) |
| returns stored apart, intermediate returns missing | 2587 | exit 0, lines good |
| single + scan angle; single without angle (no path / path) | 1925, 1921 | exit 0 rough; exit 4; exit 0 |
| no GPS time; constant GPS time | 1921 PF2, 3516 | exit 3 |
| MLS: fake EPSG:4978 + DTM without CRS; MLS + DTM; SLAM 10 files / 10 Unix paths | 1788, 2018 | exit 0 |
| TLS with stations; TLS two scans | 3950, Millyford | refused (exit 3: no GPS time; the gate refuses TLS first) |

## Product equality

- Re-cuts keep 1807 and 2014 COG and decoded COPC records identical; the CSV differs only in its format line.
- Products are identical at 2–10 slots.
- Untwine randomizes COPC overview sampling: compare decoded records in voxel-centre order.

## Not covered yet

The local matrix has no real dataset for: multi-return MLS with delivered paths, a standardized multi-file MLS
collection, a compound CRS, real ECEF coordinates, a partially covering Raster DTM on a standardized cloud, local
paths on a local cloud (synthetic tests only), TLS with GPS time, and a real (not derived) single-return ULS with
scan angle.
