# Outputs

[Back to README](../README.md)

The tool creates output directories, publishes each file by rename and writes zero-byte files at safe destinations
when it refuses.

## Trajectory CSV

```text
# format: 3dtrees-trajectory/estimated
# kind: path
# origin: estimated
# crs: EPSG:32632            (or local)
# vertical_crs: unknown      (EPSG of a compound cloud CRS when known)
# time_system: gps_sow       (or gps_adjusted)
# gps_week: unknown
# attitude: none
# source_layout: estimated_lines
# source_file: 1807.las
# converter: 3dtrees_trajectory 1.0.0
# method: lines              (lines, scan_angle, lines_scan_angle)
# quality: good              (good, rough, unreliable)
# warning: …                 (zero or more)
# delivered_path_miss_flag: false     (with --trajectory)
# delivered_path_files: 2             (with --trajectory)
# delivered_path_coverage: 1.0        (with --trajectory)
time,x,y,z,roll,pitch,yaw,name
381583.278017,596279.5609,5675269.9795,447.1477,,,,
```

One position per `--interval` (0.5 s) where the cloud has points; turns and transit legs outside the cloud have no
positions. Time and XYZ use the cloud's clock and coordinates.

## Occlusion COPC

One point per filled, occluded or unobserved voxel centre from ground to canopy top. Extra dimensions: `occlusion`
(1 filled, 2 empty, 3 occluded, 4 unobserved), `nhit`, `nmiss`, `nocc`, `occluded_fraction`
(`nocc / (nhit + nmiss + nocc)`, NaN when unobserved) and `height_above_ground`. Tags are in the `3dtrees` EVLR
(record ID 1).

## Occlusion COG

North-up, one pixel per voxel column, one band per voxel layer above ground (`occluded_0-0.5m`, `occluded_0.5-1m`,
…) and a final `unobserved_column` band. Band `k` covers `[k·s, (k+1)·s)` m above ground with step `s` = `--voxel`.

| Tag | Meaning |
| --- | --- |
| `OCCLUSION_BAND_M` | height-band step |
| `OCCLUSION_GROUND` | `dtm`, `lowest_hit` or `dtm+lowest_hit` (+ `_DTM_SHARE`, `_LOWEST_HIT_SHARE`) |
| `OCCLUSION_TRAJECTORY` | `estimated` or `delivered` |
| `OCCLUSION_PATH_METHOD`, `OCCLUSION_QUALITY` | method and quality of the path used |
| `OCCLUSION_WARNINGS` | JSON list of warnings (CRS/clock conversions, DTM assumptions, coverage) |
| `OCCLUSION_DELIVERED_*` | miss flag, median miss, path files, coverage (delivered paths) |

## Debug report (`--debug`)

`<csv>.report.json` (or `--report PATH`): tool version, files, points, CRS, time system, method, quality, median
miss, warnings, `pulse_check` (pairing counts and per-file `return_modes` with `mode`, `has_scan_angle`,
`unpaired_multi_echoes`, `line_pulses`), `occlusion` (class totals, occluded share by band, ground source, delivered
path check with segments, time windows) and `resources` (wall, CPU, peak RSS, per-stage measurements). On refusal:
`error`, `exit_code`, `pairing`. Publish it only when debugging.
