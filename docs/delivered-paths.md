# Delivered paths

[Back to README](../README.md)

A delivered path replaces the estimate as the sensor position for occlusion. The tool still estimates its own path
when `-o` is given (the estimate is optional: it may fail and occlusion still runs).

## Alignment of each file

Every `--trajectory` file is aligned on its own before anything else:

1. **CRS.** The path's `# crs:` is compared with the cloud's:
   - same EPSG, or both local: used as is;
   - different EPSG: converted to the cloud's CRS with pyproj (3D when both declare a vertical EPSG, else
     horizontal only, keeping z); warning *"The delivered path's CRS (X) was converted to the cloud's CRS (Y)."*;
   - local against EPSG (either way), or a failed conversion: occlusion is refused with the supplementary note;
   - no `# crs:` line (a non-platform file): the cloud's CRS is assumed, with a warning.
2. **Clock.** If the path's times do not overlap the cloud's GPS time, its declared `# time_system` is converted
   (`unix_utc` with leap seconds, `gps_unix`, `gps_adjusted`, or `gps_sow` with `gps_week`) to the cloud's week
   seconds or adjusted GPS time. With one path the test is the share of the cloud's span covered (≥ 50 %); with
   several paths, the share of the path's own span inside the cloud's span, so a path for one flight converts too.
   Recorded as `clock_converted` in the report and in the warnings.
3. **Duplicates.** Files with the same times and positions within 1 cm are used once (warning).

## Several paths per dataset

The tool handles one merged file with several flights, several files with one path, and several files with several
paths in the same way. The matching unit is a **pulse segment**: pulses of one input file and point source in a
contiguous GPS-time stretch, cut at pauses over 5 s and wherever any path starts or ends (each path piece widened by
1 s). For each segment:

| Paths covering ≥ 50 % of the segment | Choice |
| --- | --- |
| none | its pulses are not traced; warning *"N % of pulses had no delivered path covering their time"* |
| one | that path |
| several (duplicate export, a second platform or sensor head at the same time, the same week seconds on another day) | the path with the lowest median miss to a fixed, evenly spaced sample of up to 20,000 of the segment's multi-return pulse lines; without lines, the first path given |

Occlusion is refused when fewer than 50 % of the pulses get a path ("The point cloud and path have no overlapping GPS
times." when none does). Each pulse is then traced from its own segment's path; there is still one grid, one DTM, one
COPC and one COG per dataset.

**One path covering the whole cloud** is traced exactly as before multi-path support (same products byte for byte).

## Outputs

- CSV header: `# delivered_path_miss_flag:` (`true` when the median miss over all assigned pulses exceeds 0.1 m),
  `# delivered_path_files:`, `# delivered_path_coverage:` (share of pulses with a path), one `# warning:` per
  conversion, duplicate, uncovered share or ambiguous choice.
- COG/COPC tags: `OCCLUSION_TRAJECTORY=delivered`, `OCCLUSION_DELIVERED_MISS_FLAG`,
  `OCCLUSION_DELIVERED_MISS_MEDIAN_M`, `OCCLUSION_DELIVERED_PATH_FILES`, `OCCLUSION_DELIVERED_PATH_COVERAGE`,
  `OCCLUSION_WARNINGS`.
- Debug report `occlusion.delivered_path_check`: per file (CRS and clock conversions, duplicate of), coverage, miss
  statistics and every segment with file, point source, GPS span, candidates, chosen path and per-candidate miss.

## Examples from production data

| Dataset | Layout | Result |
| --- | --- | --- |
| 1804 | 1 cloud file, 12 paths; two heads flying the same 766 s ~10 m apart | point sources 110xx → path 19, 120xx → path 20 (6 mm vs 9.8 m); coverage 100 %, miss 6 mm |
| 1844 | 1 cloud file, 3 flights | 46 segments, coverage 100 %, miss 2 mm |
| 2019 | 2 cloud files, 3 paths | coverage 100 % |
| 2018 | 10 files, 10 Unix-time SLAM paths, single returns | coverage 100 % |
| 2585 | 25 files, 24 Unix-time DJI POS files, parallel drones | 57 overlapping segments decided (33–39 m vs 170–450 m), coverage 95 %; miss 35 m because the POS positions are rounded to 0.001° |
