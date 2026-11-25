# Inputs

[Back to README](../README.md)

## Point clouds

Pass all standardized LAS/LAZ files of one acquisition to one invocation, in consistent collection order. A
collection can hold 1 to 65,536 files. The parallel chunked reader keeps the scanner's record order within each
file, including time blocks, so keep pulse siblings and their metadata together during standardization.

- LAS/LAZ 1.2–1.4, point formats 0–10. Formats 0 and 2 have no GPS time and are refused (exit 3).
- GPS time must be finite and vary across the points (constant or missing time: exit 3). All files share one GPS
  time system (week seconds or adjusted standard GPS time).
- Coordinates in metres. Keep **ReturnNumber**, **NumberOfReturns**, **ScanAngle** (or ScanAngleRank in legacy
  formats), **PointSourceID** and the **scanner channel**: the path estimate pairs echoes of one shot by GPS time,
  point source and channel.
- A merged COPC can supply the path estimate, but occlusion needs the scanner-order LAS/LAZ files.
- Static scans (TLS) are refused (exit 2); TLS station support is planned for 1.1.

## Coordinate reference system

| Cloud CRS | Treatment |
| --- | --- |
| Projected or geographic EPSG in metres | used as is; GeoJSON in WGS84; products carry the CRS |
| Compound CRS | horizontal EPSG in `# crs:`, the vertical EPSG in `# vertical_crs:` when known |
| none | **local**: `# crs: local`, GeoJSON in cloud coordinates with `crs: local`, products without CRS |
| EPSG:4978 with the cloud centred < 1,000 km from the origin | Standardization's placeholder: **local**, as in Raster |
| EPSG:4978 further out (real ECEF) | refused (exit 2): reproject to a projected CRS in metres |

Older files that say `unknown` are read as `local`.

## Terrain raster (`--dtm`)

One GeoTIFF for the whole dataset, typically the 3Dtrees Raster `<prefix>_dtm.tif` (1 m). Each voxel column uses the
raster cell containing its centre (nearest neighbour). If both the cloud and the raster declare a real CRS, they must
agree (exit 2 otherwise). A raster without a CRS, or with the local placeholder, is sampled in the cloud's
coordinates; when only one side has a real CRS the assumption is recorded in the products and the debug report.
Nodata and uncovered columns use the lowest-hit ground estimate; a raster outside the cloud gives a warning and the
lowest-hit ground everywhere. See [occlusion.md](occlusion.md#ground).

## Delivered paths (`--trajectory`)

One or more 3Dtrees trajectory CSVs as written by the platform's converter: `#` header lines, then a header row
starting `time,x,y,z` (attitude and `name` columns may follow). Each file needs at least two positions with finite
coordinates and strictly increasing, unique GPS times.

```text
# format: 3dtrees-trajectory/1
# crs: EPSG:32632
# vertical_crs: unknown
# time_system: gps_sow
# gps_week: 2332
time,x,y,z,roll,pitch,yaw
381469.001000,596355.5460,5675393.3920,384.8050,13.187000,-8.243000,113.322000
```

The header keys `crs`, `vertical_crs`, `time_system` and `gps_week` let the tool convert the path to the cloud's CRS
and clock. Formats: `3dtrees-trajectory/1` (uploaded paths), `/2` (rows with names, e.g. scan positions),
`/estimated` (written by this tool). How several paths are matched to the cloud: [delivered-paths.md](delivered-paths.md).
