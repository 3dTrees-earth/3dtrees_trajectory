# Changelog

## 1.0.0

Re-cut of 1.0.0 (24 Nov 2026):

- Treat a geocentric CRS on a cloud centred within 1,000 km of the origin as a local placeholder. Reject Earth-centred inputs at or beyond that distance with exit 2.
- Accept terrain rasters without a CRS in the cloud coordinates. Warn when assuming a CRS for one input or falling back to lowest hits because the raster is outside the cloud.
- Return exit 127 for memory exhaustion and worker kills, with a debug refusal report and empty products at safe destinations.

First cut:

- Estimate a moving lidar trajectory from one standardized ULS or MLS collection of scanner-order LAS/LAZ files. The inputs need finite, varying GPS time in a common time system and metric coordinates. A time-aligned delivered path can drive occlusion directly. An optional terrain GeoTIFF must share the point cloud's CRS.
- Write a 3Dtrees trajectory CSV and GeoJSON, an occlusion voxel COPC and an occlusion height-band COG. Use `--debug` for a JSON report with quality, refusal reasons, and resource measurements.
- In the 3Dtrees workflow, run after Raster for ULS and MLS when each standardized file has varying GPS time and can use multiple returns, a nonzero single-return scan angle, or a standardized delivered path with timestamps. The CLI validates the clock and the selected method at runtime.
- Exit 0 indicates success. Exit 2 covers invalid inputs, static scans, and hard resource or output failures. Exit 3 means GPS time cannot support a moving path. Exit 4 means the tool could not estimate a usable path. Refusals print their reason and leave requested products empty at safe paths.
- Occlusion requires scanner-order files and uses 0.5 m voxels by default. A merged input COPC supports path estimation only. TLS station support is planned for 1.1.
- Compared with the round-10 candidate, the product version string changes to `1.0.0` in the CSV converter header and debug JSON. Product schemas and scientific results are the same.
