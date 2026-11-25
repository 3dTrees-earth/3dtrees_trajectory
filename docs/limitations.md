# Limitations and plans

[Back to README](../README.md)

## Not supported in 1.0.0

- **TLS** (static scans) is refused. A file without GPS time exits 3 before the static check; a static scan with GPS
  time exits 2. The 3Dtrees platform does not start the tool for TLS. Station support is planned for 1.1.
- **Single returns without a scan angle** need a delivered path (exit 4 otherwise), unless the file holds shots with
  first and last returns stored apart (see [path-estimation.md](path-estimation.md#pulse-pairing)).
- **Real Earth-centred coordinates** (EPSG:4978 far from the origin) are refused; reproject the cloud first.
- **Local against EPSG** delivered paths cannot be aligned; occlusion is refused for them.
- **Rounded or offset delivered paths** (e.g. DJI POS rounded to 0.001°, lever-arm offsets, clock shifts) are used as
  given; the miss flag shows the problem but the tool does not correct it.
- **Turns and transits outside the cloud** have no estimated positions.
- **Small voxels** need much more memory (about (0.5/voxel)^2.5); 0.25 m on very large flights can exceed the
  Galaxy allocation and exit 127.
- GPS week rolls and independent flights that share the same clock and point source should be split before
  processing.

## Planned for 1.1

- TLS scan stations (estimated or delivered) for occlusion.
- A DSM input as the upper bound of each column, and rays that stop at the terrain.
- `lower_threshold`: ignore voxels below a height above ground.
- An area-of-interest crop of the voxel grid.
- A Livox keyframe adapter for SLAM paths.
