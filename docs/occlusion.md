# Occlusion

[Back to README](../README.md)


```bash
3dtrees-trajectory standardized/strip1.las standardized/strip2.laz -o trajectory.csv \
  --occlusion-copc occlusion.copc.laz --occlusion-cog occlusion.tif

3dtrees-trajectory standardized/strip1.las standardized/strip2.laz --trajectory delivered.csv \
  --dtm raster/dataset_dtm.tif --occlusion-copc occlusion.copc.laz --occlusion-cog occlusion.tif
```

Supply both occlusion destinations. When `--trajectory` is supplied, occlusion uses that path. With no
`-o` the tool checks the cloud and writes the occlusion products without fitting a path.
A delivered moving-path CSV starts with `time,x,y,z`, optionally preceded by `#` metadata. It needs at least two
positions, finite coordinates, and unique timestamps in the cloud's GPS clock and coordinate system.
Resolve week seconds versus adjusted GPS time before running. Both modes require a moving sensor.

The compiled walker uses the OccPy 0.3.0 grid traversal and linearly interpolates the sensor position for each pulse. Echoes up
to one second outside the path take the nearest endpoint. More distant echoes are counted and skipped.
The tool refuses disjoint time ranges. A voxel is filled (1) if there is a registered return, empty (2)
if seen clear, occluded (3) if only traversed behind returns, or unobserved (4) in the absence of traversal.

Voxel size defaults to 0.5 m. `--voxel` accepts 0.05–10 m.

### Ground

With `--dtm PATH`, the tool uses the terrain
raster cell containing each voxel column centre (nearest-neighbour sampling). The raster tool's 1 m
`<prefix>_dtm.tif` is a supported input. If both inputs declare a real CRS, they must agree. A DTM without a
CRS or with a local placeholder is sampled in the cloud coordinates. When one input has a real CRS and the
other has no usable CRS, the products and debug report record the coordinate assumption. A DTM outside
the cloud triggers a warning and falls back to lowest hits. Nodata and uncovered columns also use lowest hits.
That estimate starts at the centre of the lowest filled voxel per column, fills empty
columns from the nearest estimate and replaces heights more than 2 m above the local 5 m median.
The highest filled voxel caps each column.

| File | Contents |
| --- | --- |
| `occlusion.copc.laz` | One point at each filled, occluded or unobserved voxel centre from ground through canopy top. Extra dimensions are `occlusion`, `nhit`, `nmiss`, `nocc`, `occluded_fraction` and `height_above_ground`. The fraction is `nocc / (nhit + nmiss + nocc)`, with NaN for unobserved voxels. |
| `occlusion.tif` | North-up COG with one pixel per voxel column and one height band per voxel layer. The default bands are `occluded_0-0.5m`, `occluded_0.5-1m` and so on. Each value is the occluded voxel count divided by the sum of filled, empty and occluded voxels, with NaN when the sum is zero. The final `unobserved_column` band gives the unobserved share of all voxels from ground through canopy top. |

With step `s`, zero-based band `k` covers `[k·s, (k+1)·s)` metres above ground. The pixel size and height-band
step both follow `--voxel` and default to 0.5 m. The COG records the step in `OCCLUSION_BAND_M` and the debug report
records it in `occlusion.cog_band_step_m`.

The COG and COPC carry `OCCLUSION_GROUND` with the value `dtm`, `lowest_hit`, or `dtm+lowest_hit`. Mixed
runs include `OCCLUSION_GROUND_DTM_SHARE` and `OCCLUSION_GROUND_LOWEST_HIT_SHARE`, measured over all
columns in the XY grid. The debug report records `occlusion.ground` and, when a DTM is supplied,
`occlusion.ground_column_shares`.

If there are no qualifying voxels, the COPC is zero bytes and the COG is valid with NaN values. On an
occlusion error both destinations are zero bytes and stderr gives the reason. A combined job can
still return its estimated trajectory after a traversal error. An invalid terrain input exits 2.
The debug report includes traversal counts, in-range class totals and the pooled occluded share per height band.

The first read builds the trajectory sample, indexes record ranges by GPS time and collects their column extrema.
The path selects which echoes overlap the occlusion run. Blocks that intersect a path boundary are decoded to
apply the exact clock limits. Separated returns are joined by GPS time, point source and scanner channel within
each file.

Occlusion uses independent workers over disjoint GPS-time windows. Workers seek the original LAS/LAZ and
trace each complete pulse once. They use the same voxel lattice and store compact counts between the ground
and canopy. The planner splits a window further in time when its working memory or echo count exceeds the
budget. Set the planning ceiling with `--grid-memory-gib` (16 GiB by default). At runtime the scheduler reserves
the estimated memory of each window within the available cgroup memory and `--threads` limit. Clean file
cache is reclaimable. Dirty and writeback pages count toward the working budget. The scheduler adjusts
subsequent reservations using each completed window's peak RSS. If a worker runs out of memory, the job exits
with code 127 so the scheduler or operator can retry it with more memory.

Workers merge counts and classify independent voxel blocks. The writer publishes COPC records and raster windows
in grid order. Repeat the same command to resume an interrupted run and reuse finished windows. Keep the
`.occlusion-*` checkpoint beside the COPC until both outputs are published. A checkpoint write failure
from exhausted storage reports exit code 2.
A merged COPC can supply the trajectory estimate. For occlusion, supply the standardized LAS/LAZ collection.

The delivered path is used for every occlusion ray when supplied. `OCCLUSION_DELIVERED_MISS_FLAG` records whether
the median perpendicular distance from that path to usable pulse lines exceeds 0.1 m. Its values are `true`,
`false`, and `unavailable` when there are no usable lines. A CSV produced in the same run records the result
as `# delivered_path_miss_flag:`. OccPy traces single-return shots independently.

Both occlusion products record the path source (`OCCLUSION_TRAJECTORY`), path method
(`OCCLUSION_PATH_METHOD`), quality (`OCCLUSION_QUALITY`), any warnings (`OCCLUSION_WARNINGS`), and the
delivered-path miss flag when applicable. A scan-angle estimate is marked `rough`. The COG stores these as
dataset tags. The COPC stores the same keys in the `3dtrees` extended variable length record (record ID 1).

COPC indexing and COG conversion use temporary files beside the product destinations. The job needs space for
the uncompressed voxel LAS, intermediate TIFF and final outputs. Wide rasters use a band-interleaved
intermediate and up to four job slots for the COG conversion. They skip intermediate compression when scratch
has room for the uncompressed raster. The image includes Untwine 1.5.1. For local
runs, install the package with `python -m pip install .` using a C++17 compiler and set
`TRAJECTORY_UNTWINE` to an executable Untwine 1.5.1 installation.

Untwine randomizes COPC overview sampling. For regression checks, compare decoded records in voxel-centre order.
[Untwine 1.5.1 sampling implementation](https://github.com/hobuinc/untwine/blob/1.5.1/bu/Processor.cpp#L99-L126)

## Voxel size

`--voxel` sets both the voxel edge and the COG height-band step (default 0.5 m; CLI range 0.05–10 m, Galaxy form
0.25–2 m). Memory grows steeply below 0.5 m: dataset 1844 (7.7 GB) peaked at 6.3 GiB with 0.5 m and 35.9 GiB with
0.25 m, about (0.5/voxel)^2.5. See [resources.md](resources.md#voxel-size).
