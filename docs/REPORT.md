# 3dtrees_trajectory: what the tool does, when, and where the time goes

Status 20 Nov 2026, `main` after round 5. Scope of 1.0.0: moving lidar (ULS and MLS). TLS is v1.1.

## 1. What the tool does

The tool takes the point cloud of one ULS or MLS dataset and does two things.

1. **Trajectory.** It estimates where the sensor was, every 0.5 s, from the points alone.
2. **Occlusion.** It traces every laser pulse from the sensor through its returns. It then marks each 0.5 m voxel as
   seen, empty, hidden behind vegetation, or never reached.

**Inputs**

| Input | Required | Notes |
| --- | --- | --- |
| One collection of standardized LAS/LAZ files per dataset | yes | The output of the standardization tool. One job takes the whole collection. Files are read in parallel chunks, in scanner order. |
| Delivered trajectory (3Dtrees trajectory CSV) | no | If given, it is always used for occlusion instead of the estimate. |
| DTM raster | no | Ground for heights above ground. Default: the lowest hits per column. |

A merged COPC is accepted for trajectory estimation. Occlusion requires the scanner-order LAS/LAZ collection.

**Outputs**

| Output | Content |
| --- | --- |
| Trajectory CSV (`3dtrees-trajectory/2`) | time, x, y, z per position, with method and quality in the header. Plus a GeoJSON line for the map. |
| JSON report | Method, quality, warnings, refusal reason, exit code, and time, CPU, RAM and scratch for every stage. |
| Occlusion COPC | One point per voxel. Class 1 filled, 2 empty, 3 occluded, 4 unobserved, plus hit, pass and hidden counts. |
| Occlusion COG | Occluded share per 1 m height band above ground, 0.5 m pixels. |

## 2. What happens when

![Decision flow](figures/flow.png)

The steps below follow `src/threedtrees_trajectory/cli.py:main`.

1. **Argument checks** (`cli.py`). Outputs must not overwrite inputs. The voxel size must be 0.05 to 10 m.
   At least one file is needed. A delivered trajectory is only accepted together with occlusion outputs. Any failure
   gives **exit 2**.
2. **Read and time-sort** (`cli.py`, `pulses.py`). Large files use independently decoded, ordered record ranges. GPS time must
   exist, be finite and vary (`pulses.py`). If not: **exit 3** (`coverage.py`). Returns are
   paired into pulses per file. With no usable pulse: **exit 4** (`cli.py`).
3. **Platform check** (`cli.py`, `estimate.py`). A static scanner is recognised when one position
   explains all pulse lines within 5 cm, or when ground points pile up 15 times or more around one spot. TLS is
   refused with **exit 2**. The TLS code is on the branch `feat/tls-stations-v1.1`.
4. **Path estimation** (`cli.py`, `estimate.py`, `solve.py`). The method is chosen per file, with independent fits started concurrently:
   - multi-return files: the pulse lines (first to last return) meet at the sensor. Method `lines`.
   - single-return files with a scan angle: cone fit on the scan angle. Method `scan_angle`.
   - mixed: lines, with gaps filled from scan angles per file. Method `lines_scan_angle`.
   - single returns **without** scan angle: no estimate is possible. **Exit 4** (`cli.py`), unless a
     delivered path is given.
5. **Quality label** (`output.py`). The median distance between pulse lines and the fitted path decides:
   - **good**: ≤ 0.1 m.
   - **rough**: 0.1 to 0.5 m, or any scan-angle position.
   - **unreliable**: > 0.5 m. The tool then retries with scan angles. Without a scan angle: **exit 4**
     (`cli.py`).
6. **Delivered path** (`cli.py`, `occlusion.py`). It always replaces the estimate for occlusion.
   The tool measures how far the pulse lines miss it. A median above 0.1 m is flagged in the report, but the path
   is still used. If the estimate fails but a delivered path exists, the job continues with the delivered path
   (`cli.py`).
7. **Occlusion grid** (`window_input.py`, `occlusion.py`). During the first read, indexed record ranges store time
   bounds and column extrema. Path coverage selects the summaries used to determine ground and canopy.
   The final vertical grid spans ground minus 1 m through canopy top plus 1 m, within the point extent.
8. **Traversal** (`time_windows.py`, `ray_walk.cpp`). Workers seek the relevant LAS/LAZ ranges, rebuild whole
   pulses and trace them in disjoint GPS windows. The window planner uses the available slots and splits
   further in time if needed. Finished count files can be reused after an interruption.
9. **Classes and products** (`occlusion.py`). Sum hit, pass and hidden counts blockwise. A hit makes a filled
   voxel, a pass makes an empty voxel, and hidden-only makes an occluded voxel. Keep an unobserved class
   for voxels with zero coverage. Then write COPC and COG.
10. **Occlusion failure** (`cli.py`). Report the error, publish zero-byte occlusion files and retain any
    successfully estimated trajectory.

### Exit codes

| Code | Meaning | Where |
| --- | --- | --- |
| 0 | Done | `cli.py` |
| 2 | Invalid arguments, unreadable input, out of memory, TLS refused, report not writable | `cli.py` |
| 3 | GPS time missing, constant or not finite | `coverage.py`, `pulses.py`, `estimate.py` |
| 4 | No path can be estimated (no pulses, single returns without scan angle, too few reliable positions) | `coverage.py`, `estimate.py`, `cli.py` |

Refusals include a reason in the JSON report. Exit 3, exit 4 and errors using a delivered path also ask the uploader to retain the path as supplementary data.

## 3. Stage resources

![Stage resources](figures/stages.png)

The local replays use ten CPU slots and a 140 GiB cap. Datasets 1807, 1844 and 2014 use local scanner-order LAZ
replays. Dataset 1788 uses its standardized MLS file. The before and after runs for each dataset use the same
source records. The stage CSV includes the spatial-spool reference and the historical 2585 measurements.

<!-- stages:start -->
| Stage | 1807 ULS | 1844 ULS | 2014 ULS | 1788 MLS |
| --- | ---: | ---: | ---: | ---: |
| Read, reconstruct pulses and index source | 4.7 s / 6.5 cores | 38.5 s / 8.2 cores | 7.6 s / 5.1 cores | 4.6 s / 3.3 cores |
| Path fitting | 1.0 s / 5.8 cores | 5.2 s / 3.4 cores | 1.0 s / 6.1 cores | 96.9 s / 2.4 cores |
| Extent and crop | 0.2 s / 1.0 cores | 0.8 s / 1.0 cores | 1.6 s / 1.0 cores | 4.8 s / 1.0 cores |
| Plan GPS windows | 0.0 s / 1.0 cores | 0.2 s / 1.0 cores | 0.1 s / 1.0 cores | 0.1 s / 1.0 cores |
| Rebuild and trace pulses | 10.7 s / 7.3 cores | 139.2 s / 7.5 cores | 24.6 s / 8.4 cores | 43.0 s / 8.0 cores |
| Merge, classify and write intermediates | 1.1 s / 1.0 cores | 9.0 s / 1.0 cores | 9.1 s / 1.0 cores | 39.3 s / 1.0 cores |
| Write voxel COPC | 3.5 s / 1.7 cores | 9.6 s / 4.0 cores | 7.7 s / 3.5 cores | 13.0 s / 5.9 cores |
| Write height-band COG | 0.5 s / 1.0 cores | 4.0 s / 2.2 cores | 3.8 s / 2.3 cores | 18.0 s / 1.6 cores |
<!-- stages:end -->

CPU use is CPU seconds divided by wall seconds. Traversal includes source seeks, decoding, pulse grouping,
the ray walk and checkpoint writes. The first read reduces per-column minima and maxima. Once the path is
known, preparation decodes only indexed ranges that straddle the accepted time limits. The measurements
separate preparation from the streaming traversal's own source reads. Extent preparation decoded one
200,000-record boundary range on 1807 and zero records on 1844, 2014 and 1788.

## 4. Streaming evaluation and speed-up

<!-- comparison:start -->
| Dataset | Spatial spool wall / CPU | GPS windows wall / CPU | Peak RSS, GiB | Peak scratch, GiB | Peak RSS per window, GiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1807 | 103 s / 1.2 cores | 22 s / 5.5 cores | 1.57 → 3.83 | 2.17 → 0.49 | 0.12 |
| 1844 | 1061 s / 1.7 cores | 208 s / 6.9 cores | 8.69 → 5.84 | 32.49 → 3.25 | 0.31 |
| 2014 | 145 s / 2.8 cores | 57 s / 5.3 cores | 16.81 → 4.66 | 6.48 → 2.36 | 0.49 |
| 1788 | 597 s / 6.6 cores | 223 s / 3.3 cores | 45.66 → 8.87 | 28.96 → 6.32 | 0.48 |
<!-- comparison:end -->

The RSS column shows the maximum of the largest process and the sampled sum of active worker RSS. Scratch
uses samples of allocated blocks. Spatial-spool scratch for the raw replays comes from archived matching-input
runs. File cache is accounted for separately by the capped job scope.

### Product equality

| Dataset | Reference | Filled / empty / occluded / unobserved | COG values differing | COPC voxel centres differing | COPC class or count fields differing |
| --- | --- | --- | ---: | ---: | ---: |
| 1807 | Round 4, one whole grid | 1,445,111 / 2,526,490 / 212,799 / 22 | 0 | 0 | 0 |
| 2014 | Round 4, one whole grid, 64 GiB planning budget | 6,951,441 / 11,182,745 / 1,001,850 / 24,589 | 0 | 0 | 0 |

The trajectory CSV, GeoJSON and COG bytes agree for both controls. The complete decoded COPC point
records also agree byte for byte when sorted by voxel centre. The compressed files differ after Untwine
indexing.

| Dataset | Round-4 spool reference | Filled / empty / occluded / unobserved change | COG cells differing | Max COG difference | Old-only / new-only COPC centres |
| --- | --- | --- | ---: | ---: | ---: |
| 1807 | 1 grid | 0 / 0 / 0 / 0 | 0 | 0 | 0 / 0 |
| 1844 | 2 spatial tiles | −1 / +1 / 0 / 0 | 0 | 0 | 1 / 0 |
| 2014 | 8 spatial tiles | 0 / 0 / 0 / 0 | 11 | 1.0 | 6 / 6 |

The prior tile cores were traced with local grid boundaries and buffered inputs. The full-grid controls
isolate equivalence from that clipping. The 1844 difference is one previously filled centre at
(366537.75, 5977455.25, 9.25). The 2014 COG differences against spatial tiles disappear when the
reference uses the same round-4 whole grid.

### GPS windows on a wide strip

The initial rectangular window grids exceeded the budget on 2014 and 1788. The compiled walker retains the
whole-grid entry and lattice calculations. Each worker stores counts only for the ground-to-canopy layers
in its indexed reach. The 1788 grid spans 2,000 × 2,000 × 1,278 voxels. Its widest planned windows each
reserved 67.1 million count cells and finished within the 16 GiB per-window budget. The largest
scanner-order ULS strip checked in local staging was 2014 at 658 × 740 m. Its per-window process
RSS and the 1788 measurements are in the comparison table. The scheduler bisects a window when
its estimated bytes or echo count exceed its budget.

Use direct GPS-window streaming for scanner-order LAS/LAZ. Pulse grouping runs inside each worker for all
job sizes. The output merge adds counts blockwise before assigning classes. Interleaved GPS times can make
several workers decode the same compressed record range. Half-open time bands assign each retained echo to
one window. The 1788 traversal decoded 41.1 million records for 24.7 million echoes traced.

### Read parallelism

<!-- read:start -->
| Dataset | Spatial-spool first read | GPS-window first read | Record-range workers |
| --- | ---: | ---: | ---: |
| 1807 | 9.6 s / 2.3 cores | 4.7 s / 6.5 cores | 10 |
| 1844 | 102.3 s / 2.5 cores | 38.5 s / 8.2 cores | 10 |
| 2014 | 11.1 s / 2.5 cores | 7.6 s / 5.1 cores | 10 |
| 1788 | 4.6 s / 2.2 cores | 4.6 s / 3.3 cores | 5 |
<!-- read:end -->

The old path used a native parallel LAZ decoder inside a sequential Python batch pipeline. The 2014
control spent 3.27 s decoding at 6.1 cores, 3.50 s reconstructing sampled pulses and 1.65 s writing
the source echo spool. Each batch waited for those operations before the next decode. Independent
record-range workers now decode, reconstruct and index their own ranges. The reader joins their pulse
columns in scanner order and carries the equal-time tail across range boundaries. The read figures
include worker startup, result collection and the ordered merge. The averages reached 8.2 cores
on 1844 and 5.1 on 2014 with ten workers.

The local 2585 COPC remains a historical control. Its previous 7.0 billion point spool run routed 28.3
billion echoes and peaked at 568 GiB scratch. A comparable scanner-order standardized 2585 collection was
unavailable in the local staging used for these replays. The 140 GiB proposal for the highest Galaxy class
still needs a full scanner-order run.

**Resource request for Galaxy (proposal, 10 cores)**

| Points | RAM | Scratch | Wall limit |
| --- | --- | --- | --- |
| ≤ 75 M | 32 GiB | 16 GiB | 1 h |
| 75–750 M | 64 GiB | 96 GiB | 2 h |
| 750 M–2.1 B | 96 GiB | 256 GiB | 5 h, estimated |
| 2.1–7.1 B | 140 GiB | 750 GiB | 10 h, estimated |

For wide grids, also size for the ground-to-canopy count cells, merge block depth and raster band count.

## 5. Limits and v1.1

- **TLS:** occlusion from provided scan positions is on branch `feat/tls-stations-v1.1`. Still open there:
  - merged clouds without scan IDs, such as 3950 (station-visibility mode);
  - a same-grid test for 2789.
- **SLAM paths** (GeoSLAM `_traj.txt`): the upload accepts them. The tool can use them only if the path time
  matches the points' GPS time. This is not yet checked on real data.
- **Livox with constant GPS time** (3516): refused, exit 3. The vendor keyframes match the points, but using them
  needs a decision.
- **Supplementary fallback:** trajectories that the upload or this tool cannot use should be kept as supplementary
  data on the dataset. They will be added when the format is supported.
- **Coverage** (production screen, 20 Nov):
  - 1,034 MLS/ULS datasets have multiple returns.
  - 57 are single-return with a scan angle.
  - 772 are refused: 408 single-return without scan angle (mostly GeoSLAM MLS), 364 with missing or constant
    GPS time.
  - 513 are unknown (no standardization metadata yet).
