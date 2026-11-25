# Resources

[Back to README](../README.md)

## Threads, memory and scratch

- The reader decodes up to five million points per batch. The job caps decoding, path fitting and window workers
  at the configured `GALAXY_SLOTS` or `--threads`, CPU affinity and cgroup quota. The default is ten slots.
  Large files decode independent record ranges and merge them in scanner order. The path solver uses threads
  for pulse-line fits and separate processes for scan-angle fits. It merges both in time-window order.
  Compressed LAS 1.4 point formats 6–10 decode only
  coordinates, GPS time, return numbers, scanner channel, point source ID and scan angle during path
  estimation. Window tracing skips scan angle, and ordinary extent rereads use only positions and GPS time
  plus the LAZ base layer. Point formats 0–5 use full decoding.
- **All pulses (default).** For plain LAS/LAZ without a delivered path, worker processes read seekable record
  ranges per GPS-time window and fit each window as soon as its pulses are complete; the solver keeps a compact
  normal system per window, so memory does not grow with the input. Echoes stored apart or overlapping files spill
  compact window records to scratch beside the job (up to about 85 GB for a 7 B-point flight such as 2585).
  COPC input, a delivered path, or a pairing fallback use the in-memory reader, which keeps a bounded sample
  (about 40 million multi-return and 10 million single-return pulses, or 40 million single-return pulses for
  single-return acquisitions) and says so in a CSV warning. `--sample` uses that sample everywhere.

  | Dataset (10 cores) | sample: wall / CPU / peak RSS | all pulses: wall / CPU / peak RSS |
  | --- | --- | --- |
  | 1844 (7.7 GB, 678 M points, 419 M pulses) | 31 s / 261 s / 1.65 GiB | 35 s / 315 s / 0.46 GiB |
  | 2585 (70 GB, 7.04 B points, 5.2 B pulses), raw files | — | 533 s, 0.50 GiB process RSS, 4.6 GiB over all workers, ~83 GB scratch |

  Before streaming (1.0.0 first cuts) all pulses on 1844 took 130 s and 34.5 GiB. Paths from all pulses and from
  the sample differ by a median of 0.1–0.2 mm (p95 2–5 mm); outputs of the streamed fit are byte-identical to the
  whole-cloud fit.
- GPS-window scratch contains durable compact counts at 12 bytes per stored voxel plus output intermediates.
  Count arrays are sparse files until the walker touches those cells. The debug report records recovered and unpaired
  echoes, executed and reused windows, per-window RSS, sampled aggregate RSS, allocated scratch, peak worker
  concurrency and memory reservations.
- `resources` records wall time, CPU, peak RSS, points read, decoder threads and cgroup measurements.
  Cgroup memory includes file cache.

- Positions cover sampled acquisition times. GPS week rolls, overlapping independent flights with the same
  clock and identifiers, and widely separated acquisitions should be split before processing.

## Voxel size

| Dataset 1844 (ULS, 7.7 GB, 678 M points), 10 cores | 0.5 m | 0.25 m |
| --- | ---: | ---: |
| wall / CPU | 137 s / 1083 s | 304 s / 1993 s |
| peak RSS (all processes) | 6.3 GiB | 35.9 GiB |
| container peak incl. page cache | 10.9 GB | 58.1 GB |
| COG bands / COPC / COG size | 97 / 132 MB / 7.1 MB | 191 / 476 MB / 45 MB |
| occluded share | 0.49 % | 0.94 % |

Memory grows about (0.5/voxel)^2.5 below 0.5 m. Very large flights at 0.25 m can exceed 140 GB and then exit 127.

## Galaxy Europe (TPV)

`usegalaxy-eu/infrastructure-playbook` `files/galaxy/tpv/tools.yml`, from the 1.0.0 CPU-scaling benchmark:

```yaml
cores: 6 if input_size <= 5 else 10
mem: min(140, 20 + (input_size * 0.55))            # GB, input = standardized LAZ collection
rules:                                             # voxels under 0.5 m
  mem: min(140, 20 + (input_size * 0.55) * (0.5 / voxel) ** 2.75)
```

Runs on Singularity, converted from the Docker image (`explicit_singularity`). Products are byte-identical at any
slot count; six slots cost little on small inputs and get scheduled sooner.

| Dataset | LAZ | 6 slots | 10 slots | Peak anon RSS |
| --- | ---: | ---: | ---: | ---: |
| 1807 ULS | 0.61 GB | 30.7 s | 25.1 s | 3.8 GiB |
| 2014 ULS | 0.48 GB | 61.5 s | 49.6 s | 4.3 GiB |
| 1844 ULS | 7.66 GB | 353.8 s | 229.8 s | 3.9 GiB |
| 1788 MLS | 0.11 GB | 113.4 s | 91.6 s | 2.8 GiB |
| 1804 ULS | 10.5 GB | 505.9 s | — | 3.7 GiB |
| 2585 ULS | 70 GB, 7.04 B points | — | 2 h 11 min | < 50 GiB cap |

## Earlier resource classes and replays

Allocate ten CPU slots for trajectory and occlusion together. Start with the total standardized point count, then
check count cells in the ground-to-canopy envelope, the extent of the merge blocks and the number of raster bands.

| Input points | RAM (GiB) | Scratch (GiB) | Wall limit |
| ---: | ---: | ---: | ---: |
| Up to 75 million | 32 | 16 | 1 h |
| >75 million to 750 million | 64 | 96 | 2 h |
| >750 million to 2.1 billion | 96 | 256 | 5 h |
| >2.1 billion to 7.1 billion | 96 | 256 | 6 h |

The standardized LAZ replay for 2027 (29.7 M points) took 78 s and peaked at 3.0 GiB of cgroup memory.
Historical scanner-order LAZ benchmarks provide headroom estimates for the larger classes:

| Dataset | Points | Read stage (s) | Full job (s) | Peak cgroup RAM (GiB) | Peak scratch (GiB) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1807 | 52.7 M | 9.0 | 101 | 2.7 | 2.2 |
| 2014, multi-return replay | 52.7 M | 10.8 | 152 | 21.8 | 6.5 |
| 2014, single-return replay | 33.8 M | 6.4 | 153 | 24.0 | 7.5 |
| 1844 | 678.1 M | 101.4 | 994 | 34.5 | 32.5 |
| 2585 | 7.04 B | 1,238.8 | 27,953 | 190.0 | 567.9 |

The 2014 single-return replay reached 24.01 GiB, so the smallest class allows 32 GiB. The 1844 replay
fits the 64 GiB class with room for the parent process and file cache. The 750 M–2.1 B class is interpolated.
The round-10 raw 2585 replay took 2 h 11 min with ten CPU slots, a 50 GiB cap and about 50 GiB of
allocated scratch. A replay from completed checkpoints measured 1 min 13 s for the final merge.
The resource class is an estimate for acquisitions with a comparable ground-to-canopy footprint.
A standardized 2585 replay is still required for Galaxy validation.

The standardized 1788 collection has 24.7 M points and a 2,000 × 2,000 × 1,278 grid at 0.5 m.
Use `--debug` for current window memory and scratch measurements. [Stage resources](REPORT.md#3-stage-resources)

A replay of 36 archived 3Dtrees ULS/MLS datasets matched the time and XYZ bytes of all 35 prior paths.
Dataset 2018 reproduced its refusal because it has single returns without a usable scan angle.
The archive and replay ran separately on the 3Dtrees compute server. Replay jobs used eight decoder slots with a ten-CPU affinity mask.

| Dataset | Points read | Archived → replay wall time | Archived → replay peak RSS |
| --- | ---: | ---: | ---: |
| 1807 YellowScan | 52.7 M | 23.14 → 11.28 s | 2.21 → 1.06 GB |
| 2018 MLS | 1.95 B | 188.48 → 161.52 s | 6.86 → 2.70 GB |
| 2585 DJI | 7.04 B | 646.63 → 712.40 s | 3.57 → 2.34 GB |

Across all 36 jobs the median wall time was 17.06 s, compared with 26.97 s in the archive. Peak replay process
RSS was 3.56 GB. The largest input charged 77.30 GB to its cgroup, which also accounted for file cache.
For 2585, the delivered-path median horizontal error remains 37.8 m. Quality labels measure pulse fit residuals.
