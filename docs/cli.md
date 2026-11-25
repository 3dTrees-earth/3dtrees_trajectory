# Command line

[Back to README](../README.md)

```text
3dtrees-trajectory CLOUD [CLOUD …] [-o CSV] [--geojson PATH] [--trajectory CSV [CSV …]] [--dtm TIF]
                   [--occlusion-copc PATH --occlusion-cog PATH] [--voxel M] [--interval S] [--sample]
                   [--threads N] [--grid-memory-gib G] [--no-platform-check] [--debug [--report PATH]] [--version]
```

| Option | Default | Range | Meaning |
| --- | --- | --- | --- |
| `CLOUD …` | | 1–65,536 files | standardized LAS/LAZ files of one acquisition, in collection order |
| `-o, --output` | | | trajectory CSV (path estimate) |
| `--geojson` | next to the CSV | | map layer |
| `--trajectory` | | 1..n files | delivered paths for occlusion ([delivered-paths.md](delivered-paths.md)) |
| `--dtm` | | one GeoTIFF | terrain for height above ground |
| `--occlusion-copc`, `--occlusion-cog` | | both together | occlusion products |
| `--voxel` | 0.5 | 0.05–10 m | voxel edge and COG band step |
| `--interval` | 0.5 | > 0 s | seconds between path positions |
| `--sample` | off | | fit the path from a bounded sample of pulses instead of all pulses; `--full` (all pulses) is the default and still accepted |
| `--threads` | `GALAXY_SLOTS` or 10 | ≥ 1 | decoder threads and GPS-time window workers |
| `--grid-memory-gib` | 16 | ≤ 64 | soft per-worker budget before further GPS-time splitting |
| `--no-platform-check` | off | | skip the static (TLS) check for a known moving acquisition |
| `--debug`, `--report` | off | | JSON quality and resource report |

The CSV, GeoJSON, report and input paths must be distinct. Without `-o` the tool checks the cloud and writes only the
occlusion products from the delivered paths.

## Galaxy wrapper

The Galaxy tool **3Dtrees: Trajectory and Occlusion** (bgruening/galaxytools, `tools/3dtrees_trajectory`) exposes:
the point cloud collection, an optional DTM, optional delivered trajectories (multiple), voxel size (0.25–2 m), path sampling interval (0.1–5 s) and
*Use all pulses* (on by default; off passes `--sample`). Threads come from `GALAXY_SLOTS`;
internals (`--grid-memory-gib`, `--debug`, `--no-platform-check`) are not exposed. Errors are detected by exit code.
