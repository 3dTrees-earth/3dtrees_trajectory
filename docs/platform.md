# 3Dtrees platform, Galaxy and releases

[Back to README](../README.md)

## Platform gate

Use the standardization metadata to check that GPS time is present and varying. Treat
`max(NumberOfReturns) > 1` as evidence of multi-returns. A single-return cloud needs usable scan angles or an
available delivered path for occlusion. Check for a scan-angle dimension and actual nonzero values, since a LAS
point format can contain an angle field filled with zero. Multi-return clouds can estimate a path without angles
when they have enough usable pulse lines. Preserve GPS time, ReturnNumber, NumberOfReturns, ScanAngle
(or ScanAngleRank for legacy formats), PointSourceID and scanner channel wherever the input format provides
them. The tool checks usable pulses, the multi-return share, return pairing,
and static-scan evidence and prints its refusal to stderr. The platform gate can route clearly unsuitable
datasets before scheduling a job.

The 3Dtrees platform (`trees_api/processing/trajectory.py`) runs the tool for ULS and MLS collections only, when every
standardized file has varying GPS time and either multiple returns, a single return with a non-zero scan angle, or a
standardized timed delivered path. All active delivered paths of the dataset are passed as one list. Paths whose
coordinates cannot be aligned with the cloud (local against EPSG, judged from the standardized summary's CRS and the
upload's EPSG) are dropped before the run with the reason `trajectory_path_crs_mismatch`. TLS is refused
(`tls_requires_trajectory_v1_1`).

Products: `trajectory` (CSV) and `occlusion` (COPC, COG) as separate generations; the occlusion COG appears
in the viewer's Raster menu and in downloads. The default workflow runs the tool after Raster with its DTM;
`refresh_trajectory` reruns it on the standardized collection.

## Galaxy and containers

- Wrapper: bgruening/galaxytools `tools/3dtrees_trajectory` (see [cli.md](cli.md#galaxy-wrapper)).
- Resources: Galaxy Europe TPV (see [resources.md](resources.md#galaxy-europe-tpv)).
- Image: `ghcr.io/3dtrees-earth/3dtrees_trajectory:<version>`; Galaxy converts it to Singularity.

```bash
make test
make image
```

The Dockerfile pins both base image digests, the Untwine and system Expat versions, and the Linux Python runtime packages.

The Docker suite runs with two CPUs and a 4 GiB limit by default. Tests generate small LAS/LAZ clouds with known positions and verify the coordinates, refusals, resource reporting,
and byte-for-byte stability across decoder thread counts.

Call `python /src/run.py` from a writable job directory. The Docker image defaults to a non-root user. For a
local bind-mounted job directory, run under the job owner's UID and GID.

```bash
docker run --rm --name trajectory-local --user "$(id -u):$(id -g)" \
  --label "owner=$(id -un)" --label created-by=trajectory-local \
  --mount "type=bind,source=$PWD,target=/data" \
  3dtrees_trajectory:dev python /src/run.py standardized/strip1.las standardized/strip2.laz \
    -o trajectory/trajectory.csv
```

The Galaxy wrapper should accept one collection of standardized LAS/LAZ files per dataset and pass all its
members to one tool job. Pass `--threads "${GALAXY_SLOTS:-4}"` from the job allocation and collect the requested
CSV and occlusion products.
Keep sampling and platform diagnostics in the CLI. File metadata and geometry refusals are checked at runtime.
CI tests the image on pull requests, pushes to `main`, and pushes of `v*` tags. After the tests pass, it
publishes `ghcr.io/3dtrees-earth/3dtrees_trajectory:sha-<short-sha>` and either `:latest` for `main` or
the package version (for example, `:1.0.0`) for a matching release tag. Galaxy converts that Docker image
to Singularity when preparing the tool container.

## Releases

A release is a `v<version>` tag on `main`; CI tests it and publishes the image tag. Before a version is installed on
Galaxy Europe a tag may be re-cut (moved) with fixes; after installation any change needs a new version, because
Galaxy caches images by tag. See [CHANGELOG.md](../CHANGELOG.md).
