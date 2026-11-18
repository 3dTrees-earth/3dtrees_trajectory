# Untwine writes the occlusion COPC. It comes from conda-forge with only the shared libraries it loads (found through
# its RPATH), so it never mixes with the Python wheels' own GDAL and PROJ.
FROM mambaorg/micromamba:1.5.10 AS untwine
ARG UNTWINE_VERSION=1.5.1
RUN micromamba create --yes --prefix /tmp/untwine-env --channel conda-forge "untwine=${UNTWINE_VERSION}" \
    && mkdir -p /tmp/untwine/bin /tmp/untwine/lib /tmp/untwine/share \
    && cp /tmp/untwine-env/bin/untwine /tmp/untwine/bin/ \
    && ldd /tmp/untwine-env/bin/untwine | awk '$3 ~ /^\/tmp\/untwine-env\// {print $3}' \
      | xargs cp -L -t /tmp/untwine/lib \
    && cp -r /tmp/untwine-env/share/proj /tmp/untwine/share/ \
    && /tmp/untwine/bin/untwine --version

FROM python:3.12.11-slim-bookworm@sha256:519591d6871b7bc437060736b9f7456b8731f1499a57e22e6c285135ae657bf7

ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
# Rasterio's bundled GDAL links against the system Expat library.
RUN apt-get update \
    && apt-get install -y --no-install-recommends libexpat1 \
    && rm -rf /var/lib/apt/lists/*
COPY requirements.txt /src/requirements.txt
RUN python -m pip install --no-cache-dir --only-binary=:all: -r /src/requirements.txt \
    # Only OccPy's compiled voxel traversal (raytr) is used; its Python layer would pull open3d, pyvista and OpenCV.
    && python -m pip install --no-cache-dir --only-binary=:all: --no-deps occpy_ls==0.3.0
COPY --from=untwine /tmp/untwine /opt/untwine
COPY src/ /src/

# Galaxy runs its job script with /bin/sh and calls /src/run.py from the job directory.
WORKDIR /data
USER 65532:65532
