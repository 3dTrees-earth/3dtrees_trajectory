# Untwine writes the occlusion COPC. It comes from conda-forge with only the shared libraries it loads (found through
# its RPATH), so it never mixes with the Python wheels' own GDAL and PROJ.
FROM mambaorg/micromamba:1.5.10@sha256:e3797091302382ea841498bc93a7b0a50f7c1448333d5e946d2d1608d0c5f43d AS untwine
RUN micromamba create --yes --prefix /tmp/untwine-env --override-channels --channel conda-forge "untwine=1.5.1" \
    && mkdir -p /tmp/untwine/bin /tmp/untwine/lib /tmp/untwine/share \
    && cp /tmp/untwine-env/bin/untwine /tmp/untwine/bin/ \
    && ldd /tmp/untwine-env/bin/untwine | awk '$3 ~ /^\/tmp\/untwine-env\// {print $3}' \
      | xargs cp -L -t /tmp/untwine/lib \
    && cp -r /tmp/untwine-env/share/proj /tmp/untwine/share/ \
    && /tmp/untwine/bin/untwine --version

FROM python:3.12.11-slim-bookworm@sha256:519591d6871b7bc437060736b9f7456b8731f1499a57e22e6c285135ae657bf7 AS ray-walker
RUN apt-get update && apt-get install -y --no-install-recommends g++
COPY src/threedtrees_trajectory/ray_walk.cpp /tmp/ray_walk.cpp
RUN g++ -std=c++17 -O3 -fno-fast-math -ffp-contract=off -fPIC -shared -static-libstdc++ -static-libgcc \
    /tmp/ray_walk.cpp -o /tmp/_ray_walk.so

FROM python:3.12.11-slim-bookworm@sha256:519591d6871b7bc437060736b9f7456b8731f1499a57e22e6c285135ae657bf7 AS runtime

ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
# Rasterio's bundled GDAL links against the system Expat library.
RUN apt-get update \
    && apt-get install -y --no-install-recommends libexpat1=2.5.0-1+deb12u5 \
    && rm -rf /var/lib/apt/lists/*
COPY requirements.txt /src/requirements.txt
RUN python -m pip install --no-cache-dir --only-binary=:all: --no-deps -r /src/requirements.txt
COPY --from=untwine /tmp/untwine /opt/untwine
COPY src/ /src/
COPY --from=ray-walker /tmp/_ray_walk.so /src/threedtrees_trajectory/_ray_walk.so

# Galaxy runs its job script with /bin/sh and calls /src/run.py from the job directory.
WORKDIR /data
USER 65532:65532

FROM runtime AS test
USER root
COPY requirements-dev.txt /src/requirements-dev.txt
RUN python -m pip install --no-cache-dir --only-binary=:all: -r /src/requirements-dev.txt \
    && python -m pip install --no-cache-dir --only-binary=:all: --no-deps occpy_ls==0.3.0
COPY tests/ /test-suite/tests/
RUN ln -s /src /test-suite/src
ENV PYTHONPATH=/src
WORKDIR /test-suite
USER 65532:65532
CMD ["python", "-m", "pytest", "-q", "-p", "no:cacheprovider", "tests"]

FROM runtime AS release
