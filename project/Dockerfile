# THCA v5.1 + v8 reproducibility container (audit fix F10)
# ---------------------------------------------------------
# Base:    python:3.12-slim (Debian bookworm-slim)
# Build:   docker build -t thyrai/dial-v5.1 .
# Run:     docker run --rm thyrai/dial-v5.1 make reproduce
#
# A reviewer running `make reproduce` inside the container should see the
# canonical results/v5/v5p1_dial_all_cancers.tsv regenerated and the
# tests/check_reproducibility.py harness print "REPRODUCIBILITY CHECK: PASS".
#
# Cache:   Set V5P1_CACHE_DIR=/data and mount a host volume there so re-runs
#          do not re-download GDC/GEO blobs (~10 GB).

FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    V5P1_CACHE_DIR=/data \
    PROJECT_ROOT=/opt/thyroid-dash/project

# Minimal build-time deps for scientific wheels (numpy/scipy/pandas BLAS,
# pillow, matplotlib fonts, plus `make` for the Makefile entrypoint).
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        make \
        libgl1 \
        libglib2.0-0 \
        libgomp1 \
        curl \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/thyroid-dash/project

# Install the FROZEN environment first (cache-friendly layer).
# tests/requirements.lock is the output of `pip freeze` against the venv
# that produced the canonical results/v5/v5p1_dial_all_cancers.tsv.
COPY project/tests/requirements.lock /tmp/requirements.lock
RUN python -m pip install --upgrade pip \
 && python -m pip install -r /tmp/requirements.lock

# Copy code, makefile, tests. Do NOT copy data_raw/ or .venv/ or logs/.
COPY project/notebooks_or_scripts /opt/thyroid-dash/project/notebooks_or_scripts
COPY project/scripts              /opt/thyroid-dash/project/scripts
COPY project/metadata             /opt/thyroid-dash/project/metadata
COPY project/results              /opt/thyroid-dash/project/results
COPY project/Makefile             /opt/thyroid-dash/project/Makefile
COPY project/tests                /opt/thyroid-dash/project/tests

# Create scratch dirs for runtime outputs and the GDC cache mount-point.
RUN mkdir -p /opt/thyroid-dash/project/data_raw \
             /opt/thyroid-dash/project/data_processed \
             /opt/thyroid-dash/project/logs \
             /opt/thyroid-dash/project/state \
             /data/thca

# Persist the GDC/GEO download cache across container restarts.
VOLUME ["/data"]

CMD ["make", "reproduce"]
