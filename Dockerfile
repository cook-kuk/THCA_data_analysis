# THCA Multi-Omics Dashboard — reproducible runtime container
# ------------------------------------------------------------
# Base:    python:3.12-slim (Debian bookworm-slim)
# Size:    ~1.4 GB after installs (torch + qiskit dominate)
# Ports:   8012 (serve_secure.py)
# Volumes: /app/project/data_raw, /app/project/data_processed
# Build:   docker build -t thca-dashboard .
# Run:     docker run --rm -p 8012:8012 thca-dashboard

FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PORT=8012 \
    BIND=0.0.0.0

# Minimal build-time deps for scientific wheels (numba/LLVM, pillow, matplotlib fonts)
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        libgl1 \
        libglib2.0-0 \
        libgomp1 \
        curl \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Python dependencies first so layer is cache-friendly.
COPY requirements.txt /app/requirements.txt
RUN python -m pip install --upgrade pip \
 && python -m pip install --index-url https://download.pytorch.org/whl/cpu torch==2.11.0+cpu \
 && python -m pip install -r /app/requirements.txt

# Copy only project sub-trees we need — NOT .venv, NOT data_raw/, NOT logs/.
COPY project/notebooks_or_scripts /app/project/notebooks_or_scripts
COPY project/scripts              /app/project/scripts
COPY project/metadata             /app/project/metadata
COPY project/reports              /app/project/reports
COPY project/results              /app/project/results
COPY project/qc                   /app/project/qc
COPY project/index.html           /app/project/index.html

# Volume mount-points for user-supplied data
RUN mkdir -p /app/project/data_raw /app/project/data_processed /app/project/logs /app/project/state
VOLUME ["/app/project/data_raw", "/app/project/data_processed"]

EXPOSE 8012

WORKDIR /app/project
ENTRYPOINT ["python", "/app/project/scripts/serve_secure.py"]
