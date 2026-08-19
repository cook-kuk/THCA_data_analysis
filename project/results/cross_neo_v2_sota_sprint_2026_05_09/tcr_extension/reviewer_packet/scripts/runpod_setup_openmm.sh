#!/usr/bin/env bash
set -euo pipefail
python -m pip install --upgrade pip
python -m pip install openmm mdtraj git+https://github.com/openmm/pdbfixer.git
python - <<'PY'
from openmm import Platform
print([Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
PY
