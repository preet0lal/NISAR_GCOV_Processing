#!/usr/bin/env bash
# MAAP DPS build step. Use as the Build Command:
# NISAR_GCOV_Processing/build_env.sh
set -euo pipefail

# MAAP invokes this script from a parent working directory. Make sibling
# Python modules importable during the Docker build as well as at runtime.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${SCRIPT_DIR}${PYTHONPATH:+:${PYTHONPATH}}"

echo "MAAP DPS dependency setup"
python --version

# --lia yes calls the gdalwarp command-line program. Rasterio's pip wheel
# provides a Python library but NOT the gdalwarp executable.
if ! command -v gdalwarp >/dev/null 2>&1; then
    echo "gdalwarp unavailable; installing GDAL CLI via conda-forge."
    if ! command -v conda >/dev/null 2>&1; then
        echo "ERROR: conda unavailable; cannot install required gdalwarp." >&2
        exit 1
    fi
    conda install -y -c conda-forge gdal
fi
command -v gdalwarp
gdalwarp --version

# Preserve preinstalled geospatial libraries where possible; only add
# packages that are missing from the selected MAAP base container.
python - <<'PY'
import importlib
import importlib.util
import subprocess
import sys

required = {
    "maap": "maap-py",
    "s3fs": "s3fs",
    "numpy": "numpy",
    "h5py": "h5py",
    "rasterio": "rasterio",
    "scipy": "scipy",
    "matplotlib": "matplotlib",
}
missing = [package for module, package in required.items()
           if importlib.util.find_spec(module) is None]
if missing:
    print("Installing missing dependencies:", ", ".join(missing), flush=True)
    subprocess.check_call(
        [sys.executable, "-m", "pip", "install", "--no-cache-dir", *missing]
    )
else:
    print("All required Python dependencies are already present.", flush=True)

for module in required:
    importlib.import_module(module)
importlib.import_module("maap.maap")
print("MAAP DPS Python dependency verification passed.", flush=True)
PY

echo "GDAL warp tool: $(command -v gdalwarp)"

# Verify CRS support inside the built image. If there is no compatible
# projection database, install pyproj's bundled database and retry.
if ! python "${SCRIPT_DIR}/proj_setup.py"; then
    echo "PROJ database validation failed; adding pyproj's bundled PROJ data."
    python -m pip install --no-cache-dir pyproj
    python "${SCRIPT_DIR}/proj_setup.py"
fi

# Check the separate GDAL CLI against the same projection environment that
# the runtime wrapper passes into its DEM/gdalwarp subprocess.
python - <<'PY'
import subprocess
from proj_setup import configure_projection_environment

env = configure_projection_environment(verbose=False)
subprocess.run(["gdalsrsinfo", "EPSG:6933"], env=env,
               check=True, capture_output=True, text=True)
print("GDAL CLI EPSG:6933 validation passed.", flush=True)
PY
