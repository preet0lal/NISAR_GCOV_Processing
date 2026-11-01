#!/usr/bin/env bash
# MAAP DPS build step. Use as the Build Command:
# NISAR_GCOV_Processing/build_env.sh
set -euo pipefail

echo "MAAP DPS dependency setup"
python --version

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

if ! command -v gdalwarp >/dev/null 2>&1; then
    echo "WARNING: gdalwarp is not on PATH; --lia yes requires the GDAL CLI." >&2
else
    echo "GDAL warp tool: $(command -v gdalwarp)"
fi
