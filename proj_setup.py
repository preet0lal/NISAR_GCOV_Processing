#!/usr/bin/env python3
"""Locate a usable PROJ database for the MAAP GDAL/Rasterio environment.

Avoid hardcoding a conda prefix: MAAP may use /opt/conda or /srv/conda.
Every candidate is tested using the actual Python/Rasterio installation.
"""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import subprocess
import sys


_CHECK = (
    "from rasterio.crs import CRS; "
    "from rasterio.warp import transform; "
    "assert CRS.from_epsg(6933).is_projected; "
    "x, y = transform('EPSG:4326', 'EPSG:6933', [78.0], [22.0]); "
    "assert len(x) == len(y) == 1; "
    "print('EPSG:6933 and WGS84 transformations verified')"
)


def _candidate_dirs(env: dict[str, str]):
    # Prioritize PROJ database provided with Rasterio to avoid wheel/conda
    # version mismatches. Test every path against the real Rasterio library.
    for module, relative in (
        ("rasterio", "proj_data"),
        ("pyproj", "proj_dir/share/proj"),
    ):
        try:
            spec = importlib.util.find_spec(module)
        except (ImportError, ModuleNotFoundError, ValueError):
            spec = None
        if spec is not None and spec.submodule_search_locations:
            package = Path(next(iter(spec.submodule_search_locations)))
            yield package / relative

    for prefix in (
        sys.prefix,
        env.get("CONDA_PREFIX"),
        env.get("VIRTUAL_ENV"),
        "/opt/conda",
        "/srv/conda",
        "/usr/local",
        "/usr",
    ):
        if prefix:
            yield Path(prefix) / "share/proj"

    for var in ("PROJ_DATA", "PROJ_LIB"):
        for value in env.get(var, "").split(os.pathsep):
            if value:
                yield Path(value)


def _attempt(env: dict[str, str]):
    try:
        proc = subprocess.run(
            [sys.executable, "-c", _CHECK],
            env=env,
            text=True,
            capture_output=True,
            timeout=30,
            check=False,
        )
        return proc.returncode == 0, (proc.stderr or proc.stdout).strip()
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, str(exc)


def configure_projection_environment(
    environ: dict[str, str] | None = None, verbose: bool = True
) -> dict[str, str]:
    """Return a copy of the environment where EPSG:6933 works in Rasterio.

    Deliberately do not modify global os.environ or guess that a database at
    some path is compatible; validate an actual CRS and coordinate transform.
    """
    base = dict(os.environ if environ is None else environ)
    # An invalid GDAL_DATA override can also disrupt a bundled GDAL wheel.
    if base.get("GDAL_DATA") and not Path(base["GDAL_DATA"]).is_dir():
        base.pop("GDAL_DATA", None)

    ok, last_error = _attempt(base)
    if ok:
        if verbose:
            print("PROJ check: EPSG:6933 already works.", flush=True)
        return base

    clean = base.copy()
    clean.pop("PROJ_DATA", None)
    clean.pop("PROJ_LIB", None)
    ok, last_error = _attempt(clean)
    if ok:
        if verbose:
            print("PROJ check: removed stale PROJ_DATA/PROJ_LIB overrides.", flush=True)
        return clean

    seen = set()
    for candidate in _candidate_dirs(base):
        candidate = candidate.resolve()
        if candidate in seen or not (candidate / "proj.db").is_file():
            continue
        seen.add(candidate)
        trial = clean.copy()
        # PROJ >= 9.1 uses PROJ_DATA; older libraries use PROJ_LIB.
        trial["PROJ_DATA"] = str(candidate)
        trial["PROJ_LIB"] = str(candidate)
        ok, last_error = _attempt(trial)
        if ok:
            if verbose:
                print(f"PROJ check: using validated database {candidate}", flush=True)
            return trial

    raise RuntimeError(
        "No compatible proj.db found for this Python/Rasterio environment; "
        "checked Rasterio and pyproj package data and standard conda/system "
        "paths. Build the MAAP image with a compatible PROJ database. "
        f"Last CRS test error: {last_error[-800:]}"
    )


if __name__ == "__main__":
    configure_projection_environment()
