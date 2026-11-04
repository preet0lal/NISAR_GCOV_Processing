#!/usr/bin/env python3

import argparse
import h5py
import os
import subprocess
import sys
from pathlib import Path

import s3fs
from maap.maap import MAAP
from proj_setup import configure_projection_environment


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--s3_url", required=True)

    parser.add_argument(
        "--bbox",
        nargs=4,
        type=float,
        metavar=("W", "S", "E", "N"),
        required=True,
    )

    args = parser.parse_args()

    basedir = Path(__file__).resolve().parent

    input_dir = Path("input").resolve()
    output_dir = Path("output").resolve()

    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Validate EPSG:6933 before acquiring credentials or downloading large HDF5.
    # Inherit the tested PROJ paths in both Python and the GDAL subprocess.
    proj_env = configure_projection_environment()

    # ---------------------------------------------------------
    # NISAR temporary S3 credentials
    # Works in MAAP ADE / DPS
    # ---------------------------------------------------------

    maap = MAAP()

    cred_url = (
        "https://nisar.asf.earthdatacloud.nasa.gov/"
        "s3credentials"
    )

    creds = maap.aws.earthdata_s3_credentials(cred_url)

    fs = s3fs.S3FileSystem(
        anon=False,
        key=creds["accessKeyId"],
        secret=creds["secretAccessKey"],
        token=creds["sessionToken"],
        client_kwargs={"region_name": "us-west-2"},
    )

    # ---------------------------------------------------------
    # Download ONE GCOV granule to this DPS worker
    # ---------------------------------------------------------

    local_h5 = input_dir / Path(args.s3_url).name

    print("Downloading:")
    print(args.s3_url)
    print(" ->", local_h5)

    fs.get(
        args.s3_url,
        str(local_h5),
    )

    print(
        f"Downloaded: "
        f"{local_h5.stat().st_size / 1024**3:.2f} GB"
    )

    # ---------------------------------------------------------
    # Credentials also available to GDAL / VSIS3
    # ---------------------------------------------------------

    env = proj_env.copy()

    env.update(
        AWS_ACCESS_KEY_ID=creds["accessKeyId"],
        AWS_SECRET_ACCESS_KEY=creds["secretAccessKey"],
        AWS_SESSION_TOKEN=creds["sessionToken"],
        AWS_REGION="us-west-2",
        AWS_DEFAULT_REGION="us-west-2",
        AWS_REQUEST_PAYER="requester",

        OMP_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        NUMEXPR_NUM_THREADS="1",
    )

    processor = basedir / "NISAR_GCOV_TO_TIFF_E.py"

    dem_vrt = (
        "/vsis3/"
        "sds-n-cumulus-prod-nisar-products/"
        "DEM/v1.2/EPSG4326/EPSG4326.vrt"
    )

    # Local to this DPS worker only
    dem_cache = input_dir / "DEM_Cache"

    w, s, e, n = args.bbox

    # Inspect available diagonal terms independently for each GCOV frequency.
    # Run each compatible pair or single channel in a separate subfolder.
    # This prevents frame-count, LIA, and quicklook filename collisions.
    modes = {
        "HH_HV": ("HHHH", "HVHV"),
        "VV_VH": ("VVVV", "VHVH"),
        "RH_RV": ("RHRH", "RVRV"),
    }
    single_modes = {
        "HHHH": "HH_ONLY", "HVHV": "HV_ONLY",
        "VVVV": "VV_ONLY", "VHVH": "VH_ONLY",
        "RHRH": "RH_ONLY", "RVRV": "RV_ONLY",
    }
    with h5py.File(local_h5, "r") as h5:
        root = "LSAR" if local_h5.name.startswith("NISAR_L") else "SSAR"
        grids_path = f"/science/{root}/GCOV/grids"
        if grids_path not in h5:
            raise RuntimeError(f"Missing GCOV grids: {grids_path}")
        grids = h5[grids_path]
        tasks = []
        for frequency in ("frequencyA", "frequencyB"):
            if frequency not in grids:
                continue
            present = {p for p in single_modes if p in grids[frequency]}
            consumed = set()
            for mode, (co, cross) in modes.items():
                if co in present and cross in present:
                    tasks.append((frequency, mode, (co, cross)))
                    consumed.update((co, cross))
            for term in sorted(present - consumed):
                tasks.append((frequency, single_modes[term], (term,)))

    if not tasks:
        raise RuntimeError("No supported diagonal GCOV intensity datasets in either frequency")

    common_cmd = [
        sys.executable, str(processor),
        "--input_dir", str(input_dir),
        "--normalization", "gamma0",
        "--filter", "enhanced_lee",
        "--filter_window", "7",
        "--grid_mode", "EASE2",
        "--bbox", str(w), str(s), str(e), str(n),
        "--lia", "yes",
        "--nisar_dem_vrt", dem_vrt,
        "--dem_cache_dir", str(dem_cache),
        "--mosaic", "no",
        "--workers", "1",
        "--tile_workers", "4",
        "--keep_intermediate", "no",
    ]

    # Nonzero exit from any frequency/polarization task fails the whole DPS job.
    # A per-task folder isolates output and proves which channels were generated.
    try:
        for frequency, mode, terms in tasks:
            destination = output_dir / frequency / mode
            destination.mkdir(parents=True, exist_ok=True)
            cmd = common_cmd + [
                "--output_dir", str(destination),
                "--frequency", frequency,
                "--pol_pair", mode,
            ]
            print("\\nRunning:", " ".join(cmd), flush=True)
            subprocess.run(cmd, check=True, env=env)
            for term in terms:
                matches = list(destination.glob(f"*_{term}.tif"))
                if not matches:
                    raise RuntimeError(
                        f"Processor returned success but {frequency}/{mode}/{term} output is missing"
                    )
            print(f"Verified {frequency}/{mode}: {', '.join(terms)}", flush=True)
    finally:
        if local_h5.exists():
            local_h5.unlink()

    print("\\nDPS processing complete.")
    print("Output:", output_dir)


if __name__ == "__main__":
    main()
