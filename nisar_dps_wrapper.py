#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
from pathlib import Path

import s3fs
from maap.maap import MAAP


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

    env = os.environ.copy()

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

    cmd = [
        sys.executable,
        str(processor),

        "--input_dir", str(input_dir),
        "--output_dir", str(output_dir),

        "--frequency", "frequencyA",
        "--pol_pair", "HH_HV",
        "--normalization", "gamma0",
        "--filter", "none",

        "--grid_mode", "EASE2",

        "--bbox",
        str(w), str(s), str(e), str(n),

        "--lia", "yes",
        "--nisar_dem_vrt", dem_vrt,
        "--dem_cache_dir", str(dem_cache),

        # ONE granule/frame per DPS job
        "--mosaic", "no",

        "--workers", "1",
        "--tile_workers", "4",

        "--keep_intermediate", "no",
    ]

    print("\nRunning:")
    print(" ".join(cmd))

    try:

        subprocess.run(
            cmd,
            check=True,
            env=env,
        )

    finally:

        # Raw 7.5-GB HDF5 is not an output product.
        # Remove it before DPS packages results.
        if local_h5.exists():
            local_h5.unlink()

    print("\nDPS processing complete.")
    print("Output:", output_dir)


if __name__ == "__main__":
    main()
