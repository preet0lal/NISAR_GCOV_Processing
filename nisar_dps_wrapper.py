#!/usr/bin/env python3

import argparse
import os
import subprocess
import shutil
import sys
from pathlib import Path
from urllib.parse import urlparse, unquote

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
        default=None,
        help="Optional clipping bounding box; omitted means native GCOV footprint",
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

    source = args.s3_url.strip()
    if not source:
        parser.error("--s3_url must not be empty")
    local_source = Path(unquote(urlparse(source).path) if source.startswith("file://") else source)
    is_s3 = source.startswith("s3://")
    local_h5 = input_dir / (Path(urlparse(source).path).name if is_s3 else local_source.name)
    if local_h5.suffix.lower() not in (".h5", ".hdf5"):
        parser.error("GCOV input must end with .h5 or .hdf5")
    created_local_copy = False
    if is_s3:
        print("Downloading:", source, "->", local_h5, flush=True)
        fs.get(source, str(local_h5))
        created_local_copy = True
    else:
        if not local_source.is_file():
            parser.error("MAAP-staged GCOV file not found: " + source)
        if local_source.resolve() != local_h5.resolve():
            shutil.copyfile(local_source, local_h5)
            created_local_copy = True
        print("Using staged GCOV:", local_source, flush=True)

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

        "--lia", "yes",
        "--nisar_dem_vrt", dem_vrt,
        "--dem_cache_dir", str(dem_cache),

        # ONE granule/frame per DPS job
        "--mosaic", "no",

        "--workers", "1",
        "--tile_workers", "4",

        "--keep_intermediate", "no",
    ]

    if args.bbox is not None:
        cmd += ["--bbox", *(str(v) for v in args.bbox)]

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
        if created_local_copy and local_h5.exists():
            local_h5.unlink()

    print("\nDPS processing complete.")
    print("Output:", output_dir)


if __name__ == "__main__":
    main()
