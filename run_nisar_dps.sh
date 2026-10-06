#!/usr/bin/env bash

set -euo pipefail

BASEDIR="$(cd "$(dirname "$0")" && pwd)"

mkdir -p input
mkdir -p output

echo "========================================"
echo "NISAR GCOV → EASE2 100 m"
echo "MAAP DPS processing"
echo "========================================"

echo "S3 URL : $1"
echo "BBOX   : $2 $3 $4 $5"

python "${BASEDIR}/nisar_dps_wrapper.py" \
    --s3_url "$1" \
    --bbox "$2" "$3" "$4" "$5"

echo "========================================"
echo "DONE"
echo "========================================"
