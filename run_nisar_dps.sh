#!/usr/bin/env bash
set -euo pipefail
BASEDIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p input output
SOURCE="${1:-}"
if [[ -z "$SOURCE" ]]; then
  shopt -s nullglob
  files=(input/*.h5 input/*.hdf5 inputs/*.h5 inputs/*.hdf5 /inputs/*.h5 /inputs/*.hdf5)
  shopt -u nullglob
  if (( ${#files[@]} == 1 )); then SOURCE="${files[0]}"; fi
fi
if [[ -z "$SOURCE" ]]; then
  echo "ERROR: Empty GCOV input. MAAP did not supply 's3-url' or stage an HDF5." >&2
  exit 2
fi
args=(--s3_url "$SOURCE")
if (( $# >= 5 )) && [[ -n "${2:-}" && -n "${3:-}" && -n "${4:-}" && -n "${5:-}" ]]; then
  args+=(--bbox "$2" "$3" "$4" "$5")
fi
python "${BASEDIR}/nisar_dps_wrapper.py" "${args[@]}"
