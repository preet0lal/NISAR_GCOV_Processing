#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'USAGE'
Usage:
  bash NISAR_GCOV_TIFF_E.sh \
    --input_dir INPUT \
    --output_dir OUTPUT \
    [options]

Purpose:
  Convert NISAR L2 GCOV HDF5 files to fixed-grid GeoTIFFs/mosaics with
  optional incidence-angle and terrain-local-incidence-angle outputs.

  SAR GeoTIFF values remain LINEAR power. Quicklook SAR panels are displayed
  in dB using 10*log10(power). With --lia yes, the quicklook becomes 2x2:
      co-pol dB | cross-pol dB
      LIA (deg) | incidence angle (deg)

Grid modes:
  --grid_mode UTM|EASE2     Default: UTM

  UTM:
    Keeps the existing fixed projected-grid workflow.
    --target_epsg auto chooses the majority native NISAR EPSG.
    --scale controls the output spacing (default 100 m).

  EASE2:
    Uses the exact EASE-Grid 2.0 M0p1 global lattice:
      EPSG:6933
      spacing = 100.089502334956 m
      global size = 347040 columns x 146160 rows
    The scene is cropped from this global lattice by exact integer global
    row/column indices. --scale, --target_epsg and grid origins are not used
    to define the EASE2 lattice. Each EASE2 output set also gets one compact
    *_EASE2_grid.mat companion with global row/col and x/y/lon/lat center vectors.

Required:
  --input_dir PATH          Folder containing NISAR GCOV .h5/.hdf5 files
  --output_dir PATH         Output folder

Grid options:
  --grid_mode MODE          UTM | EASE2 (default UTM)
  --target_epsg VALUE       UTM mode only: auto or projected EPSG code
                            Default: auto
  --scale M                 UTM mode only: output pixel size (default 100 m)
                            EASE2 always uses exact M0p1 = 100.089502334956 m
  --grid_origin_x X         UTM mode only: fixed grid origin x (default 0)
  --grid_origin_y Y         UTM mode only: fixed grid origin y (default 0)
  --ease2_mat PATH          Optional EASE2_latlon100m.mat sanity check.
                            If supplied in EASE2 mode, the Python code verifies
                            that its lat/lon centers match exact M0p1 centers.
  --bbox "W S E N"          Optional clip box in EPSG:4326

Data options:
  --frequency NAME          frequencyA | frequencyB (default frequencyA)
  --pol_pair PAIR           auto | HH_HV | VV_VH | RH_RV (default auto)
  --normalization TYPE      gamma0 | sigma0 (default gamma0)
  --mode MODE               nofilt | enhanced_lee | gamma_map | all
                            Default: nofilt
  --window N                Speckle-filter window, odd >=3 (default 7)
  --noise_floor VALUE       Linear-power threshold; 0 disables (default 0)
  --min_coverage F          Minimum valid native area per output pixel
                            Default: 0.90
  --use_exception_mask VAL  yes|no. Default: yes
                            yes: exclude pixels where inputDataExceptionMask
                                 is nonzero (e.g. NISAR sample slip/shift)
                            no:  keep those pixels instead of masking them,
                                 provided they pass every other check
  --lia yes|no              Generate incidence-angle + local-incidence-angle
                            GeoTIFFs on the same final grid and SAR-valid mask
                            (default no)
  --lia_dem_scale M         Internal terrain scale used for LIA before final
                            aggregation (default 30 m)
  --nisar_dem_vrt PATH      Optional local/VSI/URL override for NISAR DEM VRT
  --dem_cache_dir PATH      Persistent DEM cache. Existing matching DEMs are
                            reused instead of downloaded/reprojected again.
                            Default: <input_dir>/DEM_Cache

Mosaic/performance options:
  --mosaic yes|no           Mosaic adjacent compatible frames (default yes)
  --gap_threshold KM        Maximum inter-frame gap (default 50)
  --workers N               Parallel acquisition groups (default 1)
  --tile_workers N          GDAL warp threads per group (default 4)
  --tile_size_km KM         Output processing window size (default 50)
  --native_block_px N       Native HDF5 block size (default 1024)
  --keep_intermediate yes|no
                            Keep temporary native rasters (default no)

Execution options:
  --script PATH             Python processor path
  --python_exec PATH        Python executable (default python3)
  --extra_args "ARGS"       Additional whitespace-separated Python arguments
  -h, --help

Examples:

  # Exact EASE2 M0p1 (~100 m), recommended for one common global grid
  bash NISAR_GCOV_TIFF_E.sh \
    --input_dir /path/to/gcov \
    --output_dir /path/to/output_ease2 \
    --grid_mode EASE2 \
    --mode nofilt --mosaic no --lia yes

  # EASE2 with optional validation against your ancillary MAT file
  bash NISAR_GCOV_TIFF_E.sh \
    --input_dir /path/to/gcov \
    --output_dir /path/to/output_ease2 \
    --grid_mode EASE2 \
    --ease2_mat /path/to/EASE2_latlon100m.mat \
    --mode nofilt --mosaic no --lia yes

  # Existing UTM/native-projection workflow at exactly 100 m
  bash NISAR_GCOV_TIFF_E.sh \
    --input_dir /path/to/gcov \
    --output_dir /path/to/output_utm \
    --grid_mode UTM \
    --target_epsg auto \
    --scale 100 --mode nofilt --mosaic no --lia yes
USAGE
}

INPUT_DIR=""
OUTPUT_DIR=""
GRID_MODE="UTM"
TARGET_EPSG="auto"
EASE2_MAT=""
MODE="nofilt"
MOSAIC="yes"
WORKERS="1"
TILE_WORKERS="4"
TILE_SIZE_KM="50"
NATIVE_BLOCK_PX="1024"
SCALE="100"
GRID_ORIGIN_X="0"
GRID_ORIGIN_Y="0"
FREQUENCY="frequencyA"
POL_PAIR="auto"
NORMALIZATION="gamma0"
WINDOW="7"
NOISE_FLOOR="0"
MIN_COVERAGE="0.90"
USE_EXCEPTION_MASK="yes"
LIA="no"
LIA_DEM_SCALE="30"
NISAR_DEM_VRT=""
DEM_CACHE_DIR=""
EARTHDATA_COOKIE="/tmp/gdal_cookies.txt"
GAP_THRESHOLD="50"
BBOX=""
KEEP_INTERMEDIATE="no"
EXTRA_ARGS=""
PYTHON_EXEC="python3"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/NISAR_GCOV_TO_TIFF_E.py"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input_dir) INPUT_DIR="$2"; shift 2 ;;
    --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
    --grid_mode|--grid) GRID_MODE="${2^^}"; shift 2 ;;
    --target_epsg) TARGET_EPSG="$2"; shift 2 ;;
    --ease2_mat) EASE2_MAT="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --mosaic) MOSAIC="$2"; shift 2 ;;
    --workers) WORKERS="$2"; shift 2 ;;
    --tile_workers) TILE_WORKERS="$2"; shift 2 ;;
    --tile_size_km) TILE_SIZE_KM="$2"; shift 2 ;;
    --native_block_px) NATIVE_BLOCK_PX="$2"; shift 2 ;;
    --scale) SCALE="$2"; shift 2 ;;
    --grid_origin_x) GRID_ORIGIN_X="$2"; shift 2 ;;
    --grid_origin_y) GRID_ORIGIN_Y="$2"; shift 2 ;;
    --frequency) FREQUENCY="$2"; shift 2 ;;
    --pol_pair) POL_PAIR="$2"; shift 2 ;;
    --normalization) NORMALIZATION="$2"; shift 2 ;;
    --window) WINDOW="$2"; shift 2 ;;
    --noise_floor) NOISE_FLOOR="$2"; shift 2 ;;
    --min_coverage) MIN_COVERAGE="$2"; shift 2 ;;
    --use_exception_mask) USE_EXCEPTION_MASK="$2"; shift 2 ;;
    --lia) LIA="$2"; shift 2 ;;
    --lia_dem_scale) LIA_DEM_SCALE="$2"; shift 2 ;;
    --nisar_dem_vrt) NISAR_DEM_VRT="$2"; shift 2 ;;
    --dem_cache_dir) DEM_CACHE_DIR="$2"; shift 2 ;;
    --earthdata_cookie) EARTHDATA_COOKIE="$2"; shift 2 ;;
    --gap_threshold) GAP_THRESHOLD="$2"; shift 2 ;;
    --bbox) BBOX="$2"; shift 2 ;;
    --keep_intermediate) KEEP_INTERMEDIATE="$2"; shift 2 ;;
    --script) SCRIPT_PATH="$2"; shift 2 ;;
    --python_exec) PYTHON_EXEC="$2"; shift 2 ;;
    --extra_args) EXTRA_ARGS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$INPUT_DIR" ]] || { echo "--input_dir is required" >&2; exit 1; }
[[ -n "$OUTPUT_DIR" ]] || { echo "--output_dir is required" >&2; exit 1; }
[[ -d "$INPUT_DIR" ]] || { echo "Input directory not found: $INPUT_DIR" >&2; exit 1; }
[[ -f "$SCRIPT_PATH" ]] || { echo "Python script not found: $SCRIPT_PATH" >&2; exit 1; }
command -v "$PYTHON_EXEC" >/dev/null 2>&1 || { echo "Python executable not found: $PYTHON_EXEC" >&2; exit 1; }

if [[ -z "$DEM_CACHE_DIR" ]]; then
  DEM_CACHE_DIR="${INPUT_DIR}/DEM_Cache"
fi

case "$GRID_MODE" in UTM|EASE2) ;; *) echo "Invalid --grid_mode: $GRID_MODE (use UTM or EASE2)" >&2; exit 1 ;; esac
case "$MODE" in nofilt|enhanced_lee|gamma_map|all) ;; *) echo "Invalid --mode: $MODE" >&2; exit 1 ;; esac
case "$MOSAIC" in yes|no) ;; *) echo "Invalid --mosaic: $MOSAIC" >&2; exit 1 ;; esac
case "$FREQUENCY" in frequencyA|frequencyB) ;; *) echo "Invalid --frequency: $FREQUENCY" >&2; exit 1 ;; esac
case "$POL_PAIR" in auto|HH_HV|VV_VH|RH_RV) ;; *) echo "Invalid --pol_pair: $POL_PAIR" >&2; exit 1 ;; esac
case "$NORMALIZATION" in gamma0|sigma0) ;; *) echo "Invalid --normalization: $NORMALIZATION" >&2; exit 1 ;; esac
case "$KEEP_INTERMEDIATE" in yes|no) ;; *) echo "Invalid --keep_intermediate: $KEEP_INTERMEDIATE" >&2; exit 1 ;; esac
case "$USE_EXCEPTION_MASK" in yes|no) ;; *) echo "Invalid --use_exception_mask: $USE_EXCEPTION_MASK" >&2; exit 1 ;; esac
case "$LIA" in yes|no) ;; *) echo "Invalid --lia: $LIA" >&2; exit 1 ;; esac
if [[ "$GRID_MODE" == "UTM" ]]; then
  if [[ "$TARGET_EPSG" != "auto" && "$TARGET_EPSG" != "AUTO" ]]; then
    [[ "$TARGET_EPSG" =~ ^[0-9]+$ ]] || {
      echo "--target_epsg must be 'auto' or a numeric EPSG code" >&2
      exit 1
    }
  fi
else
  if [[ "$TARGET_EPSG" != "auto" && "$TARGET_EPSG" != "AUTO" && "$TARGET_EPSG" != "6933" ]]; then
    echo "EASE2 mode fixes the CRS to EPSG:6933; remove --target_epsg $TARGET_EPSG" >&2
    exit 1
  fi
fi

if [[ -n "$EASE2_MAT" && ! -f "$EASE2_MAT" ]]; then
  echo "--ease2_mat file not found: $EASE2_MAT" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"

launch_job() {
  local label="$1"
  local filter_name="$2"
  local destination="$3"
  mkdir -p "$destination"
  local log_file="${destination}/launcher_${label}_${TIMESTAMP}.log"

  local cmd=(
    "$PYTHON_EXEC" "$SCRIPT_PATH"
    --input_dir "$INPUT_DIR"
    --output_dir "$destination"
    --grid_mode "$GRID_MODE"
    --target_epsg "$TARGET_EPSG"
    --scale "$SCALE"
    --grid_origin_x "$GRID_ORIGIN_X"
    --grid_origin_y "$GRID_ORIGIN_Y"
    --frequency "$FREQUENCY"
    --pol_pair "$POL_PAIR"
    --normalization "$NORMALIZATION"
    --filter "$filter_name"
    --filter_window "$WINDOW"
    --noise_floor "$NOISE_FLOOR"
    --min_coverage "$MIN_COVERAGE"
    --use_exception_mask "$USE_EXCEPTION_MASK"
    --lia "$LIA"
    --lia_dem_scale "$LIA_DEM_SCALE"
    --dem_cache_dir "$DEM_CACHE_DIR"
    --earthdata_cookie "$EARTHDATA_COOKIE"
    --mosaic "$MOSAIC"
    --gap_threshold "$GAP_THRESHOLD"
    --workers "$WORKERS"
    --tile_workers "$TILE_WORKERS"
    --tile_size_km "$TILE_SIZE_KM"
    --native_block_px "$NATIVE_BLOCK_PX"
    --keep_intermediate "$KEEP_INTERMEDIATE"
  )

  if [[ -n "$EASE2_MAT" ]]; then
    cmd+=(--ease2_mat "$EASE2_MAT")
  fi

  if [[ -n "$NISAR_DEM_VRT" ]]; then
    cmd+=(--nisar_dem_vrt "$NISAR_DEM_VRT")
  fi

  if [[ -n "$BBOX" ]]; then
    read -r -a bbox_values <<< "$BBOX"
    [[ ${#bbox_values[@]} -eq 4 ]] || {
      echo "--bbox requires four values: W S E N" >&2
      exit 1
    }
    cmd+=(--bbox "${bbox_values[@]}")
  fi

  if [[ -n "$EXTRA_ARGS" ]]; then
    read -r -a extra_values <<< "$EXTRA_ARGS"
    cmd+=("${extra_values[@]}")
  fi

  nohup "${cmd[@]}" >"$log_file" 2>&1 &
  local pid=$!

  echo "Launched $label"
  echo "  PID             : $pid"
  echo "  Grid mode       : $GRID_MODE"
  if [[ "$GRID_MODE" == "EASE2" ]]; then
    echo "  Grid            : EASE2 M0p1 / EPSG:6933 / 100.089502334956 m"
    if [[ -n "$EASE2_MAT" ]]; then
      echo "  EASE2 MAT check : $EASE2_MAT"
    fi
  else
    echo "  EPSG selection  : $TARGET_EPSG"
    echo "  Pixel/grid      : ${SCALE} m, origin (${GRID_ORIGIN_X}, ${GRID_ORIGIN_Y})"
  fi
  echo "  Mosaic          : $MOSAIC"
  echo "  Exception mask  : $USE_EXCEPTION_MASK"
  echo "  Earthdata cookie: $EARTHDATA_COOKIE"
  echo "  LIA/incidence   : $LIA"
  if [[ "$LIA" == "yes" ]]; then
    echo "  LIA DEM scale   : ${LIA_DEM_SCALE} m"
    echo "  DEM cache       : $DEM_CACHE_DIR"
  fi
  echo "  Output          : $destination"
  echo "  Launcher log    : $log_file"
  echo "  EPSG/result log : tail -f '$log_file'"
  echo
}

if [[ "$MODE" == "all" ]]; then
  echo "WARNING: --mode all starts three independent processors."
  echo "For a scientific baseline, run --mode nofilt first."
  launch_job nofilt none "${OUTPUT_DIR}/nofilt"
  launch_job enhanced_lee enhanced_lee "${OUTPUT_DIR}/enhanced_lee"
  launch_job gamma_map gamma_map "${OUTPUT_DIR}/gamma_map"
else
  case "$MODE" in
    nofilt) FILTER_NAME="none" ;;
    enhanced_lee) FILTER_NAME="enhanced_lee" ;;
    gamma_map) FILTER_NAME="gamma_map" ;;
  esac
  launch_job "$MODE" "$FILTER_NAME" "$OUTPUT_DIR"
fi