#!/usr/bin/env python3
"""
Generic NISAR GCOV HDF5 -> fixed-grid GeoTIFF mosaics.

This version is safe for incremental processing in a directory that already
contains outputs from earlier runs.

Resume / no-overwrite behavior
------------------------------
For each expected acquisition/cluster output set:

* If the expected SAR GeoTIFFs already exist, the acquisition is SKIPPED.
* If only some expected GeoTIFFs exist, the acquisition is SKIPPED with a
  warning. Existing files are never overwritten.
* If all expected GeoTIFFs exist but metadata/grid validation reports a
  mismatch, the acquisition is still SKIPPED with a warning. Existing files
  are never overwritten.
* Only an acquisition for which none of the three expected GeoTIFFs exists is
  processed.

This makes it possible to add new NISAR HDF5 files to the input directory and
rerun the same command without recomputing existing outputs.

Destination grid modes
----------------------
--grid_mode UTM keeps the original fixed projected-grid workflow. By default,
--target_epsg auto reads the native EPSG from every selected NISAR HDF5 file
and chooses the majority EPSG.

--grid_mode EASE2 uses the exact global EASE-Grid 2.0 M0p1 lattice:
EPSG:6933, 100.089502334956 m, 347040 columns x 146160 rows. Regional outputs
are cropped from that global lattice by integer global row/column indices, so
all frames, tracks, dates, and regions are exactly aligned.

Processing architecture
-----------------------
1. Stream each GCOV frame into temporary native-grid GeoTIFFs while applying
   the product-defined validity mask:
       mask = 1..numberOfSubSwaths -> valid, fully focused samples
       mask = 0                    -> invalid/partially focused contributors
       mask = 255                  -> outside acquisition extent
   The processor also requires numberOfLooks > 0, finite covariance values,
   and no inputDataExceptionMask flags when that layer is available.
2. Open each complete native-grid raster as one WarpedVRT on the canonical
   destination grid. Processing-window size changes memory and speed, not the
   output pixel values.
3. Aggregate to the requested resolution in linear power and require a minimum
   valid source-area fraction per output pixel.
4. Mosaic adjacent compatible frames using a true per-pixel sum/count mean.
5. Save co-pol, cross-pol, and source-frame-count GeoTIFFs. In EASE2 mode,
   also save one compact *_EASE2_grid.mat companion containing the crop's global
   EASE2 row/column IDs and x/y/lon/lat pixel-center vectors.
6. Quicklooks display SAR in dB (10*log10 of linear output power).
7. With --lia yes, persist/reuse per-frame DEM subsets in DEM_Cache, then
   compute incidence/LIA from the NISAR radarGrid and the
   NISAR-modified Copernicus DEM at DEM terrain scale, aggregate them to the
   same destination grid, and mask them to the same valid SAR support.

Different dates, tracks, passes, modes, polarizations, or processor releases
are never averaged together. Single-polarization groups are skipped when a
co-pol/cross-pol pair is requested.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import fcntl
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import ExitStack
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional

import h5py
import numpy as np
import rasterio
from rasterio.transform import Affine, array_bounds
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window
from rasterio.warp import Resampling, transform_bounds
from scipy.ndimage import map_coordinates, uniform_filter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


NODATA = np.float32(-9999.0)

# ---------------------------------------------------------------------------
# EASE-Grid 2.0 Global Cylindrical Equal-Area, M0p1 (~100 m) grid.
#
# These constants match the user's smapease2forward.m definition:
#   gridid     = M0p1
#   EPSG       = 6933
#   resolution = 100.089502334956 m
#   columns    = 347040
#   rows       = 146160
#
# MATLAB/NSIDC convention uses zero-based row/column at pixel CENTERS:
#   col = 0 is the first pixel center, with its left edge at -0.5
#   row = 0 is the first pixel center, with its top edge at -0.5
#
# Therefore the global raster outer edges are exactly half the dimensions
# times the M0p1 spacing. Regional outputs are always cropped on these exact
# global cell boundaries; every output pixel can therefore be mapped back to
# one unique global M0p1 row/column.
# ---------------------------------------------------------------------------
EASE2_M0P1_EPSG = 6933
EASE2_M0P1_RESOLUTION_M = 100.089502334956
EASE2_M0P1_COLS = 347040
EASE2_M0P1_ROWS = 146160
EASE2_M0P1_LEFT = -(EASE2_M0P1_COLS / 2.0) * EASE2_M0P1_RESOLUTION_M
EASE2_M0P1_RIGHT = +(EASE2_M0P1_COLS / 2.0) * EASE2_M0P1_RESOLUTION_M
EASE2_M0P1_TOP = +(EASE2_M0P1_ROWS / 2.0) * EASE2_M0P1_RESOLUTION_M
EASE2_M0P1_BOTTOM = -(EASE2_M0P1_ROWS / 2.0) * EASE2_M0P1_RESOLUTION_M

DEFAULT_NISAR_DEM_VRT = (
    "https://nisar.asf.earthdatacloud.nasa.gov/"
    "NISAR/DEM/v1.2/EPSG4326/EPSG4326.vrt"
)

POL_PAIRS = {
    "HH_HV": ("HHHH", "HVHV"),
    "VV_VH": ("VVVV", "VHVH"),
    "RH_RV": ("RHRH", "RVRV"),
}

_FN_PATTERN = re.compile(
    r"NISAR_([LS])(\d)_([A-Z]{2})_GCOV_"
    r"(\d{3})_(\d{3})_([AD])_(\d{3})_(\d{4})_([A-Z]{4})_([AM])_"
    r"(\d{8})T(\d{6})_(\d{8})T(\d{6})_"
    r"([A-Z0-9]{6})_([PMNF])_([FP])_([A-Z])_(\d{3})"
    r"\.(?:h5|hdf5)$",
    re.IGNORECASE,
)


def setup_logger(log_path: Optional[str] = None) -> logging.Logger:
    logger = logging.getLogger(f"NISAR_GCOV_{os.getpid()}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    stream = logging.StreamHandler(sys.stdout)
    stream.setFormatter(formatter)
    logger.addHandler(stream)

    if log_path:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def parse_filename(path: str) -> Optional[dict]:
    match = _FN_PATTERN.fullmatch(os.path.basename(path))
    if not match:
        return None

    return {
        "instrument": match.group(1).upper(),
        "level": match.group(2),
        "proc_type": match.group(3).upper(),
        "cycle": match.group(4),
        "track": match.group(5),
        "pass": match.group(6).upper(),
        "frame": match.group(7),
        "bw_mode": match.group(8),
        "pole": match.group(9).upper(),
        "source": match.group(10).upper(),
        "date": match.group(11),
        "time": match.group(12),
        "end_date": match.group(13),
        "end_time": match.group(14),
        "crid": match.group(15).upper(),
        "accuracy": match.group(16).upper(),
        "coverage": match.group(17).upper(),
        "loc": match.group(18).upper(),
        "counter": match.group(19),
        "stem": os.path.basename(path).rsplit(".", 1)[0],
    }


def sar_root(info: dict) -> str:
    return "LSAR" if info["instrument"] == "L" else "SSAR"


def read_nisar_gcov_epsg(
    path: str,
    frequency: str = "frequencyA",
    pol_pair: Optional[str] = None,
) -> int:
    """Read source EPSG and optionally require a compatible polarization pair."""
    info = parse_filename(path)
    if info is None:
        raise ValueError(f"Unrecognized NISAR GCOV filename: {path}")

    root = sar_root(info)
    projection_path = f"/science/{root}/GCOV/grids/{frequency}/projection"
    group_path = f"/science/{root}/GCOV/grids/{frequency}"

    with h5py.File(path, "r") as h5:
        if group_path not in h5:
            raise KeyError(f"Missing GCOV grid group: {group_path}")

        group = h5[group_path]

        # Count only files that can actually be processed for the requested
        # polarization pair.
        if pol_pair is not None:
            select_pol_pair(group, pol_pair)

        if projection_path not in h5:
            raise KeyError(f"Missing projection dataset: {projection_path}")

        dataset = h5[projection_path]
        epsg = int(np.asarray(dataset[()]).item())

        if epsg <= 0 and "epsg_code" in dataset.attrs:
            epsg = int(np.asarray(dataset.attrs["epsg_code"]).item())

    if epsg <= 0:
        raise ValueError(f"Invalid EPSG code {epsg} in {path}")

    return epsg


def resolve_target_epsg(
    paths: list[str],
    frequency: str,
    pol_pair: str,
    requested: str,
    log: logging.Logger,
) -> int:
    """Resolve destination EPSG from an explicit value or the source majority."""
    requested_text = str(requested).strip()

    if requested_text.lower() != "auto":
        try:
            epsg = int(requested_text)
        except ValueError as error:
            raise ValueError(
                "--target_epsg must be 'auto' or a positive numeric EPSG code"
            ) from error

        if epsg <= 0:
            raise ValueError("--target_epsg must be positive")

        log.info(f"Target EPSG explicitly supplied: EPSG:{epsg}")
        return epsg

    log.info("Resolving target EPSG automatically from selected NISAR HDF5 files")
    records: list[tuple[str, int]] = []

    for path in sorted(paths):
        try:
            epsg = read_nisar_gcov_epsg(
                path,
                frequency=frequency,
                pol_pair=pol_pair,
            )
        except Exception as error:
            log.warning(
                "Skipping source EPSG candidate %s: %s",
                os.path.basename(path),
                error,
            )
            continue

        records.append((path, epsg))
        log.info("  Source CRS EPSG:%d <- %s", epsg, os.path.basename(path))

    if not records:
        raise RuntimeError("Could not read an EPSG code from any selected input file")

    counts = Counter(epsg for _, epsg in records)
    log.info("Source EPSG count summary:")
    for epsg, count in counts.most_common():
        log.info("  EPSG:%d -> %d file(s)", epsg, count)

    highest_count = max(counts.values())
    tied_majorities = {
        epsg for epsg, count in counts.items() if count == highest_count
    }

    if len(tied_majorities) == 1:
        chosen = next(iter(tied_majorities))
        log.info(
            "Automatically selected majority target CRS: EPSG:%d (%d file(s))",
            chosen,
            highest_count,
        )
        return chosen

    first_path, chosen = next(
        (path, epsg) for path, epsg in records if epsg in tied_majorities
    )
    tied_text = ", ".join(f"EPSG:{epsg}" for epsg in sorted(tied_majorities))
    log.warning(
        "Source EPSG tie among %s (%d file(s) each). Using EPSG:%d from first "
        "readable tied file: %s",
        tied_text,
        highest_count,
        chosen,
        os.path.basename(first_path),
    )
    return chosen


def compress_frames(frames: Iterable[str]) -> str:
    values = sorted({int(value) for value in frames})
    if not values:
        return "none"

    segments: list[str] = []
    start = end = values[0]

    for value in values[1:]:
        if value == end + 1:
            end = value
        else:
            segments.append(
                f"{start:03d}" if start == end else f"{start:03d}-{end:03d}"
            )
            start = end = value

    segments.append(
        f"{start:03d}" if start == end else f"{start:03d}-{end:03d}"
    )

    return "_".join(segments)


def select_pol_pair(group: h5py.Group, requested: str) -> tuple[str, str]:
    if requested != "auto":
        co_key, cross_key = POL_PAIRS[requested]
        if co_key not in group or cross_key not in group:
            raise KeyError(
                f"Requested {requested}, requiring {co_key}/{cross_key}; "
                f"available datasets: {sorted(group.keys())}"
            )
        return co_key, cross_key

    for name in ("HH_HV", "VV_VH", "RH_RV"):
        co_key, cross_key = POL_PAIRS[name]
        if co_key in group and cross_key in group:
            return co_key, cross_key

    raise KeyError(
        "No supported diagonal co/cross pair found. "
        f"Expected one of {list(POL_PAIRS.values())}; available={sorted(group.keys())}"
    )


@dataclass(frozen=True)
class FrameMeta:
    path: str
    info: dict
    root: str
    frequency: str
    co_key: str
    cross_key: str
    source_epsg: int
    width: int
    height: int
    transform: Affine
    source_bounds: tuple[float, float, float, float]
    target_bounds: tuple[float, float, float, float]


@dataclass(frozen=True)
class PreparedFrame:
    meta: FrameMeta
    co_path: str
    cross_path: str
    valid_path: str


@dataclass(frozen=True)
class PreparedGeometry:
    incidence_path: str
    lia_path: str


def inspect_frame(
    path: str,
    frequency: str,
    pol_pair: str,
    target_epsg: int,
) -> FrameMeta:
    info = parse_filename(path)
    if info is None:
        raise ValueError(f"Unrecognized NISAR GCOV filename: {path}")

    root = sar_root(info)
    group_path = f"/science/{root}/GCOV/grids/{frequency}"

    with h5py.File(path, "r") as h5:
        if group_path not in h5:
            raise KeyError(f"Missing HDF5 group: {group_path}")

        group = h5[group_path]

        required = [
            "xCoordinates",
            "yCoordinates",
            "projection",
            "mask",
            "numberOfLooks",
            "numberOfSubSwaths",
        ]
        missing = [name for name in required if name not in group]
        if missing:
            raise KeyError(f"Missing required datasets in {path}: {missing}")

        co_key, cross_key = select_pol_pair(group, pol_pair)
        x = np.asarray(group["xCoordinates"][:], dtype=np.float64)
        y = np.asarray(group["yCoordinates"][:], dtype=np.float64)
        epsg = int(np.asarray(group["projection"][()]).item())
        shape = tuple(group[co_key].shape)

        if len(shape) != 2:
            raise ValueError(f"{co_key} is not a 2-D raster in {path}")

        for dataset_name in (cross_key, "mask", "numberOfLooks"):
            if tuple(group[dataset_name].shape) != shape:
                raise ValueError(f"Shape mismatch for {dataset_name} in {path}")

    if len(x) != shape[1] or len(y) != shape[0] or len(x) < 2 or len(y) < 2:
        raise ValueError(f"Coordinate lengths do not match raster dimensions in {path}")

    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])

    if dx <= 0 or dy >= 0:
        raise ValueError(
            "Expected standard north-up GCOV orientation with increasing x and "
            f"decreasing y; found dx={dx}, dy={dy} in {path}"
        )

    atol_x = max(1e-6, abs(dx) * 1e-8)
    atol_y = max(1e-6, abs(dy) * 1e-8)

    if not np.allclose(np.diff(x), dx, rtol=0.0, atol=atol_x):
        raise ValueError(f"Irregular xCoordinates in {path}")

    if not np.allclose(np.diff(y), dy, rtol=0.0, atol=atol_y):
        raise ValueError(f"Irregular yCoordinates in {path}")

    transform = Affine(
        dx,
        0.0,
        x[0] - dx / 2.0,
        0.0,
        dy,
        y[0] - dy / 2.0,
    )

    source_bounds = array_bounds(shape[0], shape[1], transform)

    target_bounds = transform_bounds(
        rasterio.CRS.from_epsg(epsg),
        rasterio.CRS.from_epsg(target_epsg),
        *source_bounds,
        densify_pts=41,
    )

    return FrameMeta(
        path=path,
        info=info,
        root=root,
        frequency=frequency,
        co_key=co_key,
        cross_key=cross_key,
        source_epsg=epsg,
        width=shape[1],
        height=shape[0],
        transform=transform,
        source_bounds=source_bounds,
        target_bounds=target_bounds,
    )


def intersects(
    a: tuple[float, float, float, float],
    b: tuple[float, float, float, float],
) -> bool:
    return not (
        a[2] <= b[0]
        or b[2] <= a[0]
        or a[3] <= b[1]
        or b[3] <= a[1]
    )


def cluster_frames(
    frames: list[FrameMeta],
    gap_m: float,
    mosaic: bool,
) -> list[list[FrameMeta]]:
    if not mosaic:
        return [[frame] for frame in frames]

    parent = list(range(len(frames)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(i: int, j: int) -> None:
        root_i = find(i)
        root_j = find(j)
        if root_i != root_j:
            parent[root_j] = root_i

    def adjacent(a, b) -> bool:
        x_gap = max(0.0, max(a[0], b[0]) - min(a[2], b[2]))
        y_gap = max(0.0, max(a[1], b[1]) - min(a[3], b[3]))
        return x_gap <= gap_m and y_gap <= gap_m

    for i in range(len(frames)):
        for j in range(i + 1, len(frames)):
            if adjacent(frames[i].target_bounds, frames[j].target_bounds):
                union(i, j)

    grouped: dict[int, list[FrameMeta]] = defaultdict(list)
    for index, frame in enumerate(frames):
        grouped[find(index)].append(frame)

    return [
        sorted(cluster, key=lambda item: int(item.info["frame"]))
        for cluster in grouped.values()
    ]


def snap_floor(value: float, spacing: float, origin: float) -> float:
    return math.floor((value - origin) / spacing) * spacing + origin


def snap_ceil(value: float, spacing: float, origin: float) -> float:
    return math.ceil((value - origin) / spacing) * spacing + origin


def _cluster_bounds(
    frames: list[FrameMeta],
    bbox_wgs84: Optional[tuple[float, float, float, float]],
    target_epsg: int,
) -> tuple[float, float, float, float]:
    """Union frame bounds in the destination CRS and optionally intersect bbox."""
    left = min(frame.target_bounds[0] for frame in frames)
    bottom = min(frame.target_bounds[1] for frame in frames)
    right = max(frame.target_bounds[2] for frame in frames)
    top = max(frame.target_bounds[3] for frame in frames)

    if bbox_wgs84 is not None:
        clip_bounds = transform_bounds(
            "EPSG:4326",
            rasterio.CRS.from_epsg(target_epsg),
            *bbox_wgs84,
            densify_pts=41,
        )
        left = max(left, clip_bounds[0])
        bottom = max(bottom, clip_bounds[1])
        right = min(right, clip_bounds[2])
        top = min(top, clip_bounds[3])

        if left >= right or bottom >= top:
            raise ValueError("Requested --bbox does not intersect this cluster")

    return left, bottom, right, top


def _build_utm_target_grid(
    frames: list[FrameMeta],
    resolution: float,
    origin_x: float,
    origin_y: float,
    bbox_wgs84: Optional[tuple[float, float, float, float]],
    target_epsg: int,
) -> tuple[
    Affine,
    int,
    int,
    tuple[float, float, float, float],
    dict[str, str],
]:
    """
    Build the existing user-defined metric fixed grid.

    Despite the CLI label UTM, an explicit projected metre-based EPSG is still
    allowed, preserving backward compatibility. With --target_epsg auto this
    resolves to the majority native NISAR projected EPSG.
    """
    left, bottom, right, top = _cluster_bounds(
        frames,
        bbox_wgs84,
        target_epsg,
    )

    left = snap_floor(left, resolution, origin_x)
    right = snap_ceil(right, resolution, origin_x)
    bottom = snap_floor(bottom, resolution, origin_y)
    top = snap_ceil(top, resolution, origin_y)

    width = int(round((right - left) / resolution))
    height = int(round((top - bottom) / resolution))

    if width <= 0 or height <= 0:
        raise ValueError("Invalid fixed destination grid")

    transform = Affine(
        resolution,
        0.0,
        left,
        0.0,
        -resolution,
        top,
    )

    metadata = {
        "grid_mode": "UTM",
        "grid_name": "user_fixed_projected_grid",
    }
    return transform, width, height, (left, bottom, right, top), metadata


def _build_ease2_m0p1_target_grid(
    frames: list[FrameMeta],
    bbox_wgs84: Optional[tuple[float, float, float, float]],
) -> tuple[
    Affine,
    int,
    int,
    tuple[float, float, float, float],
    dict[str, str],
]:
    """
    Crop the exact global EASE-Grid 2.0 M0p1 lattice to the cluster.

    Alignment is integer-index based rather than generic floating-point
    snapping. This is the important part: all frames, tracks, dates, and
    regions use the same global cell boundaries.

    Global cell-center convention matches smapease2forward.m:
        col_global = 0..347039
        row_global = 0..146159
    """
    left, bottom, right, top = _cluster_bounds(
        frames,
        bbox_wgs84,
        EASE2_M0P1_EPSG,
    )

    # Intersect with the finite global M0p1 raster domain.
    left = max(left, EASE2_M0P1_LEFT)
    right = min(right, EASE2_M0P1_RIGHT)
    bottom = max(bottom, EASE2_M0P1_BOTTOM)
    top = min(top, EASE2_M0P1_TOP)

    if left >= right or bottom >= top:
        raise ValueError(
            "Cluster does not intersect the global EASE2 M0p1 domain"
        )

    res = EASE2_M0P1_RESOLUTION_M

    # A tiny index-space tolerance prevents a transformed bound that is
    # numerically 1e-12 beyond an exact grid edge from adding a spurious cell.
    eps = 1.0e-10

    col0 = math.floor((left - EASE2_M0P1_LEFT) / res + eps)
    col1 = math.ceil((right - EASE2_M0P1_LEFT) / res - eps)
    row0 = math.floor((EASE2_M0P1_TOP - top) / res + eps)
    row1 = math.ceil((EASE2_M0P1_TOP - bottom) / res - eps)

    col0 = max(0, min(col0, EASE2_M0P1_COLS))
    col1 = max(0, min(col1, EASE2_M0P1_COLS))
    row0 = max(0, min(row0, EASE2_M0P1_ROWS))
    row1 = max(0, min(row1, EASE2_M0P1_ROWS))

    if col1 <= col0 or row1 <= row0:
        raise ValueError("Invalid EASE2 M0p1 destination grid after clipping")

    width = col1 - col0
    height = row1 - row0

    snapped_left = EASE2_M0P1_LEFT + col0 * res
    snapped_right = EASE2_M0P1_LEFT + col1 * res
    snapped_top = EASE2_M0P1_TOP - row0 * res
    snapped_bottom = EASE2_M0P1_TOP - row1 * res

    transform = Affine(
        res,
        0.0,
        snapped_left,
        0.0,
        -res,
        snapped_top,
    )

    metadata = {
        "grid_mode": "EASE2",
        "grid_name": "EASE-Grid_2.0_M0p1",
        "ease2_grid_id": "M0p1",
        "ease2_global_rows": str(EASE2_M0P1_ROWS),
        "ease2_global_cols": str(EASE2_M0P1_COLS),
        "ease2_global_row_start": str(row0),
        "ease2_global_row_end_exclusive": str(row1),
        "ease2_global_col_start": str(col0),
        "ease2_global_col_end_exclusive": str(col1),
        "ease2_center_indexing": "zero_based_global_row_col",
    }

    return (
        transform,
        width,
        height,
        (snapped_left, snapped_bottom, snapped_right, snapped_top),
        metadata,
    )


def build_target_grid(
    frames: list[FrameMeta],
    resolution: float,
    origin_x: float,
    origin_y: float,
    bbox_wgs84: Optional[tuple[float, float, float, float]],
    target_epsg: int,
    grid_mode: str,
) -> tuple[
    Affine,
    int,
    int,
    tuple[float, float, float, float],
    dict[str, str],
]:
    grid_mode = grid_mode.upper()

    if grid_mode == "EASE2":
        if target_epsg != EASE2_M0P1_EPSG:
            raise ValueError(
                f"EASE2 mode requires EPSG:{EASE2_M0P1_EPSG}; "
                f"received EPSG:{target_epsg}"
            )
        return _build_ease2_m0p1_target_grid(
            frames,
            bbox_wgs84,
        )

    if grid_mode == "UTM":
        return _build_utm_target_grid(
            frames,
            resolution,
            origin_x,
            origin_y,
            bbox_wgs84,
            target_epsg,
        )

    raise ValueError(f"Unsupported grid mode: {grid_mode}")


def validate_ease2_mat(
    path: str,
    log: logging.Logger,
) -> None:
    """
    Optional sanity check against EASE2_latlon100m.mat.

    The MAT file is not needed to perform the reprojection; the exact M0p1
    projected lattice above defines the raster unambiguously. This check simply
    confirms that a supplied ancillary lat/lon file uses the same pixel centers.
    """
    from scipy.io import loadmat
    from rasterio.warp import transform as warp_coordinates

    mat = loadmat(path)
    if "lat" not in mat or "lon" not in mat:
        raise ValueError(
            f"{path} must contain variables named 'lat' and 'lon'"
        )

    lat = np.asarray(mat["lat"], dtype=np.float64).reshape(-1)
    lon = np.asarray(mat["lon"], dtype=np.float64).reshape(-1)

    if lat.size != EASE2_M0P1_ROWS or lon.size != EASE2_M0P1_COLS:
        raise ValueError(
            "EASE2 MAT dimensions do not match M0p1: "
            f"lat={lat.size}, lon={lon.size}; expected "
            f"{EASE2_M0P1_ROWS}, {EASE2_M0P1_COLS}"
        )

    res = EASE2_M0P1_RESOLUTION_M
    x_first = EASE2_M0P1_LEFT + 0.5 * res
    x_last = EASE2_M0P1_RIGHT - 0.5 * res
    y_first = EASE2_M0P1_TOP - 0.5 * res
    y_last = EASE2_M0P1_BOTTOM + 0.5 * res

    expected_lon, _ = warp_coordinates(
        f"EPSG:{EASE2_M0P1_EPSG}",
        "EPSG:4326",
        [x_first, x_last],
        [0.0, 0.0],
    )
    _, expected_lat = warp_coordinates(
        f"EPSG:{EASE2_M0P1_EPSG}",
        "EPSG:4326",
        [0.0, 0.0],
        [y_first, y_last],
    )

    checks = [
        ("lon_first", lon[0], expected_lon[0]),
        ("lon_last", lon[-1], expected_lon[-1]),
        ("lat_first", lat[0], expected_lat[0]),
        ("lat_last", lat[-1], expected_lat[-1]),
    ]
    tolerance_deg = 2.0e-7

    failed = [
        (name, observed, expected)
        for name, observed, expected in checks
        if not np.isfinite(observed)
        or abs(observed - expected) > tolerance_deg
    ]
    if failed:
        details = "; ".join(
            f"{name}: observed={observed:.10f}, expected={expected:.10f}"
            for name, observed, expected in failed
        )
        raise ValueError(
            "Supplied EASE2 MAT does not match the built-in M0p1 centers: "
            + details
        )

    log.info(
        "EASE2 ancillary validation PASSED: %s matches M0p1 "
        "(%d rows x %d cols, %.12f m)",
        path,
        EASE2_M0P1_ROWS,
        EASE2_M0P1_COLS,
        EASE2_M0P1_RESOLUTION_M,
    )



def ensure_ease2_grid_companion(
    path: str,
    transform: Affine,
    width: int,
    height: int,
    grid_metadata: dict[str, str],
    log: logging.Logger,
) -> None:
    """Create or validate the compact MATLAB companion for one EASE2 crop.

    The companion stores one-dimensional GLOBAL EASE2 row/column vectors rather
    than redundant two-dimensional row/column rasters.  Together with the raster
    shape, these vectors uniquely identify every GeoTIFF pixel.  Both zero-based
    (NSIDC / smapease2forward convention) and one-based (MATLAB-friendly)
    indices are saved, along with EPSG:6933 x/y and WGS84 lon/lat pixel-center
    vectors.

    Existing companion files are never overwritten.  They are validated against
    the requested crop; a mismatch is treated as an error because silently using
    a wrong EASE2 index file would break cross-scene alignment.
    """
    from scipy.io import loadmat, savemat
    from rasterio.warp import transform as warp_coordinates

    if grid_metadata.get("grid_mode") != "EASE2":
        raise ValueError("EASE2 grid companion requested for a non-EASE2 grid")

    row0 = int(grid_metadata["ease2_global_row_start"])
    row1 = int(grid_metadata["ease2_global_row_end_exclusive"])
    col0 = int(grid_metadata["ease2_global_col_start"])
    col1 = int(grid_metadata["ease2_global_col_end_exclusive"])

    if row1 - row0 != height or col1 - col0 != width:
        raise ValueError(
            "EASE2 companion index extent does not match raster dimensions: "
            f"rows {row0}:{row1} -> {row1-row0}, height={height}; "
            f"cols {col0}:{col1} -> {col1-col0}, width={width}"
        )

    res = EASE2_M0P1_RESOLUTION_M
    tol = max(1e-8, res * 1e-10)
    expected_left = EASE2_M0P1_LEFT + col0 * res
    expected_top = EASE2_M0P1_TOP - row0 * res

    if (
        abs(transform.a - res) > tol
        or abs(transform.e + res) > tol
        or abs(transform.b) > tol
        or abs(transform.d) > tol
        or abs(transform.c - expected_left) > tol
        or abs(transform.f - expected_top) > tol
    ):
        raise ValueError(
            "Target affine transform is not aligned to the exact EASE2 M0p1 "
            "global lattice"
        )

    # Compact global index vectors. These are the authoritative row/column IDs.
    row_0 = np.arange(row0, row1, dtype=np.int32)[:, None]
    col_0 = np.arange(col0, col1, dtype=np.int32)[None, :]
    row_1 = row_0 + np.int32(1)
    col_1 = col_0 + np.int32(1)

    # EPSG:6933 pixel-center coordinates, preserving raster orientation:
    # x increases left -> right; y decreases top -> bottom.
    x = (
        transform.c
        + (np.arange(width, dtype=np.float64) + 0.5) * transform.a
    )[None, :]
    y = (
        transform.f
        + (np.arange(height, dtype=np.float64) + 0.5) * transform.e
    )[:, None]

    # EPSG:6933 is cylindrical, so lon is a function of x and lat of y.
    # Transform only 1-D center vectors instead of constructing an H x W mesh.
    lon_values, _ = warp_coordinates(
        f"EPSG:{EASE2_M0P1_EPSG}",
        "EPSG:4326",
        x.reshape(-1).tolist(),
        np.zeros(width, dtype=np.float64).tolist(),
    )
    _, lat_values = warp_coordinates(
        f"EPSG:{EASE2_M0P1_EPSG}",
        "EPSG:4326",
        np.zeros(height, dtype=np.float64).tolist(),
        y.reshape(-1).tolist(),
    )
    lon = np.asarray(lon_values, dtype=np.float64)[None, :]
    lat = np.asarray(lat_values, dtype=np.float64)[:, None]

    payload = {
        # Global EASE2 indices.
        "ease2_row_0based": row_0,
        "ease2_col_0based": col_0,
        "ease2_row_1based": row_1,
        "ease2_col_1based": col_1,
        # Short aliases use the NSIDC / smapease2forward zero-based convention.
        "ease2_row": row_0,
        "ease2_col": col_0,
        # Pixel-center coordinates.
        "x": x,
        "y": y,
        "lon": lon,
        "lat": lat,
        # Crop/index metadata.
        "row_start_0based": np.int64(row0),
        "row_end_inclusive_0based": np.int64(row1 - 1),
        "row_end_exclusive_0based": np.int64(row1),
        "col_start_0based": np.int64(col0),
        "col_end_inclusive_0based": np.int64(col1 - 1),
        "col_end_exclusive_0based": np.int64(col1),
        "n_rows": np.int64(height),
        "n_cols": np.int64(width),
        "global_n_rows": np.int64(EASE2_M0P1_ROWS),
        "global_n_cols": np.int64(EASE2_M0P1_COLS),
        "resolution_m": np.float64(res),
        "epsg": np.int32(EASE2_M0P1_EPSG),
        "grid_id": "M0p1",
        "grid_name": "EASE-Grid 2.0 Global Cylindrical Equal-Area M0p1",
        "index_convention": (
            "ease2_row/ease2_col are global zero-based pixel-center indices; "
            "*_1based variables are MATLAB-friendly global indices"
        ),
        "cell_id_formula_0based": (
            f"cell_id = ease2_row_0based * {EASE2_M0P1_COLS} + "
            "ease2_col_0based"
        ),
        # GDAL-style affine terms: x = GT0 + col*GT1 + row*GT2;
        # y = GT3 + col*GT4 + row*GT5. GT0/GT3 are OUTER pixel edges.
        "geotransform_gdal": np.asarray(
            [
                transform.c,
                transform.a,
                transform.b,
                transform.f,
                transform.d,
                transform.e,
            ],
            dtype=np.float64,
        ),
    }

    def scalar(mat: dict, name: str) -> float:
        if name not in mat:
            raise ValueError(f"missing variable {name!r}")
        arr = np.asarray(mat[name]).reshape(-1)
        if arr.size != 1:
            raise ValueError(f"variable {name!r} is not scalar")
        return float(arr[0])

    def validate_existing() -> None:
        # Load only the compact fields needed to prove that the companion belongs
        # to this exact crop.  We intentionally do not read the full coordinate
        # vectors just to resume an existing run.
        variables = [
            "row_start_0based",
            "row_end_exclusive_0based",
            "col_start_0based",
            "col_end_exclusive_0based",
            "n_rows",
            "n_cols",
            "resolution_m",
            "epsg",
            "ease2_row_0based",
            "ease2_col_0based",
        ]
        mat = loadmat(path, variable_names=variables)
        checks = {
            "row_start_0based": row0,
            "row_end_exclusive_0based": row1,
            "col_start_0based": col0,
            "col_end_exclusive_0based": col1,
            "n_rows": height,
            "n_cols": width,
            "epsg": EASE2_M0P1_EPSG,
        }
        for name, expected in checks.items():
            observed = scalar(mat, name)
            if int(round(observed)) != int(expected):
                raise ValueError(
                    f"{name}: observed={observed}, expected={expected}"
                )

        observed_res = scalar(mat, "resolution_m")
        if abs(observed_res - res) > 1e-9:
            raise ValueError(
                f"resolution_m: observed={observed_res}, expected={res}"
            )

        r = np.asarray(mat.get("ease2_row_0based", [])).reshape(-1)
        c = np.asarray(mat.get("ease2_col_0based", [])).reshape(-1)
        if (
            r.size != height
            or c.size != width
            or int(r[0]) != row0
            or int(r[-1]) != row1 - 1
            or int(c[0]) != col0
            or int(c[-1]) != col1 - 1
        ):
            raise ValueError("stored EASE2 row/column vectors do not match crop")

    path = os.path.abspath(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    lock_path = path + ".lock"
    lock_handle = open(lock_path, "a+")
    try:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        if os.path.isfile(path):
            try:
                validate_existing()
            except Exception as error:
                raise RuntimeError(
                    "Existing EASE2 companion MAT does not match the requested "
                    f"grid and will not be overwritten: {path}. Details: {error}"
                ) from error
            log.info("EASE2 grid companion verified: %s", path)
            return

        tmp_path = path + f".tmp.{os.getpid()}.mat"
        try:
            savemat(
                tmp_path,
                payload,
                appendmat=False,
                do_compression=True,
                oned_as="row",
            )
            os.replace(tmp_path, path)
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

        log.info(
            "Saved EASE2 grid companion: %s (rows %d:%d, cols %d:%d)",
            path,
            row0,
            row1,
            col0,
            col1,
        )
    finally:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
        lock_handle.close()


def iter_windows(width: int, height: int, block_px: int) -> list[Window]:
    windows: list[Window] = []

    for row_off in range(0, height, block_px):
        window_height = min(block_px, height - row_off)

        for col_off in range(0, width, block_px):
            window_width = min(block_px, width - col_off)
            windows.append(
                Window(
                    col_off,
                    row_off,
                    window_width,
                    window_height,
                )
            )

    return windows


def expand_window(
    window: Window,
    halo: int,
    width: int,
    height: int,
) -> Window:
    col0 = max(0, int(window.col_off) - halo)
    row0 = max(0, int(window.row_off) - halo)
    col1 = min(width, int(window.col_off + window.width) + halo)
    row1 = min(height, int(window.row_off + window.height) + halo)

    return Window(
        col0,
        row0,
        col1 - col0,
        row1 - row0,
    )


def _windowed_stats(
    array: np.ndarray,
    size: int,
) -> tuple[np.ndarray, np.ndarray]:
    valid = np.isfinite(array).astype(np.float32)
    filled = np.where(valid > 0, array, 0.0).astype(np.float32)

    mean_numerator = uniform_filter(
        filled,
        size=size,
        mode="reflect",
    )

    square_numerator = uniform_filter(
        filled * filled,
        size=size,
        mode="reflect",
    )

    count = uniform_filter(
        valid,
        size=size,
        mode="reflect",
    )

    mean = np.full(array.shape, np.nan, dtype=np.float32)
    variance = np.full(array.shape, np.nan, dtype=np.float32)

    good = count > 0
    mean[good] = mean_numerator[good] / count[good]

    second_moment = np.zeros(array.shape, dtype=np.float32)
    second_moment[good] = square_numerator[good] / count[good]

    variance[good] = np.maximum(
        second_moment[good] - mean[good] ** 2,
        0.0,
    )

    return mean, variance


def enhanced_lee_filter(
    array: np.ndarray,
    looks: np.ndarray,
    window: int,
) -> np.ndarray:
    local_mean, local_variance = _windowed_stats(array, window)

    enl = np.where(
        np.isfinite(looks),
        np.maximum(looks, 1.0),
        np.nan,
    ).astype(np.float32)

    with np.errstate(divide="ignore", invalid="ignore"):
        ci = np.sqrt(local_variance) / local_mean
        cu = np.sqrt(1.0 / enl)
        cmax = np.sqrt(1.0 + 2.0 / enl)

    output = array.astype(np.float32).copy()

    valid = (
        np.isfinite(array)
        & np.isfinite(local_mean)
        & np.isfinite(ci)
        & np.isfinite(enl)
        & (local_mean > 0)
    )

    homogeneous = valid & (ci <= cu)
    output[homogeneous] = local_mean[homogeneous]

    textured = valid & (ci > cu) & (ci < cmax)

    if np.any(textured):
        denominator = np.maximum(
            cmax[textured] - ci[textured],
            1e-8,
        )

        weight = np.exp(
            -(ci[textured] - cu[textured]) / denominator
        ).clip(0.0, 1.0)

        output[textured] = (
            weight * local_mean[textured]
            + (1.0 - weight) * array[textured]
        ).astype(np.float32)

    output[~np.isfinite(array)] = np.nan
    return output


def gamma_map_filter(
    array: np.ndarray,
    looks: np.ndarray,
    window: int,
) -> np.ndarray:
    local_mean, local_variance = _windowed_stats(array, window)

    enl = np.where(
        np.isfinite(looks),
        np.maximum(looks, 1.0001),
        np.nan,
    ).astype(np.float32)

    with np.errstate(divide="ignore", invalid="ignore"):
        ci = np.sqrt(local_variance) / local_mean
        cu = np.sqrt(1.0 / enl)
        cmax = np.sqrt(2.0) * cu

    output = array.astype(np.float32).copy()

    valid = (
        np.isfinite(array)
        & np.isfinite(local_mean)
        & np.isfinite(ci)
        & np.isfinite(enl)
        & (local_mean > 0)
    )

    homogeneous = valid & (ci <= cu)
    output[homogeneous] = local_mean[homogeneous]

    textured = valid & (ci > cu) & (ci < cmax)

    if np.any(textured):
        ci_t = ci[textured].astype(np.float64)
        mean_t = local_mean[textured].astype(np.float64)
        observation_t = array[textured].astype(np.float64)
        enl_t = enl[textured].astype(np.float64)
        cu_t = cu[textured].astype(np.float64)

        denominator = np.maximum(
            ci_t**2 - cu_t**2,
            1e-12,
        )

        alpha = (1.0 + cu_t**2) / denominator
        b = alpha - enl_t - 1.0

        discriminant = np.maximum(
            b**2 * mean_t**2
            + 4.0 * alpha * enl_t * mean_t * observation_t,
            0.0,
        )

        output[textured] = np.maximum(
            (b * mean_t + np.sqrt(discriminant)) / (2.0 * alpha),
            0.0,
        ).astype(np.float32)

    output[~np.isfinite(array)] = np.nan
    return output


def _fill_value(dataset: h5py.Dataset, default=None):
    """Return a scalar HDF5 fill value when available."""
    for name in ("_FillValue", "fill_value"):
        if name in dataset.attrs:
            value = np.asarray(dataset.attrs[name])
            if value.size == 1:
                return value.reshape(-1)[0].item()
    return default


def _not_fill(array: np.ndarray, fill_value) -> np.ndarray:
    if fill_value is None:
        return np.ones(array.shape, dtype=bool)

    try:
        if np.isnan(fill_value):
            return ~np.isnan(array)
    except TypeError:
        pass

    return array != fill_value


def validity_profile(frame: FrameMeta) -> dict:
    """Native binary validity raster: 1=valid fully focused GCOV, 0=invalid."""
    return {
        "driver": "GTiff",
        "width": frame.width,
        "height": frame.height,
        "count": 1,
        "dtype": "uint8",
        "crs": rasterio.CRS.from_epsg(frame.source_epsg),
        "transform": frame.transform,
        "nodata": None,
        "compress": "DEFLATE",
        "predictor": 1,
        "tiled": True,
        "blockxsize": 512,
        "blockysize": 512,
        "BIGTIFF": "IF_SAFER",
    }


def native_profile(frame: FrameMeta) -> dict:
    return {
        "driver": "GTiff",
        "width": frame.width,
        "height": frame.height,
        "count": 1,
        "dtype": "float32",
        "crs": rasterio.CRS.from_epsg(frame.source_epsg),
        "transform": frame.transform,
        "nodata": float(NODATA),
        "compress": "DEFLATE",
        "predictor": 2,
        "tiled": True,
        "blockxsize": 512,
        "blockysize": 512,
        "BIGTIFF": "IF_SAFER",
    }


def prepare_native_frame(
    frame: FrameMeta,
    work_dir: str,
    normalization: str,
    filter_name: str,
    filter_window: int,
    noise_floor: float,
    use_exception_mask: bool,
    native_block_px: int,
    log: logging.Logger,
) -> PreparedFrame:
    frame_label = f"frame_{frame.info['frame']}_{os.getpid()}"

    co_path = os.path.join(
        work_dir,
        f"{frame_label}_{frame.co_key}.tif",
    )

    cross_path = os.path.join(
        work_dir,
        f"{frame_label}_{frame.cross_key}.tif",
    )

    valid_path = os.path.join(
        work_dir,
        f"{frame_label}_valid_fraction_source.tif",
    )

    windows = iter_windows(
        frame.width,
        frame.height,
        native_block_px,
    )

    halo = filter_window // 2 if filter_name != "none" else 0
    group_path = f"/science/{frame.root}/GCOV/grids/{frame.frequency}"

    log.info(
        f"  Preparing native frame {frame.info['frame']}: "
        f"{frame.height}x{frame.width}, EPSG:{frame.source_epsg}, "
        f"blocks={len(windows)}"
    )

    with h5py.File(frame.path, "r") as h5, \
         rasterio.open(co_path, "w", **native_profile(frame)) as co_dst, \
         rasterio.open(cross_path, "w", **native_profile(frame)) as cross_dst, \
         rasterio.open(valid_path, "w", **validity_profile(frame)) as valid_dst:

        group = h5[group_path]
        co_dataset = group[frame.co_key]
        cross_dataset = group[frame.cross_key]
        looks_dataset = group["numberOfLooks"]
        mask_dataset = group["mask"]
        exception_dataset = group.get("inputDataExceptionMask")

        co_fill = _fill_value(co_dataset)
        cross_fill = _fill_value(cross_dataset)
        looks_fill = _fill_value(looks_dataset)

        n_subswaths = int(
            np.asarray(group["numberOfSubSwaths"][()]).item()
        )

        if n_subswaths < 1 or n_subswaths > 254:
            raise ValueError(
                f"Invalid numberOfSubSwaths={n_subswaths} in {frame.path}"
            )

        mask_fill = _fill_value(mask_dataset, 255)
        exception_fill = (
            _fill_value(exception_dataset)
            if exception_dataset is not None
            else None
        )

        if normalization == "sigma0" and "rtcGammaToSigmaFactor" not in group:
            raise KeyError(
                "sigma0 requested, but rtcGammaToSigmaFactor is missing in "
                f"{frame.path}"
            )

        report_step = max(1, len(windows) // 5)

        total_pixels = 0
        valid_pixels = 0
        mask_zero_pixels = 0
        mask_fill_pixels = 0
        mask_other_pixels = 0

        for index, output_window in enumerate(windows, start=1):
            read_window = expand_window(
                output_window,
                halo,
                frame.width,
                frame.height,
            )

            row_slice = slice(
                int(read_window.row_off),
                int(read_window.row_off + read_window.height),
            )

            col_slice = slice(
                int(read_window.col_off),
                int(read_window.col_off + read_window.width),
            )

            co = np.asarray(
                co_dataset[row_slice, col_slice],
                dtype=np.float32,
            )

            cross = np.asarray(
                cross_dataset[row_slice, col_slice],
                dtype=np.float32,
            )

            looks = np.asarray(
                looks_dataset[row_slice, col_slice],
                dtype=np.float32,
            )

            mask = np.asarray(
                mask_dataset[row_slice, col_slice],
                dtype=np.uint8,
            )

            valid_mask = (
                (mask >= 1)
                & (mask <= n_subswaths)
            )

            valid = (
                valid_mask
                & np.isfinite(looks)
                & _not_fill(looks, looks_fill)
                & (looks > 0.0)
                & np.isfinite(co)
                & np.isfinite(cross)
                & _not_fill(co, co_fill)
                & _not_fill(cross, cross_fill)
                & (co >= 0.0)
                & (cross >= 0.0)
            )

            if use_exception_mask and exception_dataset is not None:
                exception = np.asarray(
                    exception_dataset[row_slice, col_slice],
                    dtype=np.uint8,
                )

                valid &= (
                    _not_fill(exception, exception_fill)
                    & (exception == 0)
                )

            if normalization == "sigma0":
                rtc = np.asarray(
                    group["rtcGammaToSigmaFactor"][row_slice, col_slice],
                    dtype=np.float32,
                )

                valid &= np.isfinite(rtc) & (rtc > 0)
                co = co * rtc
                cross = cross * rtc

            coverage_valid = valid.copy()

            co[~valid] = np.nan
            cross[~valid] = np.nan
            looks[~valid] = np.nan

            if noise_floor > 0:
                co[np.isfinite(co) & (co < noise_floor)] = np.nan
                cross[np.isfinite(cross) & (cross < noise_floor)] = np.nan

            if filter_name == "enhanced_lee":
                co = enhanced_lee_filter(
                    co,
                    looks,
                    filter_window,
                )
                cross = enhanced_lee_filter(
                    cross,
                    looks,
                    filter_window,
                )

            elif filter_name == "gamma_map":
                co = gamma_map_filter(
                    co,
                    looks,
                    filter_window,
                )
                cross = gamma_map_filter(
                    cross,
                    looks,
                    filter_window,
                )

            row0 = int(output_window.row_off - read_window.row_off)
            col0 = int(output_window.col_off - read_window.col_off)
            row1 = row0 + int(output_window.height)
            col1 = col0 + int(output_window.width)

            co_out = co[row0:row1, col0:col1]
            cross_out = cross[row0:row1, col0:col1]
            coverage_out = coverage_valid[row0:row1, col0:col1]
            mask_out = mask[row0:row1, col0:col1]
            valid_mask_out = valid_mask[row0:row1, col0:col1]

            total_pixels += int(mask_out.size)
            valid_pixels += int(np.count_nonzero(coverage_out))
            mask_zero_pixels += int(np.count_nonzero(mask_out == 0))
            mask_fill_pixels += int(np.count_nonzero(mask_out == mask_fill))
            mask_other_pixels += int(
                np.count_nonzero(
                    ~valid_mask_out
                    & (mask_out != 0)
                    & (mask_out != mask_fill)
                )
            )

            final_valid = (
                np.isfinite(co_out)
                & np.isfinite(cross_out)
            )

            co_out = np.where(
                final_valid,
                co_out,
                NODATA,
            ).astype(np.float32)

            cross_out = np.where(
                final_valid,
                cross_out,
                NODATA,
            ).astype(np.float32)

            valid_out = coverage_out.astype(np.uint8)

            co_dst.write(co_out, 1, window=output_window)
            cross_dst.write(cross_out, 1, window=output_window)
            valid_dst.write(valid_out, 1, window=output_window)

            if index % report_step == 0 or index == len(windows):
                log.info(
                    f"    native blocks completed: {index}/{len(windows)}"
                )

        log.info(
            f"    frame {frame.info['frame']} validity: "
            f"valid={valid_pixels:,}/{total_pixels:,} "
            f"({100.0 * valid_pixels / max(total_pixels, 1):.2f}%), "
            f"mask0={mask_zero_pixels:,}, "
            f"maskFill={mask_fill_pixels:,}, "
            f"maskOther={mask_other_pixels:,}, "
            f"nSubSwaths={n_subswaths}"
        )

        common_tags = {
            "source_hdf5": os.path.basename(frame.path),
            "source_frame": frame.info["frame"],
            "normalization": normalization,
            "filter": filter_name,
            "noise_floor_linear": str(noise_floor),
            "number_of_subswaths": str(n_subswaths),
            "native_mask_rule": (
                "retain 1..numberOfSubSwaths; reject 0 and 255"
            ),
            "input_data_exception_mask_applied": (
                "yes" if use_exception_mask else "no"
            ),
        }

        co_dst.update_tags(
            **common_tags,
            covariance_term=frame.co_key,
        )

        cross_dst.update_tags(
            **common_tags,
            covariance_term=frame.cross_key,
        )

        valid_description = (
            "1 where mask is in 1..numberOfSubSwaths, numberOfLooks>0, "
            "and both requested covariance terms are finite; 0 otherwise"
        )
        if use_exception_mask:
            valid_description += (
                ". Pixels with a nonzero inputDataExceptionMask value are "
                "also excluded (0)."
            )
        else:
            valid_description += (
                ". inputDataExceptionMask was NOT applied by request; "
                "pixels flagged by the instrument (e.g. sample slip/shift) "
                "are retained if they otherwise pass all other checks."
            )

        valid_dst.update_tags(
            **common_tags,
            description=valid_description,
        )

    return PreparedFrame(
        meta=frame,
        co_path=co_path,
        cross_path=cross_path,
        valid_path=valid_path,
    )



# ============================================================================
# OPTIONAL INCIDENCE / LOCAL-INCIDENCE-ANGLE SUPPORT
# ============================================================================


def _gdal_dem_source(value: str) -> str:
    """Return a GDAL-readable source path for local, /vsi, or HTTP DEM input."""
    text = str(value).strip()
    if text.startswith(("/vsicurl/", "/vsis3/")):
        return text
    if text.startswith(("http://", "https://")):
        return f"/vsicurl/{text}"
    return text


def _require_gdalwarp() -> str:
    executable = shutil.which("gdalwarp")
    if executable is None:
        raise RuntimeError(
            "LIA requested but gdalwarp was not found. Install GDAL, e.g. "
            "'conda install -c conda-forge gdal'."
        )
    return executable


def _dem_cache_spec(
    frame: FrameMeta,
    dem_vrt: str,
    dem_scale_m: float,
) -> tuple[str, tuple[float, float, float, float], str]:
    """Return stable cache key, buffered bounds, and GDAL source for one frame."""
    left, bottom, right, top = frame.source_bounds
    buffer_m = max(2.0 * dem_scale_m, 60.0)
    bounds = (
        float(left - buffer_m),
        float(bottom - buffer_m),
        float(right + buffer_m),
        float(top + buffer_m),
    )
    source = _gdal_dem_source(dem_vrt)
    signature = (
        f"source={source}|epsg={frame.source_epsg}|"
        f"left={bounds[0]:.6f}|bottom={bounds[1]:.6f}|"
        f"right={bounds[2]:.6f}|top={bounds[3]:.6f}|"
        f"scale={dem_scale_m:.6f}"
    )
    digest = hashlib.sha1(signature.encode("utf-8")).hexdigest()[:16]
    return digest, bounds, source


def _valid_cached_dem(
    path: str,
    target_epsg: int,
    dem_scale_m: float,
    required_bounds: tuple[float, float, float, float],
) -> bool:
    """Quick structural validation for a persistent cached DEM GeoTIFF."""
    if not os.path.isfile(path) or os.path.getsize(path) <= 0:
        return False

    try:
        with rasterio.open(path) as src:
            if src.count != 1 or src.width <= 0 or src.height <= 0:
                return False
            if src.crs != rasterio.CRS.from_epsg(target_epsg):
                return False
            tol = max(1e-6, dem_scale_m * 1e-6)
            if abs(abs(src.transform.a) - dem_scale_m) > tol:
                return False
            if abs(abs(src.transform.e) - dem_scale_m) > tol:
                return False

            left, bottom, right, top = required_bounds
            cover_tol = max(1e-6, dem_scale_m * 1.01)
            if src.bounds.left > left + cover_tol:
                return False
            if src.bounds.bottom > bottom + cover_tol:
                return False
            if src.bounds.right < right - cover_tol:
                return False
            if src.bounds.top < top - cover_tol:
                return False
    except Exception:
        return False

    return True


def _prepare_dem_for_frame(
    frame: FrameMeta,
    work_dir: str,
    dem_vrt: str,
    dem_scale_m: float,
    warp_threads: int,
    log: logging.Logger,
    dem_cache_dir: Optional[str] = None,
    earthdata_cookie: Optional[str] = None,
) -> str:
    """
    Subset/reproject the NISAR-modified Copernicus DEM to a metric terrain
    grid covering one GCOV frame.

    If dem_cache_dir is supplied, the resulting DEM is stored persistently and
    reused on future runs whenever the DEM source, CRS, buffered frame bounds,
    and terrain scale are identical. A file lock prevents parallel workers
    from downloading/building the same cached DEM simultaneously.
    """
    gdalwarp = _require_gdalwarp()
    digest, required_bounds, source = _dem_cache_spec(
        frame=frame,
        dem_vrt=dem_vrt,
        dem_scale_m=dem_scale_m,
    )
    left, bottom, right, top = required_bounds

    # Each worker gets its own writable Earthdata/GDAL cookie jar.  If a
    # known-working master cookie is supplied, seed the worker jar from it so
    # parallel gdalwarp processes inherit the authenticated Earthdata session
    # without concurrently writing the same cookie file.
    cookie_path = f"/tmp/nisar_gdal_cookie_{os.getuid()}_{os.getpid()}.txt"
    master_cookie = (
        os.path.abspath(os.path.expanduser(earthdata_cookie))
        if earthdata_cookie
        else None
    )
    if master_cookie and os.path.isfile(master_cookie):
        shutil.copyfile(master_cookie, cookie_path)
        log.info(
            "  Seeded worker GDAL cookie from %s -> %s",
            master_cookie,
            cookie_path,
        )
    else:
        Path(cookie_path).touch(exist_ok=True)
        if master_cookie:
            log.warning(
                "  Earthdata master cookie not found (%s); using a fresh worker cookie",
                master_cookie,
            )
    try:
        os.chmod(cookie_path, 0o600)
    except OSError:
        pass

    cache_path: Optional[str] = None
    lock_handle = None

    if dem_cache_dir:
        cache_dir = os.path.abspath(os.path.expanduser(dem_cache_dir))
        os.makedirs(cache_dir, exist_ok=True)
        scale_tag = f"{dem_scale_m:g}".replace(".", "p")
        cache_name = (
            f"NISAR_DEM_EPSG{frame.source_epsg}_{scale_tag}m_{digest}.tif"
        )
        cache_path = os.path.join(cache_dir, cache_name)
        lock_path = cache_path + ".lock"

        # Keep the lock file itself; flock is released automatically. Leaving
        # the tiny lock file avoids inode-race problems between waiting workers.
        lock_handle = open(lock_path, "a+")
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)

        if _valid_cached_dem(
            cache_path,
            frame.source_epsg,
            dem_scale_m,
            required_bounds,
        ):
            log.info(
                "  DEM cache HIT for frame %s -> %s",
                frame.info["frame"],
                cache_path,
            )
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
            lock_handle.close()
            return cache_path

        if os.path.exists(cache_path):
            log.warning(
                "  Cached DEM failed validation; rebuilding: %s",
                cache_path,
            )
            try:
                os.remove(cache_path)
            except OSError:
                pass

        dem_path = cache_path + f".tmp.{os.getpid()}.tif"
        log.info(
            "  DEM cache MISS for frame %s; building persistent DEM -> %s",
            frame.info["frame"],
            cache_path,
        )
    else:
        frame_label = f"frame_{frame.info['frame']}_{os.getpid()}"
        dem_path = os.path.join(
            work_dir,
            f"{frame_label}_nisar_dem_{dem_scale_m:g}m.tif",
        )
        log.info(
            "  Preparing temporary DEM for frame %s at %.1f m terrain scale",
            frame.info["frame"],
            dem_scale_m,
        )

    command = [
        gdalwarp,
        "-overwrite",
        "-of", "GTiff",
        "-t_srs", f"EPSG:{frame.source_epsg}",
        "-te", str(left), str(bottom), str(right), str(top),
        "-tr", str(dem_scale_m), str(dem_scale_m),
        "-tap",
        "-r", "bilinear",
        "-ot", "Float32",
        "-dstnodata", str(float(NODATA)),
        "-co", "TILED=YES",
        "-co", "COMPRESS=DEFLATE",
        "-co", "PREDICTOR=3",
        "-co", "BIGTIFF=IF_SAFER",
        "-wo", f"NUM_THREADS={max(1, warp_threads)}",
        "--config", "GDAL_DISABLE_READDIR_ON_OPEN", "TRUE",
        "--config", "GDAL_HTTP_NETRC", "YES",
        "--config", "GDAL_HTTP_COOKIEFILE", cookie_path,
        "--config", "GDAL_HTTP_COOKIEJAR", cookie_path,
        source,
        dem_path,
    ]

    try:
        completed = subprocess.run(
            command,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if completed.stderr.strip():
            log.debug(completed.stderr.strip())

        if not _valid_cached_dem(
            dem_path,
            frame.source_epsg,
            dem_scale_m,
            required_bounds,
        ):
            raise RuntimeError(
                f"gdalwarp completed but produced an invalid DEM: {dem_path}"
            )

        if cache_path is not None:
            os.replace(dem_path, cache_path)
            log.info(
                "  DEM cache SAVED for frame %s -> %s",
                frame.info["frame"],
                cache_path,
            )
            return cache_path

        return dem_path

    except subprocess.CalledProcessError as error:
        detail = (error.stderr or error.stdout or str(error)).strip()
        try:
            if os.path.exists(dem_path):
                os.remove(dem_path)
        except OSError:
            pass
        raise RuntimeError(
            "Could not subset the NISAR DEM. If using the default HTTPS VRT, "
            "make sure Earthdata Login credentials are configured in ~/.netrc. "
            f"gdalwarp reported: {detail}"
        ) from error
    finally:
        if lock_handle is not None and not lock_handle.closed:
            try:
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
            finally:
                lock_handle.close()


def _clean_cube(dataset: h5py.Dataset) -> np.ndarray:
    array = np.asarray(dataset[:], dtype=np.float32)
    fill = _fill_value(dataset)
    if fill is not None:
        try:
            if np.isnan(fill):
                array[~np.isfinite(array)] = np.nan
            else:
                array[array == np.float32(fill)] = np.nan
        except TypeError:
            pass
    array[~np.isfinite(array)] = np.nan
    return array


def _prepare_radar_cube(
    h_axis: np.ndarray,
    y_axis: np.ndarray,
    x_axis: np.ndarray,
    cube: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Make radar-grid coordinate axes strictly increasing for interpolation."""
    h = np.asarray(h_axis, dtype=np.float64)
    y = np.asarray(y_axis, dtype=np.float64)
    x = np.asarray(x_axis, dtype=np.float64)
    data = np.asarray(cube, dtype=np.float32)

    if data.shape != (h.size, y.size, x.size):
        raise ValueError(
            f"Radar-grid cube shape {data.shape} does not match "
            f"axes ({h.size}, {y.size}, {x.size})"
        )

    if h[0] > h[-1]:
        h = h[::-1]
        data = data[::-1, :, :]
    if y[0] > y[-1]:
        y = y[::-1]
        data = data[:, ::-1, :]
    if x[0] > x[-1]:
        x = x[::-1]
        data = data[:, :, ::-1]

    if not (
        np.all(np.diff(h) > 0)
        and np.all(np.diff(y) > 0)
        and np.all(np.diff(x) > 0)
    ):
        raise ValueError("Radar-grid h/y/x coordinate axes are not monotonic")

    return h, y, x, data


def _fractional_indices(values: np.ndarray, axis: np.ndarray) -> np.ndarray:
    flat = np.asarray(values, dtype=np.float64).ravel()
    index_axis = np.arange(axis.size, dtype=np.float64)
    result = np.interp(flat, axis, index_axis)
    outside = (
        ~np.isfinite(flat)
        | (flat < axis[0])
        | (flat > axis[-1])
    )
    result[outside] = np.nan
    return result.reshape(np.shape(values))


def _interpolate_radar_cube(
    cube: np.ndarray,
    h_axis: np.ndarray,
    y_axis: np.ndarray,
    x_axis: np.ndarray,
    heights: np.ndarray,
    y_coords: np.ndarray,
    x_coords: np.ndarray,
) -> np.ndarray:
    """Trilinear interpolation of a [height, y, x] NISAR metadata cube."""
    hi = _fractional_indices(heights, h_axis)
    yi_1d = _fractional_indices(y_coords, y_axis)
    xi_1d = _fractional_indices(x_coords, x_axis)

    yi = np.broadcast_to(yi_1d[:, None], heights.shape)
    xi = np.broadcast_to(xi_1d[None, :], heights.shape)
    valid = np.isfinite(hi) & np.isfinite(yi) & np.isfinite(xi)

    coordinates = np.vstack(
        [
            np.where(valid, hi, 0.0).ravel(),
            np.where(valid, yi, 0.0).ravel(),
            np.where(valid, xi, 0.0).ravel(),
        ]
    )

    result = map_coordinates(
        cube,
        coordinates,
        order=1,
        mode="constant",
        cval=np.nan,
        prefilter=False,
    ).reshape(heights.shape)

    result[~valid] = np.nan
    return result.astype(np.float32)


def _read_dem_halo(
    source: rasterio.DatasetReader,
    output_window: Window,
    halo: int = 1,
) -> tuple[np.ndarray, tuple[slice, slice]]:
    row0 = max(0, int(output_window.row_off) - halo)
    col0 = max(0, int(output_window.col_off) - halo)
    row1 = min(
        source.height,
        int(output_window.row_off + output_window.height) + halo,
    )
    col1 = min(
        source.width,
        int(output_window.col_off + output_window.width) + halo,
    )

    read_window = Window(
        col0,
        row0,
        col1 - col0,
        row1 - row0,
    )

    array = source.read(1, window=read_window).astype(np.float32)
    if source.nodata is not None:
        array[array == np.float32(source.nodata)] = np.nan
    array[~np.isfinite(array)] = np.nan

    inner_row0 = int(output_window.row_off) - row0
    inner_col0 = int(output_window.col_off) - col0
    inner_row1 = inner_row0 + int(output_window.height)
    inner_col1 = inner_col0 + int(output_window.width)

    return array, (
        slice(inner_row0, inner_row1),
        slice(inner_col0, inner_col1),
    )


def prepare_geometry_frame(
    prepared: PreparedFrame,
    work_dir: str,
    dem_vrt: str,
    dem_scale_m: float,
    dem_cache_dir: Optional[str],
    earthdata_cookie: Optional[str],
    geometry_block_px: int,
    warp_threads: int,
    log: logging.Logger,
) -> PreparedGeometry:
    """
    Compute ellipsoid incidence angle and terrain-local incidence angle at the
    DEM terrain scale, already limited to the GCOV SAR support footprint.

    These temporary angle rasters are later area-averaged onto the requested
    destination grid (normally 100 m). Final output masking is additionally
    tied to the same co-pol/cross-pol validity and min_coverage test as SAR.
    """
    frame = prepared.meta
    radar_group_path = f"/science/{frame.root}/GCOV/metadata/radarGrid"

    dem_path = _prepare_dem_for_frame(
        frame=frame,
        work_dir=work_dir,
        dem_vrt=dem_vrt,
        dem_scale_m=dem_scale_m,
        warp_threads=warp_threads,
        log=log,
        dem_cache_dir=dem_cache_dir,
        earthdata_cookie=earthdata_cookie,
    )

    frame_label = f"frame_{frame.info['frame']}_{os.getpid()}"
    incidence_path = os.path.join(
        work_dir,
        f"{frame_label}_incidence_demscale.tif",
    )
    lia_path = os.path.join(
        work_dir,
        f"{frame_label}_lia_demscale.tif",
    )

    with h5py.File(frame.path, "r") as h5:
        required = [
            "xCoordinates",
            "yCoordinates",
            "heightAboveEllipsoid",
            "incidenceAngle",
            "losUnitVectorX",
            "losUnitVectorY",
        ]
        missing = [
            name
            for name in required
            if f"{radar_group_path}/{name}" not in h5
        ]
        if missing:
            raise KeyError(
                f"LIA requested but radarGrid datasets are missing in "
                f"{frame.path}: {missing}"
            )

        h0 = np.asarray(
            h5[f"{radar_group_path}/heightAboveEllipsoid"][:],
            dtype=np.float64,
        )
        y0 = np.asarray(
            h5[f"{radar_group_path}/yCoordinates"][:],
            dtype=np.float64,
        )
        x0 = np.asarray(
            h5[f"{radar_group_path}/xCoordinates"][:],
            dtype=np.float64,
        )

        incidence0 = _clean_cube(
            h5[f"{radar_group_path}/incidenceAngle"]
        )
        los_x0 = _clean_cube(
            h5[f"{radar_group_path}/losUnitVectorX"]
        )
        los_y0 = _clean_cube(
            h5[f"{radar_group_path}/losUnitVectorY"]
        )

    h_axis, y_axis, x_axis, incidence_cube = _prepare_radar_cube(
        h0, y0, x0, incidence0
    )
    _, _, _, los_x_cube = _prepare_radar_cube(
        h0, y0, x0, los_x0
    )
    _, _, _, los_y_cube = _prepare_radar_cube(
        h0, y0, x0, los_y0
    )

    with rasterio.open(dem_path) as dem_source, \
         rasterio.open(prepared.valid_path) as valid_source:

        if dem_source.crs != rasterio.CRS.from_epsg(frame.source_epsg):
            raise ValueError(
                f"DEM CRS {dem_source.crs} does not match frame "
                f"EPSG:{frame.source_epsg}"
            )

        if abs(dem_source.transform.b) > 1e-12 or abs(dem_source.transform.d) > 1e-12:
            raise ValueError("DEM terrain grid must be north-up without rotation")

        profile = dem_source.profile.copy()
        profile.update(
            driver="GTiff",
            count=1,
            dtype="float32",
            nodata=float(NODATA),
            compress="DEFLATE",
            predictor=3,
            tiled=True,
            blockxsize=512,
            blockysize=512,
            BIGTIFF="IF_SAFER",
        )

        valid_on_dem = WarpedVRT(
            valid_source,
            crs=dem_source.crs,
            transform=dem_source.transform,
            width=dem_source.width,
            height=dem_source.height,
            nodata=0.0,
            resampling=Resampling.average,
            init_dest_nodata=True,
            warp_mem_limit=256,
            dtype="float32",
            NUM_THREADS=str(warp_threads),
        )

        windows = iter_windows(
            dem_source.width,
            dem_source.height,
            max(64, geometry_block_px),
        )
        report_step = max(1, len(windows) // 5)

        log.info(
            "  Computing frame %s incidence/LIA on DEM terrain grid %dx%d",
            frame.info["frame"],
            dem_source.height,
            dem_source.width,
        )

        with valid_on_dem, \
             rasterio.open(incidence_path, "w", **profile) as inc_dst, \
             rasterio.open(lia_path, "w", **profile) as lia_dst:

            inc_dst.update_tags(
                source_hdf5=os.path.basename(frame.path),
                source_frame=frame.info["frame"],
                quantity="NISAR ellipsoid-referenced incidence angle",
                units="degrees",
                source_dataset=f"{radar_group_path}/incidenceAngle",
                terrain_scale_m=str(dem_scale_m),
                sar_support_mask_applied="yes",
            )
            lia_dst.update_tags(
                source_hdf5=os.path.basename(frame.path),
                source_frame=frame.info["frame"],
                quantity="terrain local incidence angle",
                units="degrees",
                definition=(
                    "acos(LOS_target_to_sensor dot upward_terrain_normal)"
                ),
                terrain_scale_m=str(dem_scale_m),
                sar_support_mask_applied="yes",
                dem_source=str(dem_vrt),
            )

            dx = float(dem_source.transform.a)
            dy = float(dem_source.transform.e)

            for index, window in enumerate(windows, start=1):
                dem_ext, crop = _read_dem_halo(
                    dem_source,
                    window,
                    halo=1,
                )

                coverage = valid_on_dem.read(
                    1,
                    window=window,
                    out_dtype="float32",
                )

                out_shape = (int(window.height), int(window.width))
                inc_out = np.full(out_shape, NODATA, dtype=np.float32)
                lia_out = np.full(out_shape, NODATA, dtype=np.float32)

                if np.any(np.isfinite(coverage) & (coverage > 0.0)):
                    dz_dy_ext, dz_dx_ext = np.gradient(
                        dem_ext.astype(np.float64),
                        dy,
                        dx,
                        edge_order=1,
                    )

                    dem = dem_ext[crop].astype(np.float64)
                    dz_dx = dz_dx_ext[crop]
                    dz_dy = dz_dy_ext[crop]

                    normal_mag = np.sqrt(
                        1.0 + dz_dx * dz_dx + dz_dy * dz_dy
                    )
                    normal_x = -dz_dx / normal_mag
                    normal_y = -dz_dy / normal_mag
                    normal_z = 1.0 / normal_mag

                    cols = (
                        np.arange(
                            int(window.col_off),
                            int(window.col_off + window.width),
                            dtype=np.float64,
                        )
                        + 0.5
                    )
                    rows = (
                        np.arange(
                            int(window.row_off),
                            int(window.row_off + window.height),
                            dtype=np.float64,
                        )
                        + 0.5
                    )

                    x_coords = (
                        dem_source.transform.c
                        + cols * dem_source.transform.a
                    )
                    y_coords = (
                        dem_source.transform.f
                        + rows * dem_source.transform.e
                    )

                    incidence = _interpolate_radar_cube(
                        incidence_cube,
                        h_axis,
                        y_axis,
                        x_axis,
                        dem,
                        y_coords,
                        x_coords,
                    ).astype(np.float64)

                    los_x = _interpolate_radar_cube(
                        los_x_cube,
                        h_axis,
                        y_axis,
                        x_axis,
                        dem,
                        y_coords,
                        x_coords,
                    ).astype(np.float64)

                    los_y = _interpolate_radar_cube(
                        los_y_cube,
                        h_axis,
                        y_axis,
                        x_axis,
                        dem,
                        y_coords,
                        x_coords,
                    ).astype(np.float64)

                    los_z_squared = 1.0 - los_x * los_x - los_y * los_y
                    los_z_squared = np.where(
                        los_z_squared >= -1e-6,
                        np.maximum(los_z_squared, 0.0),
                        np.nan,
                    )
                    los_z = np.sqrt(los_z_squared)

                    dot = (
                        los_x * normal_x
                        + los_y * normal_y
                        + los_z * normal_z
                    )
                    dot = np.clip(dot, -1.0, 1.0)
                    lia = np.degrees(np.arccos(dot))

                    supported = (
                        np.isfinite(coverage)
                        & (coverage > 0.0)
                        & np.isfinite(dem)
                        & np.isfinite(incidence)
                        & np.isfinite(lia)
                    )

                    inc_out[supported] = incidence[supported].astype(np.float32)
                    lia_out[supported] = lia[supported].astype(np.float32)

                inc_dst.write(inc_out, 1, window=window)
                lia_dst.write(lia_out, 1, window=window)

                if index % report_step == 0 or index == len(windows):
                    log.info(
                        "    geometry blocks completed: %d/%d",
                        index,
                        len(windows),
                    )

    return PreparedGeometry(
        incidence_path=incidence_path,
        lia_path=lia_path,
    )


def destination_profile(
    width: int,
    height: int,
    transform: Affine,
    target_epsg: int,
    dtype: str,
    nodata,
) -> dict:
    return {
        "driver": "GTiff",
        "width": width,
        "height": height,
        "count": 1,
        "dtype": dtype,
        "crs": rasterio.CRS.from_epsg(target_epsg),
        "transform": transform,
        "nodata": nodata,
        "compress": "DEFLATE",
        "predictor": 2 if dtype.startswith("float") else 1,
        "tiled": True,
        "blockxsize": 512,
        "blockysize": 512,
        "BIGTIFF": "IF_SAFER",
    }


def make_quicklook(
    co_path: str,
    cross_path: str,
    out_path: str,
    co_key: str,
    cross_key: str,
    normalization: str,
    incidence_path: Optional[str] = None,
    lia_path: Optional[str] = None,
) -> None:
    """
    Create a diagnostic quicklook.

    The GeoTIFF SAR outputs remain LINEAR power.  Only the display arrays are
    converted to dB using 10*log10(power).  When incidence/LIA paths are
    supplied, the quicklook is a 2x2 panel:

        co-pol dB | cross-pol dB
        LIA (deg) | incidence angle (deg)

    Geometry panels therefore use exactly the final SAR-valid geometry rasters
    written by this processor.
    """

    geometry_requested = bool(incidence_path and lia_path)

    if geometry_requested:
        figure, axes = plt.subplots(
            2,
            2,
            figsize=(16, 13),
            dpi=140,
            squeeze=False,
        )
        sar_axes = (axes[0, 0], axes[0, 1])
        geometry_axes = (axes[1, 0], axes[1, 1])
    else:
        figure, axes_1d = plt.subplots(
            1,
            2,
            figsize=(16, 7),
            dpi=140,
        )
        sar_axes = tuple(np.atleast_1d(axes_1d))
        geometry_axes = ()

    def read_for_display(path: str) -> tuple[np.ndarray, list[float]]:
        with rasterio.open(path) as source:
            reduction = max(
                1.0,
                max(source.width, source.height) / 1800.0,
            )

            out_height = max(
                1,
                int(round(source.height / reduction)),
            )
            out_width = max(
                1,
                int(round(source.width / reduction)),
            )

            data = source.read(
                1,
                out_shape=(out_height, out_width),
                resampling=Resampling.average,
            ).astype(np.float32)

            if source.nodata is not None:
                data[data == source.nodata] = np.nan

            extent = [
                source.bounds.left,
                source.bounds.right,
                source.bounds.bottom,
                source.bounds.top,
            ]

        return data, extent

    def display_pol_name(key: str) -> str:
        # Diagonal covariance names such as HHHH/HVHV represent HH/HV power.
        return {
            "HHHH": "HH",
            "HVHV": "HV",
            "VVVV": "VV",
            "VHVH": "VH",
            "RHRH": "RH",
            "RVRV": "RV",
        }.get(key, key)

    # ------------------------------------------------------------------
    # SAR panels: linear GeoTIFF -> dB for plotting ONLY.
    # ------------------------------------------------------------------
    for axis, path, label in zip(
        sar_axes,
        (co_path, cross_path),
        (co_key, cross_key),
    ):
        linear, extent = read_for_display(path)

        db = np.full(linear.shape, np.nan, dtype=np.float32)
        positive = np.isfinite(linear) & (linear > 0.0)
        db[positive] = (10.0 * np.log10(linear[positive])).astype(np.float32)

        if not np.any(np.isfinite(db)):
            axis.text(
                0.5,
                0.5,
                "No valid data",
                ha="center",
                va="center",
                transform=axis.transAxes,
            )
            axis.set_title(f"{display_pol_name(label)} ({normalization}, dB)")
            continue

        vmin, vmax = np.nanpercentile(db, [2, 98])
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
            vmin = float(np.nanmin(db))
            vmax = float(np.nanmax(db))

        image = axis.imshow(
            db,
            extent=extent,
            origin="upper",
            cmap="gray",
            vmin=vmin,
            vmax=vmax,
        )

        figure.colorbar(
            image,
            ax=axis,
            label=f"{normalization} (dB)",
        )

        axis.set_title(f"{display_pol_name(label)} ({normalization}, dB)")
        axis.set_xlabel("Easting (m)")
        axis.set_ylabel("Northing (m)")

    # ------------------------------------------------------------------
    # Optional second row: LIA and ellipsoid incidence angle.
    # ------------------------------------------------------------------
    if geometry_requested:
        geometry_specs = (
            (lia_path, "Local Incidence Angle"),
            (incidence_path, "Incidence Angle"),
        )

        for axis, (path, title) in zip(geometry_axes, geometry_specs):
            assert path is not None
            data, extent = read_for_display(path)
            valid = np.isfinite(data)

            if not np.any(valid):
                axis.text(
                    0.5,
                    0.5,
                    "No valid data",
                    ha="center",
                    va="center",
                    transform=axis.transAxes,
                )
                axis.set_title(title)
                continue

            vmin, vmax = np.nanpercentile(data, [2, 98])
            if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
                vmin = float(np.nanmin(data))
                vmax = float(np.nanmax(data))

            image = axis.imshow(
                data,
                extent=extent,
                origin="upper",
                cmap="viridis",
                vmin=vmin,
                vmax=vmax,
            )

            figure.colorbar(
                image,
                ax=axis,
                label="degrees",
            )

            axis.set_title(title)
            axis.set_xlabel("Easting (m)")
            axis.set_ylabel("Northing (m)")

    figure.tight_layout()
    figure.savefig(
        out_path,
        bbox_inches="tight",
    )
    plt.close(figure)


# ============================================================================
# EXISTING-OUTPUT / RESUME HELPERS
# ============================================================================

def validate_existing_output(
    path: str,
    width: int,
    height: int,
    transform: Affine,
    target_epsg: int,
    expected_dtype: str,
    expected_tags: dict,
    log: logging.Logger,
) -> bool:
    """
    Check whether an existing GeoTIFF matches the current requested output.

    IMPORTANT:
    Validation failure NEVER causes the existing file to be overwritten.
    The caller will skip that acquisition and report a warning.
    """
    if not os.path.isfile(path):
        return False

    try:
        with rasterio.open(path) as src:
            if src.count != 1:
                log.warning(
                    "Existing output band-count mismatch: found=%d expected=1: %s",
                    src.count,
                    path,
                )
                return False

            if src.width != width or src.height != height:
                log.warning(
                    "Existing output dimension mismatch: found=%dx%d expected=%dx%d: %s",
                    src.width,
                    src.height,
                    width,
                    height,
                    path,
                )
                return False

            expected_crs = rasterio.CRS.from_epsg(target_epsg)
            if src.crs != expected_crs:
                log.warning(
                    "Existing output CRS mismatch: found=%s expected=EPSG:%d: %s",
                    src.crs,
                    target_epsg,
                    path,
                )
                return False

            if not src.transform.almost_equals(transform):
                log.warning(
                    "Existing output transform/grid mismatch: %s",
                    path,
                )
                return False

            if src.dtypes[0] != expected_dtype:
                log.warning(
                    "Existing output dtype mismatch: found=%s expected=%s: %s",
                    src.dtypes[0],
                    expected_dtype,
                    path,
                )
                return False

            tags = src.tags()

            for key, expected_value in expected_tags.items():
                actual_value = tags.get(key)
                expected_text = str(expected_value)

                if actual_value != expected_text:
                    log.warning(
                        "Existing output metadata mismatch for %s: "
                        "found=%r expected=%r: %s",
                        key,
                        actual_value,
                        expected_text,
                        path,
                    )
                    return False

    except Exception as error:
        log.warning(
            "Could not validate existing GeoTIFF %s: %s",
            path,
            error,
        )
        return False

    return True


def log_existing_output_set(
    log: logging.Logger,
    status: str,
    first: dict,
    frame_string: str,
    co_path: str,
    cross_path: str,
    count_path: str,
) -> None:
    log.info("=" * 80)
    log.info(status)
    log.info("  Cycle       : %s", first["cycle"])
    log.info("  Track       : %s", first["track"])
    log.info("  Pass        : %s", first["pass"])
    log.info("  Date        : %s", first["date"])
    log.info("  Frame(s)    : %s", frame_string)
    log.info("  Co-pol      : %s", co_path)
    log.info("  Cross-pol   : %s", cross_path)
    log.info("  Frame count : %s", count_path)
    log.info("=" * 80)


def process_cluster(
    cluster: list[FrameMeta],
    output_dir: str,
    cluster_index: int,
    n_clusters: int,
    grid_mode: str,
    target_epsg: int,
    resolution: float,
    origin_x: float,
    origin_y: float,
    bbox_wgs84: Optional[tuple[float, float, float, float]],
    output_tile_km: float,
    native_block_px: int,
    warp_threads: int,
    normalization: str,
    filter_name: str,
    filter_window: int,
    noise_floor: float,
    use_exception_mask: bool,
    min_coverage: float,
    make_lia: bool,
    nisar_dem_vrt: str,
    lia_dem_scale: float,
    dem_cache_dir: Optional[str],
    earthdata_cookie: Optional[str],
    keep_intermediate: bool,
    log: logging.Logger,
) -> tuple[str, Optional[tuple[str, ...]]]:
    """Process one spatial cluster with optional SAR-valid incidence/LIA outputs."""

    target_transform, width, height, bounds, grid_metadata = build_target_grid(
        cluster,
        resolution,
        origin_x,
        origin_y,
        bbox_wgs84,
        target_epsg,
        grid_mode,
    )

    output_tile_px = max(
        1,
        int(round(output_tile_km * 1000.0 / resolution)),
    )
    output_windows = iter_windows(width, height, output_tile_px)

    first = cluster[0].info
    frame_string = compress_frames(
        frame.info["frame"] for frame in cluster
    )

    filter_tag = (
        "nofilt"
        if filter_name == "none"
        else f"{filter_name}_w{filter_window}"
    )

    stem = (
        f"NISAR_{first['instrument']}_GCOV_"
        f"c{first['cycle']}_t{first['track']}_{first['pass']}_"
        f"f{frame_string}_{first['bw_mode']}_{first['pole']}_"
        f"{cluster[0].frequency}_{first['date']}_"
        f"{int(resolution)}m_{filter_tag}"
    )

    co_path = os.path.join(
        output_dir,
        f"{stem}_{cluster[0].co_key}.tif",
    )
    cross_path = os.path.join(
        output_dir,
        f"{stem}_{cluster[0].cross_key}.tif",
    )
    count_path = os.path.join(
        output_dir,
        f"{stem}_frame_count.tif",
    )
    incidence_path = os.path.join(
        output_dir,
        f"{stem}_incidence_angle_deg.tif",
    )
    lia_path = os.path.join(
        output_dir,
        f"{stem}_local_incidence_angle_deg.tif",
    )
    quicklook_path = os.path.join(
        output_dir,
        f"{stem}_quicklook.png",
    )

    ease2_grid_mat_path = (
        os.path.join(output_dir, f"{stem}_EASE2_grid.mat")
        if grid_mode.upper() == "EASE2"
        else None
    )

    if ease2_grid_mat_path is not None:
        ensure_ease2_grid_companion(
            path=ease2_grid_mat_path,
            transform=target_transform,
            width=width,
            height=height,
            grid_metadata=grid_metadata,
            log=log,
        )

    source_frames = ",".join(
        frame.info["frame"] for frame in cluster
    )
    source_files = "|".join(
        os.path.basename(frame.path) for frame in cluster
    )

    expected_common_tags = {
        **grid_metadata,
        "target_epsg": str(target_epsg),
        "resolution_m": str(resolution),
        "grid_origin_x": str(origin_x),
        "grid_origin_y": str(origin_y),
        "normalization": normalization,
        "filter": filter_name,
        "minimum_valid_source_fraction": str(min_coverage),
        "source_frames": source_frames,
        "source_files": source_files,
        "input_data_exception_mask_applied": (
            "yes" if use_exception_mask else "no"
        ),
    }

    base_paths = (co_path, cross_path, count_path)
    geometry_paths = (
        (incidence_path, lia_path)
        if make_lia
        else tuple()
    )
    companion_paths = (
        (ease2_grid_mat_path,)
        if ease2_grid_mat_path is not None
        else tuple()
    )
    expected_paths = base_paths + geometry_paths + companion_paths

    def validate_base() -> bool:
        return (
            validate_existing_output(
                path=co_path,
                width=width,
                height=height,
                transform=target_transform,
                target_epsg=target_epsg,
                expected_dtype="float32",
                expected_tags={
                    **expected_common_tags,
                    "covariance_term": cluster[0].co_key,
                },
                log=log,
            )
            and validate_existing_output(
                path=cross_path,
                width=width,
                height=height,
                transform=target_transform,
                target_epsg=target_epsg,
                expected_dtype="float32",
                expected_tags={
                    **expected_common_tags,
                    "covariance_term": cluster[0].cross_key,
                },
                log=log,
            )
            and validate_existing_output(
                path=count_path,
                width=width,
                height=height,
                transform=target_transform,
                target_epsg=target_epsg,
                expected_dtype="uint16",
                expected_tags=expected_common_tags,
                log=log,
            )
        )

    def validate_geometry() -> bool:
        if not make_lia:
            return True
        geometry_common = {
            **grid_metadata,
            "target_epsg": str(target_epsg),
            "resolution_m": str(resolution),
            "source_frames": source_frames,
            "source_files": source_files,
            "sar_valid_mask_applied": "yes",
            "lia_dem_scale_m": str(lia_dem_scale),
        }
        return (
            validate_existing_output(
                path=incidence_path,
                width=width,
                height=height,
                transform=target_transform,
                target_epsg=target_epsg,
                expected_dtype="float32",
                expected_tags={
                    **geometry_common,
                    "geometry_quantity": "incidence_angle",
                },
                log=log,
            )
            and validate_existing_output(
                path=lia_path,
                width=width,
                height=height,
                transform=target_transform,
                target_epsg=target_epsg,
                expected_dtype="float32",
                expected_tags={
                    **geometry_common,
                    "geometry_quantity": "local_incidence_angle",
                },
                log=log,
            )
        )

    base_flags = tuple(os.path.isfile(path) for path in base_paths)
    geometry_flags = tuple(
        os.path.isfile(path) for path in geometry_paths
    )
    geometry_only = False

    # Existing SAR outputs may be augmented with LIA without overwriting them.
    if all(base_flags):
        if not validate_base():
            log_existing_output_set(
                log=log,
                status=(
                    "SKIPPING EXISTING ACQUISITION: SAR OUTPUT VALIDATION "
                    "REPORTED A MISMATCH"
                ),
                first=first,
                frame_string=frame_string,
                co_path=co_path,
                cross_path=cross_path,
                count_path=count_path,
            )
            return "skipped_existing_mismatch", expected_paths

        if not make_lia:
            log_existing_output_set(
                log=log,
                status="SKIPPING ALREADY PROCESSED ACQUISITION",
                first=first,
                frame_string=frame_string,
                co_path=co_path,
                cross_path=cross_path,
                count_path=count_path,
            )
            return "skipped_existing", expected_paths

        if geometry_flags and all(geometry_flags):
            if validate_geometry():
                log_existing_output_set(
                    log=log,
                    status="SKIPPING ALREADY PROCESSED ACQUISITION + LIA",
                    first=first,
                    frame_string=frame_string,
                    co_path=co_path,
                    cross_path=cross_path,
                    count_path=count_path,
                )
                log.info("  Incidence   : %s", incidence_path)
                log.info("  Local inc.  : %s", lia_path)
                return "skipped_existing", expected_paths

            log.warning(
                "Existing incidence/LIA outputs failed validation and will NOT "
                "be overwritten."
            )
            return "skipped_existing_mismatch", expected_paths

        if any(geometry_flags):
            log.warning(
                "Only part of the requested incidence/LIA output set exists. "
                "Nothing will be overwritten. Move/remove the partial geometry "
                "output before rerunning."
            )
            return "skipped_partial", None

        geometry_only = True
        log.info("=" * 80)
        log.info("EXISTING SAR FOUND - ADDING 100 m INCIDENCE/LIA ONLY")
        log.info("  Cycle    : %s", first["cycle"])
        log.info("  Track    : %s", first["track"])
        log.info("  Pass     : %s", first["pass"])
        log.info("  Date     : %s", first["date"])
        log.info("  Frame(s) : %s", frame_string)
        log.info("=" * 80)

    elif any(base_flags):
        log.warning(
            "Partial existing SAR output set detected; acquisition skipped to "
            "preserve no-overwrite behavior."
        )
        return "skipped_partial", None

    elif make_lia and any(geometry_flags):
        log.warning(
            "Geometry outputs exist but the SAR output set is absent. Nothing "
            "will be overwritten; acquisition skipped."
        )
        return "skipped_partial", None

    if not geometry_only:
        log.info("=" * 80)
        log.info("NEW ACQUISITION - PROCESSING")
        log.info("  Cycle    : %s", first["cycle"])
        log.info("  Track    : %s", first["track"])
        log.info("  Pass     : %s", first["pass"])
        log.info("  Date     : %s", first["date"])
        log.info("  Frame(s) : %s", frame_string)
        log.info("  Filter   : %s", filter_tag)
        log.info("  LIA      : %s", "yes" if make_lia else "no")
        log.info(
            "  Exc.mask : %s",
            "applied" if use_exception_mask else "NOT applied (kept)",
        )
        log.info("=" * 80)

    work_root = os.path.join(output_dir, ".nisar_work")
    os.makedirs(work_root, exist_ok=True)
    work_dir = tempfile.mkdtemp(prefix=f"{stem}_", dir=work_root)

    log.info(
        "Cluster %d/%d: frames=%s; target grid=%dx%d; bounds=%s; "
        "output windows=%d",
        cluster_index,
        n_clusters,
        frame_string,
        height,
        width,
        tuple(round(v, 1) for v in bounds),
        len(output_windows),
    )

    if grid_mode.upper() == "EASE2":
        log.info(
            "  EASE2 M0p1 global rows [%s, %s), cols [%s, %s)",
            grid_metadata["ease2_global_row_start"],
            grid_metadata["ease2_global_row_end_exclusive"],
            grid_metadata["ease2_global_col_start"],
            grid_metadata["ease2_global_col_end_exclusive"],
        )

    try:
        # Even in geometry-only add-on mode, create temporary frame-level SAR
        # rasters. They reproduce the exact per-frame validity used by the SAR
        # output and guarantee incidence/LIA are never valid where HH/HV are not.
        prepared = [
            prepare_native_frame(
                frame=frame,
                work_dir=work_dir,
                normalization=normalization,
                filter_name=filter_name,
                filter_window=filter_window,
                noise_floor=noise_floor,
                use_exception_mask=use_exception_mask,
                native_block_px=native_block_px,
                log=log,
            )
            for frame in cluster
        ]

        prepared_geometry: list[PreparedGeometry] = []
        if make_lia:
            prepared_geometry = [
                prepare_geometry_frame(
                    prepared=item,
                    work_dir=work_dir,
                    dem_vrt=nisar_dem_vrt,
                    dem_scale_m=lia_dem_scale,
                    dem_cache_dir=dem_cache_dir,
                    earthdata_cookie=earthdata_cookie,
                    geometry_block_px=native_block_px,
                    warp_threads=warp_threads,
                    log=log,
                )
                for item in prepared
            ]

        float_profile = destination_profile(
            width,
            height,
            target_transform,
            target_epsg,
            "float32",
            float(NODATA),
        )
        count_profile = destination_profile(
            width,
            height,
            target_transform,
            target_epsg,
            "uint16",
            0,
        )

        temp_co_path = os.path.join(
            work_dir,
            f"final_{cluster[0].co_key}.tif",
        )
        temp_cross_path = os.path.join(
            work_dir,
            f"final_{cluster[0].cross_key}.tif",
        )
        temp_count_path = os.path.join(work_dir, "final_frame_count.tif")
        temp_incidence_path = os.path.join(
            work_dir,
            "final_incidence_angle_deg.tif",
        )
        temp_lia_path = os.path.join(
            work_dir,
            "final_local_incidence_angle_deg.tif",
        )

        with ExitStack() as stack:
            co_vrts = []
            cross_vrts = []
            coverage_vrts = []
            incidence_vrts = []
            lia_vrts = []

            for frame_index, item in enumerate(prepared):
                co_source = stack.enter_context(rasterio.open(item.co_path))
                cross_source = stack.enter_context(rasterio.open(item.cross_path))
                coverage_source = stack.enter_context(rasterio.open(item.valid_path))

                co_vrts.append(
                    stack.enter_context(
                        WarpedVRT(
                            co_source,
                            crs=rasterio.CRS.from_epsg(target_epsg),
                            transform=target_transform,
                            width=width,
                            height=height,
                            src_nodata=float(NODATA),
                            nodata=float(NODATA),
                            resampling=Resampling.average,
                            init_dest_nodata=True,
                            warp_mem_limit=256,
                            dtype="float32",
                            NUM_THREADS=str(warp_threads),
                        )
                    )
                )
                cross_vrts.append(
                    stack.enter_context(
                        WarpedVRT(
                            cross_source,
                            crs=rasterio.CRS.from_epsg(target_epsg),
                            transform=target_transform,
                            width=width,
                            height=height,
                            src_nodata=float(NODATA),
                            nodata=float(NODATA),
                            resampling=Resampling.average,
                            init_dest_nodata=True,
                            warp_mem_limit=256,
                            dtype="float32",
                            NUM_THREADS=str(warp_threads),
                        )
                    )
                )
                coverage_vrts.append(
                    stack.enter_context(
                        WarpedVRT(
                            coverage_source,
                            crs=rasterio.CRS.from_epsg(target_epsg),
                            transform=target_transform,
                            width=width,
                            height=height,
                            nodata=0.0,
                            resampling=Resampling.average,
                            init_dest_nodata=True,
                            warp_mem_limit=256,
                            dtype="float32",
                            NUM_THREADS=str(warp_threads),
                        )
                    )
                )

                if make_lia:
                    geometry = prepared_geometry[frame_index]
                    inc_source = stack.enter_context(
                        rasterio.open(geometry.incidence_path)
                    )
                    lia_source = stack.enter_context(
                        rasterio.open(geometry.lia_path)
                    )
                    incidence_vrts.append(
                        stack.enter_context(
                            WarpedVRT(
                                inc_source,
                                crs=rasterio.CRS.from_epsg(target_epsg),
                                transform=target_transform,
                                width=width,
                                height=height,
                                src_nodata=float(NODATA),
                                nodata=float(NODATA),
                                resampling=Resampling.average,
                                init_dest_nodata=True,
                                warp_mem_limit=256,
                                dtype="float32",
                                NUM_THREADS=str(warp_threads),
                            )
                        )
                    )
                    lia_vrts.append(
                        stack.enter_context(
                            WarpedVRT(
                                lia_source,
                                crs=rasterio.CRS.from_epsg(target_epsg),
                                transform=target_transform,
                                width=width,
                                height=height,
                                src_nodata=float(NODATA),
                                nodata=float(NODATA),
                                resampling=Resampling.average,
                                init_dest_nodata=True,
                                warp_mem_limit=256,
                                dtype="float32",
                                NUM_THREADS=str(warp_threads),
                            )
                        )
                    )

            write_sar = not geometry_only
            co_destination = (
                stack.enter_context(rasterio.open(temp_co_path, "w", **float_profile))
                if write_sar
                else None
            )
            cross_destination = (
                stack.enter_context(
                    rasterio.open(temp_cross_path, "w", **float_profile)
                )
                if write_sar
                else None
            )
            count_destination = (
                stack.enter_context(
                    rasterio.open(temp_count_path, "w", **count_profile)
                )
                if write_sar
                else None
            )
            incidence_destination = (
                stack.enter_context(
                    rasterio.open(temp_incidence_path, "w", **float_profile)
                )
                if make_lia
                else None
            )
            lia_destination = (
                stack.enter_context(
                    rasterio.open(temp_lia_path, "w", **float_profile)
                )
                if make_lia
                else None
            )

            tags = {
                **grid_metadata,
                "processing_grid": (
                    "EASE-Grid 2.0 M0p1"
                    if grid_mode.upper() == "EASE2"
                    else "canonical fixed projected grid"
                ),
                "target_epsg": str(target_epsg),
                "resolution_m": str(resolution),
                "grid_origin_x": str(origin_x),
                "grid_origin_y": str(origin_y),
                "mosaic_method": (
                    "equal-weight arithmetic mean of valid frame-level values"
                ),
                "normalization": normalization,
                "filter": filter_name,
                "filter_window": str(filter_window),
                "noise_floor_linear": str(noise_floor),
                "minimum_valid_source_fraction": str(min_coverage),
                "source_frames": source_frames,
                "source_files": source_files,
                "input_data_exception_mask_applied": (
                    "yes" if use_exception_mask else "no"
                ),
            }
            if ease2_grid_mat_path is not None:
                tags["ease2_grid_companion"] = os.path.basename(
                    ease2_grid_mat_path
                )

            if bbox_wgs84 is not None:
                tags["bbox_wgs84"] = ",".join(
                    str(value) for value in bbox_wgs84
                )

            if write_sar:
                assert co_destination is not None
                assert cross_destination is not None
                assert count_destination is not None
                co_destination.update_tags(
                    **tags,
                    covariance_term=cluster[0].co_key,
                )
                cross_destination.update_tags(
                    **tags,
                    covariance_term=cluster[0].cross_key,
                )
                count_destination.update_tags(
                    **tags,
                    description=(
                        "number of valid source frames contributing to each "
                        "output pixel"
                    ),
                )

            if make_lia:
                assert incidence_destination is not None
                assert lia_destination is not None
                geometry_tags = {
                    **grid_metadata,
                    "processing_grid": (
                        "EASE-Grid 2.0 M0p1"
                        if grid_mode.upper() == "EASE2"
                        else "canonical fixed projected grid"
                    ),
                    "target_epsg": str(target_epsg),
                    "resolution_m": str(resolution),
                    "grid_origin_x": str(origin_x),
                    "grid_origin_y": str(origin_y),
                    "source_frames": source_frames,
                    "source_files": source_files,
                    "sar_valid_mask_applied": "yes",
                    "minimum_valid_source_fraction": str(min_coverage),
                    "lia_dem_scale_m": str(lia_dem_scale),
                    "nisar_dem_source": str(nisar_dem_vrt),
                    "geometry_aggregation": (
                        "DEM-scale angles area-averaged to target grid; "
                        "then retained only where the same frame contributes "
                        "valid co-pol and cross-pol SAR"
                    ),
                }
                if ease2_grid_mat_path is not None:
                    geometry_tags["ease2_grid_companion"] = os.path.basename(
                        ease2_grid_mat_path
                    )

                incidence_destination.update_tags(
                    **geometry_tags,
                    geometry_quantity="incidence_angle",
                    units="degrees",
                    definition=(
                        "angle between LOS and ellipsoid normal at target height"
                    ),
                )
                lia_destination.update_tags(
                    **geometry_tags,
                    geometry_quantity="local_incidence_angle",
                    units="degrees",
                    definition=(
                        "angle between LOS target-to-sensor vector and "
                        "DEM-derived upward terrain normal"
                    ),
                )

            report_step = max(1, len(output_windows) // 10)

            for index, window in enumerate(output_windows, start=1):
                shape = (int(window.height), int(window.width))
                co_sum = np.zeros(shape, dtype=np.float64)
                cross_sum = np.zeros(shape, dtype=np.float64)
                co_count = np.zeros(shape, dtype=np.uint16)
                cross_count = np.zeros(shape, dtype=np.uint16)

                if make_lia:
                    incidence_sum = np.zeros(shape, dtype=np.float64)
                    lia_sum = np.zeros(shape, dtype=np.float64)
                    geometry_count = np.zeros(shape, dtype=np.uint16)

                for frame_index, (co_vrt, cross_vrt, coverage_vrt) in enumerate(
                    zip(co_vrts, cross_vrts, coverage_vrts)
                ):
                    co = co_vrt.read(1, window=window, out_dtype="float32")
                    cross = cross_vrt.read(1, window=window, out_dtype="float32")
                    coverage = coverage_vrt.read(
                        1,
                        window=window,
                        out_dtype="float32",
                    )

                    sufficiently_covered = (
                        np.isfinite(coverage)
                        & (coverage >= min_coverage)
                    )
                    valid_co = (
                        np.isfinite(co)
                        & (co != NODATA)
                        & sufficiently_covered
                    )
                    valid_cross = (
                        np.isfinite(cross)
                        & (cross != NODATA)
                        & sufficiently_covered
                    )

                    co_sum[valid_co] += co[valid_co]
                    cross_sum[valid_cross] += cross[valid_cross]
                    co_count[valid_co] += 1
                    cross_count[valid_cross] += 1

                    if make_lia:
                        incidence = incidence_vrts[frame_index].read(
                            1,
                            window=window,
                            out_dtype="float32",
                        )
                        lia = lia_vrts[frame_index].read(
                            1,
                            window=window,
                            out_dtype="float32",
                        )

                        # Critical design rule: geometry is valid only where the
                        # SAME frame actually contributes both SAR channels.
                        valid_geometry = (
                            valid_co
                            & valid_cross
                            & np.isfinite(incidence)
                            & (incidence != NODATA)
                            & np.isfinite(lia)
                            & (lia != NODATA)
                        )
                        incidence_sum[valid_geometry] += incidence[valid_geometry]
                        lia_sum[valid_geometry] += lia[valid_geometry]
                        geometry_count[valid_geometry] += 1

                common_count = np.minimum(
                    co_count,
                    cross_count,
                ).astype(np.uint16)

                if write_sar:
                    co_output = np.full(shape, NODATA, dtype=np.float32)
                    cross_output = np.full(shape, NODATA, dtype=np.float32)
                    valid_co_output = co_count > 0
                    valid_cross_output = cross_count > 0
                    co_output[valid_co_output] = (
                        co_sum[valid_co_output] / co_count[valid_co_output]
                    ).astype(np.float32)
                    cross_output[valid_cross_output] = (
                        cross_sum[valid_cross_output]
                        / cross_count[valid_cross_output]
                    ).astype(np.float32)

                    co_destination.write(co_output, 1, window=window)
                    cross_destination.write(cross_output, 1, window=window)
                    count_destination.write(common_count, 1, window=window)

                if make_lia:
                    incidence_output = np.full(shape, NODATA, dtype=np.float32)
                    lia_output = np.full(shape, NODATA, dtype=np.float32)
                    geometry_valid = (
                        (geometry_count > 0)
                        & (common_count > 0)
                    )
                    incidence_output[geometry_valid] = (
                        incidence_sum[geometry_valid]
                        / geometry_count[geometry_valid]
                    ).astype(np.float32)
                    lia_output[geometry_valid] = (
                        lia_sum[geometry_valid]
                        / geometry_count[geometry_valid]
                    ).astype(np.float32)
                    incidence_destination.write(
                        incidence_output,
                        1,
                        window=window,
                    )
                    lia_destination.write(lia_output, 1, window=window)

                if index % report_step == 0 or index == len(output_windows):
                    log.info(
                        "  fixed-grid windows completed: %d/%d",
                        index,
                        len(output_windows),
                    )

        paths_to_publish: list[tuple[str, str]] = []
        if not geometry_only:
            paths_to_publish.extend(
                [
                    (temp_co_path, co_path),
                    (temp_cross_path, cross_path),
                    (temp_count_path, count_path),
                ]
            )
        if make_lia:
            paths_to_publish.extend(
                [
                    (temp_incidence_path, incidence_path),
                    (temp_lia_path, lia_path),
                ]
            )

        if any(os.path.exists(final) for _, final in paths_to_publish):
            log.warning(
                "An expected final output appeared while this acquisition was "
                "processing. Temporary results will be discarded to preserve "
                "no-overwrite behavior."
            )
            return "skipped_race_existing", None

        for temporary, final in paths_to_publish:
            os.rename(temporary, final)

        # Quicklook is regenerated after publication.  SAR is displayed in dB
        # while the stored GeoTIFFs remain linear power.  With --lia yes the
        # second row contains LIA and ellipsoid incidence angle.  This also
        # refreshes an older 2-panel quicklook during geometry-only resume.
        try:
            make_quicklook(
                co_path,
                cross_path,
                quicklook_path,
                cluster[0].co_key,
                cluster[0].cross_key,
                normalization,
                incidence_path=incidence_path if make_lia else None,
                lia_path=lia_path if make_lia else None,
            )
        except Exception as error:
            log.warning(
                "GeoTIFFs completed, but quicklook creation failed: %s",
                error,
            )

        if not geometry_only:
            log.info("Saved: %s", co_path)
            log.info("Saved: %s", cross_path)
            log.info("Saved: %s", count_path)
        if make_lia:
            log.info("Saved: %s", incidence_path)
            log.info("Saved: %s", lia_path)
        if ease2_grid_mat_path is not None:
            log.info("EASE2 index companion: %s", ease2_grid_mat_path)

        status = "processed_geometry" if geometry_only else "processed"
        return status, expected_paths

    finally:
        if keep_intermediate:
            log.info("Keeping intermediate rasters: %s", work_dir)
        else:
            shutil.rmtree(work_dir, ignore_errors=True)
            try:
                if not os.listdir(work_root):
                    os.rmdir(work_root)
            except OSError:
                pass



def deduplicate_files(
    files: list[str],
    log: logging.Logger,
) -> list[str]:
    buckets: dict[tuple, list[tuple[str, dict]]] = defaultdict(list)

    for path in files:
        info = parse_filename(path)

        if info is None:
            log.warning(
                f"Skipping unrecognized filename: {os.path.basename(path)}"
            )
            continue

        identity = (
            info["instrument"],
            info["level"],
            info["proc_type"],
            info["cycle"],
            info["track"],
            info["pass"],
            info["frame"],
            info["bw_mode"],
            info["pole"],
            info["source"],
            info["date"],
            info["time"],
            info["end_date"],
            info["end_time"],
            info["accuracy"],
            info["coverage"],
            info["loc"],
        )

        buckets[identity].append((path, info))

    selected: list[str] = []

    for entries in buckets.values():
        entries.sort(
            key=lambda item: (
                int(item[1]["counter"]),
                item[1]["crid"],
                item[0],
            )
        )

        chosen_path, _ = entries[-1]
        selected.append(chosen_path)

        for discarded_path, _ in entries[:-1]:
            log.warning(
                f"Duplicate identity: keeping {os.path.basename(chosen_path)}; "
                f"discarding {os.path.basename(discarded_path)}"
            )

    return sorted(selected)


def group_key(info: dict) -> tuple:
    # CRID is included so products created by different processor releases are
    # not silently averaged together.
    return (
        info["instrument"],
        info["proc_type"],
        info["cycle"],
        info["track"],
        info["date"],
        info["pass"],
        info["bw_mode"],
        info["pole"],
        info["source"],
        info["crid"],
    )


def process_group_worker(
    payload: tuple,
) -> tuple[tuple, list[tuple[str, Optional[tuple[str, ...]]]]]:
    (
        key,
        file_list,
        output_dir,
        grid_mode,
        target_epsg,
        resolution,
        origin_x,
        origin_y,
        bbox_wgs84,
        frequency,
        pol_pair,
        mosaic,
        gap_m,
        output_tile_km,
        native_block_px,
        warp_threads,
        normalization,
        filter_name,
        filter_window,
        noise_floor,
        use_exception_mask,
        min_coverage,
        make_lia,
        nisar_dem_vrt,
        lia_dem_scale,
        dem_cache_dir,
        earthdata_cookie,
        keep_intermediate,
    ) = payload

    safe_key = "_".join(str(value) for value in key)
    log_path = os.path.join(
        output_dir,
        "logs",
        f"group_{safe_key}.log",
    )

    log = setup_logger(log_path)
    log.info(f"Processing group {key}: {len(file_list)} files")

    frames: list[FrameMeta] = []
    for path in file_list:
        try:
            frames.append(
                inspect_frame(
                    path,
                    frequency,
                    pol_pair,
                    target_epsg,
                )
            )
        except (KeyError, ValueError) as error:
            log.warning(
                f"Skipping incompatible file {os.path.basename(path)}: {error}"
            )

    if not frames:
        log.warning(
            f"No files in group {key} contain the requested polarization pair"
        )
        return key, []

    co_keys = {frame.co_key for frame in frames}
    cross_keys = {frame.cross_key for frame in frames}
    if len(co_keys) != 1 or len(cross_keys) != 1:
        raise ValueError(
            "Frames in one group do not contain the same covariance pair: "
            f"co={co_keys}, cross={cross_keys}"
        )

    clusters = cluster_frames(
        frames,
        gap_m=gap_m,
        mosaic=mosaic,
    )

    results: list[tuple[str, Optional[tuple[str, ...]]]] = []
    for index, cluster in enumerate(clusters, start=1):
        results.append(
            process_cluster(
                cluster=cluster,
                output_dir=output_dir,
                cluster_index=index,
                n_clusters=len(clusters),
                grid_mode=grid_mode,
                target_epsg=target_epsg,
                resolution=resolution,
                origin_x=origin_x,
                origin_y=origin_y,
                bbox_wgs84=bbox_wgs84,
                output_tile_km=output_tile_km,
                native_block_px=native_block_px,
                warp_threads=warp_threads,
                normalization=normalization,
                filter_name=filter_name,
                filter_window=filter_window,
                noise_floor=noise_floor,
                use_exception_mask=use_exception_mask,
                min_coverage=min_coverage,
                make_lia=make_lia,
                nisar_dem_vrt=nisar_dem_vrt,
                lia_dem_scale=lia_dem_scale,
                dem_cache_dir=dem_cache_dir,
                earthdata_cookie=earthdata_cookie,
                keep_intermediate=keep_intermediate,
                log=log,
            )
        )
        gc.collect()

    return key, results



def parse_bbox(
    values: Optional[list[float]],
) -> Optional[tuple[float, float, float, float]]:
    if values is None:
        return None

    west, south, east, north = map(float, values)

    if not (
        -180 <= west < east <= 180
        and -90 <= south < north <= 90
    ):
        raise ValueError(
            "--bbox must be WEST SOUTH EAST NORTH in EPSG:4326"
        )

    return west, south, east, north


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "NISAR GCOV -> fixed-grid regional GeoTIFF mosaics "
            "with safe skip/resume behavior"
        )
    )

    parser.add_argument(
        "--input_dir",
        required=True,
    )

    parser.add_argument(
        "--output_dir",
        required=True,
    )

    parser.add_argument(
        "--frequency",
        choices=["frequencyA", "frequencyB"],
        default="frequencyA",
    )

    parser.add_argument(
        "--pol_pair",
        choices=["auto", *POL_PAIRS.keys()],
        default="auto",
    )

    parser.add_argument(
        "--normalization",
        choices=["gamma0", "sigma0"],
        default="gamma0",
    )

    parser.add_argument(
        "--filter",
        choices=["none", "enhanced_lee", "gamma_map"],
        default="none",
    )

    parser.add_argument(
        "--filter_window",
        type=int,
        default=7,
    )

    parser.add_argument(
        "--noise_floor",
        type=float,
        default=0.0,
        help="Linear-power threshold; 0 disables thresholding (recommended)",
    )

    parser.add_argument(
        "--min_coverage",
        type=float,
        default=0.90,
        help=(
            "Minimum fraction of native area with mask in "
            "1..numberOfSubSwaths required in each output pixel; default 0.90"
        ),
    )

    parser.add_argument(
        "--use_exception_mask",
        choices=["yes", "no"],
        default="yes",
        help=(
            "yes (default): exclude pixels where inputDataExceptionMask is "
            "nonzero (e.g. NISAR sample slip/shift). no: keep those pixels "
            "instead of masking them, provided they pass every other "
            "validity check (subswath mask, numberOfLooks, finite/positive "
            "covariance)."
        ),
    )

    parser.add_argument(
        "--lia",
        choices=["yes", "no"],
        default="no",
        help=(
            "Generate ellipsoid incidence-angle and terrain-local-incidence-angle "
            "GeoTIFFs. Geometry is computed from the NISAR radarGrid + official "
            "NISAR-modified Copernicus DEM, aggregated to --scale, and masked to "
            "the same valid HH/HV (or selected co/cross) support. Default: no."
        ),
    )

    parser.add_argument(
        "--lia_dem_scale",
        type=float,
        default=30.0,
        help=(
            "Metric terrain-grid spacing used internally for DEM slope/LIA "
            "before aggregation to the final output grid. Default: 30 m."
        ),
    )

    parser.add_argument(
        "--nisar_dem_vrt",
        default=DEFAULT_NISAR_DEM_VRT,
        help=(
            "Official NISAR-modified Copernicus DEM global VRT URL, /vsi path, "
            "or a local DEM/VRT override."
        ),
    )

    parser.add_argument(
        "--dem_cache_dir",
        default=None,
        help=(
            "Persistent DEM cache directory. With --lia yes, an existing valid "
            "cached DEM is reused instead of downloading/reprojecting it again. "
            "Default: <input_dir>/DEM_Cache."
        ),
    )

    parser.add_argument(
        "--earthdata_cookie",
        default="/tmp/gdal_cookies.txt",
        help=(
            "Known-working Earthdata/GDAL cookie file used to seed a separate "
            "cookie jar for each parallel worker. Default: /tmp/gdal_cookies.txt. "
            "If the file is absent, workers fall back to a fresh cookie and ~/.netrc."
        ),
    )

    parser.add_argument(
        "--grid_mode",
        "--grid",
        dest="grid_mode",
        type=str.upper,
        choices=["UTM", "EASE2"],
        default="UTM",
        help=(
            "Destination grid family. UTM preserves the existing projected "
            "fixed-grid behavior using --target_epsg/--scale/--grid_origin_*. "
            "EASE2 forces the exact global EASE-Grid 2.0 M0p1 lattice "
            "(EPSG:6933, 100.089502334956 m, 347040 x 146160 globally) and "
            "crops that lattice to each scene/cluster. Default: UTM."
        ),
    )

    parser.add_argument(
        "--ease2_mat",
        default=None,
        help=(
            "Optional EASE2_latlon100m.mat file used only to validate that the "
            "ancillary lat/lon pixel centers match the built-in M0p1 lattice. "
            "Not required for processing."
        ),
    )

    parser.add_argument(
        "--target_epsg",
        default="auto",
        help=(
            "UTM-mode destination projected CRS. Default 'auto' reads all "
            "selected NISAR HDF5 source EPSGs and chooses the majority; a tie "
            "uses the first readable file among the tied EPSGs. Supply an "
            "integer such as 32611 to override. Ignored in EASE2 mode, where "
            "EPSG:6933 is fixed."
        ),
    )

    parser.add_argument(
        "--scale",
        type=float,
        default=100.0,
        help=(
            "UTM-mode output pixel size in metres (default 100). "
            "EASE2 mode always uses exact M0p1 spacing "
            "100.089502334956 m."
        ),
    )

    parser.add_argument(
        "--grid_origin_x",
        type=float,
        default=0.0,
    )

    parser.add_argument(
        "--grid_origin_y",
        type=float,
        default=0.0,
    )

    parser.add_argument(
        "--bbox",
        nargs=4,
        type=float,
        metavar=("WEST", "SOUTH", "EAST", "NORTH"),
        help="Optional clipping box in EPSG:4326",
    )

    parser.add_argument(
        "--mosaic",
        choices=["yes", "no"],
        default="yes",
    )

    parser.add_argument(
        "--gap_threshold",
        type=float,
        default=50.0,
        help="Maximum gap in km for adjacent frames to be in one mosaic",
    )

    parser.add_argument(
        "--tile_size_km",
        type=float,
        default=50.0,
        help=(
            "Fixed-grid output processing window size; does not alter values"
        ),
    )

    parser.add_argument(
        "--native_block_px",
        type=int,
        default=1024,
        help="Native-grid HDF5 read/write block size",
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel acquisition-group workers",
    )

    parser.add_argument(
        "--tile_workers",
        type=int,
        default=4,
        help="GDAL warp threads per acquisition group",
    )

    parser.add_argument(
        "--keep_intermediate",
        choices=["yes", "no"],
        default="no",
    )

    args = parser.parse_args()

    if args.filter_window < 3 or args.filter_window % 2 == 0:
        parser.error("--filter_window must be odd and >= 3")

    if args.scale <= 0 or args.tile_size_km <= 0:
        parser.error("--scale and --tile_size_km must be positive")

    if args.lia_dem_scale <= 0:
        parser.error("--lia_dem_scale must be positive")

    if args.native_block_px < 16:
        parser.error("--native_block_px must be >= 16")

    if args.workers < 1 or args.tile_workers < 1:
        parser.error("--workers and --tile_workers must be >= 1")

    if args.noise_floor < 0:
        parser.error("--noise_floor must be >= 0")

    if not (0.0 < args.min_coverage <= 1.0):
        parser.error("--min_coverage must be in (0, 1]")

    try:
        bbox_wgs84 = parse_bbox(args.bbox)
    except ValueError as error:
        parser.error(str(error))

    os.makedirs(
        args.output_dir,
        exist_ok=True,
    )

    dem_cache_dir: Optional[str] = None
    if args.lia == "yes":
        dem_cache_dir = (
            os.path.abspath(os.path.expanduser(args.dem_cache_dir))
            if args.dem_cache_dir
            else os.path.join(os.path.abspath(args.input_dir), "DEM_Cache")
        )
        os.makedirs(dem_cache_dir, exist_ok=True)

    earthdata_cookie: Optional[str] = None
    if args.lia == "yes" and args.earthdata_cookie:
        earthdata_cookie = os.path.abspath(
            os.path.expanduser(args.earthdata_cookie)
        )

    log = setup_logger(
        os.path.join(
            args.output_dir,
            "processing.log",
        )
    )

    input_dir = Path(args.input_dir)

    discovered = sorted(
        str(path)
        for path in input_dir.rglob("*")
        if path.is_file()
        and path.suffix.lower() in {".h5", ".hdf5"}
    )

    selected = deduplicate_files(
        discovered,
        log,
    )

    if not selected:
        raise SystemExit(
            "No parseable NISAR GCOV HDF5 files were found"
        )

    grid_mode = args.grid_mode.upper()

    if grid_mode == "EASE2":
        if str(args.target_epsg).strip().lower() not in {
            "auto",
            str(EASE2_M0P1_EPSG),
        }:
            parser.error(
                "--grid_mode EASE2 fixes the CRS to EPSG:6933; "
                "do not supply a different --target_epsg"
            )
        target_epsg = EASE2_M0P1_EPSG
        resolution = EASE2_M0P1_RESOLUTION_M
        origin_x = EASE2_M0P1_LEFT
        origin_y = EASE2_M0P1_BOTTOM
    else:
        try:
            target_epsg = resolve_target_epsg(
                selected,
                frequency=args.frequency,
                pol_pair=args.pol_pair,
                requested=args.target_epsg,
                log=log,
            )
        except Exception as error:
            parser.error(str(error))

        resolution = float(args.scale)
        origin_x = float(args.grid_origin_x)
        origin_y = float(args.grid_origin_y)

    try:
        target_crs = rasterio.CRS.from_epsg(target_epsg)
    except Exception as error:
        parser.error(
            f"Invalid resolved target EPSG {target_epsg}: {error}"
        )

    if not target_crs.is_projected:
        parser.error(
            f"Resolved target EPSG:{target_epsg} is not projected. "
            "Do not use EPSG:4326 with a pixel size expressed in metres."
        )

    linear_units = (
        target_crs.linear_units or ""
    ).lower()

    if "metre" not in linear_units and "meter" not in linear_units:
        parser.error(
            f"Resolved target EPSG:{target_epsg} uses units "
            f"'{target_crs.linear_units}'. Choose a projected CRS whose "
            "horizontal units are metres."
        )

    if grid_mode == "EASE2" and args.ease2_mat:
        ease2_mat_path = os.path.abspath(os.path.expanduser(args.ease2_mat))
        if not os.path.isfile(ease2_mat_path):
            parser.error(f"--ease2_mat file not found: {ease2_mat_path}")
        try:
            validate_ease2_mat(ease2_mat_path, log)
        except Exception as error:
            parser.error(str(error))

    log.info("=" * 80)
    log.info("NISAR GCOV generic fixed-grid processor")
    log.info("SAFE RESUME MODE  : existing GeoTIFFs are never overwritten")
    log.info(f"Input directory   : {args.input_dir}")
    log.info(f"Output directory  : {args.output_dir}")
    log.info(f"Grid mode         : {grid_mode}")
    if grid_mode == "EASE2":
        log.info(
            "Destination grid  : EASE-Grid 2.0 M0p1, "
            f"EPSG:{target_epsg}, {resolution:.12f} m"
        )
        log.info(
            "Global M0p1 grid  : "
            f"{EASE2_M0P1_ROWS} rows x {EASE2_M0P1_COLS} cols"
        )
        log.info(
            "Global grid edges : "
            f"left={EASE2_M0P1_LEFT:.6f}, "
            f"right={EASE2_M0P1_RIGHT:.6f}, "
            f"top={EASE2_M0P1_TOP:.6f}, "
            f"bottom={EASE2_M0P1_BOTTOM:.6f}"
        )
    else:
        log.info(f"EPSG selection    : {args.target_epsg}")
        log.info(
            f"Destination grid  : EPSG:{target_epsg}, {resolution:g} m"
        )
        log.info(
            f"Grid origin       : ({origin_x:g}, {origin_y:g})"
        )
    log.info(f"CRS definition    : {target_crs.to_string()}")
    log.info(
        f"Mosaic            : {args.mosaic}; "
        f"adjacency gap={args.gap_threshold:g} km"
    )
    log.info(
        f"Frequency/pair    : {args.frequency} / {args.pol_pair}"
    )
    log.info(
        f"Radiometry        : {args.normalization} linear power"
    )
    log.info(
        f"Filter            : {args.filter}; "
        f"window={args.filter_window}"
    )
    log.info(
        f"Noise threshold   : {args.noise_floor:g} linear power"
    )
    log.info(
        f"Min. coverage     : {args.min_coverage:.2f} of each output pixel"
    )
    log.info(
        f"Exception mask    : "
        f"{'applied (nonzero pixels excluded)' if args.use_exception_mask == 'yes' else 'NOT applied (flagged pixels kept)'}"
    )
    log.info(
        f"LIA/incidence     : {args.lia}"
    )
    if args.lia == "yes":
        log.info(
            f"LIA DEM scale     : {args.lia_dem_scale:g} m -> "
            f"final {resolution:.12g} m grid"
        )
        log.info(
            "Geometry masking  : same valid co/cross SAR support + min_coverage"
        )
        log.info(
            f"DEM cache         : {dem_cache_dir}"
        )
        log.info(
            f"Earthdata cookie  : {earthdata_cookie}"
        )
    log.info(
        f"Parallelism       : groups={args.workers}; "
        f"warp threads/group={args.tile_workers}"
    )

    if bbox_wgs84:
        log.info(
            f"Clip bbox         : {bbox_wgs84} EPSG:4326"
        )

    log.info(
        (
            "Fixed grid definition: exact EASE2 M0p1 global lattice"
            if grid_mode == "EASE2"
            else "Fixed grid definition: resolved target EPSG, scale, and grid origin"
        )
    )
    log.info(f"Selected {len(selected)} input file(s)")
    log.info("=" * 80)

    groups: dict[tuple, list[str]] = defaultdict(list)

    for path in selected:
        info = parse_filename(path)
        assert info is not None
        groups[group_key(info)].append(path)

    for key, paths in sorted(groups.items()):
        frames = compress_frames(
            parse_filename(path)["frame"]
            for path in paths
        )
        log.info(
            f"Group {key}: files={len(paths)}, frames={frames}"
        )

    payloads = [
        (
            key,
            sorted(paths),
            args.output_dir,
            grid_mode,
            target_epsg,
            resolution,
            origin_x,
            origin_y,
            bbox_wgs84,
            args.frequency,
            args.pol_pair,
            args.mosaic == "yes",
            args.gap_threshold * 1000.0,
            args.tile_size_km,
            args.native_block_px,
            args.tile_workers,
            args.normalization,
            args.filter,
            args.filter_window,
            args.noise_floor,
            args.use_exception_mask == "yes",
            args.min_coverage,
            args.lia == "yes",
            args.nisar_dem_vrt,
            args.lia_dem_scale,
            dem_cache_dir,
            earthdata_cookie,
            args.keep_intermediate == "yes",
        )
        for key, paths in sorted(groups.items())
    ]

    results: list[
        tuple[
            tuple,
            list[tuple[str, Optional[tuple[str, ...]]]],
        ]
    ] = []

    if args.workers == 1:
        for payload in payloads:
            results.append(
                process_group_worker(payload)
            )

    else:
        with ProcessPoolExecutor(
            max_workers=args.workers
        ) as pool:
            future_map = {
                pool.submit(
                    process_group_worker,
                    payload,
                ): payload[0]
                for payload in payloads
            }

            for future in as_completed(future_map):
                key = future_map[future]

                try:
                    results.append(
                        future.result()
                    )
                    log.info(
                        f"Completed group {key}"
                    )
                except Exception:
                    log.exception(
                        f"Failed group {key}"
                    )
                    raise

    status_counts = Counter()

    for _, cluster_results in results:
        for status, _ in cluster_results:
            status_counts[status] += 1

    processed_count = status_counts["processed"]
    processed_geometry_count = status_counts["processed_geometry"]
    skipped_existing_count = status_counts["skipped_existing"]
    skipped_mismatch_count = status_counts["skipped_existing_mismatch"]
    skipped_partial_count = status_counts["skipped_partial"]
    skipped_race_count = status_counts["skipped_race_existing"]

    total_clusters = sum(status_counts.values())

    log.info("=" * 80)
    log.info("PROCESSING SUMMARY")
    log.info(f"Acquisition groups examined : {len(results)}")
    log.info(f"Clusters examined           : {total_clusters}")
    log.info(f"New clusters processed      : {processed_count}")
    log.info(f"Geometry-only additions     : {processed_geometry_count}")
    log.info(f"Complete existing skipped   : {skipped_existing_count}")
    log.info(f"Existing mismatch skipped   : {skipped_mismatch_count}")
    log.info(f"Partial existing skipped    : {skipped_partial_count}")
    log.info(f"Race/external output skipped: {skipped_race_count}")
    log.info(
        f"Finished: {datetime.now().isoformat(timespec='seconds')}"
    )
    log.info("=" * 80)


if __name__ == "__main__":
    main()