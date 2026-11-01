# NISAR GCOV Processing

**Native-grid speckle filtering, 100 m EASE-Grid 2.0 export, and MAAP DPS execution**

This repository converts **NISAR L2 GCOV HDF5** granules into spatially aligned GeoTIFFs. It supports processing one granule per **NASA MAAP Distributed Processing System (DPS)** job or processing a local directory of granules. The workflow is geographically reusable: the area of interest is supplied **when a job is submitted**, not fixed during algorithm registration.

The default **MAAP DPS wrapper** applies **7 × 7 Enhanced Lee speckle filtering** at native GCOV grid spacing **before additional spatial averaging to the ~100 m output grid**. GCOV products have *already undergone multilooking during NISAR product generation*; this filter does **not** operate on pre-multilooked SLC data.

## Processing overview

```text
L2 GCOV HDF5 (one granule per MAAP DPS job)
  └─► Read native-grid co-/cross-pol covariance powers and numberOfLooks
       └─► Apply GCOV validity and exception masks
            └─► Enhanced Lee speckle filter, 7 × 7 (native grid, linear power)
                 └─► Area-average/reproject onto fixed EASE2 M0p1 (~100 m)
                      └─► Co-pol and cross-pol GeoTIFFs (linear gamma0)
                      └─► Incidence angle / local incidence angle (optional)
                      └─► Frame-count raster, EASE2 grid companion, quicklook
```

The main processor also supports **no filter** and **Gamma MAP**, and a **UTM** destination grid for local/HPC runs. Mosaicking of compatible adjacent frames is supported when enabled; the provided MAAP wrapper deliberately uses `--mosaic no` and processes **one input granule per job**.

### Default MAAP DPS settings

| Setting | Current DPS wrapper behavior |
| --- | --- |
| NISAR input | L2 GCOV HDF5 from the specified S3 URL |
| Frequency | `frequencyA` |
| Polarizations | `auto` — first available supported co-/cross-pol pair (`HH_HV`, `VV_VH`, `RH_RV`) |
| Radiometry | `gamma0`, output stored as **linear power**, not dB |
| Speckle filter | `enhanced_lee` |
| Filter window | `7` (7 × 7 native pixels) |
| Grid | Global EASE2 M0p1, **EPSG:6933**, pixel spacing **100.089502334956 m** |
| Incidence / LIA | Enabled; LIA uses the NISAR DEM workflow |
| Frame mosaic | Disabled (`--mosaic no`) |
| Intermediate rasters | Removed after successful processing |

**Important:** `--pol_pair auto` selects a pair *only if both channels are available*. It does not fabricate a missing channel or automatically fall back from `frequencyA` to `frequencyB`. Granules without a supported pair in frequency A cannot produce a co-/cross-pol output using the current DPS defaults.

## Repository files

| File | Purpose |
| --- | --- |
| `NISAR_GCOV_TO_TIFF_E.py` | Main HDF5 → GeoTIFF processor, masks, filters, grid averaging, LIA, mosaics |
| `NISAR_GCOV_TIFF_E.sh` | Command-line launcher for local or HPC processing |
| `nisar_dps_wrapper.py` | MAAP DPS wrapper: obtains credentials, downloads one granule, runs processor |
| `run_nisar_dps.sh` | MAAP runtime entry point accepting S3 URL and AOI coordinates |
| `build_env.sh` | MAAP build-time dependency check and installation |
| `.github/workflows/python-package-conda.yml` | GitHub Actions **syntax/configuration checks only**; not a MAAP container build |

## MAAP DPS registration

Register the repository in MAAP **Algorithm Catalog** using settings appropriate to your MAAP workspace:

| Registration field | Value |
| --- | --- |
| Git repository | `https://github.com/preet0lal/NISAR_GCOV_Processing.git` |
| Branch | `main` |
| Algorithm name / version | `nisar-gcov-ease2` / `main` |
| **Build Command** | `NISAR_GCOV_Processing/build_env.sh` |
| **Run Command** | `NISAR_GCOV_Processing/run_nisar_dps.sh` |
| Container image | Use a compatible MAAP base image with GDAL available |
| Resources | Select resources appropriate for the granule and LIA workload |

### Required positional inputs

Create **five positional inputs, in this exact order**, with **no geographic default values**:

| Position | Input name | Meaning |
| --- | --- | --- |
| 1 | `s3_url` | Complete NISAR GCOV HDF5 `s3://` URL |
| 2 | `bbox_w` | Western longitude (degrees east) |
| 3 | `bbox_s` | Southern latitude (degrees north) |
| 4 | `bbox_e` | Eastern longitude (degrees east) |
| 5 | `bbox_n` | Northern latitude (degrees north) |

All inputs are supplied when submitting a job, allowing the **same registration** to process India, CONUS, Australia, or another region. The longitude/latitude bounding box is an axis-aligned selection/clip envelope; it is **not** a country-boundary mask.

A conceptual MAAP submission (from an authenticated MAAP notebook) is:

```python
from maap.maap import MAAP

maap = MAAP()
# Populate these values from the granule you want to process and your AOI:
granule_s3_url = "s3://BUCKET/PREFIX/GRANULE.h5"
west, south, east, north = (-10.0, -10.0, 10.0, 10.0)  # example placeholders only

# Supply s3_url, bbox_w, bbox_s, bbox_e, bbox_n in the same
# positional order configured during MAAP algorithm registration.
# Submit via the MAAP submitJob interface used by your batch notebook.
```

**Deployment rule:** Updating GitHub does **not** update an already-built MAAP DPS image. After changing the wrapper or build script, **rebuild/re-register the algorithm**, then submit a single test granule before starting a large batch. Previously completed jobs do not change.

### Dependency installation

The MAAP **Build Command** runs `build_env.sh`, which verifies `maap-py`, `s3fs`, `numpy`, `h5py`, `rasterio`, `scipy`, and `matplotlib`, installing missing Python packages with `pip`. It warns if `gdalwarp` is unavailable; **LIA processing requires the GDAL CLI**.

An immediate failure such as:

```text
ModuleNotFoundError: No module named 's3fs'
```

means the active MAAP runtime environment does not have the download dependency. Rebuild with the Build Command set, and verify that the dependency check passes in the build log. A successful GitHub Actions check by itself does **not** validate the MAAP container.

## Local / HPC processing

For locally accessible GCOV HDF5 files, use the shell launcher directly. In contrast to the MAAP wrapper, this launcher defaults to **no filtering** unless a mode is specified.

**Enhanced Lee (7 × 7), EASE2, 100 m, single-frame processing:**

```bash
bash NISAR_GCOV_TIFF_E.sh \
  --input_dir /path/to/gcov_h5 \
  --output_dir /path/to/output_enhanced_lee \
  --grid_mode EASE2 \
  --frequency frequencyA \
  --pol_pair auto \
  --normalization gamma0 \
  --mode enhanced_lee \
  --window 7 \
  --mosaic no \
  --lia yes
```

**Unfiltered baseline:**

```bash
bash NISAR_GCOV_TIFF_E.sh \
  --input_dir /path/to/gcov_h5 \
  --output_dir /path/to/output_nofilt \
  --grid_mode EASE2 \
  --mode nofilt \
  --mosaic no
```

Other shell launcher modes: `--mode gamma_map` and `--mode all` (launches separate jobs for no filter, Enhanced Lee, and Gamma MAP). For a local UTM 100 m grid, select `--grid_mode UTM --target_epsg auto --scale 100`. Run `bash NISAR_GCOV_TIFF_E.sh --help` for the full option list.

## Outputs and QA

For Enhanced Lee with a 7 × 7 window, representative file suffixes are:

```text
*_100m_enhanced_lee_w7_HHHH.tif       # Co-pol example (HH)
*_100m_enhanced_lee_w7_HVHV.tif       # Cross-pol example (HV)
*_100m_enhanced_lee_w7_frame_count.tif
*_100m_enhanced_lee_w7_incidence_angle_deg.tif
*_100m_enhanced_lee_w7_local_incidence_angle_deg.tif
*_100m_enhanced_lee_w7_EASE2_grid.mat
*_100m_enhanced_lee_w7_quicklook.png
```

When the selected polarization pair is VV/VH or RH/RV, the channel-specific filenames differ. SAR GeoTIFF pixels are **linear power**; the quicklook displays SAR backscatter in **dB**.

**Check filtering before running a regional batch.** In the MAAP job stdout, look for:

```text
Filter            : enhanced_lee; window=7
```

and in the executed command:

```text
--filter enhanced_lee --filter_window 7
```

If the log says `Filter : none` or outputs contain `_nofilt_`, the active **MAAP deployment is still running unfiltered code**, even if the job status is `COMPLETED`.

### Filtering and validation notes

- Filtering occurs on the **native geocoded GCOV grid**, after upstream NISAR multilooking but before the processor's additional 100 m area averaging.
- The filters operate in **linear power**, not dB. The Enhanced Lee model uses the product's spatially varying `numberOfLooks`.
- The processor excludes invalid/missing GCOV samples and, by default, flags from `inputDataExceptionMask`. It requires at least **90% valid source-area coverage** by default when creating target-grid SAR outputs.
- Adjacent 100 m cells may become more correlated after native-resolution filtering. For biomass applications, compare a **no-filter baseline** with filtered results using independent validation rather than assuming greater filtering always improves AGBD accuracy.
- The MAAP wrapper's one-granule-per-job outputs are **not** automatically combined across dates, tracks, or ascending/descending passes. Mosaicking or subsequent retrieval/fusion is a separate step.

## Troubleshooting

| Symptom | Likely cause / action |
| --- | --- |
| `No module named 's3fs'` | MAAP image lacks dependency: use `build_env.sh` as Build Command, then rebuild/re-register |
| Job completes but outputs say `nofilt` | Old deployment or wrong wrapper: inspect printed `--filter` command and rebuild |
| No compatible polarizations | `frequencyA` has no supported complete co-/cross-pol pair; inspect GCOV channels/frequency |
| LIA fails at DEM stage | Check `gdalwarp`, GDAL credentials, DEM access, and the job stderr |
| Outputs already exist | Processor intentionally skips complete output sets rather than overwriting |
| GitHub Actions passed but MAAP fails | GitHub Actions checks syntax/configuration; MAAP builds/runs in a separate environment |

## License

[MIT](LICENSE). Copyright © 2026 Preet Lal.
