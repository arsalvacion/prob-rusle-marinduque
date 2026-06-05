# Probabilistic RUSLE: Multi-factor Sensitivity Analysis

Companion code for:

> Salvacion, A.R. (2026). Probabilistic soil erosion modeling with RUSLE and Monte Carlo simulation: quantifying and contextualizing rainfall-driven uncertainty in a tropical island landscape. *Submitted*.

## Overview

This repository contains the R code for the multi-factor Monte Carlo sensitivity analysis described in Sections 2.3.5, 3.5, and 4.5 of the manuscript. The script extends the primary R-only Monte Carlo analysis by simultaneously varying the rainfall-runoff erosivity (R), soil erodibility (K), vegetation cover (C), and conservation practice (P) factors to quantify total RUSLE prediction uncertainty and evaluate the robustness of hotspot classification under combined input uncertainty.

## Method

In each of 1,000 Monte Carlo iterations, the script:

1. Randomly samples a year (with replacement) from the 61-year TerraClimate record (1960–2020) to obtain the annual precipitation surface and corresponding R factor.
2. Draws independent scalar perturbation factors for K (±20%), C (±30%), and P (±20%) from uniform distributions.
3. Computes RUSLE soil loss: A = R × (K × k_mult) × LS × (C × c_mult) × (P × p_mult).
4. Updates running statistics using Welford's online algorithm (memory-efficient; does not store all 1,000 rasters).

The perturbation ranges are justified by:

- **K ±20%**: typical within-class variability in soil erodibility (Tian et al., 2025)
- **C ±30%**: uncertainty in look-up table assignment and temporal land cover change (Alewell et al., 2019)
- **P ±20%**: limited empirical basis for conservation practice values in Philippine conditions (Delgado and Canters, 2012)

## Requirements

- R ≥ 4.3
- Required packages: `terra` (≥ 1.7), `sf`, `tidyverse`, `tidyterra`, `patchwork`, `ggspatial`, `RColorBrewer`

Install all dependencies:

```r
install.packages(c("terra", "sf", "tidyverse", "tidyterra",
                    "patchwork", "ggspatial", "RColorBrewer"))
```

## Directory Structure

The script expects the following directory layout. Edit the paths in Section 1 (Configuration) of `sensitivity_analysis.R` to match your local setup.

```
project/
├── sensitivity_analysis.R
├── data/
│   ├── boundary/
│   │   └── MarProv.shp           # Province boundary
│   ├── municipalities/
│   │   └── Marmun.shp            # Municipality boundaries
│   ├── soil/
│   │   └── marinduque_K.tif      # K factor (30-m)
│   ├── topography/
│   │   └── marinduque_LS.tif     # LS factor (30-m)
│   ├── landcover/
│   │   ├── marinduque_C.tif      # C factor (30-m)
│   │   └── marinduque_P.tif      # P factor (30-m)
│   └── RUSLE Input/
│       └── A_deterministic.tif   # Deterministic RUSLE output
├── output/
│   └── rasters/
│       ├── A_mc_mean.tif         # R-only MC mean (from primary analysis)
│       ├── A_mc_sd.tif           # R-only MC standard deviation
│       ├── A_mc_cv.tif           # R-only MC coefficient of variation
│       ├── P_exceed_80.tif       # R-only exceedance probability (>80)
│       ├── ppt_ds_1960.tif       # Downscaled precipitation (1960)
│       ├── ppt_ds_1961.tif       # Downscaled precipitation (1961)
│       ├── ...                   #   ... (61 annual files)
│       └── ppt_ds_2020.tif       # Downscaled precipitation (2020)
└── output/sensitivity/           # Created by the script
    ├── rasters/
    │   ├── A_sens_mean.tif
    │   ├── A_sens_sd.tif
    │   ├── A_sens_cv.tif
    │   ├── P_sens_exceed_10.tif
    │   ├── P_sens_exceed_20.tif
    │   ├── P_sens_exceed_40.tif
    │   ├── P_sens_exceed_80.tif
    │   └── divergence_sensitivity.tif
    ├── tables/
    │   ├── sensitivity_province_summary.csv
    │   ├── sensitivity_municipal_summary.csv
    │   └── perturbation_log.csv
    └── figures/
        ├── fig_cv_comparison.png
        ├── fig_cv_histogram.png
        └── fig_divergence_sensitivity.png
```

## Input Data

| File | Description | Source |
|---|---|---|
| `MarProv.shp` | Province boundary | NAMRIA Geoportal Philippines |
| `Marmun.shp` | Municipality boundaries | NAMRIA Geoportal Philippines |
| `marinduque_K.tif` | Soil erodibility factor | BSWM soil map; values from David (1988) |
| `marinduque_LS.tif` | Slope length and steepness | ASTER GDEM 30-m; computed with SAGA GIS |
| `marinduque_C.tif` | Vegetation cover factor | NAMRIA land cover; values from David (1987) |
| `marinduque_P.tif` | Conservation practice factor | NAMRIA land cover; values from Delgado and Canters (2012) |
| `A_deterministic.tif` | Deterministic RUSLE soil loss | Computed from 61-year mean precipitation |
| `A_mc_*.tif`, `P_exceed_80.tif` | R-only Monte Carlo outputs | From primary analysis |
| `ppt_ds_YYYY.tif` | Downscaled annual precipitation (61 files) | TerraClimate; downscaled with random forest |

## Usage

1. Place input data in the directory structure shown above (or edit the paths in Section 1 of the script).
2. Run the primary R-only Monte Carlo analysis first to generate the comparison outputs (`A_mc_mean.tif`, `A_mc_sd.tif`, `A_mc_cv.tif`, `P_exceed_80.tif`).
3. Run the sensitivity analysis:

```r
source("sensitivity_analysis.R")
```

The script takes approximately 30–60 minutes depending on hardware (1,000 iterations over a ~96,000 ha raster at 30-m resolution). Progress is printed every 100 iterations.

## Outputs

### Rasters

| File | Description |
|---|---|
| `A_sens_mean.tif` | Multi-factor Monte Carlo mean soil loss (t ha⁻¹ yr⁻¹) |
| `A_sens_sd.tif` | Multi-factor Monte Carlo standard deviation |
| `A_sens_cv.tif` | Multi-factor Monte Carlo coefficient of variation (%) |
| `P_sens_exceed_{10,20,40,80}.tif` | Exceedance probabilities for erosion thresholds |
| `divergence_sensitivity.tif` | Divergence classification at >80 t ha⁻¹ yr⁻¹ threshold |

Divergence categories: 1 = confirmed severe, 2 = uncertain hotspot, 3 = hidden hotspot, 4 = confirmed not severe.

### Tables

| File | Description |
|---|---|
| `sensitivity_province_summary.csv` | Province-wide comparison of R-only vs. multi-factor statistics |
| `sensitivity_municipal_summary.csv` | Municipal-level comparison for six municipalities |
| `perturbation_log.csv` | Record of sampled years and perturbation factors for all 1,000 iterations |

### Figures

| File | Description | Manuscript reference |
|---|---|---|
| `fig_cv_comparison.png` | Side-by-side CV maps (R-only vs. multi-factor) | Fig. 9a,b |
| `fig_divergence_sensitivity.png` | Divergence map under multi-factor uncertainty | Fig. 9c |
| `fig_cv_histogram.png` | Pixel-level CV distribution comparison | Fig. 9d |

## Reproducibility

The script uses `set.seed(2024)` for full reproducibility. All 1,000 sampled years and perturbation factors are logged in `perturbation_log.csv`.

## Citation

If you use this code, please cite:

```
Salvacion, A.R. (2026). Probabilistic soil erosion modeling with RUSLE and Monte Carlo
simulation: quantifying and contextualizing rainfall-driven uncertainty in a tropical
island landscape. SUBMITTED. https://doi.org/[DOI]
```

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contact

Arnold R. Salvacion
Curtin Biometry and Agricultural Data Analytics
Centre for Crop and Disease Management
Curtin University, Bentley 6102, Australia
arnold.salvacion@curtin.edu.au
