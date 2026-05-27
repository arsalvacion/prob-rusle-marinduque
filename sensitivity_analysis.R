# ============================================================
# Multi-factor Monte Carlo sensitivity analysis for RUSLE
# ============================================================
#
# Companion code for:
#   Salvacion, A.R. (2026). Probabilistic soil erosion modeling
#   with RUSLE and Monte Carlo simulation: quantifying and
#   contextualizing rainfall-driven uncertainty in a tropical
#   island landscape. CATENA.
#
# This script extends the primary R-only Monte Carlo analysis
# by simultaneously varying K, C, P, and R to quantify total
# RUSLE prediction uncertainty and evaluate the robustness of
# hotspot classification under combined input uncertainty.
#
# Perturbation ranges:
#   R — sampled from the empirical 61-year rainfall record
#   K — +/-20% uniform (Tian et al., 2025)
#   C — +/-30% uniform (Alewell et al., 2019)
#   P — +/-20% uniform (Delgado and Canters, 2012)
#
# Perturbation design: scalar (one random multiplier per factor
# per iteration, applied uniformly across the map). This captures
# systematic uncertainty in look-up table assignment.
#
# Requirements:
#   R >= 4.3, terra >= 1.7, tidyverse, sf, patchwork, ggspatial,
#   RColorBrewer, tidyterra
#
# Author: Arnold R. Salvacion
# Contact: arsalvacion@up.edu.ph
# License: MIT
# ============================================================


# ── 0. PACKAGES ──────────────────────────────────────────────

library(terra)
library(sf)
library(tidyverse)
library(tidyterra)
library(patchwork)
library(ggspatial)
library(RColorBrewer)


# ── 1. CONFIGURATION ────────────────────────────────────────
#
# Edit the paths below to match your local directory structure.
# All other parameters can be adjusted in this section.
# ─────────────────────────────────────────────────────────────

# --- Paths ---
data_dir   <- "data"              # Root directory for input data
output_dir <- "output"            # Root directory for primary MC outputs
sens_dir   <- "output/sensitivity"  # Sensitivity analysis outputs

# Input files (relative to data_dir)
boundary_path      <- file.path(data_dir, "boundary", "MarProv.shp")
municipality_path  <- file.path(data_dir, "municipalities", "Marmun.shp")
K_path             <- file.path(data_dir, "soil", "marinduque_K.tif")
LS_path            <- file.path(data_dir, "topography", "marinduque_LS.tif")
C_path             <- file.path(data_dir, "landcover", "marinduque_C.tif")
P_path             <- file.path(data_dir, "landcover", "marinduque_P.tif")
A_det_path         <- file.path(data_dir, "RUSLE Input", "A_deterministic.tif")

# Primary analysis outputs (for comparison; from the R-only MC)
A_mc_mean_path     <- file.path(output_dir, "rasters", "A_mc_mean.tif")
A_mc_sd_path       <- file.path(output_dir, "rasters", "A_mc_sd.tif")
A_mc_cv_path       <- file.path(output_dir, "rasters", "A_mc_cv.tif")
P_exceed_80_path   <- file.path(output_dir, "rasters", "P_exceed_80.tif")

# Downscaled precipitation: individual annual TIFs (ppt_ds_YYYY.tif)
ppt_pattern        <- "ppt_ds_\\d{4}\\.tif$"
ppt_dir            <- file.path(output_dir, "rasters")

# --- Simulation parameters ---
n_iter  <- 1000L                  # Monte Carlo iterations
years   <- 1960:2020              # TerraClimate record span
seed    <- 2024L                  # Random seed for reproducibility

# --- Perturbation ranges (proportion) ---
K_range <- 0.20                   # K +/- 20%
C_range <- 0.30                   # C +/- 30%
P_range <- 0.20                   # P +/- 20%

# --- Erosion thresholds (Singh et al., 1992) ---
thresholds <- c(10, 20, 40, 80)  # t ha-1 yr-1

# --- Municipality name field in shapefile ---
mun_name_field <- "NAME_2"


# ── 2. SETUP ────────────────────────────────────────────────

set.seed(seed)

dir.create(file.path(sens_dir, "rasters"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(sens_dir, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(sens_dir, "figures"), recursive = TRUE, showWarnings = FALSE)


# ── 3. LOAD DATA ────────────────────────────────────────────

cat("Loading input data...\n")

boundary       <- vect(boundary_path)
municipalities <- vect(municipality_path)

K  <- rast(K_path)
LS <- rast(LS_path)
C  <- rast(C_path)
P  <- rast(P_path)

ref <- LS

# Align grids
if (!compareGeom(K, ref, stopOnError = FALSE)) K <- resample(K, ref, method = "near")
if (!compareGeom(C, ref, stopOnError = FALSE)) C <- resample(C, ref, method = "near")
if (!compareGeom(P, ref, stopOnError = FALSE)) P <- resample(P, ref, method = "near")

A_deterministic     <- rast(A_det_path)
A_mc_mean_primary   <- rast(A_mc_mean_path)
A_mc_sd_primary     <- rast(A_mc_sd_path)
A_mc_cv_primary     <- rast(A_mc_cv_path)
P_exceed_80_primary <- rast(P_exceed_80_path)

# Load downscaled precipitation stack
ppt_files <- list.files(ppt_dir, pattern = ppt_pattern, full.names = TRUE) |> sort()

if (length(ppt_files) != length(years)) {
  stop("Expected ", length(years), " precipitation files but found ", length(ppt_files),
       ".\nCheck ppt_dir and ppt_pattern in the configuration section.")
}

ppt_stack <- rast(ppt_files)
names(ppt_stack) <- paste0("ppt_", years)

cat("  Precipitation stack:", nlyr(ppt_stack), "layers\n")
cat("  Grid dimensions:    ", nrow(ref), "x", ncol(ref), "\n")


# ── 4. R-FACTOR STACK ───────────────────────────────────────
# R = 38.5 + 0.38 * P  (Eq. 2 in manuscript)

cat("Computing R-factor stack...\n")
R_stack <- 38.5 + 0.38 * ppt_stack
names(R_stack) <- paste0("R_", years)


# ── 5. MONTE CARLO SIMULATION ───────────────────────────────
#
# Welford's online algorithm computes running mean and variance
# without storing all 1,000 rasters in memory.

cat("Running Monte Carlo simulation (n = ", n_iter, ")...\n")

sampled_years  <- sample(years, size = n_iter, replace = TRUE)
year_to_layer  <- setNames(seq_along(years), as.character(years))

# Accumulators
mc_mean <- ref * 0
mc_M2   <- ref * 0

exceed_n <- setNames(
  lapply(thresholds, function(x) ref * 0),
  as.character(thresholds)
)

# Perturbation log
pert_log <- tibble(
  iteration = integer(),
  year      = integer(),
  k_mult    = numeric(),
  c_mult    = numeric(),
  p_mult    = numeric()
)

t0 <- Sys.time()

for (i in seq_len(n_iter)) {

  # Sample rainfall year
  yr  <- sampled_years[i]
  R_i <- R_stack[[year_to_layer[as.character(yr)]]]

  # Draw scalar perturbation factors
  k_mult <- runif(1, 1 - K_range, 1 + K_range)
  c_mult <- runif(1, 1 - C_range, 1 + C_range)
  p_mult <- runif(1, 1 - P_range, 1 + P_range)

  # RUSLE with perturbed factors
  A_i <- R_i * (K * k_mult) * LS * (C * c_mult) * (P * p_mult)

  # Welford update
  delta   <- A_i - mc_mean
  mc_mean <- mc_mean + delta / i
  delta2  <- A_i - mc_mean
  mc_M2   <- mc_M2 + delta * delta2

  # Exceedance counts
  for (thr in thresholds) {
    exceed_n[[as.character(thr)]] <- exceed_n[[as.character(thr)]] +
      ifel(A_i > thr, 1, 0)
  }

  # Log
  pert_log <- bind_rows(pert_log, tibble(
    iteration = i, year = yr,
    k_mult = k_mult, c_mult = c_mult, p_mult = p_mult
  ))

  # Progress
  if (i %% 100 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("  %4d / %d  [%.1f min elapsed, ~%.1f min remaining]\n",
                i, n_iter, elapsed, (n_iter - i) / (i / elapsed)))
  }
}

cat("Simulation complete.\n")


# ── 6. SUMMARY RASTERS ──────────────────────────────────────

mc_var <- mc_M2 / (n_iter - 1)
mc_sd  <- sqrt(mc_var)
mc_cv  <- ifel(mc_mean > 0, (mc_sd / mc_mean) * 100, NA)

exceed_prob <- lapply(exceed_n, function(x) x / n_iter)

names(mc_mean) <- "mean"
names(mc_sd)   <- "sd"
names(mc_cv)   <- "cv"


# ── 7. DIVERGENCE ANALYSIS ──────────────────────────────────
# Compare deterministic classification against multi-factor
# probabilistic at the Very Severe threshold (>80 t/ha/yr).
#
# Categories:
#   1 = Confirmed severe     (det >= 80 AND P(A>80) >= 0.50)
#   2 = Uncertain hotspot    (det >= 80 AND P(A>80) <  0.50)
#   3 = Hidden hotspot       (det <  80 AND P(A>80) >= 0.50)
#   4 = Confirmed not severe (det <  80 AND P(A>80) <  0.50)

det_class <- classify(
  A_deterministic,
  rcl = matrix(c(-Inf, 80, 0, 80, Inf, 1), ncol = 3, byrow = TRUE)
)

P80 <- exceed_prob[["80"]]

divergence <- ifel(
  det_class == 1 & P80 >= 0.50, 1,
  ifel(det_class == 1 & P80 < 0.50, 2,
       ifel(det_class == 0 & P80 >= 0.50, 3, 4))
)


# ── 8. SAVE OUTPUTS ─────────────────────────────────────────

cat("Saving outputs...\n")

writeRaster(mc_mean, file.path(sens_dir, "rasters", "A_sens_mean.tif"), overwrite = TRUE)
writeRaster(mc_sd,   file.path(sens_dir, "rasters", "A_sens_sd.tif"),   overwrite = TRUE)
writeRaster(mc_cv,   file.path(sens_dir, "rasters", "A_sens_cv.tif"),   overwrite = TRUE)

for (thr in thresholds) {
  writeRaster(exceed_prob[[as.character(thr)]],
              file.path(sens_dir, "rasters", paste0("P_sens_exceed_", thr, ".tif")),
              overwrite = TRUE)
}

writeRaster(divergence,
            file.path(sens_dir, "rasters", "divergence_sensitivity.tif"),
            overwrite = TRUE)

write.csv(pert_log,
          file.path(sens_dir, "tables", "perturbation_log.csv"),
          row.names = FALSE)


# ── 9. SUMMARY STATISTICS ───────────────────────────────────

# Helper: extract the last column from zonal() regardless of
# terra version (column layout changed in terra >= 1.7).
zonal_val <- function(raster, zones, ...) {
  z <- zonal(raster, zones, ...)
  z[[ncol(z)]]
}

# -- Province-wide --

primary_cv_vals <- values(A_mc_cv_primary, na.rm = TRUE)
sens_cv_vals    <- values(mc_cv, na.rm = TRUE)
div_vals        <- values(divergence, na.rm = TRUE)
n_valid         <- length(div_vals)

# Primary divergence (for comparison)
P80_prim <- P_exceed_80_primary
div_prim <- ifel(
  det_class == 1 & P80_prim >= 0.50, 1,
  ifel(det_class == 1 & P80_prim < 0.50, 2,
       ifel(det_class == 0 & P80_prim >= 0.50, 3, 4))
)
div_prim_vals <- values(div_prim, na.rm = TRUE)
n_valid_p     <- length(div_prim_vals)

province_stats <- tibble(
  Metric = c(
    "Province-wide mean (t/ha/yr)",
    "Province-wide SD (t/ha/yr)",
    "Province-wide CV (%)",
    "CV 5th percentile (%)",
    "CV 95th percentile (%)",
    "P(A > 80): province mean",
    "Confirmed severe (%)",
    "Uncertain hotspot (%)",
    "Hidden hotspot (%)",
    "Confirmed not severe (%)"
  ),
  R_only = c(
    round(global(A_mc_mean_primary, "mean", na.rm = TRUE)[[1]], 1),
    round(global(A_mc_sd_primary,   "mean", na.rm = TRUE)[[1]], 1),
    round(mean(primary_cv_vals), 1),
    round(quantile(primary_cv_vals, 0.05), 1),
    round(quantile(primary_cv_vals, 0.95), 1),
    round(global(P_exceed_80_primary, "mean", na.rm = TRUE)[[1]], 3),
    round(sum(div_prim_vals == 1) / n_valid_p * 100, 1),
    round(sum(div_prim_vals == 2) / n_valid_p * 100, 1),
    round(sum(div_prim_vals == 3) / n_valid_p * 100, 1),
    round(sum(div_prim_vals == 4) / n_valid_p * 100, 1)
  ),
  Multi_factor = c(
    round(global(mc_mean, "mean", na.rm = TRUE)[[1]], 1),
    round(global(mc_sd,   "mean", na.rm = TRUE)[[1]], 1),
    round(mean(sens_cv_vals), 1),
    round(quantile(sens_cv_vals, 0.05), 1),
    round(quantile(sens_cv_vals, 0.95), 1),
    round(global(P80, "mean", na.rm = TRUE)[[1]], 3),
    round(sum(div_vals == 1) / n_valid * 100, 1),
    round(sum(div_vals == 2) / n_valid * 100, 1),
    round(sum(div_vals == 3) / n_valid * 100, 1),
    round(sum(div_vals == 4) / n_valid * 100, 1)
  )
)

write.csv(province_stats,
          file.path(sens_dir, "tables", "sensitivity_province_summary.csv"),
          row.names = FALSE)

# -- Municipal-level --

mun_proj <- project(municipalities, mc_mean)

muni_stats <- tibble(
  Municipality     = mun_proj[[mun_name_field]],
  Det_Mean         = round(zonal_val(A_deterministic, mun_proj, fun = "mean", na.rm = TRUE), 1),
  Sens_MC_Mean     = round(zonal_val(mc_mean,         mun_proj, fun = "mean", na.rm = TRUE), 1),
  Sens_MC_SD       = round(zonal_val(mc_sd,           mun_proj, fun = "mean", na.rm = TRUE), 1),
  Sens_MC_CV       = round(zonal_val(mc_cv,           mun_proj, fun = "mean", na.rm = TRUE), 1),
  Primary_MC_CV    = round(zonal_val(A_mc_cv_primary, mun_proj, fun = "mean", na.rm = TRUE), 1),
  P_exceed_80_sens = round(zonal_val(P80,             mun_proj, fun = "mean", na.rm = TRUE), 3),
  P_exceed_80_prim = round(zonal_val(P_exceed_80_primary, mun_proj, fun = "mean", na.rm = TRUE), 3)
)

write.csv(muni_stats,
          file.path(sens_dir, "tables", "sensitivity_municipal_summary.csv"),
          row.names = FALSE)

cat("\n=== Province-wide summary ===\n")
print(province_stats, n = 20)
cat("\n=== Municipal summary ===\n")
print(muni_stats)


# ── 10. FIGURES ──────────────────────────────────────────────

cat("Generating figures...\n")

# Shared map theme
map_theme <- theme_bw() +
  theme(panel.background = element_rect(fill = "#d2eeff"))

map_extras <- function(loc = "br") {
  list(
    annotation_north_arrow(location = loc, which_north = "true",
                           style = north_arrow_minimal(),
                           pad_x = unit(0.1, "cm"),
                           pad_y = unit(0.8, "cm")),
    annotation_scale(location = loc, width_hint = 0.3, text_cex = 0.8)
  )
}

# -- Fig. 9a,b: CV comparison maps --

cv_pal <- brewer.pal(9, "YlOrRd")
cv_scale <- scale_fill_stepsn(
  colours = cv_pal, name = "CV (%)",
  breaks = seq(0, 40, 5), limits = c(0, 40),
  na.value = "transparent"
)

p_cv_ronly <- ggplot() +
  geom_spatraster(data = A_mc_cv_primary) +
  cv_scale +
  geom_spatvector(data = boundary, fill = NA, color = "grey30", linewidth = 0.4) +
  geom_spatvector(data = municipalities, fill = NA, color = "grey50", linewidth = 0.2) +
  map_theme +
  labs(x = "Longitude", y = "Latitude", title = "(a) R-only") +
  map_extras()

p_cv_multi <- ggplot() +
  geom_spatraster(data = mc_cv) +
  cv_scale +
  geom_spatvector(data = boundary, fill = NA, color = "grey30", linewidth = 0.4) +
  geom_spatvector(data = municipalities, fill = NA, color = "grey50", linewidth = 0.2) +
  map_theme +
  labs(x = "Longitude", y = "Latitude", title = "(b) Multi-factor") +
  map_extras()

fig_cv <- p_cv_ronly + p_cv_multi +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(sens_dir, "figures", "fig_cv_comparison.png"),
       fig_cv, width = 14, height = 7, dpi = 300)

# -- Fig. 9c: CV histogram --

cv_df <- bind_rows(
  tibble(CV = as.numeric(primary_cv_vals), Analysis = "R-only"),
  tibble(CV = as.numeric(sens_cv_vals),    Analysis = "Multi-factor")
)

fig_hist <- ggplot(cv_df, aes(x = CV, fill = Analysis)) +
  geom_histogram(alpha = 0.6, bins = 80, position = "identity") +
  scale_fill_manual(values = c("R-only" = "#2171b5", "Multi-factor" = "#d73027")) +
  labs(x = "Coefficient of variation (%)", y = "Number of pixels",
       title = "(c) Distribution of pixel-level CV") +
  theme_minimal(base_size = 12) +
  theme(legend.position = c(0.8, 0.8))

ggsave(file.path(sens_dir, "figures", "fig_cv_histogram.png"),
       fig_hist, width = 8, height = 5, dpi = 300)

# -- Fig. 9d: Divergence map --

levels(divergence) <- data.frame(
  id    = 1:4,
  label = c("Confirmed severe", "Uncertain hotspot",
            "Hidden hotspot", "Confirmed not severe")
)

fig_div <- ggplot() +
  geom_spatraster(data = divergence) +
  scale_fill_manual(
    values = c("Confirmed severe"     = "#d73027",
               "Uncertain hotspot"    = "#fc8d59",
               "Hidden hotspot"       = "#91bfdb",
               "Confirmed not severe" = "#e0e0e0"),
    name = "Classification", na.translate = FALSE
  ) +
  geom_spatvector(data = boundary, fill = NA, color = "grey30", linewidth = 0.4) +
  geom_spatvector(data = municipalities, fill = NA, color = "grey50", linewidth = 0.2) +
  map_theme +
  theme(legend.position = "inside",
        legend.position.inside = c(0.80, 0.22)) +
  labs(x = "Longitude", y = "Latitude",
       title = "(d) Divergence: deterministic vs. multi-factor (>80 t/ha/yr)") +
  map_extras()

ggsave(file.path(sens_dir, "figures", "fig_divergence_sensitivity.png"),
       fig_div, width = 8, height = 8, dpi = 300)


# ── 11. CONSOLE SUMMARY ─────────────────────────────────────

cat("\n")
cat("============================================================\n")
cat("  SENSITIVITY ANALYSIS COMPLETE\n")
cat("============================================================\n\n")
cat(sprintf("  R-only CV (province mean):        %.1f%%\n", mean(primary_cv_vals)))
cat(sprintf("  Multi-factor CV (province mean):   %.1f%%\n", mean(sens_cv_vals)))
cat(sprintf("  Multi-factor CV range (5-95%%):     %.1f%% - %.1f%%\n",
            quantile(sens_cv_vals, 0.05), quantile(sens_cv_vals, 0.95)))
cat(sprintf("  CV increase factor:                %.1fx\n",
            mean(sens_cv_vals) / mean(primary_cv_vals)))
cat(sprintf("  Rainfall share of total CV:        %.0f%%\n",
            mean(primary_cv_vals) / mean(sens_cv_vals) * 100))

n_det_severe <- sum(div_vals %in% c(1, 2))
n_uncertain  <- sum(div_vals == 2)
n_hidden     <- sum(div_vals == 3)

cat(sprintf("\n  Uncertain hotspots:  %d pixels (%.1f%% of det. Very Severe)\n",
            n_uncertain, ifelse(n_det_severe > 0, n_uncertain / n_det_severe * 100, 0)))
cat(sprintf("  Hidden hotspots:     %d pixels (%.1f%% of all valid pixels)\n",
            n_hidden, n_hidden / n_valid * 100))
cat(sprintf("  Total reclassified:  %.1f%% of all pixels\n",
            (n_uncertain + n_hidden) / n_valid * 100))
cat("\nOutputs saved to:", sens_dir, "\n")
cat("============================================================\n")
