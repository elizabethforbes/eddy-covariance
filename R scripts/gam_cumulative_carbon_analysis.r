# =============================================================================
# Cumulative Annual Carbon Flux Analysis
#
# Computes annual and running-cumulative NEE, GPP, and Reco for EF01 (organic)
# and EF02 (conventional) using observed EC tower data gap-filled with fitted
# GAM models. Also calculates Net Ecosystem Carbon Balance (NECB) by adding
# harvest carbon export to annual NEE.
#
# Models used for gap-filling:
#   bam_NEE              — bam() with AR(1), rho = 0.41
#   bam_GPP              — bam() with AR(1), rho = 0.61
#   gamm_Reco_smoothed   — gamm() with corARMA(p=1,q=1), 7-day smoothed climate
#
# These models must already be fitted and present in the R environment.
# Source the relevant scripts first:
#   source(here::here("R scripts", "gam_avgdaily_cropstage_NEE.r"))
#   source(here::here("R scripts", "gam_avgdaily_cropstage_GPP.r"))
#   source(here::here("R scripts", "gam_avgdaily_cropstage_Reco.r"))
#
# Also required: ec_daily_mngmnt3 (full dataset, all years, all crops)
#
# Analysis period: 2018–2021 (years with management + yield data)
# Units assumed: gC m⁻² d⁻¹ for all flux variables (NEE, GPP, Reco)
#   → verify with str(ec_daily_mngmnt) before running Section 6
#
# Kyle Hemes gap-filled dataset: available but not used here for methodological
# coherence — same model developed for mechanistic analysis used for gap-filling.
# (See email correspondence, March 2026.)
# =============================================================================

library(mgcv)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(zoo)      # rollmean for Reco lagged climate predictors

# =============================================================================
# SECTION 1: Data preparation
# =============================================================================

# Build full 2018–2021 dataset from the un-trimmed source.
# Unlike data1 (2018–2020, no alfalfa), data_cumul includes all crops and 2021,
# so that cumulative totals reflect actual field carbon exchange across the
# entire observation period.

data_cumul <- ec_daily_mngmnt3 %>%
  filter(year %in% 2018:2021) %>%
  mutate(
    management        = factor(as.character(management)),
    crop_stage_simple = relevel(factor(as.character(crop_stage_simple)), ref = "mature"),
    year_f            = factor(year)
  ) %>%
  arrange(management, year, doy) %>%
  group_by(management, year) %>%
  mutate(AR.start = row_number() == 1) %>%
  ungroup()


# Compute 7-day rolling means for Reco model predictors.
# Computed per management field (not reset at year boundaries), matching the
# training data behaviour in gam_avgdaily_cropstage_Reco.r.
# NA values appear for the first 6 rows of each management group after sorting.
data_cumul <- data_cumul %>%
  arrange(management, year, doy) %>%
  group_by(management) %>%
  mutate(
    air_temperature_7 = zoo::rollmean(air_temperature, 7, fill = NA, align = "right"),
    VPD_7             = zoo::rollmean(VPD,             7, fill = NA, align = "right")
  ) %>%
  ungroup()

cat("data_cumul dimensions:", nrow(data_cumul), "rows,", ncol(data_cumul), "columns\n")
cat("Years present:", sort(unique(data_cumul$year)), "\n")
cat("Management levels:", levels(data_cumul$management), "\n")
cat("Crop stage levels:", levels(data_cumul$crop_stage_simple), "\n")

# Quick check: flag rows where Reco gap-filling won't be possible
# (rolling mean NAs at start of each management group)
n_rollmean_na <- sum(is.na(data_cumul$air_temperature_7) | is.na(data_cumul$VPD_7))
cat("Rows with NA rolling means (unfillable for Reco):", n_rollmean_na, "\n")
# 882 total

# =============================================================================
# SECTION 1b: Bring in ERA5-filled climate predictors
# =============================================================================
#
# gapfill_VPDandAirT_forcumulative.r (run separately before this script) fits
# bias-correction regressions for ERA5 air temperature and VPD against tower
# observations and stores results in `df` with columns:
#   air_temperature_filled  — tower air_temperature; ERA5-corrected where NA
#   VPD_filled              — tower VPD (Pa); ERA5-corrected where NA
#
# We join those two columns onto data_cumul by management + date, then overwrite
# air_temperature and VPD so that downstream rolling-mean and gap-fill steps
# have continuous climate predictors.
#
# NOTE: df$VPD_filled is already in Pa (same units as data_cumul$VPD).
#       df is built with select() which drops year_f — do NOT use df directly
#       as a substitute for data_cumul.

if (!exists("df") || !all(c("VPD_filled", "air_temperature_filled") %in% names(df))) {
  stop(paste(
    "ERA5-filled climate columns not found.\n",
    "Please run gapfill_VPDandAirT_forcumulative.r first, then re-run this script."
  ))
}

climate_filled <- df %>%
  select(management, date, air_temperature_filled, VPD_filled)

data_cumul <- data_cumul %>%
  left_join(climate_filled, by = c("management", "date")) %>%
  mutate(
    airt_era5_filled = is.na(air_temperature) & !is.na(air_temperature_filled),
    VPD_era5_filled  = is.na(VPD)             & !is.na(VPD_filled),
    air_temperature  = coalesce(air_temperature, air_temperature_filled),
    VPD              = coalesce(VPD,             VPD_filled)
  ) %>%
  select(-air_temperature_filled, -VPD_filled)

era5_fill_summary <- data_cumul %>%
  group_by(management, year) %>%
  summarise(
    n_airt_era5_filled = sum(airt_era5_filled, na.rm = TRUE),
    n_VPD_era5_filled  = sum(VPD_era5_filled,  na.rm = TRUE),
    n_airt_still_na    = sum(is.na(air_temperature)),
    n_VPD_still_na     = sum(is.na(VPD)),
    .groups = "drop"
  )
cat("\nERA5 climate gap-fill coverage:\n")
print(era5_fill_summary)

# Re-compute 7-day rolling means now that climate predictors have fewer NAs.
# (Reco model uses air_temperature_7 and VPD_7; more filled days = more
#  gap-rows that can receive Reco predictions.)
data_cumul <- data_cumul %>%
  arrange(management, year, doy) %>%
  group_by(management) %>%
  mutate(
    air_temperature_7 = zoo::rollmean(air_temperature, 7, fill = NA, align = "right"),
    VPD_7             = zoo::rollmean(VPD,             7, fill = NA, align = "right")
  ) %>%
  ungroup()

n_rollmean_na_after <- sum(is.na(data_cumul$air_temperature_7) | is.na(data_cumul$VPD_7))
cat("Rolling-mean NAs after ERA5 fill:", n_rollmean_na_after,
    "(was", n_rollmean_na, "before)\n")
# Rolling-mean NAs after ERA5 fill: 12 (was 882 before)
# six at the start of each time series (conv and org)

# =============================================================================
# SECTION 2: Identify year random effect term names
# =============================================================================

# Helper: extract smooth term labels containing "year_f" from a gam/bam model.
# These labels are passed to the `exclude` argument of predict() when predicting
# for years not in the training data (2021), so that year-to-year variation is
# marginalised out and predictions reflect "average year" conditions.

get_re_terms <- function(model) {
  labels <- sapply(model$smooth, function(s) s$label)
  labels[grepl("year_f", labels)]
}

re_terms_nee  <- get_re_terms(bam_NEE)
re_terms_gpp  <- get_re_terms(bam_GPP)
re_terms_reco <- get_re_terms(gamm_Reco_smoothed$gam)

cat("\nYear RE terms — NEE:  ", paste(re_terms_nee,  collapse = ", "), "\n")
cat("Year RE terms — GPP:  ", paste(re_terms_gpp,  collapse = ", "), "\n")
cat("Year RE terms — Reco: ", paste(re_terms_reco, collapse = ", "), "\n")

# =============================================================================
# SECTION 3: Gap-filling function
# =============================================================================

#' Gap-fill a flux column using GAM/BAM model predictions
#'
#' Strategy:
#'   - Observed (non-NA) rows: returned as-is; SE = 0.
#'   - Gap rows in training years (2018–2020): predict with year random effects
#'     included. year_f factor set to match model's fitted levels.
#'   - Gap rows in extrapolation year (2021): predict with year RE excluded
#'     (average-year prediction). year_f set to "2018" as a valid placeholder
#'     (value does not affect predictions when the term is excluded).
#'
#' @param data         data_cumul or a subset thereof
#' @param model        fitted bam or gam object (for gamm, pass model$gam)
#' @param flux_col     character; column name of observed flux (e.g. "NEE")
#' @param re_terms     character vector; year RE term labels from get_re_terms()
#' @param is_gamm      logical; TRUE if model is the $gam component of a gamm()
#' @param training_yrs integer vector; years the model was fitted on
#'
#' @return data.frame with columns:
#'   flux_filled   — observed where available, GAM-predicted in gaps
#'   flux_se       — 0 for observed; predict() SE for gap-filled rows
#'   is_gapfilled  — logical; TRUE where a prediction was used

gapfill_flux <- function(data, model, flux_col, re_terms,
                         is_gamm = FALSE, training_yrs = 2018:2020) {

  flux_obs     <- data[[flux_col]]
  flux_filled  <- flux_obs
  flux_se      <- rep(0.0, nrow(data))
  is_gapfilled <- is.na(flux_obs)

  if (!any(is_gapfilled)) {
    return(data.frame(flux_filled = flux_filled, flux_se = flux_se,
                      is_gapfilled = is_gapfilled))
  }

  gap_idx          <- which(is_gapfilled)
  training_gap_idx <- gap_idx[data$year[gap_idx] %in% training_yrs]
  extrap_gap_idx   <- gap_idx[!data$year[gap_idx] %in% training_yrs]

  # Year factor levels: derived directly from training_yrs rather than
  # model$model$year_f, because bam() does not reliably store the full model
  # frame and model$model$year_f may be NULL.
  model_year_levels <- as.character(sort(training_yrs))
  fallback_year_f   <- model_year_levels[1]   # placeholder for extrap rows

  # Crop stage levels: model$xlevels is reliably stored by both bam() and gam()
  # for factors appearing in the formula (parametric or smooth terms).
  # Falls back to unique levels in the training-period data if xlevels is absent.
  model_stage_levels <- if (!is.null(model$xlevels$crop_stage_simple)) {
    model$xlevels$crop_stage_simple
  } else {
    as.character(unique(data$crop_stage_simple[data$year %in% training_yrs &
                                                !is.na(data$crop_stage_simple)]))
  }

  predict_fn <- if (is_gamm) {
    function(nd, excl) predict.gam(model, newdata = nd, se.fit = TRUE, exclude = excl)
  } else {
    function(nd, excl) predict.bam(model, newdata = nd, se.fit = TRUE, exclude = excl)
  }

  # ---- Training-period gaps (include year RE) --------------------------------
  if (length(training_gap_idx) > 0) {
    nd <- data.frame(data[training_gap_idx, ])   # coerce to plain data.frame
    # Safeguard: if year_f is absent (e.g. data_cumul rebuilt without it in the
    # current session), create it from year before coercing to factor.
    if (!"year_f" %in% names(nd)) nd$year_f <- as.character(nd$year)
    nd$year_f <- factor(as.character(nd$year_f), levels = model_year_levels)
    nd$crop_stage_simple <- factor(as.character(nd$crop_stage_simple),
                                   levels = model_stage_levels)
    pred <- predict_fn(nd, excl = NULL)
    flux_filled[training_gap_idx] <- pred$fit
    flux_se[training_gap_idx]     <- pred$se.fit
  }

  # ---- Extrapolation gaps (exclude year RE, average-year prediction) ---------
  if (length(extrap_gap_idx) > 0) {
    nd <- data.frame(data[extrap_gap_idx, ])     # coerce to plain data.frame
    nd$year_f <- factor(fallback_year_f, levels = model_year_levels)
    nd$crop_stage_simple <- factor(as.character(nd$crop_stage_simple),
                                   levels = model_stage_levels)
    pred <- predict_fn(nd, excl = re_terms)
    flux_filled[extrap_gap_idx] <- pred$fit
    flux_se[extrap_gap_idx]     <- pred$se.fit
  }

  data.frame(flux_filled = flux_filled, flux_se = flux_se,
             is_gapfilled = is_gapfilled)
}

# =============================================================================
# SECTION 4: Apply gap-filling
# =============================================================================

cat("\nGap-filling NEE...\n")
nee_gf <- gapfill_flux(data_cumul, bam_NEE, "NEE", re_terms_nee)
data_cumul$NEE_filled    <- nee_gf$flux_filled
data_cumul$NEE_se        <- nee_gf$flux_se
data_cumul$NEE_gapfilled <- nee_gf$is_gapfilled

cat("Gap-filling GPP...\n")
gpp_gf <- gapfill_flux(data_cumul, bam_GPP, "GPP", re_terms_gpp)
data_cumul$GPP_filled    <- gpp_gf$flux_filled
data_cumul$GPP_se        <- gpp_gf$flux_se
data_cumul$GPP_gapfilled <- gpp_gf$is_gapfilled

cat("Gap-filling Reco...\n")
# gamm_Reco_smoothed uses air_temperature_7 and VPD_7.
# Rows where rolling means are NA cannot be gap-filled; flag them separately.
reco_gf <- gapfill_flux(data_cumul, gamm_Reco_smoothed$gam,
                        "Reco", re_terms_reco, is_gamm = TRUE)
data_cumul$Reco_filled     <- reco_gf$flux_filled
data_cumul$Reco_se         <- reco_gf$flux_se
data_cumul$Reco_gapfilled  <- reco_gf$is_gapfilled
# Rows that were already NA in Reco AND have NA rolling means remain NA after
# gap-filling (gapfill_flux skips rows with NA predictors inside predict()).
# Flag these explicitly for reporting:
data_cumul$Reco_unfillable <- data_cumul$Reco_gapfilled &
  (is.na(data_cumul$air_temperature_7) | is.na(data_cumul$VPD_7))

# Gap-fill summary (print before proceeding)
gap_summary <- data_cumul %>%
  group_by(management, year) %>%
  summarise(
    n_days         = n(),
    n_gap_NEE      = sum(is.na(NEE)),
    n_gap_GPP      = sum(is.na(GPP)),
    n_gap_Reco     = sum(is.na(Reco)),
    pct_gap_NEE    = round(100 * n_gap_NEE  / n_days, 1),
    pct_gap_GPP    = round(100 * n_gap_GPP  / n_days, 1),
    pct_gap_Reco   = round(100 * n_gap_Reco / n_days, 1),
    n_unfillable_Reco = sum(Reco_unfillable),
    .groups        = "drop"
  )

cat("\n--- Gap-fill summary ---\n")
print(gap_summary)

# =============================================================================
# SECTION 5: Harvest carbon export
# =============================================================================

# Carbon fraction and test weight conversion factors.
# In EC sign convention: carbon leaving the ecosystem is a POSITIVE flux (loss),
# so NECB = NEE_annual + C_harvested. Fields with NECB < 0 are net carbon sinks.
#
# Conversion formula:
#   C_export (gC m⁻²) = yield × kg_per_unit × moisture_corr × c_fraction
#                       / 4046.86 (m² per acre) × 1000 (kg → g)
#
# Sources:
#   Test weights, standard moisture:
#     USDA Federal Grain Inspection Service (2013). Grain Inspection Handbook,
#     Book II. USDA GIPSA, Washington, DC.
#   Carbon fractions (cereal grains, legumes):
#     IPCC (2006). 2006 IPCC Guidelines for National GHG Inventories, Volume 4:
#     Agriculture, Forestry and Other Land Use, Chapter 11, Table 11.2. IGES, Japan.
#   Carbon fraction — soybean (elevated vs. cereals due to ~20% lipid content):
#     Watanabe I., 1976. Transformation Factor from CO2 Net Assimilation to Dry Matter in Crop Plants
#     AGricultural and Food Sciences, Jarq-japan Agricultural Research Quarterly
#   Carbon fraction — alfalfa hay (DM basis):
#     Bolinder, M.A., et al. (2007). An approach for estimating net primary
#     productivity and annual carbon inputs to soil for common agricultural
#     crops in Canada. Agriculture, Ecosystems & Environment 118(1-4), 29–42.

c_conv <- tribble(
  ~crop,         ~yield_units, ~kg_per_unit, ~moisture_corr, ~c_fraction, ~note,
  "corn",        "bu/ac",      25.401,       0.845,          0.44,        
  "56 lb/bu (page 1-27, table 1.9: USDA FGIS) at 15.5% std moisture (USDA FGIS 2013); C: IPCC (2006) Vol.4 Ch.11 Table 11.2",
  # https://grains.org/corn_report/corn-harvest-quality-report-2018-2019/4/ confirms ~15-16% moisture content nationally
  # https://corn.ces.ncsu.edu/news/harvesting-corn-what-grain-moisture-should-i-harvest-corn-at/ confirms that <= 15% is ideal
  # https://www.smallfarmcanada.ca/resources/standard-weights-per-bushel-for-agricultural-commodities confirms standard weights in kg to lb
  # https://courses.ecampus.oregonstate.edu/ans312/two/corn_trans.htm for general percentages of macronutrients that guide conversion to C content
  # and https://www.feedtables.com/content/maize assuming that kCal/kg divided by 100 is approx. the carbon content of corn grain
  
  "soybeans",    "bu/ac",      27.216,       0.870,          0.51,        
  "60 lb/bu (page 1-27, table 1.9: USDA FGIS) at 13.0% std moisture (USDA FGIS 2013); C elevated due to lipid content (Monfreda et al. 2008)",
  # https://www.smallfarmcanada.ca/resources/standard-weights-per-bushel-for-agricultural-commodities standard conversion, lb to kg
  # https://cropwatch.unl.edu/2018/enhancing-soybean-storage-starts-harvest-moisture/ confirms 13% market moisture content for storage
  # Watanabe, 1976: soybeans themselves have approx. 51% carbon content at ripeness
  
  "barley",      "bu/ac",      21.772,       0.860,          0.42,        
  "48 lb/bu (page 1-27, table 1.9: USDA FGIS) at 14.0% std moisture (USDA FGIS 2013); C: IPCC (2006) Vol.4 Ch.11 Table 11.2",
  # https://www.smallfarmcanada.ca/resources/standard-weights-per-bushel-for-agricultural-commodities standard conversion, lb to kg
  # https://extension.umn.edu/small-grains-harvest-and-storage/drying-wheat-and-barley confirms 13-14% moisture content for barley
  # Rogers et al., 2025 indicates 42.2%C (from 422g/kg grain) (table 2)
  
  "alfalfa hay", "ton/ac",     907.185,      1.000,          0.45,        
  "Already on DM basis per management notes; C: Bolinder et al. (2007)"
  # assumes 0.45% on average of all aboveground plant parts
  
) %>%
  mutate(
    # gC exported per unit yield per m²
    gC_per_unit_per_m2 = kg_per_unit * moisture_corr * c_fraction / 4046.86 * 1000,
    kgC_per_unit_per_m2 = kg_per_unit *moisture_corr *c_fraction / 4046.86
    # reminder: 4046.86 (m² per acre) × 1000 (kg → g)
  )

# Harvest events: one row per harvest date
# Dates marked "estimated" in management records are flagged in note column
harvest_events <- tribble(
  ~management,    ~date,         ~year, ~crop,         ~yield,  ~units,    ~date_estimated,
  "organic",      "2018-10-26",  2018,  "corn",         91.443, "bu/ac",   FALSE,
  "organic",      "2019-11-06",  2019,  "soybeans",     26.712, "bu/ac",   FALSE,
  "organic",      "2020-07-09",  2020,  "barley",       54.407, "bu/ac",   FALSE,
  "organic",      "2021-06-05",  2021,  "alfalfa hay",   1.000, "ton/ac",  FALSE,
  "organic",      "2021-08-10",  2021,  "alfalfa hay",   1.000, "ton/ac",  FALSE,
  "conventional", "2018-10-04",  2018,  "corn",        194.000, "bu/ac",   FALSE,
  "conventional", "2019-10-01",  2019,  "soybeans",     56.000, "bu/ac",   TRUE,
  "conventional", "2020-07-01",  2020,  "barley",       43.000, "bu/ac",   TRUE,
  "conventional", "2021-10-01",  2021,  "soybeans",     77.000, "bu/ac",   TRUE
) %>%
  mutate(
    date = as.Date(date),
    doy  = as.integer(format(date, "%j"))
  ) %>%
  left_join(c_conv %>% select(crop, gC_per_unit_per_m2, kgC_per_unit_per_m2), by = "crop") %>% 
  mutate(C_exported_gC_m2 = yield * gC_per_unit_per_m2,
         C_exported_kgC_m2 = yield * kgC_per_unit_per_m2)

# Annual harvest C export (sum of cuts within each field-year)
harvest_by_year <- harvest_events %>%
  group_by(management, year) %>%
  summarise(
    C_harvested_gC_m2 = sum(C_exported_gC_m2),
    C_harvested_kgC_m2 = sum(C_exported_kgC_m2),
    harvest_dates     = paste(format(date, "%b %d"), collapse = " + "),
    crops             = paste(unique(crop), collapse = " + "),
    any_date_estimated = any(date_estimated),
    .groups           = "drop"
  )

cat("\n--- Harvest C export (gC m⁻²) ---\n")
print(harvest_by_year %>% select(management, year, crops, C_harvested_gC_m2,
                                 C_harvested_kgC_m2, harvest_dates, any_date_estimated))

# =============================================================================
# SECTION 6: Unit conversion — µmol CO₂ m⁻² s⁻¹  →  kgC m⁻² d⁻¹ →  gC m⁻² d⁻¹
# =============================================================================
# EC tower fluxes are daily averages in µmol CO₂ m⁻² s⁻¹.
# To convert to gC m⁻² d⁻¹ (the unit needed for annual accumulation):
#
#   gC m⁻² d⁻¹ = µmol CO₂ m⁻² s⁻¹
#                 × 86,400  s d⁻¹
#                × 12.011   g C mol⁻¹
#                × 1×10⁻⁶   mol µmol⁻¹
#               = µmol m⁻² s⁻¹ × 1.037750

# To convert to kgC m⁻² d⁻¹ 
#   kgC m⁻² d⁻¹ = µmol CO₂ m⁻² s⁻¹
#                 × 86,400  s d⁻¹
#                × 12.011   g C mol⁻¹
#                × 1×10⁻⁶   mol µmol⁻¹
#               = µmol m⁻² s⁻¹ × 1.037750
#               x 0.001 (g to kg) = 0.00103775
#
# Applied to both observed columns (NEE, GPP, Reco) AND gap-filled columns
# (NEE_filled, GPP_filled, Reco_filled). SEs scale by the same factor.
# Harvest C export table (Section 5) is already in gC m⁻² and unaffected.

UMOL_TO_GC_PER_DAY <- 86400 * 12.011 * 1e-6   # = 1.037750
UMOL_TO_KGC_PER_DAY <- 86400 * 12.011 * 1e-6 * 0.001 # = 0.00103775

cat("\n--- Raw flux range (µmol CO₂ m⁻² s⁻¹, pre-conversion) ---\n")
cat("NEE  range:", range(data_cumul$NEE,  na.rm = TRUE), "\n")
cat("GPP  range:", range(data_cumul$GPP,  na.rm = TRUE), "\n")
cat("Reco range:", range(data_cumul$Reco, na.rm = TRUE), "\n")

data_cumul <- data_cumul %>%
  mutate(
    # Convert observed columns (for reference / diagnostics)
    NEE_g  = NEE  * UMOL_TO_GC_PER_DAY,
    GPP_g  = GPP  * UMOL_TO_GC_PER_DAY,
    Reco_g = Reco * UMOL_TO_GC_PER_DAY,
    NEE_kg = NEE  * UMOL_TO_KGC_PER_DAY,
    GPP_kg = GPP  * UMOL_TO_KGC_PER_DAY,
    Reco_kg = Reco* UMOL_TO_KGC_PER_DAY,
    # Convert gap-filled columns and their SEs
    NEE_filled_g  = NEE_filled  * UMOL_TO_GC_PER_DAY,
    GPP_filled_g  = GPP_filled  * UMOL_TO_GC_PER_DAY,
    Reco_filled_g = Reco_filled * UMOL_TO_GC_PER_DAY,
    NEE_g_se      = NEE_se      * UMOL_TO_GC_PER_DAY,
    GPP_g_se      = GPP_se      * UMOL_TO_GC_PER_DAY,
    Reco_g_se     = Reco_se     * UMOL_TO_GC_PER_DAY,
    NEE_filled_kg  = NEE_filled  * UMOL_TO_KGC_PER_DAY,
    GPP_filled_kg  = GPP_filled  * UMOL_TO_KGC_PER_DAY,
    Reco_filled_kg = Reco_filled * UMOL_TO_KGC_PER_DAY,
    NEE_kg_se      = NEE_se      * UMOL_TO_KGC_PER_DAY,
    GPP_kg_se      = GPP_se      * UMOL_TO_KGC_PER_DAY,
    Reco_kg_se     = Reco_se     * UMOL_TO_KGC_PER_DAY
  )

cat("\n--- Converted flux range (gC m⁻² d⁻¹, post-conversion) ---\n")
cat("NEE  range:", range(data_cumul$NEE_g,  na.rm = TRUE), "\n")
cat("GPP  range:", range(data_cumul$GPP_g,  na.rm = TRUE), "\n")
cat("Reco range:", range(data_cumul$Reco_g, na.rm = TRUE), "\n")

# =============================================================================
# SECTION 6b: Partial-year flag
# =============================================================================
# Both 2018 and 2021 are partial-year observations:
#   2018: tower data starts late in the year (does not begin on Jan 1)
#   2021: management records end at the final harvest date for each field
#           Organic (EF01):      data ends 2021-08-10 (second alfalfa cut)
#           Conventional (EF02): data ends 2021-10-01 (soybean harvest)
#
# Annual totals for these years should NOT be directly compared to full-year
# totals without a clear caveat. However, because 2021 data ends AT the final
# harvest, NECB for 2021 does capture the complete growing-season carbon
# balance including final C export — which is the ecologically meaningful
# quantity.
#
# partial_year is flagged data-driven (first_doy > 14 OR last_doy < 351)
# rather than hardcoded by year. The n_days column documents coverage.

partial_year_flag <- data_cumul %>%
  group_by(management, year) %>%
  summarise(
    first_doy = min(doy),
    last_doy  = max(doy),
    n_days    = n(),
    # Flag any year where data does not cover (roughly) the full calendar year.
    # Threshold: starts after Jan 14 (doy > 14) OR ends before Dec 17 (doy < 351).
    # This catches both 2018 (late start) and 2021 (early end at final harvest).
    partial_year = (min(doy) > 14) | (max(doy) < 351),
    .groups   = "drop"
  )

cat("\n--- Observation coverage by field-year ---\n")
print(partial_year_flag)

# =============================================================================
# SECTION 7: Cumulative sums
# =============================================================================

# Sum daily filled flux values within each management × year.
# Uncertainty is propagated as: SE_cumulative(t) = sqrt( Σ SE_i² for i = 1..t )
# This treats gap-filled days as having independent prediction errors and
# observed days as error-free. As a result, uncertainty bounds reflect
# gap-filling uncertainty only (a conservative lower bound on total uncertainty).

# update, April 2026: use only the 2018 data from Aug. 13, 2018 which is the earliest record we have for org.
# update 2, April 2026: exclude the 2021 data from main analysis, potentially include in supplemental:

data_cumul_2 <- data_cumul %>% 
  filter(date >= "2018-08-13") # removed 105 observations from conventional field

cumul_daily <- data_cumul_2 %>%
  filter(year != 2021) %>% 
  arrange(management, year, doy) %>%
  group_by(management, year) %>%
  mutate(
    # Running cumulative sums
    cumNEE_g   = cumsum(replace_na(NEE_filled_g,  0)),
    cumGPP_g   = cumsum(replace_na(GPP_filled_g,  0)),
    cumReco_g  = cumsum(replace_na(Reco_filled_g, 0)),
    cumNEE_kg  = cumsum(replace_na(NEE_filled_kg,  0)),
    cumGPP_kg  = cumsum(replace_na(GPP_filled_kg,  0)),
    cumReco_kg = cumsum(replace_na(Reco_filled_kg, 0)),
    # Propagated uncertainty (SE of cumulative sum).
    # replace_na(..., 0) prevents a single NA predictor in one gap row from
    # propagating through cumsum() and wiping out all subsequent SE values.
    # Rows where predict() returned NA (unfillable gaps) contribute 0 variance.
    cumNEE_g_se  = sqrt(cumsum(replace_na(NEE_g_se^2,  0))),
    cumGPP_g_se  = sqrt(cumsum(replace_na(GPP_g_se^2,  0))),
    cumReco_g_se = sqrt(cumsum(replace_na(Reco_g_se^2, 0))),
    cumNEE_kg_se  = sqrt(cumsum(replace_na(NEE_kg_se^2,  0))),
    cumGPP_kg_se  = sqrt(cumsum(replace_na(GPP_kg_se^2,  0))),
    cumReco_kg_se = sqrt(cumsum(replace_na(Reco_kg_se^2, 0))),
    # Running gap-fill count
    n_days         = row_number(),
    cum_gapfill_NEE  = cumsum(as.integer(NEE_gapfilled)),
    cum_gapfill_GPP  = cumsum(as.integer(GPP_gapfilled)),
    cum_gapfill_Reco = cumsum(as.integer(Reco_gapfilled))
  ) %>%
  ungroup()

# Annual totals: take the last row of each management-year
annual_totals <- cumul_daily %>%
  group_by(management, year) %>%
  summarise(
    n_days          = n(),
    first_doy       = min(doy),
    last_doy        = max(doy),
    partial_year    = (min(doy) > 14) | (max(doy) < 351),  # late start OR early end
    n_gap_NEE       = sum(NEE_gapfilled),
    n_gap_GPP       = sum(GPP_gapfilled),
    n_gap_Reco      = sum(Reco_gapfilled),
    pct_gap_NEE     = round(100 * n_gap_NEE  / n_days, 1),
    # in g:
    NEE_g_annual      = sum(replace_na(NEE_filled_g,  0)),
    GPP_g_annual      = sum(replace_na(GPP_filled_g,  0)),
    Reco_g_annual     = sum(replace_na(Reco_filled_g, 0)),
    NEE_annual_g_se   = sqrt(sum(replace_na(NEE_g_se^2, 0))), # there are 15 NAs in this row
    GPP_annual_g_se   = sqrt(sum(GPP_g_se^2)),
    Reco_annual_g_se  = sqrt(sum(Reco_g_se^2)),
    # in kg:
    NEE_kg_annual      = sum(replace_na(NEE_filled_kg,  0)),
    GPP_kg_annual      = sum(replace_na(GPP_filled_kg,  0)),
    Reco_kg_annual     = sum(replace_na(Reco_filled_kg, 0)),
    NEE_annual_kg_se   = sqrt(sum(NEE_kg_se^2)),
    GPP_annual_kg_se   = sqrt(sum(GPP_kg_se^2)),
    Reco_annual_kg_se  = sqrt(sum(Reco_kg_se^2)),
    .groups         = "drop"
  )

# =============================================================================
# SECTION 7b: Mass-balance diagnostics
# =============================================================================
# NEE, GPP, and Reco are gap-filled with INDEPENDENT models, so there is no
# algebraic guarantee that NEE_filled = Reco_filled - GPP_filled.  Before
# trusting the annual totals, two checks are run:
#
#   Check 1 — Sign convention on OBSERVED days only:
#     Computes NEE_implied = Reco_obs - GPP_obs and NEE_implied_alt = GPP_obs - Reco_obs
#     and regresses each against observed NEE. A slope near +1 with low RMSE
#     confirms the sign convention (positive NEE = source in standard EC notation).
#
#   Check 2 — Annual mass-balance closure:
#     Compares NEE_annual against (Reco_annual - GPP_annual).
#     Large discrepancies (> ~50 gC m⁻²) indicate independent model drift and
#     should be noted in the methods; NEE_filled is the primary EC measurement
#     and should be treated as authoritative over partitioned GPP/Reco.

# --- Check 1: sign convention on fully-observed days -------------------------
obs_days <- data_cumul %>%
  filter(!is.na(NEE_g) & !is.na(GPP_g) & !is.na(Reco_g)) %>%
  mutate(
    NEE_g_implied     = Reco_g - GPP_g,   # standard: NEE = Reco - GPP
    NEE_g_implied_alt = GPP_g  - Reco_g   # alternative sign convention
  )

lm_std <- lm(NEE ~ NEE_g_implied,     data = obs_days)
lm_alt <- lm(NEE ~ NEE_g_implied_alt, data = obs_days)

rmse <- function(m) round(sqrt(mean(resid(m)^2)), 4)

cat("\n=== Mass-balance check 1: sign convention (observed days only) ===\n")
cat(sprintf("  NEE ~ (Reco - GPP):  slope = %.3f, R2 = %.3f, RMSE = %.4f\n",
            coef(lm_std)[2], summary(lm_std)$r.squared, rmse(lm_std)))
cat(sprintf("  NEE ~ (GPP - Reco):  slope = %.3f, R2 = %.3f, RMSE = %.4f\n",
            coef(lm_alt)[2], summary(lm_alt)$r.squared, rmse(lm_alt)))
cat("  >> Convention closest to slope=1 / lowest RMSE is the correct one.\n")

# --- Check 2: annual mass-balance closure ------------------------------------
cat("\n=== Mass-balance check 2: annual closure (gap-filled totals) ===\n")
cat("  Expected: NEE_annual ≈ Reco_annual - GPP_annual  (or GPP - Reco; see Check 1)\n\n")

mb_check <- annual_totals %>%
  mutate(
    implied_NEE_std = Reco_g_annual - GPP_g_annual,
    implied_NEE_alt = GPP_g_annual  - Reco_g_annual,
    closure_std     = round(NEE_g_annual - implied_NEE_std, 1),
    closure_alt     = round(NEE_g_annual - implied_NEE_alt, 1)
  ) %>%
  select(management, year, n_days, pct_gap_NEE,
         NEE_g_annual, GPP_g_annual, Reco_g_annual,
         implied_NEE_std, closure_std,
         implied_NEE_alt, closure_alt)

print(as.data.frame(mb_check))
cat("\n  closure_std = NEE_annual - (Reco - GPP): near 0 = good closure\n")
cat("  closure_alt = NEE_annual - (GPP - Reco): near 0 = good closure\n")
cat("  Large values indicate independent model drift over gap-filled periods.\n")

# --- Check 3: daily residuals split by gap-filled vs observed ----------------
cat("\n=== Mass-balance check 3: daily closure by data type ===\n")
data_cumul %>%
  mutate(
    NEE_implied = Reco_filled - GPP_filled,
    daily_closure = NEE_filled - NEE_implied
  ) %>%
  group_by(management, year, gap_type = if_else(NEE_gapfilled, "gap_filled", "observed")) %>%
  summarise(
    n            = n(),
    mean_closure = round(mean(daily_closure, na.rm = TRUE), 3),
    sd_closure   = round(sd(daily_closure,   na.rm = TRUE), 3),
    .groups      = "drop"
  ) %>%
  print()
cat("\n  mean_closure near 0 = models agree on average\n")
cat("  Systematic non-zero mean_closure on gap_filled rows = model drift\n")

# --- Check 4: daily flux inspection plot for gap-filled days -----------------
# Plots three model-predicted daily values on the same axis for gap-filled rows:
#   NEE_filled        (from bam_NEE)
#   Reco_filled       (from gamm_Reco_smoothed)
#   GPP_filled        (from bam_GPP, shown as negative to match NEE sign space)
#   Reco - GPP        (implied NEE from the two component models)
#
# Divergence between NEE_filled and (Reco_filled - GPP_filled) on gap days
# indicates which model is driving the discrepancy and whether it looks
# physically plausible for the crop/season.

gap_daily_inspect <- data_cumul %>%
  filter(NEE_gapfilled) %>%
  select(management, year, doy, crop_stage_simple,
         NEE_filled, GPP_filled, Reco_filled) %>%
  mutate(implied_NEE = Reco_filled - GPP_filled,
         closure     = NEE_filled - implied_NEE) %>%
  tidyr::pivot_longer(
    cols      = c(NEE_filled, Reco_filled, GPP_filled, implied_NEE),
    names_to  = "flux",
    values_to = "value"
  ) %>%
  mutate(flux = factor(flux,
                       levels = c("NEE_filled", "implied_NEE",
                                  "Reco_filled", "GPP_filled"),
                       labels = c("NEE (bam_NEE)",
                                  "Reco \u2212 GPP (implied)",
                                  "Reco (gamm)",
                                  "GPP (bam_GPP)")))

p_gapfill_inspect <- ggplot(gap_daily_inspect,
                            aes(x = doy, y = value, colour = flux)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  facet_grid(year ~ management, scales = "free_y") +
  scale_colour_manual(
    values = c("NEE (bam_NEE)"       = "black",
               "Reco \u2212 GPP (implied)" = "firebrick",
               "Reco (gamm)"         = "steelblue",
               "GPP (bam_GPP)"       = "darkgreen")
  ) +
  labs(
    title   = "Gap-filled days: daily model predictions",
    subtitle = "Black = NEE model | Red = Reco\u2212GPP implied NEE | Blue = Reco | Green = GPP\nLarge black\u2013red divergence indicates model drift on gap days",
    x       = "Day of year",
    y       = expression("Flux (gC m"^{-2}*" d"^{-1}*")"),
    colour  = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

print(p_gapfill_inspect)

# =============================================================================
# SECTION 8: NECB
# =============================================================================

# NECB = NEE_annual + C_harvested
# Sign convention maintained from EC measurements:
#   Negative NECB = net carbon sink (ecosystem accumulating carbon)
#   Positive NECB = net carbon source (ecosystem losing carbon overall)
#
# Note: NECB here accounts for grain/forage export only (above-ground removal).
# Root and residue C inputs to soil partially offset harvest losses but are not
# estimated here due to data limitations; reported values represent a minimum
# (most positive) NECB estimate. See Chapin et al. (2006) for NECB framework.
#
# Reference: Chapin, F.S. III, et al. (2006). Reconciling carbon-cycle concepts,
# terminology, and methods. Ecosystems 9, 1041–1050.

annual_totals <- annual_totals %>%
  left_join(harvest_by_year %>%
              select(management, year, C_harvested_gC_m2, crops,
                     harvest_dates, any_date_estimated),
            by = c("management", "year")) %>%
  mutate(
    # Add C export; propagate NEE SE only (yield treated as fixed/measured)
    NECB    = NEE_g_annual + C_harvested_gC_m2,
    NECB_se = NEE_annual_g_se
  )

cat("\n--- Annual flux summary (gC m⁻² [partial yr for 2021]) ---\n")
cat("Note: 2021 organic ends doy",
    partial_year_flag$last_doy[partial_year_flag$management == "organic"   & partial_year_flag$year == 2021],
    "(Aug 10); conventional ends doy",
    partial_year_flag$last_doy[partial_year_flag$management == "conventional" & partial_year_flag$year == 2021],
    "(Oct 1)\n\n")
print(annual_totals %>%
        select(management, year, n_days, partial_year, NEE_g_annual, GPP_g_annual,
               Reco_g_annual, C_harvested_gC_m2, NECB, pct_gap_NEE, any_date_estimated))

# =============================================================================
# SECTION 9: Visualizations
# =============================================================================

# Colour palette: consistent with main analysis
mgmt_colors <- c("organic" = "#E69F00", "conventional" = "tomato")
mgmt_labels <- c("organic" = "Organic (EF01)", "conventional" = "Conventional (EF02)")

# --- 9a: Running cumulative NEE, by year -------------------------------------
# Harvest event markers (per management): vertical dotted lines
# update, April 2026: remove 2021 data from main analysis (consider including in supplement)
harvest_doy_by_year <- harvest_events %>%
  filter(year != 2021) %>% 
  select(management, year, doy, crop) %>%
  mutate(management = as.character(management))

# Diagnostic: print max cumulative SE to assess ribbon visibility.
# If max SE is small relative to plot range, ribbon will be near-invisible
# at any alpha. In that case, omit ribbon and note in caption instead.
cat("\n--- Cumulative NEE SE range (max per field-year) ---\n")
print(
  cumul_daily %>%
    group_by(management, year) %>%
    summarise(max_cumNEE_se    = round(max(cumNEE_g_se), 2),
              ci_95_halfwidth  = round(1.96 * max(cumNEE_g_se), 2),
              .groups          = "drop")
)
# If ci_95_halfwidth < ~15 gC m⁻² the ribbon will be very narrow on this scale.
# In that case: remove geom_ribbon and add to caption:
#   "Gap-filling uncertainty (95% CI from propagated prediction SE) was
#    <X gC m⁻² at all time points and is not shown."

p_cumNEE <- ggplot(cumul_daily,
                   aes(x = doy, y = cumNEE_g, colour = management, fill = management)) +
  geom_ribbon(aes(ymin = cumNEE_g - 1.96 * cumNEE_g_se,
                  ymax = cumNEE_g + 1.96 * cumNEE_g_se),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.45) +
  geom_vline(
    data = harvest_doy_by_year,
    aes(xintercept = doy, colour = management),
    linetype = "dotted", linewidth = 0.65, alpha = 0.75,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ year, 
             # nrow = 2, 
             scales = "fixed") +
  scale_colour_manual(values = mgmt_colors, labels = mgmt_labels) +
  scale_fill_manual(values = mgmt_colors, labels = mgmt_labels) +
  labs(
    # title    = "Running cumulative NEE by field and year",
    # subtitle = "Shaded band = 95% CI (gap-filling uncertainty) | Dotted lines = harvest dates",
    x        = "Day of year",
    y        = expression("Cumulative NEE (gC m"^{-2}*")"),
    colour   = NULL, fill = NULL
  ) +
  ylim(-410, 610) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

p_cumNEE
# --- 9b: Running cumulative NEE vs NECB (step-function harvest jumps) --------
# Build a step-function version of NECB for each field-year by adding harvest
# C export at the harvest doy. Before the first harvest: NECB = cumNEE.
# After each harvest: NECB = cumNEE + cumulative C exported to that point.

cumul_with_harvest <- cumul_daily %>%
  select(management, year, doy, cumNEE_g, cumNEE_g_se) %>%
  left_join(harvest_events %>%
              select(management, year, doy_harvest = doy, C_exported_gC_m2),
            by = c("management", "year"),
            relationship = "many-to-many") %>%
  group_by(management, year, doy) %>%
  summarise(
    cumNEE     = first(cumNEE_g),
    cumNEE_se  = first(cumNEE_g_se),
    cum_C_exp  = sum(C_exported_gC_m2[doy_harvest <= doy], na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(cumNECB = cumNEE + cum_C_exp)

p_cumNECB <- ggplot(cumul_with_harvest,
                    aes(x = doy, colour = management, fill = management)) +
  # 95% CI ribbon on NECB (harvest C treated as fixed, so NECB SE = NEE SE)
  geom_ribbon(aes(ymin = cumNECB - 1.96 * cumNEE_se,
                  ymax = cumNECB + 1.96 * cumNEE_se),
              alpha = 0.15, colour = NA) +
  geom_line(aes(y = cumNEE),  linetype = "solid",  linewidth = 0.8, alpha = 0.25) +
  geom_line(aes(y = cumNECB), linetype = "solid",  linewidth = 0.8) +
  geom_vline(
    data = harvest_doy_by_year,
    aes(xintercept = doy, colour = management),
    linetype = "dotted", linewidth = 0.65, alpha = 0.75,
    inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.45) +
  facet_wrap(~ year, 
             # nrow = 2, 
             scales = "fixed") +
  scale_colour_manual(values = mgmt_colors, labels = mgmt_labels) +
  scale_fill_manual(values = mgmt_colors, labels = mgmt_labels) +
  labs(
    # title    = "Cumulative NEE vs. NECB by field and year",
    # subtitle = "Faint = NEE only | Bold = NECB (NEE + harvest C export) | Shaded = 95% CI | Dotted = harvest date",
    x        = "Day of year",
    y        = expression("Cumulative NECB (NEE + harvest C) (gC m"^{-2}*")"),
    colour   = NULL, fill = NULL
  ) +
  ylim(-410, 610) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

p_cumNECB
# --- 9c: Annual NEE bar chart -------------------------------------------------
# 2021 bars are hatched (alpha = 0.5) to signal partial-year coverage.
# A footnote annotation is added instead of overloading the legend.
p_annual_nee <- ggplot(annual_totals,
                       aes(x = factor(year), y = NEE_g_annual, fill = management,
                           alpha = ifelse(partial_year, 0.45, 1.0))) +
  geom_col(position = position_dodge(0.72), width = 0.65, colour = "grey30",
           linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = NEE_g_annual - 1.96 * NEE_annual_g_se,
        ymax = NEE_g_annual + 1.96 * NEE_annual_g_se),
    position = position_dodge(0.72), width = 0.25, linewidth = 0.6, alpha = 1
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = mgmt_colors, labels = mgmt_labels) +
  scale_alpha_identity() +
  labs(
    # title    = "Annual NEE",
    # subtitle = "Faded bars (2021) = partial-year observation (ends at final harvest date)",
    subtitle = "Faded bars (2018) = partial-year | Uncertainty from gap-filling SE only",
    x = "Year",
    y = expression("NEE (gC m"^{-2}*")"),
    fill  = NULL
  ) +
  # ylim(-250, 610) +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# --- 9d: Annual NECB bar chart -----------------------------------------------
p_annual_necb <- ggplot(annual_totals,
                        aes(x = factor(year), y = NECB, fill = management,
                            alpha = ifelse(partial_year, 0.45, 1.0))) +
  geom_col(position = position_dodge(0.72), width = 0.65, colour = "grey30",
           linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = NECB - 1.96 * NECB_se,
        ymax = NECB + 1.96 * NECB_se),
    position = position_dodge(0.72), width = 0.25, linewidth = 0.6, alpha = 1
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = mgmt_colors, labels = mgmt_labels) +
  scale_alpha_identity() +
  labs(
    # title    = "Annual NECB (NEE + harvest C export)",
    # subtitle = "Faded bars (2021) = partial-year | Uncertainty from gap-filling SE only",
    subtitle = "Faded bars (2018) = partial-year | Uncertainty from gap-filling SE only",
    x = element_blank(),
    y = expression("NECB (gC m"^{-2}*")"),
    fill  = NULL
  ) +
  # ylim(-250, 610) +
  theme_bw(base_size = 8) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

# --- 9e: Annual GPP and Reco -------------------------------------------------
annual_gpp_reco <- annual_totals %>%
  select(management, year,
         GPP_g_annual, GPP_annual_g_se,
         Reco_g_annual, Reco_annual_g_se) %>%
  pivot_longer(
    cols      = c(GPP_g_annual, Reco_g_annual),
    names_to  = "flux",
    values_to = "value"
  ) %>%
  mutate(
    se   = ifelse(flux == "GPP_g_annual", GPP_annual_g_se, Reco_annual_g_se),
    flux = recode(flux, "GPP_g_annual" = "GPP", "Reco_g_annual" = "Reco")
  )

p_annual_gpp_reco <- ggplot(annual_gpp_reco,
                            aes(x = factor(year), y = value, fill = management)) +
  geom_col(position = position_dodge(0.72), width = 0.65) +
  geom_errorbar(
    aes(ymin = value - 1.96 * se, ymax = value + 1.96 * se),
    position = position_dodge(0.72), width = 0.25, linewidth = 0.6
  ) +
  facet_wrap(~ flux, scales = "free_y") +
  scale_fill_manual(values = mgmt_colors, labels = mgmt_labels) +
  labs(
    title = "Annual GPP and Ecosystem Respiration",
    x = element_blank(),
    y = expression("Flux (gC m"^{-2}*" yr"^{-1}*")"),
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# --- 9f: Data availability figure --------------------------------------------
# Calendar-style raster showing data status for every day × field × year.
# Four categories:
#   Observed              — tower NEE accepted (no gap-filling applied)
#   Gap-filled (met OK)   — NEE rejected by QC; met sensors still online;
#                           GAM gap-filled from tower climate predictors
#   Gap-filled (ERA5 met) — complete tower outage; ERA5 climate used to
#                           enable GAM gap-filling
#   Unfilled              — gap remained after all gap-filling attempts
#                           (rare; typically first days of year before
#                           7-day rolling mean is available for Reco)
#
# Harvest events are marked as triangles on the lower edge of each year-row.
# Observation period for partial years (2018, 2021) is visible from where
# the tiles begin and end.

gap_vis <- data_cumul %>%
  mutate(
    # Classify each day into data-status categories
    status = case_when(
      !NEE_gapfilled                                      ~ "Observed",
      # NEE_gapfilled & airt_era5_filled & !is.na(NEE_filled) ~ "Gap-filled (ERA5 met)",
      # NEE_gapfilled & !airt_era5_filled & !is.na(NEE_filled) ~ "Gap-filled (met OK)",
      NEE_gapfilled & !is.na(NEE_filled) ~ "Gap-filled",
      NEE_gapfilled & is.na(NEE_filled)                   ~ "Unfilled",
      TRUE                                                ~ "Unfilled"
    ),
    status = factor(status,
                    levels = c("Observed",
                               # "Gap-filled (met OK)",
                               # "Gap-filled (ERA5 met)",
                               "Gap-filled",
                               "Unfilled")),
    field_label = factor(
      dplyr::recode(as.character(management),
                    "organic"      = "Organic",
                    "conventional" = "Conventional"),
      levels = c("Conventional", "Organic")
    ),
    # Reverse year order so 2018 is at top of each panel
    year_f = factor(year, levels = rev(sort(unique(year))))
  )

# Colour palette: accessible, prints legibly in greyscale
status_colours <- c(
  "Observed"              = "#2166AC",   # dark blue
  # "Gap-filled (met OK)"   = "#92C5DE",   # light blue
  # "Gap-filled (ERA5 met)" = "#F4A582",   # salmon
  "Gap-filled"            = "#92C5DE",
  "Unfilled"              = "#CFCFCF"    # light grey
)

p_data_avail <- gap_vis %>% 
  filter(year_f != 2021) %>% 
  ggplot(aes(x = doy, y = year_f, fill = status)) +
  geom_tile(height = 0.85, linewidth = 0) +
  # # Harvest triangles just below each year strip
  # geom_point(
  #   data = harvest_markers,
  #   aes(x = doy, y = year_f, shape = "Harvest"),
  #   inherit.aes = FALSE,
  #   colour = "grey20", size = 2, stroke = 0.4,
  #   position = position_nudge(y = -0.48)
  # ) +
  scale_fill_manual(values = status_colours,
                    name   = "Data status") +
  # scale_shape_manual(values = c("Harvest" = 25),   # filled downward triangle
  #                    name   = NULL) +
  scale_x_continuous(
    breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
    labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
    expand = c(0, 0)
  ) +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  facet_wrap(~ field_label, ncol = 1) +
  labs(
    title   = "NEE data availability by field and year",
    x       = NULL,
    y       = element_blank()
  ) +
  guides(
    fill  = guide_legend(order = 1, nrow = 2,
                         override.aes = list(size = 4, height = 0.7)),
    shape = guide_legend(order = 2,
                         override.aes = list(size = 3, colour = "grey20",
                                             fill  = "grey20"))
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position   = "bottom",
    # legend.box        = "vertical",
    legend.margin     = margin(t = 1),
    panel.grid        = element_blank(),
    strip.text        = element_text(face = "bold", size = 11),
    axis.ticks.x      = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1))

print(p_data_avail)

# --- Render all plots --------------------------------------------------------
print(p_cumNEE)
print(p_cumNECB)
print(p_annual_nee + p_annual_necb)   # side by side via patchwork
print(p_annual_gpp_reco)
print(p_gapfill_inspect)
print(p_data_avail)

# =============================================================================
# SECTION 10: Publication-ready summary table
# =============================================================================
#
# Columns:
#   Field            — management system (organic / conventional)
#   Year             — calendar year († = partial-year observation)
#   Crop             — crop(s) grown that year
#   Period (DOY)     — first and last day of observation
#   Days             — number of observation days
#   Gap (%)          — percentage of NEE days gap-filled
#   NEE              — annual cumulative NEE ± 95% CI (gC m⁻² yr⁻¹)
#   Harvest C export — carbon removed at harvest (gC m⁻² yr⁻¹)
#   NECB             — net ecosystem carbon balance ± 95% CI (gC m⁻² yr⁻¹)
#   C balance        — qualitative designation (sink / near-neutral / source)
#
# GPP and Reco are excluded from this table
#
# Uncertainty: ±1.96 × propagated gap-filling prediction SE (observed days
# contribute zero variance; harvest C treated as fixed/measured).
# Converted from µmol CO₂ m⁻² s⁻¹ × 86400 × 12.011 × 10⁻⁶ = gC m⁻² d⁻¹.

library(flextable)
library(officer)   # for Word export

# Qualitative C-balance designation based on annual NEE sign and magnitude.
# Threshold: |NEE| < 20 gC m⁻² = near-neutral (within typical measurement
# uncertainty for a partial-year or low-gap-fraction field-year).
designate_balance <- function(nee) {
  dplyr::case_when(
    nee < -20  ~ "Sink",
    nee >  20  ~ "Source",
    TRUE       ~ "Near-neutral"
  )
}

# Build observation-period dates from data_cumul for the period column
obs_period <- data_cumul %>%
  group_by(management, year) %>%
  summarise(
    date_start = format(min(date, na.rm = TRUE), "%d %b"),
    date_end   = format(max(date, na.rm = TRUE), "%d %b"),
    .groups = "drop"
  ) %>%
  mutate(period = paste0(date_start, "\u2013", date_end))   # en-dash

table_data <- annual_totals %>%
  left_join(obs_period, by = c("management", "year")) %>%
  mutate(
    # Field label
    Field = dplyr::recode(as.character(management),
                          "organic"       = "Organic",
                          "conventional"  = "Conventional"),
    # Year with dagger for partial years
    Year  = if_else(partial_year,
                    paste0(year, "\u2020"),   # †
                    as.character(year)),
    # Formatted NEE and NECB
    NEE_fmt  = sprintf("%.0f \u00b1 %.0f", NEE_g_annual, 1.96 * NEE_annual_g_se),
    NECB_fmt = sprintf("%.0f \u00b1 %.0f", NECB,       1.96 * NECB_se),
    C_exp    = sprintf("%.0f", C_harvested_gC_m2),
    Gap      = sprintf("%.0f", pct_gap_NEE),
    Days     = as.character(n_days),
    Balance  = designate_balance(NEE_g_annual)
  ) %>%
  arrange(Field, year) %>%
  select(
    Field, Year, Crop = crops, 
    # Period = period, 
    Days, `Gap-filled (%)` = Gap,
    `NEE` = NEE_fmt, `Harvest C` = C_exp, `NECB` = NECB_fmt,
    `C balance` = Balance
  )

table_data_NEEonly <- annual_totals %>%
  left_join(obs_period, by = c("management", "year")) %>%
  mutate(
    # Field label
    Field = dplyr::recode(as.character(management),
                          "organic"       = "Organic",
                          "conventional"  = "Conventional"),
    # Year with dagger for partial years
    Year  = if_else(partial_year,
                    paste0(year, "\u2020"),   # †
                    as.character(year)),
    # Formatted NEE and NECB
    NEE_fmt  = sprintf("%.0f \u00b1 %.0f", NEE_g_annual, 1.96 * NEE_annual_g_se),
    # NECB_fmt = sprintf("%.0f \u00b1 %.0f", NECB,       1.96 * NECB_se),
    # C_exp    = sprintf("%.0f", C_harvested_gC_m2),
    Gap      = sprintf("%.0f", pct_gap_NEE),
    Days     = as.character(n_days),
    Balance  = designate_balance(NEE_g_annual)
  ) %>%
  arrange(Field, year) %>%
  select(
    Field, Year, Crop = crops, 
    # Period = period, 
    Days, `Gap-filled (%)` = Gap,
    `NEE` = NEE_fmt, 
    # `Harvest C` = C_exp, 
    # `NECB` = NECB_fmt,
    `C balance` = Balance
  )

table_data_with_NECB <- annual_totals %>%
  left_join(obs_period, by = c("management", "year")) %>%
  mutate(
    # Field label
    Field = dplyr::recode(as.character(management),
                          "organic"       = "Organic",
                          "conventional"  = "Conventional"),
    # Year with dagger for partial years
    Year  = if_else(partial_year,
                    paste0(year, "\u2020"),   # †
                    as.character(year)),
    # Formatted NEE and NECB
    NEE_fmt  = sprintf("%.0f \u00b1 %.0f", NEE_g_annual, 1.96 * NEE_annual_g_se),
    NECB_fmt = sprintf("%.0f \u00b1 %.0f", NECB,       1.96 * NECB_se),
    C_exp    = sprintf("%.0f", C_harvested_gC_m2),
    Gap      = sprintf("%.0f", pct_gap_NEE),
    Days     = as.character(n_days),
    Balance  = designate_balance(NECB)
  ) %>%
  arrange(Field, year) %>%
  select(
    Field, Year, Crop = crops, 
    # Period = period, 
    Days, `Gap-filled (%)` = Gap,
    `NEE` = NEE_fmt, `Harvest C` = C_exp, `NECB` = NECB_fmt,
    `C balance with export` = Balance
  )

# Print plain version to console
cat("\n=== ANNUAL CARBON BALANCE SUMMARY TABLE ===\n")
cat("Units: gC m\u207b\u00b2. Uncertainty = \u00b11.96 \u00d7 propagated gap-fill SE.\n")
cat("\u2020 Partial-year observation (2018: late start; 2021: ends at final harvest).\n\n")
print(as.data.frame(table_data))

# --- Formatted flextable -----------------------------------------------------
ft <- flextable(table_data) %>%

  # Header labels
  set_header_labels(
    Field      = "Field",
    Year       = "Year",
    Crop       = "Crop",
    # Period     = "Observation period",
    Days       = "Days",
    `Gap-filled (%)`  = "Gap-filled (%)",
    NEE        = "NEE\n(gC m\u207b\u00b2)",
    `Harvest C`= "Harvest C\n(gC m\u207b\u00b2)",
    NECB       = "NECB\n(gC m\u207b\u00b2)",
    `C balance`= "C balance"
  ) %>%

  # Merge repeated Field cells (organic / conventional grouping)
  merge_v(j = "Field") %>%
  valign(j = "Field", valign = "top") %>%

  # Right-align numeric columns, left-align text
  align(j = c("NEE", "Harvest C", "NECB", "Days", "Gap-filled (%)"),
        align = "right", part = "all") %>%
  align(j = c("Field", "Year", "Crop", 
              # "Period", 
              "C balance"),
        align = "left", part = "all") %>%

  # Colour-code C balance column
  color(i = ~ `C balance` == "Sink",         j = "C balance", color = "#2E7D32") %>%
  color(i = ~ `C balance` == "Source",       j = "C balance", color = "#C62828") %>%
  color(i = ~ `C balance` == "Near-neutral", j = "C balance", color = "#757575") %>%
  bold(i = ~ `C balance` == "Sink",   j = "C balance") %>%
  bold(i = ~ `C balance` == "Source", j = "C balance") %>%

  # Light shading for partial-year rows
  # bg(i = ~ grepl("\u2018", Year), bg = "#F5F5F5") %>%
  
  # light shading for each management:
  bg(i = ~ grepl("Conventional", Field), bg = "#ffe5e1") %>% 
  bg(i = ~ grepl("Organic", Field), bg = "#fff0ce") %>% 

  # Column widths (inches, for Word)
  width(j = "Field",       width = 1.4) %>%
  width(j = "Year",        width = 0.45) %>%
  width(j = "Crop",        width = 1.3) %>%
  # width(j = "Period",      width = 1.1) %>%
  width(j = "Days",        width = 0.4) %>%
  width(j = "Gap-filled (%)",     width = 0.5) %>%
  width(j = "NEE",         width = 0.95) %>%
  width(j = "Harvest C",   width = 0.85) %>%
  width(j = "NECB",        width = 0.95) %>%
  width(j = "C balance",   width = 0.85) %>%

  # Borders
  border_outer(part = "all",    border = officer::fp_border(width = 1.2)) %>%
  border_inner_h(part = "body", border = officer::fp_border(width = 2,
                                                             color = "grey70")) %>%
  hline(part = "header",        border = officer::fp_border(width = 1.2)) %>%

  # Font
  # font(family = "serif") %>%
  fontsize(size = 10, part = "body") %>%
  fontsize(size = 10, part = "header") %>%
  bold(part = "header") %>%

  # Footnote
  add_footer_lines(paste0(
    "\u2020 Partial-year observations: 2018 data begins after tower installation; ",
    # "2021 data ends at final harvest date (organic: 10 Aug; conventional: 1 Oct). ",
    "Sink = net carbon uptake (NEE < \u221220 gC m\u207b\u00b2); ",
    "Source = net carbon loss (NEE > +20 gC m\u207b\u00b2); ",
    "Near-neutral = |\u200bNEE\u200b| \u2264 20 gC m\u207b\u00b2. ",
    "Uncertainty = \u00b11.96 \u00d7 propagated gap-filling prediction SE. ",
    "Harvest C export treated as fixed (no uncertainty). ",
    "NEE sign convention: negative = net carbon sink."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  italic(part = "footer") %>%
  color(part = "footer", color = "grey30") %>%

  set_table_properties(layout = "autofit")

# Define the text format with "Times New Roman"
# text_format_times <- fp_text(font.family = "Times New Roman", font.size = 10)

# Apply the font to the entire table
# ft <- style(ft, pr_t = text_format_times, part = "all")
ft

# NEE only: 
ft2 <- flextable(table_data_NEEonly) %>%
  
  # Header labels
  set_header_labels(
    Field      = "Field",
    Year       = "Year",
    Crop       = "Crop",
    # Period     = "Observation period",
    Days       = "Days",
    `Gap-filled (%)`  = "Gap-filled (%)",
    NEE        = "NEE\n(gC m\u207b\u00b2)",
    # `Harvest C`= "Harvest C\n(gC m\u207b\u00b2)",
    # NECB       = "NECB\n(gC m\u207b\u00b2)",
    `C balance`= "C balance"
  ) %>%
  
  # Merge repeated Field cells (organic / conventional grouping)
  merge_v(j = "Field") %>%
  valign(j = "Field", valign = "top") %>%
  
  # Right-align numeric columns, left-align text
  align(j = c("NEE", "Days", "Gap-filled (%)"),
        align = "right", part = "all") %>%
  align(j = c("Field", "Year", "Crop", "C balance"),
        align = "left", part = "all") %>%
  
  # Colour-code C balance column
  color(i = ~ `C balance` == "Sink",         j = "C balance", color = "#2E7D32") %>%
  color(i = ~ `C balance` == "Source",       j = "C balance", color = "#C62828") %>%
  color(i = ~ `C balance` == "Near-neutral", j = "C balance", color = "#757575") %>%
  bold(i = ~ `C balance` == "Sink",   j = "C balance") %>%
  bold(i = ~ `C balance` == "Source", j = "C balance") %>%
  
  # Light shading for partial-year rows
  # bg(i = ~ grepl("\u2020", Year), bg = "#F5F5F5") %>%
  
  # light shading for each management:
  bg(i = ~ grepl("Conventional", Field), bg = "#ffe5e1") %>% 
  bg(i = ~ grepl("Organic", Field), bg = "#fff0ce") %>% 
  
  # Column widths (inches, for Word)
  width(j = "Field",       width = 1.4) %>%
  width(j = "Year",        width = 0.45) %>%
  width(j = "Crop",        width = 1.3) %>%
  # width(j = "Period",      width = 1.1) %>%
  width(j = "Days",        width = 0.4) %>%
  width(j = "Gap-filled (%)",     width = 0.5) %>%
  width(j = "NEE",         width = 0.95) %>%
  # width(j = "Harvest C",   width = 0.85) %>%
  # width(j = "NECB",        width = 0.95) %>%
  width(j = "C balance",   width = 0.85) %>%
  
  # Borders
  border_outer(part = "all",    border = officer::fp_border(width = 1.2)) %>%
  border_inner_h(part = "body", border = officer::fp_border(width = 0.4,
                                                            color = "grey70")) %>%
  hline(part = "header",        border = officer::fp_border(width = 1.2)) %>%
  
  # Font
  # font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "body") %>%
  fontsize(size = 10, part = "header") %>%
  bold(part = "header") %>%
  
  # Footnote
  add_footer_lines(paste0(
    "\u2020 Partial-year observations: 2018 data begins after tower installation; ",
    # "2021 data ends at final harvest date (organic: 10 Aug; conventional: 1 Oct). ",
    "Sink = net carbon uptake (NEE < \u221220 gC m\u207b\u00b2); ",
    "Source = net carbon loss (NEE > +20 gC m\u207b\u00b2); ",
    "Near-neutral = |\u200bNEE\u200b| \u2264 20 gC m\u207b\u00b2. ",
    "Uncertainty = \u00b11.96 \u00d7 propagated gap-filling prediction SE. ",
    "NEE sign convention: negative = net carbon sink."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  italic(part = "footer") %>%
  color(part = "footer", color = "grey30") %>%
  
  set_table_properties(layout = "autofit")

print(ft2)

# with NECB: 
ft3 <- flextable(table_data_with_NECB) %>%
  
  # Header labels
  set_header_labels(
    Field      = "Field",
    Year       = "Year",
    Crop       = "Crop",
    Period     = "Observation period",
    Days       = "Days",
    `Gap-filled (%)`  = "Gap-filled (%)",
    NEE        = "NEE\n(gC m\u207b\u00b2)",
    `Harvest C`= "Harvest C\n(gC m\u207b\u00b2)",
    NECB       = "NECB\n(gC m\u207b\u00b2)",
    `C balance with export`= "C balance with export"
  ) %>%
  
  # Merge repeated Field cells (organic / conventional grouping)
  merge_v(j = "Field") %>%
  valign(j = "Field", valign = "top") %>%
  
  # Right-align numeric columns, left-align text
  align(j = c("NEE", "Harvest C", "NECB", "Days", "Gap-filled (%)"),
        align = "right", part = "all") %>%
  align(j = c("Field", "Year", "Crop", "Period", "C balance with export"),
        align = "left", part = "all") %>%
  
  # Colour-code C balance column
  color(i = ~ `C balance with export` == "Sink",         j = "C balance with export", color = "#2E7D32") %>%
  color(i = ~ `C balance with export` == "Source",       j = "C balance with export", color = "#C62828") %>%
  color(i = ~ `C balance with export` == "Near-neutral", j = "C balance with export", color = "#757575") %>%
  bold(i = ~ `C balance with export` == "Sink",   j = "C balance with export") %>%
  bold(i = ~ `C balance with export` == "Source", j = "C balance with export") %>%
  
  # Light shading for partial-year rows
  bg(i = ~ grepl("\u2020", Year), bg = "#F5F5F5") %>%
  
  # Column widths (inches, for Word)
  width(j = "Field",       width = 1.4) %>%
  width(j = "Year",        width = 0.45) %>%
  width(j = "Crop",        width = 1.3) %>%
  width(j = "Period",      width = 1.1) %>%
  width(j = "Days",        width = 0.4) %>%
  width(j = "Gap-filled (%)",     width = 0.5) %>%
  width(j = "NEE",         width = 0.95) %>%
  width(j = "Harvest C",   width = 0.85) %>%
  width(j = "NECB",        width = 0.95) %>%
  width(j = "C balance with export",   width = 0.85) %>%
  
  # Borders
  border_outer(part = "all",    border = officer::fp_border(width = 1.2)) %>%
  border_inner_h(part = "body", border = officer::fp_border(width = 0.4,
                                                            color = "grey70")) %>%
  hline(part = "header",        border = officer::fp_border(width = 1.2)) %>%
  
  # Font
  # font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "body") %>%
  fontsize(size = 10, part = "header") %>%
  bold(part = "header") %>%
  
  # Footnote
  add_footer_lines(paste0(
    "Sink, integrating harvest = net carbon uptake (NECB < \u221220 gC m\u207b\u00b2); ",
    "Source, integrating harvest = net carbon loss (NECB > +20 gC m\u207b\u00b2); ",
    "Near-neutral = |\u200bNEE\u200b| \u2264 20 gC m\u207b\u00b2. ",
    "Uncertainty = \u00b11.96 \u00d7 propagated gap-filling prediction SE. "
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  italic(part = "footer") %>%
  color(part = "footer", color = "grey30") %>%
  
  set_table_properties(layout = "autofit")

print(ft3)

# --- Export to Word ----------------------------------------------------------
# doc <- officer::read_docx() %>%
#   officer::body_add_par("Table 1. Annual carbon balance summary for organic (EF01) and conventional (EF02) fields, 2018\u20132021.",
#                          style = "Normal") %>%
#   officer::body_add_par("", style = "Normal") %>%   # spacer
#   flextable::body_add_flextable(ft)
# 
# output_path <- here::here("R scripts", "annual_carbon_balance_table.docx")
# print(doc, target = output_path)
# cat("\nTable exported to:", output_path, "\n")
