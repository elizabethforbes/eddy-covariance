# this script is a comparison of the model outcomes with and without the addition of
# end-of-2020 empirical mean daily flux data from the organic field. In the initial 
# round of analyses, we had excluded the data labeled as "alfalfa hay", though in reality
# this was a grass mix dominated by orchardgrass but also containing brome and alfalfa.
# This grass mix was undersown in the barley crop, such that it grew out and "popped"
# after barley was harvested in August 2020. In mid-September 2020, the field was 
# additionally drill-seeded with orchardgrass to fill it out. The field went dormant in winter, 
# and in 2021 the cash crop was alfalfa-dominated hay (harvested twice, though these data
# were not included in the original model comparison because by 2021 the cash crop rotations
# across the two fields no longer matched, as the conventional field went back to soy).

library(tidyverse)
library(mgcv)
library(dplyr)
library(gratia)
library(DHARMa)
library(zoo) 

################################################################################
#  define data
################################################################################

# dat1: the data without anything labeled alfalfa
dat1 <- ec_daily_mngmnt3 %>% 
  filter(year == 2018 | year == 2019 | year == 2020) %>% # restrict years to 2018, 2019, and 2020
  filter(is.na(current_crop) | current_crop != "alfalfa hay")
# 1742 individual mean daily observations of flux; broken down as follows
# 2018 conventional   246
# 2018 organic        141
# 2019 conventional   365
# 2019 organic        365
# 2020 conventional   366
# 2020 organic        259

# make sure management is in factor form
dat1$management <- factor(as.character(dat1$management))
# same for crop stage:
dat1$crop_stage_simple <- factor(dat1$crop_stage_simple)
# set reference level crop stage as "mature":
dat1$crop_stage_simple <- relevel(dat1$crop_stage_simple, ref = "mature")

# dat2: the data with the alfalfa 2020 rows
dat2 <- ec_daily_mngmnt3 %>% 
  filter(year == 2018 | year == 2019 | year == 2020) # restrict years to 2018, 2019, and 2020
  # keep alfalfa rows from 2020 organic
# 1849 individual mean daily observations of flux; broken down as follows
# 2018 conventional   246
# 2018 organic        141
# 2019 conventional   365
# 2019 organic        365
# 2020 conventional   366
# 2020 organic        366 # the only difference, with addition of 107 days from mid-Sept to EOY

dat2$management <- factor(as.character(dat2$management))
dat2$crop_stage_simple <- factor(dat2$crop_stage_simple)
dat2$crop_stage_simple <- relevel(dat2$crop_stage_simple, ref = "mature")

################################################################################
# model building
################################################################################

# mod1: the accepted model structure without the 2020 organic alfalfa rows
# Mark the first row of each management-year block (resets AR process)
dat1 <- dat1 %>%
  dplyr::arrange(management, year, doy) %>%
  dplyr::group_by(management, year) %>%
  dplyr::mutate(AR.start = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# convert year to factor to avoid linear function of year:
dat1$year_f <- factor(dat1$year)

mod1 <- 
  bam(
  NEE ~
    management * crop_stage_simple +
    s(year_f, bs = "re", by = as.factor(management)) +
    s(doy, bs = "cc", k = 20) +
    te(air_temperature, VPD,
       by = as.factor(management),
       k  = c(8, 8)) +
    s(days_since_tillage, by = as.factor(management), k = 8),
  data     = dat1,
  method   = "fREML",
  rho      = 0.41,           # phi from gamm() corAR1 estimate
  AR.start = dat1$AR.start
)

appraise(mod1)

# mod2: the accepted model structure WITH the 2020 organic alfalfa rows
# Mark the first row of each management-year block (resets AR process)
dat2 <- dat2 %>%
  dplyr::arrange(management, year, doy) %>%
  dplyr::group_by(management, year) %>%
  dplyr::mutate(AR.start = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# convert year to factor to avoid linear function of year:
dat2$year_f <- factor(dat2$year)

mod2 <- 
  bam(
    NEE ~
      management * crop_stage_simple +
      s(year_f, bs = "re", by = as.factor(management)) +
      s(doy, bs = "cc", k = 20) +
      te(air_temperature, VPD,
         by = as.factor(management),
         k  = c(8, 8)) +
      s(days_since_tillage, by = as.factor(management), k = 8),
    data     = dat2,
    method   = "fREML",
    rho      = 0.41,           # phi from gamm() corAR1 estimate
    AR.start = dat2$AR.start
  )

appraise(mod2)

################################################################################
# gap-filling, air temp and VPD
################################################################################

# prep data:
dat1_cumul <- dat1 %>%
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

dat2_cumul <- dat2 %>%
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

# introduce VPD calculation
# VPD calculation functions (in kPa)
calc_es <- function(temp_C) {
  0.6108 * exp((17.27 * temp_C) / (temp_C + 237.3))
}

calc_ea <- function(es, RH) {
  (RH / 100) * es
}

calc_vpd <- function(temp_C, RH) {
  es <- calc_es(temp_C)
  ea <- calc_ea(es, RH)
  es - ea
}

# calculate corrected ERA5 air temp and VPD values:
dat1_filled <- dat1 %>% 
  mutate(
    # ERA5-derived VPD (kPa):
    VPD_era5_calc = calc_vpd(airt_era5, RH_era5),
    
    # bias-corrected ERA5 values for both air temperature and VPD:
    airt_era5_corrected   = coef(era5_correction_airt)[1]+ # intercept
      coef(era5_correction_airt)[2]*airt_era5, # slope*value
    VPD_era5_corrected    = coef(era5_correction_VPD)[1]+
      coef(era5_correction_VPD)[2]*VPD_era5_calc,
    # convert VPD in kPa back to Pa to match VPD in model structure:
    VPD_era5_corrected_Pa = VPD_era5_corrected * 1000 ,
    
    # gap filling:
    air_temperature_filled = if_else(is.na(air_temperature),
                                     airt_era5_corrected,
                                     air_temperature), # if air_temperature is NA, fill with corrected ERA5 temp
    VPD_filled = if_else(is.na(VPD),
                         VPD_era5_corrected_Pa,
                         VPD)
  )
dat2_filled <- dat2 %>% 
  mutate(
    # ERA5-derived VPD (kPa):
    VPD_era5_calc = calc_vpd(airt_era5, RH_era5),
    
    # bias-corrected ERA5 values for both air temperature and VPD:
    airt_era5_corrected   = coef(era5_correction_airt)[1]+ # intercept
      coef(era5_correction_airt)[2]*airt_era5, # slope*value
    VPD_era5_corrected    = coef(era5_correction_VPD)[1]+
      coef(era5_correction_VPD)[2]*VPD_era5_calc,
    # convert VPD in kPa back to Pa to match VPD in model structure:
    VPD_era5_corrected_Pa = VPD_era5_corrected * 1000 ,
    
    # gap filling:
    air_temperature_filled = if_else(is.na(air_temperature),
                                     airt_era5_corrected,
                                     air_temperature), # if air_temperature is NA, fill with corrected ERA5 temp
    VPD_filled = if_else(is.na(VPD),
                         VPD_era5_corrected_Pa,
                         VPD)
  )
dat1_filled <- dat1_filled %>%
  select(management, date, air_temperature_filled, VPD_filled)
dat2_filled <- dat2_filled %>%
  select(management, date, air_temperature_filled, VPD_filled)

# merge:
dat1_cumul <- dat1_cumul %>%
  left_join(dat1_filled, by = c("management", "date")) %>%
  mutate(
    airt_era5_filled = is.na(air_temperature) & !is.na(air_temperature_filled),
    VPD_era5_filled  = is.na(VPD)             & !is.na(VPD_filled),
    air_temperature  = coalesce(air_temperature, air_temperature_filled),
    VPD              = coalesce(VPD,             VPD_filled)
  ) %>%
  select(-air_temperature_filled, -VPD_filled)
dat2_cumul <- dat2_cumul %>%
  left_join(dat2_filled, by = c("management", "date")) %>%
  mutate(
    airt_era5_filled = is.na(air_temperature) & !is.na(air_temperature_filled),
    VPD_era5_filled  = is.na(VPD)             & !is.na(VPD_filled),
    air_temperature  = coalesce(air_temperature, air_temperature_filled),
    VPD              = coalesce(VPD,             VPD_filled)
  ) %>%
  select(-air_temperature_filled, -VPD_filled)

# Re-compute 7-day rolling means now that climate predictors have fewer NAs.
# (Reco model uses air_temperature_7 and VPD_7; more filled days = more
#  gap-rows that can receive Reco predictions.)
dat1_cumul <- dat1_cumul %>%
  arrange(management, year, doy) %>%
  group_by(management) %>%
  mutate(
    air_temperature_7 = zoo::rollmean(air_temperature, 7, fill = NA, align = "right"),
    VPD_7             = zoo::rollmean(VPD,             7, fill = NA, align = "right")
  ) %>%
  ungroup()
dat2_cumul <- dat2_cumul %>%
  arrange(management, year, doy) %>%
  group_by(management) %>%
  mutate(
    air_temperature_7 = zoo::rollmean(air_temperature, 7, fill = NA, align = "right"),
    VPD_7             = zoo::rollmean(VPD,             7, fill = NA, align = "right")
  ) %>%
  ungroup()

# random effect terms of year:
get_re_terms <- function(model) {
  labels <- sapply(model$smooth, function(s) s$label)
  labels[grepl("year_f", labels)]
}

re_terms_mod1  <- get_re_terms(mod1)
re_terms_mod2 <- get_re_terms(mod2)

################################################################################
# gap-filling, NEE
###############################################################################

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

# apply function to both models:
nee_mod1_gf <- gapfill_flux(dat1_cumul, mod1, "NEE", re_terms_mod1)
nee_mod2_gf <- gapfill_flux(dat2_cumul, mod2, "NEE", re_terms_mod2)

# append to dataframe:
dat1_cumul$NEE_filled    <- nee_mod1_gf$flux_filled
dat1_cumul$NEE_se        <- nee_mod1_gf$flux_se
dat1_cumul$NEE_gapfilled <- nee_mod1_gf$is_gapfilled

dat2_cumul$NEE_filled    <- nee_mod2_gf$flux_filled
dat2_cumul$NEE_se        <- nee_mod2_gf$flux_se
dat2_cumul$NEE_gapfilled <- nee_mod2_gf$is_gapfilled

# =============================================================================
# Harvest carbon export
# =============================================================================
# Conversion formula:
#   C_export (gC m⁻²) = yield × kg_per_unit × moisture_corr × c_fraction
#                       / 4046.86 (m² per acre) × 1000 (kg → g)

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

UMOL_TO_GC_PER_DAY <- 86400 * 12.011 * 1e-6   # = 1.037750
UMOL_TO_KGC_PER_DAY <- 86400 * 12.011 * 1e-6 * 0.001 # = 0.00103775

dat1_cumul <- dat1_cumul %>%
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
    NEE_g_se      = NEE_se      * UMOL_TO_GC_PER_DAY,
    NEE_filled_kg  = NEE_filled  * UMOL_TO_KGC_PER_DAY,
    NEE_kg_se      = NEE_se      * UMOL_TO_KGC_PER_DAY,
  )
dat2_cumul <- dat2_cumul %>%
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
    NEE_g_se      = NEE_se      * UMOL_TO_GC_PER_DAY,
    NEE_filled_kg  = NEE_filled  * UMOL_TO_KGC_PER_DAY,
    NEE_kg_se      = NEE_se      * UMOL_TO_KGC_PER_DAY,
  )

dat1_cumul <- dat1_cumul %>% 
  filter(date >= "2018-08-13") # removed 105 observations from conventional field
dat2_cumul <- dat2_cumul %>% 
  filter(date >= "2018-08-13") # removed 105 observations from conventional field
# so the two partial year start dates match (for conv. and org.) for accumulation calc

# cumulative sums:
dat1_cumul_daily <- dat1_cumul %>%
  arrange(management, year, doy) %>%
  group_by(management, year) %>%
  mutate(
    # Running cumulative sums
    cumNEE_g   = cumsum(replace_na(NEE_filled_g,  0)),
    cumNEE_kg  = cumsum(replace_na(NEE_filled_kg,  0)),
    # Propagated uncertainty (SE of cumulative sum).
    # replace_na(..., 0) prevents a single NA predictor in one gap row from
    # propagating through cumsum() and wiping out all subsequent SE values.
    # Rows where predict() returned NA (unfillable gaps) contribute 0 variance.
    cumNEE_g_se  = sqrt(cumsum(replace_na(NEE_g_se^2,  0))),
    cumNEE_kg_se  = sqrt(cumsum(replace_na(NEE_kg_se^2,  0))),
    # Running gap-fill count
    n_days         = row_number(),
    cum_gapfill_NEE  = cumsum(as.integer(NEE_gapfilled)),
  ) %>%
  ungroup()

dat2_cumul_daily <- dat2_cumul %>%
  arrange(management, year, doy) %>%
  group_by(management, year) %>%
  mutate(
    # Running cumulative sums
    cumNEE_g   = cumsum(replace_na(NEE_filled_g,  0)),
    cumNEE_kg  = cumsum(replace_na(NEE_filled_kg,  0)),
    # Propagated uncertainty (SE of cumulative sum).
    # replace_na(..., 0) prevents a single NA predictor in one gap row from
    # propagating through cumsum() and wiping out all subsequent SE values.
    # Rows where predict() returned NA (unfillable gaps) contribute 0 variance.
    cumNEE_g_se  = sqrt(cumsum(replace_na(NEE_g_se^2,  0))),
    cumNEE_kg_se  = sqrt(cumsum(replace_na(NEE_kg_se^2,  0))),
    # Running gap-fill count
    n_days         = row_number(),
    cum_gapfill_NEE  = cumsum(as.integer(NEE_gapfilled)),
  ) %>%
  ungroup()

# Annual totals: take the last row of each management-year
annual_totals_dat1 <- dat1_cumul_daily %>%
  group_by(management, year) %>%
  summarise(
    n_days          = n(),
    first_doy       = min(doy),
    last_doy        = max(doy),
    partial_year    = (min(doy) > 14) | (max(doy) < 351),  # late start OR early end
    n_gap_NEE       = sum(NEE_gapfilled),
    pct_gap_NEE     = round(100 * n_gap_NEE  / n_days, 1),
    # in g:
    NEE_g_annual      = sum(replace_na(NEE_filled_g,  0)),
    NEE_annual_g_se   = sqrt(sum(replace_na(NEE_g_se^2, 0))), # there are 15 NAs in this row
    # in kg:
    NEE_kg_annual      = sum(replace_na(NEE_filled_kg,  0)),
    NEE_annual_kg_se   = sqrt(sum(NEE_kg_se^2)),
    .groups         = "drop"
  )
annual_totals_dat2 <- dat2_cumul_daily %>%
  group_by(management, year) %>%
  summarise(
    n_days          = n(),
    first_doy       = min(doy),
    last_doy        = max(doy),
    partial_year    = (min(doy) > 14) | (max(doy) < 351),  # late start OR early end
    n_gap_NEE       = sum(NEE_gapfilled),
    pct_gap_NEE     = round(100 * n_gap_NEE  / n_days, 1),
    # in g:
    NEE_g_annual      = sum(replace_na(NEE_filled_g,  0)),
    NEE_annual_g_se   = sqrt(sum(replace_na(NEE_g_se^2, 0))), # there are 15 NAs in this row
    # in kg:
    NEE_kg_annual      = sum(replace_na(NEE_filled_kg,  0)),
    NEE_annual_kg_se   = sqrt(sum(NEE_kg_se^2)),
    .groups         = "drop"
  )

# NECB totals:
annual_totals_dat1 <- annual_totals_dat1 %>%
  left_join(harvest_by_year %>%
              select(management, year, C_harvested_gC_m2, crops,
                     harvest_dates, any_date_estimated),
            by = c("management", "year")) %>%
  mutate(
    # Add C export; propagate NEE SE only (yield treated as fixed/measured)
    NECB    = NEE_g_annual + C_harvested_gC_m2,
    NECB_se = NEE_annual_g_se
  )
annual_totals_dat2 <- annual_totals_dat2 %>% 
  left_join(harvest_by_year %>%
              select(management, year, C_harvested_gC_m2, crops,
                     harvest_dates, any_date_estimated),
            by = c("management", "year")) %>%
  mutate(
    # Add C export; propagate NEE SE only (yield treated as fixed/measured)
    NECB    = NEE_g_annual + C_harvested_gC_m2,
    NECB_se = NEE_annual_g_se
  )

# for comparison: write to csv
write.csv(dat1_cumul_daily, "no_alfalfa_dat1.csv")
write.csv(dat2_cumul_daily, "yes_alfalfa_dat2.csv")
write.csv(annual_totals_dat1, "no_alfalfa_annualtotals.csv")
write.csv(annual_totals_dat2, "yes_alfalfa_annualtotals.csv")

# TAKEAWAYS FROM COMPARISON BETWEEN THE OLD-GAP-FILLED DATASET AND THE NEW GAP-FILLED
# DATASET: new one is parameterized using a model that included empirical training
# data from 2020 in the organic field from mid-September through the end of the year,
# and reflecting a grass mix field that was undersown in the barley field and grew up before
# entering dormancy and transitioning into a majority alfalfa hay field in 2021.

# The resulting gap-filling changed the 2019 predictions substantially, such that
# the subsequent *2019* organic cumulative NEE was much more negative than it had been prior
# to the inclusion of those 2020 data in the organic field. This shift (-208.7 gC m⁻²) is
# real, but requires careful framing:

# The new model is predicting substantially less source behavior during non-
# growing-season periods in the organic field. The shift is concentrated in 
# dormant (-74.6 gC m⁻²), mature (-70.3 gC m⁻²), and vegetative (-47.2 gC m⁻²) 
# stages, and by month it's heaviest in April (-46 gC m⁻²), October (-50 gC m⁻²), 
# and November-December (-60 gC m⁻²) — exactly the off-season periods where the
# new fall 2020 organic observations now provide training data that previously 
# didn't exist.

# the old model was extrapolating the organic field's fall/winter 
# behavior from conventional field data and the shared seasonal smooth, with no 
# actual organic field observations from those periods. It predicted larger 
# source conditions (old dormant: +1.59 µmol m⁻² s⁻¹, old mature: +0.58) than 
# what the new observations show those conditions actually look like (new 
# dormant: +0.87, new mature: −0.55). The new training data is correcting a 
# systematic over-prediction of source conditions in the organic field during 
# non-peak-season periods.

# IMPORTANT FRAMING: the new training data comes from a grass mix 
# establishing after barley harvest — which is not the same as a rye cover crop 
# establishing after soybean harvest (what actually happened in fall 2019 in the organic field). 
# The model can't distinguish those contexts; it sees "organic, vegetative/dormant, 
# November" and predicts accordingly. So the 2019 correction is directionally 
# defensible but grounded in an imperfect analog. It is likely necessary to note 
# in the methods that the organic field's off-season periods before 2020 rely on 
# predictions that were informed by the 2020 grass-mix establishment phase as the closest
# analog available from the empirical dataset.

# The change in cumulative NEE for 2020 is clean: 106 days of additional real observations 
# averaging +0.57 µmol m⁻² s⁻¹ added directly to the cumulative sum.

# table comparing empirical data, mod1 data, and mod2 data across growth stages:
################################################################################
# SUMMARY: Mean daily NEE by management × year × crop stage
# empirical = observed values only (gap-filled days excluded)
# mod1      = NEE_filled from model without 2020 fall organic data
# mod2      = NEE_filled from model with 2020 fall organic data
# units: µmol m⁻² s⁻¹
# Note: rows present only in dat2 (2020 organic Sep17–Dec31) will have
#       NA for mod1 — these are flagged via mod1_only_na
################################################################################

nee_stage_comparison <- dat2_cumul_daily %>%
  left_join(
    dat1_cumul_daily %>%
      select(management, year, doy, NEE_filled) %>%
      rename(NEE_filled_mod1 = NEE_filled),
    by = c("management", "year", "doy")
  ) %>%
  group_by(management, year, crop_stage_simple) %>%
  summarise(
    n_days             = n(),
    n_observed         = sum(!is.na(NEE)),
    pct_gapfilled_mod2 = round(100 * sum(NEE_gapfilled) / n(), 1),
    n_mod1_na          = sum(is.na(NEE_filled_mod1)),  # days in mod2 but not mod1
    mean_NEE_empirical = round(mean(NEE,              na.rm = TRUE), 3),
    mean_NEE_mod1      = round(mean(NEE_filled_mod1,  na.rm = TRUE), 3),
    mean_NEE_mod2      = round(mean(NEE_filled,       na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  mutate(
    mod2_minus_mod1 = round(mean_NEE_mod2 - mean_NEE_mod1, 3)
  ) %>%
  arrange(management, year, crop_stage_simple)

print(nee_stage_comparison, n = Inf)

# write.csv(nee_stage_comparison, "nee_stage_model_comparison.csv", row.names = FALSE)