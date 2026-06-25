# =============================================================================
# GAM models: Management effects on NEE (EC tower data)
#
# Two model structures compared:
#   Model 1 (additive):  separate smooths for air_temperature and VPD
#   Model 2 (tensor):    ti(air_temperature, VPD) interaction term added
#
# Dropped from previous version:
#   - soilt_avg (ERA5-modeled, not empirical)
#   - vwc_avg   (ERA5-modeled, not empirical)
#   - ti(soilt_avg, vwc_avg, ...) tensor interaction
#
# Retained:
#   - management * crop_stage_simple (fixed interaction)
#   - s(year, bs = "re", by = management) (random effect of year per management)
#   - s(doy, bs = "cc", k = 20) (shared seasonal cycle, cyclical cubic spline)
#   - s(air_temperature, by = management, k = 10)
#   - s(VPD, by = management, k = 10)
#   - s(days_since_tillage, by = management, k = 8)
#
# Note on concurvity: air_temperature and VPD covary strongly (VPD is a
# function of temp and RH). Check post-fit with concrvity(model, full = TRUE).
# Values > 0.8 for any term indicate a potential problem. RH is excluded
# entirely for this reason -- VPD already integrates humidity information.
# 
# Model comparison: use AIC and REML score. Note that REML scores are only
# directly comparable between models with identical fixed effect structures;
# use AIC (ML-fitted) for comparing models 1 vs 2.
# =============================================================================

library(mgcv)
library(dplyr)
library(gratia)
library(DHARMa)

# make sure management is in factor form
ec_daily_mngmnt3$management <- factor(as.character(ec_daily_mngmnt3$management))
# same for crop stage:
ec_daily_mngmnt3$crop_stage_simple <- factor(ec_daily_mngmnt3$crop_stage_simple)
# set reference level crop stage as "mature":
ec_daily_mngmnt3$crop_stage_simple <- relevel(ec_daily_mngmnt3$crop_stage_simple, ref = "mature")

# -----------------------------------------------------------------------------
# Mechanisms of difference in fluxes: matched data across 2018-2020
# -----------------------------------------------------------------------------

# Because 2021 represented a transition from annual grains to perennial in the organic
# field, I am removing these data from both fields for the direct analysis of average
# daily fluxes by management. The question: how do the systems respond differently to climate drivers?
# this question is better answered by temporally and biophysically-matched datasets.

data1 <- ec_daily_mngmnt3 %>%
  filter(year == 2018 | year == 2019 | year == 2020) # restrict years to 2018, 2019, and 2020
# data1 <- data1 %>% 
#   filter(is.na(current_crop) | current_crop != "alfalfa hay")

# Note: rather than use NDVI or some other indicator of crop development (esp
# across different crops), crop_stage_simple is used as a categorical predictor.

# I have chosen to include an interaction between management type and crop stage
# in the linear predictor variables as it is reasonable to expect that management
# type will result in different outcomes at each crop stage IRT metabolism and thus fluxes.

# -----------------------------------------------------------------------------
# Model 1: Additive -- separate smooths for air_temperature and VPD
# -----------------------------------------------------------------------------

gam_NEE_additive <- gam(
  NEE ~
    management *
    crop_stage_simple +
    
    # smoothing factors:
    s(year, bs = "re", by = as.factor(management)) +   # year random effect, per management
    s(doy, bs = "cc", k = 20) +                        # shared seasonal cycle
    s(air_temperature, by = as.factor(management), k = 10) +
    s(VPD,             by = as.factor(management), k = 10) +
    s(days_since_tillage, by = as.factor(management), k = 8),

  data   = data1,
  method = "REML"
)

# -----------------------------------------------------------------------------
# Model 2: Tensor interaction -- ti(air_temperature, VPD) added
#
# Rationale: hot+dry days (high temp, high VPD) may suppress NEE differently
# than cool+dry or hot+humid conditions. The ti() terms decompose the
# interaction cleanly, leaving the main effects (s() terms) interpretable.
# -----------------------------------------------------------------------------

gam_NEE_tensor <- gam(
  NEE ~
    management *
    crop_stage_simple +

    s(year, bs = "re", by = as.factor(management)) +
    s(doy, bs = "cc", k = 20) +

    s(air_temperature, by = as.factor(management), k = 10) +
    s(VPD,             by = as.factor(management), k = 10) +

    ti(air_temperature, VPD,
       by = as.factor(management),
       k  = c(8, 8)) +

    s(days_since_tillage, by = as.factor(management), k = 8),

  data   = data1,
  method = "REML"
)

# -----------------------------------------------------------------------------
# Model 3: FULL tensor product -- te(air_temperature, VPD) added, and removal
# of the independent air temp and VPD terms; thus the marginal smooths are 
# encompassed by the tensor product smooth.
#
# Rationale: hot+dry days (high temp, high VPD) may suppress NEE differently
# than cool+dry or hot+humid conditions. The te()) term does NOT decompose the
# interaction, so the main effects of temp and VPD are not interpretable, but 
# this structure reduces the risk of concurvity (aka when a smooth term can be
# approximated by a combo of others).
# -----------------------------------------------------------------------------

gam_NEE_tensor2 <- gam(
  NEE ~
    management *
    crop_stage_simple +
    s(year, bs = "re", by = as.factor(management)) +
    s(doy, bs = "cc", k = 20) +
    te(air_temperature, VPD,
       by = as.factor(management),
       k  = c(8, 8)) +
    s(days_since_tillage, by = as.factor(management), k = 8),
  data   = data1,
  method = "REML"
)

# -----------------------------------------------------------------------------
# Model comparison
# -----------------------------------------------------------------------------

# Refit with ML for valid AIC comparison (REML scores not comparable across
# models with different smooth/fixed structures). We will use the REML-fitted 
# models for interpretation, but direct comparison of goodness-of-fit will
# not be meaningful, hence re-fitting here.
gam_NEE_additive_ml <- update(gam_NEE_additive, method = "ML")
gam_NEE_tensor_ml   <- update(gam_NEE_tensor,   method = "ML")
gam_NEE_tensor2_ml <- update(gam_NEE_tensor2, method = "ML")

AIC(gam_NEE_additive_ml, gam_NEE_tensor_ml, gam_NEE_tensor2_ml)
#                         df      AIC
# gam_NEE_additive_ml 51.15988 4725.101
# gam_NEE_tensor_ml   60.10648 4657.745
# gam_NEE_tensor2_ml  65.01893 4655.843

# check concurvity of ti() model: basically, can one smooth term be approximated by a combination of the other smooths in the model?
concurvity(gam_NEE_tensor, full = TRUE)
# Concurvity check: ti() model rejected in favour of te()
# Air temperature and VPD show high concurvity in the decomposed ti() model
# (observed > 0.99 for marginal smooths), consistent with their strong
# covariance in the Hudson Valley. Marginal effects are not independently
# identifiable; joint effect should be modeled with te() instead.

# Summary and diagnostics (gam_NEE with full tensor interaction — pre-autocorrelation correction)
summary(gam_NEE_tensor2)

# checking residuals:
appraise(gam_NEE_tensor2)
simulateResiduals(gam_NEE_tensor2, plot = TRUE)
# possible quantile deviations detected; possible issue with autocorrelation
acf(residuals(gam_NEE_tensor2), main = "ACF of GAM residuals") 
# significantly high ACF at lag-1
# daily EC is inherently autocorrelated, i.e., subsequent days' fluxes determined partly
# by previous days' fluxes.

# add AR(1) correlation structure with gamm():
gamm_NEE <- gamm(
  NEE ~
    management * crop_stage_simple +
    s(year, bs = "re", by = as.factor(management)) +
    s(doy, bs = "cc", k = 20) +
    te(air_temperature, VPD, by = as.factor(management), k = c(8, 8)) +
    s(days_since_tillage, by = as.factor(management), k = 8),
  data        = data1,
  method      = "REML",
  correlation = corAR1(form = ~ 1 | management/year) # AR1 autocorr structure that resets by year within each field
)
# investigate phi (the AR(1) autocorrelation structure value)
coef(gamm_NEE$lme$modelStruct$corStruct, unconstrained = FALSE)
# Phi 
# 0.4110967  
# interpret: 41% of today's residual is carried into the next day's. 
# Previous (gam()) model was underestimating SE without autocorrelation accounting.
acf(residuals(gamm_NEE$lme, type = "normalized"), main = "ACF of gamm LME normalized residuals") 
# lag is absorbed, indicating better model fit

# -----------------------------------------------------------------------------
# Autocorrelation correction: bam() with AR(1)
#
# ACF (autocorrelation function) of gam_NEE_tensor2 residuals showed significant lag-1 
# autocorrelation (>0.30), confirmed by gamm() corAR1 estimate of phi = 0.41. 

# Because the data contains gaps (tower outages), corAR1 within gamm() does not account
# for unequal time spacing within years. bam() with rho + AR.start handles
# this explicitly: the AR(1) is reset at the start of each management-year
# block, so correlation does not bleed across year boundaries or gaps.
# -----------------------------------------------------------------------------

# Mark the first row of each management-year block (resets AR process)
data1 <- data1 %>%
  dplyr::arrange(management, year, doy) %>%
  dplyr::group_by(management, year) %>%
  dplyr::mutate(AR.start = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# convert year to factor to avoid linear function of year:
data1$year_f <- factor(data1$year)

bam_NEE <- bam(
  NEE ~
    management * crop_stage_simple +
    s(year_f, bs = "re", by = as.factor(management)) +
    s(doy, bs = "cc", k = 20) +
    te(air_temperature, VPD,
       by = as.factor(management),
       k  = c(8, 8)) +
    s(days_since_tillage, by = as.factor(management), k = 8),
  data     = data1,
  method   = "fREML",
  rho      = 0.41,           # phi from gamm() corAR1 estimate
  AR.start = data1$AR.start
)

appraise(bam_NEE)

# Check ACF of bam residuals -- these DO NOT reflect the AR(1) correction directly, so we whiten manually
acf(residuals(bam_NEE), main = "ACF of bam_NEE residuals") # correlated residuals; as expected, high
# plot the decorrelated residuals as addressed in the model structure:
r <- residuals(bam_NEE)

r_white <- numeric(length(r))
r_white[1] <- r[1]
for (i in 2:length(r)) {
  if (data1$AR.start[i]) {
    r_white[i] <- r[i]        # reset at start of each management-year block
  } else {
    r_white[i] <- r[i] - 0.42 * r[i - 1]
  }
}
acf(r_white, main = "ACF of whitened bam_NEE residuals") # decorrelated residuals: within the range

summary(bam_NEE)

# -----------------------------------------------------------------------------
# Summary tables
# -----------------------------------------------------------------------------

library(gt)

# -- Parametric coefficients (fixed effects) ----------------------------------
# Note: tbl_regression() labels by variable name, not coefficient name, so
# we use broom::tidy() + gt() directly for full control over interaction labels.
tbl_nee_parametric <- broom::tidy(bam_NEE, parametric = TRUE) %>%
  dplyr::mutate(
    term = dplyr::case_when(
      term == "(Intercept)"                                      ~ "Intercept (conventional, mature)",
      term == "managementorganic"                                ~ "Organic vs. Conventional",
      term == "crop_stage_simpledormant"                         ~ "Dormant",
      term == "crop_stage_simpleearly"                           ~ "Early growth",
      term == "crop_stage_simplegrain_fill"                      ~ "Grain fill",
      term == "crop_stage_simplereproductive"                    ~ "Reproductive",
      term == "crop_stage_simplevegetative"                      ~ "Vegetative",
      term == "crop_stage_simplefallow"                          ~ "Fallow",
      term == "managementorganic:crop_stage_simpledormant"       ~ "Organic \u00d7 Dormant",
      term == "managementorganic:crop_stage_simpleearly"         ~ "Organic \u00d7 Early growth",      
      term == "managementorganic:crop_stage_simplefallow"         ~ "Organic \u00d7 Fallow",      
      term == "managementorganic:crop_stage_simplegrain_fill"    ~ "Organic \u00d7 Grain fill",
      term == "managementorganic:crop_stage_simplereproductive"  ~ "Organic \u00d7 Reproductive",
      term == "managementorganic:crop_stage_simplevegetative"    ~ "Organic \u00d7 Vegetative",
      TRUE ~ term
    ),
    sig = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE            ~ ""
    )
  ) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value, sig) %>%
  gt() %>%
  cols_label(
    term      = "Term",
    estimate  = "Estimate",
    std.error = "Std. Error",
    statistic = "t",
    p.value   = "p-value",
    sig       = ""
  ) %>%
  fmt_number(columns = c(estimate, std.error, statistic), decimals = 3) %>%
  fmt(
    columns = p.value,
    fns = function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  ) %>%
  tab_header(
    title    = md("**Table 1. Parametric coefficients: NEE GAM (bam_NEE)**"),
    subtitle = md("Reference: conventional management, mature crop stage")
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = everything(), rows = p.value < 0.05)
  )

# -- Smooth terms (approximate significance) ----------------------------------
tbl_nee_smooths <- broom::tidy(bam_NEE, parametric = FALSE) %>%
  dplyr::mutate(
    term = dplyr::case_when(
      grepl("s\\(year_f\\).*conventional", term) ~ "s(year): conventional",
      grepl("s\\(year_f\\).*organic",      term) ~ "s(year): organic",
      grepl("s\\(doy\\)",                term) ~ "s(doy) \u2014 shared seasonal cycle",
      grepl("te.*conventional",          term) ~ "te(Temp \u00d7 VPD): conventional",
      grepl("te.*organic",               term) ~ "te(Temp \u00d7 VPD): organic",
      grepl("days_since.*conventional",  term) ~ "s(days since tillage): conventional",
      grepl("days_since.*organic",       term) ~ "s(days since tillage): organic",
      TRUE ~ term
    )
  ) %>%
  dplyr::select(term, edf, ref.df, statistic, p.value) %>%
  gt() %>%
  cols_label(
    term      = "Smooth term",
    edf       = "EDF",
    ref.df    = "Ref. df",
    statistic = "F",
    p.value   = "p-value"
  ) %>%
  fmt_number(columns = c(edf, ref.df, statistic), decimals = 2) %>%
  fmt(
    columns = p.value,
    fns = function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  ) %>%
  tab_header(
    title    = md("**Table 2. Smooth terms: NEE GAM (bam_NEE)**"),
    subtitle = md("Approximate significance; EDF = effective degrees of freedom")
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = p.value, rows = p.value < 0.05)
  )

tbl_nee_parametric
tbl_nee_smooths

# Residual diagnostics (gratia)
appraise(bam_NEE)

# DHARMa residuals: tests for overdispersion, uniformity, outliers
simulationOutput <- DHARMa::simulateResiduals(bam_NEE, n = 500)
plot(simulationOutput)
# some issues here: indications of heteroscedasticity

gam.check(bam_NEE)
# however, the histogram is clean as are the response vs. fitted values.
# the distribution of residuals is very peaked, but that's ok; and, 
# the model converged quickly (12 iterations) and the k checks indicate good
# (aka not overfitting) the smooths, and the te() terms are only using 14 of 63 EDFs
# (aka, also not overfitting)

# -----------------------------------------------------------------------------
# Visualizations
# -----------------------------------------------------------------------------

# source(here::here("R scripts", "gam_NEE_vizfunctions.r"))
source(here::here("R scripts", "plot_databymanagement.r"))
source(here::here("R scripts", "plot_te_difference.r"))

# Main finding figure ---------------------------------------------------------
# p1 <- plot_coef_plot(bam_NEE, flux_var = "NEE")
# p1

# Raw observed means, by management x stage — no climate correction -----------
p2 <- plot_management_comparison(data1, flux_var = "NEE")
p2 
# + coord_cartesian(ylim = c(-10, 10))    # before saving / returning

# difference plot for te() surface: effects of climate on VARIABLE as indicated
# by the difference in organic - conventional
n_te <- plot_te_difference(
  model           = bam_NEE,
  data            = data1,
  flux_var        = "NEE",
  temp_var        = "air_temperature",
  climate_var     = "VPD",
  better_direction = "negative",
  layout          = "temp_vpd"
)
n_te
