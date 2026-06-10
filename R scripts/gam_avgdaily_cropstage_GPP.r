# =============================================================================
# GAM models: Management effects on GPP (EC tower data)

library(mgcv)
library(dplyr)
library(gratia)
library(DHARMa)

# if data is not already prepped (e.g., data1 produced in gam_avgdaily_cropstage_NEE.r),
# put management and crop stage in factor form, and set "mature" as the reference stage
# trim dataset to only the years 2018-2020, and remove any alfalfa crop rows

# (see "gam_avgdaily_cropstage_NEE.r" for reference)

# -----------------------------------------------------------------------------
# Model 1: Additive -- separate smooths for air_temperature and VPD
# -----------------------------------------------------------------------------

gam_GPP_additive <- gam(
  GPP ~
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
# Rationale: hot+dry days (high temp, high VPD) may suppress GPP differently
# than cool+dry or hot+humid conditions. The ti() terms decompose the
# interaction cleanly, leaving the main effects (s() terms) interpretable.
# -----------------------------------------------------------------------------

gam_GPP_tensor <- gam(
  GPP ~
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

gam_GPP_tensor2 <- gam(
  GPP ~
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
gam_GPP_additive_ml <- update(gam_GPP_additive, method = "ML")
gam_GPP_tensor_ml   <- update(gam_GPP_tensor,   method = "ML")
gam_GPP_tensor2_ml <- update(gam_GPP_tensor2, method = "ML")

AIC(gam_GPP_additive_ml, gam_GPP_tensor_ml, gam_GPP_tensor2_ml)
#                         df      AIC
# gam_GPP_additive_ml 48.94357 5082.339
# gam_GPP_tensor_ml   63.80241 5043.723
# gam_GPP_tensor2_ml  64.04496 5046.048

# Likelihood ratio test: tensor2 (te) NOT significantly better than tensor (ti)
anova(gam_GPP_tensor, gam_GPP_tensor2)
# Resid. Df Resid. Dev     Df Deviance      F Pr(>F)
# 1    1060.8     5027.6                              
# 2    1059.0     5024.0 1.7726   3.6889 0.4452 0.6169

# check concurvity of ti() model:
concurvity(gam_GPP_tensor, full = TRUE)
# Concurvity check: ti() model rejected in favour of te()
# Air temperature and VPD show high concurvity in the decomposed ti() model
# (observed > 0.98 for marginal smooths), consistent with their strong
# covariance in the Hudson Valley. Marginal effects are not independently
# identifiable; joint effect modelled with te() instead.

# Summary and diagnostics (gam_NEE with full tensor interaction — pre-autocorrelation correction)
summary(gam_GPP_tensor2)

# checking residuals:
appraise(gam_GPP_tensor2)
simulateResiduals(gam_GPP_tensor2, plot = TRUE)
# possible quantile deviations detected; possible issue with autocorrelation
acf(residuals(gam_GPP_tensor2), main = "ACF of GAM residuals") 
# significantly high ACF at lag-1: >0.45
# daily EC is inherently autocorrelated, i.e., subsequent days' fluxes determined partly
# by previous days' fluxes.

# -----------------------------------------------------------------------------
# Investigate temporal autocorrelation:
# -----------------------------------------------------------------------------

# add AR(1) correlation structure with gamm():
gamm_GPP <- gamm(
  GPP ~
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
coef(gamm_GPP$lme$modelStruct$corStruct, unconstrained = FALSE)
# Phi 
# 0.6006473  
# interpret: ~60% of today's residual is carried into the next day's. Previous (gam()) model 
# was underestimating SE without autocorrelation accounting.
acf(residuals(gamm_GPP$lme, type = "normalized"), main = "ACF of gamm LME normalized residuals") 
# lag is absorbed, indicating better model fit

# -----------------------------------------------------------------------------
# Autocorrelation correction: bam() with AR(1)
#
# ACF (autocorrelation function) of gam_GPP_tensor2 residuals showed significant lag-1 
# autocorrelation (>0.45), confirmed by gamm() corAR1 estimate of phi = 0.60. 

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

bam_GPP <- bam(
  GPP ~
    management * crop_stage_simple +
    s(year_f, bs = "re", by = as.factor(management)) +
    s(doy, bs = "cc", k = 20) +
    te(air_temperature, VPD,
       by = as.factor(management),
       k  = c(8, 8)) +
    s(days_since_tillage, by = as.factor(management), k = 8),
  data     = data1,
  method   = "fREML",
  rho      = 0.61,           # phi from gamm() corAR1 estimate
  AR.start = data1$AR.start
)

appraise(bam_GPP)

# Check ACF of bam residuals -- these DO NOT reflect the AR(1) correction directly, so we whiten manually
acf(residuals(bam_GPP), main = "ACF of bam_GPP residuals") # correlated residuals; as expected, high
# plot the decorrelated residuals as addressed in the model structure:
r <- residuals(bam_GPP)

r_white <- numeric(length(r))
r_white[1] <- r[1]
for (i in 2:length(r)) {
  if (data1$AR.start[i]) {
    r_white[i] <- r[i]        # reset at start of each management-year block
  } else {
    r_white[i] <- r[i] - 0.42 * r[i - 1]
  }
}
acf(r_white, main = "ACF of whitened bam_GPP residuals") # decorrelated residuals: within the range

summary(bam_GPP)

# -----------------------------------------------------------------------------
# Summary tables: GPP
# -----------------------------------------------------------------------------

library(gt)

# -- Parametric coefficients (fixed effects) ----------------------------------
# Note: tbl_regression() labels by variable name, not coefficient name, so
# we use broom::tidy() + gt() directly for full control over interaction labels.
tbl_gpp_parametric <- broom::tidy(bam_GPP, parametric = TRUE) %>%
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
    title    = md("**Table 1. Parametric coefficients: GPP GAM (bam_GPP)**"),
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
tbl_gpp_smooths <- broom::tidy(bam_GPP, parametric = FALSE) %>%
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
    title    = md("**Table 2. Smooth terms: GPP GAM (bam_GPP)**"),
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

tbl_gpp_parametric
tbl_gpp_smooths

# -----------------------------------------------------------------------------
# Visualizations
# -----------------------------------------------------------------------------

# source(here::here("R scripts", "gam_NEE_vizfunctions.r"))
source(here::here("R scripts", "plot_databymanagement.r"))
source(here::here("R scripts", "plot_te_difference.r"))


# Main finding figure ---------------------------------------------------------
# g1 <- plot_coef_plot(bam_GPP, flux_var = "GPP")

# Raw observed means, by management x stage — no climate correction -----------
g2 <- plot_management_comparison(data1, flux_var = "GPP")
g2  + coord_cartesian(ylim = c(-5, 20))    # before saving / returning

# difference plot for te() surface: effects of climate on VARIABLE as indicated
# by the difference in organic - conventional

g_te <- plot_te_difference(model           = bam_GPP,
                           data            = data1,
                           flux_var        = "GPP",
                           temp_var        = "air_temperature",
                           climate_var     = "VPD",
                           better_direction = "positive",
                           layout          = "temp_vpd")
g_te
