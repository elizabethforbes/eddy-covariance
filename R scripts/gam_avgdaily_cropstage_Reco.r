# =============================================================================
# GAM models: Management effects on Ecosystem Respiration (Reco) (EC tower data)

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

gam_Reco_additive <- gam(
  Reco ~
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
# Rationale: hot+dry days (high temp, high VPD) may suppress Reco differently
# than cool+dry or hot+humid conditions. The ti() terms decompose the
# interaction cleanly, leaving the main effects (s() terms) interpretable.
# -----------------------------------------------------------------------------

gam_Reco_tensor <- gam(
  Reco ~
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
# Rationale: hot+dry days (high temp, high VPD) may suppress Reco differently
# than cool+dry or hot+humid conditions. The te()) term does NOT decompose the
# interaction, so the main effects of temp and VPD are not interpretable, but 
# this structure reduces the risk of concurvity (aka when a smooth term can be
# approximated by a combo of others).
# -----------------------------------------------------------------------------

gam_Reco_tensor2 <- gam(
  Reco ~
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
gam_Reco_additive_ml <- update(gam_Reco_additive, method = "ML")
gam_Reco_tensor_ml   <- update(gam_Reco_tensor,   method = "ML")
gam_Reco_tensor2_ml <- update(gam_Reco_tensor2, method = "ML")

AIC(gam_GPP_additive_ml, gam_Reco_tensor_ml, gam_Reco_tensor2_ml)
#                         df      AIC
# gam_GPP_additive_ml 49.39249 4651.388
# gam_Reco_tensor_ml  64.18081 3188.126
# gam_Reco_tensor2_ml 63.88245 3188.036

# Likelihood ratio test: tensor2 (te) significantly better than tensor (ti)
anova(gam_Reco_tensor, gam_Reco_tensor2)
# no p-value due to negative deviance aka no significant improvement in fit

# check concurvity of ti() model:
concurvity(gam_Reco_tensor, full = TRUE)
# Concurvity check: ti() model rejected in favour of te()
# Air temperature and VPD show high concurvity in the decomposed ti() model
# (observed > 0.97 for marginal smooths), consistent with their strong
# covariance in the Hudson Valley. Marginal effects are not independently
# identifiable; joint effect modelled with te() instead.

# Summary and diagnostics (gam_NEE with full tensor interaction — pre-autocorrelation correction)
summary(gam_Reco_tensor2)

# checking residuals:
appraise(gam_Reco_tensor2)
simulateResiduals(gam_Reco_tensor2, plot = TRUE)
# possible quantile deviations detected; possible issue with autocorrelation
acf(residuals(gam_Reco_tensor2), main = "ACF of GAM residuals") 
# significantly high ACF at lag-1: >0.80
# daily EC is inherently autocorrelated, i.e., subsequent days' fluxes determined partly
# by previous days' fluxes.

# -----------------------------------------------------------------------------
# Investigate temporal autocorrelation:
# -----------------------------------------------------------------------------

# add AR(1) correlation structure with gamm():
gamm_Reco <- gamm(
  Reco ~
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
coef(gamm_Reco$lme$modelStruct$corStruct, unconstrained = FALSE)
# Phi 
# 0.9947313  
# interpret: ~99% of today's residual is carried into the next day's. Previous (gam()) model 
# was underestimating SE without autocorrelation accounting.
acf(residuals(gamm_Reco$lme, type = "normalized"), main = "ACF of gamm LME normalized residuals") 
# lag is not totally absorbed, indicating better model fit but not perfect. ACF > 0.50 now.

# -----------------------------------------------------------------------------
# Autocorrelation correction: bam() with AR(1)
#
# ACF (autocorrelation function) of gam_GPP_tensor2 residuals showed significant lag-1 
# autocorrelation (>0.80), confirmed by gamm() corAR1 estimate of phi = 0.99. 

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

bam_Reco <- bam(
  Reco ~
    management * crop_stage_simple +
    s(year_f, bs = "re", by = as.factor(management)) +
    s(doy, bs = "cc", k = 20) +
    te(air_temperature, VPD,
       by = as.factor(management),
       k  = c(8, 8)) +
    s(days_since_tillage, by = as.factor(management), k = 8),
  data     = data1,
  method   = "fREML",
  rho      = 0.99,           # phi from gamm() corAR1 estimate
  AR.start = data1$AR.start
)

appraise(bam_Reco) # WOOF: this is way worse!

# Check ACF of bam residuals -- these DO NOT reflect the AR(1) correction directly, so we whiten manually
acf(residuals(bam_Reco), main = "ACF of bam_Reco residuals") # correlated residuals; as expected, high
# plot the decorrelated residuals as addressed in the model structure:
r <- residuals(bam_Reco)

r_white <- numeric(length(r))
r_white[1] <- r[1]
for (i in 2:length(r)) {
  if (data1$AR.start[i]) {
    r_white[i] <- r[i]        # reset at start of each management-year block
  } else {
    r_white[i] <- r[i] - 0.42 * r[i - 1]
  }
}
acf(r_white, main = "ACF of whitened bam_Reco residuals") # decorrelated residuals: within the range
# ACF still really high at lag-1; around .99 still.

# -----------------------------------------------------------------------------
# Dealing with autocorrelation: Ecosystem Respiration
# Reco is more susceptible to past environmental conditions: months or weeks
# current autocorrelation structure assumes day-to-day autocorrelation, so no
# lagged or cumulative variables
# therefore: residuals remain highly autocorrelated
# -----------------------------------------------------------------------------

# attempt 1: calculate lagged air temp and VPD variables (7 days) to assess longer autocorrelation
data1$air_temperature_7 <- zoo::rollmean(data1$air_temperature, 7, fill = NA, align = "right")
data1$VPD_7 <- zoo::rollmean(data1$VPD, 7, fill = NA, align = "right")

bam_Reco2 <- bam(
  Reco ~
    management * crop_stage_simple +
    s(year_f, bs = "re", by = as.factor(management)) +
    s(doy, bs = "cc", k = 20) +
    te(air_temperature, VPD,
       by = as.factor(management),
       k  = c(8, 8)) +
    te(air_temperature_7, VPD_7,
       by = as.factor(management),
       k  = c(8, 8)) +
    s(days_since_tillage, by = as.factor(management), k = 8),
  data     = data1,
  method   = "fREML",
  AR.start = data1$AR.start
)

appraise(bam_Reco2) # better observed vs. fitted and residuals vs. linear, more outliers in QQ

# Check ACF of bam residuals -- these DO NOT reflect the AR(1) correction directly, so we whiten manually
acf(residuals(bam_Reco2), main = "ACF of bam_Reco2 (with climate lags) residuals") # correlated residuals; lag-1 = ~.78
# plot the decorrelated residuals as addressed in the model structure:
r <- residuals(bam_Reco2)
r_white[1] <- r[1]
for (i in 2:length(r)) {
  if (data1$AR.start[i]) {
    r_white[i] <- r[i]        # reset at start of each management-year block
  } else {
    r_white[i] <- r[i] - 0.42 * r[i - 1]
  }
}
acf(r_white, main = "ACF of whitened bam_Reco2 (with climate lags) residuals") # decorrelated residuals: within the range
# ACF still really high at lag-1; around .7 still.

# attempt 2: use more complex, higher-order structure for ACF, like corARMA
# AR(1) assumes correlation only between adjacent residuals decaying exponentially.
# ARMA(1,1) adds a moving average term that can model short-term shocks or noise, 
# improving fit for residuals with more complex temporal patterns.
# Ecosystem respiration residuals often have both long memory (AR) and 
# short-term noise (MA) components due to biological and measurement variability.
library(nlme)
gamm_Reco2 <- gamm(
# gamm_Reco3 <- gamm(
# gamm_Reco4 <- gamm(
  Reco ~
    management * crop_stage_simple +
    s(year_f, bs = "re", by = as.factor(management)) +
    s(doy, bs = "cc", k = 20) +
    te(air_temperature, VPD,
       by = as.factor(management),
       k  = c(8, 8)) +
    s(days_since_tillage, by = as.factor(management), k = 8),
  data = data1,
  correlation = corARMA(form = ~ 1 | management/year_f, p = 1, q = 1),
  # correlation = corARMA(form = ~ 1 | year_f, p = 1, q = 1),
  method = "REML"
  # method = "ML"
)

# gamm_Reco3 = same model with interaction in corARMA, just ML structure
# gamm_Reco4 = also ML structure, without interaction in corARMA
# AIC(gamm_Reco3, gamm_Reco4)
# # df      AIC
# # gamm_Reco3$lme 34 434.1785
# # gamm_Reco4$lme 34 431.4870 indicates slightly better fit without interaction

coef(gamm_Reco2$lme$modelStruct$corStruct, unconstrained = FALSE)
# Phi 
# 0.9934581  
# interpret: ~99% of today's residual is carried into the next day's. Previous (gam()) model 
# was underestimating SE without autocorrelation accounting.
acf(residuals(gamm_Reco2$lme, type = "normalized"), main = "ACF of gamm LME normalized residuals") 
# lag-1 is mostly absorbed, however!

# attempt 3: Ecosystem respiration responds to slow-changing drivers (smoothed 
# climate variables) but also has inherent biological memory and process noise 
# that cause residual autocorrelation. Use both corARMA structure and 7-day
# smoothed air temp and VPD variables.

gamm_Reco_smoothed <- gamm(
  Reco ~
    management * crop_stage_simple +
    s(year_f, bs = "re", by = as.factor(management)) +
    s(doy, bs = "cc", k = 20) +
    te(air_temperature_7, VPD_7,
       by = as.factor(management),
       k  = c(8, 8)) +
    s(days_since_tillage, by = as.factor(management), k = 8),
  data = data1,
  correlation = corARMA(form = ~ 1 | management/year_f, p = 1, q = 1),
  method = "REML"
  # method = "ML"
)
# investigate phi (the AR(1) autocorrelation structure value)
coef(gamm_Reco_smoothed$lme$modelStruct$corStruct, unconstrained = FALSE)
# Phi 
# 0.989817  
# interpret: ~99% of today's residual is carried into the next day's. Previous (gam()) model 
# was underestimating SE without autocorrelation accounting.
acf(residuals(gamm_Reco_smoothed$lme, type = "normalized"), main = "ACF of gamm LME normalized residuals") 
# lag is absorbed! no significant remaining temporal dependence.

# assess model fit:
resid_norm <- residuals(gamm_Reco_smoothed$lme, type = "normalized")
plot(fitted(gamm_Reco_smoothed$lme), resid_norm,
     xlab = "fitted values", ylab = "normalized residuals", main = "Residuals vs. Fitted")
abline(h = 0, col = "red", lty = 5) # some patterns of scatter around 7 but generally okay

qqnorm(resid_norm) 
qqline(resid_norm, col = "red", lty = 5) # a bit off, indicating outliers
hist(resid_norm, breaks = 30, main = "Histogram of Normalized Residuals",
     xlab = "Residual Value") 
# however the residuals are fairly normally distributed
# outliers (see above, qq plot) are apparent at -6 and 5

# final autocorrelation estimated parameters with ARMA (1,1) structure:
coef(gamm_Reco_smoothed$lme$modelStruct$corStruct, unconstrained = FALSE)
# Phi1    Theta1 
# 0.989817 0.482016  

# phi = influence of previous time step's values on this one. 0.99 indicates strong persistence/memory
# in Reco from day to day, suggesting gradual changes in Reco over time and general temporal stability

# theta = moving average, indicates relationship between this value and past error terms, aka random shocks
# unexpected deviations in Reco have moderate influence on current Reco value; short term environmental
# disturbances will affect it but impact will fade pretty fast

# calculate a pseudo R2 for this model since R2 is not useful here as it reflects two different
# correlation structures (lm and gam) PLUS the corARMA structure makes it impossible to 
# account for autocorrelation in predicted values

library(DescTools)
# use CCC (concordance correlation coefficient) bc it's not possible to calculate R2
CCC(obs, fits_lme)$rho.c
#    est    lwr.ci    upr.ci
# 0.6307068 0.5946755 0.6642096

# -----------------------------------------------------------------------------
# Summary tables: Reco (gamm model with corARMA structure and smoothed climate)
# -----------------------------------------------------------------------------

library(gt)

# -- Parametric coefficients (fixed effects) ----------------------------------
# Note: tbl_regression() labels by variable name, not coefficient name, so
# we use broom::tidy() + gt() directly for full control over interaction labels.
tbl_reco_parametric <- broom::tidy(gamm_Reco_smoothed$gam, parametric = TRUE) %>%
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
    title    = md("**Table 1. Parametric coefficients: Reco GAMM with 7day averaged climate variables**"),
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
tbl_reco_smooths <- broom::tidy(gamm_Reco_smoothed$gam, parametric = FALSE) %>%
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
    title    = md("**Table 2. Smooth terms: Reco GAMM with 7day averaged climate variables**"),
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

tbl_reco_parametric
tbl_reco_smooths

# -----------------------------------------------------------------------------
# Visualizations
# -----------------------------------------------------------------------------

# source(here::here("R scripts", "gam_NEE_vizfunctions.r"))
# source(here::here("R scripts", "gam_NEE_vizfunctions.r"))
source(here::here("R scripts", "plot_databymanagement.r"))
source(here::here("R scripts", "plot_te_difference.r"))


# Main finding figure ---------------------------------------------------------
# g1 <- plot_coef_plot(bam_GPP, flux_var = "GPP")

# Main finding figure ---------------------------------------------------------
# r1 <- plot_coef_plot(gamm_Reco_smoothed$gam, flux_var = "Reco")

# Raw observed means, by management x stage — no climate correction -----------
r2 <- plot_management_comparison(data1, flux_var = "Reco")
r2 + labs( 
  y = (bquote(respiration ~ "(µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")")),
  title = "Respiration by Management and Crop Stage")

# difference plot for te() surface: effects of climate on VARIABLE as indicated
# by the difference in organic - conventional
r_te <- plot_te_difference(
  model           = gamm_Reco_smoothed,
  data            = data1,
  flux_var        = "R[eco]",
  temp_var        = "air_temperature_7",
  climate_var     = "VPD_7",
  better_direction = "neutral",
  layout          = "temp_vpd"
)
r_te


# Predicted by stage ------------------------------------------------------
# r3 <- plot_predicted_by_stage(gamm_Reco_smoothed$gam, data1, flux_var = "Reco")
# corrected for climate (VPD and air temperature)

# Plot the te() surface: effects of VPD*temperature on NEE across management --
# r4 <- plot_te_surface(gamm_Reco_smoothed$gam,
                      # temp_var = "air_temperature_7", vpd_var = "VPD_7",
                      # data1, flux_var = "Reco")
# organic maintains higher atmospheric fixation through primary productivity
# in hot, dry conditions, more so than the conventional field. Its atmospheric 
# fixation potential is much higher at cold, humid levels than conventional; probably
# reflects the lack of fallow conditions in winter. in extremely hot, fairly dry conditions
# (temp > 20C, VPD < 1000) the organic field maintains a positive productivity level
# while the conventional field is non-productive.

# # add corner labels:
# r4 + 
#   annotate("label", x = -4, y = .25,  label = "cool & humid air\n(low demand)",
#               color = "grey20", size = 2.6, label.size = 0, fill = alpha("white", 0.75)) +
#   annotate("label", x = 20,  y = .25,  label = "warm & humid air\n(productive)",
#            color = "grey20", size = 2.6, label.size = 0, fill = alpha("white", 0.75)) +
#   annotate("label", x = -4, y = 1.3, label = "cool & dry air\n(uncommon)",
#            color = "grey20", size = 2.6, label.size = 0, fill = alpha("white", 0.75)) +
#   annotate("label", x = 20,  y = 1.3, label = "hot & dry air\n(drought stress)",
#            color = "grey20", size = 2.6, label.size = 0, fill = alpha("white", 0.75))

# -- Model performance --------------------------------------------------------
# r5 <- plot_gam_observed_vs_predicted(gamm_Reco_smoothed$gam, data1, flux_var = "Reco")

# Residuals by crop stage: check for systematic stage-level bias
# r6 <- plot_gam_residuals_by_stage(gamm_Reco_smoothed$gam, data1, flux_var = "Reco")
