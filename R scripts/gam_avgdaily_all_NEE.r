# =============================================================================
# GAM models: Management effects on NEE (EC tower data) for all crop stages,
# with management as the only driver
#
# Two model structures compared:
#   Model 1 (additive):  separate smooths for air_temperature and VPD
#   Model 2 (tensor):    ti(air_temperature, VPD) interaction term added
#
# Retained:
#   - management (fixed interaction)
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
ec_daily_mngmnt$management <- factor(as.character(ec_daily_mngmnt$management))
# make sure year is in factor form
ec_daily_mngmnt$year <- factor(ec_daily_mngmnt$year)

# -----------------------------------------------------------------------------
# Model 1: Additive -- separate smooths for air_temperature and VPD
# -----------------------------------------------------------------------------

gam_NEE_additive <- gam(
  NEE ~
    management +

    # smoothing factors:
    s(year, bs = "re", by = as.factor(management)) +   # year random effect, per management
    s(doy, bs = "cc", k = 20) +                        # shared seasonal cycle
    s(air_temperature, by = as.factor(management), k = 10) +
    s(VPD,             by = as.factor(management), k = 10) +
    s(days_since_tillage, by = as.factor(management), k = 8),

  data   = ec_daily_mngmnt,
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
    management +
    
    s(year, bs = "re", by = as.factor(management)) +
    s(doy, bs = "cc", k = 20) +

    s(air_temperature, by = as.factor(management), k = 10) +
    s(VPD,             by = as.factor(management), k = 10) +

    ti(air_temperature, VPD,
       by = as.factor(management),
       k  = c(8, 8)) +

    s(days_since_tillage, by = as.factor(management), k = 8),

  data   = ec_daily_mngmnt,
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
    management +
    s(year, bs = "re", by = as.factor(management)) +
    s(doy, bs = "cc", k = 20) +
    te(air_temperature, VPD,
       by = as.factor(management),
       k  = c(8, 8)) +
    s(days_since_tillage, by = as.factor(management), k = 8),
  data   = ec_daily_mngmnt,
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
# gam_NEE_additive_ml 51.25718 7127.791
# gam_NEE_tensor_ml   63.40903 7085.167
# gam_NEE_tensor2_ml  62.35945 7080.500

# check concurvity of ti() model:
concurvity(gam_NEE_tensor, full = TRUE)
# para s(year):as.factor(management)conventional s(year):as.factor(management)organic    s(doy)
# worst       1                                 1.0000000                            1.0000000 0.9648129
# observed    1                                 0.6687036                            0.7211138 0.8011951
# estimate    1                                 0.8079585                            0.7995640 0.3519846
# s(air_temperature):as.factor(management)conventional s(air_temperature):as.factor(management)organic
# worst                                               0.9999993                                       0.9999995
# observed                                            0.9705408                                       0.9678078
# estimate                                            0.9746772                                       0.9864982
# s(VPD):as.factor(management)conventional s(VPD):as.factor(management)organic
# worst                                   0.9999997                           0.9999993
# observed                                0.9953478                           0.9993279
# estimate                                0.9882211                           0.9917520
# ti(air_temperature,VPD):as.factor(management)conventional ti(air_temperature,VPD):as.factor(management)organic
# worst                                                    0.9999997                                            0.9999997
# observed                                                 0.9314694                                            0.9375607
# estimate                                                 0.3902418                                            0.3813842
# s(days_since_tillage):as.factor(management)conventional s(days_since_tillage):as.factor(management)organic
# worst                                                  0.9386042                                          0.9999244
# observed                                               0.6193592                                          0.8223721
# estimate                                               0.6093005                                          0.8325202
# VPD and air temperature show high concurvity in the decomposed ti() model (observed > 0.97 for marginal smooths)

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
    management +
    s(year, bs = "re", by = as.factor(management)) +
    s(doy, bs = "cc", k = 20) +
    te(air_temperature, VPD, by = as.factor(management), k = c(8, 8)) +
    s(days_since_tillage, by = as.factor(management), k = 8),
  data        = ec_daily_mngmnt,
  method      = "REML",
  correlation = corAR1(form = ~ 1 | management/year) # AR1 autocorr structure that resets by year within each field
)
# investigate phi (the AR(1) autocorrelation structure value)
coef(gamm_NEE$lme$modelStruct$corStruct, unconstrained = FALSE)
# Phi 
# 0.6280484  
# interpret: 63% of today's residual is carried into the next day's. Previous (gam()) model was underestimating SE without autocorrelation accounting.
acf(residuals(gamm_NEE$lme, type = "normalized"), main = "ACF of gamm LME normalized residuals") 
# lag is absorbed, indicating better model fit

# -----------------------------------------------------------------------------
# Model 4 and 5: Autocorrelation correction: bam() with AR(1)
#
# ACF (autocorrelation function) of gam_NEE_tensor2 residuals showed significant lag-1 
# autocorrelation (>0.30), confirmed by gamm() corAR1 estimate of phi = 0.63 

# Because the data contains gaps (tower outages), corAR1 within gamm() does not account
# for unequal time spacing within years. bam() with rho + AR.start handles
# this explicitly: the AR(1) is reset at the start of each management-year
# block, so correlation does not bleed across year boundaries or gaps.
# -----------------------------------------------------------------------------

# Mark the first row of each management-year block (resets AR process)
ec_daily_mngmnt <- ec_daily_mngmnt %>%
  dplyr::arrange(management, year, doy) %>%
  dplyr::group_by(management, year) %>%
  dplyr::mutate(AR.start = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

bam_NEE <- bam(
  NEE ~
    management +
    s(year, bs = "re", by = as.factor(management)) +
    s(doy, bs = "cc", k = 20) +
    te(air_temperature, VPD,
       by = as.factor(management),
       k  = c(8, 8)) +
    s(days_since_tillage, by = as.factor(management), k = 8),
  data     = ec_daily_mngmnt,
  method   = "fREML",
  rho      = 0.63,           # phi from gamm() corAR1 estimate
  AR.start = ec_daily_mngmnt$AR.start
)

appraise(bam_NEE)

# Check ACF of bam residuals -- these DO NOT reflect the AR(1) correction directly, so we whiten manually
acf(residuals(bam_NEE), main = "ACF of bam_NEE residuals") # correlated residuals; as expected, high
# plot the decorrelated residuals as addressed in the model structure:
r <- residuals(bam_NEE)

r_white <- numeric(length(r))
r_white[1] <- r[1]
for (i in 2:length(r)) {
  if (ec_daily_mngmnt$AR.start[i]) {
    r_white[i] <- r[i]        # reset at start of each management-year block
  } else {
    r_white[i] <- r[i] - 0.42 * r[i - 1]
  }
}
acf(r_white, main = "ACF of whitened bam_NEE residuals") # decorrelated residuals: within the range

summary(bam_NEE)


