# =============================================================================
# GAM models: Management effects on NPP, respiration (chamber survey data)

library(mgcv)
library(dplyr)
library(gratia)
library(DHARMa)
library(lubridate)
library(emmeans)
library(patchwork)

# =============================================================================
# upload data: use "datacleaning_chambers.r" if not yet loaded
# =============================================================================
# source(here("r scripts", "datacleaning_chambers.r"))

# make sure management is in factor form
chamber_daily_mngmnt$management <- factor(as.character(chamber_daily_mngmnt$management))
# same for crop stage:
chamber_daily_mngmnt$crop_stage_simple <- factor(chamber_daily_mngmnt$crop_stage_simple)
# set reference level crop stage as "mature":
chamber_daily_mngmnt$crop_stage_simple <- relevel(chamber_daily_mngmnt$crop_stage_simple, ref = "mature")
# replace any NaNs with real NAs:
# Replace NaN with NA
chamber_daily_mngmnt <- chamber_daily_mngmnt %>%
  mutate(across(where(is.numeric), ~ na_if(., NaN)))
# calculate DOY:
chamber_daily_mngmnt$doy <- yday(chamber_daily_mngmnt$date)

# The models for chamber data will be much simpler, as there are fewer observations
# and thus fewer degrees of freedom to use up with complex interactions.

# -----------------------------------------------------------------------------
# Model development: contemplating crop stage, and how (or if) to incorporate
# -----------------------------------------------------------------------------

# summary with crop stage simple reveals possible issues with singularities --
# both levels of management need to appear for each level of crop_stage_simple
# for the models to be symmetrical.

table(chamber_daily_mngmnt$management, chamber_daily_mngmnt$crop_stage_simple)
#                 mature dormant early fallow grain_fill reproductive vegetative
# conventional      5       0     7      8          9            4         12
# organic          14       2     1      0          9            6         21

# ok so fallow is a 0 for organic in this data, and dormant is a 0 in conventional
# which means it never even makes it to the model comparison
# there is a zero in the *interaction* between organic and vegetative, because organic has 
# 21 veg obs versus conventional's 12...but conventional has 8 obs in fallow. 
# A model including crop stage will try to 'difference' against the 'mature' reference level.

# for this dataset, which is much smaller than the EC data of average fluxes daily,
# there is asymmetry in the crop stage distributions, so an interaction is not supported.

# subbing in NDVI is not terribly useful with ground-level observations, as NDVI
# is directly correlated with NEE as measured by hand; inclusion as a proxy for 
# crop stage may instead absorb actually interesting trends in NEE.

# Because of these covarying factors, and the relatively low sample size compared
# to the EC tower data, not including a crop stage into the model and will instead
# rely on day-of-year to encompass seasonality.
# it's not useful to include crop stage or NDVI as a proxy for these models, given
# the relatively low sample size. What we're really interested in is the effect
# of management on NEE and respiration, period, while accounting for seasonal 
# variance and single-year occurrence of drought (2020).

# -----------------------------------------------------------------------------
# Modeling: respiration, considering inter-year variation, drought, cropping
# -----------------------------------------------------------------------------

gam_Rs_simple <- gam(
  Rs_umolm2sec_mean ~
    management +
    s(doy, k = 25) +
    s(airtempC_mean, by = management, k = 12) +
    s(soilm_perc_mean, by = management, k = 11),
  data = chamber_daily_mngmnt, method = "REML")

gam.check(gam_Rs_simple)
simulateResiduals(gam_Rs_simple, plot = TRUE) # fit issues:

# initial fit of Gaussian GAM (first model attempts here) revealed heteroscedasticity in
# residuals, as well as a flipped sign for the organic effect (negative, indicating lower
# respiration than conventional on average). Re-fitting with a Gamma family and log link resolved
# the residuals diagnostics issue, and produced the positive (and biologically plausible/
# meaningful) estimate of organic management on Reco. this is hopefully to capture 
# the near-exponentially positive relationship between
# Rs and temperature, which may be causing the misfit in distribution in the 
# Gaussian version of events, which resulted in heteroskedasticity.

# The previous model-building has revealed some 
# intricacies that we need to deal with: a gamma family for the Rs model,
# disuse of the tensor interaction between soil moisture and air temp because
# they are less tightly related than air temp and VPD (as in EC models), and
# dropping of air temp from the model because it is encompassed in DOY (versus
# soil moisture, which can vary with soil biology, which may also be dependent
# on management and thus needs to stay in the model.

# remove two negative values bc gamma family cannot deal with negatives (and it's
# biologically realistic; they're close to zero anyway so likely within the range
# of noise around zero)
dat <- chamber_daily_mngmnt %>% filter(Rs_umolm2sec_mean > 0)

# 1. Baseline: management only
gam_Rs_0 <- gam(Rs_umolm2sec_mean ~ management,
                family = Gamma(link = "log"),
                data = dat, method = "ML")

# 2. Add seasonal trajectory
gam_Rs_1 <- gam(Rs_umolm2sec_mean ~ management + s(doy, k = 25),
                family = Gamma(link = "log"),
                data = dat, method = "ML")

# 3. Add soil moisture -- shared smooth
gam_Rs_2 <- gam(Rs_umolm2sec_mean ~ management + s(doy, k = 25) +
                  s(soilm_perc_mean, k = 8),
                family = Gamma(link = "log"),
                data = dat, method = "ML")

# 4. Management-varying soil moisture -- your current final
gam_Rs_final <- gam(Rs_umolm2sec_mean ~ management + s(doy, k = 25) +
                      s(soilm_perc_mean, by = management, k = 8),
                    family = Gamma(link = "log"),
                    data = dat, method = "ML")

# 5. Add drought resilience test by having management interact with year id
gam_Rs_drought <- gam(Rs_umolm2sec_mean ~ management * factor(year) +
                        s(doy, k = 25) +
                        s(soilm_perc_mean, by = management, k = 8),
                      family = Gamma(link = "log"),
                      data = dat, method = "ML")

AIC(gam_Rs_0, gam_Rs_1, gam_Rs_2, gam_Rs_final, gam_Rs_drought)
#                    df      AIC
# gam_Rs_0        3.000000 584.8685
# gam_Rs_1        8.992322 499.6187
# gam_Rs_2        9.869282 500.9810
# gam_Rs_final   10.690155 498.3722
# gam_Rs_drought 13.421699 487.8272

# drought has the best model fit: refit with REML, evaluate with DHARMa
gam_Rs_drought <- update(gam_Rs_drought, REML = TRUE)
simulateResiduals(gam_Rs_drought, plot = TRUE) # model diagnostics look good

# summarize:
summary(gam_Rs_drought)
# organic itself does not have a significant effect on Rs compared to conventional,
# but it's interaction with year does (p = 0.01) as does yaer itself (p < 0.0001)

# residuals vs. fitted colored by year -- to identify systematic discrepancies
# 2019 in red, 2020 (drought year) in blue
plot(fitted(gam_Rs_drought), residuals(gam_Rs_drought),
     col = ifelse(dat$year == 2019, "red", "blue"),
     pch = 16); abline(h = 0, lty = 2)
# residuals are well mixed across fitted range despite year, indicating that the
# management * year interaction accounts for differences in year (drought in 2020)
# which was a better fit in AIC comparison

# ACF -- 
acf(residuals(gam_Rs_drought))
# unfortunately there is a lag-1 autocorrelation that is concerning: need to figure
# out how to deal with this because biweekly measurements may not be independent from
# each other due to the biological 'memory' of soil (slow moving change in response
# to things like substrate availability and ambient climate)

# however, there aren't that many observations so using a corAR1 adjustment in gamm
# will be too clunky. try using a gls() model with a correlation structure, applied
# to a parametric version of the model.

#library(nlme)

# fit with AR(1) correlation structure, reset by management/year block
gls_Rs_drought <- gls(
  Rs_umolm2sec_mean ~ management * factor(year) + doy + soilm_perc_mean,
  correlation = corAR1(form = ~ 1 | management/year),
  data = chamber_daily_mngmnt %>% 
    arrange(management, year, doy),
  method = "REML"
)

summary(gls_Rs_drought)

# gamm with AR1 -- will estimate phi but with limited precision given block size
gamm_Rs_drought <- gamm(
  Rs_umolm2sec_mean ~ management * factor(year) +
    s(doy, k = 25) +
    s(soilm_perc_mean, by = management, k = 8),
  family = Gamma(link = "log"),
  correlation = corAR1(form = ~ 1 | management/year),
  data = dat %>% arrange(management, year, doy),
  method = "REML"
)

# check phi
coef(gamm_Rs_drought$lme$modelStruct$corStruct, unconstrained = FALSE)
# phi = 0.7841415

# check whether ACF is absorbed
acf(residuals(gamm_Rs_drought$lme, type = "normalized")) # successful absorption
summary(gamm_Rs_drought$gam)
# parametric interaction:             Estimate Std. Error t value Pr(>|t|)    
# managementorganic:factor(year)2020  0.58883    0.58996   0.998    0.321    
# not significant, but indicates same magnitude of effect as in uncorrected model

# results:
# the uncorrected model (gam_Rs_drought) shows significant, ecologically reasonable and
# coherent patterns
# correcting for autocorrelation results in loss of significance, but preservation of
# effect direction and magnitude
# the correction is probably overcorrecting, due to short block lengths (n = 12
# observations per block (block = management * year))
# phi = 0.78 for the autocorrelation corrected gamm model, which is very high
# and indicates that the AR(1) estimate is overshooting bc the time series is too
# short to accommodate all that complexity, resulting in inflated phi.
# 

# USE THE UNCORRECTED MODEL and state clearly: "We observed residual autocorrelation
# in the management*year blocks (0.35-0.45 ACF lag-1), indicating that SE's from the 
# final GAM were optimistic. an Ar(1)-corrected model (gamm with corAR1) yielded
# directionally-consistent estimates (interaction beta = 0.59), but largely attenuated
# significance (p = 0.321), likely reflecting phi overestimation given limited
# within-block sample size (n = 12 per block). Given this uncertainty, chamber-based
# Rs results are interpreted as directionally supportive of the EC analysis rather 
# than independently conclusive. Primary model results here are uncorrected for temporal
# autocorrelation for model parsimony, with the above caveat."

# Additionally, the normalized residuals showed no significant autocorrelation after correction
# (e.g., the gamm normalized model), with a phi = 0.78 that nonetheless whitened the
# residuals. This indicates that even tho phi is pretty high and thus imprecise, the model
# is internally consistent with the uncorrected one, validating the model structure even if point
# estimates are uncertain.

summary(gam_Rs_drought)
# Parametric coefficients:
#                                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          1.7707     0.1096  16.159  < 2e-16 ***
# managementorganic                    0.1670     0.1571   1.063   0.2907    
# factor(year)2020                    -0.7012     0.1673  -4.190 6.82e-05 ***
# managementorganic:factor(year)2020   0.5811     0.2274   2.555   0.0124 *   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#                                               edf Ref.df      F p-value    
#   s(doy)                                    5.375  6.706 11.825 < 2e-16 ***
#   s(soilm_perc_mean):managementconventional 1.000  1.000 11.321 0.00115 ** 
#   s(soilm_perc_mean):managementorganic      1.000  1.000  0.095 0.75866    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.626   Deviance explained =   68%
# -ML = 241.25  Scale est. = 0.27137   n = 96

# extract the outputs from the drought model (uncorrected, with notes as above)
summary(gam_Rs_drought)$p.table
#                                     Estimate Std. Error   t value     Pr(>|t|)
# (Intercept)                         1.7706918  0.1095818 16.158629 1.359556e-27
# managementorganic                   0.1670246  0.1570924  1.063225 2.907045e-01
# factor(year)2020                   -0.7011688  0.1673466 -4.189918 6.819819e-05
# managementorganic:factor(year)2020  0.5810679  0.2274230  2.555010 1.240733e-02

summary(gam_Rs_drought)$s.table
#                                               edf   Ref.df           F     p-value
# s(doy)                                    5.374720 6.705837 11.82474219 0.000000000
# s(soilm_perc_mean):managementconventional 1.000032 1.000062 11.32073302 0.001151907
# s(soilm_perc_mean):managementorganic      1.000006 1.000012  0.09501121 0.758660454

# according to the above, the year 2020 resulted in lower respiration rates (likely
# in response to drought) but that the interaction with organic management
# attenuated this effect of drought

# -----------------------------------------------------------------------------
# Modeling: NEE, considering inter-year variation, drought, cropping
# -----------------------------------------------------------------------------

# going to start with a similar approach as for respiration, though I suspect
# there will be less of a risk of temporal autocorrelation given that NEE
# can change a bit more on-a-dime than respiration can in response to changing
# climate conditions (see notes on EC version of this analysis).

# so, a stepwise walk through model complexity to keep it as parsimonious as
# possible while also integrating important variables like the year*climate interaction
# (which is critical here bc NEE is the primary carbon balance metric in question
# for resilience in this study)

# 1. Baseline
gam_NEE_0 <- gam(NEE_umolm2sec_mean ~ management,
                 data = chamber_daily_mngmnt, method = "ML")

# 2. Seasonal trajectory
gam_NEE_1 <- gam(NEE_umolm2sec_mean ~ management + s(doy, k = 25),
                 data = chamber_daily_mngmnt, method = "ML")

# 3. Soil moisture shared
gam_NEE_2 <- gam(NEE_umolm2sec_mean ~ management + s(doy, k = 25) +
                   s(soilm_perc_mean, k = 8),
                 data = chamber_daily_mngmnt, method = "ML")

# 4. Management-varying soil moisture
gam_NEE_3 <- gam(NEE_umolm2sec_mean ~ management + s(doy, k = 25) +
                   s(soilm_perc_mean, by = management, k = 8),
                 data = chamber_daily_mngmnt, method = "ML")

# 5. Drought resilience -- key hypothesis test
gam_NEE_drought <- gam(NEE_umolm2sec_mean ~ management * factor(year) +
                         s(doy, k = 25) +
                         s(soilm_perc_mean, by = management, k = 8),
                       data = chamber_daily_mngmnt, method = "ML")

AIC(gam_NEE_0, gam_NEE_1, gam_NEE_2, gam_NEE_3, gam_NEE_drought)
#                     df      AIC
# gam_NEE_0        3.000000 596.7869
# gam_NEE_1        8.028393 580.7682
# gam_NEE_2        8.926209 582.7505
# gam_NEE_3       10.083677 583.2977
# gam_NEE_drought 11.875568 581.1716
# indicates good fit for the drought model, and includes the ecologically-
# needed variable of soil moisture (where model 1 doesn't). Fits better than just 
# management-varying soil moisture alone.

# update to REML fit
gam_NEE_drought <- update(gam_NEE_drought, REML = TRUE)
simulateResiduals(gam_NEE_drought, plot = TRUE) # TONS of issues with residual/predicted
# quantile deviations detected and combined adjusted quantile test significant
# inverted u-shape means variance is highest at intermediate predicted values,
# and compressed at high and low predicted values.

# summarize:
summary(gam_NEE_drought)

# residuals vs. fitted colored by year -- to identify systematic discrepancies
# 2019 in red, 2020 (drought year) in blue
plot(fitted(gam_NEE_drought), residuals(gam_NEE_drought),
     col = ifelse(dat$year == 2019, "red", "blue"),
     pch = 16); abline(h = 0, lty = 2)
# residuals are well mixed across fitted range despite year, indicating that the
# management * year interaction accounts for differences in year (drought in 2020)

# ACF -- potentially some temporal autocorrelation to deal with
acf(residuals(gam_NEE_drought))

# troubleshooting:
# check ACF separately per management x year block
chamber_daily_mngmnt %>%
  group_by(management, year) %>%
  arrange(doy, .by_group = TRUE) %>%
  group_map(~ {
    r <- residuals(gam_NEE_drought)[as.numeric(rownames(.x))]
    acf(r, main = paste(.y$management, .y$year))
  })
# because all the within-block variance does NOT have evidence of autocorrelation,
# we do not need to account for it in the model. the elevated lag-1 in the acf of the
# model itself was likely a dataframe ordering artifact, and the blocks are internally
# temporally independent. TL;DR -- NEE 'memory' between visits is minimal or nonexistent.

# the heteroscedasticity in residuals for the uncorrected final model 
# may be driven by the fact that NEE variance scales with flux magnitude in 
# both directions. Try including observation-level weights using SE column:
gam_NEE_drought_wtd <- gam(
  NEE_umolm2sec_mean ~
    management * factor(year) +
    s(doy, k = 25) +
    s(soilm_perc_mean, by = management, k = 8),
  weights = 1 / (NEE_daily_se^2), # observation-level weighting, as SE captures real 
  # measurement uncertainty that scales with flux magnitude as exhibited by orig.
  # model's heteroscedasticity; here, we are down-weighting high-variance observations
  family = gaussian,
  data = chamber_daily_mngmnt,
  method = "REML"
)

simulateResiduals(gam_NEE_drought_wtd, plot = TRUE) # no weird heteroscedasticity,
# minor overdispersion but not major.

summary(gam_NEE_drought_wtd)
# Parametric coefficients:
#                                       Estimate Std. Error t value Pr(>|t|)  
#   (Intercept)                          0.2763     0.3715   0.744   0.4592  
#   managementorganic                   -1.9186     0.8004  -2.397   0.0188 *
#   factor(year)2020                    -0.5121     0.4247  -1.206   0.2313  
#   managementorganic:factor(year)2020   0.1708     1.0810   0.158   0.8749  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#                                               edf   Ref.df  F  p-value    
#   s(doy)                                    5.579  6.833 2.668 0.015032 *  
#   s(soilm_perc_mean):managementconventional 1.000  1.001 0.747 0.390020    
#   s(soilm_perc_mean):managementorganic      3.934  4.743 4.799 0.000894 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
#   R-sq.(adj) =  0.407   Deviance explained = 49.8%
#   -REML = 266.77  Scale est. = 23.382    n = 98

# resolves the residuals issue, and in doing so it fixes the underlying
# issue of collar-level variance in NEE that scaled with flux magnitude. Here,
# we're weighting each observation by the inverse of its squared daily SE to account
# for that pattern of variance in a Gaussian gam.

# takeaways:
# 1) organic results in significantly greater uptake (more negative NEE) by 
# almost 2umol/m2/sec when compared to conventional management

# 2) nonlinear relationship between organic management and soil moisture, indicating
# a greater, more complex relationship with moisture availability -- likely
# because of greater SOM water retention, more active root/microbial community, etc.

# 3) drought year * management is not significant

# extract the outputs from the drought model (uncorrected, with notes as above)
summary(gam_NEE_drought_wtd)$p.table
#                                     Estimate Std. Error     t value   Pr(>|t|)
# (Intercept)                         0.2762619  0.3715222  0.7435947 0.45920932
# managementorganic                  -1.9185649  0.8004339 -2.3969059 0.01876681
# factor(year)2020                   -0.5121292  0.4247424 -1.2057409 0.23132357
# managementorganic:factor(year)2020  0.1707521  1.0809858  0.1579596 0.87487006

summary(gam_NEE_drought_wtd)$s.table
#                                               edf   Ref.df         F      p-value
# s(doy)                                    5.578774 6.833272 2.6684455 0.0150317212
# s(soilm_perc_mean):managementconventional 1.000281 1.000528 0.7472386 0.3900199833
# s(soilm_perc_mean):managementorganic      3.934387 4.743116 4.7994285 0.0008943208

# -----------------------------------------------------------------------------
# Plotting survey data
# -----------------------------------------------------------------------------

chamber_daily_mngmnt <- chamber_daily_mngmnt %>%
  mutate(period = if_else(date < as.Date("2020-01-01"), "2019", "2020"))

# extract doy smooth predictions for each model
library(lubridate)

# create prediction grid -- one row per doy x management x year combination
pred_grid <- expand.grid(
  doy            = seq(min(chamber_daily_mngmnt$doy), 
                       max(chamber_daily_mngmnt$doy), by = 1),
  management     = c("conventional", "organic"),
  year           = c("2019", "2020"),
  soilm_perc_mean = median(chamber_daily_mngmnt$soilm_perc_mean, na.rm = TRUE)
) %>%
  mutate(
    # reconstruct a date for plotting on your existing x-axis
    date = as.POSIXct(strptime(paste(year, doy), "%Y %j")),
    year = factor(year)
  )

# NEE predictions
nee_pred <- predict(gam_NEE_drought_wtd, newdata = pred_grid, se.fit = TRUE)
pred_grid$NEE_fit    <- nee_pred$fit
pred_grid$NEE_se     <- nee_pred$se.fit
pred_grid$NEE_upper  <- nee_pred$fit + 1.96 * nee_pred$se.fit
pred_grid$NEE_lower  <- nee_pred$fit - 1.96 * nee_pred$se.fit

# Rs predictions (on response scale -- Gamma log link)
rs_pred <- predict(gam_Rs_drought, newdata = pred_grid, 
                   se.fit = TRUE, type = "response")
pred_grid$Rs_fit    <- rs_pred$fit
pred_grid$Rs_se     <- rs_pred$se.fit
pred_grid$Rs_upper  <- rs_pred$fit + 1.96 * rs_pred$se.fit
pred_grid$Rs_lower  <- rs_pred$fit - 1.96 * rs_pred$se.fit

# get observed doy range per year
doy_range <- chamber_daily_mngmnt %>%
  group_by(year) %>%
  summarise(doy_min = min(doy), doy_max = max(doy))

pred_grid <- pred_grid %>%
  left_join(doy_range, by = "year") %>%
  filter(doy >= doy_min & doy <= doy_max)

# Plot NEE over time for both fields
p1 <- ggplot() +
  geom_point(data = chamber_daily_mngmnt, 
             aes(x = date, y = NEE_umolm2sec_mean, color = management),
             alpha = 0.6, size = 1.5) +
  # GAM ribbon
  geom_ribbon(data = pred_grid,
              aes(x = date, ymin = NEE_lower, ymax = NEE_upper,
                  fill = management),
              alpha = 0.2) +
  # GAM fitted line
  geom_line(data = pred_grid,
            aes(x = date, y = NEE_fit, color = management),
            linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "darkred", linetype = "dashed")+
  labs(x = element_blank(),
    y = expression(paste("NEE (µmol", " m"^-2, " s"^-1, ")")),
    color = "management") +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00"))+
  theme_minimal() +
  theme(legend.position = "top")+
  scale_x_datetime(
    date_labels = "%b %d, %Y",  # Format: "Jan 01, 2023"
    date_breaks = "1 month"             # Set tick marks every month
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate text for readability
p1

# Plot soil respiration over time for both fields
p2 <- ggplot() +  # Convert NEE to respiration
  geom_point(data=chamber_daily_mngmnt, 
             aes(x = date, y = Rs_umolm2sec_mean, color = management),
             alpha = 0.6, size = 1.5) +
  geom_ribbon(data = pred_grid,
              aes(x = date, ymin = Rs_lower, ymax = Rs_upper,
                  fill = management),
              alpha = 0.2) +
  geom_line(data = pred_grid,
            aes(x = date, y = Rs_fit, color = management),
            linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "darkred", linetype = "dashed")+
  labs(
    x = element_blank(),
    y = expression(paste("Respiration (µmol", " m"^-2, " s"^-1, ")")),
    color = "management", fill = "management"
  ) +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00"))+
  theme_minimal() +
  theme(legend.position = "none")+
  scale_x_datetime(
    date_labels = "%b %d, %Y",  # Format: "Jan 01, 2023"
    date_breaks = "1 month"             # Set tick marks every month
  )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate text for readability
p2

# Combine plots
p1+p2+plot_layout(ncol = 1)
# p1 + p2 + p3 + plot_layout(ncol = 1)

#################################
# Soil moisture response curves
#################################

# extract smooth estimates for soil moisture terms
sm_NEE <- smooth_estimates(gam_NEE_drought_wtd) %>% 
  filter(grepl("soilm_perc_mean", .smooth))

sm_Rs <- smooth_estimates(gam_Rs_drought) %>%
  filter(grepl("soilm_perc_mean", .smooth))

# plot with ggplot, partial effect on y-axis, 
# ribbon for CI, rug for data density
p1 <- ggplot(sm_NEE, aes(x = soilm_perc_mean, y = .estimate, 
                   ymin = .estimate - 2*.se, ymax = .estimate + 2*.se,
                   fill = management, color = management)) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  geom_rug(data = chamber_daily_mngmnt, 
           aes(x = soilm_perc_mean), inherit.aes = FALSE,
           alpha = 0.3) +
  facet_wrap(~ management) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Partial effect on NEE", x = "Soil moisture (%)") +
  theme_bw() +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  theme(legend.position = "none")

p2 <- ggplot(sm_Rs, aes(x = soilm_perc_mean, y = .estimate, 
                   ymin = .estimate - 2*.se, ymax = .estimate + 2*.se,
                   fill = management, color = management)) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  geom_rug(data = chamber_daily_mngmnt, 
           aes(x = soilm_perc_mean), inherit.aes = FALSE,
           alpha = 0.3) +
  facet_wrap(~ management) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Partial effect on respiration", x = "Soil moisture (%)") +
  theme_bw() +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  theme(legend.position = "none")

p1+p2+plot_layout(ncol = 1)

##################################################################
# Management x Year marginal means (drought contrast)
##################################################################

# extract marginal means, holding smooths at median
emm_NEE <- emmeans(gam_NEE_drought_wtd, ~ management * factor(year),
                   at = list(soilm_perc_mean = median(chamber_daily_mngmnt$soilm_perc_mean,
                                                      na.rm = TRUE),
                             doy = 180),
                   data = chamber_daily_mngmnt)

emm_Rs <- emmeans(gam_Rs_drought, ~ management * factor(year),
                  at = list(soilm_perc_mean = median(chamber_daily_mngmnt$soilm_perc_mean,
                                                     na.rm = TRUE),
                            doy = 180),
                  data = chamber_daily_mngmnt,
                  type = "response") # important bc the model uses a log link

emm_NEE_df <- as.data.frame(emm_NEE)
emm_Rs_df  <- as.data.frame(emm_Rs)

# check column names -- emmeans names the CI columns differently 
# depending on family (asymmetric for Gamma)
names(emm_NEE_df)
names(emm_Rs_df)

# standardize column names for combining
emm_NEE_df <- emm_NEE_df %>%
  mutate(response_var = "NEE",
         estimate = emmean)

emm_Rs_df <- emm_Rs_df %>%
  mutate(response_var = "Ecosystem Respiration",
         estimate = response)

emm_combined <- bind_rows(emm_NEE_df, emm_Rs_df) %>%
  mutate(
    year = as.character(year),
    year_label = ifelse(year == "2019", "2019\n(soybean)", "2020\n(barley/drought)"),
    management = factor(management, levels = c("conventional", "organic"),
                        labels = c("Conventional", "Organic"))
  )

# plot

emm_combined$response_var <- as.factor(emm_combined$response_var)

# Update custom_labels to match
custom_labels <- c(
  'Ecosystem Respiration' = 'Respiration (\u03BCmol CO m\u207B\u00B2 s\u207B\u00B9)',
  'NEE' = 'NEE (\u03BCmol CO m\u207B\u00B2 s\u207B\u00B9)'
)

# Convert to factor
custom_labels <- as.factor(custom_labels)

p_emm <- ggplot(emm_combined, 
                aes(x = year_label, y = estimate, 
                    color = management, group = management)) +
  geom_point(size = 3, position = position_dodge(0.15)) +
  geom_line(position = position_dodge(0.15)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.1,
                position = position_dodge(0.15)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~ response_var,
             scales = "free_y",
             # labeller = labeller(response_var = c(
               # "Ecosystem Respiration" = expression(paste('Respiration ('*mu~ 'mol' ~CO~ m^-2~s^-1*')')),
               # "NEE" = expression(paste('NEE ('*mu~ 'mol' ~CO~ m^-2~s^-1*')'))
             # ))) +
             # labeller = label_parsed)+
             labeller = labeller(custom_labels))+
  scale_color_manual(values = c("Conventional" = "#E05C5C", 
                                "Organic"      = "#C49A00")) +
  labs(x = NULL, y = "Estimated marginal mean flux",
       color = "Management") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

p_emm

##################################################################
# TABLES
##################################################################

library(gtsummary)
library(broom)
library(flextable)
# set theme to compact: reduces row padding, font size
set_gtsummary_theme(theme_gtsummary_compact(set_theme = TRUE))

# table of results for Reco analysis:
ch_eco_tab <- 
  tbl_regression(
    gam_Rs_drought,
    exponentiate = FALSE,
    conf.int = TRUE) %>% 
  add_significance_stars(
    hide_ci = TRUE, hide_se = FALSE,
    hide_p = FALSE,
    pattern = "{p.value}{stars}"
  ) %>% 
  remove_row_type(type = "reference") %>% 
  italicize_levels() %>% bold_labels() %>% 
  modify_header(estimate~"**Estimate**") %>% 
  # update level names for the interactions
  modify_table_body(~ .x %>% 
                      dplyr::mutate(
                        label = dplyr::recode(label, 
                                              "factor(year)" = "year",
                                              "management * factor(year)" = "management * year",
                                              "s(soilm_perc_mean):managementconventional" = "s(% soil moisture * conventional)",
                                              "s(soilm_perc_mean):managementorganic" = "s(% soil moisture * organic)")
                      ))
# remove footnotes:
ch_eco_tab$table_styling$abbreviation <- 
  ch_eco_tab$table_styling$abbreviation %>% 
  dplyr::filter(column != c("conf.low", "std.error"))
# ch_eco_tab <- 
#   ch_eco_tab %>% 
#   modify_caption("**Predictors of respiration (chamber-based)**") %>% 
#   as_gt()

# NEE:
ch_nee_tab <- 
  tbl_regression(
    gam_NEE_drought_wtd,
    exponentiate = FALSE,
    conf.int = TRUE) %>% 
  add_significance_stars(
    hide_ci = TRUE, hide_se = FALSE,
    hide_p = FALSE,
    pattern = "{p.value}{stars}"
  ) %>% 
  remove_row_type(type = "reference") %>% 
  italicize_levels() %>% bold_labels() %>% 
  modify_header(estimate~"**Estimate**") %>% 
  # update level names for the interactions
  modify_table_body(~ .x %>% 
                      dplyr::mutate(
                        label = dplyr::recode(label, 
                                              "factor(year)" = "year",
                                              "management * factor(year)" = "management * year",
                                              "s(soilm_perc_mean):managementconventional" = "s(% soil moisture * conventional)",
                                              "s(soilm_perc_mean):managementorganic" = "s(% soil moisture * organic)")
                      ))
# remove footnotes:
ch_nee_tab$table_styling$abbreviation <- 
  ch_nee_tab$table_styling$abbreviation %>% 
  dplyr::filter(column != c("conf.low", "std.error"))
# ch_nee_tab <- 
  # ch_nee_tab %>% 
  # modify_caption("**Predictors of NEE (chamber-based)**") %>% 
  # as_gt()

# combine:
tbl_merge(
  tbls = list(ch_eco_tab, ch_nee_tab),
  tab_spanner = c(
    "**Respiration**", "**NEE**")
) %>% 
  bold_labels() %>% 
  modify_caption("**Predictors of fluxes (chamber-based)**") %>% 
  as_gt()
