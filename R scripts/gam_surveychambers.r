# =============================================================================
# GAM models: Management effects on NPP, respiration (chamber survey data)

library(mgcv)
library(dplyr)
library(gratia)
library(DHARMa)
library(lubridate)

# =============================================================================
# upload data: use "datacleaning_chambers.r" if not yet loaded
# =============================================================================
source(here("r scripts", "datacleaning_chambers.r"))

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

# I have chosen to include an interaction between management type and crop stage
# in the linear predictor variables as it is reasonable to expect that management
# type will result in different outcomes at each crop stage IRT metabolism and thus fluxes.

# The models for chamber data will be much simpler, as there are fewer observations
# and thus fewer degrees of freedom to use up with complex interactions.

# -----------------------------------------------------------------------------
# Model development: starting simple with a baseline model, no crop stage
# -----------------------------------------------------------------------------

# linear baseline -- are the climate vars significant at all by management?
lm_climate <- lm(NEE_umolm2sec_mean ~ management +
                   airtempC_mean + soilm_perc_mean,
                 data = chamber_daily_mngmnt)
summary(lm_climate)
# Estimate Std. Error t value Pr(>|t|)  
#   (Intercept)        4.92379    3.06231   1.608   0.1112  
#   managementorganic -1.80042    0.98511  -1.828   0.0708 .
#   airtempC_mean     -0.27203    0.10660  -2.552   0.0123 *
#   soilm_perc_mean    0.04330    0.08947   0.484   0.6295  

# there does seem to be an effect of management on NEE as measured by survey chambers,
# and as expected there is a significant effect of air temperature

# management-varying responses -- does each site respond differently 
# to the same temperature or moisture conditions?
gam_climate_mgmt <- gam(
  NEE_umolm2sec_mean ~
    management +
    s(doy, k = 25) + # day of year
    s(airtempC_mean, by = as.factor(management), k = 5)+
    s(soilm_perc_mean, by = as.factor(management), k = 6),
  data   = chamber_daily_mngmnt,
  method = "REML"
)
gam.check(gam_climate_mgmt) # adjusted k for soil moisture to k = 6 
summary(gam_climate_mgmt)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)        -0.2346     0.7256  -0.323   0.7472  
# managementorganic  -1.7278     0.9757  -1.771   0.0799 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#                                                         edf   Ref.df  F   p-value   
#   s(airtempC_mean):as.factor(management)conventional   1.000  1.000 0.109  0.7426   
#   s(airtempC_mean):as.factor(management)organic        1.551  1.922 5.205  0.0059 **
#   s(soilm_perc_mean):as.factor(management)conventional 1.000  1.000 0.474  0.4929   
#   s(soilm_perc_mean):as.factor(management)organic      1.000  1.001 0.069  0.7935 

# -----------------------------------------------------------------------------
# Model development: contemplating crop stage, and how (or if) to incorporate
# -----------------------------------------------------------------------------

# summary with crop stage simple reveals confusing lack of effect for organic and the vegetative phase, indicating
# possible issues with singularities -- both levels of management need to appear for each
# level of crop_stage_simple.
table(chamber_daily_mngmnt$management, chamber_daily_mngmnt$crop_stage_simple)
#               dormant early fallow grain_fill mature reproductive vegetative
# conventional       0     7      8          9      5            4         12
# organic            2     2      0          9     14            6         20

# ok so fallow is a 0 for organic in this data, and dormant is a 0 in conventional
# which means it never even makes it to the model comparison
# there is a zero in the *interaction* between organic and vegetative, because organic has 
# 20 veg obs versus conventional's 12...but conventional has 8 obs in fallow. 
# This model is trying to 'difference' against the 'mature' reference level.

# for this dataset, which is much smaller than the EC data of average fluxes daily,
# there is asymmetry in the crop stage distributions, so the interaction is not supported.
# I have 7 stages, 2 management levels, and a handful of empty or low values in their interactions.

# next steps, option 1: collapse the simple stages more aggressively so their are no
# zeros in the management * crop stage interaction.

chamber_daily_mngmnt <- chamber_daily_mngmnt %>%
  mutate(crop_stage_3 = case_when(
    crop_stage_simple %in% c("fallow", "early", "dormant")           ~ "low_canopy",
    crop_stage_simple %in% c("vegetative", "reproductive")           ~ "active_growth",
    crop_stage_simple %in% c("grain_fill", "mature")                 ~ "late_season"
  ))

table(chamber_daily_mngmnt$management, chamber_daily_mngmnt$crop_stage_3)

# add this adjusted crop stage:
gam_climate_mgmt_stage <- gam(
  NEE_umolm2sec_mean ~
    management + 
    crop_stage_3 +
    s(airtempC_mean, by = as.factor(management), k = 5) +
    s(soilm_perc_mean, by = as.factor(management), k = 6),
  data   = chamber_daily_mngmnt,
  method = "REML"
)
summary(gam_climate_mgmt_stage)
# Estimate Std. Error t value Pr(>|t|)
# (Intercept)              -0.6556     0.9933  -0.660    0.511
# managementorganic        -1.1297     1.0633  -1.062    0.291
# crop_stage_3late_season  -0.6269     1.0612  -0.591    0.556
# crop_stage_3low_canopy    1.5193     1.3866   1.096    0.276
# 
# Approximate significance of smooth terms:
#                                                         edf Ref.df     F p-value   
#   s(airtempC_mean):as.factor(management)conventional   1.00  1.000 0.201 0.65476   
#   s(airtempC_mean):as.factor(management)organic        1.54  1.906 4.839 0.00829 **
#   s(soilm_perc_mean):as.factor(management)conventional 1.00  1.000 0.053 0.81911   
#   s(soilm_perc_mean):as.factor(management)organic      1.00  1.000 0.116 0.73495   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.149   Deviance explained = 21.6%
# -REML = 279.41  Scale est. = 21.296    n = 98

# the model plus crop stage doesn't add much to explaining the variance; the R2 value
# for this model remains similar.
AIC(gam_climate_mgmt, gam_climate_mgmt_stage)
#                           df      AIC
# gam_climate_mgmt       7.922997 587.1128
# gam_climate_mgmt_stage 9.906927 588.7236
# confirmed: let's not worry about crop stage at the moment, for these supporting models.
# there are not enough observations to be splitting them up into crop stages with 
# low explanatory power for each one.

# CONCLUSION: NOT INCLUDING CROP STAGE IN SURVEY CHAMBER MODELS.

# -----------------------------------------------------------------------------
# Model development: necessary to include NDVI? and how to deal with co-variance?
# -----------------------------------------------------------------------------

# Model with doy only
gam_doy <- gam(NEE_umolm2sec_mean ~ 
                 management + 
                 s(doy, k = 25), 
               data = chamber_daily_mngmnt, 
               method = "REML")

# Model with climate variables only
gam_climate <- gam(NEE_umolm2sec_mean ~
                     management +
                     s(airtempC_mean, k = 12) +
                     s(soilm_perc_mean, k = 11),
                   data = chamber_daily_mngmnt, method = "REML")

# Model with doy and climate variables
gam_both <- gam(NEE_umolm2sec_mean ~ 
                  management +
                  s(doy, k = 25) +
                  s(airtempC_mean, k = 12) +
                  s(soilm_perc_mean, k = 11),
                data = chamber_daily_mngmnt, method = "REML")

# Compare AIC
AIC(gam_doy, gam_climate, gam_both)
#                 df      AIC
# gam_doy      8.764667 580.1626
# gam_climate  5.972523 587.5808
# gam_both    11.700213 584.0419

# continuous NDVI rather than crop stage to indicate where plants are at in terms of phenology:
# including in the model will essentially cover up any effects of drivers, bc
# NDVI tends to covary with NEE (one is measure of greenness; one is measure of the CO2 uptake from that green)

# final investigation: tensor interaction between climate variables

gam_tensor <- gam(NEE_umolm2sec_mean ~
                    management + 
                    s(doy, k = 25) + 
                    te(airtempC_mean, soilm_perc_mean, by = management, k = c(8, 8)),
                  data = chamber_daily_mngmnt, method = "REML")
AIC(gam_doy, gam_climate, gam_both, gam_tensor)
# gam_doy      8.764667 580.1626
# gam_climate  5.972523 587.5808
# gam_both    11.700213 584.0419
# gam_tensor  15.513535 581.1968

# tensor interaction results in best fit.
summary(gam_tensor)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)        0.02201    0.72977   0.030   0.9760  
# managementorganic -1.93515    0.97136  -1.992   0.0496 *
# 
# Approximate significance of smooth terms:
#                                                               edf Ref.df     F p-value  
#   s(doy)                                                   5.162  6.511 2.284  0.0395 *
#   te(airtempC_mean,soilm_perc_mean):managementconventional 3.000  3.001 1.607  0.1938  
#   te(airtempC_mean,soilm_perc_mean):managementorganic      3.001  3.002 0.351  0.7888  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.259   Deviance explained = 35.2%
# -REML = 263.12  Scale est. = 18.548    n = 98
simulateResiduals(gam_tensor, plot = TRUE) # residuals look good, no significant issues

# split the tensor interaction and remove the air temp variable (which covaries with DOY)

# Model with doy and soil moisture only:
gam_final <- gam(NEE_umolm2sec_mean ~ 
                  management +
                  s(doy, k = 25) +
                  s(soilm_perc_mean, by = management, k = 11),
                data = chamber_daily_mngmnt, method = "REML")
gam.check(gam_final)
simulateResiduals(gam_final, plot = TRUE)

AIC(gam_tensor, gam_final)



# -----------------------------------------------------------------------------
# Modeling: respiration
# -----------------------------------------------------------------------------

# start with similar model structure as for NEE, and again not using NDVI to avoid covariance:
gam_Rs_simple <- gam(
  Rs_umolm2sec_mean ~
    management +
    s(doy, k = 25) +
    s(airtempC_mean, by = management, k = 12) +
    s(soilm_perc_mean, by = management, k = 11),
  data = chamber_daily_mngmnt, method = "REML")

# add tensor interaction:
gam_Rs_tensor <- gam(
  Rs_umolm2sec_mean ~
    management + 
    s(doy, k = 25) +
    te(airtempC_mean, soilm_perc_mean, by = management, k = c(8, 8)),
  data = chamber_daily_mngmnt, method = "REML")

AIC(gam_Rs_simple, gam_Rs_tensor)
#                   df      AIC
# gam_Rs_simple 14.96618 619.2105
# gam_Rs_tensor 19.75492 624.7611

gam.check(gam_Rs_simple)
gam.check(gam_Rs_tensor) # both have normally fitted residuals
simulateResiduals(gam_Rs_simple, plot = TRUE) # both models have some fit issues, however
simulateResiduals(gam_Rs_tensor, plot = TRUE) # both models have some fit issues, however

# the residual checks indicate some issues, however.  I don't think it's an issue of 
# needing to integrate rolling averages for climate data (as for EC analysis) since
# these data are not continuous or even daily.

# re-run with gamma family, with log link
# this is hopefully to capture the near-exponentially positive relationship between
# Rs and temperature, which may be causing the misfit in distribution
# remove two negatives:
dat <- chamber_daily_mngmnt %>% filter(Rs_umolm2sec_mean > 0)
gam_Rs_gamma <- gam(
  Rs_umolm2sec_mean ~
    management + 
    s(doy, k = 25) +
    te(airtempC_mean, soilm_perc_mean, by = management, k = c(8, 8)),
  family = Gamma(link = "log"),
  data = dat, method = "REML")
# check model
gam.check(gam_Rs_gamma) # k is adequate for all parameters
simulateResiduals(gam_Rs_gamma, plot = TRUE)
# model outputs:
summary(gam_Rs_gamma)

gam_Rs_gamma_alt <- gam(
  Rs_umolm2sec_mean ~
    management + 
    s(doy, k = 25) +
    s(airtempC_mean, by = management, k = 8) +
    s(soilm_perc_mean, by = management, k = 8),
  family = Gamma(link = "log"),
  data = dat, method = "REML")
gam.check(gam_Rs_gamma_alt)
simulateResiduals(gam_Rs_gamma_alt, plot = TRUE)
summary(gam_Rs_gamma_alt)

AIC(gam_Rs_gamma, gam_Rs_gamma_alt)
# df      AIC
# gam_Rs_gamma     17.33596 507.3735
# gam_Rs_gamma_alt 14.83139 503.8567

# The best fit is separate smooths of air temp and moisture

# coefficient for organic parameter is +0.416, which is equivalent to 1.52 in non-log scale
# indicates that organic respires 1.5x more than conventional, on average, which tracks
# with the EC data analysis and makes sense if we consider biology: more SOM, more 
# biological activity, etc.
# initial fit of Gaussian GAM (first two model attempts here) revealed heteroscedasticity in
# residuals, as well as a flipped sign for the organic effect (negative, indicating lower
# respiration than conventional on average). Re-fitting with a Gamma family and log link resolved
# the residuals diagnostics issue, and produced the positive (and biologically plausible/
# meaningful) estimate of organic management on Reco.
# e.g.: "Organic management was associated with 52% greater ecosystem respiration (β = 0.416 on 
# log scale, exp(β) = 1.52, 95% CI on response scale [exp(0.416 − 1.96×0.124), 
# exp(0.416 + 1.96×0.124)] = [1.19, 1.93], p = 0.001)."

# TAKEAWAY:
# NEE: organic = −1.94 µmol m⁻² s⁻¹ vs conventional → organic takes up more CO2
# R_eco: organic = exp(0.416) = 1.52× conventional → organic respires more CO2
# GPP (= R_eco − NEE): organic must therefore have substantially higher GPP

# -----------------------------------------------------------------------------
# Modeling: respiration, considering inter-year variation, drought, cropping
# -----------------------------------------------------------------------------

# it's not useful to include crop stage or NDVI as a proxy for these models, given
# the relatively low sample size. What we're really interested in is the effect
# of management on NEE and respiration, period, while accounting for seasonal 
# variance and drought (2020). The previous model-building has revealed some 
# intricacies that we need to deal with: a gamma family for the Rs model,
# disuse of the tensor interaction between soil moisture and air temp because
# they are less tightly related than air temp and VPD (as in EC models), and
# dropping of air temp from the model because it is encompassed in DOY (versus
# soil moisture, which can vary with soil biology, which may also be dependent
# on management and thus needs to stay in the model.

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

# residuals vs. fitted colored by year -- to identify systematic discrepancies
# 2019 in red, 2020 (drought year) in blue
plot(fitted(gam_Rs_drought), residuals(gam_Rs_drought),
     col = ifelse(dat$year == 2019, "red", "blue"),
     pch = 16); abline(h = 0, lty = 2)

# ACF -- confirm no autocorrelation correction needed
acf(residuals(gam_Rs_drought))

# -----------------------------------------------------------------------------
# Plotting survey data
# -----------------------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(dplyr)

chamber_daily_mngmnt <- chamber_daily_mngmnt %>%
  mutate(period = if_else(date < as.Date("2020-01-01"), "2019", "2020"))

# Plot NEE over time for both fields
p1 <- ggplot(chamber_daily_mngmnt, aes(x = date, y = NEE_umolm2sec_mean, color = management)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(
    aes(group = interaction(management, period)),
    method = "loess",
    se = TRUE,
    span = 0.5
  ) +  geom_hline(yintercept = 0, color = "darkred", linetype = "dashed")+
  labs(
    x = element_blank(),
    y = expression(paste("NEE (µmol CO", " m"^-2, " s"^-1, ")")),
    color = "Management"
  ) +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  theme_minimal() +
  theme(legend.position = "top")+
  scale_x_datetime(
    date_labels = "%b %d, %Y",  # Format: "Jan 01, 2023"
    date_breaks = "1 month"             # Set tick marks every month
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate text for readability

# Plot soil respiration over time for both fields
p2 <- ggplot(chamber_daily_mngmnt, aes(x = date, y = Rs_umolm2sec_mean, color = management)) +  # Convert NEE to respiration
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(
    aes(group = interaction(management, period)),
    method = "loess",
    se = TRUE,
    span = 0.5
  ) +  geom_hline(yintercept = 0, color = "darkred", linetype = "dashed")+
  labs(
    x = element_blank(),
    y = expression(paste("Respiration (µmol CO", " m"^-2, " s"^-1, ")")),
    color = "Management"
  ) +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  theme_minimal() +
  theme(legend.position = "none")+
  scale_x_datetime(
    date_labels = "%b %d, %Y",  # Format: "Jan 01, 2023"
    date_breaks = "1 month"             # Set tick marks every month
  )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate text for readability

# Calculate and plot GPP over time for both fields
p3 <- ggplot(chamber_daily_mngmnt, aes(x = date, y = Rs_umolm2sec_mean-NEE_umolm2sec_mean, color = management)) +  # Convert NEE to respiration
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(
    aes(group = interaction(management, period)),
    method = "loess",
    se = TRUE,
    span = 0.5
  ) +  geom_hline(yintercept = 0, color = "darkred", linetype = "dashed")+
  labs(
    x = element_blank(),
    y = expression(paste("GPP (resp. - NEE) (µmol CO", " m"^-2, " s"^-1, ")")),
    color = "Management"
  ) +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  theme_minimal() +
  theme(legend.position = "none")+
  scale_x_datetime(
    date_labels = "%b %d, %Y",  # Format: "Jan 01, 2023"
    date_breaks = "1 month"             # Set tick marks every month
  )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate text for readability

# Combine plots
p1 + p2 + p3 + plot_layout(ncol = 1)

# GPP by soil moisture:
ggplot(chamber_daily_mngmnt, aes(x = Rs_umolm2sec_mean-NEE_umolm2sec_mean, 
                                 y = soilm_perc_mean, 
                                 color = management)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE) +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  labs(
    x = "Soil Moisture, VWC (%)",
    y = expression(paste("GPP (µmol CO", " m"^-2, " s"^-1, ")")),
    color = "Management"
  ) +
  theme_minimal()+
  facet_wrap(~year)

# against tensor interaction:
grid <- expand.grid(
  airtempC_mean = seq(min(chamber_daily_mngmnt$airtempC_mean, na.rm=TRUE),
                      max(chamber_daily_mngmnt$airtempC_mean, na.rm=TRUE), length = 60),
  soilm_perc_mean = seq(min(chamber_daily_mngmnt$soilm_perc_mean, na.rm=TRUE),
                        max(chamber_daily_mngmnt$soilm_perc_mean, na.rm=TRUE), length = 60),
  doy = median(chamber_daily_mngmnt$doy, na.rm=TRUE)
)

pred_conv <- predict(gam_tensor,
                     newdata = mutate(grid, management = "conventional"), se.fit = TRUE)
pred_org  <- predict(gam_tensor,
                     newdata = mutate(grid, management = "organic"), se.fit = TRUE)

grid$diff   <- pred_org$fit - pred_conv$fit
grid$diff_se <- sqrt(pred_org$se.fit^2 + pred_conv$se.fit^2)
grid$org <- pred_org$fit
grid$conv <- pred_conv$fit

common_limits <- range(c(grid$org, grid$conv), na.rm = TRUE)
common_limits <- range(grid$diff, na.rm = TRUE)

scale_common_NEE <- scale_fill_gradient2(
  low = "darkgreen", mid = "white", high = "brown",
  midpoint = 0,
  limits = common_limits,
  name = expression(paste(Delta * NEE, " (", mu, "mol/m²/s", ")"))
)

p1 <- ggplot(grid, aes(airtempC_mean, soilm_perc_mean, fill = conv)) +
       geom_raster() +
      scale_common_NEE +
       geom_contour(aes(z = conv), color = "white", alpha = 0.4) +
       # scale_fill_gradient2(low = "brown", mid = "white", high = "darkgreen", midpoint = 0,
                            # name = "Conventional\nNEE (µmol/m²/s¹)") +
       labs(x = "Air temperature (°C)", y = "Soil moisture (% VWC)") +
       theme_minimal()+ geom_point(data = chamber_daily_mngmnt,
                                   aes(airtempC_mean, soilm_perc_mean),
                                   inherit.aes = FALSE, size = 0.5, alpha = 0.3) +
  annotate("label", x = 8, y = 3, label = "cool & dry\n(dormant)",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75)) +
  annotate("label", x = 28, y = 3, label = "hot & dry\n(drought stress)",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75)) +
  annotate("label", x = 8, y = 25, label = "cool & wet",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75)) +
  annotate("label", x = 28, y = 25, label = "hot & wet",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75))+
  ggtitle("Conventional")



p2 <- ggplot(grid, aes(airtempC_mean, soilm_perc_mean, fill = org)) +
  geom_raster() +
  scale_common +
  geom_contour(aes(z = org), color = "white", alpha = 0.4) +
  # scale_fill_gradient2(low = "brown", mid = "white", high = "darkgreen", midpoint = 0,
                       # name = "Organic \nNEE (µmol/m²/s)") +
  labs(x = "Air temperature (°C)", y = "Soil moisture (% VWC)") +
  theme_minimal()+ geom_point(data = chamber_daily_mngmnt,
                              aes(airtempC_mean, soilm_perc_mean),
                              inherit.aes = FALSE, size = 0.5, alpha = 0.3) +
  annotate("label", x = 8, y = 3, label = "cool & dry\n(dormant)",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75)) +
  annotate("label", x = 28, y = 3, label = "hot & dry\n(drought stress)",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75)) +
  annotate("label", x = 8, y = 25, label = "cool & wet",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75)) +
  annotate("label", x = 28, y = 25, label = "hot & wet",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75))+
  ggtitle("Organic")

p2 + p1

# difference plot:
p3 <- ggplot(grid, aes(airtempC_mean, soilm_perc_mean, fill = grid$diff)) +
  geom_raster() +
  scale_common_NEE +
  geom_contour(aes(z = org), color = "white", alpha = 0.4) +
  labs(x = "Air temperature (°C)", y = "Soil moisture (% VWC)") +
  theme_minimal()+ 
    geom_point(data = chamber_daily_mngmnt,
                              aes(airtempC_mean, soilm_perc_mean),
                              inherit.aes = FALSE, size = 0.5, alpha = 0.3) +
  annotate("label", x = 8, y = 3, label = "cool & dry\n(dormant)",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75)) +
  annotate("label", x = 28, y = 3, label = "hot & dry\n(drought stress)",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75)) +
  annotate("label", x = 8, y = 25, label = "cool & wet",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75)) +
  annotate("label", x = 28, y = 25, label = "hot & wet",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75))+
  ggtitle("NEE: Organic - Conventional")
p3

# respiration version:

# against tensor interaction:
pred_conv <- predict(gam_Rs_gamma,
                     newdata = mutate(grid, management = "conventional"), se.fit = TRUE)
pred_org  <- predict(gam_Rs_gamma,
                     newdata = mutate(grid, management = "organic"), se.fit = TRUE)

grid$diff   <- pred_org$fit - pred_conv$fit
grid$diff_se <- sqrt(pred_org$se.fit^2 + pred_conv$se.fit^2)
grid$org <- pred_org$fit
grid$conv <- pred_conv$fit

common_limits <- range(c(grid$org, grid$conv), na.rm = TRUE)
common_limits <- range(grid$diff, na.rm = TRUE)

scale_common_ER <- scale_fill_gradient2(
  low = "lightblue", mid = "white", high = "navyblue",
  midpoint = 0,
  limits = common_limits,
  name = "soil respiration\n(µmol/m²/s)"
)

p4 <- ggplot(grid, aes(airtempC_mean, soilm_perc_mean, fill = grid$diff)) +
  geom_raster() +
  scale_common_ER +
  geom_contour(aes(z = org), color = "white", alpha = 0.4) +
  labs(x = "Air temperature (°C)", y = "Soil moisture (% VWC)") +
  theme_minimal()+ 
  geom_point(data = chamber_daily_mngmnt,
             aes(airtempC_mean, soilm_perc_mean),
             inherit.aes = FALSE, size = 0.5, alpha = 0.3) +
  annotate("label", x = 8, y = 3, label = "cool & dry\n(dormant)",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75)) +
  annotate("label", x = 28, y = 3, label = "hot & dry\n(drought stress)",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75)) +
  annotate("label", x = 8, y = 25, label = "cool & wet",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75)) +
  annotate("label", x = 28, y = 25, label = "hot & wet",
           color = "grey20", size = 2.8, label.size = 0,
           fill = alpha("white", 0.75))+
  ggtitle("Soil respiration: Organic - Conventional")
p4
