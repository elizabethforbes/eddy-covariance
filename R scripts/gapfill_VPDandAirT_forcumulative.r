# =============================================================================
# Cumulative Annual Carbon Flux Analysis: gap-filling climate data
#
# to compute annual and running-cumulative NEE, GPP, and Reco for EF01 (organic)
# and EF02 (conventional) using observed EC tower data plus gap-filled with fitted
# GAM models, we need continuous (i.e., gap-filled) air temperature and VPD data.
# Use ERA5 modeled air temperature and RH to calculate VPD for those dates that 
# VPD does not exist due to tower outages; compare empirical air temp/VPD to 
# those generated with ERA5 data.
#
# required: ec_daily_mngmnt (full dataset, all years, all crops)
#
# Analysis period: 2018–2021 (years with management + yield data)

# =============================================================================

library(mgcv)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(zoo)      # rollmean for Reco lagged climate predictors

# =============================================================================
# Data preparation
# =============================================================================

# Build full 2018–2021 dataset from the un-trimmed source.
# Unlike data1 (2018–2020, no alfalfa), data_cumul includes all crops and 2021,
# so that cumulative totals reflect actual field carbon exchange across the
# entire observation period.

df <- ec_daily_mngmnt3 %>%
  filter(year %in% 2018:2021) %>%
  mutate(
    management        = factor(as.character(management)),
    crop_stage_simple = relevel(factor(as.character(crop_stage_simple)), ref = "mature"),
    year_f            = factor(year)
  ) %>%
  arrange(management, year, doy) %>%
  group_by(management, year) %>%
  mutate(AR.start = row_number() == 1) %>%
  ungroup() %>% 
  select(management, date, year, month, doy, NEE, GPP, Reco, NEE_se, Reco_se, 
         # climate variables
         air_temperature, RH, VPD, RH_era5, airt_era5,
         # crop system variables
         current_crop, days_since_fertilizer, days_since_tillage, days_since_herbicide,
         days_since_planting, crop_stage, crop_stage_simple)

# calculate VPD and air temperature using the ERA5 values in the df to fill 
# gaps in the tower data (which there are for each row there is no flux data)

# temp_C (air temperature in °C), RH (relative humidity in %), 
# and observed_vpd (empirical VPD in kPa), measured by towers when online

# Define the VPD calculation functions, which uses saturation vapor pressure (es),
# and then actual vapor pressure (ea) based on relative humidity

# Use dplyr to calculate VPD and flag differences

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

# Tolerance thresholds
vpd_tolerance <- 0.05  # kPa
temp_tolerance <- 0.5  # °C

# air_temperature (empirical air temp °C)
# RH (empirical relative humidity %)
# VPD (empirical VPD in Pa)
# airt_era5 (ERA5 air temp °C)
# RH_era5 (ERA5 relative humidity %)

library(tidyr)
library(ggpmisc)

df <- df %>%
  # Convert observed VPD from Pa to kPa
  mutate(observed_VPD_kPa = VPD / 1000) %>%
  
  # Calculate VPD from empirical temperature and RH
  mutate(calculated_vpd_empirical = calc_vpd(air_temperature, RH)) %>%
  
  # Calculate VPD from ERA5 temperature and RH
  mutate(calculated_vpd_era5 = calc_vpd(airt_era5, RH_era5)) %>%
  
  # Flag when calculated empirical VPD differs from observed VPD beyond tolerance
  mutate(vpd_diff_flag = abs(calculated_vpd_empirical - observed_VPD_kPa) > vpd_tolerance) %>%
  
  # Flag when ERA5 VPD differs from empirical calculated VPD beyond tolerance (for validation)
  mutate(era5_vpd_diff_flag = abs(calculated_vpd_era5 - calculated_vpd_empirical) > vpd_tolerance) %>%
  
  # Flag temperature differences between ERA5 and empirical beyond tolerance
  mutate(temp_diff_flag = abs(airt_era5 - air_temperature) > temp_tolerance)

# Optional: View flagged rows for review
flagged <- df %>%
  filter(vpd_diff_flag | era5_vpd_diff_flag | temp_diff_flag)

# View(flagged)
# summary table of flags:
flagged_summary <- df %>%
  group_by(management, year) %>%
  summarise(
    n_days                = n(),
    n_VPD_obs_vs_calc     = sum(vpd_diff_flag, na.rm = TRUE),
    n_VPD_obs_vs_ERA5calc = sum(era5_vpd_diff_flag, na.rm = TRUE),
    n_temp_obs_vs_ERA5    = sum(temp_diff_flag, na.rm = TRUE),
    pct_flagged_vpd_calc  = round(100 * n_VPD_obs_vs_calc / n_days, 1),
    pct_flagged_vpd_ERA5  = round(100 * n_VPD_obs_vs_ERA5calc / n_days, 1),
    pct_flagged_temp_ERA5 = round(100 * n_temp_obs_vs_ERA5 / n_days, 1),
    .groups               = "drop"
  )

# plot the linear relationships between observed VPD and calculated: 
plot <- ggplot(df, aes(x = observed_VPD_kPa, y = calculated_vpd_empirical)) +
  
  # first dataset: those VPD's calculated with tower-based RH and air temperature data:
  geom_point(alpha = 0.25, color = "steelblue") +
  
  # second dataset: those VPD's calculated with ERA5-based RH and air temperature data:
  geom_point(aes(x=observed_VPD_kPa, calculated_vpd_era5),
             alpha = 0.15, color = "lightpink")+
 
  # reference line: y = x, perfect 1:1 relationship
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 1, color = "tomato") +
  
  # first linear regression (tower-based data calculation)
  geom_smooth(method = "lm")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x, parse = TRUE, color = "steelblue",

    label.x = 0.05, label.y = 0.95
  ) +
  
  # second linear regression (ERA5-based data calculation):
  geom_smooth(aes(y = calculated_vpd_era5), method = "lm", se = FALSE, color = "maroon") +
  stat_poly_eq(
    aes(y = calculated_vpd_era5, label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x, parse = TRUE, color = "maroon",
    label.x = 0.05, label.y = 1
  ) +
  
  labs(
  x = "Observed VPD (kPa)",
  y = "ERA5 VPD (maroon), calculated with observed climate data (blue) (kPa)"
  ) +
  theme_minimal()

plot

# linear relationship between observed and calculated-with-observed-data VPD = 0.816, R2 = 0.98
# linear relationship between observed and calculated-with-ERA5-data VPD = 0.723, R2 = 0.75

# fit correction model: tower as response variable, ERA5 as predictor
# use only rows where both values are non-NA
era5_correction_VPD <- lm(observed_VPD_kPa ~ calculated_vpd_era5,
                          data = df, na.action = na.omit)
era5_correction_airt <- lm(air_temperature ~ airt_era5,
                           data = df, na.action = na.omit)

# =============================================================================
# Interpretation
# =============================================================================

# It is probable that the difference in VPD between observed (daily average) and calculated using tower RH and air temperature (daily averages)
# is due to a Jensen's inequality -- an averaging order error. The tower's observations are all originally at 30min increments; above, I'm 
# computing VPD from already-averaged-to-daily air temperature and relative humidity. Since Teten's equation is exponential/nonlinear, averaging first and 
# THEN computing VPD underestimates relative to computing VPD first, then averaging. The slope (0.816) reflects that underestimate. No need to fix it,
# as this sanity check and the above explanation of a lack of 1:1 relationship does exactly what we want -- sanity checks.

# as for the difference observed in the ERA5-derived VPD calculation and the observed VPD: this is very likely due to differences of scale and 
# land-surface feedbacks. ERA5 is from 30km gridded modules and represent a spatial average; the tower measures a specific ag surface with more dynamic
# relationships with ecology/biology/geology right beneath. As a result, there is significant scatter especially at higher VPD (hotter, dryer).

# NEXT STEPS: interpolate VPD and air temperature for gaps in tower data, using ERA5 data and using a linear regression correction to correct 
# based on relationship with observed values

# =============================================================================
# calculate and correct: VPD and air temperature where there are tower gaps
# =============================================================================

df <- df %>%
  # remove flag columns for clarity
  select(!c(vpd_diff_flag, era5_vpd_diff_flag, temp_diff_flag, calculated_vpd_empirical, calculated_vpd_era5))

summary(era5_correction_VPD)
summary(era5_correction_airt)

df <- df %>% 
  mutate(
    # ERA5-derived VPD (kPa):
    VPD_era5_calc = calc_vpd(airt_era5, RH_era5),
    
    # bias-corrected ERA5 values for both air temperature and VPD:
    airt_era5_corrected   = coef(era5_correction_airt)[1]+ # intercept
                            coef(era5_correction_airt)[2]*airt_era5, # slope*value
    VPD_era5_corrected    = coef(era5_correction_VPD)[1]+
                            coef(era5_correction_VPD)[2]*VPD_era5_calc,
    VPD_era5_corrected_Pa = VPD_era5_corrected * 1000 , # back to Pa to match VPD in model
    
    # gap filling:
    air_temperature_filled = if_else(is.na(air_temperature),
                                     airt_era5_corrected,
                                     air_temperature), # if air_temperature is NA, fill with corrected ERA5 temp
    VPD_filled = if_else(is.na(VPD),
                         VPD_era5_corrected_Pa,
                         VPD)
  )

# Verify: gap coverage after climate fill
cat("\n--- Climate gap-fill coverage ---\n")
df %>%
  group_by(management, year) %>%
  summarise(
    n_airt_still_na = sum(is.na(air_temperature_filled)),
    n_vpd_still_na  = sum(is.na(VPD_filled)),
    .groups = "drop"
  ) %>% print()


