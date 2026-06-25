# script to merge management data with EC data, Rs data (survey chamber)
library(tidyverse)
library(mgcv)
library(mgcViz)
library(gratia)
library(dplyr)
library(lubridate)
library(readxl)

# ============================================================================
# 1. UPLOAD MANAGEMENT DATA
# ============================================================================

mngmnt <- read_xlsx("management_data_collated.xlsx", sheet = 1)

# ============================================================================
# 2. CREATE MANAGEMENT VARIABLES and APPEND TO EC HOURLY AVERAGE FLUXES DATA
# ============================================================================

# use custom function that calculates days since last management event
source("merge_management_data_timeseries.R")

# select just EFO1, EFO2
mng_2 <- mngmnt %>% 
  filter(management != "dairy pasture")
ec_hrly_2 <- ec_hrly_avg %>% 
  filter(management != "dairy pasture")

# merge management data with the hourly average EC data:
ec_hrly_wmngmnt <- create_management_vars(ec_hrly_2, mng_2)

# includes variables like days since last [fill in the blank] as well as general decay values like
# a 20-day decay for tillage, 15-day decay for fertilizer, etc., to represent the pulses that may occur

# merge these data with the BWS meteorological data. for modeling, use tower-specific data when possible; but where none exist, use the met-tower data
# this is because there is just one met tower, versus site-specific data at each EC tower.
ec_hrly_wmngmnt <- ec_hrly_wmngmnt %>% 
  left_join(bws_avghrly, by = join_by(date, hour))
# add day of year (doy) with lubridate's yday()
ec_hrly_wmngmnt$doy <- yday(ec_hrly_wmngmnt$date)

# cumulative daily NEE and Reco and GPP, plus average for all other vars, and max for PAR: Shahan et al., 2022
ec_daily_cum <- ec_hrly_wmngmnt %>% 
  rowwise() %>%
  group_by(management, date) %>% 
  dplyr::summarize(
            NEE_cum = sum(co2_flux.c, na.rm = TRUE), # daily cumulative NEE
            Reco_cum = sum(Reco, na.rm = TRUE), # daily cumulative Reco
            GPP_cum = sum(GPP_f, na.rm = TRUE), # daily cumulative GPP
            VPD = mean(VPD, na.rm = TRUE),
            NDVI = mean(ndvi, na.rm = TRUE),
            airt_C = mean(airt, na.rm = TRUE),
            RH = mean(rh.x, na.rm = TRUE),
            soilt_C = mean(soilt_avg, na.rm = TRUE),
            VWC = mean(vwc_avg, na.rm = TRUE),
            # management data:
            # last_tillage = last_tillage,
            days_since_tillage = mean(days_since_tillage),
            # last_planting = last_planting,
            days_since_planting = mean(days_since_planting),
            # last_fertilizer = last_fertilizer,
            days_since_fertilizer = mean(days_since_fertilizer),
            # last_herbicide = last_herbicide,
            days_since_herbicide = mean(days_since_herbicide),
            # last_harvest = last_harvest,
            days_since_harvest = mean(days_since_harvest),
            tillage_pulse = mean(tillage_pulse),
            # crop_stage = crop_stage,
            fert_effect = mean(fert_effect),
            herbicide_recent = mean(herbicide_recent),
            post_harvest = mean(post_harvest),
            # BWS tower data:
            pressure_hg = mean(pressure_hg, na.rm = TRUE),
            rain_in = mean(rain_in, na.rm = TRUE),
            PAR = max(par_umolm2sec), # maximum daily PAR
            # extraneous variables:
            month = mean(month.x),
            year = mean(year.x),
            doy = mean(doy))

# ============================================================================
# 2. Experiment with models describing fluxes using this merged hourly dataset
# ============================================================================

# GAMs are weird with categorical variables; convert management to factor
ec_hrly_wmngmnt$management <- factor(ec_hrly_wmngmnt$management)
ec_daily_cum$management <- factor (ec_daily_cum$management)

# Basic model structure: hourly fluxes
model <- gam(co2_flux.c ~ 
               # Field type
               management +
               # Smooth terms for continuous drivers
               s(par_umolm2sec, by = management) +  # photosynthetically active radiation; from met tower
               s(airt, by = management) +  # air temperature
               # s(VPD, by = management) +   # vapor pressure deficit
               s(vwc_avg, by = management) +   # soil water content
               s(soilt_avg, by = management) + # soil temp on average
               # Temporal structure
               s(doy, bs = "cc", by = management) +  # day of year (cyclic)
               s(year.x, bs = "re"),  # random effect for year
             # +
               # Interactions if needed
               # ti(par_umolm2sec, airt, by = management),
             data = ec_hrly_wmngmnt,
             method = "REML",
             family = gaussian())


summary(model)
simulateResiduals(model, plot = TRUE)

# Basic model structure: daily cumulative fluxes
model <- gam(
  # NEE_cum ~ 
  log(Reco_cum) ~
               # Field type
               management +
               # Smooth terms for continuous drivers
               s(PAR, by = management) +  # photosynthetically active radiation; from met tower
               s(airt_C, by = management) +  # air temperature
               # s(VPD, by = management) +   # vapor pressure deficit
               s(VWC, by = management) +   # soil water content
               s(soilt_C, by = management) + # soil temp on average
               # Temporal structure
               s(doy, bs = "cc", by = management) +  # day of year (cyclic)
               s(year, bs = "re"),  # random effect for year
             # +
             # Interactions if needed
             data = ec_daily_cum,
             method = "REML",
             family = gaussian())

summary(model)
simulateResiduals(model, n=1000, plot = TRUE)
testOutliers(model)

# test for temporal autocorrelation:
res <- simulateResiduals(model, plot = F)
# aggregate residuals by subgroup (location), can't test autocorrelation with more than one residual per time
res1 <- recalculateResiduals(res, sel = ec_daily_cum$management == "conventional")
res2 <- recalculateResiduals(res, sel = ec_daily_cum$management == "organic")

testTemporalAutocorrelation(res1, time = 1:length(res1$scaledResiduals))
# DW = 0.32199, p-value < 2.2e-16
# alternative hypothesis: true autocorrelation is not 0
testTemporalAutocorrelation(res2, time = 1:length(res1$scaledResiduals))


# ============================================================================
# 3. Add management data:
# ============================================================================
# Full model with management integration
m1_full <- gam(
  # NEE_cum ~ 
  log(Reco_cum) ~
                 # Fixed effects
                 management +

                 # Environmental responses by field type
                 s(PAR, by = management, k = 20) +
                 s(airt_C, by = management, k = 15) +
                 s(VPD, by = management, k = 10) +
                 s(VWC, by = management, k = 10) +
                 
                 # Management pulses (smooth the decay functions)
                 s(tillage_pulse, by = management, k = 5) +
                 # s(fert_effect, by = management, k = 5) +
                 s(days_since_planting, by = management, k = 5) +
                 # s(herbicide_recent, by = management) +  # linear interaction
                 # post_harvest * management +
                 
                 # Temporal structure
                 s(doy, bs = "cc", k = 20, by = management) +
                 s(year, bs = "re"),
                 # s(hour, bs = "cc", k = 12) +  # diurnal
                 
                 # Key interaction
                 # ti(PAR, airt_C, by = management, k = c(8, 8)),
               
               data = ec_daily_cum,
               method = "REML",
               family = gaussian())

summary(m1_full)
simulateResiduals(m1_full, plot = TRUE)

# ============================================================================
# 4. visualize full model comparing Reco across sites with management data:
# ============================================================================
library(gratia)
library(ggplot2)

# 1. Compare environmental responses
draw(m1_full, select = c("s(VWC):managementorganic", "s(VWC):managementconventional"))

# 2. Management pulse effects
plot_management <- function(model) {
  new_data <- expand.grid(
    management = c("organic", "conventional"),
    tillage_pulse = seq(0, 1, length = 100),
    # Set other variables to median/reference
    PAR = 500,
    airt_C = 20,
    VPD = 1,
    VWC = 0.3,
    soilt_C = 7.8,
    # crop_stage = "vegetative_late",
    fert_effect = 0,
    days_since_planting = 120,
    herbicide_recent = 0,
    post_harvest = 0,
    doy = 180,
    year = "2022",
    hour = 12
  )
  
  pred <- predict(model, newdata = new_data, se.fit = TRUE)
  new_data$fit <- pred$fit
  new_data$se <- pred$se.fit
  
  ggplot(new_data, aes(x = tillage_pulse, y = fit, color = management)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = fit - 1.96*se, ymax = fit + 1.96*se, fill = management), 
                alpha = 0.3) +
    labs(x = "Tillage pulse effect", y = "NEE (µmol m-2 s-1)",
         title = "Effect of tillage on NEE by field type")
}

plot_management(m1_full)

# 3. Phenological trajectories
ggplot(hourly_data_mgmt %>% 
         group_by(management, days_since_planting) %>%
         summarise(NEE = mean(NEE, na.rm = TRUE)),
       aes(x = days_since_planting, y = NEE, color = management)) +
  geom_line() +
  geom_vline(data = management_events %>% filter(event_type == "fertilizer"),
             aes(xintercept = days_since_planting, color = management),
             linetype = "dashed", alpha = 0.5) +
  facet_wrap(~year)