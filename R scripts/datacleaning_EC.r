library(tidyverse)
library(readxl)
library(purrr)
library(lubridate)
library(hms)
library(here)

# ============================================================================
# 1. upload eddy covariance data from project folder:
# ============================================================================

# eddy covariance final data:
ec_conv <-  read_csv(here("eddy_covariance_fluxdata", "conv_final.csv"))
ec_conv <- ec_conv %>%
  mutate("management" = "conventional")

ec_org <- read_csv(here("eddy_covariance_fluxdata", "organic_final.csv"))
ec_org <- ec_org %>%
  mutate("management" = "organic")

all_ec <- rbind(ec_conv, ec_org) # 215,064 rows of 155 variables

# select just those variables you'll use in analyses:
all_ec_filtered <- all_ec %>%
  select(date, month, hour, management,
         Reco, co2_flux.c, GPP_f, VPD,
         ndvi, rh, vwc1, vwc2, precip, airt, soilt1, soilt2) %>%
  #add year
  mutate(year = year(date)) %>% 
  # rename a few variables
  rename(NEE = co2_flux.c)

# ============================================================================
# 2. create hourly...
# ============================================================================

ec_hrly_avg <- all_ec_filtered %>% 
  rowwise() %>%
  mutate(soilt_avg = rowMeans(cbind(soilt1, soilt2), na.rm = TRUE),
         vwc_avg = rowMeans(cbind(vwc1, vwc2), na.rm = TRUE)) %>% 
  select(!c(soilt1, soilt2, vwc1, vwc2)) %>% 
  group_by(management, date, hour) %>% 
  summarize(across(is.numeric, mean, na.rm = TRUE), .groups = 'drop') %>% 
  ungroup()

# ============================================================================
# 3. daily...
# ============================================================================

ec_daily_avg <- all_ec_filtered %>% 
  rowwise() %>%
  mutate(soilt_avg = rowMeans(cbind(soilt1, soilt2), na.rm = TRUE),
         vwc_avg = rowMeans(cbind(vwc1, vwc2), na.rm = TRUE)) %>% 
  select(!c(soilt1, soilt2, vwc1, vwc2)) %>% 
  group_by(management, date) %>% 
  summarize(
    Reco_se = sd(Reco, na.rm = TRUE) / sqrt(sum(!is.na(Reco))),
    NEE_se = sd(NEE, na.rm = TRUE) / sqrt(sum(!is.na(NEE))),
    across(is.numeric, mean, na.rm = TRUE),
    .groups = 'drop'
  ) %>% 
  ungroup() %>% 
  select(!c(hour)) %>% 
  mutate(doy = yday(date))
 
# ============================================================================
# 4. and weekly averages of the EC data for further analyses:
# ============================================================================

ec_weekly_avg <- all_ec_filtered %>% 
  rowwise() %>% 
  mutate(soilt_avg = rowMeans(cbind(soilt1, soilt2), na.rm = TRUE),
         vwc_avg = rowMeans(cbind(vwc1, vwc2), na.rm = TRUE)) %>% 
  select(!c(soilt1, soilt2, vwc1, vwc2)) %>% 
  group_by(management, week = floor_date(date, "week")) %>%
  summarize(
    Reco_se = sd(Reco, na.rm = TRUE) / sqrt(sum(!is.na(Reco))),
    NEE_se = sd(NEE, na.rm = TRUE) / sqrt(sum(!is.na(NEE))),
    across(is.numeric, mean, na.rm = TRUE),
    .groups = 'drop'
  ) %>% 
  ungroup() %>% 
  select(!c(hour))

# ============================================================================
# 5. remove unnecessary items from global environment:
# ============================================================================

rm(ec_org, ec_conv, all_ec, all_ec_filtered)
