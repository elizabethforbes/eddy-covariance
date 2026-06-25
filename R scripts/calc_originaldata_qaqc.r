library(tidyverse)
library(readr)
library(here)

org <- read_csv(here("eddy_covariance_fluxdata", "organic_final.csv"))
conv <- read_csv(here("eddy_covariance_fluxdata", "conv_final.csv"))

# organic: all
org %>% filter(!is.na(co2_flux)) %>% summarize(n_obs = n())
# 40793
# organic: cleaned
org %>% filter(!is.na(co2_flux.c)) %>% summarize(n_obs = n())
# 33496

# percent of org half-hourly obs retained after QAQC:
(33496/40793)*100
# 82.11213

# conventional: all
conv %>% filter(!is.na(co2_flux)) %>% summarize(n_obs = n())
# 46660
# conventional: cleaned
conv %>% filter(!is.na(co2_flux.c)) %>% summarize(n_obs = n())
# 37622

# percent of conv half-hourly obs retained after QAQC:
(37622/46660)*100
# 80.63009

# total percentage, org and conv:
(33496+37622)/(40793+46660)*100
# 81.3214
