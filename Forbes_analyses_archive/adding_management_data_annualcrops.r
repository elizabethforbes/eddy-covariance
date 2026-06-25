# this script uses the custom merge_management_data_timeseries.R script to merge the management data we have (individual dates with actions) to the GHG exchange data
# (continuous over several years) using a range of phenology logic, etc.

# Load necessary libraries
library(dplyr)
library(lubridate)
library(zoo)  # for rolling functions
library(readxl)
library(lubridate)

source(here("R scripts", "merge_management_data_timeseries_FINAL.r"))

source(here("R scripts", "datacleaning_management.r"))

source(here("R scripts", "datacleaning_EC.r"))

# merge daily average EC data with management:
ec_daily_mngmnt <- create_management_vars_multicrop(
  ec_daily_avg,
  mngmnt) 

# same for chamber data, and convert Date to as.Date:
chamber_daily_mngmnt <- chamber_daily_site %>% 
  mutate(date = as.Date(date, tz = "EST")) %>% 
  create_management_vars_multicrop(mngmnt)


