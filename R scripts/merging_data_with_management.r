library(tidyverse)
library(readxl)
library(readr)
library(here)
library(purrr)
library(lubridate)
library(hms)

source(here("R scripts", "datacleaning_EC.r"))
source(here("R scripts", "datacleaning_chambers.r"))
source(here("R scripts", "datacleaning_management.r"))

# Using a custom function, calculate continuous variables for each management action. 
# The custom function is informed by characteristics like the expected 
# crop-specific phenology, and days since actions like tillage.
source(here("R scripts", "merge_management_data_covercropsandperennials_v3_6_0_6.r"))

# merge continuous daily average flux data with management data using above function:
ec_daily_mngmnt <- create_management_vars_multicrop(
  ec_daily_avg,
  mngmnt)

# manually adjust post-9/15/2020 crop stages to all 'mature' before fall hardening, 
# to reflect that the alfalfa (or more likely, orchardgrass) was drill-planted 
# into an already-mature grass field that was a mix of brome, orchardgrass, and 
# alfalfa. In spring, the alfalfa will outcompete these grasses and the hay field 
# will be harvested as primarily an alfalfa field.
write_csv(ec_daily_mngmnt, "ec_with_mngngment.csv")

# edited manually and re-read:
ec_daily_mngmnt <- read_csv("ec_with_mngngment_final.csv")
# make sure date is in date class:
ec_daily_mngmnt <- ec_daily_mngmnt %>% 
  mutate(date = as.Date(date, tz = "EST", format = "%m/%d/%y"))

# pull the assigned crop stages from the above to the daily chamber data so that the crop stages match across the two datasets:
chamber_daily_mngmnt <- chamber_daily_site %>%
  # convert date from POSIXct to Date:
  mutate(date = as.Date(date, tz = "EST")) %>%
  left_join(ec_daily_mngmnt %>% 
              select(management, date, current_crop, crop_stage, crop_stage_simple),
            by = c("management", "date"))