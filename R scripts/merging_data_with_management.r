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

# ec_daily_mngmnt <- create_management_vars_multicrop(
#   ec_daily_avg,
#   mngmnt)

ec_daily_mngmnt2 <- create_management_vars_multicrop(
  ec_daily_avg, mngmnt
)

# adjust the 2020, organic post-barley-harvest crop stage to reflect that the 
# brome/orchardgrass/alfalfa grass field was undersown (then also drill-planted in mid-Sept)
# In spring, the alfalfa will outcompete these grasses and the hay field 
# will be harvested as primarily an alfalfa field. However, the undersown crop will be
# in early vegetative phase (not germination nor mature) after barley is harvested, 
# and will "pop" in a way akin to rapid spring growth before hardening in fall and
# going into dormancy. We want to be able to compare this grass field to the winter 
# rye cover crop going into the conventional field at the same time (end of 2020).

# categorize as "vegetative", and then transition to "mature" for fall hardening,
# before entering dormancy using the logic cues detailed in the supplement:
ec_daily_mngmnt3 <- ec_daily_mngmnt2 %>% 
  mutate(
    crop_stage = case_when(
      management == "organic" &
        date > as.Date("2020-07-08") & # date the barley crop is harvested
        date < as.Date("2020-10-31") # date the grass mix crop goes into fall hardening
      ~ "vegetative",
      TRUE ~ crop_stage),
    crop_stage_simple = case_when(
      management == "organic" &
        date > as.Date("2020-07-08") &
        date < as.Date("2020-10-31")
      ~ "vegetative",
      management == "organic" & crop_stage == "fall hardening"
      ~ "mature",
      TRUE ~ crop_stage_simple
    )
  )
#   

# pull the assigned crop stages from the above to the daily chamber data so that the crop stages match across the two datasets:
chamber_daily_mngmnt <- chamber_daily_site %>%
  # convert date from POSIXct to Date:
  mutate(date = as.Date(date, tz = "EST")) %>%
  left_join(ec_daily_mngmnt3 %>% 
              select(management, date, current_crop, crop_stage, crop_stage_simple),
            by = c("management", "date"))
