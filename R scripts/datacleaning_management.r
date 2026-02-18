library(tidyverse)
library(readxl)
library(purrr)
library(lubridate)
library(hms)

# ============================================================================
# 1. upload, clean management data:
# ============================================================================

mngmnt <- read_xlsx("management_data_collated.xlsx", sheet = 1)

# Convert date to proper format
mngmnt$date <- as.Date(mngmnt$date, format = "%m/%d/%y")
mngmnt$duration <- as.numeric(mngmnt$duration)
mngmnt <- mngmnt %>% 
  mutate(management = case_when(field == 'EF01' ~ 'organic',
                          field == 'EF02' ~ 'conventional',
                          field == 'EF03' ~ 'dairy pasture'))

# remove dairy pasture data for these analyses, unnecessary columns:
mngmnt <- mngmnt %>% 
  filter(management != "dairy pasture") %>% 
  select(!c(field, EFO3_S2))
