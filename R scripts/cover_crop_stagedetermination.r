# BASE CODE CO-GENERATED WITH CLAUDE.AI, adjusted as needed

# Load necessary libraries
library(dplyr)
library(lubridate)
library(zoo)  # for rolling functions
library(readxl)
library(lubridate)

# Read data
# df <- read_xlsx("covercrops_notes.xlsx", sheet = 1)
# df <- read_xlsx("covercrops_notes.xlsx", sheet = 2)
df <- read_xlsx("covercrops_notes.xlsx", sheet = 3)
# df <- read_xlsx("covercrops_notes.xlsx", sheet = 4)

# Ensure date column is POSIXct
df$date <- as.POSIXct(df$date)

# Planting and termination dates
# planting_date <- as.POSIXct("2018-10-27") # organic 1
# termination_date <- as.POSIXct("2019-05-08") # organic 1
# planting_date <- as.POSIXct("2019-11-11") # organic 2
# termination_date <- as.POSIXct("2020-03-27") # organic 2
planting_date <- as.POSIXct("2019-10-15") # conventional 1
termination_date <- as.POSIXct("2020-04-01") # conventional 1
# planting_date <- as.POSIXct("2020-07-20") # conventional 2
# termination_date <- as.POSIXct("2021-05-31") # conventional 2

# Smooth soil temp and NDVI before loop
df$soil_temp_smooth <- rollmean(df$soilt_avg, k = 5, fill = NA, align = "right")
df$ndvi_smooth <- rollmean(df$ndvi, k = 5, fill = NA, align = "right")

# Define ordered phases (all lowercase)
phase_order <- c("pre-planting", "germination/emergence", "fall tillering", 
                 "dormant", "spring growth", "rapid growth", "terminated", "transition", "insufficient data")

# Helper function to get phase index (case-insensitive)
get_phase_index <- function(phase) {
  match(tolower(phase), phase_order)
}

current_phase <- NA
df$stage <- NA_character_

for(i in seq_len(nrow(df))) {
  current_date <- df$date[i]
  ndvi <- df$ndvi[i]
  soil_temp <- df$soilt_avg[i]
  soil_temp_smooth <- df$soil_temp_smooth[i]
  
  # Skip if data missing
  if(is.na(ndvi) | is.na(soil_temp) | is.na(current_date)) {
    df$stage[i] <- "insufficient data"
    next
  }
  
  # Determine new phase
  new_phase <- NA
  
  if(current_date < planting_date - days(1)) {
    new_phase <- "pre-planting"
  } else if(current_date >= termination_date) {
    new_phase <- "terminated"
  } else if(ndvi < 0.20 & soil_temp < 0.55) {
    new_phase <- "dormant"
  } else if(ndvi >= 0.25 & current_date <= planting_date + days(10)) {
    new_phase <- "germination/emergence"
  } else if(ndvi >= 0.20 & soil_temp >= 0.55) {   # Raised threshold to 34°F, or 0.55C according to HV seed company literature
    new_phase <- "fall tillering"
  } else {
    if(soil_temp_smooth >= 1 & 
       ndvi >= 0.2 & 
       !(month(current_date) %in% c(9, 10, 11, 12))) {
      new_phase <- "spring growth"
    } else if(ndvi > 0.25 & 
              soil_temp_smooth > 2) {
      new_phase <- "rapid growth"
    } else {
      new_phase <- "transition"
    }
  }
  
  # Update logic with sticky transition and no backward move dormant->fall tillering
  if(is.na(current_phase)) {
    current_phase <- new_phase
  } else if(new_phase == "transition") {
    if(current_phase == "transition" || is.na(current_phase)) {
      current_phase <- new_phase
    }
  } else {
    if(current_phase == "dormant" && new_phase == "fall tillering") {
      # Ignore backward move, keep dormant
      cat(sprintf("[Info] Ignored backward transition from dormant to fall tillering on %s\n", current_date))
    } else if(get_phase_index(new_phase) > get_phase_index(current_phase)) {
      cat(sprintf("[Phase change] %s -> %s on %s\n", current_phase, new_phase, current_date))
      current_phase <- new_phase
    } else if(new_phase == current_phase) {
      # Optional: update current phase if same phase detected again
      current_phase <- new_phase
    }
  }
  
  df$stage[i] <- current_phase
}


# Optional: Save results to CSV
# write.csv(df, "winter_rye_stages_efo1_1.csv", row.names = FALSE)
# write.csv(df, "winter_rye_stages_efo1_2.csv", row.names = FALSE)
write.csv(df, "winter_rye_stages_efo2_1.csv", row.names = FALSE)
# write.csv(df, "winter_rye_stages_efo2_2.csv", row.names = FALSE)



