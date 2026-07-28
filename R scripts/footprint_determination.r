library(tidyverse)
library(purrr)
library(here)

read_eddypro <- function(path) {
  headers <- names(read.csv(path, skip = 1, nrow = 0))
  read.csv(path, skip = 3, header = FALSE, col.names = headers,
           na.strings = c("-9999", "")) %>%
    mutate(across(c(qc_co2_flux, u., x_90.), as.numeric))
}

data_dir <- here("eddy_covariance_fluxdata")

fp_files <- list(
  organic      = file.path(data_dir, c("org_2018_fulloutput.csv",
                                       "org_2019_fulloutput.csv",
                                       "org_2020_fulloutput.csv")),
  conventional = file.path(data_dir, c("conv_2018_fulloutput.csv",
                                       "conv_2019_fulloutput_1.csv",
                                       "conv_2019_fulloutput_2.csv",
                                       "conv_2020_fulloutput.csv"))
)

# read data:
fp_data <- imap_dfr(fp_files, function(files, mgmt) {
  map_dfr(files, read_eddypro) %>% 
    mutate(management = mgmt)
})

# Filter to valid flux periods: see Kyle's notes (methods section)
fp_filtered <- fp_data %>%
  filter(qc_co2_flux <= 1,  # remove observations with a QAQC of 2, keeping 0 (excellent) and 1 (good)
         u. >= 0.2,         # keep only observations with friction velocity >= 0.2
         x_90. > 0,         # filter out the -9999s (errors, gaps in the data)
         !is.na(x_90.)) %>% # filter out instances where the 90% distance is NA
  mutate(date = as.Date(date),
         year = year(date))

# Summary: by management × year
fp_summary <- fp_filtered %>%
  filter(year != 2021) %>% 
  group_by(management, year) %>%
  summarise(
    n        = n(),
    median   = median(x_90.),
    q25      = quantile(x_90., 0.25),
    q75      = quantile(x_90., 0.75),
    p95      = quantile(x_90., 0.95),
    .groups  = "drop"
  )

# Summary: by management only
fp_combined <- fp_filtered %>%
  group_by(management) %>%
  summarise(
    n        = n(),
    median   = median(x_90.),
    q25      = quantile(x_90., 0.25),
    q75      = quantile(x_90., 0.75),
    p95      = quantile(x_90., 0.95),
    .groups  = "drop"
  )

print(fp_summary)
print(fp_combined)
