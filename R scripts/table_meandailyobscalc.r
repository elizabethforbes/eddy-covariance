library(tidyverse)
library(gt)

##############################################################################
# table with number of half-hourly observations in each daily average calc
##############################################################################

# summarize number of daily averages calculated with each number of half-hourly
# observations: total of 1849 obs
# define your observations, binned in groups of 10:
custom_breaks <- c(0, 1, 11, 21, 31, 41, Inf)
custom_labels <- c("0", "1-10", "11-20", "21-30", "31-40", ">40")
table <- data1 %>% 
  group_by(NEE_n) %>% 
  summarise(dailyavgs_obs_n = n(),
            .groups = "drop") %>% 
  mutate(NEE_n_bin = case_when(
    NEE_n == 0 ~ "0",
    is.na(NEE_n) ~ "NA",
    TRUE ~ as.character(cut(NEE_n,
                            breaks = custom_breaks,
                            labels = custom_labels,
                            right = FALSE,
                            include.lowest = TRUE))
  ))

table_binned <- table %>%
  group_by(NEE_n_bin) %>%
  summarise(
    total_days = sum(dailyavgs_obs_n),
    .groups = "drop"
  ) %>%
  mutate(
    percentage = (total_days / sum(total_days))
  ) %>%
  as.data.frame()

table_binned %>%
  gt() %>%
  tab_header(
    title = "Number of Half-hourly Observations per Day",
  ) %>%
  cols_label(
    NEE_n_bin = "Num. of obs.",
    total_days = "Days",
    percentage = "% of Total Days"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  fmt_percent(columns = percentage, decimals = 1)

##############################################################################
# table with proportions of observations from each field, each stage:
##############################################################################

summary_table <- data1 %>%
  drop_na(NEE) %>%  # Remove rows with NA in NEE
  group_by(management, crop_stage_simple) %>%
  summarize(
    count = n(),  # Count rows for each group
    .groups = "drop"  # Drop grouping structure after summarizing
  ) %>%
  group_by(management) %>%  # Group again by management to calculate proportions
  mutate(
    proportion = count / sum(count),  # Calculate proportion within each management group
    proportion = round(proportion, 3)  # Round to 3 decimal places
  ) %>%
  ungroup()  # Remove grouping structure

summary_table %>%
  gt() %>%
  tab_header(
    title = "Count, Proportion of Empirical Observations by Management and Crop Stage"
  ) %>%
  fmt_number(
    columns = proportion,
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_body()
  ) %>%
  cols_label(
    management = "Management",
    crop_stage_simple = "Crop Stage",
    count = "Count",
    proportion = "Proportion"
  ) %>%
  tab_spanner(
    label = "Group Details",
    columns = c(management, crop_stage_simple)
  ) %>%
  tab_spanner(
    label = "Statistics",
    columns = c(count, proportion)
  ) %>% 
  fmt_percent(columns = proportion, decimals = 1)

##############################################################################
# table of mean fluxes by management, crop stage
##############################################################################

library(dplyr)

mean_summary <- data1 %>%
  drop_na(NEE) %>%  # Remove rows with NA in any of these columns
  group_by(management, crop_stage_simple) %>%
  summarize(
    mean_NEE = mean(NEE, na.rm = TRUE),  # Mean NEE
    mean_GPP = mean(GPP, na.rm = TRUE),  # Mean GPP
    mean_Reco = mean(Reco, na.rm = TRUE),    # Mean Re (Ecosystem Respiration)
    .groups = "drop"  # Drop grouping structure after summarizing
  ) %>%
  mutate(
    across(c(mean_NEE, mean_GPP, mean_Reco), ~ round(.x, 3))  # Round to 3 decimal places
  )

library(gt)

# crop_stage_colors <- c(
#   "mature"    = "#FF6B6B",  # Light red
#   "dormant"      = "#4ECDC4",  # Teal
#   "late"     = "#FFE66D",  # Light yellow
#   "senescing" = "#A78BFA", # Light purple
#   "flowering" = "#6BCB77"  # Light green
# )


mean_summary %>%
  gt() %>%
  tab_header(
    title = "Mean NEE, GPP, and Reco by Management and Crop Stage"
    ) %>%
  fmt_number(
    columns = c(mean_NEE, mean_GPP, mean_Reco),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_body()
  ) %>%
  cols_label(
    management = "Management",
    crop_stage_simple = "Crop Stage",
    mean_NEE = "Mean NEE",
    mean_GPP = "Mean GPP",
    mean_Reco = "Mean Reco"
  ) %>%
  tab_spanner(
    label = "Group Details",
    columns = c(management, crop_stage_simple)
  ) %>%
  tab_spanner(
    label = "Mean Fluxes (all μmol CO₂ m⁻² s⁻¹)",
    columns = c(mean_NEE, mean_GPP, mean_Reco)
  )
