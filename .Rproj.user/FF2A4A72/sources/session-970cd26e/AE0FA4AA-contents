# ============================================================================
# Complete Gap-Filling Workflow for EC Daily Data
# Fills missing dates and values using GAM predictions
# ============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)  # for interpolation
library(ggplot2)

# Assuming you have:
# - ec_daily_mngmnt2: your original data
# - g_NEE, g_GPP, g_Reco: your fitted GAM models

# ============================================================================
# STEP 1: Create Complete Date Sequence
# ============================================================================

cat("=== STEP 1: Creating Complete Date Sequence ===\n\n")

# Get date range for each management system
date_range <- ec_daily_mngmnt2 %>%
  group_by(management) %>%
  summarise(
    min_date = min(date, na.rm = TRUE),
    max_date = max(date, na.rm = TRUE),
    n_original = n(),
    .groups = "drop"
  )

print(date_range)

# Create complete daily sequence for each management
complete_dates <- date_range %>%
  rowwise() %>%
  mutate(
    date_seq = list(seq.Date(min_date, max_date, by = "day"))
  ) %>%
  unnest(date_seq) %>%
  select(management, date = date_seq)

cat("\nOriginal dataset rows:", nrow(ec_daily_mngmnt2), "\n")
cat("Complete sequence rows:", nrow(complete_dates), "\n")
cat("Missing dates to fill:", nrow(complete_dates) - nrow(ec_daily_mngmnt2), "\n\n")

# Merge with original data
ec_complete <- complete_dates %>%
  left_join(ec_daily_mngmnt2, by = c("management", "date")) %>%
  mutate(
    # Flag row type
    data_exists = !is.na(NEE) | !is.na(GPP_f) | !is.na(Reco),
    row_type = ifelse(data_exists, "original_data", "added_missing_date")
  )

# Summary
gap_summary <- ec_complete %>%
  group_by(management, row_type) %>%
  summarise(n = n(), .groups = "drop")

cat("Breakdown by management:\n")
print(gap_summary)
cat("\n")

# ============================================================================
# STEP 2A: OPTION 1 - Interpolate Predictor Variables from Surrounding Days
# ============================================================================

cat("=== STEP 2A: Interpolating Predictor Variables ===\n\n")

# Identify which predictors need interpolation
predictors_to_interpolate <- c("VPD", "soilt_avg", "vwc_avg", "airt", 
                                "precip", "rh", "ndvi")

# Check missingness before interpolation
cat("Missingness BEFORE interpolation:\n")
missing_before <- ec_complete %>%
  summarise(across(all_of(predictors_to_interpolate), 
                   ~sum(is.na(.))))
print(t(missing_before))

# Interpolate continuous variables
ec_complete <- ec_complete %>%
  arrange(management, date) %>%
  group_by(management) %>%
  mutate(
    # Linear interpolation (na.approx)
    # maxgap parameter limits how many consecutive NAs to interpolate
    VPD = na.approx(VPD, na.rm = FALSE, maxgap = 7),  # max 7 day gap
    soilt_avg = na.approx(soilt_avg, na.rm = FALSE, maxgap = 7),
    vwc_avg = na.approx(vwc_avg, na.rm = FALSE, maxgap = 7),
    airt = na.approx(airt, na.rm = FALSE, maxgap = 7),
    precip = na.approx(precip, na.rm = FALSE, maxgap = 7),
    rh = na.approx(rh, na.rm = FALSE, maxgap = 7),
    ndvi = na.approx(ndvi, na.rm = FALSE, maxgap = 14),  # NDVI changes slower
  ) %>%
  ungroup()

# Check missingness after interpolation
cat("\nMissingness AFTER interpolation:\n")
missing_after <- ec_complete %>%
  summarise(across(all_of(predictors_to_interpolate), 
                   ~sum(is.na(.))))
print(t(missing_after))

# Fill categorical variables (forward then backward fill)
ec_complete <- ec_complete %>%
  arrange(management, date) %>%
  group_by(management) %>%
  fill(crop_stage_simple, .direction = "down") %>%
  fill(crop_stage_simple, .direction = "up") %>%
  fill(current_crop, .direction = "down") %>%
  fill(current_crop, .direction = "up") %>%
  fill(crop_type, .direction = "down") %>%
  fill(crop_type, .direction = "up") %>%
  fill(in_season, .direction = "down") %>%
  fill(in_season, .direction = "up") %>%
  ungroup()

# Recalculate time-based variables
ec_complete <- ec_complete %>%
  mutate(
    year = as.factor(year(date)),
    month = month(date),
    doy = yday(date)
  )

# Handle tillage timing (if applicable)
# For new rows, estimate based on last known tillage
ec_complete <- ec_complete %>%
  arrange(management, date) %>%
  group_by(management) %>%
  fill(last_tillage, .direction = "down") %>%
  mutate(
    last_tillage = as.Date(last_tillage),
    days_since_tillage = as.numeric(date - last_tillage)
  ) %>%
  ungroup()

# Flag which variables were interpolated
ec_complete <- ec_complete %>%
  mutate(
    interpolated_predictors = case_when(
      row_type == "added_missing_date" ~ "all",
      !data_exists ~ "some",
      TRUE ~ "none"
    )
  )

cat("\nInterpolation complete!\n\n")

# Save checkpoint
write.csv(ec_complete %>% select(management, date, VPD, soilt_avg, vwc_avg, 
                                  crop_stage_simple, row_type),
          "/home/claude/checkpoint_interpolated_predictors.csv",
          row.names = FALSE)

# ============================================================================
# STEP 2B: OPTION 2 - Load External Weather Data (if available)
# ============================================================================

cat("=== STEP 2B: External Data Option (PLACEHOLDER) ===\n\n")

# If you have external data, load it here and merge
# Example structure:

# external_weather <- read.csv("path/to/external_weather.csv")
# external_weather$date <- as.Date(external_weather$date)
# 
# ec_complete <- ec_complete %>%
#   left_join(external_weather %>% 
#               select(date, VPD_external, temp_external, precip_external),
#             by = "date") %>%
#   mutate(
#     # Use external data to fill gaps, prefer on-site when available
#     VPD = coalesce(VPD, VPD_external),
#     airt = coalesce(airt, temp_external),
#     precip = coalesce(precip, precip_external)
#   )

cat("To use external data:\n")
cat("  1. Load external weather data as a dataframe with 'date' column\n")
cat("  2. Uncomment and modify the code above\n")
cat("  3. Use coalesce() to prefer on-site data, fill with external\n\n")

# ============================================================================
# STEP 3: Quality Check on Predictor Variables
# ============================================================================

cat("=== STEP 3: Checking Predictor Quality ===\n\n")

# Check for remaining NAs in key predictors
final_nas <- ec_complete %>%
  summarise(
    n_total = n(),
    across(c(VPD, soilt_avg, vwc_avg, crop_stage_simple, doy, year),
           ~sum(is.na(.)), .names = "na_{.col}")
  )

print(final_nas)

# If still have NAs, identify problematic periods
if(sum(final_nas[,grep("^na_", names(final_nas))]) > 0) {
  cat("\nWARNING: Still have missing predictor values!\n")
  cat("Identifying problematic periods:\n\n")
  
  problem_rows <- ec_complete %>%
    filter(if_any(c(VPD, soilt_avg, vwc_avg, crop_stage_simple), is.na)) %>%
    select(management, date, VPD, soilt_avg, vwc_avg, crop_stage_simple, row_type)
  
  print(head(problem_rows, 20))
  
  cat("\nOptions:\n")
  cat("  1. Use longer maxgap in interpolation\n")
  cat("  2. Get external weather data for these periods\n")
  cat("  3. Exclude these dates from gap-filling\n\n")
}

# Visualize interpolated vs original values
p_interp_check <- ec_complete %>%
  filter(year == 2019) %>%  # Example year
  ggplot(aes(x = date, y = VPD)) +
  geom_line(color = "gray70") +
  geom_point(data = ec_complete %>% 
               filter(year == 2019, row_type == "original_data"),
             color = "black", size = 1) +
  geom_point(data = ec_complete %>%
               filter(year == 2019, row_type == "added_missing_date"),
             color = "red", size = 1, alpha = 0.5) +
  facet_wrap(~management, ncol = 1) +
  labs(title = "VPD: Original (black) vs Interpolated (red) - 2019 Example",
       x = "Date", y = "VPD (Pa)") +
  theme_bw()

ggsave("/home/claude/interpolation_check_VPD.png", p_interp_check,
       width = 12, height = 6, dpi = 300)

cat("Saved: interpolation_check_VPD.png\n\n")

# ============================================================================
# STEP 4: Generate GAM Predictions for ALL Rows
# ============================================================================

cat("=== STEP 4: Generating GAM Predictions ===\n\n")

# Check that all predictors needed by models are present
# Get model formula to see what's needed
# formula(g_NEE)  # Uncomment to check

cat("Predicting NEE...\n")
pred_NEE <- predict(g_NEE, newdata = ec_complete, se.fit = TRUE)
ec_complete$NEE_predicted <- pred_NEE$fit
ec_complete$NEE_se <- pred_NEE$se.fit

cat("Predicting GPP...\n")
pred_GPP <- predict(g_GPP, newdata = ec_complete, se.fit = TRUE)
ec_complete$GPP_predicted <- pred_GPP$fit
ec_complete$GPP_se <- pred_GPP$se.fit

cat("Predicting Reco...\n")
pred_Reco <- predict(g_Reco, newdata = ec_complete, se.fit = TRUE)
ec_complete$Reco_predicted <- pred_Reco$fit
ec_complete$Reco_se <- pred_Reco$se.fit

cat("Predictions complete!\n\n")

# Check for prediction failures (NAs)
pred_nas <- ec_complete %>%
  summarise(
    na_NEE_pred = sum(is.na(NEE_predicted)),
    na_GPP_pred = sum(is.na(GPP_predicted)),
    na_Reco_pred = sum(is.na(Reco_predicted))
  )

if(sum(pred_nas) > 0) {
  cat("WARNING: Some predictions failed (returned NA):\n")
  print(pred_nas)
  cat("\nThis likely means predictor values were missing or out of range\n\n")
}

# ============================================================================
# STEP 5: Create Gap-Filled Columns with Sanity Check
# ============================================================================

cat("=== STEP 5: Creating Gap-Filled Dataset ===\n\n")

ec_complete <- ec_complete %>%
  mutate(
    # Keep original measured values
    NEE_original = NEE,
    GPP_original = GPP_f,
    Reco_original = Reco,
    
    # Gap-filled: use measured when available, predicted when missing
    NEE_gapfilled = ifelse(!is.na(NEE), NEE, NEE_predicted),
    GPP_gapfilled = ifelse(!is.na(GPP_f), GPP_f, GPP_predicted),
    Reco_gapfilled = ifelse(!is.na(Reco), Reco, Reco_predicted),
    
    # SANITY CHECK: Calculate NEE from Reco - GPP
    NEE_from_GppReco = Reco_gapfilled - GPP_gapfilled,
    
    # Difference between the two estimates
    sanity_check_diff = NEE_gapfilled - NEE_from_GppReco,
    sanity_check_pct = 100 * abs(sanity_check_diff) / (abs(NEE_gapfilled) + 0.1),
    
    # Flag data source
    NEE_source = case_when(
      !is.na(NEE) ~ "measured",
      is.na(NEE) & row_type == "original_data" ~ "gapfilled_missing_value",
      row_type == "added_missing_date" ~ "gapfilled_missing_date",
      TRUE ~ "unknown"
    ),
    
    GPP_source = case_when(
      !is.na(GPP_f) ~ "measured",
      is.na(GPP_f) & row_type == "original_data" ~ "gapfilled_missing_value",
      row_type == "added_missing_date" ~ "gapfilled_missing_date",
      TRUE ~ "unknown"
    ),
    
    Reco_source = case_when(
      !is.na(Reco) ~ "measured",
      is.na(Reco) & row_type == "original_data" ~ "gapfilled_missing_value",
      row_type == "added_missing_date" ~ "gapfilled_missing_date",
      TRUE ~ "unknown"
    )
  )

cat("Gap-filled dataset created!\n\n")

# ============================================================================
# STEP 6: Sanity Check Analysis
# ============================================================================

cat("=== STEP 6: Sanity Check - NEE_predicted vs NEE_from_GppReco ===\n\n")

# Overall statistics
sanity_stats <- ec_complete %>%
  filter(!is.na(NEE_gapfilled) & !is.na(NEE_from_GppReco)) %>%
  summarise(
    n_total = n(),
    mean_abs_diff = mean(abs(sanity_check_diff), na.rm = TRUE),
    median_abs_diff = median(abs(sanity_check_diff), na.rm = TRUE),
    max_abs_diff = max(abs(sanity_check_diff), na.rm = TRUE),
    mean_pct_diff = mean(sanity_check_pct, na.rm = TRUE),
    n_large_diff = sum(abs(sanity_check_diff) > 5, na.rm = TRUE),
    pct_large_diff = round(100 * n_large_diff / n_total, 1)
  )

print(sanity_stats)
cat("\n")

# By management
sanity_by_mgmt <- ec_complete %>%
  filter(!is.na(NEE_gapfilled)) %>%
  group_by(management) %>%
  summarise(
    mean_abs_diff = mean(abs(sanity_check_diff), na.rm = TRUE),
    max_abs_diff = max(abs(sanity_check_diff), na.rm = TRUE),
    mean_pct_diff = mean(sanity_check_pct, na.rm = TRUE),
    .groups = "drop"
  )

cat("By management:\n")
print(sanity_by_mgmt)
cat("\n")

# Identify outliers (large discrepancy)
outliers <- ec_complete %>%
  filter(abs(sanity_check_diff) > 5) %>%  # >5 g C m⁻² d⁻¹ difference
  arrange(desc(abs(sanity_check_diff))) %>%
  select(management, date, NEE_gapfilled, NEE_from_GppReco, sanity_check_diff,
         NEE_source, crop_stage_simple, VPD, soilt_avg)

if(nrow(outliers) > 0) {
  cat("WARNING:", nrow(outliers), 
      "days with large discrepancy (|diff| > 5 g C m⁻² d⁻¹):\n")
  print(head(outliers, 15))
  cat("\nConsider investigating these dates for:\n")
  cat("  - Extreme environmental conditions\n")
  cat("  - Model extrapolation beyond training range\n")
  cat("  - Predictor interpolation issues\n\n")
  
  write.csv(outliers, "/home/claude/sanity_check_outliers.csv", 
            row.names = FALSE)
}

# ============================================================================
# STEP 7: Quality Checks - Prediction Ranges
# ============================================================================

cat("=== STEP 7: Checking Prediction Ranges ===\n\n")

# Are predictions within reasonable bounds?
range_check <- ec_complete %>%
  summarise(
    # NEE
    NEE_meas_min = min(NEE_original, na.rm = TRUE),
    NEE_meas_max = max(NEE_original, na.rm = TRUE),
    NEE_pred_min = min(NEE_predicted, na.rm = TRUE),
    NEE_pred_max = max(NEE_predicted, na.rm = TRUE),
    
    # GPP
    GPP_meas_min = min(GPP_original, na.rm = TRUE),
    GPP_meas_max = max(GPP_original, na.rm = TRUE),
    GPP_pred_min = min(GPP_predicted, na.rm = TRUE),
    GPP_pred_max = max(GPP_predicted, na.rm = TRUE),
    
    # Reco
    Reco_meas_min = min(Reco_original, na.rm = TRUE),
    Reco_meas_max = max(Reco_original, na.rm = TRUE),
    Reco_pred_min = min(Reco_predicted, na.rm = TRUE),
    Reco_pred_max = max(Reco_predicted, na.rm = TRUE)
  )

cat("Range comparison (Measured vs Predicted):\n")
print(t(range_check))
cat("\n")

# Flag if predictions go beyond measured range
if(range_check$NEE_pred_min < range_check$NEE_meas_min - 5 ||
   range_check$NEE_pred_max > range_check$NEE_meas_max + 5) {
  cat("WARNING: NEE predictions extend beyond measured range!\n")
  cat("This suggests model extrapolation - check these periods\n\n")
}

# Distribution comparison
p_distribution <- ec_complete %>%
  select(management, NEE_original, NEE_gapfilled, NEE_from_GppReco) %>%
  pivot_longer(cols = c(NEE_original, NEE_gapfilled, NEE_from_GppReco),
               names_to = "type", values_to = "NEE") %>%
  mutate(type = factor(type, levels = c("NEE_original", "NEE_gapfilled", 
                                         "NEE_from_GppReco"),
                       labels = c("Measured", "Gap-filled (GAM)", 
                                  "Gap-filled (Reco-GPP)"))) %>%
  ggplot(aes(x = NEE, fill = type)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~management) +
  scale_fill_manual(values = c("Measured" = "black",
                                "Gap-filled (GAM)" = "blue",
                                "Gap-filled (Reco-GPP)" = "red")) +
  labs(title = "Distribution Comparison: Measured vs Gap-Filled NEE",
       x = "NEE (g C m⁻² d⁻¹)",
       fill = "Data Type") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("/home/claude/distribution_comparison.png", p_distribution,
       width = 10, height = 5, dpi = 300)

cat("Saved: distribution_comparison.png\n\n")

# ============================================================================
# STEP 8: Visualizations
# ============================================================================

cat("=== STEP 8: Creating Visualizations ===\n\n")

# PLOT 1: Sanity check scatter
p_sanity <- ec_complete %>%
  filter(NEE_source != "measured") %>%  # Only gap-filled
  ggplot(aes(x = NEE_gapfilled, y = NEE_from_GppReco)) +
  geom_point(aes(color = NEE_source), alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~management) +
  scale_color_manual(values = c("gapfilled_missing_value" = "blue",
                                 "gapfilled_missing_date" = "orange")) +
  labs(title = "Sanity Check: NEE_predicted vs NEE_from_GppReco",
       subtitle = "Gap-filled values only - should fall near 1:1 line",
       x = "NEE from GAM prediction (g C m⁻² d⁻¹)",
       y = "NEE from Reco - GPP (g C m⁻² d⁻¹)",
       color = "Gap-fill Type") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("/home/claude/sanity_check_scatter.png", p_sanity,
       width = 10, height = 5, dpi = 300)
cat("Saved: sanity_check_scatter.png\n")

# PLOT 2: Time series with gap-filled sections
plot_year <- 2019  # Choose a representative year

p_timeseries <- ec_complete %>%
  filter(year(date) == plot_year) %>%
  ggplot(aes(x = date)) +
  # Ribbon showing uncertainty
  geom_ribbon(aes(ymin = NEE_gapfilled - 1.96*NEE_se,
                  ymax = NEE_gapfilled + 1.96*NEE_se,
                  fill = NEE_source), alpha = 0.3) +
  # Gap-filled line (GAM predictions)
  geom_line(aes(y = NEE_gapfilled), color = "blue", size = 0.8) +
  # Sanity check line (Reco - GPP)
  geom_line(aes(y = NEE_from_GppReco), color = "red", 
            size = 0.5, linetype = "dashed", alpha = 0.7) +
  # Original measured points on top
  geom_point(aes(y = NEE_original), size = 1, color = "black") +
  facet_wrap(~management, ncol = 1) +
  scale_fill_manual(values = c("measured" = "transparent",
                                "gapfilled_missing_value" = "lightblue",
                                "gapfilled_missing_date" = "pink"),
                    labels = c("Measured", "Gap-filled (missing value)",
                               "Gap-filled (missing date)")) +
  labs(title = paste("Gap-Filled NEE Time Series:", plot_year),
       subtitle = "Black points = measured | Blue line = GAM predictions | Red dashed = Reco-GPP check",
       x = "Date",
       y = "NEE (g C m⁻² d⁻¹)",
       fill = "Data Source") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("/home/claude/timeseries_gapfilled_2019.png", p_timeseries,
       width = 14, height = 8, dpi = 300)
cat("Saved: timeseries_gapfilled_2019.png\n")

# PLOT 3: Sanity check difference over time
p_diff_time <- ec_complete %>%
  filter(year(date) == plot_year, NEE_source != "measured") %>%
  ggplot(aes(x = date, y = sanity_check_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-5, 5), linetype = "dotted", color = "red") +
  geom_line(aes(color = NEE_source)) +
  geom_point(aes(color = NEE_source), size = 1, alpha = 0.5) +
  facet_wrap(~management, ncol = 1) +
  scale_color_manual(values = c("gapfilled_missing_value" = "blue",
                                 "gapfilled_missing_date" = "orange")) +
  labs(title = paste("Sanity Check Difference Over Time:", plot_year),
       subtitle = "Difference = NEE_predicted - NEE_from_GppReco | Red lines = ±5 threshold",
       x = "Date",
       y = "Difference (g C m⁻² d⁻¹)",
       color = "Gap-fill Type") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("/home/claude/sanity_check_timeseries.png", p_diff_time,
       width = 14, height = 6, dpi = 300)
cat("Saved: sanity_check_timeseries.png\n\n")

# ============================================================================
# STEP 9: Summary Statistics
# ============================================================================

cat("=== STEP 9: Summary Statistics ===\n\n")

# Gap-filling summary by year
gapfill_summary <- ec_complete %>%
  group_by(management, year(date), NEE_source) %>%
  summarise(n_days = n(), .groups = "drop") %>%
  rename(year = `year(date)`) %>%
  pivot_wider(names_from = NEE_source, values_from = n_days, values_fill = 0) %>%
  mutate(
    total_days = measured + gapfilled_missing_value + gapfilled_missing_date,
    pct_measured = round(100 * measured / total_days, 1),
    pct_gapfilled_total = round(100 * (gapfilled_missing_value + gapfilled_missing_date) / total_days, 1)
  )

print(gapfill_summary)

write.csv(gapfill_summary, "/home/claude/gapfill_summary_by_year.csv",
          row.names = FALSE)
cat("\nSaved: gapfill_summary_by_year.csv\n\n")

# ============================================================================
# STEP 10: Save Gap-Filled Datasets
# ============================================================================

cat("=== STEP 10: Saving Gap-Filled Datasets ===\n\n")

# Full dataset with all columns
write.csv(ec_complete, "/home/claude/ec_daily_COMPLETE_gapfilled_FULL.csv",
          row.names = FALSE)
cat("Saved: ec_daily_COMPLETE_gapfilled_FULL.csv\n")

# Clean version with essential columns
ec_gapfilled_clean <- ec_complete %>%
  select(
    # Identifiers
    management, date, year, month, doy,
    
    # Original measured values
    NEE_original, GPP_original, Reco_original,
    
    # Gap-filled values (PRIMARY - use these for analyses!)
    NEE_gapfilled, GPP_gapfilled, Reco_gapfilled,
    
    # Sanity check
    NEE_from_GppReco, sanity_check_diff,
    
    # Uncertainty estimates
    NEE_se, GPP_se, Reco_se,
    
    # Metadata
    NEE_source, GPP_source, Reco_source,
    row_type, interpolated_predictors,
    
    # Key predictors (for reference)
    crop_stage_simple, VPD, soilt_avg, vwc_avg
  )

write.csv(ec_gapfilled_clean, "/home/claude/ec_daily_gapfilled_CLEAN.csv",
          row.names = FALSE)
cat("Saved: ec_daily_gapfilled_CLEAN.csv\n\n")

# ============================================================================
# STEP 11: Create Summary Report
# ============================================================================

sink("/home/claude/gapfilling_report.txt")

cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("GAP-FILLING SUMMARY REPORT\n")
cat("Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

cat("DATASET OVERVIEW:\n")
cat("-" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("Total rows in complete dataset:", nrow(ec_complete), "\n")
cat("Date range:", min(ec_complete$date), "to", max(ec_complete$date), "\n")
cat("Management systems:", paste(unique(ec_complete$management), collapse = ", "), "\n\n")

cat("GAP-FILLING SUMMARY:\n")
cat("-" %>% rep(70) %>% paste(collapse = ""), "\n")
print(gapfill_summary)

cat("\n\nSANITY CHECK RESULTS:\n")
cat("-" %>% rep(70) %>% paste(collapse = ""), "\n")
print(sanity_stats)

cat("\n\nPREDICTION RANGE CHECK:\n")
cat("-" %>% rep(70) %>% paste(collapse = ""), "\n")
print(t(range_check))

cat("\n\nFILES CREATED:\n")
cat("-" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("1. ec_daily_COMPLETE_gapfilled_FULL.csv - Complete dataset with all columns\n")
cat("2. ec_daily_gapfilled_CLEAN.csv - Essential columns only (USE THIS!)\n")
cat("3. gapfill_summary_by_year.csv - Summary statistics\n")
cat("4. sanity_check_outliers.csv - Days with large discrepancy\n")
cat("5. checkpoint_interpolated_predictors.csv - Predictor interpolation check\n")
cat("6. Multiple diagnostic plots (.png)\n\n")

cat("RECOMMENDATIONS:\n")
cat("-" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("1. Review sanity_check_outliers.csv for problematic predictions\n")
cat("2. Examine diagnostic plots for quality assessment\n")
cat("3. Use NEE_gapfilled (not NEE_from_GppReco) for analyses\n")
cat("4. Report gap-filling methodology in methods section\n")
cat("5. Include uncertainty (NEE_se) in cumulative calculations\n\n")

cat("NEXT STEPS:\n")
cat("-" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("- Calculate cumulative fluxes using gap-filled data\n")
cat("- Repeat for other crop stages if needed\n")
cat("- Compare results with/without gap-filling\n")

sink()

cat("\n" %>% rep(70) %>% paste(collapse = "="), "\n")
cat("GAP-FILLING COMPLETE!\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

cat("Output files created in /home/claude/\n")
cat("Review the report and diagnostic plots before proceeding.\n\n")
