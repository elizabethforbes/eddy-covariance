# ============================================================================
# Year-Specific Resilience Analysis
# Examining how management × VPD responses vary by year (drought vs. normal)
# ============================================================================

library(mgcv)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Assuming you have already:
# 1. Loaded your data as ec_daily_mngmnt2
# 2. Created drought_stress variable
# 3. Fitted g_resilience model

# ============================================================================
# STEP 1: Identify Drought vs. Normal Years
# ============================================================================

cat("=== STEP 1: Characterizing Years ===\n\n")

year_characteristics <- ec_daily_mngmnt2 %>%
  group_by(year) %>%
  summarise(
    # VPD metrics
    mean_VPD = mean(VPD, na.rm = TRUE),
    median_VPD = median(VPD, na.rm = TRUE),
    max_VPD = max(VPD, na.rm = TRUE),
    sd_VPD = sd(VPD, na.rm = TRUE),
    
    # Drought stress metrics
    n_high_vpd_days = sum(VPD > quantile(VPD, 0.75, na.rm = TRUE), na.rm = TRUE),
    pct_high_vpd = round(100 * n_high_vpd_days / n(), 1),
    
    # Soil moisture
    mean_vwc = mean(vwc_avg, na.rm = TRUE),
    median_vwc = median(vwc_avg, na.rm = TRUE),
    
    # Precipitation
    total_precip = sum(precip, na.rm = TRUE),
    mean_precip = mean(precip, na.rm = TRUE),
    
    # Temperature
    mean_soilt = mean(soilt_avg, na.rm = TRUE),
    mean_airt = mean(airt, na.rm = TRUE),
    
    # Data quality
    n_days = n(),
    n_NEE_obs = sum(!is.na(NEE)),
    
    .groups = "drop"
  ) %>%
  arrange(desc(mean_VPD)) %>%
  mutate(
    drought_year = mean_VPD > median(mean_VPD),  # Simple classification
    year_type = ifelse(drought_year, "Drought", "Normal")
  )

print(year_characteristics)
cat("\n")

# Save year characteristics
write.csv(year_characteristics, 
          "/home/claude/year_characteristics.csv", 
          row.names = FALSE)

# Identify the drought year(s)
drought_years <- year_characteristics %>% 
  filter(drought_year) %>% 
  pull(year)

cat("Drought year(s):", paste(drought_years, collapse = ", "), "\n\n")

# ============================================================================
# STEP 2: Create Predictions Across VPD Gradient for Each Year
# ============================================================================

cat("=== STEP 2: Generating Predictions ===\n\n")

# Approach A: Hold all other variables at overall medians
pred_by_year_simple <- expand.grid(
  VPD = seq(min(ec_daily_mngmnt2$VPD, na.rm = TRUE),
            max(ec_daily_mngmnt2$VPD, na.rm = TRUE),
            length.out = 100),
  management = c("conventional", "organic"),
  year = factor(unique(ec_daily_mngmnt2$year)),
  crop_stage_simple = "vegetative",  # Representative stage
  drought_stress = FALSE,
  # Overall medians
  soilt_avg = median(ec_daily_mngmnt2$soilt_avg, na.rm = TRUE),
  vwc_avg = median(ec_daily_mngmnt2$vwc_avg, na.rm = TRUE),
  doy = 200,  # Mid-season
  days_since_tillage = 30
)

# Get predictions
cat("Predicting NEE for simple approach (overall medians)...\n")
preds_simple <- predict(g_resilience, newdata = pred_by_year_simple, se.fit = TRUE)
pred_by_year_simple$NEE <- preds_simple$fit
pred_by_year_simple$se <- preds_simple$se.fit
pred_by_year_simple$lower <- pred_by_year_simple$NEE - 1.96 * pred_by_year_simple$se
pred_by_year_simple$upper <- pred_by_year_simple$NEE + 1.96 * pred_by_year_simple$se

# Approach B: Use year-specific environmental conditions (BETTER)
year_medians <- ec_daily_mngmnt2 %>%
  group_by(year, management) %>%
  summarise(
    median_soilt = median(soilt_avg, na.rm = TRUE),
    median_vwc = median(vwc_avg, na.rm = TRUE),
    median_doy = median(doy, na.rm = TRUE),
    .groups = "drop"
  )

pred_by_year_specific <- expand.grid(
  VPD = seq(min(ec_daily_mngmnt2$VPD, na.rm = TRUE),
            max(ec_daily_mngmnt2$VPD, na.rm = TRUE),
            length.out = 100),
  year = factor(unique(ec_daily_mngmnt2$year)),
  management = c("conventional", "organic")
) %>%
  left_join(year_medians, by = c("year", "management")) %>%
  mutate(
    crop_stage_simple = factor("vegetative"),
    drought_stress = FALSE,
    soilt_avg = median_soilt,
    vwc_avg = median_vwc,
    doy = median_doy,
    days_since_tillage = 30
  )

# Get predictions with year-specific conditions
cat("Predicting NEE with year-specific conditions...\n")
preds_specific <- predict(g_resilience, newdata = pred_by_year_specific, se.fit = TRUE)
pred_by_year_specific$NEE <- preds_specific$fit
pred_by_year_specific$se <- preds_specific$se.fit
pred_by_year_specific$lower <- pred_by_year_specific$NEE - 1.96 * pred_by_year_specific$se
pred_by_year_specific$upper <- pred_by_year_specific$NEE + 1.96 * pred_by_year_specific$se

cat("Predictions complete!\n\n")

# ============================================================================
# STEP 3: Calculate Management Differences for Each Year
# ============================================================================

cat("=== STEP 3: Calculating Management Differences ===\n\n")

# Simple approach
difference_by_year_simple <- pred_by_year_simple %>%
  select(VPD, year, management, NEE, se) %>%
  pivot_wider(names_from = management, 
              values_from = c(NEE, se)) %>%
  mutate(
    difference = NEE_organic - NEE_conventional,
    difference_flipped = -1 * difference,  # Positive = organic better at uptake
    se_diff = sqrt(se_organic^2 + se_conventional^2),
    lower_diff = difference_flipped - 1.96 * se_diff,
    upper_diff = difference_flipped + 1.96 * se_diff
  ) %>%
  left_join(year_characteristics %>% select(year, year_type, mean_VPD), 
            by = "year")

# Year-specific approach
difference_by_year_specific <- pred_by_year_specific %>%
  select(VPD, year, management, NEE, se) %>%
  pivot_wider(names_from = management, 
              values_from = c(NEE, se)) %>%
  mutate(
    difference = NEE_organic - NEE_conventional,
    difference_flipped = -1 * difference,
    se_diff = sqrt(se_organic^2 + se_conventional^2),
    lower_diff = difference_flipped - 1.96 * se_diff,
    upper_diff = difference_flipped + 1.96 * se_diff
  ) %>%
  left_join(year_characteristics %>% select(year, year_type, mean_VPD), 
            by = "year")

# ============================================================================
# STEP 4: Summary Statistics by Year
# ============================================================================

year_effects_summary <- difference_by_year_specific %>%
  group_by(year, year_type, mean_VPD) %>%
  summarise(
    max_org_advantage = max(difference_flipped, na.rm = TRUE),
    min_org_advantage = min(difference_flipped, na.rm = TRUE),
    avg_org_advantage = mean(difference_flipped, na.rm = TRUE),
    vpd_at_max_advantage = VPD[which.max(difference_flipped)],
    vpd_at_min_advantage = VPD[which.min(difference_flipped)],
    range_org_advantage = max_org_advantage - min_org_advantage,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_VPD))

print(year_effects_summary)

write.csv(year_effects_summary, 
          "/home/claude/year_effects_summary.csv", 
          row.names = FALSE)

# ============================================================================
# STEP 5: VISUALIZATIONS
# ============================================================================

cat("\n=== STEP 5: Creating Visualizations ===\n\n")

# PLOT 1: NEE Response Curves by Year (Year-Specific Conditions)
# Shows actual predicted NEE for each management system in each year

p1 <- ggplot(pred_by_year_specific, aes(x = VPD, y = NEE, color = management)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = management),
              alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~ year, ncol = 2) +
  scale_color_manual(values = c("conventional" = "#E69F00", "organic" = "#009E73")) +
  scale_fill_manual(values = c("conventional" = "#E69F00", "organic" = "#009E73")) +
  labs(title = "NEE Response to VPD by Year",
       subtitle = "Year-specific soil conditions | More negative = greater carbon uptake",
       x = "VPD (Pa)",
       y = "Predicted NEE (g C m⁻² d⁻¹)",
       color = "Management",
       fill = "Management") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("/home/claude/Fig_NEE_by_year_VPD.png", p1, 
       width = 10, height = 8, dpi = 300)
cat("Saved: Fig_NEE_by_year_VPD.png\n")

# PLOT 2: Management Difference by Year
# Shows organic advantage/disadvantage across VPD for each year

p2 <- ggplot(difference_by_year_specific, 
             aes(x = VPD, y = difference_flipped)) +
  geom_ribbon(aes(ymin = lower_diff, ymax = upper_diff),
              alpha = 0.3, fill = "steelblue") +
  geom_line(size = 1.2, color = "darkblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ year, ncol = 2) +
  labs(title = "Organic Carbon Uptake Advantage by Year",
       subtitle = "Positive = organic has greater carbon uptake | Shaded = 95% CI",
       x = "VPD (Pa)",
       y = "Organic Advantage\n(sign-reversed ΔNEE, g C m⁻² d⁻¹)") +
  theme_bw(base_size = 12)

ggsave("/home/claude/Fig_management_difference_by_year.png", p2, 
       width = 10, height = 8, dpi = 300)
cat("Saved: Fig_management_difference_by_year.png\n")

# PLOT 3: Year Comparison Overlay - All Years on One Plot
# Shows if drought years have different patterns

# Add year labels with VPD info
difference_by_year_specific <- difference_by_year_specific %>%
  mutate(year_label = paste0(year, " (", year_type, ")"))

p3 <- ggplot(difference_by_year_specific, 
             aes(x = VPD, y = difference_flipped, 
                 color = year_label, group = year)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_viridis_d(option = "plasma", end = 0.9) +
  labs(title = "Does Organic Advantage Vary by Year?",
       subtitle = "Comparing drought vs. normal years",
       x = "VPD (Pa)",
       y = "Organic Carbon Uptake Advantage (g C m⁻² d⁻¹)",
       color = "Year") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

ggsave("/home/claude/Fig_year_comparison_overlay.png", p3, 
       width = 10, height = 6, dpi = 300)
cat("Saved: Fig_year_comparison_overlay.png\n")

# PLOT 4: Comparison of Drought vs. Normal Years (Averaged)
# Groups years into drought vs. normal and shows average pattern

avg_by_year_type <- difference_by_year_specific %>%
  group_by(VPD, year_type) %>%
  summarise(
    mean_diff = mean(difference_flipped, na.rm = TRUE),
    se_diff = sd(difference_flipped, na.rm = TRUE) / sqrt(n()),
    lower = mean_diff - 1.96 * se_diff,
    upper = mean_diff + 1.96 * se_diff,
    .groups = "drop"
  )

p4 <- ggplot(avg_by_year_type, aes(x = VPD, y = mean_diff, color = year_type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = year_type),
              alpha = 0.2, color = NA) +
  geom_line(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Drought" = "#D55E00", "Normal" = "#0072B2")) +
  scale_fill_manual(values = c("Drought" = "#D55E00", "Normal" = "#0072B2")) +
  labs(title = "Organic Advantage: Drought vs. Normal Years",
       subtitle = "Average pattern across year types",
       x = "VPD (Pa)",
       y = "Organic Carbon Uptake Advantage (g C m⁻² d⁻¹)",
       color = "Year Type",
       fill = "Year Type") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("/home/claude/Fig_drought_vs_normal_years.png", p4, 
       width = 10, height = 6, dpi = 300)
cat("Saved: Fig_drought_vs_normal_years.png\n")

# PLOT 5: Add Observed Data Points
# Show where each year actually fell on the VPD gradient

observed_summary <- ec_daily_mngmnt2 %>%
  filter(!is.na(NEE)) %>%
  group_by(year, management) %>%
  summarise(
    mean_VPD = mean(VPD, na.rm = TRUE),
    mean_NEE = mean(NEE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(year_characteristics %>% select(year, year_type), by = "year")

p5 <- ggplot(pred_by_year_specific, aes(x = VPD, y = NEE, color = management)) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(data = observed_summary, 
             aes(x = mean_VPD, y = mean_NEE, shape = year_type),
             size = 5, stroke = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~ year, ncol = 2) +
  scale_color_manual(values = c("conventional" = "#E69F00", "organic" = "#009E73")) +
  scale_shape_manual(values = c("Drought" = 17, "Normal" = 16)) +
  labs(title = "NEE Response to VPD with Observed Year Means",
       subtitle = "Lines = model predictions | Points = actual year mean",
       x = "VPD (Pa)",
       y = "NEE (g C m⁻² d⁻¹)",
       color = "Management",
       shape = "Year Type") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("/home/claude/Fig_NEE_by_year_with_observed.png", p5, 
       width = 10, height = 8, dpi = 300)
cat("Saved: Fig_NEE_by_year_with_observed.png\n")

# PLOT 6: Summary Figure - Key Metrics by Year
# Bar plot showing max advantage, min advantage, etc.

p6 <- year_effects_summary %>%
  select(year, year_type, max_org_advantage, min_org_advantage, avg_org_advantage) %>%
  pivot_longer(cols = c(max_org_advantage, min_org_advantage, avg_org_advantage),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, 
                         levels = c("max_org_advantage", "avg_org_advantage", "min_org_advantage"),
                         labels = c("Max Advantage", "Average", "Max Disadvantage"))) %>%
  ggplot(aes(x = year, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ year_type, scales = "free_x") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Summary of Organic Advantage Across VPD Range by Year",
       x = "Year",
       y = "Organic Advantage (g C m⁻² d⁻¹)",
       fill = "Metric") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("/home/claude/Fig_year_summary_metrics.png", p6, 
       width = 10, height = 6, dpi = 300)
cat("Saved: Fig_year_summary_metrics.png\n")

# ============================================================================
# STEP 6: Statistical Tests
# ============================================================================

cat("\n=== STEP 6: Statistical Comparisons ===\n\n")

# Test: Is organic advantage different in drought vs. normal years?
# At high VPD (>75th percentile)

high_vpd_threshold <- quantile(ec_daily_mngmnt2$VPD, 0.75, na.rm = TRUE)

high_vpd_comparison <- difference_by_year_specific %>%
  filter(VPD > high_vpd_threshold) %>%
  group_by(year, year_type) %>%
  summarise(
    mean_advantage = mean(difference_flipped, na.rm = TRUE),
    .groups = "drop"
  )

cat("Organic advantage at high VPD by year:\n")
print(high_vpd_comparison)

# Simple t-test
if(length(unique(high_vpd_comparison$year_type)) == 2) {
  t_result <- t.test(mean_advantage ~ year_type, data = high_vpd_comparison)
  cat("\nT-test: Drought vs Normal years at high VPD\n")
  cat("Mean advantage in Drought years:", 
      mean(high_vpd_comparison$mean_advantage[high_vpd_comparison$year_type == "Drought"]), "\n")
  cat("Mean advantage in Normal years:", 
      mean(high_vpd_comparison$mean_advantage[high_vpd_comparison$year_type == "Normal"]), "\n")
  cat("p-value:", t_result$p.value, "\n")
}

# ============================================================================
# STEP 7: Create Combined Summary Report
# ============================================================================

cat("\n=== STEP 7: Creating Summary Report ===\n\n")

# Combine key findings
summary_report <- list(
  year_characteristics = year_characteristics,
  year_effects = year_effects_summary,
  high_vpd_comparison = high_vpd_comparison
)

# Save as RDS for later use
saveRDS(summary_report, "/home/claude/resilience_analysis_summary.rds")

# Create text summary
sink("/home/claude/resilience_analysis_summary.txt")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("YEAR-SPECIFIC RESILIENCE ANALYSIS SUMMARY\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

cat("YEAR CHARACTERISTICS:\n")
cat("-" %>% rep(70) %>% paste(collapse = ""), "\n")
print(year_characteristics)

cat("\n\nYEAR EFFECTS ON ORGANIC ADVANTAGE:\n")
cat("-" %>% rep(70) %>% paste(collapse = ""), "\n")
print(year_effects_summary)

cat("\n\nCOMPARISON AT HIGH VPD:\n")
cat("-" %>% rep(70) %>% paste(collapse = ""), "\n")
print(high_vpd_comparison)

cat("\n\nKEY FINDINGS:\n")
cat("-" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("1. Drought year(s):", paste(drought_years, collapse = ", "), "\n")
cat("2. Number of years analyzed:", nrow(year_characteristics), "\n")
cat("3. VPD range:", round(min(ec_daily_mngmnt2$VPD, na.rm = TRUE), 1), 
    "to", round(max(ec_daily_mngmnt2$VPD, na.rm = TRUE), 1), "Pa\n")

sink()

cat("\n=== ANALYSIS COMPLETE ===\n\n")
cat("Files created:\n")
cat("  1. year_characteristics.csv\n")
cat("  2. year_effects_summary.csv\n")
cat("  3. resilience_analysis_summary.rds\n")
cat("  4. resilience_analysis_summary.txt\n")
cat("  5. Fig_NEE_by_year_VPD.png\n")
cat("  6. Fig_management_difference_by_year.png\n")
cat("  7. Fig_year_comparison_overlay.png\n")
cat("  8. Fig_drought_vs_normal_years.png\n")
cat("  9. Fig_NEE_by_year_with_observed.png\n")
cat(" 10. Fig_year_summary_metrics.png\n\n")

cat("Next steps:\n")
cat("  - Examine the figures to see if drought years show different patterns\n")
cat("  - Repeat analysis for GPP and Reco models\n")
cat("  - Consider adding crop stage as another facet\n")
cat("  - Test specific VPD thresholds if you see interesting patterns\n")
