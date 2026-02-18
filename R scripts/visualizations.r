library(tidyverse)

dormancy_break_date <- as.Date("2019-02-17")

df %>% 
  # filter(management == "organic") %>%
  # filter(date < "2019-06-15") %>% 
  ggplot(aes(x = date)) +
  geom_line(aes(y = ndvi), color = "forestgreen", size = 1) +
  geom_line(aes(y = soilt_avg / 15), color = "steelblue", size = 1) +
  # geom_vline(xintercept = dormancy_break_date, color = "purple", linetype = "dotted", size = 1) + #indicating break in dormancy for 2019
  geom_hline(yintercept = 0 / 15, linetype = "dashed", color = "blue", size = 0.7) + # indicates OC == 32F
  geom_hline(yintercept = 4 / 15, linetype = "dashed", color = "red", size = 0.7) + # indicates 4C == 39F
  scale_y_continuous(
    name = "NDVI",
    limits = c(min(df$ndvi), max(df$ndvi)),
    sec.axis = sec_axis(~ . * 15, name = "Soil Temperature (°C)")
  ) +
  # scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  theme_minimal() +
  labs(
    title = "Winter Rye NDVI and Soil Temperature",
    x = "",
    # subtitle = "Purple dotted line = Dormancy Break (~Feb 17, 2019)"
  ) +
  theme(
    axis.title.y.left = element_text(color = "forestgreen", size = 12),
    axis.title.y.right = element_text(color = "steelblue", size = 12),
    axis.text.y.left = element_text(color = "forestgreen"),
    axis.text.y.right = element_text(color = "steelblue"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )


#### plotting interactions, NEE or Reco
library(mgcv)
library(ggplot2)
library(dplyr)
library(mgcViz)

# Assuming your fitted model is called gam_model

# Create new data frame for prediction
newdata <- expand.grid(
  management = c("conventional", "organic"),
  crop_stage_simple = levels(ec_daily_mngmnt2$crop_stage_simple),
  year = mean(ec_daily_mngmnt2$year),          # hold other variables constant at mean or typical values
  doy = median(ec_daily_mngmnt2$doy),
  soilt_avg = mean(ec_daily_mngmnt2$soilt_avg, na.rm = TRUE),
  vwc_avg = mean(ec_daily_mngmnt2$vwc_avg, na.rm = TRUE),
  VPD = mean(ec_daily_mngmnt2$VPD, na.rm = TRUE),
  days_since_tillage = mean(ec_daily_mngmnt2$days_since_tillage, na.rm = TRUE)
)

# Predict fitted values with standard errors
# pred <- predict(g3.3, newdata, se.fit = TRUE, type = "response")
pred <- predict(r1.3, newdata, se.fit = TRUE, type = "response")

newdata <- newdata %>%
  mutate(
    fit = pred$fit,
    se = pred$se.fit,
    lower = fit - 1.96 * se,
    upper = fit + 1.96 * se
  )

# Plot predicted NEE by crop stage and management
ggplot(newdata, aes(x = crop_stage_simple, y = fit, color = management, group = management)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_line(position = position_dodge(width = 0.5), size = 1) +
  labs(
    # title = "Predicted NEE by Crop Stage and Management Type",
    title = "Predicted Ecosystem Respiration by Crop Stage and Management Type",
    x = "Crop Stage",
    # y = "Predicted NEE",
    y = "Predicted R(eco)",
    color = "Management"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### plotting smooth terms, NEE

# Convert to mgcViz object
viz <- getViz(g3.3)

# Plot smooth of day of year (seasonality), third smooth in the model
plot(sm(viz, 3)) + l_fitLine() + l_ciLine() + l_rug() + ggtitle("Seasonal effect (Day of Year)")

# Plot smooth of soil temp by management
plot(sm(viz, 6)) + l_fitLine() + l_ciLine() + l_rug() + ggtitle("Soil Temperature effect -- Conventional")
plot(sm(viz, 7)) + l_fitLine() + l_ciLine() + l_rug() + ggtitle("Soil Temperature effect -- Organic")

# Plot tensor interaction smooth of soil temp and moisture by management
plot(sm(viz, 4)) + l_fitContour() + ggtitle("Soil Temp & Moisture interaction -- Conventional")
plot(sm(viz, 5)) + l_fitContour() + ggtitle("Soil Temp & Moisture interaction -- Organic")

# Plot smooth of days since tillage by management
plot(sm(viz, 12)) + l_fitLine() + l_ciLine() + l_rug() + ggtitle("Days since tillage effect -- Conventional")
plot(sm(viz, 13)) + l_fitLine() + l_ciLine() + l_rug() + ggtitle("Days since tillage effect -- Organic")

#### plotting effect size for year-to-year variation, NEE
# Predict terms with SE
terms_pred <- predict(g3.3, 
                      newdata = ec_daily_mngmnt2,
                      type = "terms", se.fit = TRUE)

term_names <- colnames(terms_pred$fit)

# Indices for year random effects by management
idx_conventional <- grep("s\\(year\\):as.factor\\(management\\)conventional", term_names)
idx_organic <- grep("s\\(year\\):as.factor\\(management\\)organic", term_names)

# Add predictions to original data
pred_df <- ec_daily_mngmnt2 %>%
  mutate(
    re_conventional = terms_pred$fit[, idx_conventional],
    se_conventional = terms_pred$se.fit[, idx_conventional],
    re_organic = terms_pred$fit[, idx_organic],
    se_organic = terms_pred$se.fit[, idx_organic]
  )

# Summarize to get one value per year and management
df_conventional <- pred_df %>%
  filter(management == "conventional") %>%
  group_by(year) %>%
  summarize(
    random_effect = mean(re_conventional, na.rm = TRUE),
    se = mean(se_conventional, na.rm = TRUE)
  ) %>%
  mutate(
    management = "conventional",
    lower = random_effect - 1.96 * se,
    upper = random_effect + 1.96 * se
  )

df_organic <- pred_df %>%
  filter(management == "organic") %>%
  group_by(year) %>%
  summarize(
    random_effect = mean(re_organic, na.rm = TRUE),
    se = mean(se_organic, na.rm = TRUE)
  ) %>%
  mutate(
    management = "organic",
    lower = random_effect - 1.96 * se,
    upper = random_effect + 1.96 * se
  )

# Combine
df_year_re <- bind_rows(df_conventional, df_organic)
df_year_re <- df_year_re %>% drop_na()

# Plot
ggplot(df_year_re, aes(x = year, y = random_effect, color = management, group = management)) +
  geom_line() +
  geom_point() +
  # geom_ribbon(aes(ymin = lower, ymax = upper, fill = management), alpha = 0.2, color = NA) +
  labs(
    title = "Year-specific Random Effects on NEE by Management",
    x = "Year",
    y = "Estimated Random Effect (NEE deviation)",
    color = "Management",
    fill = "Management"
  ) +
  theme_minimal()

