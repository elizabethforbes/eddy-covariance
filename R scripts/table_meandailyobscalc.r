library(tidyverse)

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