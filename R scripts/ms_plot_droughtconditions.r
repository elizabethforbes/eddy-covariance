library(tidyverse)
library(here)

# plotting a time series of drought conditions

library(patchwork)

# shared x-axis theme -- suppress x-axis labels on top two panels
theme_top <- theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   panel.grid.minor = element_blank())
theme_bottom <- theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      panel.grid.minor = element_blank())

# summarize to calculate general average across the empirical data:
droughtconditions <- ec_daily_mngmnt3 %>% 
  filter(year %in% c(2018, 2019, 2020)) %>%
  group_by(date) %>%
  # average across fields
  summarize(VPD_mean = mean(VPD, na.rm = TRUE), 
            airt_mean = mean(air_temperature, na.rm = TRUE),
            vwc_mean = mean(vwc_avg_era5),
            .groups = "drop")

# identify mean summertime air temp:

  
# Panel 1: air temperature
p_temp <- droughtconditions %>%
  ggplot(aes(x = date, y = airt_mean)) +
  geom_line(colour = "grey70", linewidth = 0.3) +
  geom_smooth(method = "loess", span = 0.1,
              colour = "steelblue", se = FALSE, linewidth = 0.8) +
  # geom_hline(yintercept = 0, linetype = "dashed",
             # colour = "grey40", linewidth = 0.4) +
  # average monthly temp over the study period: 
  # https://www.ncei.noaa.gov/access/monitoring/climate-at-a-glance/county/time-series/NY-021/tavg/1/0/2018-2020
  # versus historical 'summer' baseline: https://nyaspubs.onlinelibrary.wiley.com/doi/full/10.1111/nyas.15240
  # geom_hline(yintercept = 19, linetype = "dashed",
             # colour = "darkblue", linewidth = 0.4) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  labs(y = "Air temperature (°C)") +
  theme_bw(base_size = 11) + theme_top

# Panel 2: VPD
p_vpd <- droughtconditions %>%
  ggplot(aes(x = date, y = VPD_mean/1000)) +
  geom_line(colour = "grey70", linewidth = 0.3) +
  geom_smooth(method = "loess", span = 0.1,
              colour = "darkred", se = FALSE, linewidth = 0.8) +
  # geom_hline(yintercept = 1.0, linetype = "dashed",
             # colour = "tan3", linewidth = 0.4) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  labs(y = "VPD (kPa)") +
  theme_bw(base_size = 11) + theme_top

# Panel 3: ERA5 VWC
p_vwc <- droughtconditions %>%
  ggplot(aes(x = date, y = vwc_mean)) +
  geom_line(colour = "grey70", linewidth = 0.3) +
  geom_smooth(method = "loess", span = 0.1,
              colour = "sienna", se = FALSE, linewidth = 0.8) +
  # https://www.nature.com/articles/s41467-024-49244-7 referencing this global average for humid environments
  # geom_hline(yintercept = 0.26, linetype = "dashed",
             # colour = "orange", linewidth = 0.4)+
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  labs(y = "Soil moisture\n(ERA5 VWC)") +
  theme_bw(base_size = 11) + theme_top

# panel 4: NYS drought conditions data
ny_drought <- read_csv(here("eddy_covariance_fluxdata", "USDM-columbia-county-ny.csv"))
ny_drought <- ny_drought %>% 
  mutate(date_start = as.Date(ValidStart),
         date_end = as.Date(ValidEnd)) %>%
  filter(year(date_start) %in% c(2018, 2019, 2020)) %>%
  # pivot:
  select(date_start, date_end, None, D0, D1, D2, D3, D4) %>% 
  rowwise() %>% 
  mutate(date = list(seq(date_start, date_end, by = "day"))) %>% 
  unnest(date) %>% 
  ungroup() %>% 
  select(date, None, D0, D1, D2, D3, D4) %>%
  # filter to match my data:
  filter(date >= as.Date("2018-08-13"),
         date <= as.Date("2020-12-31")) %>% 
  pivot_longer(-date, names_to = "category", values_to = "pct") %>% 
  mutate(category = factor(category,
                           levels = c("D4", "D3", "D2", "D1", "D0", "None"))) %>% 
  # select only those categories that appear in the data: D1, and D0
  filter(category %in% c("D0", "D1")) %>% 
  # relevel so that D1 stacks on top of D0 in the plot
  mutate(category = factor(category, levels = c("D1", "D0")))

# define so we can trim down the county data to the dates of our EC data:
ec_date_range <- range(droughtconditions$date)

ny_drought_plot <- ny_drought %>%
  filter(category %in% c("D0", "D1"),
         date >= ec_date_range[1],
         date <= ec_date_range[2]) %>%
  pivot_wider(names_from = category,
              values_from = pct) %>% 
  # replace NAs with 0
  mutate(across(c(D0, D1), ~replace_na(., 0)))

p_drought <- ny_drought_plot %>%
  ggplot(aes(x = date)) +
  geom_col(aes(y = D0), fill = "#FFFF00", width = 1) +
  geom_col(aes(y = D1), fill = "#FCD37F", width = 1) +
  # manual legend as annotated colored rectangles
  annotate("rect", 
           xmin = as.Date("2018-09-01"), xmax = as.Date("2018-11-15"),
           ymin = 78, ymax = 90, fill = "#FFFF00", color = "grey50", linewidth = 0.3) +
  annotate("rect",
           xmin = as.Date("2018-09-01"), xmax = as.Date("2018-11-15"),
           ymin = 63, ymax = 75, fill = "#FCD37F", color = "grey50", linewidth = 0.3) +
  annotate("text",
           x = as.Date("2018-11-30"), y = 84,
           label = "D0: Abnormally dry", hjust = 0, size = 3, color = "grey20") +
  annotate("text",
           x = as.Date("2018-11-30"), y = 69,
           label = "D1: Moderate drought", hjust = 0, size = 3, color = "grey20") +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  scale_y_continuous(limits = c(0, 100),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Area in drought\n(% of county)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position = "none")

# Stack with patchwork
p_temp / p_vpd / p_vwc / p_drought +
  plot_layout(heights = c(1, 1, 1, 1)) +
  plot_annotation(
    caption = "Air temperature, VPD averaged across both tower sites. ERA5 VWC = volumetric water content."
  )

