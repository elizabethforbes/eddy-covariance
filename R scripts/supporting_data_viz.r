# ============================================================================
# VISUALIZATION of SUPPORTING DATA, eddy covariance study
# ============================================================================
# 
# Purpose: Publication-ready visualizations of ERA5-modeled soil moisture,
# empirical SOM across fields
#
# Requirements:
#   - ggplot2, dplyr, tidyr, viridis, patchwork, scales
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)
library(scales)
library(readxl)
library(here)
library(DHARMa)
library(lmtest)
library(nlme)
library(mgcv)
library(lubridate)
library(patchwork)

my_custom_theme <- function() {
  theme_minimal(base_size = 12, base_family = "Helvetica") %+replace%
    theme(
      panel.grid.major = element_line(color = "gray85", size = 0.15),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(color = "gray20"),
      axis.ticks = element_line(color = "gray40"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.background = element_blank(),
      strip.background = element_rect(fill = "gray90", color = NA),
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
      plot.background = element_blank()
    )
}


# ============================================================================
# Visualize soil moisture collected using handheld chambers for surveys
# ============================================================================

chamber_daily_mngmnt <- chamber_daily_mngmnt %>% 
  mutate(year = as.factor(year),
         month_num = month(date),
         month = factor(month_num, levels = 1:12, labels = month.abb),
         management = as.factor(management))

ggplot(chamber_daily_mngmnt, aes(y = soilm_perc_mean, color = management, fill = management)) +
  geom_boxplot(aes(x = month, group = interaction(month, management)),
               position = position_dodge(width = 0.8),
               alpha = 0.5, outlier.shape = NA, color = "black") +
  geom_jitter(aes(x = month, group = interaction(month, management)),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.5, size = 1.5) +
  geom_smooth(aes(x = month, group = management), method = "loess", 
              se = FALSE, size = 1) +
  # add horizontal line at 5% VWC (low end of permanent wilting point for crops in sandy soils)
  geom_hline(yintercept = 5, linetype = "dashed", color = "darkred", size = 1)+
  scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) + 
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  scale_x_discrete(drop = FALSE, labels = month.abb) +
  labs(x = NULL, y = "Soil moisture (%)", color = "Management", fill = "Management") +
  my_custom_theme()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~year)

# ============================================================================
# Visualize VPD collected from EC towers
# ============================================================================

dat <- ec_daily_mngmnt %>% 
  mutate(year = as.factor(year),
         # month_num = month(date),
         month = factor(month, levels = 1:12, labels = month.abb),
         management = as.factor(management))

ggplot(dat, aes(y = VPD/1000))+
                # , color = management, fill = management)) +
  # horizontal band indicating "ideal" VPD range for plant growth:
  geom_rect(data = dat[1,],
    aes(
      xmin = -Inf,  # Extend to the left edge
      xmax = Inf,   # Extend to the right edge
      ymin = 0.8,     # Lower y-bound of the band
      ymax = 1.2      # Upper y-bound of the band
    ),
    fill = "olivedrab",
    alpha = 0.25  # Transparency
  ) +
  geom_boxplot(aes(x = month),
                   # group = interaction(month, management)),
               position = position_dodge(width = 0.8),
               alpha = 0.5, outlier.shape = NA, color = "black") +
  geom_jitter(aes(x = month), 
                  # group = interaction(month, management)),
              # geom_point(aes(color = above_avg), size = 2, alpha = 0.7)
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.5, size = 1.5) +
  geom_smooth(aes(x = month), 
                  # group = management), 
              method = "loess", 
              se = FALSE, size = 1) +
  # scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) + 
  # scale_color_manual(values = ) +
  scale_x_discrete(drop = FALSE, labels = month.abb) +
  labs(x = NULL, y = "VPD (kPa)", color = "Management", fill = "Management") +
  my_custom_theme()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # facet_wrap(~year)

# ============================================================================
# table of crop yields per acre across 2018-2020
# ============================================================================

library(flextable)

yield_data <- matrix(c("corn (2018)", "91.4", "194",
                       "soy (2019)", "26.7", "56",
                       "barley (2020)", "54.4", "43"), ncol = 3, byrow = TRUE)
colnames(yield_data) <- c("crop", "organic", "conventional")
# rownames(yield_data) <- c("corn (2018)", "soybeans (2019)", "barley (2020)")
yield_data <- as.data.frame(yield_data)

# add converted units of carbon data
# yield_data <- yield_data %>% 
  # mutate(Freq2 = "453", "167", "84", "213", "80", "106")

y_table <- flextable(yield_data) %>%
  # Header labels
  set_header_labels(
    crop = "Crop",
    organic      = "Organic Yield (bu/ac)",
    conventional       = "Conventional Yield (bu/ac)") %>% 
  
  # # Column widths (inches, for Word)
  width(j = "crop",       width = 1.2) %>%
  width(j = "organic",        width = 1.5) %>%
  width(j = "conventional",        width = 1.5) %>%
  
  
  # light shading for each management:
  bg(i = ~ grepl("conventional", organic), bg = "#ffe5e1") %>% 
  bg(i = ~ grepl("organic", conventional), bg = "#fff0ce") %>% 
  
  # Borders
  border_outer(part = "all",    border = officer::fp_border(width = 1.5)) %>%
  border_inner_h(part = "body", border = officer::fp_border(width = 0.4,
                                                            color = "grey70")) %>%
  hline(part = "header",        border = officer::fp_border(width = 1.5))
  
y_table

