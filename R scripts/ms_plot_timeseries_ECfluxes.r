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

# ============================================================================
# plot time series of daily flux estimates:
# ============================================================================

stage_colors <- c("fallow" = "tan4", "early" = "olivedrab2", "vegetative" = "green4", "reproductive" = "orchid3", 
                  "grain_fill" = "#D94801", "mature" = "goldenrod1", "dormant" = "#7B6888")


stage_order <- c("fallow", "early", "vegetative", "reproductive", 
                 "grain_fill", "mature", "dormant")

plot_timeseries_with_stages <- function(data, flux_var = "NEE", year = NULL) {
  
  if (!is.null(year)) {
    data <- data %>% filter(.data$year == !!year)
  }
  
  data$crop_stage_simple <- factor(data$crop_stage_simple, levels = stage_order)
  
  # Calculate ymin and ymax for flux variable, to include background temp data
  ymin <- min(data[[flux_var]], na.rm = TRUE)
  ymax <- max(data[[flux_var]], na.rm = TRUE)
  
  ggplot(data, aes(x = date, y = .data[[flux_var]])) +
    geom_point(aes(colour = crop_stage_simple), alpha = 0.6, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
    
    # add background tiles colored by air temperature:
    geom_tile(aes(y = (ymin + ymax)/2, height = ymax - ymin, fill = data1$airt_era5), alpha = 0.15) +
    
    # geom_smooth(colour = "darkgreen", 
    #             method = gam, 
    #             formula = y~s(x, bs = "cs"),
    #             linewidth = 0.8, se = TRUE) +
    
    scale_colour_manual(values = stage_colors, 
                        name = "Crop stage",
                        na.translate = FALSE) +  # This hides NA from the legend
    
    # color for tiles:
    scale_fill_viridis_c(name = "Air Temp (°C)",
                         alpha = 0.2) +  # continuous color scale for background
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
    facet_wrap(~ management, ncol = 1, labeller = labeller(
      management = c("conventional" = "Conventional", "organic" = "Organic")
    ),
    drop = TRUE) +
    labs(
      x        = NULL,
      y        = bquote(.(flux_var) ~ "(g C m"^-2 ~ "d"^-1 ~ ")")
    ) +
    theme_flux() +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.position  = "bottom",
      legend.box = "vertical", # stacks legends vertically
      legend.box.just = "center",
      legend.direction = "horizontal",
      strip.text       = element_text(face = "bold")
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "top", 
        title.hjust = 0.5, 
        order = 1,
      ),
      colour = guide_legend(
        nrow = 1, 
        title.position = "top", 
        title.hjust = 0.5, 
        order = 2,
        override.aes = list(size = 3)
      ))
}

nt <- plot_timeseries_with_stages(data1, flux_var = "NEE_gapfilled")
gt <- plot_timeseries_with_stages(data1, flux_var = "GPP_gapfilled")
rt <- plot_timeseries_with_stages(data1, flux_var = "Reco_gapfilled")

# plot all three together:
nt1 <- plot_timeseries_with_stages(data1, flux_var = "NEE_gapfilled")+theme(legend.position = "none")
gt1 <- plot_timeseries_with_stages(data1, flux_var = "GPP_gapfilled")+theme(legend.position = "none")
nt1/gt1/rt
# save in size 1200 x 1200 to get labels, etc. correct