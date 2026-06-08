# # Color palettes for crop stages
# stage_colors <- c(
#   "fallow" = "#D3D3D3",        # Light gray
#   "recovery" = "#FF6B6B",      # Red (high respiration)
#   "early" = "#95E1D3",         # Light teal
#   "regrowth" = "#4ECDC4",      # Teal
#   "vegetative" = "#45B7D1",    # Blue
#   "reproductive" = "#F38181",   # Pink
#   "grain_fill" = "#FFA07A",    # Light coral
#   "mature" = "#DAA520",        # Goldenrod
#   "dormant" = "#9B9B9B"        # Gray
# )

# Color palette for management
management_colors <- c(
  "conventional" = "tomato",  # Red
  "organic" = "#E69F00"        # Green
)

plot_management_comparison <- function(data, flux_var = "NEE") {
  
  stage_order <- c("fallow", "recovery", "early", "regrowth", "vegetative",
                   "reproductive", "grain_fill", "mature", "dormant")
  data$crop_stage_simple <- factor(
    data$crop_stage_simple,
    levels = intersect(stage_order, unique(data$crop_stage_simple))
  ) 
  
  # Filter raw data to groups with n >= 5 (so points match the bars shown)
  data_filtered <- data %>%
    dplyr::group_by(management, crop_stage_simple) %>%
    dplyr::filter(sum(!is.na(.data[[flux_var]])) >= 5) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(crop_stage_simple = droplevels(crop_stage_simple))   # <-- new
  
  # Summary for bars and error bars
  summary_data <- data_filtered %>%
    dplyr::group_by(management, crop_stage_simple) %>%
    dplyr::summarise(
      mean_flux = mean(.data[[flux_var]], na.rm = TRUE),
      se_flux   = sd(.data[[flux_var]], na.rm = TRUE) / sqrt(dplyr::n()),
      n         = dplyr::n(),
      .groups   = "drop"
    )
  
  ggplot() +
    # Bars first (so points overlay them)
    geom_col(
      data = summary_data,
      aes(x = crop_stage_simple, y = mean_flux, fill = management),
      position = position_dodge(0.9), alpha = 0.6
    ) +
    # Raw points — note explicit data argument and jitterdodge position
    geom_jitter(
      data = data_filtered,
      aes(x = crop_stage_simple, y = .data[[flux_var]],
          color = management, group = management),
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
      alpha = 0.25, size = 0.8, show.legend = FALSE
    ) +
    # Error bars on top of everything
    geom_errorbar(
      data = summary_data,
      aes(x = crop_stage_simple,
          ymin = mean_flux - se_flux,
          ymax = mean_flux + se_flux,
          group = management),
      position = position_dodge(0.9), width = 0.2
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = management_colors,  name = "Management") +
    scale_color_manual(values = management_colors, guide = "none") +
    labs(
      title    = paste0(flux_var, " by Management and Crop Stage"),
      subtitle = "mean ± SE, points = daily observations (min. n=5)",
      x        = element_blank(),
      y        = bquote(.(flux_var) ~ "(µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")")
    ) +
    theme_flux() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
}

