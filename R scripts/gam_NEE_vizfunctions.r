# ============================================================================
# GAM VISUALIZATION SCRIPTS FOR AGRICULTURAL CARBON FLUX ANALYSIS
# ============================================================================
# 
# Purpose: Publication-ready visualizations for GAM model outcomes analyzing
#          carbon fluxes (NEE, GPP, Reco) across crop phenology stages
#
# Author: Generated for eddy covariance flux analysis
# Date: March 2026
#
# Requirements:
#   - ggplot2, dplyr, tidyr, viridis, patchwork, scales
#   - GAM model outputs (mgcv package)
#   - Daily EC flux data with crop_stage_simple categories
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)
library(scales)
library(mgcv)
library(gratia)
library(stringr)

# ============================================================================
# THEME SETTINGS - Publication Quality
# ============================================================================

# Custom theme for all plots
theme_flux <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      # Text
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = "gray30", hjust = 0),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      
      # Grid
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      
      # Background
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "gray80"),
      legend.background = element_rect(fill = "white", color = NA),
      
      # Margins
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Color palettes for crop stages
stage_colors <- c(
  "fallow" = "#D3D3D3",        # Light gray
  "recovery" = "#FF6B6B",      # Red (high respiration)
  "early" = "#95E1D3",         # Light teal
  "regrowth" = "#4ECDC4",      # Teal
  "vegetative" = "#45B7D1",    # Blue
  "reproductive" = "#F38181",   # Pink
  "grain_fill" = "#FFA07A",    # Light coral
  "mature" = "#DAA520",        # Goldenrod
  "dormant" = "#9B9B9B"        # Gray
)

# Color palette for management
management_colors <- c(
  "conventional" = "tomato",  # Red
  "organic" = "#E69F00"        # Green
)

# ============================================================================
# 1. OBSERVED VS PREDICTED - Basic GAM Performance
# ============================================================================

plot_gam_observed_vs_predicted <- function(model, data, flux_var = "NEE", 
                                           title = NULL) {
  #' Plot observed vs predicted values with 1:1 line
  #' 
  #' @param model GAM model object
  #' @param data Data frame with observed values (not used directly, for compatibility)
  #' @param flux_var Character, name of flux variable (NEE, GPP, Reco)
  #' @param title Character, plot title (auto-generated if NULL)
  
  # Get model data (only rows actually used in fitting)
  model_data <- model$model
  
  # Get predictions and observed values
  model_data$predicted <- predict(model, type = "response")
  model_data$observed <- model_data[[1]]  # Response variable is first column
  
  # Calculate R² and RMSE
  r_squared <- cor(model_data$observed, model_data$predicted, use = "complete.obs")^2
  rmse <- sqrt(mean((model_data$observed - model_data$predicted)^2, na.rm = TRUE))
  
  # Auto title
  if (is.null(title)) {
    title <- paste0("GAM Performance: ", flux_var)
  }
  
  # Plot
  p <- ggplot(model_data, aes(x = observed, y = predicted)) +
    # Hexbin for density
    geom_hex(bins = 50, alpha = 0.7) +
    scale_fill_viridis_c(option = "plasma", name = "Count") +
    
    # 1:1 line
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "red", linewidth = 1) +
    
    # Smooth
    geom_smooth(method = "lm", color = "blue", linewidth = 1, alpha = 0.2) +
    
    # Labels and theme
    labs(
      title = title,
      subtitle = sprintf("R² = %.3f | RMSE = %.2f | n = %d", r_squared, rmse, nrow(model_data)),
      x = bquote("Observed" ~ .(flux_var) ~ "(µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")"),
      y = bquote("Predicted" ~ .(flux_var) ~ "(µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")")
    ) +
    theme_flux() +
    coord_equal() +
    theme(legend.position = "right")
  
  return(p)
}

# ============================================================================
# 2. RESIDUALS BY CROP STAGE - Check for bias
# ============================================================================

plot_gam_residuals_by_stage <- function(model, data, flux_var = "NEE") {
  #' Plot residuals by crop stage to check for systematic bias
  #' 
  #' @param model GAM model object
  #' @param data Data frame with crop_stage_simple
  #' @param flux_var Character, name of flux variable
  
  # Get model data (only rows actually used in fitting)
  # This handles NA values that were automatically excluded
  model_data <- model$model
  model_data$residuals <- residuals(model, type = "response")
  
  # If crop_stage_simple not in model_data, get it from original data
  if (!"crop_stage_simple" %in% names(model_data)) {
    # Match rows using row names
    model_rows <- as.numeric(rownames(model_data))
    model_data$crop_stage_simple <- data$crop_stage_simple[model_rows]
  }
  
  # Order stages
  stage_order <- c("fallow", "recovery", "early", "regrowth", "vegetative",
                   "reproductive", "grain_fill", "mature", "dormant")
  model_data$crop_stage_simple <- factor(model_data$crop_stage_simple, 
                                         levels = stage_order)
  
  # Plot
  p <- ggplot(model_data, aes(x = crop_stage_simple, y = residuals, 
                              fill = crop_stage_simple)) +
    # Violin + boxplot
    geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    
    # Zero line
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    
    # Colors
    scale_fill_manual(values = stage_colors, guide = "none") +
    
    # Labels
    labs(
      title = paste0("Model Residuals by Crop Stage: ", flux_var),
      subtitle = "Check for systematic bias across phenology stages",
      x = "Crop Stage",
      y = bquote("Residuals (µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")")
    ) +
    theme_flux() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# ============================================================================
# 3. SMOOTH EFFECTS - Visualize GAM smooth terms
# ============================================================================

plot_gam_smooth_effects <- function(model, term, data = NULL, 
                                    rug = TRUE, title = NULL) {
  #' Plot GAM smooth term effects
  #' 
  #' @param model GAM model object
  #' @param term Character, name of smooth term (e.g., "s(air_temperature)")
  #' @param data Data frame (optional, for rug plot)
  #' @param rug Logical, add rug plot showing data distribution
  #' @param title Character, plot title
  
  # Extract smooth
  smooth_data <- gratia::smooth_estimates(model, smooth = term)
  
  # Auto title
  if (is.null(title)) {
    title <- paste0("Smooth Effect: ", gsub("s\\(|\\)", "", term))
  }
  
  # Plot
  p <- ggplot(smooth_data, aes(x = .data[[names(smooth_data)[1]]], y = .estimate)) +
    # Confidence interval
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), 
                alpha = 0.3, fill = "steelblue") +
    
    # Smooth line
    geom_line(color = "steelblue", linewidth = 1.5) +
    
    # Zero line
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    
    # Rug if requested
    {if (rug && !is.null(data)) {
      var_name <- gsub("s\\(|\\)", "", term)
      geom_rug(data = data, aes(x = .data[[var_name]], y = NULL), 
               sides = "b", alpha = 0.1)
    }} +
    
    labs(
      title = title,
      x = gsub("s\\(|\\)", "", term),
      y = "Effect on Response"
    ) +
    theme_flux()
  
  return(p)
}

# ============================================================================
# 4. PARTIAL EFFECTS BY CROP STAGE - Stage-specific responses
# ============================================================================

plot_partial_effects_by_stage <- function(model, data, var, flux_var = "NEE",
                                          stages = NULL) {
  #' Plot partial effects of a variable across different crop stages
  #' 
  #' @param model GAM model object
  #' @param data Data frame (optional, for compatibility)
  #' @param var Character, predictor variable name
  #' @param flux_var Character, response variable name
  #' @param stages Character vector, stages to plot (NULL = all)
  
  # Get model data (actual data used in fitting)
  model_data <- model$model
  
  # Select stages
  if (is.null(stages)) {
    # Get stages that are actually in the model
    available_stages <- unique(model_data$crop_stage_simple)
    stages <- intersect(c("vegetative", "reproductive", "grain_fill", "mature"), 
                        available_stages)
  }
  
  # Create prediction grid
  var_range <- range(model_data[[var]], na.rm = TRUE)
  pred_grid <- expand.grid(
    var_val = seq(var_range[1], var_range[2], length.out = 100),
    crop_stage_simple = stages
  )
  names(pred_grid)[1] <- var
  
  # Get all predictor variable names from the model
  model_vars <- all.vars(formula(model))
  model_vars <- model_vars[model_vars != names(model_data)[1]]  # Remove response
  
  # Add all other predictors at their mean/mode values
  for (col in model_vars) {
    if (col %in% names(pred_grid)) {
      next  # Already added (var or crop_stage_simple)
    }
    
    if (col %in% names(model_data)) {
      if (is.numeric(model_data[[col]])) {
        # Numeric: use mean
        pred_grid[[col]] <- mean(model_data[[col]], na.rm = TRUE)
      } else {
        # Factor/character: use mode (most common value)
        pred_grid[[col]] <- names(sort(table(model_data[[col]]), decreasing = TRUE))[1]
      }
    }
  }
  
  # Get predictions
  pred_grid$predicted <- predict(model, newdata = pred_grid, type = "response")
  pred_grid$se <- predict(model, newdata = pred_grid, type = "response", se.fit = TRUE)$se.fit
  pred_grid$lower <- pred_grid$predicted - 1.96 * pred_grid$se
  pred_grid$upper <- pred_grid$predicted + 1.96 * pred_grid$se
  
  # Plot
  p <- ggplot(pred_grid, aes(x = .data[[var]], y = predicted, 
                             color = crop_stage_simple, fill = crop_stage_simple)) +
    # Confidence ribbons
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    
    # Lines
    geom_line(linewidth = 1.2) +
    
    # Colors
    scale_color_manual(values = stage_colors, name = "Crop Stage") +
    scale_fill_manual(values = stage_colors, name = "Crop Stage") +
    
    # Labels
    labs(
      title = paste0("Partial Effect: ", var, " on ", flux_var),
      subtitle = "Stage-specific responses with 95% CI",
      x = var,
      y = bquote(.(flux_var) ~ "(µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")")
    ) +
    theme_flux() +
    theme(legend.position = "right")
  
  return(p)
}

# ============================================================================
# 5. STAGE COMPARISON - Boxplots with significance
# ============================================================================

plot_flux_by_stage <- function(data, flux_var = "NEE", 
                               add_stats = TRUE, facet_by = NULL) {
  #' Compare flux values across crop stages
  #' 
  #' @param data Data frame with flux data
  #' @param flux_var Character, flux variable to plot
  #' @param add_stats Logical, add sample sizes and medians
  #' @param facet_by Character, variable to facet by (e.g., "management")
  
  # Order stages
  stage_order <- c("fallow", "recovery", "early", "regrowth", "vegetative",
                   "reproductive", "grain_fill", "mature", "dormant")
  data$crop_stage_simple <- factor(data$crop_stage_simple, levels = stage_order)
  
  # Calculate stats if requested
  if (add_stats) {
    stats <- data %>%
      group_by(crop_stage_simple) %>%
      summarise(
        n = n(),
        median_val = median(.data[[flux_var]], na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # Base plot
  p <- ggplot(data, aes(x = crop_stage_simple, y = .data[[flux_var]], 
                        fill = crop_stage_simple)) +
    # Violin + boxplot
    geom_violin(alpha = 0.6, draw_quantiles = 0.5) +
    geom_boxplot(width = 0.2, outlier.alpha = 0.3, alpha = 0.8) +
    
    # Colors
    scale_fill_manual(values = stage_colors, guide = "none") +
    
    # Zero line
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
    
    # Stats labels
    {if (add_stats) {
      geom_text(data = stats, 
                aes(x = crop_stage_simple, label = paste0("n=", n), y = Inf), 
                vjust = 1.5, size = 3, inherit.aes = FALSE)
    }} +
    
    # Facet if requested
    {if (!is.null(facet_by)) facet_wrap(as.formula(paste("~", facet_by)), ncol = 1)} +
    
    # Labels
    labs(
      title = paste0(flux_var, " by Crop Phenology Stage"),
      subtitle = "Distribution comparison across growth stages",
      x = "Crop Stage",
      y = bquote(.(flux_var) ~ "(µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")")
    ) +
    theme_flux() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# ============================================================================
# 6. TIME SERIES WITH STAGES - Show seasonal patterns
# ============================================================================

plot_timeseries_with_stages <- function(data, flux_var = "NEE", 
                                        management = NULL, year = NULL) {
  #' Plot flux time series with crop stage coloring
  #' 
  #' @param data Data frame with date, flux, crop_stage_simple
  #' @param flux_var Character, flux variable name
  #' @param management Character, filter to specific management (optional)
  #' @param year Integer, filter to specific year (optional)
  
  # Filter data
  if (!is.null(management)) {
    data <- data %>% filter(.data$management == !!management)
  }
  if (!is.null(year)) {
    data <- data %>% filter(year(.data$date) == !!year)
  }
  
  # Order stages
  stage_order <- c("fallow", "recovery", "early", "regrowth", "vegetative",
                   "reproductive", "grain_fill", "mature", "dormant")
  data$crop_stage_simple <- factor(data$crop_stage_simple, levels = stage_order)
  
  # Plot
  p <- ggplot(data, aes(x = date, y = .data[[flux_var]])) +
    # Points colored by stage
    geom_point(aes(color = crop_stage_simple), alpha = 0.6, size = 1) +
    
    # Smooth trend
    geom_smooth(color = "black", linewidth = 1, se = TRUE, alpha = 0.2) +
    
    # Zero line
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    
    # Colors
    scale_color_manual(values = stage_colors, name = "Crop Stage") +
    
    # Date scale
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
    
    # Labels
    labs(
      title = paste0(flux_var, " Time Series"),
      subtitle = ifelse(!is.null(management), 
                        paste("Management:", management), 
                        "All fields"),
      x = "Date",
      y = bquote(.(flux_var) ~ "(µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")")
    ) +
    theme_flux() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) +
    guides(color = guide_legend(nrow = 2, override.aes = list(size = 3)))
  
  return(p)
}

# ============================================================================
# 7. MANAGEMENT COMPARISON - Conventional vs Organic
# ============================================================================

plot_management_comparison <- function(data, flux_var = "NEE") {
  #' Compare flux across management systems and crop stages
  #' 
  #' @param data Data frame with management and crop_stage_simple
  #' @param flux_var Character, flux variable
  
  # Order stages
  stage_order <- c("fallow", "recovery", "early", "regrowth", "vegetative",
                   "reproductive", "grain_fill", "mature", "dormant")
  data$crop_stage_simple <- factor(data$crop_stage_simple, levels = stage_order)
  
  # Calculate means and SE
  summary_data <- data %>%
    group_by(management, crop_stage_simple) %>%
    summarise(
      mean_flux = mean(.data[[flux_var]], na.rm = TRUE),
      se_flux = sd(.data[[flux_var]], na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(n >= 5)  # Only show stages with at least 5 observations
  
  # Plot
  p <- ggplot(summary_data, aes(x = crop_stage_simple, y = mean_flux, 
                                fill = management)) +
    # Bars with error bars
    geom_col(position = position_dodge(0.9), alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_flux - se_flux, ymax = mean_flux + se_flux),
                  position = position_dodge(0.9), width = 0.2) +
    geom_jitter() +
    
    # Zero line
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    
    # Colors
    scale_fill_manual(values = management_colors, name = "Management") +
    
    # Labels
    labs(
      title = paste0(flux_var, " by Management System and Crop Stage"),
      subtitle = "Mean ± SE (minimum n=5)",
      x = "Crop Stage",
      y = bquote(.(flux_var) ~ "(µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")")
    ) +
    theme_flux() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
  
  return(p)
}

# ============================================================================
# 8. HEATMAP - Flux patterns across time and stages
# ============================================================================

plot_flux_heatmap <- function(data, flux_var = "NEE", 
                              time_unit = "month", management = NULL) {
  #' Create heatmap of flux values by time and crop stage
  #' 
  #' @param data Data frame
  #' @param flux_var Character, flux variable
  #' @param time_unit Character, "month" or "doy" (day of year)
  #' @param management Character, filter to management (optional)
  
  # Filter
  if (!is.null(management)) {
    data <- data %>% filter(.data$management == !!management)
  }
  
  # Create time variable
  if (time_unit == "month") {
    data$time_var <- month(data$date, label = TRUE)
  } else {
    data$time_var <- yday(data$date)
  }
  
  # Aggregate
  heatmap_data <- data %>%
    group_by(time_var, crop_stage_simple) %>%
    summarise(
      mean_flux = mean(.data[[flux_var]], na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(n >= 3)  # Require at least 3 observations
  
  # Plot
  p <- ggplot(heatmap_data, aes(x = time_var, y = crop_stage_simple, 
                                fill = mean_flux)) +
    geom_tile(color = "white", linewidth = 0.5) +
    
    # Color scale
    scale_fill_gradient2(
      low = "tan3", mid = "white", high = "forestgreen", midpoint = 0,
      name = bquote(.(flux_var) ~ "(µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")")
    ) +
    
    # Labels
    labs(
      title = paste0(flux_var, " Patterns: Time × Crop Stage"),
      subtitle = ifelse(!is.null(management), paste("Management:", management), "All fields"),
      x = ifelse(time_unit == "month", "Month", "Day of Year"),
      y = "Crop Stage"
    ) +
    theme_flux() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  return(p)
}

# ============================================================================
# 9. MULTI-PANEL SUMMARY - All flux variables
# ============================================================================

plot_flux_summary <- function(data, stages = NULL, management = NULL) {
  #' Create multi-panel summary of NEE, GPP, and Reco
  #' 
  #' @param data Data frame with all flux variables
  #' @param stages Character vector, stages to include (NULL = all)
  #' @param management Character, filter to management (optional)
  
  # Filter
  if (!is.null(management)) {
    data <- data %>% filter(.data$management == !!management)
  }
  
  if (!is.null(stages)) {
    data <- data %>% filter(.data$crop_stage_simple %in% stages)
  }
  
  # Create individual plots
  p_nee <- plot_flux_by_stage(data, "NEE", add_stats = FALSE) +
    labs(title = "Net Ecosystem Exchange (NEE)") +
    theme(axis.title.x = element_blank())
  
  p_gpp <- plot_flux_by_stage(data, "GPP_gapfilled", add_stats = FALSE) +
    labs(title = "Gross Primary Productivity (GPP)") +
    theme(axis.title.x = element_blank())
  
  p_reco <- plot_flux_by_stage(data, "Reco_gapfilled", add_stats = FALSE) +
    labs(title = "Ecosystem Respiration (Reco)")
  
  # Combine
  combined <- p_nee / p_gpp / p_reco +
    plot_annotation(
      title = "Carbon Flux Components by Crop Stage",
      subtitle = ifelse(!is.null(management), paste("Management:", management), "All fields"),
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
  
  return(combined)
}

# ============================================================================
# 10. TE() SURFACE - Joint temperature × VPD response by management
# ============================================================================

plot_te_surface <- function(model, data, flux_var = "NEE",
                            temp_var = "air_temperature", vpd_var = "VPD") {
  #' Plot the te(air_temperature, VPD) partial effect surface for each
  #' management type as side-by-side heatmaps. Observed data overlaid as
  #' points to show where the surface is well-supported.
  #'
  #' @param model  GAM model object (must have te(..., by = management) term)
  #' @param data   Data frame used to fit model (for observed-data overlay)
  #' @param flux_var Character, flux variable label for axis/title
  #' @param temp_var Character, temperature column name
  #' @param vpd_var  Character, VPD column name

  # Extract 2D smooth estimates for te() terms
  te_est <- gratia::smooth_estimates(model, unconditional = TRUE) %>%
    gratia::add_confint() %>%
    dplyr::filter(grepl("^te\\(", .smooth)) %>%
    dplyr::mutate(
      Management = stringr::str_to_title(
        as.character(.data[["as.factor(management)"]])
      )
    )

  # Observed data for overlay (density context)
  obs <- data %>%
    dplyr::select(dplyr::all_of(c(temp_var, vpd_var, "management"))) %>%
    dplyr::mutate(Management = stringr::str_to_title(management))

  p <- ggplot(te_est,
              aes(x = .data[[temp_var]], y = .data[[vpd_var]]/1000,
                  fill = .estimate)) +
    geom_tile() +
    # Observed data as translucent points
    geom_point(data = obs,
               aes(x = .data[[temp_var]], y = .data[[vpd_var]]/1000),
               inherit.aes = FALSE,
               colour = "black", alpha = 0.15, size = 0.75) +
    scale_fill_gradient2(
      low = "tan3", mid = "white", high = "forestgreen", midpoint = 0,
      name = bquote(Delta ~ .(flux_var) ~ "\n(g C m"^-2 ~ "d"^-1 ~ ")")
    ) +
    facet_wrap(~ Management, ncol = 2) +
    labs(
      title    = bquote("Combined effect of climate variables on" ~ .(flux_var) ~ "relative to average"),
      subtitle = "dark points = observed data",
      x        = "Air temperature (\u00b0C)",
      y        = "VPD (kPa)"
    ) +
    theme_flux() +
    theme(
      panel.grid  = element_blank(),
      strip.text  = element_text(face = "bold")
    )

  return(p)
}

# ============================================================================
# 11. TILLAGE SMOOTH - s(days_since_tillage) by management
# ============================================================================

plot_tillage_smooth <- function(model, data, flux_var = "NEE") {
  #' Plot the s(days_since_tillage) partial effect smooth for each management
  #' type. Faceted so the flat conventional line and the nonlinear organic
  #' response are directly comparable.
  #'
  #' @param model    GAM model object
  #' @param data     Data frame used to fit model (for rug plot)
  #' @param flux_var Character, flux variable label

  # Extract smooth estimates for the tillage terms
  till_est <- gratia::smooth_estimates(model, unconditional = TRUE) %>%
    gratia::add_confint() %>%
    dplyr::filter(grepl("days_since_tillage", .smooth)) %>%
    dplyr::mutate(
      Management = stringr::str_to_title(
        as.character(.data[["as.factor(management)"]])
      )
    )

  # Rug data
  rug_data <- data %>%
    dplyr::select(days_since_tillage, management) %>%
    dplyr::mutate(Management = stringr::str_to_title(management)) %>%
    dplyr::filter(!is.na(days_since_tillage))

  p <- ggplot(till_est,
              aes(x = days_since_tillage, y = .estimate,
                  colour = Management, fill = Management)) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci),
                alpha = 0.25, colour = NA) +
    geom_line(linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "gray50", linewidth = 0.6) +
    geom_rug(data = rug_data,
             aes(x = days_since_tillage),
             sides = "b", alpha = 0.15,
             inherit.aes = FALSE) +
    scale_colour_manual(values = management_colors, guide = "none") +
    scale_fill_manual(values   = management_colors, guide = "none") +
    facet_wrap(~ Management, ncol = 2, scales = "free_y") +
    labs(
      title    = paste0("Effect of days since tillage on ", flux_var),
      subtitle = "Conventional: non-significant (p = 0.962) | Organic: nonlinear, p = 0.040 (EDF = 3.03)",
      x        = "Days since tillage",
      # y        = bquote("Partial effect on" ~ .(flux_var) ~ "(g C m"^-2 ~ "d"^-1 ~ ")")
      y        = bquote("Partial effect on" ~ .(flux_var) ~ "(µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")")
    ) +
    theme_flux() +
    theme(strip.text = element_text(face = "bold"))

  return(p)
}

# ============================================================================
# 12. COEFFICIENT PLOT - Parametric interaction terms (main finding figure)
# ============================================================================

plot_coef_plot <- function(model, flux_var = "NEE",
                           exclude_intercept = TRUE,
                           ci_level = 0.95) {
  #' Dot-and-whisker (forest) plot of parametric fixed-effect coefficients.
  #'
  #' The main story: crop stage MAIN EFFECTS (red) show how conventional NEE
  #' shifts relative to conventional-at-mature. Organic x stage INTERACTIONS
  #' (green) show how organic departs from that same conventional pattern at
  #' each stage. Faded points = p >= 0.05.
  #'
  #' @param model           GAM/bam model object
  #' @param flux_var        Character, flux variable label for axis/title
  #' @param exclude_intercept Logical, omit the intercept row (default TRUE)
  #' @param ci_level        Numeric, confidence level for intervals (default 0.95)

  z <- qnorm(1 - (1 - ci_level) / 2)

  coef_df <- broom::tidy(model, parametric = TRUE) %>%
    dplyr::mutate(
      label = dplyr::case_when(
        term == "(Intercept)"                                      ~ "Intercept (conventional, mature)",
        term == "managementorganic"                                ~ "Organic vs. conventional",
        term == "crop_stage_simpledormant"                         ~ "Dormant",
        term == "crop_stage_simpleearly"                           ~ "Early growth",
        term == "crop_stage_simplefallow"                          ~ "Fallow",
        term == "crop_stage_simplegrain_fill"                      ~ "Grain fill",
        term == "crop_stage_simplereproductive"                    ~ "Reproductive",
        term == "crop_stage_simplevegetative"                      ~ "Vegetative",
        term == "managementorganic:crop_stage_simpledormant"       ~ "Organic \u00d7 Dormant",
        term == "managementorganic:crop_stage_simpleearly"         ~ "Organic \u00d7 Early growth",
        term == "managementorganic:crop_stage_simplefallow"        ~ "Organic \u00d7 Fallow",
        term == "managementorganic:crop_stage_simplegrain_fill"    ~ "Organic \u00d7 Grain fill",
        term == "managementorganic:crop_stage_simplereproductive"  ~ "Organic \u00d7 Reproductive",
        term == "managementorganic:crop_stage_simplevegetative"    ~ "Organic \u00d7 Vegetative",
        TRUE ~ term
      ),
      term_type = dplyr::case_when(
        grepl("Intercept|^Organic vs\\.", label) ~ "Reference / management main effect",
        grepl("Organic \u00d7", label)           ~ "Organic \u00d7 stage interaction",
        TRUE                                     ~ "Conventional crop stage effect"
      ),
      lower = estimate - z * std.error,
      upper = estimate + z * std.error,
      sig   = p.value < 0.05
    )

  if (exclude_intercept) {
    coef_df <- coef_df %>% dplyr::filter(term != "(Intercept)")
  }

  # Fix display order: management main effect, then stage effects, then interactions
  label_order <- c(
    "Organic vs. conventional",
    "Dormant", "Early growth", "Fallow", "Vegetative", "Reproductive", "Grain fill",
    "Organic \u00d7 Dormant", "Organic \u00d7 Early growth", "Organic \u00d7 Fallow",
    "Organic \u00d7 Vegetative", "Organic \u00d7 Reproductive", "Organic \u00d7 Grain fill"
  )
  # Keep only levels actually present
  label_order  <- label_order[label_order %in% coef_df$label]
  coef_df$label <- factor(coef_df$label, levels = rev(label_order))

  type_colors <- c(
    "Reference / management main effect" = "#555555",
    "Conventional crop stage effect"     = "#E74C3C",
    "Organic \u00d7 stage interaction"   = "#27AE60"
  )

  # Separator line between main effects and interactions
  n_main <- sum(coef_df$term_type != "Organic \u00d7 stage interaction")

  p <- ggplot(coef_df,
              aes(x = estimate, y = label,
                  color = term_type, alpha = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "gray50", linewidth = 0.7) +
    geom_pointrange(aes(xmin = lower, xmax = upper),
                    linewidth = 0.8, size = 0.55,
                    position = position_dodge(0)) +
    scale_color_manual(values = type_colors, name = NULL) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35), guide = "none") +
    labs(
      title    = paste0("Parametric fixed effects: ", flux_var, " GAM"),
      subtitle = paste0(
        round(ci_level * 100), "% CI shown  |  faded = p \u2265 0.05  |  ",
        "Reference: conventional \u00d7 mature"
      ),
      # x = bquote("Coefficient estimate (g C m"^-2 ~ "d"^-1 ~ ")"),
      x = bquote("Coefficient estimate (µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")"),
      y = NULL
    ) +
    theme_flux() +
    theme(
      legend.position    = "bottom",
      legend.direction   = "horizontal",
      panel.grid.major.y = element_line(color = "gray92", linewidth = 0.3),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3)
    ) +
    guides(color = guide_legend(override.aes = list(size = 0.8, linewidth = 1)))

  return(p)
}


# ============================================================================
# 13. PREDICTED NEE BY MANAGEMENT x CROP STAGE
# ============================================================================

plot_predicted_by_stage <- function(model, data, flux_var = "NEE",
                                    ci_level   = 0.95,
                                    exclude_re = TRUE) {
  z <- qnorm(1 - (1 - ci_level) / 2)
  
  stage_levels <- levels(droplevels(model$model$crop_stage_simple))
  mgmt_levels  <- levels(droplevels(model$model$management))
  
  # Check the MODEL'S variables, not data, to avoid ambiguity when both exist
  model_vars <- names(model$model)
  
  if (all(c("air_temperature_7", "VPD_7") %in% model_vars)) {
    temp_mean     <- mean(data$air_temperature_7, na.rm = TRUE)
    vpd_mean      <- mean(data$VPD_7, na.rm = TRUE)
    temp_var_name <- "air_temperature_7"
    vpd_var_name  <- "VPD_7"
  } else if (all(c("air_temperature", "VPD") %in% model_vars)) {
    temp_mean     <- mean(data$air_temperature, na.rm = TRUE)
    vpd_mean      <- mean(data$VPD, na.rm = TRUE)
    temp_var_name <- "air_temperature"
    vpd_var_name  <- "VPD"
  } else {
    stop("Neither air_temperature/VPD nor air_temperature_7/VPD_7 found in model.")
  }
  
  newdata <- expand.grid(
    management        = mgmt_levels,
    crop_stage_simple = stage_levels,
    stringsAsFactors  = FALSE
  ) %>%
    dplyr::mutate(
      management        = factor(management, levels = mgmt_levels),
      crop_stage_simple = factor(crop_stage_simple, levels = stage_levels),
      doy               = as.integer(median(data$doy, na.rm = TRUE)),
      days_since_tillage = mean(data$days_since_tillage, na.rm = TRUE),
      year_f            = factor(levels(droplevels(model$model$year_f))[2],
                                 levels = levels(droplevels(model$model$year_f))),
      !!sym(temp_var_name) := temp_mean,
      !!sym(vpd_var_name)  := vpd_mean
    )
  
  re_terms <- if (exclude_re) {
    re_labels <- sapply(model$smooth, function(s) s$label)
    re_labels[grepl("year_f", re_labels)]
  } else {
    character(0)
  }
  
  preds <- predict(model, newdata = newdata, type = "response",
                   se.fit = TRUE, exclude = re_terms)
  
  newdata$fit   <- preds$fit
  newdata$lower <- preds$fit - z * preds$se.fit
  newdata$upper <- preds$fit + z * preds$se.fit
  
  xorder <- c("dormant", "early", "vegetative", "reproductive",
              "grain_fill", "mature", "fallow")
  xorder <- xorder[xorder %in% stage_levels]
  newdata$crop_stage_simple <- factor(newdata$crop_stage_simple, levels = xorder)
  
  stage_labels <- c(
    dormant      = "Dormant",
    early        = "Early\ngrowth",
    vegetative   = "Vegetative",
    reproductive = "Reproductive",
    grain_fill   = "Grain\nfill",
    mature       = "Mature",
    fallow       = "Fallow"
  )
  
  p <- ggplot(newdata,
              aes(x = crop_stage_simple, y = fit,
                  color = management, group = management)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "gray50", linewidth = 0.6) +
    geom_line(position = position_dodge(0.45),
              linewidth = 0.8, alpha = 0.55) +
    geom_pointrange(aes(ymin = lower, ymax = upper),
                    position = position_dodge(0.45),
                    linewidth = 0.9, size = 0.65) +
    scale_color_manual(
      values = management_colors,
      name   = "Management",
      labels = c("conventional" = "Conventional", "organic" = "Organic")
    ) +
    scale_x_discrete(labels = stage_labels) +
    labs(
      title    = paste0("Predicted ", flux_var, " by management \u00d7 crop stage"),
      subtitle = paste0(
        round(ci_level * 100),
        "% CI  |  Climate predictors held at means  |  Year random effects excluded"
      ),
      x = "Crop stage",
      
      y = bquote("Predicted" ~ .(flux_var) ~ "(µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")")
      # y = bquote("Predicted" ~ .(flux_var) ~ "(g C m"^-2 ~ "d"^-1 ~ ")")
    ) +
    theme_flux() +
    theme(
      legend.position = "top",
      axis.text.x     = element_text(size = 10)
    )
  
  return(p)
}


# ============================================================================
# 14. EXAMPLE USAGE
# ============================================================================

# Example function to generate all plots
generate_all_plots <- function(data, model, output_dir = "./plots") {
  #' Generate all standard plots and save to directory
  #' 
  #' @param data Data frame with flux and predictor data
  #' @param model Fitted GAM model
  #' @param output_dir Character, directory to save plots
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. Model performance
  p1 <- plot_gam_observed_vs_predicted(model, data, "NEE")
  ggsave(file.path(output_dir, "01_observed_vs_predicted.png"), 
         p1, width = 8, height = 6, dpi = 300)
  
  # 2. Residuals
  p2 <- plot_gam_residuals_by_stage(model, data, "NEE")
  ggsave(file.path(output_dir, "02_residuals_by_stage.png"), 
         p2, width = 10, height = 6, dpi = 300)
  
  # 3. Flux by stage
  p3 <- plot_flux_by_stage(data, "NEE", add_stats = TRUE)
  ggsave(file.path(output_dir, "03_nee_by_stage.png"), 
         p3, width = 10, height = 6, dpi = 300)
  
  # 4. Management comparison
  p4 <- plot_management_comparison(data, "NEE")
  ggsave(file.path(output_dir, "04_management_comparison.png"), 
         p4, width = 12, height = 6, dpi = 300)
  
  # 5. Time series (2019 as representative mid-study year)
  p5 <- plot_timeseries_with_stages(data, "NEE", year = 2019)
  ggsave(file.path(output_dir, "05_timeseries_2019.png"),
         p5, width = 12, height = 6, dpi = 300)

  # 6. Heatmap
  p6 <- plot_flux_heatmap(data, "NEE", time_unit = "month")
  ggsave(file.path(output_dir, "06_heatmap.png"),
         p6, width = 10, height = 6, dpi = 300)

  # 7. Multi-panel summary
  p7 <- plot_flux_summary(data)
  ggsave(file.path(output_dir, "07_flux_summary.png"),
         p7, width = 10, height = 12, dpi = 300)

  # 8. te() response surface: joint temp x VPD effect by management
  p8 <- plot_te_surface(model, data, "NEE")
  ggsave(file.path(output_dir, "08_te_surface_NEE.png"),
         p8, width = 12, height = 5, dpi = 300)

  # 9. Days since tillage smooth by management
  p9 <- plot_tillage_smooth(model, data, "NEE")
  ggsave(file.path(output_dir, "09_tillage_smooth_NEE.png"),
         p9, width = 10, height = 5, dpi = 300)

  # 10. Coefficient plot: parametric fixed effects (main finding figure)
  p10 <- plot_coef_plot(model, flux_var = "NEE")
  ggsave(file.path(output_dir, "10_coef_plot_NEE.png"),
         p10, width = 9, height = 7, dpi = 300)

  # 11. Predicted NEE by management x crop stage (population-level marginal means)
  p11 <- plot_predicted_by_stage(model, data, flux_var = "NEE")
  ggsave(file.path(output_dir, "11_predicted_by_stage_NEE.png"),
         p11, width = 9, height = 6, dpi = 300)

  cat("All plots saved to:", output_dir, "\n")
}

# ============================================================================
# Example GAM model specification for reference
# ============================================================================

# Example model (adjust based on your actual predictors):
# library(mgcv)
# 
# model_nee <- gam(
#   NEE ~ 
#     crop_stage_simple +
#     s(air_temperature, by = crop_stage_simple, k = 5) +
#     s(VPD, k = 5) +
#     s(PAR, k = 5) +
#     management,
#   data = data,
#   method = "REML"
# )

cat("GAM visualization scripts loaded successfully!\n")
cat("Available functions:\n")
cat("  - plot_gam_observed_vs_predicted()  # observed vs predicted, R\u00b2 + RMSE\n")
cat("  - plot_gam_residuals_by_stage()     # residual violin/box by crop stage\n")
cat("  - plot_gam_smooth_effects()         # 1-D smooth effects (doy, etc.)\n")
cat("  - plot_partial_effects_by_stage()   # partial effects across stages\n")
cat("  - plot_flux_by_stage()              # observed flux distributions by stage\n")
cat("  - plot_timeseries_with_stages()     # time series coloured by crop stage\n")
cat("  - plot_management_comparison()      # mean \u00b1 SE by management x stage\n")
cat("  - plot_flux_heatmap()               # time x stage heatmap\n")
cat("  - plot_flux_summary()               # multi-panel NEE/GPP/Reco\n")
cat("  - plot_te_surface()                 # te(Temp,VPD) 2-D surface by management\n")
cat("  - plot_tillage_smooth()             # s(days_since_tillage) by management\n")
cat("  - plot_coef_plot()                  # dot-and-whisker: parametric fixed effects\n")
cat("  - plot_predicted_by_stage()         # model-estimated marginal means by mgmt x stage\n")
cat("  - generate_all_plots()              # save all standard plots to directory\n")
