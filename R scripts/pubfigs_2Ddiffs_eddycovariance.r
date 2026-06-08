# function for 2d VPD/air temperature plots on each flux:

plot_te_difference <- function(model, data, flux_var = "NEE",
                               temp_var = "air_temperature", vpd_var = "VPD") {
  #' Plot the difference (Organic - Conventional) in the te(temp, VPD) partial
  #' effect surface. Isolates the management contrast in joint climate response,
  #' independent of parametric management offsets or other smooths.
  
  # Extract te() estimates for each management level
  te_est <- gratia::smooth_estimates(model, unconditional = TRUE) %>%
    gratia::add_confint() %>%
    dplyr::filter(grepl("^te\\(", .smooth)) %>%
    dplyr::mutate(
      Management = stringr::str_to_title(
        as.character(.data[["as.factor(management)"]])
      )
    )
  
  # Pivot to get Organic and Conventional in side-by-side columns
  te_diff <- te_est %>%
    dplyr::select(dplyr::all_of(c(temp_var, vpd_var)),
                  Management, .estimate, .se) %>%
    tidyr::pivot_wider(names_from = Management,
                       values_from = c(.estimate, .se)) %>%
    dplyr::mutate(
      diff    = .estimate_Organic - .estimate_Conventional,
      diff_se = sqrt(.se_Organic^2 + .se_Conventional^2)
    )
  
  # Mask grid cells far from observed data (too.far-style)
  obs <- data %>%
    dplyr::select(dplyr::all_of(c(temp_var, vpd_var))) %>%
    na.omit()
  nn <- fields::rdist(te_diff[, c(temp_var, vpd_var)], as.matrix(obs))
  te_diff$min_dist <- apply(nn, 1, min)
  te_diff$diff_masked <- ifelse(
    te_diff$min_dist < quantile(te_diff$min_dist, 0.8, na.rm = TRUE),
    te_diff$diff, NA_real_
  )
  
  # Auto-place corner labels based on data ranges
  xr <- range(te_diff[[temp_var]], na.rm = TRUE)
  yr <- range(te_diff[[vpd_var]], na.rm = TRUE)
  x_lo <- xr[1] + 0.10 * diff(xr); x_hi <- xr[2] - 0.10 * diff(xr)
  y_lo <- yr[1] + 0.07 * diff(yr); y_hi <- yr[2] - 0.07 * diff(yr)
  
  ggplot(te_diff, aes(x = .data[[temp_var]], y = .data[[vpd_var]],
                      fill = diff_masked)) +
    geom_tile() +
    geom_point(data = data,
               aes(x = .data[[temp_var]], y = .data[[vpd_var]]),
               inherit.aes = FALSE,
               colour = "black", alpha = 0.15, size = 0.7) +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#d73027", midpoint = 0,
      name = bquote("Organic - Conventional\n" * Delta ~ .(flux_var) ~
                      "(g C m"^-2 ~ "d"^-1 * ")"),
      na.value = "grey92"
    ) +
    annotate("label", x = x_lo, y = y_lo, label = "cool & humid\n(low growth)",
             color = "grey20", size = 2.6, label.size = 0,
             fill = alpha("white", 0.6)) +
    annotate("label", x = x_hi, y = y_lo, label = "warm & humid\n(productive)",
             color = "grey20", size = 2.6, label.size = 0,
             fill = alpha("white", 0.6)) +
    annotate("label", x = x_lo, y = y_hi, label = "cool & dry air\n(uncommon)",
             color = "grey20", size = 2.6, label.size = 0,
             fill = alpha("white", 0.6)) +
    annotate("label", x = x_hi, y = y_hi,
             label = "hot & dry air\n(drought stress)",
             color = "grey20", size = 2.6, label.size = 0,
             fill = alpha("white", 0.6)) +
    # annotate("segment",
    #          x = x_lo, y = y_lo + 0.12 * diff(yr),
    #          xend = x_hi, yend = y_hi - 0.12 * diff(yr),
    #          arrow = arrow(type = "closed", length = unit(3, "mm")),
    #          color = "grey20", linewidth = 0.6) +
    labs(
      title    = bquote("Management contrast in" ~ .(flux_var) ~
                          "across the temp x VPD envelope"),
      subtitle = "Blue = organic takes up more CO2; gray = unobserved",
      x        = "Air temperature (\u00b0C)",
      y        = "VPD (Pa)"
    ) +
    theme_flux() +
    theme(panel.grid = element_blank())
}

plot_te_difference(bam_NEE, data1, flux_var = "NEE", temp_var = "air_temperature", vpd_var = "VPD")

plot_te_difference_2 <- function(model, data, flux_var = "NEE",
                               temp_var = "air_temperature", vpd_var = "VPD") {
  #' Plot the difference (Organic - Conventional) in the te(temp, VPD) partial
  #' effect surface. Isolates the management contrast in joint climate response,
  #' independent of parametric management offsets or other smooths.
  
  # Extract te() estimates for each management level
  te_est <- gratia::smooth_estimates(model, unconditional = TRUE) %>%
    gratia::add_confint() %>%
    dplyr::filter(grepl("^te\\(", .smooth)) %>%
    dplyr::mutate(
      Management = stringr::str_to_title(
        as.character(.data[["as.factor(management)"]])
      )
    )
  
  # Pivot to get Organic and Conventional in side-by-side columns
  te_diff <- te_est %>%
    dplyr::select(dplyr::all_of(c(temp_var, vpd_var)),
                  Management, .estimate, .se) %>%
    tidyr::pivot_wider(names_from = Management,
                       values_from = c(.estimate, .se)) %>%
    dplyr::mutate(
      diff    = .estimate_Organic - .estimate_Conventional,
      diff_se = sqrt(.se_Organic^2 + .se_Conventional^2)
    )
  
  # Mask grid cells far from observed data (too.far-style)
  obs <- data %>%
    dplyr::select(dplyr::all_of(c(temp_var, vpd_var))) %>%
    na.omit()
  nn <- fields::rdist(te_diff[, c(temp_var, vpd_var)], as.matrix(obs))
  te_diff$min_dist <- apply(nn, 1, min)
  te_diff$diff_masked <- ifelse(
    te_diff$min_dist < quantile(te_diff$min_dist, 0.8, na.rm = TRUE),
    te_diff$diff, NA_real_
  )
  
  # Auto-place corner labels based on data ranges
  xr <- range(te_diff[[temp_var]], na.rm = TRUE)
  yr <- range(te_diff[[vpd_var]], na.rm = TRUE)
  x_lo <- xr[1] + 0.10 * diff(xr); x_hi <- xr[2] - 0.10 * diff(xr)
  y_lo <- yr[1] + 0.07 * diff(yr); y_hi <- yr[2] - 0.07 * diff(yr)
  
  ggplot(te_diff, aes(x = .data[[temp_var]], y = .data[[vpd_var]],
                      fill = diff_masked)) +
    geom_tile() +
    geom_point(data = data,
               aes(x = .data[[temp_var]], y = .data[[vpd_var]]),
               inherit.aes = FALSE,
               colour = "black", alpha = 0.15, size = 0.7) +
    # switch the color scale around to align with GPP sign convention:
    scale_fill_gradient2(
      high = "#2166ac", mid = "white", low = "#d73027", midpoint = 0,
      name = bquote("Organic - Conventional\n" * Delta ~ .(flux_var) ~
                      "(g C m"^-2 ~ "d"^-1 * ")"),
      na.value = "grey92"
    ) +
    annotate("label", x = x_lo, y = y_lo, label = "cool & humid\n(low growth)",
             color = "grey20", size = 2.6, label.size = 0,
             fill = alpha("white", 0.6)) +
    annotate("label", x = x_hi, y = y_lo, label = "warm & humid\n(productive)",
             color = "grey20", size = 2.6, label.size = 0,
             fill = alpha("white", 0.6)) +
    annotate("label", x = x_lo, y = y_hi, label = "cool & dry air\n(uncommon)",
             color = "grey20", size = 2.6, label.size = 0,
             fill = alpha("white", 0.6)) +
    annotate("label", x = x_hi, y = y_hi,
             label = "hot & dry air\n(drought stress)",
             color = "grey20", size = 2.6, label.size = 0,
             fill = alpha("white", 0.6)) +
    # annotate("segment",
    #          x = x_lo, y = y_lo + 0.12 * diff(yr),
    #          xend = x_hi, yend = y_hi - 0.12 * diff(yr),
    #          arrow = arrow(type = "closed", length = unit(3, "mm")),
    #          color = "grey20", linewidth = 0.6) +
    labs(
      title    = bquote("Management contrast in" ~ .(flux_var) ~
                          "across the temp x VPD envelope"),
      subtitle = "Blue = organic takes up more CO2; gray = unobserved",
      x        = "Air temperature (\u00b0C)",
      y        = "VPD (Pa)"
    ) +
    theme_flux() +
    theme(panel.grid = element_blank())
}

plot_te_difference_2(bam_GPP, data1, flux_var = "GPP", temp_var = "air_temperature", vpd_var = "VPD")

plot_te_difference_3 <- function(model, data, flux_var = "NEE",
                               temp_var = "air_temperature", vpd_var = "VPD",
                               mask_quantile = 0.85) {
  
  # Accept either a gam/bam object or a gamm result; pull out the gam piece
  if (inherits(model, "gamm") || (is.list(model) && !is.null(model$gam))) {
    model <- model$gam
  }
  
  # Get smooth estimates WITHOUT unconditional (which gamm doesn't support)
  te_est <- gratia::smooth_estimates(model) %>%
    gratia::add_confint() %>%
    dplyr::filter(grepl("^te\\(", .smooth))
  
  # Be tolerant of either as.factor(management) or management as the column name
  mgmt_col <- intersect(c("as.factor(management)", "management"), names(te_est))
  if (length(mgmt_col) == 0) {
    stop("Couldn't find management column in smooth_estimates output. Available: ",
         paste(names(te_est), collapse = ", "))
  }
  mgmt_col <- mgmt_col[1]
  
  te_est <- te_est %>%
    dplyr::mutate(Management = stringr::str_to_title(as.character(.data[[mgmt_col]])))
  
  # Pivot to compute organic - conventional
  te_diff <- te_est %>%
    dplyr::select(dplyr::all_of(c(temp_var, vpd_var)),
                  Management, .estimate, .se) %>%
    tidyr::pivot_wider(names_from = Management,
                       values_from = c(.estimate, .se)) %>%
    dplyr::mutate(
      diff    = .estimate_Organic - .estimate_Conventional,
      diff_se = sqrt(.se_Organic^2 + .se_Conventional^2)
    )
  
  # Mask far-from-data
  obs <- data %>%
    dplyr::select(dplyr::all_of(c(temp_var, vpd_var))) %>%
    na.omit()
  nn <- fields::rdist(te_diff[, c(temp_var, vpd_var)], as.matrix(obs))
  te_diff$min_dist <- apply(nn, 1, min)
  te_diff$diff_masked <- ifelse(
    te_diff$min_dist < quantile(te_diff$min_dist, mask_quantile, na.rm = TRUE),
    te_diff$diff, NA_real_
  )
  
  # Auto-place corner labels
  xr <- range(te_diff[[temp_var]], na.rm = TRUE)
  yr <- range(te_diff[[vpd_var]], na.rm = TRUE)
  x_lo <- xr[1] + 0.10 * diff(xr); x_hi <- xr[2] - 0.10 * diff(xr)
  y_lo <- yr[1] + 0.07 * diff(yr); y_hi <- yr[2] - 0.07 * diff(yr)
  
  ggplot(te_diff, aes(x = .data[[temp_var]], y = .data[[vpd_var]],
                      fill = diff_masked)) +
    geom_tile() +
    geom_point(data = data,
               aes(x = .data[[temp_var]], y = .data[[vpd_var]]),
               inherit.aes = FALSE,
               colour = "black", alpha = 0.15, size = 0.7) +
    scale_fill_gradient2(
      high = "#2166ac", mid = "white", low = "#d73027", midpoint = 0,
      name = bquote("Organic - Conventional\n" * Delta ~ .(flux_var) ~
                      "(g C m"^-2 ~ "d"^-1 * ")"),
      na.value = "grey92"
    ) +
    annotate("label", x = x_lo, y = y_lo, label = "cool & humid\n(low metabolic forcing)",
             color = "grey20", size = 2.6, label.size = 0,
             fill = alpha("white", 0.8)) +
    annotate("label", x = x_hi, y = y_lo, label = "warm & humid\n(active soils)",
             color = "grey20", size = 2.6, label.size = 0,
             fill = alpha("white", 0.8)) +
    annotate("label", x = x_lo, y = y_hi, label = "cool & dry air\n(uncommon)",
             color = "grey20", size = 2.6, label.size = 0,
             fill = alpha("white", 0.8)) +
    annotate("label", x = x_hi, y = y_hi,
             label = "hot & dry air\n(atmospheric drought)",
             color = "grey20", size = 2.6, label.size = 0,
             fill = alpha("white", 0.8)) +
    # annotate("segment",
    #          x = x_lo, y = y_lo + 0.12 * diff(yr),
    #          xend = x_hi, yend = y_hi - 0.12 * diff(yr),
    #          arrow = arrow(type = "closed", length = unit(3, "mm")),
    #          color = "grey20", linewidth = 0.6) +
    labs(
      title    = bquote("Management contrast in" ~ .(flux_var) ~
                          "across the temp x VPD envelope"),
      subtitle = "Color = organic minus conventional; gray = unobserved",
      x        = "Air temperature (\u00b0C)",
      y        = "VPD (Pa)"
    ) +
    theme_flux() +
    theme(panel.grid = element_blank())
}
plot_te_difference_3(gamm_Reco_smoothed$gam, data1, flux_var = "Reco", temp_var = "air_temperature_7", vpd_var = "VPD_7")

# FOR RECO: 
"The temperature × VPD response of ecosystem respiration differed in sensitivity 
# rather than direction between management systems, with organic R_eco showing 
# a steeper increase along the stress gradient than conventional R_eco. This 
# greater respiratory responsiveness coincided with substantially greater 
# photosynthetic uptake under the same conditions, yielding a net management 
# advantage for organic systems despite their elevated respiratory cost (Fig. X)."

# FOR NEE:
# In your color convention (red = positive ΔNEE = organic releases more / takes up less; 
# blue = negative ΔNEE = organic takes up more), the resilience signal lives in the 
# top-right "atmospheric drought" corner. If organic is more resilient to high VPD 
# and high temperature, that corner should be blue — organic maintains more uptake 
# than conventional when the air is hot and dry. The deeper the blue, the stronger 
# the resilience signal at peak atmospheric stress.

# The bottom-left "cool & humid" corner should be near zero (white) — under benign 
# atmospheric conditions there's nothing to differentiate the two systems, so the 
# management effect should fade. If you see strong color there, it's likely either 
# the parametric management offset showing through (if the systems differ on 
# average across all conditions) or extrapolation from sparser regions, both worth investigating.

# The bottom-right "warm & humid" corner is the productive optimum and may also 
# approach zero if both systems perform similarly under good conditions. A persistent 
# blue here would suggest organic generally outperforms even in benign conditions 
# — a "more productive on average" claim rather than a resilience claim. 
# Worth being able to distinguish.

# The top-left "cool & dry air" corner is going to be sparse and the masking 
# should mostly gray it out. Good — that corner often has so little data that 
# any pattern there is unreliable, and showing it in gray rather than as a colored 
# prediction signals to the reader "we don't have enough data to say."

# "Across the joint envelope of air temperature and atmospheric water demand,
# the organic system maintained substantially greater net CO2 uptake than the 
# conventional system specifically under conditions of co-occurring high 
# temperature and high vapor pressure deficit, with the management contrast 
# intensifying along the gradient from benign to atmospherically stressed conditions (Fig. X)."

