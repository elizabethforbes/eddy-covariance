# =============================================================================
# plot_te_difference.R

# written in collab. with Claude, May 11th 2026
#
# Single function for producing (Contrast - Reference) tensor-product difference
# plots from mgcv GAM / BAM / GAMM models with by-management te() smooths.
#
# Designed to give a CONSISTENT color reading across NEE, GPP, and R_eco:
#
#   * NEE  -> better_direction = "negative"
#            (more negative diff = contrast level takes up more CO2 = favorable)
#   * GPP  -> better_direction = "positive"
#            (more positive diff = contrast level produces more = favorable)
#   * Reco -> better_direction = "neutral"
#            (more respiration is not unambiguously good or bad; no flip)
#
# In all three cases, the function ensures BLUE = contrast-favorable (or
# simply "contrast minus reference" for the neutral case) and RED = the
# opposite direction. The legend subtitle adapts so readers don't have to
# decode the color logic separately for each flux.
#
# Supports both temp x VPD layouts (EC tower data) and temp x soil moisture
# layouts (chamber data), with appropriate corner labels and stress-gradient
# arrow direction for each.
#
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(gratia)
library(stringr)
library(fields)


#' Plot tensor-product (contrast - reference) climate response surface
#'
#' Extracts te() partial effects from a fitted mgcv model, computes the
#' difference between two management levels at matched (temp, climate) grid
#' points, masks regions far from observed data, and renders a heat-map style
#' difference plot with corner labels, an optional stress-gradient arrow, and
#' a consistent blue=favorable color logic.
#'
#' @param model           A fitted mgcv gam/bam object, or a gamm result list.
#'                        The model must contain a te(x, y, by = factor) term
#'                        with `by` indicating management.
#' @param data            The data frame used to fit the model (for too-far
#'                        masking and optional point overlay).
#' @param flux_var        Character. Flux name for axis/title (e.g. "NEE",
#'                        "GPP", "R[eco]"). Parsed as a plotmath expression so
#'                        you can use subscripts.
#' @param temp_var        Column name for the temperature predictor.
#' @param climate_var     Column name for the second climate predictor
#'                        (VPD for EC layout, soil moisture for chamber).
#' @param better_direction One of "negative", "positive", or "neutral".
#'                        Controls the color flip and legend wording.
#' @param reference_level Title-case label of the reference management level
#'                        (subtracted FROM). Default "Conventional".
#' @param contrast_level  Title-case label of the contrast management level
#'                        (subtracted INTO). Default "Organic".
#' @param mgmt_col_candidates Candidate column names where gratia stores the
#'                        management factor. Typically includes both with and
#'                        without the as.factor() wrapper.
#' @param layout          "temp_vpd" or "temp_sm". Drives corner labels and
#'                        stress-arrow direction.
#' @param climate_units   Character used in the y-axis label. Default depends
#'                        on layout.
#' @param mask_quantile   Quantile of nearest-neighbor distance above which
#'                        grid cells are masked as unobserved (default 0.85).
#' @param show_points     Logical. Overlay observed (x, y) points?
#' @param show_arrow      Logical. Draw the diagonal stress-gradient arrow?
#' @param show_corner_labels Logical. Draw the four corner labels?
#' @param corner_labels   Optional named list with elements 'bl', 'br', 'tl',
#'                        'tr' (bottom-left, bottom-right, top-left, top-right)
#'                        to override defaults. Each value is a string; use
#'                        "\n" for line breaks.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # NEE: organic-favorable = more negative diff
#' plot_te_difference(bam_NEE, data_NEE,
#'                    flux_var = "NEE",
#'                    temp_var = "air_temperature", climate_var = "VPD",
#'                    better_direction = "negative", layout = "temp_vpd")
#'
#' # GPP: organic-favorable = more positive diff
#' plot_te_difference(bam_GPP, data_GPP,
#'                    flux_var = "GPP",
#'                    temp_var = "air_temperature", climate_var = "VPD",
#'                    better_direction = "positive", layout = "temp_vpd")
#'
#' # R_eco: no clear "favorable" direction; use neutral
#' plot_te_difference(gamm_Reco_smoothed, data1,
#'                    flux_var = "R[eco]",
#'                    temp_var = "air_temperature_7", climate_var = "VPD_7",
#'                    better_direction = "neutral", layout = "temp_vpd")
#' }
plot_te_difference <- function(
    model,
    data,
    flux_var          = "NEE",
    temp_var          = "air_temperature",
    climate_var       = "VPD",
    better_direction  = c("negative", "positive", "neutral"),
    reference_level   = "Conventional",
    contrast_level    = "Organic",
    mgmt_col_candidates = c("as.factor(management)", "management"),
    layout            = c("temp_vpd", "temp_sm"),
    climate_units     = NULL,
    mask_quantile     = 0.85,
    show_points       = TRUE,
    show_arrow        = TRUE,
    show_corner_labels = TRUE,
    corner_labels     = NULL
) {

  better_direction <- match.arg(better_direction)
  layout           <- match.arg(layout)

  # ---- 1. Unwrap gamm objects --------------------------------------------
  # IMPORTANT: use [["gam"]] (not $gam) to avoid R's partial-matching, which
  # would otherwise grab model$gamma on bam/gam objects and silently break
  # everything downstream. Only unwrap if the object is NOT already a gam/bam
  # and has a 'gam' component (which is the gamm() return-list pattern).
  if (!inherits(model, "gam") && is.list(model) && "gam" %in% names(model)) {
    model <- model[["gam"]]
  }

  # ---- 2. Extract te() partial effects -----------------------------------
  # NOTE: unconditional = TRUE doesn't work for gamm fits (no Vc matrix),
  # so we drop it for portability.
  te_est <- gratia::smooth_estimates(model) |>
    gratia::add_confint() |>
    dplyr::filter(grepl("^te\\(", .smooth))

  if (nrow(te_est) == 0) {
    stop("No te() smooth terms found in the model.")
  }

  # ---- 3. Locate management column ---------------------------------------
  mgmt_col <- intersect(mgmt_col_candidates, names(te_est))
  if (length(mgmt_col) == 0) {
    stop("Management column not found in smooth_estimates output.\n",
         "  Available columns: ", paste(names(te_est), collapse = ", "), "\n",
         "  Pass mgmt_col_candidates explicitly if your column has a different name.")
  }
  mgmt_col <- mgmt_col[1]

  te_est <- te_est |>
    dplyr::mutate(Management = stringr::str_to_title(as.character(.data[[mgmt_col]])))

  # ---- 4. Pivot wide and compute difference ------------------------------
  te_diff <- te_est |>
    dplyr::select(dplyr::all_of(c(temp_var, climate_var)),
                  Management, .estimate, .se) |>
    tidyr::pivot_wider(names_from  = Management,
                       values_from = c(.estimate, .se))

  est_c <- paste0(".estimate_", contrast_level)
  est_r <- paste0(".estimate_", reference_level)
  se_c  <- paste0(".se_",       contrast_level)
  se_r  <- paste0(".se_",       reference_level)

  if (!all(c(est_c, est_r) %in% names(te_diff))) {
    stop("Expected management columns not found after pivot.\n",
         "  Got: ", paste(names(te_diff), collapse = ", "), "\n",
         "  Looking for reference_level = '", reference_level,
         "', contrast_level = '", contrast_level, "'.")
  }

  te_diff <- te_diff |>
    dplyr::mutate(
      diff    = .data[[est_c]] - .data[[est_r]],
      diff_se = sqrt(.data[[se_c]]^2 + .data[[se_r]]^2)
    )

  # ---- 5. Mask grid cells far from observed data -------------------------
  obs <- data |>
    dplyr::select(dplyr::all_of(c(temp_var, climate_var))) |>
    na.omit()

  nn <- fields::rdist(te_diff[, c(temp_var, climate_var)], as.matrix(obs))
  te_diff$min_dist    <- apply(nn, 1, min)
  te_diff$diff_masked <- ifelse(
    te_diff$min_dist < quantile(te_diff$min_dist, mask_quantile, na.rm = TRUE),
    te_diff$diff, NA_real_
  )
  
  # this will gray out space representing combos of temperature and VPD that do
  # not occur in the dataset, so we can see where climate conditions are sparse
  # IRL and thus not make assumptions about the effect of management in conditions
  # that don't really occur (e.g., cool and dry conditions aren't common bc
  # cold air can't hold much vapor deficit)

  # ---- 6. Color logic ----------------------------------------------------
  # Goal: green always means "contrast-level favorable" (or simply "contrast
  # minus reference" for neutral framing). tan is the opposite.
  if (better_direction == "negative") {
    # Negative diff = contrast does more uptake = favorable -> green at low (AKA uptake) end
    fill_low  <- "forestgreen"
    fill_high <- "tan3"
    legend_sub <- paste0(contrast_level, " minus ", reference_level, "; Green = ", contrast_level,
                         " has greater uptake")
  } else if (better_direction == "positive") {
    # Positive diff = contrast does more = favorable -> green at HIGH (aka productivity) end (flip)
    fill_low  <- "tan3"
    fill_high <- "forestgreen"
    legend_sub <- paste0(contrast_level, " minus ", reference_level, "; Green = ", contrast_level,
                         " has higher productivity")
  } else { # neutral
    # No value loading; use standard diverging palette
    fill_low  <- "tan3"
    fill_high <- "forestgreen"
    legend_sub <- paste0(contrast_level, " minus ", reference_level, "; Green = ", contrast_level,
                         " has higher respiration rates")
  }

  # ---- 7. Auto-place corner labels and stress arrow ----------------------
  xr <- range(te_diff[[temp_var]],    na.rm = TRUE)
  yr <- range(te_diff[[climate_var]], na.rm = TRUE)
  x_lo <- xr[1] + 0.10 * diff(xr); x_hi <- xr[2] - 0.10 * diff(xr)
  y_lo <- yr[1] + 0.07 * diff(yr); y_hi <- yr[2] - 0.07 * diff(yr)

  default_corners_vpd <- list(
    bl = "cool & humid\n(low demand)",
    br = "warm & humid\n(productive)",
    tl = "cool & dry air\n(uncommon)",
    tr = "hot & dry air\n(atmospheric drought)"
  )
  default_corners_sm <- list(
    bl = "cool & dry\n(dormant)",
    br = "hot & dry\n(drought stress)",
    tl = "cool & wet\n(benign)",
    tr = "hot & wet\n(productive)"
  )

  corners <- if (layout == "temp_vpd") default_corners_vpd else default_corners_sm
  if (!is.null(corner_labels)) {
    for (k in intersect(names(corner_labels), names(corners))) {
      corners[[k]] <- corner_labels[[k]]
    }
  }

  corner_df <- data.frame(
    x = c(x_lo, x_hi, x_lo, x_hi),
    y = c(y_lo, y_lo, y_hi, y_hi),
    label = c(corners$bl, corners$br, corners$tl, corners$tr)
  )

  # Stress arrow: for VPD layout points from bottom-left (benign) to top-right
  # (hot+dry-air); for SM layout points from top-left (cool+wet, benign) to
  # # bottom-right (hot+dry-soil).
  # if (layout == "temp_vpd") {
  #   arrow_xy <- list(x    = x_lo,
  #                    y    = y_lo + 0.12 * diff(yr),
  #                    xend = x_hi,
  #                    yend = y_hi - 0.12 * diff(yr))
  # } else {
  #   arrow_xy <- list(x    = x_lo,
  #                    y    = y_hi - 0.12 * diff(yr),
  #                    xend = x_hi,
  #                    yend = y_lo + 0.12 * diff(yr))
  # }

  # ---- 8. Axis label defaults --------------------------------------------
  if (is.null(climate_units)) {
    climate_units <- if (layout == "temp_vpd") "VPD (kPa)" else "Soil moisture (% VWC)"
  }

  # ---- 9. Legend / title text --------------------------------------------
  flux_expr <- tryCatch(parse(text = flux_var)[[1]], error = function(e) flux_var)

  legend_title <- bquote(
    # .(contrast_level) - .(reference_level) ~ "\n" *
      Delta * .(flux_expr) ~ "(g C m"^-2 ~ "d"^-1 * ")"
  )


  plot_title <- bquote(
    "Contrast in" ~ .(flux_expr) ~
      "across temperature × " *
      .(if (layout == "temp_vpd") "VPD" else "soil moisture")
  )

  # ---- 10. Build the plot ------------------------------------------------
  p <- ggplot(te_diff,
              aes(x = .data[[temp_var]], y = .data[[climate_var]]/1000,
                  fill = diff_masked)) +
    geom_tile() +
    scale_fill_gradient2(
      low      = fill_low,
      mid      = "white",
      high     = fill_high,
      midpoint = 0,
      name     = legend_title,
      na.value = "grey92"
    )

  if (show_points) {
    p <- p + geom_point(data = data,
                        aes(x = .data[[temp_var]], y = .data[[climate_var]]/1000),
                        inherit.aes = FALSE,
                        colour = "black", alpha = 0.15, size = 0.7)
  }

  if (show_corner_labels) {
    p <- p + geom_label(data = corner_df,
                        aes(x = x, y = y/1000, label = label),
                        inherit.aes = FALSE,
                        color = "grey20", size = 2.6, label.size = 0,
                        fill = alpha("white", 0.8))
  }

  # if (show_arrow) {
    # p <- p + annotate("segment",
                      # x = arrow_xy$x,    y    = arrow_xy$y,
                      # xend = arrow_xy$xend, yend = arrow_xy$yend,
                      # arrow = arrow(type = "closed", length = unit(3, "mm")),
                      # color = "grey20", linewidth = 0.6)
  # }

  p <- p +
    labs(
      title    = plot_title,
      subtitle = legend_sub,
      x        = "Air temperature (°C)",
      y        = climate_units
    ) +
    theme_minimal() +
    theme(
      panel.grid      = element_blank(),
      legend.position = "right",
      plot.title      = element_text(size = 11),
      plot.subtitle   = element_text(size = 9, color = "grey30")
    )

  return(p)
}


# =============================================================================
# Example usage (uncomment to run)
# =============================================================================
#
# # ---- NEE: blue = organic uptakes more ----
# p_nee <- plot_te_difference(
#   model           = bam_NEE,
#   data            = data1,
#   flux_var        = "NEE",
#   temp_var        = "air_temperature",
#   climate_var     = "VPD",
#   better_direction = "negative",
#   layout          = "temp_vpd"
# )
# 
# p_nee
# #
# # # ---- GPP: blue = organic produces more (color flipped relative to NEE) ----
# p_gpp <- plot_te_difference(
#   model           = bam_GPP,
#   data            = data1,
#   flux_var        = "GPP",
#   temp_var        = "air_temperature",
#   climate_var     = "VPD",
#   better_direction = "positive",
#   layout          = "temp_vpd"
# )
# 
# p_gpp
# #
# # # ---- Reco (GAMM): neutral framing, no value loading on colors ----
# p_reco <- plot_te_difference(
#   model           = gamm_Reco_smoothed,
#   data            = data1,
#   flux_var        = "R[eco]",
#   temp_var        = "air_temperature_7",
#   climate_var     = "VPD_7",
#   better_direction = "neutral",
#   layout          = "temp_vpd"
# )
# 
# p_reco
# #
# # # ---- Chamber NEE temp x soil moisture layout ----
# p_chamber <- plot_te_difference(
#   model            = gam_simplified,
#   data             = chamber_daily_mngmnt,
#   flux_var         = "NEE",
#   temp_var         = "airtempC_mean",
#   climate_var      = "soilm_perc_mean",
#   better_direction = "negative",
#   layout           = "temp_sm",
#   climate_units    = "Soil moisture (% VWC)"
# )
# 
# p_chamber
# #
# # # ---- Combine all three EC plots with patchwork for a paper figure ----
# library(patchwork)
# (p_nee / p_gpp / p_reco) + plot_layout(guides = "keep")
# 
# (p_nee + theme(legend.position = "bottom")) + 
#   (p_gpp + ylab(NULL) + theme(legend.position = "bottom")) + 
#   (p_reco + ylab(NULL) + theme(legend.position = "bottom")) + 
#   plot_layout(guides = "keep")
# #
