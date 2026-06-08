# ============================================================================
# TILLAGE SMOOTH - s(days_since_tillage) by management
# ============================================================================

# Color palette for management
management_colors <- c(
  "Conventional" = "tomato",  # Red
  "Organic" = "#E69F00"        # Green
)

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
    scale_fill_manual(values = management_colors, guide = "none") +
    facet_wrap(~ Management, ncol = 2, scales = "free_y") +
    labs(
      # title    = paste0("Effect of days since tillage on ", flux_var),
      # subtitle = "Conventional: non-significant (p = 0.962) | Organic: nonlinear, p = 0.040 (EDF = 3.03)",
      x        = "Days since tillage",
      y        = bquote("Partial effect on" ~ .(flux_var) ~ "(µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ ")")
    ) +
    theme_flux() +
    theme(strip.text = element_text(face = "bold"))
  
  return(p)
}

t1 <- plot_tillage_smooth(bam_NEE, data = data1, flux_var = "NEE")
t2 <- plot_tillage_smooth(bam_GPP, data = data1, flux_var = "GPP")
t3 <- plot_tillage_smooth(gamm_Reco_smoothed, data = data1, flux_var = "Reco")

library(patchwork)
t1+t2+t3+plot_layout(ncol = 2)
