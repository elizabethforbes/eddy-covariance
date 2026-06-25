# -----------------------------------------------------------------------------
# Summary tables
# -----------------------------------------------------------------------------

library(gt)
library(tidyverse)

################################################################################
# NEE
################################################################################

# -- Parametric coefficients (fixed effects) ----------------------------------
# Note: tbl_regression() labels by variable name, not coefficient name, so
# we use broom::tidy() + gt() directly for full control over interaction labels.
tbl_nee_parametric <- broom::tidy(bam_NEE, parametric = TRUE) %>%
  dplyr::mutate(
    term = dplyr::case_when(
      term == "(Intercept)"                                      ~ "Intercept (conventional, mature)",
      term == "managementorganic"                                ~ "Organic vs. Conventional",
      term == "crop_stage_simpledormant"                         ~ "Dormant",
      term == "crop_stage_simpleearly"                           ~ "Early growth",
      term == "crop_stage_simplegrain_fill"                      ~ "Grain fill",
      term == "crop_stage_simplereproductive"                    ~ "Reproductive",
      term == "crop_stage_simplevegetative"                      ~ "Vegetative",
      term == "crop_stage_simplefallow"                          ~ "Fallow",
      term == "managementorganic:crop_stage_simpledormant"       ~ "Organic \u00d7 Dormant",
      term == "managementorganic:crop_stage_simpleearly"         ~ "Organic \u00d7 Early growth",      
      term == "managementorganic:crop_stage_simplefallow"         ~ "Organic \u00d7 Fallow",      
      term == "managementorganic:crop_stage_simplegrain_fill"    ~ "Organic \u00d7 Grain fill",
      term == "managementorganic:crop_stage_simplereproductive"  ~ "Organic \u00d7 Reproductive",
      term == "managementorganic:crop_stage_simplevegetative"    ~ "Organic \u00d7 Vegetative",
      TRUE ~ term
    ),
    sig = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE            ~ ""
    )
  ) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value, sig) %>%
  gt() %>%
  cols_label(
    term      = "Term",
    estimate  = "Estimate",
    std.error = "Std. Error",
    statistic = "t",
    p.value   = "p-value",
    sig       = ""
  ) %>%
  fmt_number(columns = c(estimate, std.error, statistic), decimals = 3) %>%
  fmt(
    columns = p.value,
    fns = function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  ) %>%
  tab_header(
    title    = md("**Table 1. Parametric coefficients: NEE GAM (bam_NEE)**"),
    subtitle = md("Reference: conventional management, mature crop stage")
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = everything(), rows = p.value < 0.05)
  )

# -- Smooth terms (approximate significance) ----------------------------------
tbl_nee_smooths <- broom::tidy(bam_NEE, parametric = FALSE) %>%
  dplyr::mutate(
    term = dplyr::case_when(
      grepl("s\\(year_f\\).*conventional", term) ~ "s(year): conventional",
      grepl("s\\(year_f\\).*organic",      term) ~ "s(year): organic",
      grepl("s\\(doy\\)",                term) ~ "s(doy) \u2014 shared seasonal cycle",
      grepl("te.*conventional",          term) ~ "te(Temp \u00d7 VPD): conventional",
      grepl("te.*organic",               term) ~ "te(Temp \u00d7 VPD): organic",
      grepl("days_since.*conventional",  term) ~ "s(days since tillage): conventional",
      grepl("days_since.*organic",       term) ~ "s(days since tillage): organic",
      TRUE ~ term
    )
  ) %>%
  dplyr::select(term, edf, ref.df, statistic, p.value) %>%
  gt() %>%
  cols_label(
    term      = "Smooth term",
    edf       = "EDF",
    ref.df    = "Ref. df",
    statistic = "F",
    p.value   = "p-value"
  ) %>%
  fmt_number(columns = c(edf, ref.df, statistic), decimals = 2) %>%
  fmt(
    columns = p.value,
    fns = function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  ) %>%
  tab_header(
    title    = md("**Table 2. Smooth terms: NEE GAM (bam_NEE)**"),
    subtitle = md("Approximate significance; EDF = effective degrees of freedom")
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = p.value, rows = p.value < 0.05)
  )

tbl_nee_parametric
tbl_nee_smooths

################################################################################
# GPP
################################################################################
tbl_gpp_parametric <- broom::tidy(bam_GPP, parametric = TRUE) %>%
  dplyr::mutate(
    term = dplyr::case_when(
      term == "(Intercept)"                                      ~ "Intercept (conventional, mature)",
      term == "managementorganic"                                ~ "Organic vs. Conventional",
      term == "crop_stage_simpledormant"                         ~ "Dormant",
      term == "crop_stage_simpleearly"                           ~ "Early growth",
      term == "crop_stage_simplegrain_fill"                      ~ "Grain fill",
      term == "crop_stage_simplereproductive"                    ~ "Reproductive",
      term == "crop_stage_simplevegetative"                      ~ "Vegetative",
      term == "crop_stage_simplefallow"                          ~ "Fallow",
      term == "managementorganic:crop_stage_simpledormant"       ~ "Organic \u00d7 Dormant",
      term == "managementorganic:crop_stage_simpleearly"         ~ "Organic \u00d7 Early growth",      
      term == "managementorganic:crop_stage_simplefallow"         ~ "Organic \u00d7 Fallow",      
      term == "managementorganic:crop_stage_simplegrain_fill"    ~ "Organic \u00d7 Grain fill",
      term == "managementorganic:crop_stage_simplereproductive"  ~ "Organic \u00d7 Reproductive",
      term == "managementorganic:crop_stage_simplevegetative"    ~ "Organic \u00d7 Vegetative",
      TRUE ~ term
    ),
    sig = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE            ~ ""
    )
  ) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value, sig) %>%
  gt() %>%
  cols_label(
    term      = "Term",
    estimate  = "Estimate",
    std.error = "Std. Error",
    statistic = "t",
    p.value   = "p-value",
    sig       = ""
  ) %>%
  fmt_number(columns = c(estimate, std.error, statistic), decimals = 3) %>%
  fmt(
    columns = p.value,
    fns = function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  ) %>%
  tab_header(
    title    = md("**Table 1. Parametric coefficients: GPP GAM (bam_GPP)**"),
    subtitle = md("Reference: conventional management, mature crop stage")
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = everything(), rows = p.value < 0.05)
  )

# -- Smooth terms (approximate significance) ----------------------------------
tbl_gpp_smooths <- broom::tidy(bam_GPP, parametric = FALSE) %>%
  dplyr::mutate(
    term = dplyr::case_when(
      grepl("s\\(year_f\\).*conventional", term) ~ "s(year): conventional",
      grepl("s\\(year_f\\).*organic",      term) ~ "s(year): organic",
      grepl("s\\(doy\\)",                term) ~ "s(doy) \u2014 shared seasonal cycle",
      grepl("te.*conventional",          term) ~ "te(Temp \u00d7 VPD): conventional",
      grepl("te.*organic",               term) ~ "te(Temp \u00d7 VPD): organic",
      grepl("days_since.*conventional",  term) ~ "s(days since tillage): conventional",
      grepl("days_since.*organic",       term) ~ "s(days since tillage): organic",
      TRUE ~ term
    )
  ) %>%
  dplyr::select(term, edf, ref.df, statistic, p.value) %>%
  gt() %>%
  cols_label(
    term      = "Smooth term",
    edf       = "EDF",
    ref.df    = "Ref. df",
    statistic = "F",
    p.value   = "p-value"
  ) %>%
  fmt_number(columns = c(edf, ref.df, statistic), decimals = 2) %>%
  fmt(
    columns = p.value,
    fns = function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  ) %>%
  tab_header(
    title    = md("**Table 2. Smooth terms: GPP GAM (bam_GPP)**"),
    subtitle = md("Approximate significance; EDF = effective degrees of freedom")
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = p.value, rows = p.value < 0.05)
  )

tbl_gpp_parametric
tbl_gpp_smooths

################################################################################
# Reco
################################################################################
tbl_reco_parametric <- broom::tidy(gamm_Reco_smoothed$gam, parametric = TRUE) %>%
  dplyr::mutate(
    term = dplyr::case_when(
      term == "(Intercept)"                                      ~ "Intercept (conventional, mature)",
      term == "managementorganic"                                ~ "Organic vs. Conventional",
      term == "crop_stage_simpledormant"                         ~ "Dormant",
      term == "crop_stage_simpleearly"                           ~ "Early growth",
      term == "crop_stage_simplegrain_fill"                      ~ "Grain fill",
      term == "crop_stage_simplereproductive"                    ~ "Reproductive",
      term == "crop_stage_simplevegetative"                      ~ "Vegetative",
      term == "crop_stage_simplefallow"                          ~ "Fallow",
      term == "managementorganic:crop_stage_simpledormant"       ~ "Organic \u00d7 Dormant",
      term == "managementorganic:crop_stage_simpleearly"         ~ "Organic \u00d7 Early growth",      
      term == "managementorganic:crop_stage_simplefallow"         ~ "Organic \u00d7 Fallow",      
      term == "managementorganic:crop_stage_simplegrain_fill"    ~ "Organic \u00d7 Grain fill",
      term == "managementorganic:crop_stage_simplereproductive"  ~ "Organic \u00d7 Reproductive",
      term == "managementorganic:crop_stage_simplevegetative"    ~ "Organic \u00d7 Vegetative",
      TRUE ~ term
    ),
    sig = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE            ~ ""
    )
  ) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value, sig) %>%
  gt() %>%
  cols_label(
    term      = "Term",
    estimate  = "Estimate",
    std.error = "Std. Error",
    statistic = "t",
    p.value   = "p-value",
    sig       = ""
  ) %>%
  fmt_number(columns = c(estimate, std.error, statistic), decimals = 3) %>%
  fmt(
    columns = p.value,
    fns = function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  ) %>%
  tab_header(
    title    = md("**Table 1. Parametric coefficients: Reco GAMM with 7day averaged climate variables**"),
    subtitle = md("Reference: conventional management, mature crop stage")
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = everything(), rows = p.value < 0.05)
  )

# -- Smooth terms (approximate significance) ----------------------------------
tbl_reco_smooths <- broom::tidy(gamm_Reco_smoothed$gam, parametric = FALSE) %>%
  dplyr::mutate(
    term = dplyr::case_when(
      grepl("s\\(year_f\\).*conventional", term) ~ "s(year): conventional",
      grepl("s\\(year_f\\).*organic",      term) ~ "s(year): organic",
      grepl("s\\(doy\\)",                term) ~ "s(doy) \u2014 shared seasonal cycle",
      grepl("te.*conventional",          term) ~ "te(Temp \u00d7 VPD): conventional",
      grepl("te.*organic",               term) ~ "te(Temp \u00d7 VPD): organic",
      grepl("days_since.*conventional",  term) ~ "s(days since tillage): conventional",
      grepl("days_since.*organic",       term) ~ "s(days since tillage): organic",
      TRUE ~ term
    )
  ) %>%
  dplyr::select(term, edf, ref.df, statistic, p.value) %>%
  gt() %>%
  cols_label(
    term      = "Smooth term",
    edf       = "EDF",
    ref.df    = "Ref. df",
    statistic = "F",
    p.value   = "p-value"
  ) %>%
  fmt_number(columns = c(edf, ref.df, statistic), decimals = 2) %>%
  fmt(
    columns = p.value,
    fns = function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  ) %>%
  tab_header(
    title    = md("**Table 2. Smooth terms: Reco GAMM with 7day averaged climate variables**"),
    subtitle = md("Approximate significance; EDF = effective degrees of freedom")
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = p.value, rows = p.value < 0.05)
  )

tbl_reco_parametric
tbl_reco_smooths
