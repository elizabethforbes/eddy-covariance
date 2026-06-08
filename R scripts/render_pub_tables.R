# =============================================================================
# render_pub_tables.R
#
# Renders publication-ready Word tables from the tibbles produced by
# make_pub_tables.R, using the flextable package.
#
# Produces two .docx outputs:
#   - tables_main_paper.docx   : selected coefficients + smooth terms per model
#                                + a combined fit-statistics table
#   - tables_supplement.docx   : full coefficient tables for all three models
#
# Each model gets two flextables (parametric + smooth), each placed on its own
# page in the resulting docx. Significant terms (p < 0.05) are bolded; bold can
# be tuned with `sig_bold_threshold`.
#
# Dependencies: flextable, officer, dplyr
#
# Author: Bee Forbes
# =============================================================================

library(flextable)
library(officer)
library(dplyr)

# Source the data layer (adjust path as needed)
# source("make_pub_tables.R")


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

#' Build one or more footer lines describing a model's fit statistics
.fit_stats_footer <- function(stats_row) {
  parts <- c(
    paste0("Family: ", stats_row$family, " (", stats_row$link, " link)"),
    paste0("n = ", stats_row$n)
  )

  if (!is.na(stats_row$R2_adj)) {
    parts <- c(parts, paste0("R² (adj.) = ", stats_row$R2_adj))
  }
  if (!is.na(stats_row$deviance_explained)) {
    parts <- c(parts, paste0("Deviance explained = ",
                             stats_row$deviance_explained, "%"))
  }
  if (!is.na(stats_row$REML_or_GCV)) {
    parts <- c(parts, paste0("REML / GCV = ", stats_row$REML_or_GCV))
  }

  base_line <- paste(parts, collapse = "; ")

  if (!is.na(stats_row$note) && nzchar(stats_row$note)) {
    c(base_line, paste0("Note: ", stats_row$note))
  } else {
    base_line
  }
}


# -----------------------------------------------------------------------------
# Render functions: parametric, smooth, combined-stats
# -----------------------------------------------------------------------------

#' Render a parametric coefficients tibble as a publication-ready flextable
render_parametric_ft <- function(df,
                                 caption = NULL,
                                 fit_stats_lines = NULL,
                                 sig_bold_threshold = 0.05) {

  if ("term_raw" %in% names(df)) df <- dplyr::select(df, -term_raw)

  ft <- flextable::flextable(df) |>
    flextable::theme_booktabs() |>
    flextable::set_header_labels(
      term     = "Term",
      Estimate = "β",      # Greek beta
      SE       = "SE",
      t        = "t",
      p        = "p-value",
      sig      = ""
    ) |>
    flextable::align(part = "all", j = "term", align = "left") |>
    flextable::align(part = "all",
                     j = c("Estimate", "SE", "t", "p"),
                     align = "right") |>
    flextable::align(part = "all", j = "sig", align = "center")

  # Bold significant rows
  sig_idx <- which(!is.na(df$p) & df$p < sig_bold_threshold)
  if (length(sig_idx) > 0) {
    ft <- flextable::bold(ft, i = sig_idx, part = "body")
  }

  if (!is.null(caption)) {
    ft <- flextable::set_caption(ft, caption)
  }

  if (!is.null(fit_stats_lines)) {
    ft <- flextable::add_footer_lines(ft, values = fit_stats_lines)
  }

  ft |>
    flextable::set_table_properties(width = 1, layout = "autofit") |>
    flextable::fontsize(size = 10, part = "body") |>
    flextable::fontsize(size = 10, part = "header")
}


#' Render a smooth-terms tibble as a publication-ready flextable
render_smooth_ft <- function(df,
                             caption = NULL,
                             fit_stats_lines = NULL,
                             sig_bold_threshold = 0.05) {

  if ("term_raw" %in% names(df)) df <- dplyr::select(df, -term_raw)

  ft <- flextable::flextable(df) |>
    flextable::theme_booktabs() |>
    flextable::set_header_labels(
      term     = "Smooth term",
      EDF      = "EDF",
      Ref.df   = "Ref. df",
      F        = "F",
      p        = "p-value",
      sig      = ""
    ) |>
    flextable::align(part = "all", j = "term", align = "left") |>
    flextable::align(part = "all",
                     j = c("EDF", "Ref.df", "F", "p"),
                     align = "right") |>
    flextable::align(part = "all", j = "sig", align = "center")

  sig_idx <- which(!is.na(df$p) & df$p < sig_bold_threshold)
  if (length(sig_idx) > 0) {
    ft <- flextable::bold(ft, i = sig_idx, part = "body")
  }

  if (!is.null(caption)) {
    ft <- flextable::set_caption(ft, caption)
  }

  if (!is.null(fit_stats_lines)) {
    ft <- flextable::add_footer_lines(ft, values = fit_stats_lines)
  }

  ft |>
    flextable::set_table_properties(width = 1, layout = "autofit") |>
    flextable::fontsize(size = 10, part = "body") |>
    flextable::fontsize(size = 10, part = "header")
}


#' Render the combined model-fit-statistics table
render_combined_stats_ft <- function(combined_stats, caption = NULL) {

  ft <- flextable::flextable(combined_stats) |>
    flextable::theme_booktabs() |>
    flextable::set_header_labels(
      model              = "Model",
      family             = "Family",
      link               = "Link",
      n                  = "n",
      R2_adj             = "R² (adj.)",
      deviance_explained = "Dev. expl. (%)",
      REML_or_GCV        = "REML / GCV",
      note               = "Note"
    ) |>
    flextable::align(part = "all", align = "left") |>
    flextable::set_table_properties(width = 1, layout = "autofit") |>
    flextable::fontsize(size = 9, part = "all")

  if (!is.null(caption)) {
    ft <- flextable::set_caption(ft, caption)
  }

  ft
}


# -----------------------------------------------------------------------------
# Top-level driver: produce both docx files
# -----------------------------------------------------------------------------

#' Render and save all main-paper and supplement tables to .docx files
#'
#' @param results            List returned by make_all_pub_tables().
#' @param output_main_docx   Path for main-paper tables.
#' @param output_supp_docx   Path for supplement tables.
#' @param sig_bold_threshold p threshold below which rows are bolded.
#' @return Invisible list with both sets of flextables.
render_all_to_docx <- function(results,
                               output_main_docx   = "tables_main_paper.docx",
                               output_supp_docx   = "tables_supplement.docx",
                               sig_bold_threshold = 0.05) {

  per_model      <- results$per_model
  combined_stats <- results$combined_stats
  model_names    <- names(per_model)

  # ---- Main paper tables ----
  main_tables <- list()

  for (i in seq_along(model_names)) {
    nm        <- model_names[i]
    fit_lines <- .fit_stats_footer(per_model[[nm]]$model_stats)

    main_tables[[paste0(nm, "_parametric")]] <- render_parametric_ft(
      per_model[[nm]]$main_parametric,
      caption = paste0("Table ", i,
                       "A. Selected parametric coefficients: ", nm,
                       " (reference: conventional management, mature crop stage)."),
      fit_stats_lines    = fit_lines,
      sig_bold_threshold = sig_bold_threshold
    )

    main_tables[[paste0(nm, "_smooth")]] <- render_smooth_ft(
      per_model[[nm]]$smooth_terms,
      caption = paste0("Table ", i, "B. Smooth-term significance: ", nm, "."),
      sig_bold_threshold = sig_bold_threshold
    )
  }

  main_tables[["combined_stats"]] <- render_combined_stats_ft(
    combined_stats,
    caption = "Table X. Model fit summary across flux models."
  )

  do.call(flextable::save_as_docx,
          c(main_tables, list(path = output_main_docx)))
  message("Wrote ", normalizePath(output_main_docx))

  # ---- Supplement tables ----
  supp_tables <- list()

  for (i in seq_along(model_names)) {
    nm        <- model_names[i]
    fit_lines <- .fit_stats_footer(per_model[[nm]]$model_stats)

    supp_tables[[paste0(nm, "_supp_parametric")]] <- render_parametric_ft(
      per_model[[nm]]$supp_parametric,
      caption = paste0("Table S", i,
                       "A. Full parametric coefficients: ", nm,
                       " (reference: conventional management, mature crop stage)."),
      fit_stats_lines    = fit_lines,
      sig_bold_threshold = sig_bold_threshold
    )

    supp_tables[[paste0(nm, "_supp_smooth")]] <- render_smooth_ft(
      per_model[[nm]]$smooth_terms,
      caption = paste0("Table S", i, "B. Full smooth-term significance: ", nm, "."),
      sig_bold_threshold = sig_bold_threshold
    )
  }

  do.call(flextable::save_as_docx,
          c(supp_tables, list(path = output_supp_docx)))
  message("Wrote ", normalizePath(output_supp_docx))

  invisible(list(main = main_tables, supp = supp_tables))
}


# =============================================================================
# Example usage (uncomment to run)
# =============================================================================
#
# source("make_pub_tables.R")
#
# # Build all tibbles
results <- make_all_pub_tables(
  models = list(
    NEE  = bam_NEE,
    GPP  = bam_GPP,
    Reco = gamm_Reco_smoothed
  ),
  keep_terms     = c("Intercept", "Organic"),
  sig_threshold  = 0.10,
  decimals       = 3,
  r2_notes       = c(
    Reco = paste0("GAMM R² inflated by ARMA(1,1) correlation structure ",
                  "(fitted-value variance > observed); use CCC = 0.631 ",
                  "(95% CI 0.595-0.664) instead.")
  )
)
#
# # Render to two .docx files
render_all_to_docx(
  results,
  output_main_docx = "tables/tables_main_paper.docx",
  output_supp_docx = "tables/tables_supplement.docx",
  sig_bold_threshold = 0.05
)
#
# # Or render a single table on its own (e.g., to drop into a slide deck):

######
# NEE
######
ft_nee_main <- render_parametric_ft(
  results$per_model$NEE$main_parametric,
  caption = "Selected coefficients: NEE",
  fit_stats_lines = .fit_stats_footer(results$per_model$NEE$model_stats)
)
ft_nee_main   # preview in RStudio Viewer, make 450 x 415
# # flextable::save_as_image(ft_nee_main, path = "ft_nee_main.png")

ft_nee_smooth <- render_smooth_ft(
  results$per_model$NEE$smooth_terms,
  caption = "Selected smooth terms: NEE",
  fit_stats_lines = .fit_stats_footer(results$per_model$NEE$model_stats)
)
ft_nee_smooth # make 450 x 415

######
# GPP
######
ft_gpp_main <- render_parametric_ft(
  results$per_model$GPP$main_parametric,
  caption = "Selected coefficients: GPP",
  fit_stats_lines = .fit_stats_footer(results$per_model$GPP$model_stats)
)
ft_gpp_main   # preview in RStudio Viewer, make 450 x 415
# # flextable::save_as_image(ft_nee_main, path = "ft_nee_main.png")

ft_gpp_smooth <- render_smooth_ft(
  results$per_model$GPP$smooth_terms,
  caption = "Selected smooth terms: GPP",
  fit_stats_lines = .fit_stats_footer(results$per_model$GPP$model_stats)
)
ft_gpp_smooth # make 450 x 415

######
# resp
######
ft_reco_main <- render_parametric_ft(
  results$per_model$Reco$main_parametric,
  caption = "Selected coefficients: Respiration",
  fit_stats_lines = .fit_stats_footer(results$per_model$NEE$model_stats)
)
ft_reco_main   # preview in RStudio Viewer, make 450 x 415
# # flextable::save_as_image(ft_nee_main, path = "ft_nee_main.png")

ft_reco_smooth <- render_smooth_ft(
  results$per_model$Reco$smooth_terms,
  caption = "Selected smooth terms: Respiration",
  fit_stats_lines = .fit_stats_footer(results$per_model$Reco$model_stats)
)
ft_reco_smooth # make 450 x 415

######
# all
######

ft_all_main <- render_combined_stats_ft(results$combined_stats)
ft_all_main

#########
# chamber
#########

results_ch <- make_all_pub_tables(
  models = list(
    NEE  = gam_NEE_drought_wtd,
    # GPP  = bam_GPP,
    Reco = gam_Rs_drought
  ),
  keep_terms     = c("Intercept", "Organic"),
  sig_threshold  = 0.10,
  decimals       = 3,
  r2_notes       = c(
    Reco = paste0("GAMM R² inflated by ARMA(1,1) correlation structure ",
                  "(fitted-value variance > observed); use CCC = 0.631 ",
                  "(95% CI 0.595-0.664) instead.")
  )
)

# NEE
ft_nee_ch <- render_parametric_ft(
  results_ch$per_model$NEE$main_parametric,
  caption = "Selected coefficients: NEE, chamber",
  fit_stats_lines = .fit_stats_footer(results_ch$per_model$NEE$model_stats)
)
ft_nee_ch   # preview in RStudio Viewer, make 450 x 415
# # flextable::save_as_image(ft_nee_main, path = "ft_nee_main.png")

ft_nee_smooth <- render_smooth_ft(
  results_ch$per_model$NEE$smooth_terms,
  caption = "Selected smooth terms: NEE, chamber",
  fit_stats_lines = .fit_stats_footer(results_ch$per_model$NEE$model_stats)
)
ft_nee_smooth # make 450 x 415

# respiration:
ft_reco_ch <- render_parametric_ft(
  results_ch$per_model$Reco$main_parametric,
  caption = "Selected coefficients: respiration, chamber",
  fit_stats_lines = .fit_stats_footer(results_ch$per_model$Reco$model_stats)
)
ft_reco_ch   # preview in RStudio Viewer, make 450 x 415
# # flextable::save_as_image(ft_nee_main, path = "ft_nee_main.png")

ft_reco_smooth <- render_smooth_ft(
  results_ch$per_model$Reco$smooth_terms,
  caption = "Selected smooth terms: respiration, chamber",
  fit_stats_lines = .fit_stats_footer(results_ch$per_model$Reco$model_stats)
)
ft_reco_smooth # make 450 x 415
