# =============================================================================
# make_pub_tables.R
#
# Extracts publication-ready summary tables from mgcv gam / bam / gamm models.
# For each model it returns:
#   * main_parametric : filtered parametric coefficients (significant terms
#                       + always-keep terms; intended for main paper)
#   * supp_parametric : full parametric coefficient table (for supplement)
#   * smooth_terms    : smooth-term significance (one table, used everywhere)
#   * model_stats     : R^2, deviance explained, n, family, link, REML
#
# Returns tibbles only. No rendering. Use gt, flextable, kableExtra, or just
# write_csv() downstream to produce final formatted output.
#
# Author: Bee Forbes
# =============================================================================

library(dplyr)
library(tibble)
library(broom)
library(stringr)
library(purrr)
library(readr)


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

#' Add significance-star column based on p.value
.add_sig_stars <- function(df, p_col = "p.value") {
  df %>%
    dplyr::mutate(
      sig = dplyr::case_when(
        .data[[p_col]] < 0.001 ~ "***",
        .data[[p_col]] < 0.01  ~ "**",
        .data[[p_col]] < 0.05  ~ "*",
        .data[[p_col]] < 0.10  ~ ".",
        TRUE                   ~ ""
      )
    )
}

#' Clean ugly mgcv term names into human-readable labels
#' Returns the *original* in a `term_raw` column and the cleaned label in `term`
.clean_terms <- function(df,
                         management_col_pattern = "management",
                         crop_stage_col_pattern = "crop_stage_simple") {
  df %>%
    dplyr::mutate(
      term_raw = term,
      term = term %>%
        stringr::str_replace_all(paste0(management_col_pattern, "organic"),
                                 "Organic") %>%
        stringr::str_replace_all(paste0(management_col_pattern, "conventional"),
                                 "Conventional") %>%
        stringr::str_replace_all(crop_stage_col_pattern, "") %>%
        stringr::str_replace_all("as\\.factor\\(management\\)", "") %>%
        stringr::str_replace_all(":", " × ") %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_replace_all("grain fill", "Grain fill") %>%
        stringr::str_replace_all("\\bearly\\b", "Early growth") %>%
        stringr::str_replace_all("\\bvegetative\\b", "Vegetative") %>%
        stringr::str_replace_all("\\breproductive\\b", "Reproductive") %>%
        stringr::str_replace_all("\\bmature\\b", "Mature") %>%
        stringr::str_replace_all("\\bdormant\\b", "Dormant") %>%
        stringr::str_replace_all("\\bfallow\\b", "Fallow") %>%
        stringr::str_replace_all("\\bregrowth\\b", "Regrowth") %>%
        stringr::str_replace_all("\\(Intercept\\)", "Intercept") %>%
        stringr::str_squish()
    )
}

#' Unwrap gamm objects to get the $gam component
.unwrap_gamm <- function(model) {
  if (!inherits(model, "gam") && is.list(model) && "gam" %in% names(model)) {
    return(model[["gam"]])
  }
  model
}


# -----------------------------------------------------------------------------
# Main function: tables for a single model
# -----------------------------------------------------------------------------

#' Build pub-ready tables for one model
#'
#' @param model       A fitted gam / bam / gamm object.
#' @param model_name  Character label used in the model_stats row.
#' @param keep_terms  Character vector of term *labels* (post-cleaning) that
#'                    should always appear in main_parametric even if n.s.
#'                    Default keeps the management main effect and intercept.
#' @param sig_threshold p threshold for inclusion in main_parametric.
#' @param decimals    Decimal places for numeric rounding.
#' @param r2_note     Optional string to attach to model_stats (e.g. for
#'                    flagging the GAMM R^2 issue). NA to skip.
#' @return List of four tibbles.
make_pub_tables <- function(model,
                            model_name      = "model",
                            keep_terms      = c("Intercept", "Organic"),
                            sig_threshold   = 0.10,
                            decimals        = 3,
                            r2_note         = NA_character_) {

  model <- .unwrap_gamm(model)
  sm    <- summary(model)

  # ---- Parametric coefficients --------------------------------------------
  para_full <- broom::tidy(model, parametric = TRUE) %>%
    .clean_terms() %>%
    .add_sig_stars(p_col = "p.value") %>%
    dplyr::select(term, term_raw,
                  Estimate  = estimate,
                  SE        = std.error,
                  t         = statistic,
                  p         = p.value,
                  sig) %>%
    dplyr::mutate(dplyr::across(c(Estimate, SE, t, p),
                                \(x) round(x, decimals)))

  always_keep_mask <- para_full$term %in% keep_terms
  sig_mask         <- !is.na(para_full$p) & para_full$p < sig_threshold

  para_main <- para_full %>%
    dplyr::filter(always_keep_mask | sig_mask)

  # ---- Smooth terms -------------------------------------------------------
  smooth_full <- broom::tidy(model, parametric = FALSE) %>%
    .clean_terms() %>%
    .add_sig_stars(p_col = "p.value") %>%
    dplyr::select(term, term_raw,
                  EDF       = edf,
                  Ref.df    = ref.df,
                  F         = statistic,
                  p         = p.value,
                  sig) %>%
    dplyr::mutate(dplyr::across(c(EDF, `Ref.df`, F, p),
                                \(x) round(x, decimals)))

  # ---- Model fit statistics -----------------------------------------------
  stats <- tibble::tibble(
    model               = model_name,
    family              = sm$family$family,
    link                = sm$family$link,
    n                   = sm$n,
    R2_adj              = if (!is.null(sm$r.sq)) round(sm$r.sq, decimals) else NA_real_,
    deviance_explained  = if (!is.null(sm$dev.expl)) round(sm$dev.expl * 100, 1) else NA_real_,
    REML_or_GCV         = if (!is.null(sm$sp.criterion)) round(sm$sp.criterion, 2) else NA_real_,
    note                = r2_note
  )

  list(
    main_parametric = para_main,
    supp_parametric = para_full,
    smooth_terms    = smooth_full,
    model_stats     = stats
  )
}


# -----------------------------------------------------------------------------
# Convenience: build tables for several models and optionally write CSVs
# -----------------------------------------------------------------------------

#' Build pub-ready tables for multiple models
#'
#' @param models           Named list of fitted gam / bam / gamm objects.
#' @param keep_terms       Vector of term labels always kept in main_parametric.
#' @param sig_threshold    p threshold for main_parametric.
#' @param decimals         Decimal rounding.
#' @param r2_notes         Optional named character vector of per-model R^2
#'                         notes (e.g. for the gamm Reco model). Names should
#'                         match `names(models)`.
#' @param output_dir       If non-NULL, CSVs are written to this directory:
#'                         one main_parametric, one supp_parametric, one
#'                         smooth_terms file per model, plus a combined
#'                         all_model_stats.csv.
#'
#' @return Named list of per-model table lists, plus a combined `model_stats`
#'         tibble with one row per model.
make_all_pub_tables <- function(models,
                                keep_terms     = c("Intercept", "Organic"),
                                sig_threshold  = 0.10,
                                decimals       = 3,
                                r2_notes       = NULL,
                                output_dir     = NULL) {

  if (is.null(names(models))) {
    stop("`models` must be a named list (names become file prefixes).")
  }

  per_model <- purrr::map2(
    .x = models, .y = names(models),
    .f = \(mod, nm) make_pub_tables(
      mod,
      model_name    = nm,
      keep_terms    = keep_terms,
      sig_threshold = sig_threshold,
      decimals      = decimals,
      r2_note       = if (!is.null(r2_notes) && nm %in% names(r2_notes))
                        r2_notes[[nm]] else NA_character_
    )
  )

  combined_stats <- purrr::map_dfr(per_model, "model_stats")

  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    for (nm in names(per_model)) {
      readr::write_csv(per_model[[nm]]$main_parametric,
                       file.path(output_dir, paste0(nm, "_main_parametric.csv")))
      readr::write_csv(per_model[[nm]]$supp_parametric,
                       file.path(output_dir, paste0(nm, "_supp_parametric.csv")))
      readr::write_csv(per_model[[nm]]$smooth_terms,
                       file.path(output_dir, paste0(nm, "_smooth_terms.csv")))
    }
    readr::write_csv(combined_stats,
                     file.path(output_dir, "all_model_stats.csv"))
    message("Wrote tables to ", normalizePath(output_dir))
  }

  list(
    per_model = per_model,
    combined_stats = combined_stats
  )
}


# =============================================================================
# Example usage (uncomment to run)
# =============================================================================
#
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
    Reco = paste0("GAMM R^2 inflated by ARMA(1,1) correlation structure ",
                  "(fitted-value variance > observed). Use CCC = 0.631 ",
                  "(95% CI 0.595-0.664) instead.")
  ),
  output_dir = "tables/"
)

# Access any table directly:
results$per_model$NEE$main_parametric
results$per_model$NEE$smooth_terms
results$combined_stats
#
# # Or pass to gt / flextable for rendering later:
# gt::gt(results$per_model$NEE$main_parametric)
