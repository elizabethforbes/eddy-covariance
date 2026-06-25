library(broom)
library(gt)   # or kableExtra, flextable, etc.

################################################################################
### NEE
################################################################################

# Parametric coefficients
para_tbl <- broom::tidy(bam_NEE, parametric = TRUE) %>%
  dplyr::filter(p.value < 0.10 | term == "managementorganic") %>%
  # ^ keep significant + the management main effect even if n.s.
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

# Smooth terms  
smooth_tbl <- broom::tidy(bam_NEE, parametric = FALSE) %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

# Render with gt -- NEE
gt::gt(para_tbl) %>%
  gt::tab_header(title = "Selected parametric coefficients: NEE") %>%
  gt::fmt_number(decimals = 3)
gt::gt(smooth_tbl) %>% 
  gt::tab_header(title = "Selected smooth coefficients: NEE") %>% 
  gt::fmt_number(decimals = 3)

################################################################################
### GPP
################################################################################

# Parametric coefficients
para_tbl <- broom::tidy(bam_GPP, parametric = TRUE) %>%
  dplyr::filter(p.value < 0.10 | term == "managementorganic") %>%
  # ^ keep significant + the management main effect even if n.s.
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

# Smooth terms  
smooth_tbl <- broom::tidy(bam_GPP, parametric = FALSE) %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

# Render with gt -- GPP
gt::gt(para_tbl) %>%
  gt::tab_header(title = "Selected parametric coefficients: GPP") %>%
  gt::fmt_number(decimals = 3)
gt::gt(smooth_tbl) %>% 
  gt::tab_header(title = "Selected smooth coefficients: GPP") %>% 
  gt::fmt_number(decimals = 3)

################################################################################
### respiration
################################################################################

# Parametric coefficients
para_tbl <- broom::tidy(gamm_Reco_smoothed, parametric = TRUE) %>%
  dplyr::filter(p.value < 0.10 | term == "managementorganic") %>%
  # ^ keep significant + the management main effect even if n.s.
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

# Smooth terms  
smooth_tbl <- broom::tidy(gamm_Reco_smoothed, parametric = FALSE) %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

# Render with gt -- Reco
gt::gt(para_tbl) %>%
  gt::tab_header(title = "Selected parametric coefficients: Reco") %>%
  gt::fmt_number(decimals = 3)
gt::gt(smooth_tbl) %>% 
  gt::tab_header(title = "Selected smooth coefficients: Reco") %>% 
  gt::fmt_number(decimals = 3)
