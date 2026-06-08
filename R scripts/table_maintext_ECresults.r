library(broom)
library(gt)   # or kableExtra, flextable, etc.

# Parametric coefficients
para_tbl <- broom::tidy(bam_NEE, parametric = TRUE) %>%
  dplyr::filter(p.value < 0.10 | term == "managementorganic") %>%
  # ^ keep significant + the management main effect even if n.s.
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

# Smooth terms  
smooth_tbl <- broom::tidy(bam_NEE, parametric = FALSE) %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

# Render with gt
gt::gt(para_tbl) %>%
  gt::tab_header(title = "Selected parametric coefficients") %>%
  gt::fmt_number(decimals = 3)

para_tbl
