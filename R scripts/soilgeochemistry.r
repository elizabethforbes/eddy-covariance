# ============================================================================
# MODELING and VISUALIZATION of SOIL DATA: EC study, 2020 core harvest
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)
library(scales)
library(readxl)
library(here)
library(DHARMa)
library(lmtest)
library(nlme)
library(mgcv)
library(lubridate)
library(patchwork)

my_custom_theme <- function() {
  theme_minimal(base_size = 12, base_family = "Helvetica") %+replace%
    theme(
      panel.grid.major = element_line(color = "gray85", size = 0.15),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(color = "gray20"),
      axis.ticks = element_line(color = "gray40"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.background = element_blank(),
      strip.background = element_rect(fill = "gray90", color = NA),
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
      plot.background = element_blank()
    )
}

# ============================================================================
# SOM: plot, model
# ============================================================================

# soil type: Knickerbocker fine sandy loam
som <- readxl::read_xlsx(here("eddy_covariance_fluxdata", 
                              "EF0_LOI (2020, EF01 & EF02).xlsx"), sheet = 2)

# rename "sample_depth_cm" to "section_depth_cm" for table merging later
som <- som %>% 
  rename(section_depth_cm = sample_depth_cm)

# summarize: average the 4 LOI replicate subsamples per depth, per core
som_dat <- som %>% 
  filter(section_depth_cm <= 70) %>% 
  group_by(management, field_replicate, sample_location, section_depth_cm) %>% 
  summarize(
    OM_perc = mean(OM_perc, na.rm = TRUE),
    .groups = "drop") %>% 
  mutate(
    management = factor(management, levels = c("conventional", "organic")),
    core_id = paste(management, field_replicate, sample_location, sep = "_"))

# initial plot:
p_som <- som_dat %>% 
  ggplot(aes(x = factor(section_depth_cm), y = OM_perc, fill = management)) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", alpha = 0.5) +
  geom_jitter(aes(color = management), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) + 
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  labs(x = "Sample Depth (cm)", y = "OM (%)") +
  my_custom_theme()

p_som
# theme_minimal()

# model differences by treatment, 10cm increments to 70cm
# GAMM with management-specific smooths
# random effect for core (replicate/location) to account for nesting
som_gamm <- gamm(OM_perc ~ management +
                   s(section_depth_cm, by = management, k = 5),
                 random = list(core_id = ~1),
                 data = som_dat)
summary(som_gamm$gam)
# Parametric coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)         1.6273     0.1224  13.296   <2e-16 ***
#   managementorganic   0.3483     0.1731   2.012   0.0476 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#                                               edf Ref.df     F p-value    
#   s(section_depth_cm):managementconventional 2.003  2.003 153.6  <2e-16 ***
#   s(section_depth_cm):managementorganic      1.000  1.000 197.8  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.751   
# Scale est. = 0.071954  n = 84

# plot the smooth terms:
gratia::draw(som_gamm$gam) & theme_classic() 

# ============================================================================
# Bulk density: 2mm (fine fraction)
# ============================================================================

# soil type: Knickerbocker fine sandy loam
# bd <- readxl::read_xlsx(here("eddy_covariance_fluxdata", "Bulk Density Fall 2020.xlsx"), sheet = 3)
sn <- readxl::read_xlsx(here("eddy_covariance_fluxdata", "2015-2020_Nitrogen_Results.xlsx"), sheet = 9)

# summarize:
bd_summ <- sn %>% 
  filter(section_depth_cm <= 70) %>% # filter to the deepest consistently-sampled depth: 70cm
  group_by(management, replicate, location) %>% 
  summarise(
    finemass_g_70cm = sum(bd_under2mm_gcm3 * section_vol_cm3, na.rm = TRUE), # calculate mass of each section's fine fraction
    volume_core = sum(section_vol_cm3, na.rm = TRUE), # calculate the volume of the core to 70cm; should be identical
    max_depth = n() * 10,
    .groups = "drop_last") %>% 
  mutate(
    bd_70cmcore_fine_gcm3 = finemass_g_70cm / volume_core,
    management = factor(management, levels = c("conventional", "organic"))
  )

# initial plot:
bd_summ %>%
  # indicate whether a core is less than 1m in depth:
  # mutate(depth_category = if_else(max_depth_cm < 100, "<1m", "=1m")) %>%
  ggplot(aes(x = management, y = bd_70cmcore_fine_gcm3, fill = management)) +
  geom_violin(alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  # indicate depth disparity with shape:
  geom_jitter(aes(color = management),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 4, alpha = 0.8) +
  # scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  # scale_shape_manual(values = c("<1m" = 17, "=1m" = 16)) +  # triangles for <1m, circles for ≥1m
  labs(x = NULL, y = "Bulk Density to 70cm (g/cm³)") +
  my_custom_theme()+
  theme(legend.position = "none")

# plot with 10cm increments:
sn %>%
  # filter out the 75cm depth section:
  filter(section_depth_cm != 75) %>% 
  ggplot(aes(x = factor(section_depth_cm), y = bd_under2mm_gcm3, fill = management)) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", alpha = 0.5) +
  geom_jitter(aes(color = management), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  labs(x = "Sample Depth (cm)", y = "bulk density (g/cm³)", ) +
  my_custom_theme()

# ============================================================================
# Carbon, nitrogen concentrations:
# ============================================================================

# percent nitrogen, descending by depth:
p_n <- sn %>%
  # filter out the 75cm depth section:
  filter(section_depth_cm != 75) %>% 
  ggplot(aes(x = factor(section_depth_cm), y = perc_N, fill = management)) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", alpha = 0.5) +
  geom_jitter(aes(color = management), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  labs(x = "Sample Depth (cm)", y = "percent nitrogen content", ) +
  my_custom_theme()

# percent carbon, descending by depth:
p_c <- sn %>%
  # filter out the 75cm depth section:
  filter(section_depth_cm != 75) %>% 
  ggplot(aes(x = factor(section_depth_cm), y = perc_C, fill = management)) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", alpha = 0.5) +
  geom_jitter(aes(color = management), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  labs(x = "Sample Depth (cm)", y = "percent carbon content", ) +
  my_custom_theme()

p_n
p_c

# ============================================================================
# calculate total carbon, nitrogen stocks to max depth:
# ============================================================================

# summarize to 70cm depth (the shallowest depth across all cores)
stocks_summ <- sn %>% 
  filter(section_depth_cm <= 70) %>%
  group_by(management, replicate, location) %>% 
  summarize(
    Cstock_kgm2_70cm = sum(stockC_kgm2, na.rm = TRUE),
    Nstock_kgm2_70cm = sum(stockN_kgm2, na.rm = TRUE),
    .groups = "drop_last")

# Profile plot: stocks by depth and management
p_cdepths <- sn %>% 
  # filter out the 75cm depth section:
  filter(section_depth_cm <= 70) %>% 
  ggplot(aes(y = stockC_kgm2, x = as.factor(section_depth_cm), fill = management)) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", alpha = 0.5) +
  geom_jitter(aes(color = management), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  # scale_y_reverse() +
  labs(y = "Carbon stock (kg/m² per 10cm)", x = "Depth (cm)") +
  theme(legend.position = "bottom") +
  my_custom_theme() + theme(legend.position = "none")

p_ndepths <- sn %>% 
  # filter out the 75cm depth section:
  filter(section_depth_cm <= 70) %>% 
  ggplot(aes(y = stockN_kgm2, x = as.factor(section_depth_cm), fill = management)) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", alpha = 0.5) +
  geom_jitter(aes(color = management), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  # scale_y_reverse() +
  labs(y = "Nitrogen stock (kg/m² per 10cm)", x = "Depth (cm)") +
  theme(legend.position = "bottom") +
  my_custom_theme()

p_cdepths + p_ndepths

# ============================================================================
# model: comparison of C, N stocks by management
# ============================================================================

# model differences in total stock (to 80cm) by treatment:
cs_mod <- lm(Cstock_kgm2_70cm ~ management, data = stocks_summ)
anova(cs_mod)
# Analysis of Variance Table
# 
# Response: Cstock_kgm2_70cm
# Df  Sum Sq Mean Sq F value Pr(>F)
# management  1  4.4133  4.4133  3.2724 0.1006
# Residuals  10 13.4862  1.3486    

DHARMa::plotQQunif(cs_mod)
DHARMa::plotResiduals(cs_mod)

summary(cs_mod)
# Min      1Q  Median      3Q     Max 
# -1.6169 -0.7132 -0.1032  0.4521  2.2214 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         5.0673     0.4741  10.688 8.61e-07 ***
#   managementorganic   1.2129     0.6705   1.809    0.101    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.161 on 10 degrees of freedom
# Multiple R-squared:  0.2466,	Adjusted R-squared:  0.1712 
# F-statistic: 3.272 on 1 and 10 DF,  p-value: 0.1006

ns_mod <- lm(Nstock_kgm2_70cm ~ management, data = stocks_summ)
anova(ns_mod)
# Analysis of Variance Table
# 
# Response: Nstock_kgm2_70cm
# Df   Sum Sq  Mean Sq F value  Pr(>F)  
# management  1 0.035967 0.035967  6.0144 0.03412 *
#   Residuals  10 0.059801 0.005980                       

plotQQunif(ns_mod)
plotResiduals(ns_mod)

summary(ns_mod)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.08911 -0.05076 -0.00853  0.03713  0.15568 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)        0.67518    0.03157  21.387 1.11e-09 ***
#   managementorganic  0.10949    0.04465   2.452   0.0341 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07733 on 10 degrees of freedom
# Multiple R-squared:  0.3756,	Adjusted R-squared:  0.3131 
# F-statistic: 6.014 on 1 and 10 DF,  p-value: 0.03412

# ============================================================================
# model: comparison of C, N stocks by management AND by sample depth
# ============================================================================

sn_mod <- sn %>%
  filter(section_depth_cm <= 70) %>%
  mutate(
    management = factor(management, levels = c("conventional", "organic")),
    core_id = paste(management, replicate, location, sep = "_")
  )

# verify
class(sn_mod$management)  # should be "factor"
levels(sn_mod$management) # should show both levels

# carbon stocks by depth:
c_gamm <- gamm(stockC_kgm2 ~ management +
                 s(section_depth_cm, by = management, k = 5),
               random = list(core_id = ~1),
               data = sn_mod)

summary(c_gamm$gam)
#   Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)        0.72390    0.06258  11.568   <2e-16 ***
#   managementorganic  0.17327    0.08850   1.958   0.0539 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#                                               edf   Ref.df    F   p-value    
#   s(section_depth_cm):managementconventional 2.681  2.681  95.54  <2e-16 ***
#   s(section_depth_cm):managementorganic      3.174  3.174 104.53  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.818   
# Scale est. = 0.032831  n = 84

gratia::draw(c_gamm$gam) & theme_bw()
gam.check(c_gamm$gam)

# nitrogen stocks by depth:
n_gamm <- gamm(stockN_kgm2 ~ management +
                 s(section_depth_cm, by = management, k = 5),
               random = list(core_id = ~1),
               data = sn_mod)

summary(n_gamm$gam)
#   Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)       0.096455   0.004167  23.147  < 2e-16 ***
#   managementorganic 0.015642   0.005893   2.654  0.00966 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#                                               edf   Ref.df     F p-value    
#   s(section_depth_cm):managementconventional 1.739  1.739 91.53  <2e-16 ***
#   s(section_depth_cm):managementorganic      3.597  3.597 94.95  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.826   
# Scale est. = 0.00023259  n = 84

gratia::draw(n_gamm$gam) & theme_bw()
gam.check(n_gamm$gam)

# ============================================================================
# TABLES
# ============================================================================

library(gtsummary)
library(flextable)
# set theme to compact: reduces row padding, font size
set_gtsummary_theme(theme_gtsummary_compact(set_theme = TRUE))

# table of results for SOM analysis:
SOM_tab <- 
  tbl_regression(
  som_gamm$gam,
  exponentiate = FALSE,
  conf.int = TRUE) %>% 
  add_significance_stars(
    hide_ci = TRUE, hide_se = FALSE,
    hide_p = FALSE,
    pattern = "{p.value}{stars}"
  ) %>% 
  italicize_levels() %>% bold_labels() %>% 
  modify_header(estimate~"**Estimate**") %>% 
  # update level names for the interactions
  modify_table_body(~ .x %>% 
                      dplyr::mutate(
                        label = dplyr::recode(label, "s(section_depth_cm):managementconventional" = "sample depth * conventional",
                                              "s(section_depth_cm):managementorganic" = "sample depth * organic")
                      ))
# remove footnotes:
SOM_tab$table_styling$abbreviation <- 
  SOM_tab$table_styling$abbreviation %>% 
  dplyr::filter(column != c("conf.low", "std.error"))
SOM_tab %>% 
  modify_caption("**SOM concentration by management and depth (10cm increments to 70cm)**") %>% 
  as_gt()

# make side-by-side table for C and N stocks
# helper function:
make_tbl <- function(model) {
  tbl_regression(
    model,
    exponentiate = FALSE,
    conf.int     = TRUE
  ) |>
    add_significance_stars(
      hide_ci  = TRUE,
      hide_se  = FALSE,
      hide_p = FALSE,
      pattern  = "{p.value}{stars}"
    ) |>
    italicize_levels() %>% 
    remove_row_type(type = "reference") |>
    modify_header(estimate ~ "**Estimate**", std.error ~ "**SE**")
}

Cstock_tab <- make_tbl(c_gamm$gam)
Nstock_tab <- make_tbl(n_gamm$gam)


# update level names for the interactions
Cstock_tab <- 
  Cstock_tab %>% 
  modify_table_body(~ .x %>% 
      dplyr::mutate(
        label = dplyr::recode(label, "s(section_depth_cm):managementconventional" = "sample depth * conventional",
                          "s(section_depth_cm):managementorganic" = "sample depth * organic")))
Nstock_tab <- 
  Nstock_tab %>% 
  modify_table_body(~ .x %>% 
                      dplyr::mutate(
                        label = dplyr::recode(label, "s(section_depth_cm):managementconventional" = "sample depth * conventional",
                                              "s(section_depth_cm):managementorganic" = "sample depth * organic")))
# Remove footnotes
Cstock_tab$table_styling$abbreviation <- 
  Cstock_tab$table_styling$abbreviation %>% 
  dplyr::filter(column != c("conf.low", "std.error"))
Nstock_tab$table_styling$abbreviation <- 
  Nstock_tab$table_styling$abbreviation %>% 
  dplyr::filter(column != c("conf.low", "std.error"))

# combine:
tbl_merge(
  tbls = list(Cstock_tab, Nstock_tab),
  tab_spanner = c(
    "**Carbon stocks**", "**Nitrogen stocks**")
) %>% 
  bold_labels() %>% 
  modify_caption("**Predictors of C, N stocks by field management and sample depth**") %>% 
  as_gt()

# combine with SOM table, too:
tbl_merge(
  tbls = list(SOM_tab, Cstock_tab, Nstock_tab),
  tab_spanner = c("**SOM content**", "**C stocks**", "**N stocks**")
  ) %>% 
  bold_labels() %>% 
  modify_caption("**Predictors of SOM content, carbon and nitrogen stocks by management and depth**") %>% 
  as_gt()
