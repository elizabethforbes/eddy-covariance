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

# summarize:
som_summ <- som %>% 
  group_by(management, field_replicate, sample_location, sample_depth_cm) %>% 
  summarize(across(where(is.numeric), mean, na.rm = TRUE)) %>% 
  select(!c(sample_replicate)) %>% 
  mutate(core_id = paste(julian_date, management, field_replicate, sample_location),
         management = factor(management, levels = c("conventional", "organic")))

# summarize to 1m using bulk density:
som_summ <- som_summ %>% 
  left_join(bd, by = join_by(management, 
                             field_replicate == replicate,
                             sample_location == location,
                             sample_depth_cm == depth_cm)) 
som_summ <- som_summ %>% 
  mutate(mass_somsample = S_BD_gcm3.x * volume_layer.x, # volume in cm3 for each 3" by 10cm segment, bulk density of sieved soil
         weighted_som = (OM_perc/100) * mass_somsample) # get SOM in mass by multiplying %som against mass of sample in that layer

weighted_avg_som <- som_summ %>% 
  # filter out depths > 80cm
  filter(sample_depth_cm < 80) %>% 
  group_by(management, field_replicate, sample_location) %>% 
  summarise(
    total_weighted_som = sum(weighted_som, na.rm = TRUE),
    total_mass = sum(mass_somsample, na.rm = TRUE),
    somperc_avg_80cm = (total_weighted_som/total_mass)*100
  ) %>% 
  select(total_weighted_som, total_mass, somperc_avg_80cm, management)

# initial plot:
p_som <- som_summ %>% 
  ggplot(aes(x = factor(sample_depth_cm), y = OM_perc, fill = management)) +
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

# model differences by treatment, 80cm OM
som_mod <- lm(somperc_avg_80cm ~ management, data = weighted_avg_som)
anova(som_mod)
# Response: somperc_avg_80cm
#             Df  Sum Sq Mean Sq F value Pr(>F)  
# management  1 0.43064 0.43064  4.6119 0.0573 .
# Residuals  10 0.93376 0.09338                        

summary(som_mod)# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.35268 -0.15632 -0.06239  0.08847  0.59234 
# 
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)         1.4612     0.1248  11.713 3.67e-07 ***
#   managementorganic   0.3789     0.1764   2.148   0.0573 .  
# 
# Residual standard error: 0.3056 on 10 degrees of freedom
# Multiple R-squared:  0.3156,	Adjusted R-squared:  0.2472 
# F-statistic: 4.612 on 1 and 10 DF,  p-value: 0.0573

DHARMa::plotQQunif(som_mod)
DHARMa::plotResiduals(som_mod)

# plot of core-level OM concentration:
weighted_avg_som %>%
  ggplot(aes(x = management, y = somperc_avg_80cm, fill = management)) +
  # geom_violin(alpha = 0.5, color = NA) +
  geom_boxplot(color = "black", alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = management),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 4, alpha = 0.8) +
  scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  labs(x = NULL, y = "Organic matter content to 80cm (%)") +
  my_custom_theme()+
  theme(legend.position = "none")

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
  filter(section_depth_cm != 75) %>% 
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
  filter(section_depth_cm != 75) %>% 
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

# Total stocks boxplot by management
p_cstock <- 
  ggplot(stocks_summ, aes(x = management, y = Cstock_kgm2_70cm, fill = management)) +
  geom_boxplot(color = "black", alpha = 0.5) +
  geom_jitter(aes(color = management))+
  scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  labs(y = "Total Carbon Stock (kg/m² to 0.7m)",
       x = element_blank()) +
  my_custom_theme()+
  theme(legend.position = "none")

p_nstock <- 
  ggplot(stocks_summ, aes(x = management, y = Nstock_kgm2_70cm, fill = management)) +
  geom_boxplot(color = "black", alpha = 0.5) +
  geom_jitter(aes(color = management))+
  scale_fill_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  scale_color_manual(values = c("conventional" = "tomato", "organic" = "#E69F00")) +
  labs(y = "Total Nitrogen Stock (kg/m² to 0.7m)",
       x = element_blank()) +
  my_custom_theme()+
  theme(legend.position = "none")

p_cstock
p_nstock

p_cstock + p_nstock 

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
