# =============================================================================
# Generating cumulative calculations of NEE, GPP, Reco using gap-filled data
# =============================================================================

library(dplyr)

# Standard Errors (SE) for each NEE value propagate as the square root of 
# summed variances (assuming independent errors).
nee_yearly_cumu <- ec_daily_mngmnt %>% 
  select(management, date, year,
         NEE, NEE_gapfilled, NEE_gf_se,
         air_temperature, VPD, current_crop,
         days_since_tillage, days_since_fertilizer, days_since_herbicide,
         crop_stage, crop_stage_simple) %>% 
  arrange(date) %>% 
  group_by(management, year) %>% 
  mutate(nee_cumu = cumsum(NEE_gapfilled),
         variance_cumu = cumsum(NEE_gf_se^2),          # Sum of variances (SE^2)
         nee_se_cumu = sqrt(variance_cumu)    # Propagated SE as sqrt of variance
         ) %>% 
  ungroup()

# repeat for GPP, Reco
gpp_yearly_cumu <- ec_daily_mngmnt %>% 
  select(management, date, year,
         GPP, GPP_gapfilled,             # swap out variables here
         air_temperature, VPD, current_crop,
         days_since_tillage, days_since_fertilizer, days_since_herbicide,
         crop_stage, crop_stage_simple) %>% 
  arrange(date) %>% 
  group_by(management, year) %>% 
  mutate(gpp_cumu = cumsum(GPP_gapfilled)
         # variance_cumu = cumsum(GPP_gf_se^2),          # Sum of variances (SE^2)
         # gpp_se_cumu = sqrt(variance_cumu)    # Propagated SE as sqrt of variance
  ) %>% 
  ungroup()

reco_yearly_cumu <- ec_daily_mngmnt %>% 
  select(management, date, year,
         Reco, Reco_gapfilled, Reco_gf_se,             # swap out variables here
         air_temperature, VPD, current_crop,
         days_since_tillage, days_since_fertilizer, days_since_herbicide,
         crop_stage, crop_stage_simple) %>% 
  arrange(date) %>% 
  group_by(management, year) %>% 
  mutate(reco_cumu = cumsum(Reco_gapfilled),
         variance_cumu = cumsum(Reco_gf_se^2),          # Sum of variances (SE^2)
         reco_se_cumu = sqrt(variance_cumu)    # Propagated SE as sqrt of variance
  ) %>% 
  ungroup()
