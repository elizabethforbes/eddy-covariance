create_management_vars_multicrop <- function(data, events) {
  
  result <- data
  
  # =========================================================================
  # STEP 1: Identify current crop for each date/field
  # =========================================================================
  
  # Get planting and harvest events with crop info
  planting_events <- events %>%
    filter(category == "planting") %>%
    select(management, plant_date = date, crop) %>%
    arrange(management, plant_date)
  
  harvest_events <- events %>%
    filter(category == "harvest") %>%
    select(management, harvest_date = date, crop) %>%
    arrange(management, harvest_date)
  
  # For each row, find current crop
  current_crop_info <- lapply(1:nrow(result), function(i) {
    row_field <- result$management[i]
    row_date <- result$date[i]
    
    # Most recent planting
    recent_plants <- planting_events %>%
      filter(management == row_field, plant_date <= row_date) %>%
      arrange(desc(plant_date))
    
    if(nrow(recent_plants) == 0) {
      return(list(
        crop = NA_character_,
        plant_date = as.Date(NA),
        harvest_date = as.Date(NA),
        in_season = 0,
        crop_type = NA_character_
      ))
    }
    
    # Most recent harvest
    recent_harvests <- harvest_events %>%
      filter(management == row_field, harvest_date <= row_date) %>%
      arrange(desc(harvest_date))
    
    # Check each recent planting to see if it's still active
    for(j in 1:nrow(recent_plants)) {
      this_plant <- recent_plants[j, ]
      
      # Find harvest for this specific crop planting
      matching_harvest <- recent_harvests %>%
        filter(crop == this_plant$crop, harvest_date > this_plant$plant_date) %>%
        slice(1)
      
      # If no harvest yet, or harvest is in future, this is current crop
      if(nrow(matching_harvest) == 0) {
        # Special case: alfalfa is perennial (no annual harvest end date)
        crop_type <- if_else(this_plant$crop == "alfalfa hay", "perennial", "annual")
        
        return(list(
          crop = this_plant$crop,
          plant_date = this_plant$plant_date,
          harvest_date = as.Date(NA),
          in_season = 1,
          crop_type = crop_type
        ))
      }
      
      if(matching_harvest$harvest_date > row_date) {
        crop_type <- if_else(this_plant$crop == "alfalfa", "perennial", "annual")
        
        return(list(
          crop = this_plant$crop,
          plant_date = this_plant$plant_date,
          harvest_date = matching_harvest$harvest_date,
          in_season = 1,
          crop_type = crop_type
        ))
      }
    }
    
    # Between crops (most recent was harvested)
    last_plant <- recent_plants[1, ]
    last_harvest <- recent_harvests %>%
      filter(crop == last_plant$crop) %>%
      slice(1)
    
    return(list(
      crop = NA_character_,
      plant_date = last_plant$plant_date,
      harvest_date = if(nrow(last_harvest) > 0) last_harvest$harvest_date else as.Date(NA),
      in_season = 0,
      crop_type = "fallow"
    ))
  })
  
  # Add current crop info to result
  result$current_crop <- sapply(current_crop_info, function(x) x$crop)
  result$current_plant_date <- as.Date(sapply(current_crop_info, function(x) x$plant_date), 
                                       origin = "1970-01-01")
  # result$current_harvest_date <- as.Date(sapply(current_crop_info, function(x) x$harvest_date), 
                                         # origin = "1970-01-01")
  result$in_season <- sapply(current_crop_info, function(x) x$in_season)
  result$crop_type <- sapply(current_crop_info, function(x) x$crop_type)
  
  # =========================================================================
  # STEP 2: Calculate other management events (tillage, fert, herbicide) that aren't necessarily crop-specific
  # =========================================================================
  
  for(event_type_name in c("tillage", 
                           "fertilizer", 
                           "herbicide" 
                           # "harvest"
                           )) {
    
    type_events <- events %>%
      filter(category == event_type_name) %>%
      select(management, event_date = date)
    
    if(nrow(type_events) == 0) {
      result[[paste0("last_", event_type_name)]] <- as.Date(NA)
      result[[paste0("days_since_", event_type_name)]] <- NA_real_
      next
    }
    
    last_event_dates <- sapply(1:nrow(result), function(i) {
      row_field <- result$management[i]
      row_date <- result$date[i]
      
      relevant <- type_events %>%
        filter(management == row_field, event_date <= row_date)
      
      if(nrow(relevant) == 0) {
        return(as.Date(NA))
      } else {
        return(max(relevant$event_date))
      }
    })
    
    last_event_dates <- as.Date(last_event_dates, origin = "1970-01-01")
    
    result[[paste0("last_", event_type_name)]] <- last_event_dates
    result[[paste0("days_since_", event_type_name)]] <- 
      as.numeric(result$date - last_event_dates)
  }
  
  # =========================================================================
  # STEP 3: Calculate days since planting (for current crop), 
  # days since harvest (for current crop)
  # =========================================================================
  
  result <- result %>%
    group_by(crop_type)
    mutate(
      days_since_planting = as.numeric(date - current_plant_date),
      # days_to_harvest = as.numeric(current_harvest_date - date),
      # days_since_any_harvest = as.numeric(date - current_plant_date),
      days_since_thiscrop_harvest = as.numeric(date - harvest_date)
    )
  
  # =========================================================================
  # STEP 4: CROP-SPECIFIC PHENOLOGY STAGES
  # =========================================================================
  
  result <- result %>%
    mutate(
      crop_stage = case_when(
        # No crop / fallow
        in_season == 0 | is.na(current_crop) ~ "fallow",
        
        # info from -- https://corn.agronomy.wisc.edu/AA/pdfs/A102.pdf (relative maturity = 105 days from planting)
        # other resources https://www.agry.purdue.edu/ext/corn/news/timeless/grainfill.html, https://nebraskacorn.gov/cornstalk/corn101/corn-growth-stages-explained/
        # ---- CORN (C4, ~110 days) ----
        current_crop == "corn" & days_since_planting < 0 ~ "freshly_planted",
        current_crop == "corn" & days_since_planting <= 10 ~ "emergence",
        current_crop == "corn" & days_since_planting <= 30 ~ "vegetative_early",  # V1-V5
        current_crop == "corn" & days_since_planting <= 50 ~ "vegetative_rapid",  # V6-V15
        current_crop == "corn" & days_since_planting <= 58 ~ "tasseling",         # VT
        current_crop == "corn" & days_since_planting <= 60 ~ "silking",           # R1-R2
        current_crop == "corn" & days_since_planting <= 90 ~ "grain_fill",       # R3-R5
        current_crop == "corn" & days_since_planting <= 110 ~ "maturity",         # R6
        current_crop == "corn" ~ "senescence",
        
        # ---- SOYBEAN (C3, ~120-150 days) ----
        # info and estimates on timing from https://agricbusiness.com.ng/soybean-growth-stages-days/
        current_crop == "soybean" & days_since_planting < 0 ~ "freshly_planted",
        current_crop == "soybean" & days_since_planting <= 7 ~ "emergence",      # VE
        current_crop == "soybean" & days_since_planting <= 45 ~ "vegetative",     # V1-V4
        current_crop == "soybean" & days_since_planting <= 55 ~ "flowering",      # R1-R2
        current_crop == "soybean" & days_since_planting <= 70 ~ "pod_development", # R3-R4
        current_crop == "soybean" & days_since_planting <= 100 ~ "seed_fill",     # R5-R6
        current_crop == "soybean" & days_since_planting <= 120 ~ "maturity",      # R7-R8
        current_crop == "soybean" ~ "senescence",
        
        # ---- BARLEY (C3, spring barley ~90-120 days) ----
        # https://cals.cornell.edu/field-crops/small-grains/growth-stages
        current_crop == "barley" & days_since_planting < 0 ~ "freshly_planted",
        current_crop == "barley" & days_since_planting <= 14 ~ "emergence",
        current_crop == "barley" & days_since_planting <= 35 ~ "tillering",
        current_crop == "barley" & days_since_planting <= 60 ~ "stem_elongation",
        current_crop == "barley" & days_since_planting <= 75 ~ "heading",
        current_crop == "barley" & days_since_planting <= 95 ~ "grain_fill",
        current_crop == "barley" & days_since_planting <= 115 ~ "maturity",
        current_crop == "barley" ~ "senescence",
        
        # ---- ALFALFA (perennial, multiple cuts per year) ----
        # For alfalfa, consider days since last harvest more relevant
        # So code by days since planting first, then days since harvest after the initial planting
        
        #  # No harvest yet or days_since_harvest is NA or 0: Use days_since_planting thresholds
        #     (is.na(days_since_harvest) | days_since_harvest == 0) & days_since_planting <= 14  ~ "Early Vegetative",
        #     (is.na(days_since_harvest) | days_since_harvest == 0) & days_since_planting <= 30  ~ "Vegetative Growth",
        #     (is.na(days_since_harvest) | days_since_harvest == 0) & days_since_planting <= 40  ~ "Bud Stage",
        #     (is.na(days_since_harvest) | days_since_harvest == 0) & days_since_planting <= 50  ~ "Early Bloom",
        #     (is.na(days_since_harvest) | days_since_harvest == 0) & days_since_planting <= 60  ~ "Full Bloom",
        #     (is.na(days_since_harvest) | days_since_harvest == 0) & days_since_planting > 60   ~ "Maturity or Beyond First Harvest",
        #     
        #     # After first harvest: use days_since_harvest thresholds
        #     days_since_harvest <= 7   ~ "Dormant / Recovery",
        #     days_since_harvest <= 14  ~ "Early Vegetative",
        #     days_since_harvest <= 21  ~ "Vegetative Growth",
        #     days_since_harvest <= 28  ~ "Bud Stage",
        #     days_since_harvest <= 35  ~ "Early Bloom",
        #     days_since_harvest > 35   ~ "Full Bloom or Maturity",
        #     
        #     # Default fallback
        #     TRUE ~ NA_character_
        #   )
        # }
        # 
        # 
        
        
        current_crop == "alfalfa hay" & days_since_planting <= 30 ~ "seedling",
        current_crop == "alfalfa hay" & is.na(last_harvest) ~ "vegetative",
        current_crop == "alfalfa hay" & !is.na(last_harvest) & 
          as.numeric(date - last_harvest) <= 7 ~ "post_cut_recovery",
        current_crop == "alfalfa hay" & !is.na(last_harvest) & 
          as.numeric(date - last_harvest) <= 28 ~ "vegetative",
        current_crop == "alfalfa hay" & !is.na(last_harvest) &
          as.numeric(date - last_harvest) <= 40 ~ "bud",
        current_crop == "alfalfa hay" & !is.na(last_harvest) & 
          as.numeric(date - last_harvest) <= 60 ~ "bloom",
        current_crop == "alfalfa hay" ~ "maturity",
        
        # ---- COVER CROPS (variable, ~60-120 days) ----
        current_crop == "cover_crop" & days_since_planting < 0 ~ "freshly_planted",
        current_crop == "cover_crop" & days_since_planting <= 14 ~ "emergence",
        current_crop == "cover_crop" & days_since_planting <= 45 ~ "establishment",
        current_crop == "cover_crop" & days_since_planting <= 90 ~ "vegetative",
        current_crop == "cover_crop" ~ "mature",
        
        # Default
        TRUE ~ "unknown"
      ),
      
      # Simplified stage grouping for modeling
      crop_stage_simple = case_when(
        crop_stage == "fallow" ~ "fallow",
        crop_stage %in% c("freshly_planted", "emergence", "establishment") ~ "early",
        crop_stage %in% c("vegetative", "vegetative_early", "vegetative_rapid", 
                          "tillering", "growth", "regrowth") ~ "vegetative",
        crop_stage %in% c("tasseling", "silking", "flowering", "heading", 
                          "stem_elongation", "pre_flower") ~ "reproductive",
        crop_stage %in% c("grain_fill", "seed_fill", "pod_development") ~ "grain_fill",
        crop_stage %in% c("maturity", "mature_growth", "mature", "senescence", 
                          "post_cut_recovery") ~ "mature",
        TRUE ~ "unknown"
      )
    )
  
  # =========================================================================
  # STEP 5: MANAGEMENT EFFECT VARIABLES (crop-specific)
  # =========================================================================
  
  
  # result <- result %>%
  #   mutate(
  #     # ---- Tillage effects ----
  #     tillage_pulse = case_when(
  #       is.na(days_since_tillage) ~ 0,
  #       days_since_tillage > 60 ~ 0,
  #       TRUE ~ exp(-days_since_tillage / 20)
  #     ),
  #     
  #     # ---- Fertilizer effects ----
  #     fert_effect = case_when(
  #       is.na(days_since_fertilizer) ~ 0,
  #       days_since_fertilizer > 60 ~ 0,
  #       days_since_fertilizer <= 45 ~ exp(-days_since_fertilizer / 15),
  #       TRUE ~ 0
  #     ),
  #     
  #     # ---- Herbicide effects ----
  #     herbicide_recent = case_when(
  #       is.na(days_since_herbicide) ~ 0,
  #       days_since_herbicide <= 14 ~ 1,
  #       TRUE ~ 0
  #     ),
  #     
  #     # ---- Post-harvest effects ----
  #     post_harvest = case_when(
  #       current_crop == "alfalfa" & !is.na(last_harvest) ~ 
  #         if_else(as.numeric(date - last_harvest) <= 14, 1, 0),
  #       
  #       in_season == 0 & !is.na(current_harvest_date) ~ 
  #         if_else(as.numeric(date - current_harvest_date) <= 30, 1, 0),
  #       
  #       TRUE ~ 0
  #     ),
  #     
  #     # ---- Crop functional type ----
  #     crop_functional_type = case_when(
  #       is.na(current_crop) ~ "none",
  #       current_crop == "corn" ~ "C4_annual",
  #       current_crop %in% c("soybean", "barley") ~ "C3_annual",
  #       current_crop == "alfalfa" ~ "C3_perennial",
  #       current_crop == "cover_crop" ~ "C3_annual",
  #       TRUE ~ "unknown"
  #     )
  #   )
  
  return(result)
}
