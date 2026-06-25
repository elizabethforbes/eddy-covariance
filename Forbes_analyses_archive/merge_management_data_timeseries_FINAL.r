create_management_vars_multicrop <- function(data, events) {
  
  # Ensure the column name matches between data and events
  # The function expects 'management' but data may have 'field'
  if("field" %in% names(events) && !"management" %in% names(events)) {
    events <- events %>% rename(management = field)
  }
  if("field" %in% names(data) && !"management" %in% names(data)) {
    data <- data %>% rename(management = field)
  }
  
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
      
      # SPECIAL HANDLING FOR PERENNIAL CROPS
      # Perennial crops (like alfalfa) are never "ended" by harvest events
      # They continue until explicitly terminated or replaced by another planting
      if(this_plant$crop == "alfalfa hay") {
        return(list(
          crop = this_plant$crop,
          plant_date = this_plant$plant_date,
          harvest_date = as.Date(NA),  # Perennials don't have a final harvest date
          in_season = 1,
          crop_type = "perennial"
        ))
      }
      
      # ANNUAL CROPS - check if harvest has ended the season
      # Find harvest for this specific crop planting
      matching_harvest <- recent_harvests %>%
        filter(crop == this_plant$crop, harvest_date > this_plant$plant_date) %>%
        slice(1)
      
      # If no harvest yet, or harvest is in future, this is current crop
      if(nrow(matching_harvest) == 0) {
        crop_type <- "annual"
        
        return(list(
          crop = this_plant$crop,
          plant_date = this_plant$plant_date,
          harvest_date = as.Date(NA),
          in_season = 1,
          crop_type = crop_type
        ))
      }
      
      if(matching_harvest$harvest_date > row_date) {
        crop_type <- "annual"
        
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
  result$current_harvest_date <- as.Date(sapply(current_crop_info, function(x) x$harvest_date), 
                                         origin = "1970-01-01")
  result$in_season <- sapply(current_crop_info, function(x) x$in_season)
  result$crop_type <- sapply(current_crop_info, function(x) x$crop_type)
  
  # =========================================================================
  # STEP 2: Calculate other management events (tillage, fert, herbicide) that aren't necessarily crop-specific
  # =========================================================================
  
  # for(event_type_name in c("tillage", 
  #                          "fertilizer", 
  #                          "herbicide" 
  #                          # "harvest"
  #                          )) {
  #   
  #   type_events <- events %>%
  #     filter(category == event_type_name) %>%
  #     select(management, event_date = date)
  #   
  #   if(nrow(type_events) == 0) {
  #     result[[paste0("last_", event_type_name)]] <- as.Date(NA)
  #     result[[paste0("days_since_", event_type_name)]] <- NA_real_
  #     next
  #   }
  #   
  #   last_event_dates <- sapply(1:nrow(result), function(i) {
  #     row_field <- result$management[i]
  #     row_date <- result$date[i]
  #     
  #     relevant <- type_events %>%
  #       filter(management == row_field, event_date <= row_date)
  #     
  #     if(nrow(relevant) == 0) {
  #       return(as.Date(NA))
  #     } else {
  #       return(max(relevant$event_date))
  #     }
  #   })
  #   
  #   last_event_dates <- as.Date(last_event_dates, origin = "1970-01-01")
  #   
  #   result[[paste0("last_", event_type_name)]] <- last_event_dates
  #   result[[paste0("days_since_", event_type_name)]] <- 
  #     as.numeric(result$date - last_event_dates)
  # }
  
  for(event_type_name in c("tillage", "fertilizer", "herbicide")) {
    
    type_events <- events %>%
      filter(category == event_type_name) %>%
      select(management, event_date = date)
    
    # Create event indicator column name
    event_col_name <- paste0(event_type_name, "_event")
    
    # Initialize the event indicator column with NA
    result[[event_col_name]] <- NA_integer_
    
    if(nrow(type_events) == 0) {
      result[[paste0("last_", event_type_name)]] <- as.Date(NA)
      result[[paste0("days_since_", event_type_name)]] <- NA_real_
      next
    }
    
    # Identify rows in result where an event occurred on that date and management
    event_occurrences <- sapply(1:nrow(result), function(i) {
      any(type_events$management == result$management[i] & 
            type_events$event_date == result$date[i])
    })
    
    # Assign 1 where event occurred, else NA stays
    result[[event_col_name]][event_occurrences] <- 1L
    
    # Existing code to find last event dates and days since
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
    result[[paste0("days_since_", event_type_name)]] <- as.numeric(result$date - last_event_dates)
  }
  
  
  # =========================================================================
  # STEP 2.5: Calculate last harvest date for all harvest types
  # This is needed for alfalfa phenology stages
  # =========================================================================
  
  # Get ALL harvest events (not just crop-specific ones)
  all_harvest_events <- events %>%
    filter(category %in% c("harvest", "haying")) %>%
    select(management, event_date = date)
  
  if(nrow(all_harvest_events) > 0) {
    last_harvest_dates <- sapply(1:nrow(result), function(i) {
      row_field <- result$management[i]
      row_date <- result$date[i]
      
      relevant <- all_harvest_events %>%
        filter(management == row_field, event_date <= row_date)
      
      if(nrow(relevant) == 0) {
        return(as.Date(NA))
      } else {
        return(max(relevant$event_date))
      }
    })
    
    result$last_harvest <- as.Date(last_harvest_dates, origin = "1970-01-01")
  } else {
    result$last_harvest <- as.Date(NA)
  }
  
  # =========================================================================
  # STEP 2.6: Calculate last harvest date for the SAME crop (important for perennials)
  # =========================================================================
  
  # Get all harvest events with crop info and management
  all_harvest_events_crop <- events %>%
    filter(category %in% c("harvest", "haying")) %>%
    select(management, crop, event_date = date)
  
  # Calculate last harvest date of the current crop for each row
  last_harvest_same_crop_dates <- sapply(1:nrow(result), function(i) {
    row_field <- result$management[i]
    row_crop <- result$current_crop[i]
    row_date <- result$date[i]
    
    relevant <- all_harvest_events_crop %>%
      filter(management == row_field,
             crop == row_crop,
             event_date <= row_date)
    
    if(nrow(relevant) == 0) {
      return(as.Date(NA))
    } else {
      return(max(relevant$event_date))
    }
  })
  
  result$last_harvest_same_crop <- as.Date(last_harvest_same_crop_dates, origin = "1970-01-01")
  
  # Calculate days since last harvest for alfalfa only
  result <- result %>%
    mutate(
      days_since_last_harvest_same_crop = ifelse(
        current_crop == "alfalfa hay",
        as.numeric(date - last_harvest_same_crop),
        NA_real_
      )
    )
  
  
  
  # =========================================================================
  # STEP 3: Calculate days since planting (for current crop), 
  # days since harvest (for current crop)
  # =========================================================================
  
  result <- result %>%
    mutate(
      days_since_planting = as.numeric(date - current_plant_date),
      days_to_harvest = as.numeric(current_harvest_date - date),
      # For annual crops: days since the harvest that ended the season
      # For perennial crops: this will be NA (since current_harvest_date is NA)
      days_since_season_end_harvest = as.numeric(date - current_harvest_date),
      # For all crops: days since the most recent harvest event of any kind
      days_since_last_harvest = as.numeric(date - last_harvest)
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
        
        # # ---- ALFALFA (perennial, multiple cuts per year) ----
        # # For alfalfa, consider days since last harvest more relevant, but post-planting, it's the same as for annual crops:
        # # First check if it's early establishment (within 30 days of planting)
        # # https://greg.app/alfalfa-lifecycle/, https://library.ndsu.edu/server/api/core/bitstreams/4de5d19c-5d02-4522-b6e4-0291bbcc6d61/content, 
        # # https://alfalfa.ucdavis.edu/sites/g/files/dgvnsk12586/files/media/documents/UCAlfalfa8289GrowthDev-reg.pdf, https://extension.psu.edu/harvest-management-of-alfalfa,
        # # https://greenlifelawnandlandscape.com/alfalfa-plant-stages-of-growth/
        # 
        # # AFTER FIRST HARVEST: https://eos.com/blog/growing-alfalfa/, https://legacy.research.agrilife.org/wp-content/uploads/sites/3/2011/10/nmsugrowthstages_11.pdf,
        # # https://extension.psu.edu/harvest-management-of-alfalfa, https://bookstore.ksre.ksu.edu/pubs/alfalfa-growth-and-development-poster-20-x-30_MF3348.pdf
        # # After harvest, determine stage based on days since harvest, since regrowth is different than initial emergence and growth:

        # # Alfalfa, Pre-harvest stages (no harvest yet)
        # current_crop == "alfalfa hay" & is.na(days_since_last_harvest_same_crop) & days_since_planting >= 0 & days_since_planting <= 10 ~ "germination",
        # current_crop == "alfalfa hay" & is.na(days_since_last_harvest_same_crop) & days_since_planting > 10 & days_since_planting <= 30 ~ "seedling",
        # current_crop == "alfalfa hay" & is.na(days_since_last_harvest_same_crop) & days_since_planting > 30 & days_since_planting <= 72 ~ "vegetative",
        # current_crop == "alfalfa hay" & is.na(days_since_last_harvest_same_crop) & days_since_planting > 72 & days_since_planting <= 93 ~ "flowering",
        # current_crop == "alfalfa hay" & is.na(days_since_last_harvest_same_crop) & days_since_planting > 93 & days_since_planting <= 128 ~ "seed production",
        # current_crop == "alfalfa hay" & is.na(days_since_last_harvest_same_crop) & days_since_planting > 128 ~ "maturity",
        # 
        # # Alfalfa, Post-harvest stages (harvests have occurred)
        # current_crop == "alfalfa hay" & !is.na(days_since_last_harvest_same_crop) & days_since_last_harvest_same_crop >= 0 & days_since_last_harvest_same_crop <= 5 ~ "post_cut_recovery",
        # current_crop == "alfalfa hay" & !is.na(days_since_last_harvest_same_crop) & days_since_last_harvest_same_crop > 5 & days_since_last_harvest_same_crop <= 10 ~ "seedling_regrowth",
        # current_crop == "alfalfa hay" & !is.na(days_since_last_harvest_same_crop) & days_since_last_harvest_same_crop > 10 & days_since_last_harvest_same_crop <= 30 ~ "vegetative",
        # current_crop == "alfalfa hay" & !is.na(days_since_last_harvest_same_crop) & days_since_last_harvest_same_crop > 30 & days_since_last_harvest_same_crop <= 35 ~ "bud",
        # current_crop == "alfalfa hay" & !is.na(days_since_last_harvest_same_crop) & days_since_last_harvest_same_crop > 35 & days_since_last_harvest_same_crop <= 65 ~ "flowering",
        # current_crop == "alfalfa hay" & !is.na(days_since_last_harvest_same_crop) & days_since_last_harvest_same_crop > 65 ~ "maturity",
        # 
        # # ---- COVER CROPS (variable, ~60-120 days) ----
        # current_crop == "cover crop" & days_since_planting < 0 ~ "freshly_planted",
        # current_crop == "cover crop" & days_since_planting <= 14 ~ "emergence",
        # current_crop == "cover crop" & days_since_planting <= 45 ~ "establishment",
        # current_crop == "cover crop" & days_since_planting <= 90 ~ "vegetative",
        # current_crop == "cover crop" ~ "mature",
        
        # update: assign cover crops and alfalfa as NAs and use separate function to categorize them because they overwinter
        current_crop == "cover crop" ~ NA,
        current_crop == "alfalfa hay" ~ NA,
        
        # Default
        TRUE ~ "unknown"
      ),
      
      # Simplified stage grouping for modeling
      crop_stage_simple = case_when(
        crop_stage %in% c("fallow", "post_cut_recovery") ~ # categorizing post cut as fallow for phenology/gas exchange purposes
          "fallow",
        
        crop_stage %in% c("freshly_planted", "emergence", "establishment",
                          "germination", "seedling", "seedling_regrowth") ~ 
          "early",
        
        crop_stage %in% c("vegetative", "vegetative_early", "vegetative_rapid", 
                          "tillering", "growth", "regrowth") ~ 
          "vegetative",
        
        crop_stage %in% c("tasseling", "silking", "flowering", "heading", 
                          "stem_elongation", "pre_flower", "bud") ~ 
          "reproductive",
        
        crop_stage %in% c("grain_fill", "seed_fill", "pod_development", "bloom", "seed production") ~ 
          "grain_fill",
        
        crop_stage %in% c("maturity", "mature_growth", "mature", "senescence") ~ 
          "mature", # putting senescent stage here because it represents a slowdown, not a cessation, of metabolism like dormant does
        
        crop_stage %in% NA ~ NA,
        
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
