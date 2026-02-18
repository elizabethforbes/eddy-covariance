library(tidyverse)
library(readxl)
library(purrr)
library(lubridate)
library(hms)
library(here)

# ============================================================================
# 1. upload and summarize tower-based PAR data: 2020
# ============================================================================

# organic, 2020
parorg_2020 <- read_csv(here("2020 PAR etc. data", "2020 EF01 Par Data.csv"))
parorg_2020 <- parorg_2020 %>% 
  select(!`...6`) %>% 
  rename(datetime = `Date Time, GMT-04:00`,
         airt_F = `Temp, °F (LGR S/N: 20205284, SEN S/N: 20168088, LBL: Air Temperature)`,
         rh = `RH, % (LGR S/N: 20205284, SEN S/N: 20168088)`,
         PAR = `PAR, µmol/m²/s (LGR S/N: 20205284, SEN S/N: 20189138)`,
         soilt_F = `Temp, °F (LGR S/N: 20205284, SEN S/N: 20245643, LBL: Soil Temperature)`) %>% 
  # convert datettime from character to datetime
  mutate(datetime = as_datetime(datetime, 
                                format = "%m/%d/%Y %H:%M",
                                tz = NULL)) %>% 
  mutate(date = date(datetime),
         hour = hour(datetime)) %>% 
  # convert F to C for both temps:
  mutate(soilt_C = (soilt_F - 32) * 5/9,
         airt_C = (airt_F - 32)*5/9) %>% 
  select(!c(soilt_F, airt_F)) %>% 
  # average by the hour
  group_by(date, hour) %>% 
  summarise(across(is.numeric, mean, na.rm = TRUE)) %>% 
  ungroup() %>% 
  # add a "management" column (here, organic)
  mutate(management = "organic")

# conventional, 2020
parconv_2020 <- read_csv(here("2020 PAR etc. data", "2020 EF02 Par Data.csv"))
parconv_2020 <- parconv_2020 %>% 
  # select(!`...6`) %>% 
  rename(datetime = `Date Time, GMT-04:00`,
         airt_F = `Temp, °F (LGR S/N: 20205284, SEN S/N: 20168088, LBL: Air Temperature)`,
         rh = `RH, % (LGR S/N: 20205284, SEN S/N: 20168088)`,
         PAR = `PAR, µmol/m²/s (LGR S/N: 20205284, SEN S/N: 20189138)`,
         soilt_F = `Temp, °F (LGR S/N: 20205284, SEN S/N: 20245643, LBL: Soil Temperature)`) %>% 
  # convert datettime from character to datetime
  mutate(datetime = as_datetime(datetime, 
                                format = "%m/%d/%Y %H:%M",
                                tz = NULL)) %>% 
  mutate(date = date(datetime),
         hour = hour(datetime)) %>% 
  # convert F to C for both temps:
  mutate(soilt_C = (soilt_F - 32) * 5/9,
         airt_C = (airt_F - 32)*5/9) %>% 
  select(!c(soilt_F, airt_F)) %>% 
  # average by the hour
  group_by(date, hour) %>% 
  summarise(across(is.numeric, mean, na.rm = TRUE)) %>% 
  ungroup() %>% 
  # add a "management" column (here, organic)
  mutate(management = "conventional")

# ============================================================================
# 2. upload and summarize tower-based PAR data: 2019 (very different formatting)
# ============================================================================

# lots of individual files with different names; remove those with confusing names

read_filtered_csvs <- function(subfolder, 
                               required_strings, 
                               excluded_string,
                               case_sensitive = FALSE) {
  # Construct the full path to the subfolder
  folder_path <- here(subfolder)
  
  # Check if folder exists
  if (!dir.exists(folder_path)) {
    stop("Subfolder does not exist: ", folder_path)
  }
  
  # Get all CSV files in the subfolder
  csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    warning("No CSV files found in ", folder_path)
    return(NULL)
  }
  
  # Filter files based on criteria and group by required string
  file_groups <- lapply(required_strings, function(req_string) {
    matching_files <- sapply(csv_files, function(file) {
      filename <- basename(file)
      
      # Check if excluded string is present
      if (grepl(excluded_string, filename, ignore.case = !case_sensitive)) {
        return(FALSE)
      }
      
      # Check if this required string is present
      has_this_string <- grepl(req_string, filename, ignore.case = !case_sensitive)
      
      # Count how many OTHER required strings are present
      other_strings <- required_strings[required_strings != req_string]
      has_other_strings <- any(sapply(other_strings, function(s) {
        grepl(s, filename, ignore.case = !case_sensitive)
      }))
      
      # Must have this string and NOT have any other required strings
      return(has_this_string && !has_other_strings)
    })
    
    csv_files[matching_files]
  })
  
  names(file_groups) <- required_strings
  
  # Read and combine files for each group
  result_list <- lapply(required_strings, function(req_string) {
    files <- file_groups[[req_string]]
    
    if (length(files) == 0) {
      warning("No files found for string: ", req_string)
      return(NULL)
    }
    
    # Read all files in this group
    data_list <- lapply(files, read.csv)
    
    # Find common columns across all data frames in this group
    all_columns <- lapply(data_list, names)
    common_columns <- Reduce(intersect, all_columns)
    
    if (length(common_columns) == 0) {
      warning("No common columns found for files with string: ", req_string)
      return(NULL)
    }
    
    # Select only common columns from each data frame
    data_list <- lapply(data_list, function(df) {
      df[, common_columns, drop = FALSE]
    })
    
    # Combine using rbind
    combined_data <- do.call(rbind, data_list)
    
    # Add column indicating which required string this came from
    combined_data$category <- req_string
    
    # Reset row names
    rownames(combined_data) <- NULL
    
    return(combined_data)
  })
  
  names(result_list) <- required_strings
  
  # Remove NULL entries (groups with no files)
  result_list <- result_list[!sapply(result_list, is.null)]
  
  return(result_list)
}

# use custom function:
test <- read_filtered_csvs(
  subfolder = "pre2020 par etc. data",
  required_strings = c("EF01", "EF02"),
  excluded_string = c("EF03", "dairy"),
  case_sensitive = FALSE)
list2env(test, .GlobalEnv)

# rename colnames, dataframes:
parorg_2019 <- EF01 %>% 
  rename(datetime = datetime,
         airt_F = airtemp_F,
         PAR = par_umolm2sec,
         soilt_F = soiltemp_F) %>% 
  # convert datettime from character to datetime
  mutate(datetime = as_datetime(datetime, 
                                format = "%m/%d/%y %H:%M",
                                tz = NULL)) %>% 
  mutate(date = date(datetime),
         hour = hour(datetime)) %>% 
  # convert F to C for both temps:
  mutate(soilt_C = (soilt_F - 32) * 5/9,
         airt_C = (airt_F - 32)*5/9) %>% 
  select(!c(soilt_F, airt_F)) %>% 
  # average by the hour
  group_by(date, hour) %>% 
  summarise(across(is.numeric, mean, na.rm = TRUE)) %>% 
  mutate(across(everything(), ~ifelse(is.nan(.), NA, .))) %>% 
  ungroup() %>% 
  # add a "management" column (here, organic)
  mutate(management = "organic")

parconv_2019 <- EF02 %>% 
  rename(datetime = datetime,
         airt_F = airtemp_F,
         PAR = par_umolm2sec,
         soilt_F = soiltemp_F) %>% 
  # convert datettime from character to datetime
  mutate(datetime = as_datetime(datetime, 
                                format = "%m/%d/%y %H:%M",
                                tz = NULL)) %>% 
  mutate(date = date(datetime),
         hour = hour(datetime)) %>% 
  # convert F to C for both temps:
  mutate(soilt_C = (soilt_F - 32) * 5/9,
         airt_C = (airt_F - 32)*5/9) %>% 
  select(!c(soilt_F, airt_F)) %>% 
  # average by the hour
  group_by(date, hour) %>% 
  summarise(across(is.numeric, mean, na.rm = TRUE)) %>% 
  mutate(across(everything(), ~ifelse(is.nan(.), NA, .))) %>% 
  ungroup() %>% 
  # add a "management" column (here, conventional)
  mutate(management = "conventional")

# ============================================================================
# 3. merge all hourly averages of PAR using rbind (long form), calculate daily max
# ============================================================================

pardat_hourly <- rbind(parorg_2019, parconv_2019, parorg_2020, parconv_2020)

pardat_dailymax <- pardat_hourly %>% 
  group_by(management, date) %>% 
    summarize(
      across(is.numeric, max, na.rm = TRUE),
      .groups = 'drop'
    ) %>% 
      ungroup() %>% 
  select(!c(rh, hour, soilt_C, airt_C))

# ============================================================================
# 4. remove unnecessary items from global environment
# ============================================================================

rm(parorg_2019, parorg_2020, parconv_2019, parconv_2020, test,
   EF01, EF02)