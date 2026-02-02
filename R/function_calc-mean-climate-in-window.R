## Function to calculate mean of climate variable in a specified window ####

# Author:           Cherine C. Jantzen
# Created:          2025-05-26
# Last updated:     2025-05-26

################################
## Content: This function is used to easily calculate the mean or sum of daily values of a climate variable for a selected window of the year.
################################

# load packages
library(tibble)
library(dplyr)
library(lubridate)


# Arguments:
# - winStart_DOY: integer, specifying day of year the window starts
# - winEnd_DOY: integer, specifying the day of year the window ends
# - clim_var: string, specifying the climate variable included in df_climate for which the mean should be calculated
# - method: one of "mean" or "sum", defines whether to calculate the mean of the daily values in the window or the sum


# Function:
calc_mean.climate_in_window <- function(winStart_DOY,
                                        winEnd_DOY,
                                        clim_var,
                                        method = c("mean", "sum")){
  
  # create dates out of DOY (based on a non-leap year)
  win_dates <- tibble::tibble(date = c(as.Date(winStart_DOY-1, origin = "2025-01-01"), 
                                       as.Date(winEnd_DOY-1, origin = "2025-01-01"))) %>% 
    dplyr::mutate(month = lubridate::month(date), 
                  day = lubridate::day(date),
                  month_day_ID = month * 100 + day)
  
  if (method == "mean") {
    
    df_climate %>%
      dplyr::mutate(month_day_ID = month * 100 + day) %>% 
      dplyr::filter(dplyr::between(month_day_ID, min(win_dates$month_day_ID), max(win_dates$month_day_ID))) %>% 
      dplyr::rename("variable" = all_of(clim_var)) %>% 
      dplyr::summarise(mean_clim_var = mean(variable, na.rm = TRUE), .by = year) %>% 
      dplyr::mutate(yearPlus1 = year + 1,
                    yearPlus2 = year + 2)
    
  } else if(method == "sum") {
    
    df_climate %>%
      dplyr::mutate(month_day_ID = month * 100 + day) %>% 
      dplyr::filter(dplyr::between(month_day_ID, min(win_dates$month_day_ID), max(win_dates$month_day_ID))) %>% 
      dplyr::rename("variable" = all_of(clim_var)) %>% 
      dplyr::summarise(sum_clim_var = sum(variable, na.rm = TRUE), .by = year) %>% 
      dplyr::mutate(yearPlus1 = year + 1,
                    yearPlus2 = year + 2)
    
  }
}
