## Data preparation for analysis - DwC-A of beechnuts and climate variables ####

# Author:           Cherine C. Jantzen
# Created:          2024-11-20
# Last updated:     2026-02-02

################################
## Content: This script prepares the climate data for analysis and combines it with the annual beechnut counts per tree 
## into the final input file that is used for all further analyses.
################################

# load libraries
library(dplyr)
library(here)
library(stringr)

## get function to calculate climate variables
source(here::here("R", "function_calc-mean-climate-in-window.R"))

## get climate data from KNMI (deBilt_1974_2025-20251104.txt)
knmi <- read.csv(here::here("data", "deBilt_1974_2025-20251104.txt"))

## load individual tree input data (number of nuts per category per tree per year)
nuts <- read.csv(here::here("data", "beechnut-count-data.csv"))


# Climate variables -------------------------------------------------------

# prepare climate (i.e., knmi) data for analysis
## all climate values are stored as 0.1 degrees (temp) or mm (precipitation) and have to be multiplied by 10 to get degrees and mm

df_climate <- knmi %>% 
  dplyr::mutate(year = as.integer(substring(YYYYMMDD, 1, 4)),
                month = as.integer(substring(YYYYMMDD, 5, 6)),
                day = as.integer(substring(YYYYMMDD, 7, 8)),
                meanDayTemp = (TN + TX)/20, # mean daily temperature in degree Celsius
                minDayTemp = TN/10, # minimum daily temperature in degree Celsius
                maxDayTemp = TX/10, # maximum daily temperature in degree Celsius
                raw_sumDayPrec = RH/10, # daily sum of precipitation in millimetre
  ) %>% 
  # if daily precipitation sum <0.05 mm, data is encoded as -1 by data provider. To avoid negative values when calculating sums, we set those values to zero
  dplyr::mutate(sumDayPrec = dplyr::case_when(raw_sumDayPrec == -0.1 ~ 0,
                                              TRUE ~ raw_sumDayPrec)) %>% 
  dplyr::select("year", "month", "day", "meanDayTemp", "minDayTemp", "maxDayTemp", "sumDayPrec")

# save processed climate data for later analyses
write.csv(df_climate, file = here::here("data", "df_climate-processed-knmi-data.csv"), row.names = FALSE)


## Calculate temperatures and precipitation in certain windows of the year, as used in the literature ####
### Note that all values are z-scored to put them on a similar scale for later modelling


## Temperatures ####

# calculate DOY of summer window from literature (01/06 to 31/07)
startDOY_summer <- lubridate::yday("2025-06-01")
endDOY_summer <- lubridate::yday("2025-07-31")

# calculate mean maximum summer temperatures (June & July)
temp_meanMaxSummer <- calc_mean.climate_in_window(winStart_DOY = startDOY_summer,
                                                  winEnd_DOY = endDOY_summer,
                                                  method = "mean",
                                                  clim_var = "maxDayTemp") %>%
  dplyr::rename("temp_meanMaxSummer" = "mean_clim_var") %>% 
  dplyr::mutate(z_maxTempSummer = scale(temp_meanMaxSummer, center = TRUE, scale = TRUE))


# calculate mean daily temperatures (June & July) 
temp_meanDaySummer <- calc_mean.climate_in_window(winStart_DOY = startDOY_summer,
                                                  winEnd_DOY = endDOY_summer,
                                                  method = "mean",
                                                  clim_var = "meanDayTemp") %>%
  dplyr::rename("temp_meanDaySummer" = "mean_clim_var") %>% 
  dplyr::mutate(z_meanTempSummer = scale(temp_meanDaySummer, center = TRUE, scale = TRUE))


# calculate mean minimum temperature (June & July)
temp_meanMinSummer <- calc_mean.climate_in_window(winStart_DOY = startDOY_summer,
                                                  winEnd_DOY = endDOY_summer,
                                                  method = "mean",
                                                  clim_var = "minDayTemp") %>%
  dplyr::rename("temp_meanMinSummer" = "mean_clim_var") %>% 
  dplyr::mutate(z_minTempSummer = scale(temp_meanMinSummer, center = TRUE, scale = TRUE))



## Precipitation ####

# calculate the mean precipitation in summer (June & July)
prec_Summer <- calc_mean.climate_in_window(winStart_DOY = startDOY_summer,
                                           winEnd_DOY = endDOY_summer,
                                           method = "sum",
                                           clim_var = "sumDayPrec") %>%
  dplyr::rename("prec_Summer" = "sum_clim_var") %>% 
  dplyr::mutate(z_prec_Summer = scale(prec_Summer, center = TRUE, scale = TRUE))



## Growing season ####

# calculate temperature in the growing season (May to August)
temp_meanGrow <- calc_mean.climate_in_window(winStart_DOY = lubridate::yday("2025-05-01"),
                                             winEnd_DOY = lubridate::yday("2025-08-31"),
                                             method = "mean",
                                             clim_var = "meanDayTemp") %>%
  dplyr::rename("temp_meanGrow" = "mean_clim_var") %>% 
  dplyr::mutate(z_temp_meanGrow = scale(temp_meanGrow, center = TRUE, scale = TRUE))


# calculate the mean precipitation in growing season (March to April)
prec_Grow <- calc_mean.climate_in_window(winStart_DOY = lubridate::yday("2025-03-01"),
                                         winEnd_DOY = lubridate::yday("2025-04-30"),
                                         method = "sum",
                                         clim_var = "sumDayPrec") %>%
  dplyr::rename("prec_Grow" = "sum_clim_var") %>% 
  dplyr::mutate(z_prec_Grow = scale(prec_Grow, center = TRUE, scale = TRUE))


# combine all climate variables in one data frame
climate_variables <- temp_meanDaySummer %>% 
  dplyr::select("year") %>% 
  dplyr::left_join(temp_meanDaySummer %>% 
                     dplyr::select("z_meanTempSummer_T1" = "z_meanTempSummer", "yearPlus1"), 
                   by = c("year" = "yearPlus1")) %>% 
  dplyr::left_join(temp_meanDaySummer %>% 
                     dplyr::select("z_meanTempSummer_T2" = "z_meanTempSummer", "yearPlus2"), 
                   by = c("year" = "yearPlus2")) %>% 
  dplyr::left_join(temp_meanMinSummer %>% 
                     dplyr::select("z_minTempSummer_T1" = "z_minTempSummer", "yearPlus1"), 
                   by = c("year" = "yearPlus1")) %>% 
  dplyr::left_join(temp_meanMinSummer %>% 
                     dplyr::select("z_minTempSummer_T2" = "z_minTempSummer", "yearPlus2"), 
                   by = c("year" = "yearPlus2")) %>% 
  dplyr::left_join(temp_meanMaxSummer %>% 
                     dplyr::select("z_maxTempSummer_T1" = "z_maxTempSummer", "yearPlus1"), 
                   by = c("year" = "yearPlus1")) %>% 
  dplyr::left_join(temp_meanMaxSummer %>% 
                     dplyr::select("z_maxTempSummer_T2" = "z_maxTempSummer", "yearPlus2"), 
                   by = c("year" = "yearPlus2")) %>% 
  dplyr::left_join(temp_meanGrow %>% 
                     dplyr::select("year", "z_temp_meanGrow"), by = "year") %>% 
  dplyr::left_join(prec_Summer %>% 
                     dplyr::select("z_prec_Summer_T1" = "z_prec_Summer", "yearPlus1"), 
                   by = c("year" = "yearPlus1")) %>% 
  dplyr::left_join(prec_Summer %>% 
                     dplyr::select("z_prec_Summer_T2" = "z_prec_Summer", "yearPlus2"), 
                   by = c("year" = "yearPlus2")) %>% 
  dplyr::left_join(prec_Grow %>% 
                     dplyr::select("year", "z_prec_Grow"), by = "year") 


# Combine seed production & climate data ----------------------------

# bind individual tree data to climate variables
di <- nuts %>% 
  dplyr::left_join(climate_variables, by = c("WinterYear" = "year")) %>% 
  # floor all counts to get integers needed for models (some values in the data are calculated, causing non-integer values)
  dplyr::mutate(TotalWhole = floor(TotalWhole),
                TotalEmpty = floor(TotalEmpty),
                TotalNuts = floor(TotalNuts),
                TotalPredated = floor(TotalPredated))

# save final input data frame used for all following analysis

## select directory to store data frame
write.csv(di, file = here::here("data", "di_input-data.csv"), row.names = FALSE)
