# Sensitivity windows: Iteratiive sliding window approach ####

# Author:         Cherine C. Jantzen
# Created:        2025-10-24
# Last updated:   2025-11-08

########################################
## Content: This script finds the weather cues of beechnut production by using an iterative sliding window approach.
##
## NOTE: This script was run on an external server and runs models in parallel to reduce running time. As running each call of "prepare_and_run_model"
## takes around 14 hours, results of each batch and the the summary of all models are written to file. All necessary files for further processing the results are also
## written to file and further used in another script.
########################################


# load packages
library(dplyr)
library(purrr)
library(here)
library(glmmTMB)
library(stringr)
library(future)
library(future.apply)
library(future.callr)


# Data preparation --------------------------------------------------------

# read in data
df_climate <- read.csv(here::here("data", "df_climate-processed-knmi-data.csv"))
di <- read.csv(here::here("data", "di_input-data.csv"))

# source functions
source(here::here("R", "functions_sensitivity-window-analyses.R"))

# create folders to store results
if (!dir.exists(here::here("plots"))) {dir.create(here::here("plots"))}
if (!dir.exists(here::here("results"))) {dir.create(here::here("results"))}

# get seed production from previous year
bci_t1 <- di %>% 
  dplyr::select(TreeID, WinterYear, TotalNuts) %>% 
  dplyr::mutate(year_plus1 = WinterYear + 1) %>% 
  dplyr::rename("TotalNuts_T1" = "TotalNuts") %>% 
  dplyr::select(!WinterYear)

# subset di to only necessary columns for analysis
di_sub <- di %>%
  dplyr::left_join(bci_t1, by = c("TreeID" = "TreeID", "WinterYear" = "year_plus1")) %>% 
  dplyr::arrange(WinterYear) %>%
  dplyr::mutate(TreeID = as.factor(TreeID)) %>%
  dplyr::select(!c("TotalWhole", "TotalEmpty", "TotalPredated"))

# Set up batching parameters for running in parallel
no_batches <- 450
batch_size <- 35


## Create climate windows -----------------------------------------------

# calculate climate variables for all possible windows between spring equinox and autumn equinox (DOY 79 to DOY 265)
# commented out as the output file is already read in in l.34 to speed up the script (necessary to run, when running script for the first time)

all_windows <- purrr::map(.x = c(7:140), # minimum window length 7 days, maximum length 140 days (1 week to 20 weeks)
                          .f = ~{

                            window_length <- .x


                            # calculate the mean daily climate variables per window per year
                            temp_windows <- purrr::map(.x = (79:((265 + 1) - window_length)),
                                                       .f = ~{

                                                         df_climate %>%
                                                           dplyr::group_by(year) %>%
                                                           # slice data in respective window
                                                           dplyr::slice(.x:(.x + (window_length - 1))) %>%
                                                           dplyr::summarise(meanMeanWinTemp = mean(meanDayTemp, na.rm = TRUE),
                                                                            meanMaxWinTemp = mean(maxDayTemp, na.rm = TRUE),
                                                                            meanMinWinTemp = mean(minDayTemp, na.rm = TRUE),
                                                                            sumWinPrec = sum(sumDayPrec, na.rm = TRUE)) %>%
                                                           dplyr::mutate(Start_DOY = .x,
                                                                         End_DOY = .x + (window_length - 1),
                                                                         windowID = paste(Start_DOY, End_DOY, sep = "_"))

                                                       }
                            ) %>%  dplyr::bind_rows()
                          }, .progress = TRUE) %>%
  dplyr::bind_rows()

# write to file to avoid re-running
write.csv(all_windows, file = here::here("data", "all_windows_climate.csv"), row.names = FALSE)

# z-score the variables
z_all_windows <- all_windows %>% 
  dplyr::mutate(z_meanMeanWinTemp = scale(meanMeanWinTemp, center = TRUE, scale = TRUE)[,1], 
                z_meanMaxWinTemp = scale(meanMaxWinTemp, center = TRUE, scale = TRUE)[,1],
                z_meanMinWinTemp = scale(meanMinWinTemp, center = TRUE, scale = TRUE)[,1],
                z_sumWinPrec = scale(sumWinPrec, center = TRUE, scale = TRUE)[,1],
                .by = windowID) %>% 
  dplyr::select(!c("meanMeanWinTemp", "meanMaxWinTemp", "meanMinWinTemp", "sumWinPrec"))



# Model temperatures for T1 - First iteration -----------------------------

# join beech crop data with climate windows of previous year
di_win <- di_sub %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T1 = year + 1),
                   by = c("WinterYear" = "year_T1"), relationship = "many-to-many")

## Create model functions ####

## model function for max temperature in T1 ----
model_maxTempT1 <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanGrow + z_meanMaxWinTemp + z_maxTempSummer_T2 +
                               z_prec_Grow + z_prec_Summer_T1  + z_prec_Summer_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}


## model function for min temperature in T1 ----
model_minTempT1 <- function(datfr) {

  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanGrow + z_meanMinWinTemp + z_maxTempSummer_T2 +
                               z_prec_Grow + z_prec_Summer_T1  + z_prec_Summer_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## model function for mean teamperature in T1 ----
model_meanTempT1 <- function(datfr) {

  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanGrow + z_meanMeanWinTemp + z_maxTempSummer_T2 +
                               z_prec_Grow + z_prec_Summer_T1  + z_prec_Summer_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run models --------------------------------------------------------------

# First iteration: decide which temperature variable we use for the remaining models. (min, max, or mean temperature)

prepare_and_run_model(model = model_maxTempT1,
                      variable_name = "z_maxTempT1")

prepare_and_run_model(model = model_minTempT1,
                      variable_name = "z_minTempT1")

prepare_and_run_model(model = model_meanTempT1,
                      variable_name = "z_meanTempT1")

## Look at model outputs ####

### max temperature ####

summary_maxTempT1 <- get(load(here::here("results", "z_maxTempT1/summary_z_maxTempT1.rda")))

# # plot results
plots <- plot_results(dataframe = summary_maxTempT1,
             plot_title = "z_maxTempT1")

# heat map
plots[1]

# AIC peak
plots[2]

# ### min temperature ####
summary_minTempT1 <- get(load(here::here("results", "z_minTempT1/summary_z_minTempT1.rda")))

# plot results
plots <- plot_results(dataframe = summary_minTempT1,
             plot_title = "z_minTempT1")

# heat map
plots[1]

# AIC peak
plots[2]

# ### mean temperature ####
summary_meanTempT1 <- get(load(here::here("results", "z_meanTempT1/summary_z_meanTempT1.rda")))

# plot results
plots <- plot_results(dataframe = summary_meanTempT1,
             plot_title = "z_meanTempT1")

# heat map
plots[1]

# AIC peak
plots[2]

# compare model outcomes for three different temperature variables
summary_maxTempT1 %>% dplyr::filter(deltaAIC == 0) %>% dplyr::pull(AIC, variable)
summary_minTempT1 %>% dplyr::filter(deltaAIC == 0) %>% dplyr::pull(AIC, variable)
summary_meanTempT1 %>% dplyr::filter(deltaAIC == 0) %>% dplyr::pull(AIC, variable)
# max temperature has clearly lowest AIC (deltaAIC 40.73 to mean and 5.38 to min temperature)

# get window days of best model
best_wind_maxTempT1 <- summary_maxTempT1 %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate temperature for maxTemp in T1 for selected window
maxTempT1_firstIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_maxTempT1$Start_DOY, End_DOY == best_wind_maxTempT1$End_DOY) %>%
  dplyr::mutate(yearPlus1 = year + 1) %>%
  dplyr::select("yearPlus1", "z_temp_meanMax_firstIt_T1" = "z_meanMaxWinTemp")

# add max temperature of T1 for selected window
di_win_new <- di_sub %>%
  dplyr::left_join(maxTempT1_firstIt, by = c("WinterYear" = "yearPlus1"))


# # Model temperatures for T2 - First iteration -----------------------------

# combine beechnut data with windows in year - 2
di_win <- di_win_new %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T2 = year + 2),
                   by = c("WinterYear" = "year_T2"), relationship = "many-to-many") %>%
  dplyr::select(!c("year"))

## model function for max temperature in T2 ----------------------------------------

model_maxTempT2 <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanGrow + z_temp_meanMax_firstIt_T1 + z_meanMaxWinTemp +
                               z_prec_Grow + z_prec_Summer_T1  + z_prec_Summer_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## model function for min temperature in T2 ----------------------------------------

model_minTempT2 <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanGrow + z_temp_meanMax_firstIt_T1 + z_meanMinWinTemp +
                               z_prec_Grow + z_prec_Summer_T1  + z_prec_Summer_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## model function for mean temperature in T2 ----------------------------------------

model_meanTempT2 <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanGrow + z_temp_meanMax_firstIt_T1 + z_meanMeanWinTemp +
                               z_prec_Grow + z_prec_Summer_T1  + z_prec_Summer_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run models ####

# #max temp
prepare_and_run_model(model = model_maxTempT2,
                      variable_name = "z_maxTempT2")

# # min temp
prepare_and_run_model(model = model_minTempT2,
                      variable_name = "z_minTempT2")

# # mean temp
prepare_and_run_model(model = model_meanTempT2,
                      variable_name = "z_meanTempT2")

## Look at model output ####

### max temperature ####
summary_maxTempT2 <- get(load(here::here("results", "z_maxTempT2/summary_z_maxTempT2.rda")))

# plot results
plots <- plot_results(dataframe = summary_maxTempT2,
             plot_title = "z_maxTempT2")

# heat map
plots[1]

# AIC peak
plots[2]

### min temperature ####
summary_minTempT2 <- get(load(here::here("results", "z_minTempT2/summary_z_minTempT2.rda")))

# plot results
plots <- plot_results(dataframe = summary_minTempT2,
             plot_title = "z_minTempT2")

# heat map
plots[1]

# AIC peak
plots[2]

## mean temperature ####
summary_meanTempT2 <- get(load(here::here("results", "z_meanTempT2/summary_z_meanTempT2.rda")))

# plot results
plots <- plot_results(dataframe = summary_meanTempT2,
             plot_title = "z_meanTempT2")

# heat map
plots[1]

# AIC peak
plots[2]

# compare model outcomes for three different temperature variables
summary_maxTempT2 %>% dplyr::filter(deltaAIC == 0) %>% dplyr::pull(AIC, variable)
summary_minTempT2 %>% dplyr::filter(deltaAIC == 0) %>% dplyr::pull(AIC, variable)
summary_meanTempT2 %>% dplyr::filter(deltaAIC == 0) %>% dplyr::pull(AIC, variable)
# max temperature has lowest AIC (deltaAIC 33.95 (mean) and 15.06 (min))

## decide for max temperature and add to data frame
best_wind_maxTempT2 <- summary_maxTempT2 %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate temperature for max Temp in T2 for selected window
maxTempT2_firstIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_maxTempT2$Start_DOY, End_DOY == best_wind_maxTempT2$End_DOY) %>%
  dplyr::mutate(yearPlus2 = year + 2) %>%
  dplyr::select("yearPlus2", "z_temp_meanMax_firstIt_T2" = "z_meanMaxWinTemp")

# add T2 temperatures
di_win_new <- di_win_new %>%
  dplyr::left_join(maxTempT2_firstIt, by = c("WinterYear" = "yearPlus2"))


# Model temperatures for T0 - First iteration -----------------------------

# Combine beechnut data with climate windows of year of seed fall
di_win <- di_win_new  %>%
  dplyr::left_join(z_all_windows,
                   by = c("WinterYear" = "year"), relationship = "many-to-many")

## model function for mean temperature in T0 ----------------------------

# mean temperature
model_meanTempT0 <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_meanMeanWinTemp + z_temp_meanMax_firstIt_T1 + z_temp_meanMax_firstIt_T2 +
                               z_prec_Grow + z_prec_Summer_T1  + z_prec_Summer_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## model function for max temperature in T0 ----------------------------

# mean temperature
model_maxTempT0 <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_meanMaxWinTemp + z_temp_meanMax_firstIt_T1 + z_temp_meanMax_firstIt_T2 +
                               z_prec_Grow + z_prec_Summer_T1  + z_prec_Summer_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## model function for min temperature in T0 ----------------------------

# mean temperature
model_minTempT0 <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_meanMinWinTemp + z_temp_meanMax_firstIt_T1 + z_temp_meanMax_firstIt_T2 +
                               z_prec_Grow + z_prec_Summer_T1  + z_prec_Summer_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run models ####

# mean temperature
prepare_and_run_model(model = model_meanTempT0,
                      variable_name = "z_meanTempT0")

# max temperature
prepare_and_run_model(model = model_maxTempT0,
                      variable_name = "z_maxTempT0")

# min temperature
prepare_and_run_model(model = model_minTempT0,
                      variable_name = "z_minTempT0")

## Look at model output ####

### mean temperature  ####
summary_meanTempT0 <- get(load(here::here("results", "z_meanTempT0/summary_z_meanTempT0.rda")))

# plot results
plots <- plot_results(dataframe = summary_meanTempT0,
             plot_title = "z_meanTempT0")

# heat map
plots[1]

# AIC peak
plots[2]

### max temperature  ####
summary_maxTempT0 <- get(load(here::here("results", "z_maxTempT0/summary_z_maxTempT0.rda")))

# plot results
plots <- plot_results(dataframe = summary_maxTempT0,
             plot_title = "z_maxTempT0")

# heat map
plots[1]

# AIC peak
plots[2]

### min temperature  ####
summary_minTempT0 <- get(load(here::here("results", "z_minTempT0/summary_z_minTempT0.rda")))

# plot results
plots <- plot_results(dataframe = summary_minTempT0,
             plot_title = "z_minTempT0")

# heat map
plots[1]

# AIC peak
plots[2]

# compare model outcomes for three different temperature variables
summary_maxTempT0 %>% dplyr::filter(deltaAIC == 0) %>% dplyr::pull(AIC, variable)
summary_minTempT0 %>% dplyr::filter(deltaAIC == 0) %>% dplyr::pull(AIC, variable)
summary_meanTempT0 %>% dplyr::filter(deltaAIC == 0) %>% dplyr::pull(AIC, variable)
# mean temperature is best suited; delta AIC: 19.49 (max) & 51.26 (min)

# decide for best variable and add temperatures
best_wind_meanTempT0 <- summary_meanTempT0 %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate temperature for meanTemp in T0 for selected window
meanTempT0_firstIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_meanTempT0$Start_DOY, End_DOY == best_wind_meanTempT0$End_DOY) %>%
  dplyr::select("year", "z_temp_meanMean_firstIt_T0" = "z_meanMeanWinTemp")

# add T0 temperatures
di_win_new <- di_win_new %>%
  dplyr::left_join(meanTempT0_firstIt, by = c("WinterYear" = "year"))


# Model precipitation for T1 - First iteration ----------------------------

# join beech crop data with climate windows of one years prior to seed fall
di_win <- di_win_new %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T1 = year + 1),
                   by = c("WinterYear" = "year_T1"), relationship = "many-to-many")

## model function for precipitation in T1 ####

model_PrecT1 <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_firstIt_T0 + z_temp_meanMax_firstIt_T1 + z_temp_meanMax_firstIt_T2 +
                               z_prec_Grow + z_sumWinPrec + z_prec_Summer_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

# ## Run model ####
prepare_and_run_model(model = model_PrecT1,
                      variable_name = "z_PrecT1")

## Look at model output ####
summary_PrecT1 <- get(load(here::here("results", "z_PrecT1/summary_z_PrecT1.rda")))

# plot results
plots <- plot_results(dataframe = summary_PrecT1,
             plot_title = "z_PrecT1")

# heat map
plots[1]

# AIC peak
plots[2]

# add values for best model to data
best_wind_PrecT1 <- summary_PrecT1 %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate precipitation sum in T1 for selected window
precT1_firstIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_PrecT1$Start_DOY, End_DOY == best_wind_PrecT1$End_DOY) %>%
  dplyr::mutate(yearPlus1 = year + 1) %>%
  dplyr::select("yearPlus1", "z_prec_firstIt_T1" = "z_sumWinPrec")

# add T1 precipitations
di_win_new <- di_win_new %>%
  dplyr::left_join(precT1_firstIt, by = c("WinterYear" = "yearPlus1"))


# Model precipitation for T2 - First iteration ----------------------------

# join beech crop data with climate windows of two years prior to seed fall
di_win <- di_win_new  %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T2 = year + 2),
                   by = c("WinterYear" = "year_T2"), relationship = "many-to-many")

## model function for precipitation in T2  ----------------------------

model_PrecT2 <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_firstIt_T0 + z_temp_meanMax_firstIt_T1 + z_temp_meanMax_firstIt_T2 +
                               z_prec_Grow + z_prec_firstIt_T1 + z_sumWinPrec +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run model ####
prepare_and_run_model(model = model_PrecT2,
                      variable_name = "z_PrecT2")

## Look at output ####
summary_PrecT2 <- get(load(here::here("results", "z_PrecT2/summary_z_PrecT2.rda")))

# plot results
plots <- plot_results(dataframe = summary_PrecT2,
             plot_title = "z_precT2")

# heat map
plots[1]

# AIC peak
plots[2]

# add precipitation of best window to data
best_wind_PrecT2 <- summary_PrecT2 %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate precipitation in T2 for selected window
precT2_firstIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_PrecT2$Start_DOY, End_DOY == best_wind_PrecT2$End_DOY) %>%
  dplyr::mutate(yearPlus2 = year + 2) %>%
  dplyr::select("yearPlus2", "z_prec_firstIt_T2" = "z_sumWinPrec")

# add T2 precipitations
di_win_new <- di_win_new %>%
  dplyr::left_join(precT2_firstIt, by = c("WinterYear" = "yearPlus2"))

# Model precipitation in T0 - First iteration -----------------------------

# join beech crop data with climate windows of year of seed fall
di_win <- di_win_new %>%
  dplyr::left_join(z_all_windows,
                   by = c("WinterYear" = "year"), relationship = "many-to-many")

## model precipitation in T0 ####

model_PrecT0 <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_firstIt_T0 + z_temp_meanMax_firstIt_T1 + z_temp_meanMax_firstIt_T2 +
                               z_sumWinPrec + z_prec_firstIt_T1 + z_prec_firstIt_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run model ####
prepare_and_run_model(model = model_PrecT0,
                      variable_name = "z_PrecT0")

## Look at model output ####
summary_PrecT0 <- get(load(here::here("results", "z_PrecT0/summary_z_PrecT0.rda")))

# plot results
plots <- plot_results(dataframe = summary_PrecT0,
             plot_title = "z_precT0")

# heat map
plots[1]

# AIC peak
plots[2]

# add variable from best window to data
best_wind_PrecT0 <- summary_PrecT0 %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate precipitation for T0 for selected window
precT0_firstIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_PrecT0$Start_DOY, End_DOY == best_wind_PrecT0$End_DOY) %>%
  dplyr::select("year", "z_prec_firstIt_T0" = "z_sumWinPrec")

# add T0 precipitations
di_win_new <- di_win_new %>%
  dplyr::left_join(precT0_firstIt, by = c("WinterYear" = "year"))



# II. Second iteration ----------------------------------------------------

# delete old window data from input data frame for second iteration to save size
di_win_secIt <- di_win_new %>%
  dplyr::select(!c("z_temp_meanGrow", "z_maxTempSummer_T1", "z_maxTempSummer_T2",
                   "z_prec_Grow", "z_prec_Summer_T1", "z_prec_Summer_T2"))

# Model temperature in T1 - Second iteration -----------------------------------

di_win <- di_win_secIt %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T1 = year + 1),
                   by = c("WinterYear" = "year_T1"), relationship = "many-to-many")

## model function for max temperature in T1 ----
model_maxTempT1_secIt <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_firstIt_T0 + z_meanMaxWinTemp + z_temp_meanMax_firstIt_T2 +
                               z_prec_firstIt_T0 + z_prec_firstIt_T1 + z_prec_firstIt_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

# ## Run model ####
prepare_and_run_model(model = model_maxTempT1_secIt,
                      variable_name = "z_maxTempT1_secIt")

## Look at model output ####
summary_maxTempT1_secIt <- get(load(here::here("results", "z_maxTempT1_secIt/summary_z_maxTempT1_secIt.rda")))

# plot results
plots <- plot_results(dataframe = summary_maxTempT1_secIt,
             plot_title = "z_maxTempT1_secIt")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_maxTempT1_secIt <- summary_maxTempT1_secIt %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate temperature for maxTemp in T1 for selected window
maxTempT1_secIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_maxTempT1_secIt$Start_DOY, End_DOY == best_wind_maxTempT1_secIt$End_DOY) %>%
  dplyr::mutate(yearPlus1 = year + 1) %>%
  dplyr::select("yearPlus1", "z_temp_meanMax_secIt_T1" = "z_meanMaxWinTemp")

# add T1 temperatures
di_win_secIt <- di_win_secIt %>%
  dplyr::left_join(maxTempT1_secIt, by = c("WinterYear" = "yearPlus1"))


# Model temperature in T2 - Second iteration ------------------------------
di_win <- di_win_secIt %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T2 = year + 2),
                   by = c("WinterYear" = "year_T2"), relationship = "many-to-many")

## model function for temperature in T2 ####
model_maxTempT2_secIt <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_firstIt_T0 + z_temp_meanMax_secIt_T1 + z_meanMaxWinTemp +
                               z_prec_firstIt_T0 + z_prec_firstIt_T1 + z_prec_firstIt_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

# ## Run model ####
prepare_and_run_model(model = model_maxTempT2_secIt,
                      variable_name = "z_maxTempT2_secIt")

## Look at model output ####
summary_maxTempT2_secIt <- get(load(here::here("results", "z_maxTempT2_secIt/summary_maxTempT2_secIt.rda")))

# plot results
plots <- plot_results(dataframe = summary_maxTempT2_secIt,
             plot_title = "z_maxTempT2_secIt")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_maxTempT2_secIt <- summary_maxTempT2_secIt %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate temperature for maxTemp in T2 for selected window
maxTempT2_secIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_maxTempT2_secIt$Start_DOY, End_DOY == best_wind_maxTempT2_secIt$End_DOY) %>%
  dplyr::mutate(yearPlus2 = year + 2) %>%
  dplyr::select("yearPlus2", "z_temp_meanMax_secIt_T2" = "z_meanMaxWinTemp")

# add T2 temperatures
di_win_secIt <- di_win_secIt %>%
  dplyr::left_join(maxTempT2_secIt, by = c("WinterYear" = "yearPlus2"))


# Model temperature in T0 - Second iteration ------------------------------
di_win <- di_win_secIt %>%
  dplyr::left_join(z_all_windows,
                   by = c("WinterYear" = "year"), relationship = "many-to-many")

## Model function ----
model_meanTempT0_secIt <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_meanMeanWinTemp + z_temp_meanMax_secIt_T1 + z_temp_meanMax_secIt_T2 +
                               z_prec_firstIt_T0 + z_prec_firstIt_T1 + z_prec_firstIt_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run model ####
prepare_and_run_model(model = model_meanTempT0_secIt,
                      variable_name = "z_meanTempT0_secIt")

## Look at model output ####
summary_meanTempT0_secIt <- get(load(here::here("results", "z_meanTempT0_secIt/summary_z_meanTempT0_secIt.rda")))

# plot results
plots <- plot_results(dataframe = summary_meanTempT0_secIt,
             plot_title = "z_meanTempT0_secIt")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_meanTempT0_secIt <- summary_meanTempT0_secIt %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate temperature for meanTemp in T0 for selected window
meanTempT0_secIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_meanTempT0_secIt$Start_DOY, End_DOY == best_wind_meanTempT0_secIt$End_DOY) %>%
  dplyr::select("year", "z_temp_meanMean_secIt_T0" = "z_meanMeanWinTemp")

# add T0 temperatures
di_win_secIt <- di_win_secIt %>%
  dplyr::left_join(meanTempT0_secIt, by = c("WinterYear" = "year"))


# Model precipitation in T1 - Second iteration ----------------------------
di_win <- di_win_secIt %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T1 = year + 1),
                   by = c("WinterYear" = "year_T1"), relationship = "many-to-many")

## model function for precipitation in T1 -----
model_PrecT1_secIt <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_secIt_T0 + z_temp_meanMax_secIt_T1 + z_temp_meanMax_secIt_T2 +
                               z_prec_firstIt_T0 + z_sumWinPrec  + z_prec_firstIt_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

# Run model ####
prepare_and_run_model(model = model_PrecT1_secIt,
                      variable_name = "z_PrecT1_secIt")

## Look at model output ####
summary_PrecT1_secIt <- get(load(here::here("results", "z_PrecT1_secIt/summary_z_PrecT1_secIt.rda")))

# plot results
plots <- plot_results(dataframe = summary_PrecT1_secIt,
             plot_title = "z_precT1_secIt")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_PrecT1_secIt <- summary_PrecT1_secIt %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate temperature for precipitation in T1 for selected window
PrecT1_secIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_PrecT1_secIt$Start_DOY, End_DOY == best_wind_PrecT1_secIt$End_DOY) %>%
  dplyr::mutate(yearPlus1 = year + 1) %>%
  dplyr::select("yearPlus1", "z_prec_secIt_T1" = "z_sumWinPrec")

# add T1 temperatures
di_win_secIt <- di_win_secIt %>%
  dplyr::left_join(PrecT1_secIt, by = c("WinterYear" = "yearPlus1"))


# Model precipitation in T2 - Second iteration ----------------------------
di_win <- di_win_secIt %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T2 = year + 2),
                   by = c("WinterYear" = "year_T2"), relationship = "many-to-many")

## model function precipitation in T2 ----
model_PrecT2_secIt <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_secIt_T0 + z_temp_meanMax_secIt_T1 + z_temp_meanMax_secIt_T2 +
                               z_prec_firstIt_T0  + z_prec_secIt_T1 +  z_sumWinPrec +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run model ####
prepare_and_run_model(model = model_PrecT2_secIt,
                      variable_name = "z_PrecT2_secIt")

## Look at model output ####
summary_PrecT2_secIt <- get(load(here::here("results", "z_PrecT2_secIt/summary_z_PrecT2_secIt.rda")))

# plot results
plots <- plot_results(dataframe = summary_PrecT2_secIt,
             plot_title = "z_precT2_secIt")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_PrecT2_secIt <- summary_PrecT2_secIt %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate precipitation in T2 for selected window
PrecT2_secIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_PrecT2_secIt$Start_DOY, End_DOY == best_wind_PrecT2_secIt$End_DOY) %>%
  dplyr::mutate(yearPlus2 = year + 2) %>%
  dplyr::select("yearPlus2", "z_prec_secIt_T2" = "z_sumWinPrec")

# add T2 temperatures
di_win_secIt <- di_win_secIt %>%
  dplyr::left_join(PrecT2_secIt, by = c("WinterYear" = "yearPlus2"))


# Model precipitation in T0 - Second iteration ----------------------------
di_win <- di_win_secIt %>%
  dplyr::left_join(z_all_windows,
                   by = c("WinterYear" = "year"), relationship = "many-to-many")

## model precipitation in T0 ----
model_PrecT0_secIt <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_secIt_T0 + z_temp_meanMax_secIt_T1 + z_temp_meanMax_secIt_T2 +
                               z_sumWinPrec + z_prec_secIt_T1 +  z_prec_secIt_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run model ####
prepare_and_run_model(model = model_PrecT0_secIt,
                      variable_name = "z_PrecT0_secIt")

## Look at model output ####
summary_PrecT0_secIt <- get(load(here::here("results", "z_PrecT0_secIt/summary_z_PrecT0_secIt.rda")))

# # plot results
plots <- plot_results(dataframe = summary_meanPrecT0_secIt,
             plot_title = "z_precT0_secIt")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_PrecT0_secIt <- summary_PrecT0_secIt %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate precipitation in T0 for selected window
PrecT0_secIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_PrecT0_secIt$Start_DOY, End_DOY == best_wind_PrecT0_secIt$End_DOY) %>%
  dplyr::select("year", "z_prec_secIt_T0" = "z_sumWinPrec")

# add T0 precipitation
di_win_secIt <- di_win_secIt %>%
  dplyr::left_join(PrecT0_secIt, by = c("WinterYear" = "year"))


# III. Third iteration ----------------------------------------------------

# delete old window data from input data frame for second iteration to save size
di_win_thirdIt <- di_win_secIt %>%
  dplyr::select(!c("z_temp_meanMean_firstIt_T0", "z_temp_meanMax_firstIt_T1", "z_temp_meanMax_firstIt_T2", 
                   "z_prec_firstIt_T0", "z_prec_firstIt_T1", "z_prec_firstIt_T2"))

# Model temperature in T1 - Third iteration -----------------------------------
di_win <- di_win_thirdIt %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T1 = year + 1),
                   by = c("WinterYear" = "year_T1"), relationship = "many-to-many")


## Model function ----
model_maxTempT1_thirdIt <- function(datfr) {
  
  
  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_secIt_T0 + z_meanMaxWinTemp + z_temp_meanMax_secIt_T2 +
                               z_prec_secIt_T0 + z_prec_secIt_T1 + z_prec_secIt_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

# ## Run model ####
prepare_and_run_model(model = model_maxTempT1_thirdIt,
                      variable_name = "z_maxTempT1_thirdIt")

## Look at model output ####
summary_maxTempT1_thirdIt <- get(load(here::here("results", "z_maxTempT1_thirdIt/summary_z_maxTempT1_thirdIt.rda")))

# plot results
plots <- plot_results(dataframe = summary_maxTempT1_thirdIt,
                      plot_title = "z_maxTempT1_thirdIt")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_maxTempT1_thirdIt <- summary_maxTempT1_thirdIt %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate temperature for maxTemp in T1 for selected window
maxTempT1_thirdIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_maxTempT1_thirdIt$Start_DOY, End_DOY == best_wind_maxTempT1_thirdIt$End_DOY) %>%
  dplyr::mutate(yearPlus1 = year + 1) %>%
  dplyr::select("yearPlus1", "z_temp_meanMax_thirdIt_T1" = "z_meanMaxWinTemp")

# add T1 temperatures
di_win_thirdIt <- di_win_thirdIt %>%
  dplyr::left_join(maxTempT1_thirdIt, by = c("WinterYear" = "yearPlus1"))


# Model temperature in T2 - Third iteration ------------------------------
di_win <- di_win_thirdIt %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T2 = year + 2),
                   by = c("WinterYear" = "year_T2"), relationship = "many-to-many")

## Model function ####
model_maxTempT2_thirdIt <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_secIt_T0 + z_temp_meanMax_thirdIt_T1 + z_meanMaxWinTemp +
                               z_prec_secIt_T0 + z_prec_secIt_T1 + z_prec_secIt_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run model ####
prepare_and_run_model(model = model_maxTempT2_thirdIt,
                      variable_name = "z_maxTempT2_thirdIt")

## Look at model output ####
summary_maxTempT2_thirdIt <- get(load(here::here("results", "z_maxTempT2_thirdIt/summary_z_maxTempT2_thirdIt.rda")))

# plot results
plots <- plot_results(dataframe = summary_maxTempT2_thirdIt,
             plot_title = "z_maxTempT2_thirdIt")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_maxTemp2_thirdIt <- summary_maxTempT2_thirdIt %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate temperature for maxTemp in T2 for selected window
maxTempT2_thirdIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_maxTemp2_thirdIt$Start_DOY, End_DOY == best_wind_maxTemp2_thirdIt$End_DOY) %>%
  dplyr::mutate(yearPlus2 = year + 2) %>%
  dplyr::select("yearPlus2", "z_temp_meanMax_thirdIt_T2" = "z_meanMaxWinTemp")

# add T2 temperatures
di_win_thirdIt <- di_win_thirdIt %>%
  dplyr::left_join(maxTempT2_thirdIt, by = c("WinterYear" = "yearPlus2"))


# Model temperature in T0 - Third iteration ------------------------------
di_win <- di_win_thirdIt %>%
  dplyr::left_join(z_all_windows,
                   by = c("WinterYear" = "year"), relationship = "many-to-many")

## Model function ----
model_meanTempT0_thirdIt <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_meanMeanWinTemp + z_temp_meanMax_thirdIt_T1 + z_temp_meanMax_thirdIt_T2 +
                               z_prec_secIt_T0 + z_prec_secIt_T1 + z_prec_secIt_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run model ####
prepare_and_run_model(model = model_meanTempT0_thirdIt,
                      variable_name = "z_meanTempT0_thirdIt")

## Look at model putput ####
summary_meanTempT0_thirdIt <- get(load(here::here("results", "z_meanTempT0_thirdIt/summary_z_meanTempT0_thirdIt.rda")))

# plot results
plots <- plot_results(dataframe = summary_meanTempT0_thirdIt,
             plot_title = "z_meanTempT0_thirdIt")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_meanTempT0_thirdIt <- summary_meanTempT0_thirdIt %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate temperature for meanTemp in T0 for selected window
meanTempT0_thirdIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_meanTempT0_thirdIt$Start_DOY, End_DOY == best_wind_meanTempT0_thirdIt$End_DOY) %>%
  dplyr::select("year", "z_temp_meanMean_thirdIt_T0" = "z_meanMeanWinTemp")

# add T0 temperatures
di_win_thirdIt <- di_win_thirdIt %>%
  dplyr::left_join(meanTempT0_thirdIt, by = c("WinterYear" = "year"))


# Model precipitation in T1 - Third iteration ----------------------------
di_win <- di_win_thirdIt %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T1 = year + 1),
                   by = c("WinterYear" = "year_T1"), relationship = "many-to-many")

## Model function -----
model_PrecT1_thirdIt <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_thirdIt_T0 + z_temp_meanMax_thirdIt_T1 + z_temp_meanMax_thirdIt_T2 +
                               z_prec_secIt_T0 + z_sumWinPrec  + z_prec_secIt_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

# Run model ####
prepare_and_run_model(model = model_PrecT1_thirdIt,
                      variable_name = "z_PrecT1_thirdIt")

## Look at model output ####
summary_PrecT1_thirdIt <- get(load(here::here("results", "z_PrecT1_thirdIt/summary_z_PrecT1_thirdIt.rda")))

# plot results
plots <- plot_results(dataframe = summary_PrecT1_thirdIt,
             plot_title = "z_precT1_thirdIt")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_PrecT1_thirdIt <- summary_PrecT1_thirdIt %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate temperature for precipitation in T1 for selected window
PrecT1_thirdIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_PrecT1_thirdIt$Start_DOY, End_DOY == best_wind_PrecT1_thirdIt$End_DOY) %>%
  dplyr::mutate(yearPlus1 = year + 1) %>%
  dplyr::select("yearPlus1", "z_prec_thirdIt_T1" = "z_sumWinPrec")

# add T1 temperatures
di_win_thirdIt <- di_win_thirdIt %>%
  dplyr::left_join(PrecT1_thirdIt, by = c("WinterYear" = "yearPlus1"))


# Model precipitation in T2 - Third iteration ----------------------------
di_win <- di_win_thirdIt %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T2 = year + 2),
                   by = c("WinterYear" = "year_T2"), relationship = "many-to-many")

## Model function ----
model_PrecT2_thirdIt <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_thirdIt_T0 + z_temp_meanMax_thirdIt_T1 + z_temp_meanMax_thirdIt_T2 +
                               z_prec_secIt_T0  + z_prec_thirdIt_T1 +  z_sumWinPrec +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run model ####
prepare_and_run_model(model = model_PrecT2_thirdIt,
                      variable_name = "z_PrecT2_thirdIt")

## Look at model output ####
summary_PrecT2_thirdIt <- get(load(here::here("results", "z_PrecT2_thirdIt/summary_z_PrecT2_thirdIt.rda")))

# plot results
plots <- plot_results(dataframe = summary_PrecT2_thirdIt,
             plot_title = "z_precT2_thirdIt")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_PrecT2_thirdIt <- summary_PrecT2_thirdIt %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate precipitation in T2 for selected window
PrecT2_thirdIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_PrecT2_thirdIt$Start_DOY, End_DOY == best_wind_PrecT2_thirdIt$End_DOY) %>%
  dplyr::mutate(yearPlus2 = year + 2) %>%
  dplyr::select("yearPlus2", "z_prec_thirdIt_T2" = "z_sumWinPrec")

# add T2 temperatures
di_win_thirdIt <- di_win_thirdIt %>%
  dplyr::left_join(PrecT2_thirdIt, by = c("WinterYear" = "yearPlus2"))


# Model precipitation in T0 - Third iteration ----------------------------
di_win <- di_win_thirdIt %>%
  dplyr::left_join(z_all_windows,
                   by = c("WinterYear" = "year"), relationship = "many-to-many")

## Model function----
model_PrecT0_thirdIt <- function(datfr) {


  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_thirdIt_T0 + z_temp_meanMax_thirdIt_T1 + z_temp_meanMax_thirdIt_T2 +
                               z_sumWinPrec + z_prec_thirdIt_T1 +  z_prec_thirdIt_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run model ####
prepare_and_run_model(model = model_PrecT0_thirdIt,
                      variable_name = "z_PrecT0_thirdIt")

## Look at model output ####
summary_PrecT0_thirdIt <- get(load(here::here("results", "z_PrecT0_thirdIt/summary_z_PrecT0_thirdIt.rda")))

# plot results
plots <- plot_results(dataframe = summary_meanPrecT0_thirdIt,
             plot_title = "z_precT0_thirdIt")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_PrecT0_thirdIt <- summary_PrecT0_thirdIt %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate precipitation in T0 for selected window
PrecT0_thirdIt <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_PrecT0_thirdIt$Start_DOY, End_DOY == best_wind_PrecT0_thirdIt$End_DOY) %>%
  dplyr::select("year", "z_prec_thirdIt_T0" = "z_sumWinPrec")

# add T0 precipitation
di_win_thirdIt <- di_win_thirdIt %>%
  dplyr::left_join(PrecT0_thirdIt, by = c("WinterYear" = "year"))



# Fourth iteration --------------------------------------------------------

# delete old window data from input data frame for second iteration to save size
di_win_4It <- di_win_thirdIt %>%
  dplyr::select(!c("z_temp_meanMean_secIt_T0", "z_temp_meanMax_secIt_T1", "z_temp_meanMax_secIt_T2", 
                   "z_prec_secIt_T0", "z_prec_secIt_T1", "z_prec_secIt_T2"))

# Model temperature in T1 - Third iteration -----------------------------------
di_win <- di_win_4It %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T1 = year + 1),
                   by = c("WinterYear" = "year_T1"), relationship = "many-to-many")


## Model function ----
model_maxTempT1_4It <- function(datfr) {
  
  
  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_thirdIt_T0 + z_meanMaxWinTemp + z_temp_meanMax_thirdIt_T2 +
                               z_prec_thirdIt_T0 + z_prec_thirdIt_T1 + z_prec_thirdIt_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

# ## Run model ####
prepare_and_run_model(model = model_maxTempT1_4It,
                      variable_name = "z_maxTempT1_4It")

## Look at model output ####
summary_maxTempT1_4It <- get(load(here::here("results", "z_maxTempT1_4It/summary_z_maxTempT1_4It.rda")))

# plot results
plots <- plot_results(dataframe = summary_maxTempT1_4It,
                      plot_title = "z_maxTempT1_4It")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_maxTempT1_4It <- summary_maxTempT1_4It %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate temperature for maxTemp in T1 for selected window
maxTempT1_4It <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_maxTempT1_4It$Start_DOY, End_DOY == best_wind_maxTempT1_4It$End_DOY) %>%
  dplyr::mutate(yearPlus1 = year + 1) %>%
  dplyr::select("yearPlus1", "z_temp_meanMax_4It_T1" = "z_meanMaxWinTemp")

# add T1 temperatures
di_win_4It <- di_win_4It %>%
  dplyr::left_join(maxTempT1_4It, by = c("WinterYear" = "yearPlus1"))


# Model temperature in T2 - Fourth iteration ------------------------------
di_win <- di_win_4It %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T2 = year + 2),
                   by = c("WinterYear" = "year_T2"), relationship = "many-to-many")

## Model function ####
model_maxTempT2_4It <- function(datfr) {
  
  
  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_thirdIt_T0 + z_temp_meanMax_4It_T1 + z_meanMaxWinTemp +
                               z_prec_thirdIt_T0 + z_prec_thirdIt_T1 + z_prec_thirdIt_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run model ####
prepare_and_run_model(model = model_maxTempT2_4It,
                      variable_name = "z_maxTempT2_4It")

## Look at model output ####
summary_maxTempT2_4It <- get(load(here::here("results", "z_maxTempT2_4It/summary_z_maxTempT2_4It.rda")))

# plot results
plots <- plot_results(dataframe = summary_maxTempT2_4It,
                      plot_title = "z_maxTempT2_4It")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_maxTemp2_4It <- summary_maxTempT2_4It %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate temperature for maxTemp in T2 for selected window
maxTempT2_4It <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_maxTemp2_4It$Start_DOY, End_DOY == best_wind_maxTemp2_4It$End_DOY) %>%
  dplyr::mutate(yearPlus2 = year + 2) %>%
  dplyr::select("yearPlus2", "z_temp_meanMax_4It_T2" = "z_meanMaxWinTemp")

# add T2 temperatures
di_win_4It <- di_win_4It %>%
  dplyr::left_join(maxTempT2_4It, by = c("WinterYear" = "yearPlus2"))


# Model temperature in T0 - Fourth iteration ------------------------------
di_win <- di_win_4It %>%
  dplyr::left_join(z_all_windows,
                   by = c("WinterYear" = "year"), relationship = "many-to-many")

## Model function ----
model_meanTempT0_4It <- function(datfr) {
  
  
  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_meanMeanWinTemp + z_temp_meanMax_4It_T1 + z_temp_meanMax_4It_T2 +
                               z_prec_thirdIt_T0 + z_prec_thirdIt_T1 + z_prec_thirdIt_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run model ####
prepare_and_run_model(model = model_meanTempT0_4It,
                      variable_name = "z_meanTempT0_4It")

## Look at model output ####
summary_meanTempT0_4It <- get(load(here::here("results", "z_meanTempT0_4It/summary_z_meanTempT0_4It.rda")))

# plot results
plots <- plot_results(dataframe = summary_meanTempT0_4It,
                      plot_title = "z_meanTempT0_4It")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_meanTempT0_4It <- summary_meanTempT0_4It %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate temperature for meanTemp in T0 for selected window
meanTempT0_4It <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_meanTempT0_4It$Start_DOY, End_DOY == best_wind_meanTempT0_4It$End_DOY) %>%
  dplyr::select("year", "z_temp_meanMean_4It_T0" = "z_meanMeanWinTemp")

# add T0 temperatures
di_win_4It <- di_win_4It %>%
  dplyr::left_join(meanTempT0_4It, by = c("WinterYear" = "year"))


# Model precipitation in T1 - Fourth iteration ----------------------------
di_win <- di_win_4It %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T1 = year + 1),
                   by = c("WinterYear" = "year_T1"), relationship = "many-to-many")

## Model function -----
model_PrecT1_4It <- function(datfr) {
  
  
  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_4It_T0 + z_temp_meanMax_4It_T1 + z_temp_meanMax_4It_T2 +
                               z_prec_thirdIt_T0 + z_sumWinPrec + z_prec_thirdIt_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

# Run model ####
prepare_and_run_model(model = model_PrecT1_4It,
                      variable_name = "z_PrecT1_4It")

## Look at model output ####
summary_PrecT1_4It <- get(load(here::here("results", "z_PrecT1_4It/summary_z_PrecT1_4It.rda")))

# plot results
plots <- plot_results(dataframe = summary_PrecT1_4It,
                      plot_title = "z_precT1_4It")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_PrecT1_4It <- summary_PrecT1_4It %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate temperature for precipitation in T1 for selected window
PrecT1_4It <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_PrecT1_4It$Start_DOY, End_DOY == best_wind_PrecT1_4It$End_DOY) %>%
  dplyr::mutate(yearPlus1 = year + 1) %>%
  dplyr::select("yearPlus1", "z_prec_4It_T1" = "z_sumWinPrec")

# add T1 temperatures
di_win_4It <- di_win_4It %>%
  dplyr::left_join(PrecT1_4It, by = c("WinterYear" = "yearPlus1"))


# Model precipitation in T2 - Fourth iteration ----------------------------
di_win <- di_win_4It %>%
  dplyr::left_join(z_all_windows %>%
                     dplyr::mutate(year_T2 = year + 2),
                   by = c("WinterYear" = "year_T2"), relationship = "many-to-many")

## Model function ----
model_PrecT2_4It <- function(datfr) {
  
  
  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_4It_T0 + z_temp_meanMax_4It_T1 + z_temp_meanMax_4It_T2 +
                               z_prec_thirdIt_T0  + z_prec_4It_T1 +  z_sumWinPrec +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run model ####
prepare_and_run_model(model = model_PrecT2_4It,
                      variable_name = "z_PrecT2_4It")

## Look at model output ####
summary_PrecT2_4It <- get(load(here::here("results", "z_PrecT2_4It/summary_z_PrecT2_4It.rda")))

# plot results
plots <- plot_results(dataframe = summary_PrecT2_4It,
                      plot_title = "z_precT2_4It")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_PrecT2_4It <- summary_PrecT2_4It %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate precipitation in T2 for selected window
PrecT2_4It <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_PrecT2_4It$Start_DOY, End_DOY == best_wind_PrecT2_4It$End_DOY) %>%
  dplyr::mutate(yearPlus2 = year + 2) %>%
  dplyr::select("yearPlus2", "z_prec_4It_T2" = "z_sumWinPrec")

# add T2 temperatures
di_win_4It <- di_win_4It %>%
  dplyr::left_join(PrecT2_4It, by = c("WinterYear" = "yearPlus2"))


# Model precipitation in T0 - Fourth iteration ----------------------------
di_win <- di_win_4It %>%
  dplyr::left_join(z_all_windows,
                   by = c("WinterYear" = "year"), relationship = "many-to-many")

## Model function----
model_PrecT0_4It <- function(datfr) {
  
  
  # model
  output <- glmmTMB::glmmTMB(TotalNuts ~
                               z_temp_meanMean_4It_T0 + z_temp_meanMax_4It_T1 + z_temp_meanMax_4It_T2 +
                               z_sumWinPrec + z_prec_4It_T1 +  z_prec_4It_T2 +
                               TotalNuts_T1 + (1|TreeID),
                             data = datfr,
                             family = nbinom2(link = "log"),
                             ziformula = ~ .)
}

## Run model ####
prepare_and_run_model(model = model_PrecT0_4It,
                      variable_name = "z_PrecT0_4It")

## Look at model output ####
summary_PrecT0_4It <- get(load(here::here("results", "z_PrecT0_4It/summary_z_PrecT0_4It.rda")))

# plot results
plots <- plot_results(dataframe = summary_PrecT0_4It,
                      plot_title = "z_precT0_4It")

# heat map
plots[1]

# AIC peak
plots[2]

# filter for lowest AIC window
best_wind_PrecT0_4It <- summary_PrecT0_4It %>%
  dplyr::filter(AIC == min(AIC, na.rm = TRUE))

# calculate precipitation in T0 for selected window
PrecT0_4It <- z_all_windows %>%
  dplyr::filter(Start_DOY == best_wind_PrecT0_4It$Start_DOY, End_DOY == best_wind_PrecT0_4It$End_DOY) %>%
  dplyr::select("year", "z_prec_4It_T0" = "z_sumWinPrec")

# add T0 precipitation
di_win_4It <- di_win_4It %>%
  dplyr::left_join(PrecT0_4It, by = c("WinterYear" = "year"))




# Build final data file and export main results--------------------------------

# help file to re-assign WinterYear back to ordinal years
ord_year <- tibble::tibble(WinterYear = c(1976:2025),
                           ord_year = c(1:50))

# show names of final data frame
names(di_win_4It)

# delete climate variables from previous iterations
di_final <- di_win_4It %>%
  dplyr::select(!c("z_temp_meanMax_thirdIt_T1", "z_temp_meanMax_thirdIt_T2", "z_temp_meanMean_thirdIt_T0", "z_prec_thirdIt_T1",
                   "z_prec_thirdIt_T2", "z_prec_thirdIt_T0", "z_meanTempSummer_T1", "z_meanTempSummer_T2", "z_minTempSummer_T1",
                   "z_minTempSummer_T2")) %>%
  dplyr::left_join(ord_year, by = "WinterYear")

# write the final data frame to file
write.csv(di_final, file = here::here("data", "di_after4thIt.csv"), row.names = FALSE)


## Prepare files for processing (requires to have all summary and best window objects loaded in the environment)

## best window objects ##

# get names of all best_wind objects in the environment
names_best_wind <- ls()[grepl("^best_wind_", ls())]

# collect all objects in a list
best_wind_list <- mget(names_best_wind)

# convert list into data frame
sens_windows <- do.call(rbind, best_wind_list) %>%
  dplyr::arrange(desc(AIC)) %>%
  dplyr::mutate("Window_start" = paste0(Start_DOY, " (",
                                        substring(as.Date(Start_DOY, origin = "2025-01-01"), first = 9), "/",
                                        substring(as.Date(Start_DOY, origin = "2025-01-01"), first = 6, last = 7), ")"),
                "Window_end" = paste0(End_DOY, " (",
                                      substring(as.Date(End_DOY, origin = "2025-01-01"), first = 9), "/",
                                      substring(as.Date(End_DOY, origin = "2025-01-01"), first = 6, last = 7), ")")) %>%
  dplyr::select(!c("windowID", "deltaAIC", "batch_file_name")) %>%
  dplyr::rename("Window_length" = "DOY_duration")

write.csv(sens_windows, file = here::here("data", "table_sens_windows.csv"), row.names = FALSE)


## summary objects ##

# get names of all summary files
names_summary <- ls()[grepl("^summary_", ls())]

# collect all objects in a list
summary_list <- mget(names_summary)

summary_all <- do.call(rbind, summary_list)

write.csv(summary_all, file = here::here("data", "summary_all_models.csv"), row.names = FALSE)