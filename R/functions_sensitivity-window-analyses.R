## Functions for sensitivity window analyses ####

# Author:           Cherine C. Jantzen
# Created:          2025-05-20
# Last updated:     2025-08-15

################################
## Content: This script contains three functions that are used to automate the process of finding sensitivity windows. These functions partly 
## run processes in parallel and are tailored to the server used here and the analysis at hand, rather than meant to be used outside of this analysis. They are simply intended
## to reduce redundancy in the script, as similar code is run several times.
################################

# load libraries
library(future)
library(future.apply)
library(future.callr)
library(dplyr)
library(tidyr)
library(parallelly)
library(purrr)
library(here)
library(ggplot2)
library(data.table)
library(stringr)

# Function I: Parallel modelling of sensitivity windows  ------------------

## Description ----
### This function uses the future package to calculate models in parallel. For an input list containing a nested data frame for each window in each list item and a defined number of batches and their size, it subsets the list into these batches, runs a specified model on each list item and saves the models of each batch to file. Afterwards it re-reads all files, gives the summary of each model and saves the AIC and p-values for each window in a summary output, which is also written to file. Batching and writing to file is done because the input data can be very large and secures that the model outputs are partially saved in case the function run gets disrupted. 

## Arguments ----
### clim_var: string, specifying the name of the climate variable for which the windows are tested and that is assigned to the output data
### model_function: function that contains only the model which should be used 
### input_list: list with nested data tables (or data frames). Each list item should be named with the window name and contain a data frame with all model varaiables for this window
### batch_size: integer, specifying the size of each batch
### no_batches: integer, specifying the number of batches that the data should be split in


sens_windows <- function(clim_var,
                         model_function,
                         input_list,
                         batch_size,
                         no_batches){
  
  future::plan(future.callr::callr(workers = 10))
  
  options(future.globals.maxSize = +Inf)
  
  if (!dir.exists(here::here(paste0("results/", clim_var)))) {dir.create(here::here(paste0("results/", clim_var)))}
  
  # loop model over all batches
  purrr::walk(.x = c(0 : (no_batches - 1)),
              .f = ~{
                
                ls_subset <- input_list[(.x * batch_size + 1) : ((.x + 1) * batch_size)]
                
                if (.x == (no_batches - 1)) {
                  
                  ls_subset <- input_list[(.x * batch_size + 1) : length(input_list)]
                  
                }
                
                ls_models <- future.apply::future_lapply(ls_subset, purrr::possibly(model_function, otherwise = NA),  future.stdout = NA)
                
                save(ls_models, file =  here::here(paste0("results/", clim_var, "/", clim_var, "-batch", (.x + 1), ".rda")))   
                
              }, .progress = TRUE) 
  
  # End the parallel session
  future::plan(sequential)
  
  # get file names of all batches to load and process them
  file_names <- list.files(
    path = here::here(paste0("results/", clim_var)),
    pattern = paste0(clim_var, "-batch", ".*\\.rda$"),  
    full.names = TRUE)
  
  
  summary_df <- purrr::map(.x = file_names,
                           .f = ~{
                             
                             # loop through all batches and load each R object (model output of one window)
                             load(.x)
                             
                             # save batch name to make models easier findable later
                             batch_name <- stringr::str_extract(.x, pattern = "(?<=-).*")
                             
                             purrr::map(.x = 1:length(ls_models),
                                        .f = ~{
                                          
                                          summary <- summary(ls_models[[.x]])
                                          
                                          if (is.list(summary) == FALSE) {
                                            
                                            AIC_value <- NA
                                            
                                          } else {
                                            
                                            AIC_value <- summary$AICtab["AIC"]
                                            
                                          }
                                          
                                          output <- tibble::tibble(variable = clim_var, 
                                                                   AIC = AIC_value,
                                                                   windowID = names(ls_models[.x]),
                                                                   Start_DOY = as.numeric(stringr::str_extract(string = windowID, pattern = "[^_]+")),
                                                                   End_DOY = as.numeric(stringr::str_extract(string = windowID, pattern = "(?<=_).*")),
                                                                   batch_file_name = batch_name)
                                          
                                          
                                          return(output)
                                          
                                        }) 
                             
                           }, .progress = TRUE) %>% bind_rows()
  
  summary_df <- summary_df %>% 
    dplyr::mutate(DOY_duration = (End_DOY - Start_DOY) + 1,
                  deltaAIC = AIC - min(AIC, na.rm = TRUE))
  
  save(summary_df, 
       file = here::here(paste0("results/", clim_var, "/summary_", clim_var, ".rda")))
  
}



# Function II: Prepare model input files and run models  -----------

## Description: ----
### This function prepares the inout file for the model function, and prepares the remaining 
## input parameters. These steps are the same for all climate variables, while the model function has to 
## be specified individually outside of this function. The model function is called and the models run.


## Arguments: ----
### model: model function for the climate variable in focus (has to be specified beforehand as it varies between variables)
### variable_name: string, specifying the name of the focal cliamte variale that is input for "sens_windows" function and requivalent to its clim_var argument

prepare_and_run_model <- function(model,
                                  variable_name){
  
  # create nested data frames per window for easier processing
  by_window <- di_win %>%
    dplyr::group_by(windowID) %>%
    tidyr::nest()
  
  # create input list for models
  ls_by_window <- by_window %>%
    dplyr::rowwise() %>%
    dplyr::pull(data, name = windowID)
  
  # convert data frames to data tables to speed up computation time
  ls_by_window_dt <- purrr::map(.x = names(ls_by_window),
                                .f = ~{
                                  
                                  ls_by_window %>%
                                    purrr::pluck(.x) %>%
                                    data.table::setDT()
                                  
                                }
  ) %>% setNames(names(ls_by_window))
  
  # remove unneeded large objects from environment
  rm(ls_by_window, by_window)
  
  # set batch parameters
  batch_size <- floor(length(ls_by_window_dt) / no_batches)
  
  ## Run model ####
  sens_windows(clim_var = variable_name,
               model_function = model,
               batch_size = batch_size,
               no_batches = no_batches,
               input_list = ls_by_window_dt)
}


# Function III: Visualise results ----------------------------------------

## Description: ----
### This function creates two plots oout of the summary of the model results. First, a heat map of the AIC values and
### secon, a plot showing AIC values spread over window opening days to visualize where peaks are. Both plots are stored 
### in a list and saved to the plots folder. 

## Arguments: ----
### dataframe:  a dataframe, the summary output of the "sens_window" function
### plot title: string, specifying the climate variable name as the title of the plots and of the saved plot objects


plot_results <- function(dataframe,
                         plot_title){
  
  p1 <- ggplot2::ggplot(dataframe,
                        ggplot2::aes(x = Start_DOY, y = End_DOY, fill = AIC)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = "viridis") +
    ggplot2::labs(x = "Start DOY of window",
                  y = "End DOY of window", 
                  title = plot_title)
  
  ggplot2::ggsave(plot = ggplot2::last_plot(), filename = here::here("plots", paste0("heat-map_", plot_title, ".png")), width = 30, height = 20, unit = "cm")
  
  
  p2 <- ggplot2::ggplot(dataframe, ggplot2::aes(x = Start_DOY, y = AIC, colour = AIC)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "viridis") +
    ggplot2::labs(x = "Start DOY of window", 
                  title = plot_title)
  
  ggplot2::ggsave(plot = ggplot2::last_plot(), filename = here::here("plots", paste0(plot_title, ".png")), width = 30, height = 20, unit = "cm")
  
  return(list(p1, p2))
}

