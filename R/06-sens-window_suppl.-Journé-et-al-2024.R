# Analysis supplementary material - Sensitivity windows based on appraoch used in Journé et al., 2024 ####

# Author:       Cherine C. Jantzen
# Created:      2025-01-28
# Last updated: 2025-12-19


################################
## Content: Here we reproduce the analysis of Journé et al., (2024) Nat Pla (10.1038/s41477-024-01651-w) to compare the results
## of their method to find the best window of sensitivity with the results of the approach we used. This analysis is part of the supplementary 
## material and based on the textual description of the methods in Journé et al.
################################

# load packages
library(dplyr)
library(purrr)
library(here)
library(ggplot2)
library(glmmTMB)
library(tidyr)

# load data
di <- read.csv(here::here("data", "di_input-data.csv"))
df_climate <- read.csv(here::here("data", "df_climate-processed-knmi-data.csv"))

# create directory for plots if not existing
if (!dir.exists(here::here("plots"))) {dir.create(here::here("plots"))}


# I. Slice data in windows and calculate climate variables per window --------

# set window length 
window_length <- 7

# calculate the mean daily temperature of every possible window of the given window length of the year
temp_windows <- purrr::map(.x = (1:((366 + 1) - window_length)), 
                           .f = ~{
                             
                             df_climate %>% 
                               dplyr::group_by(year) %>% 
                               # slice data in respective window
                               dplyr::slice(.x:(.x + (window_length - 1))) %>% 
                               dplyr::summarise(meanMeanWinTemp = mean(meanDayTemp, na.rm = TRUE), 
                                                meanMaxWinTemp = mean(maxDayTemp, na.rm = TRUE),
                                                meanMinWinTemp = mean(minDayTemp, na.rm = TRUE),
                                                sumWinPrec = sum(sumDayPrec, na.rm = TRUE)) %>% 
                               dplyr::mutate(DOY_windowOpen = paste0("DOY_", .x),
                                             DOY_windowClosed = paste0("DOY_", (.x + (window_length - 1))))
                             
                           },
                           
                           .progress = TRUE
) %>%  dplyr::bind_rows()


# II. Find sensitivity window with Spearman correlation -----------------------

# get population level nuts and ln transform the annual seed output
pop_nuts <- di %>% 
  dplyr::summarise(sum_TotalNuts = sum(TotalNuts, na.rm = TRUE),
                   n_trees = dplyr::n(),
                   .by = WinterYear) %>% 
  dplyr::mutate(pop_TotalNuts = sum_TotalNuts/n_trees,
                ln_pop_TotalNuts = log(pop_TotalNuts)) %>% 
  dplyr::select("WinterYear", "ln_pop_TotalNuts")

get_sensWind_spearman <- function(rel_year = c("T0", "T1", "T2"),
                                  climate_variable) {
  
  
  # pivot temperatures per window to wider format  
  temp_win_pivot <- temp_windows %>% 
    dplyr::select("variable" = all_of(climate_variable), "DOY_windowClosed", "year") %>% 
    tidyr::pivot_wider(names_from = DOY_windowClosed, values_from = variable) %>% 
    mutate(year_T1 = year + 1,
           year_T2 = year + 2) 
  
  
  if (rel_year == "T1") {
    
    # Calculate correlations for T1
    cor_matrix <- pop_nuts %>% 
      dplyr::left_join(temp_win_pivot %>% 
                         dplyr::select(!c("year", "year_T2")), 
                       by = c("WinterYear" = "year_T1")) %>% 
      dplyr::select(!"WinterYear") %>% 
      cor(method = "spearman", use = "pairwise")
    
  } else if (rel_year == "T2") {
    
    # Calculate correlations for T2 
    cor_matrix <- pop_nuts %>% 
      dplyr::left_join(temp_win_pivot %>% 
                         dplyr::select(!c("year", "year_T1")), 
                       by = c("WinterYear" = "year_T2")) %>% 
      dplyr::select(!"WinterYear") %>% 
      cor(method = "spearman", use = "pairwise")
    
  } else {
    
    # Calculate correlations for T0
    cor_matrix <- pop_nuts %>% 
      dplyr::left_join(temp_win_pivot %>% 
                         dplyr::select(!c("year_T2", "year_T1")), 
                       by = c("WinterYear" = "year")) %>% 
      dplyr::select(!"WinterYear") %>% 
      cor(method = "spearman", use = "pairwise")
    
  }
  
  # set diagonal of matrix to NA (i.e., self-correlation)
  diag(cor_matrix) <- NA
  
  # convert matrix to data frame with numeric DOY when window opens and correlation coefficient
  cor_df <- as.data.frame(cor_matrix) %>% 
    tibble::rownames_to_column(var = "DOY_closed") %>% 
    dplyr::filter(DOY_closed != "ln_pop_TotalNuts") %>% 
    dplyr::mutate(DOY_closed = as.numeric(stringr::str_remove(DOY_closed, pattern = "DOY_"))) %>% 
    dplyr::select("DOY_closed", "Rho" = "ln_pop_TotalNuts") %>% 
    dplyr::mutate(seRho = sqrt((1 - Rho^2)^2 * ((1 + (Rho^2 / 2)) / (nrow(pop_nuts) - 3))))
  
  # find the end day of the window with the highest correlation & convert to calendar date
  end_date <- cor_df %>% 
    dplyr::mutate(abs_Rho = abs(Rho)) %>% 
    dplyr::filter(abs_Rho == max(abs_Rho)) %>% 
    dplyr::pull(DOY_closed) %>% 
    lubridate::as_date(origin = "2013-01-01") # chose random year that's not a leap year
  
  return(list(cor_df, end_date))
  
}


# Run function for all possible combinations of relative year, climate variable and period

## expand grid to get all possible combinations
combinations <- expand.grid(climate_var = c("meanMeanWinTemp", "sumWinPrec", "meanMaxWinTemp"), 
                            rel_year = c("T1", "T2", "T0"))

## run function over all combinations
spearman_windows <- purrr::map2(.x = combinations$climate_var,
                                .y = combinations$rel_year,
                                .f = ~{
                                  
                                  
                                  correlations <- get_sensWind_spearman(rel_year = .y, 
                                                                        climate_variable = .x) 
                                  
                                  cor_df <- correlations[[1]] %>% 
                                    dplyr::mutate(rel_year = .y, 
                                                  climate_variable = .x)
                                  
                                  max_cor <- tibble::tibble(rel_year = .y, 
                                                            climate_variable = .x,
                                                            end_dateWin = correlations[[2]])
                                  
                                  
                                  
                                  output_df <- dplyr::bind_rows(cor_df, max_cor)
                                  
                                  return(output_df)
                                  
                                }, .progress = "spearman correlations"
) %>% dplyr::bind_rows()


## Figure S3: Spearman correlations for mean and max daily temperature and seed production ####
spearman_windows %>% 
  dplyr::filter(is.na(end_dateWin), climate_variable != "sumWinPrec") %>% 
  dplyr::mutate(newDOY_closed = dplyr::case_when(rel_year == "T0" ~ DOY_closed + 723,
                                                 rel_year == "T1" ~ DOY_closed + 358,
                                                 TRUE ~ DOY_closed)) %>% 
  ggplot2::ggplot(ggplot2::aes(x = newDOY_closed, y = Rho)) + 
  ggplot2::geom_point(size = 2, ggplot2::aes(colour = Rho > 0), alpha = 0.6) + 
  ggplot2::geom_vline(xintercept = 172, linetype = "dashed") + 
  ggplot2::geom_vline(xintercept = 537, linetype = "dashed") + 
  ggplot2::geom_vline(xintercept = 902, linetype = "dashed") +
  # ggplot2::geom_line(ggplot2::aes(y = Rho + seRho), alpha = 0.3) +
  # ggplot2::geom_line(ggplot2::aes(y = Rho - seRho), alpha = 0.3) +
  ggplot2::geom_hline(yintercept = 0) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::labs(x = "Day the window closes", y = "Spearman correlation (Rho)") +
  ggplot2::scale_colour_manual(values = c("firebrick", "dodgerblue")) +
  ggplot2::scale_x_continuous(breaks = c(0, 172, 365, 537, 730, 902)) +
  ggplot2::theme_minimal(base_size = 17) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::lims(y = c(-0.6, 0.6)) + 
  ggplot2::facet_wrap(~ climate_variable, nrow = 2, ncol = 1, labeller = labeller(climate_variable = c("meanMeanWinTemp" = "mean daily temperature", "meanMaxWinTemp" = "maximum daily temperature")))

ggplot2::ggsave(plot = ggplot2::last_plot(), 
                file = here::here("plots", "Figure_S3.png"), 
       units = "cm", width = 30, height = 25, dpi = 600)
