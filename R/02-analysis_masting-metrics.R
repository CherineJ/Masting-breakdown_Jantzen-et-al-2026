## Analysis beech crop data - Masting metrics ####

# Author:       Cherine C. Jantzen
# Created:      2025-01-13
# Last updated: 2026-02-02

################################
## Content: In this script, we analyse the classical masting metrics, i.e., synchrony and inter-annual variation in seed production, 
## as well as the temporal changes in beechnut production.
################################

# I. Preparation -------------------------------------------------------

# load packages
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(mgcv)
library(gratia)
library(here)

# get colour palette to re-create published plots
source(here::here("R", "colour_palette.R"))

# create a folder for plots, if not existing already
if (!dir.exists(here::here("plots"))) {dir.create(here::here("plots"))}

# load input data 
di <- read.csv(here::here("data", "di_input-data.csv"))

# set length of sliding windows (in years)
win_length <- 5

# set nut variable for which synchrony is calculated
nut_variable <- "TotalNuts"

# get seed production from previous year
bci_t1 <- di %>% 
  dplyr::select("TreeID", "WinterYear", "TotalNuts") %>% 
  dplyr::mutate(year_plus1 = WinterYear + 1) %>% 
  dplyr::rename("TotalNuts_T1" = "TotalNuts") %>% 
  dplyr::select(!WinterYear)

# floor counts to get integers for the calculated positions & prepare factor variables for models
di <- di %>%
  dplyr::left_join(bci_t1, by = c("TreeID" = "TreeID", "WinterYear" = "year_plus1")) %>% 
  dplyr::arrange(WinterYear) %>%
  dplyr::mutate(TotalNuts = floor(TotalNuts),
                TotalWhole = floor(TotalWhole),
                TreeID = as.factor(TreeID),
                # transform WinterYear into ordinal year to avoid problems in the models
                ord_year = WinterYear - (min(WinterYear) - 1),
                # assign "probability" of a zero to each record (1 means probability of 100% for values that are indeed zero, 
                # and 0 mean probability of 0% that a value is zero for all values >0) (done for plotting purposes)
                prob_zi = dplyr::case_when(TotalNuts == 0 ~ 1,
                                           TRUE ~ 0))

# II. Trend in seed production ----------------------------------------------

# GLMM seed production over time
m_trend <- glmmTMB::glmmTMB(TotalNuts ~ ord_year + (1 | TreeID) + TotalNuts_T1, 
                            data = di, 
                            family = glmmTMB::nbinom2(link = "log"), 
                            ziformula = ~ .)

summary(m_trend)

# model diagnostics
plot(DHARMa::simulateResiduals(m_trend))

# get model predictions for zero-inflation model (Probability of producing a zero-year/no nuts)
pred_zi <- ggeffects::predict_response(m_trend, terms = c("ord_year [all]"), type = "zi_prob")

# add WinterYear back to predictions for plotting
pred_zi <- pred_zi %>% 
  dplyr::left_join(di %>% 
                     dplyr::select("WinterYear", "ord_year"), 
                   by = c("x" = "ord_year"))

# get model predictions for conditional model (How many nuts are there when there are nuts)
## note: gives warning of Jensen's inequality, setting bias_correction = TRUE does however not change the results and just takes very long
pred_cond <- ggeffects::predict_response(m_trend, terms = c("ord_year [all]"), type = "fixed")

# add WinterYear back to predictions for plotting
pred_cond <- pred_cond %>% 
  dplyr::left_join(di %>% 
                     dplyr::select("WinterYear", "ord_year"), 
                   by = c("x" = "ord_year"))

# predict combined model results (conditional model taking into account zero-inflation)
pred_trend <- ggeffects::predict_response(m_trend, terms = c("ord_year [all]"),
                                          type = "zero_inflated")


# add WinterYear back to predictions for plotting
pred_trend <- pred_trend %>% 
  dplyr::left_join(di %>% 
                     dplyr::select("WinterYear", "ord_year"), 
                   by = c("x" = "ord_year"))


### 2. Figure 1, Panel a, b, c ####

# calculate the population annual mean of beechnuts
pop_mean_nuts <- di %>%  
  dplyr::summarise(pop_mean = mean(TotalNuts), .by = WinterYear)

# Panel a: beech crop over time
p_bci <- ggplot2::ggplot(pop_mean_nuts, ggplot2::aes(x = WinterYear, y = pop_mean)) +
  ggplot2::geom_line(data = di, ggplot2::aes(x = WinterYear, y = TotalNuts, group = as.factor(TreeID)), 
                     colour = col_per_tree, alpha = 0.2) +
  ggplot2::geom_line(colour = col_pop_mean, linewidth = 1.5) +
  ggplot2::geom_point(size = 3, colour = col_pop_mean) +
  ggplot2::theme_classic(base_size = 17) +
  ggplot2::labs(y = "Total number of beechnuts per m²", x = "Year") +
  ggplot2::scale_x_continuous(breaks = c(seq(1975, 2025, by = 5))) +
  ggplot2::scale_y_continuous(breaks = c(seq(0, max(di$TotalNuts), by = 500)))


# Panel b: probability of producing no beechnuts at all in a year (population-level)
p_zi_prob <- ggplot2::ggplot(pred_zi, ggplot2::aes(x = WinterYear, y = predicted)) +
  ggplot2::geom_point(data = di, ggplot2::aes(x = WinterYear, y = prob_zi), size = 2, alpha = 0.1, colour = col_per_tree) +
  ggplot2::geom_ribbon(ggplot2::aes(x = WinterYear, y = predicted, ymin = conf.low, ymax = conf.high), 
                       fill = nut_trend_split, alpha = 0.2) +
  ggplot2::geom_line(linewidth = 1.5, colour = nut_trend_split, alpha = 0.7) +
  ggplot2::theme_classic(base_size = 17) +
  ggplot2::labs(y = "p of zero-years", x = "Year") +
  ggplot2::scale_x_continuous(breaks = c(seq(1975, 2025, by = 10)))


# Panel c: beechnut trend conditional model (how many nuts are produced, if they are produced; population-level)
p_nut_trends <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(data = pred_cond, ggplot2::aes(x = WinterYear, y = predicted, ymin = conf.low, ymax = conf.high), 
                       fill = nut_trend_split, alpha = 0.2) +
  ggplot2::geom_line(data = pred_cond, ggplot2::aes(x = WinterYear, y = predicted), 
                     linewidth = 1.5, colour = nut_trend_split, linetype = "dashed", alpha = 0.7) +
  ggplot2::geom_point(data = pop_mean_nuts, ggplot2::aes(x = WinterYear, y = pop_mean), 
                      size = 2, colour = col_per_tree, alpha = 0.1) +
  ggplot2::geom_line(data = pred_trend, ggplot2::aes(x = WinterYear, y = predicted), linewidth = 1.5) +
  ggplot2::geom_ribbon(data = pred_trend, ggplot2::aes(x = WinterYear, y = predicted, ymin = conf.low, ymax = conf.high), 
                       alpha = 0.3) +
  ggplot2::theme_classic(base_size = 17) +
  ggplot2::labs(y = "Beechnuts per m²", x = "Year") +
  ggplot2::scale_x_continuous(breaks = c(seq(1975, 2025, by = 10)))


# combine three panels
plots_trends <- ggpubr::ggarrange(
  
  ggpubr::ggarrange(p_bci, labels = "a", font.label = list(size = 18, face = "bold")),
  ggpubr::ggarrange(p_zi_prob, p_nut_trends, ncol = 1, nrow = 2, labels = c("b", "c"), 
                    font.label = list(size = 18, face = "bold")),
  
  widths = c(1.5, 1)
)


# III. Synchrony analysis ------------------------------------------------------

## Calculate synchrony per tree per year 

# prepare the individual-tree beech crop data for pairwise correlation by pivoting to wide format
di_pivot <- di %>%
  dplyr::arrange(WinterYear) %>% 
  dplyr::mutate(TreeID = paste0("Tree_", TreeID)) %>% 
  dplyr::select("TreeID", "WinterYear", "variable" = all_of(nut_variable)) %>%
  tidyr::pivot_wider(names_from = TreeID, values_from = variable)


# calculate the pairwise pearson correlation between every pair of trees per window 

## Note: despite filtering our trees that have less than 3 observations in a window, there are some warnings 
## remaining that the standard deviation of the correlation is zero. This is mainly caused by trees that have 
## always the same value and are "perfectly synced"

synchrony <- purrr::map(.x = min(di$WinterYear):(max(di$WinterYear) - (win_length - 1)),
                        .f = ~{
                          
                          # slice data per window
                          syn_window <- di_pivot %>% 
                            dplyr::filter(WinterYear >= .x, WinterYear <= (.x + (win_length - 1))) %>% 
                            dplyr::select(!WinterYear) %>% 
                            # filter out trees that in this window have less than 3 observations
                            dplyr::select_if(!colSums(is.na(.)) > 2) %>% 
                            dplyr::select_if(!(apply(., 2, function(x) length(unique(x[!is.na(x)]))) == 1))
                          
                          # calculate pairwise pearson rank correlation for window across all trees
                          rho_matrix <- cor(syn_window, method = "pearson", use = "pairwise")
                          
                          # set diagonal of matrix to NA (self-correlation)
                          diag(rho_matrix) <- NA
                          
                          # create output matrix
                          rho_df <- as.data.frame(rho_matrix) %>% 
                            # assign the year in which the window opens to correlation coefficients
                            dplyr::mutate(year_windowOpen = .x) %>% 
                            tibble::rownames_to_column(var = "TreeID")
                          
                          
                          return(rho_df)
                          
                        }, 
                        
) %>%  dplyr::bind_rows()


# calculate the mean correlation coefficient (synchrony) per window of one tree to all other trees
df_synchrony <- tibble::tibble(year_windowOpen = synchrony$year_windowOpen, 
                               TreeID = as.factor(synchrony$TreeID),
                               meanRho = synchrony %>% 
                                 dplyr::select(!c("year_windowOpen", "TreeID")) %>% 
                                 rowMeans(na.rm = TRUE)) %>% 
  # filter out tree-year combinations with NaN (there are ideally none)
  dplyr::filter(!is.na(meanRho))


# calculate population mean of synchrony per window
pop_synchrony <- df_synchrony %>%
  dplyr::summarise(pop_meanRho = mean(meanRho),
                   pop_sdRho = sd(meanRho),
                   .by = year_windowOpen)

# model the time trend of population means of synchrony with GAM
mp_pop_sync <- mgcv::gam(pop_meanRho ~ s(year_windowOpen), data = pop_synchrony, method = "REML")

summary(mp_pop_sync)

# model diagnostics
mgcv::gam.check(mp_pop_sync)
gratia::appraise(mp_pop_sync, point_col = "steelblue", point_alpha = 0.4, n_bin = "scott")

# predict on model
sync_predict <- mgcv::predict.gam(mp_pop_sync, newdata = pop_synchrony, se.fit = TRUE) %>%
  dplyr::bind_cols()

# create data frame for plotting with model predictions
pop_synchrony <- pop_synchrony %>%
  dplyr::mutate(pred.sync = sync_predict$fit,
                pred.sync.se = sync_predict$se.fit,
                # set ymax to maximum of 1 to create correct error bars in figures (synchrony cannot be bigger than 1)
                ymax_popMean = dplyr::case_when((pop_meanRho + pop_sdRho) > 1 ~ 1,
                                                TRUE ~ (pop_meanRho + pop_sdRho)))

## Figure 1, panel d: plot the population mean of synchrony ####
p_sync <- ggplot2::ggplot(pop_synchrony, ggplot2::aes(x = year_windowOpen, y = pop_meanRho)) +
  ggplot2::geom_point(size = 1, data = df_synchrony, ggplot2::aes(y = meanRho), alpha = 0.4, colour = col_per_tree) +
  ggplot2::geom_errorbar(ggplot2::aes(ymax = ymax_popMean, ymin = (pop_meanRho - pop_sdRho)), alpha = 0.8) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_line(ggplot2::aes(y = pred.sync), linewidth = 1.5, colour = col_pop_mean) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = pred.sync - pred.sync.se, ymax =  pred.sync + pred.sync.se), alpha = 0.3) +
  ggplot2::labs(x = "Start year of 5-year window",
                y = "Between-tree synchrony") +
  ggplot2::theme_classic(base_size = 17) +
  ggplot2::scale_x_continuous(breaks = c(seq(1975, 2025, by = 5))) + 
  ggplot2::geom_hline(yintercept = 0, linetype = "dotted") 



# IV. Inter annual variation (CVi) -----------------------------------------

# calculate CVi for every window and tree
df_CVi_win <- purrr::map(.x = min(di$WinterYear):(max(di$WinterYear) - (win_length - 1)),
                         .f = ~{
                           
                           df <-  di %>% 
                             dplyr::filter(WinterYear >= .x, WinterYear <= (.x + (win_length - 1))) 
                           
                           include_trees <- df %>%
                             dplyr::group_by(TreeID) %>% 
                             dplyr::count() %>% 
                             dplyr::filter(n > 2) %>% 
                             dplyr::pull(TreeID)
                           
                           df %>% 
                             dplyr::filter(TreeID %in% include_trees) %>% 
                             dplyr::summarise(mean_Total = mean(TotalNuts, na.rm = TRUE),
                                              sd_Total = sd(TotalNuts, na.rm = TRUE),
                                              .by = TreeID) %>% 
                             dplyr::mutate(p.CVi = sd_Total / mean_Total, # Pearson CV
                                           k.CVi = sqrt((p.CVi^2)/(1 + p.CVi^2)), # Kvålseth's CV
                                           year_windowOpen = .x)
                           
                         }
) %>%  dplyr::bind_rows()

# to get one value for CVi (as done in literature) take the mean CVi per window 
df_CVi <- df_CVi_win %>% 
  dplyr::filter(!(is.na(p.CVi) | is.na(k.CVi))) %>% 
  dplyr::mutate(TreeID = as.factor(TreeID)) %>% 
  dplyr::summarise(mean_p.CVi = mean(p.CVi), 
                   sd_p.CVi = sd(p.CVi), 
                   mean_k.CVi = mean(k.CVi), 
                   sd_k.CVi = sd(k.CVi),
                   .by = year_windowOpen)


### Pearson CVi (reported in supplements) ####

# model the temporal change in mean CVi
gam_p.CVi <- mgcv::gam(mean_p.CVi ~ s(year_windowOpen), data = df_CVi, method = "REML")

summary(gam_p.CVi)

# model diagnostics
gratia::draw(gam_p.CVi, residuals = TRUE)
mgcv::gam.check(gam_p.CVi)
gratia::appraise(gam_p.CVi, point_col = "steelblue", point_alpha = 0.4, n_bin = "scott")

# prediction based on gam
mean_p.CVi_predict <- mgcv::predict.gam(gam_p.CVi, newdata = df_CVi, se.fit = TRUE) %>%
  dplyr::bind_cols()

df_CVi <- df_CVi %>%
  dplyr::mutate(pred_p.CVi = mean_p.CVi_predict$fit,
                pred_p.CVi.se = mean_p.CVi_predict$se.fit)

## Figure S1 - Pearson CVi ####
p_p.CVi <- ggplot2::ggplot(df_CVi, ggplot2::aes(x = year_windowOpen)) +
  ggplot2::geom_point(data = df_CVi_win, ggplot2::aes(y = p.CVi), alpha = 0.4, colour = col_per_tree, size = 1) +
  ggplot2::geom_errorbar(ggplot2::aes(ymax = mean_p.CVi + sd_p.CVi, ymin = (mean_p.CVi - sd_p.CVi)), alpha = 0.8) +
  ggplot2::geom_point(ggplot2::aes(y = mean_p.CVi), size = 3) +
  ggplot2::geom_line(ggplot2::aes(y = pred_p.CVi), colour = col_pop_mean, linewidth = 1.5) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = pred_p.CVi - pred_p.CVi.se, ymax =  pred_p.CVi + pred_p.CVi.se), alpha = 0.3) +
  ggplot2::scale_x_continuous(breaks = c(seq(1975, 2025, by = 5))) +
  ggplot2::labs(x = "Start year of 5-year window",
                y = "Pearson CVi (mean)") +
  ggplot2::theme_classic(base_size = 17) 

# save figure for supplementary material (Figure S1)
ggplot2::ggsave(plot = p_p.CVi,  filename = here::here("plots", "Figure_S1.png"), 
                width = 20, height = 12, unit = "cm")



### Kvålseth CVi (reported in main text) ####

# model the temporal change in mean CVi
gam_k.CVi <- mgcv::gam(mean_k.CVi ~ s(year_windowOpen), data = df_CVi, method = "REML")

summary(gam_k.CVi)

# model diagnostics
gratia::draw(gam_k.CVi, residuals = TRUE)
mgcv::gam.check(gam_k.CVi)
gratia::appraise(gam_k.CVi, point_col = "steelblue", point_alpha = 0.4, n_bin = "scott")

# prediction based on gam
mean_k.CVi_predict <- mgcv::predict.gam(gam_k.CVi, newdata = df_CVi, se.fit = TRUE) %>%
  dplyr::bind_cols()

# add predictions to data frame
df_CVi <- df_CVi %>%
  dplyr::mutate(pred_k.CVi = mean_k.CVi_predict$fit,
                pred_k.CVi.se = mean_k.CVi_predict$se.fit,
                # manually assign maximum value of 1 to SD to create correct error bars in plot (k.CVi cannot be bigger than 1)
                ymax_sd_k.CVi = dplyr::case_when((mean_k.CVi + sd_k.CVi) > 1 ~ 1,
                                                 TRUE ~ (mean_k.CVi + sd_k.CVi)))

## Figure 1 - panel e: Kvålseth CVi ####
p_k.CVi <- ggplot2::ggplot(df_CVi, ggplot2::aes(x = year_windowOpen)) +
  ggplot2::geom_point(data = df_CVi_win, ggplot2::aes(y = k.CVi), alpha = 0.4, colour = col_per_tree, size = 1) +
  ggplot2::geom_errorbar(ggplot2::aes(ymax = ymax_sd_k.CVi, ymin = (mean_k.CVi - sd_k.CVi)), alpha = 0.8) +
  ggplot2::geom_point(ggplot2::aes(y = mean_k.CVi), size = 3) +
  ggplot2::geom_line(ggplot2::aes(y = pred_k.CVi), colour = col_pop_mean, linewidth = 1.5) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = pred_k.CVi - pred_k.CVi.se, ymax =  pred_k.CVi + pred_k.CVi.se), alpha = 0.3) +
  ggplot2::scale_x_continuous(breaks = c(seq(1975, 2025, by = 5))) +
  ggplot2::labs(x = "Start year of 5-year window",
                y = "Kvålseth CVi (mean)") +
  ggplot2::theme_classic(base_size = 17) 


# V. Make final figure 1 -------------------------------------------------------

# Combine beech crop & trend, synchrony and CVi plots into final figure
ggpubr::ggarrange(
  
  plots_trends, 
  ggpubr::ggarrange(p_sync, p_k.CVi, labels = list("d", "e"), 
                    font.label = list(size = 18, face = "bold")),
  
  nrow = 2, ncol = 1, heights = c(1.5, 1)
  
)

# save Figure 1 of main text
ggplot2::ggsave(plot = ggplot2::last_plot(), filename = here::here("plots", "Figure_1.png"), 
                width = 45, height = 30, unit = "cm")
