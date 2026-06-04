## Processing of results from sensitivity window analysis ####

# Author:           Cherine C. Jantzen 
# Created:          2025-09-11
# Last updated:     2026-02-02


################################
## Content: This script processes the results of the weather cue (i.e., sensitivity windows) identification for each climate variable. It also analyses 
## temporal trends in the weather cues,the effects of the weather cues on annual seed production and the role of resource reserves.
################################

# I. Preparation ----------------------------------------------------------

# load packages
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(rstudioapi)
library(glmmTMB)
library(lme4)
library(corrplot)
library(here)
library(performance)

# get colour palette for plots
source(here::here("R", "colour_palette.R"))

# get summary table of results of all tested windows in all iterations
summary_allModels <- read.csv(here::here("data", "summary_all_models.csv"))

# get table with best windows of each variable of all iterations
sens_wind <- read.csv(here::here("data", "table_sens_windows.csv"))

# get data frame with all climate variables set to their respective best windows of 4th iteration
di_win_4It <- read.csv(here::here("data", "di_after4thIt.csv")) %>% 
  dplyr::mutate(TreeID = as.factor(TreeID))

# load non-transformed climate data for each window
all_windows <- read.csv(here::here("data", "all_windows_climate.csv"))

# load beechnut input data
di <- read.csv(here::here("data", "di_input-data.csv"))

# create a folder for plots, if not existing already
if (!dir.exists(here::here("plots"))) {dir.create(here::here("plots"))}


# Visualize sliding window results ----------------------------------------

## Figure 2 - panel a: Position of best windows in the year ####

# position an length of best windows per variable
p_windows <- sens_wind %>% 
  # include only results of fourth iteration
  dplyr::filter(stringr::str_detect(variable, pattern = "4It")) %>% 
  tidyr::pivot_longer(cols = c(Start_DOY, End_DOY), names_to = "DOY") %>% 
  # rename climate variables for easier plotting
  dplyr::mutate(variable = dplyr::case_when(variable == "z_meanTempT0_4It" ~ "meanTempT0",
                                            variable == "z_maxTempT1_4It" ~ "maxTempT1",
                                            variable == "z_maxTempT2_4It" ~ "maxTempT2",
                                            variable == "z_PrecT0_4It" ~ "PrecT0", 
                                            variable == "z_PrecT1_4It" ~ "PrecT1",  
                                            variable == "z_PrecT2_4It" ~ "PrecT2"), 
                # convert variables to factor for desired order in plot
                fct_variable = factor(variable, 
                                      levels = c("PrecT2", "PrecT1", "PrecT0", "maxTempT2","maxTempT1", "meanTempT0")),
                # add calendar dates for start and end of window for labelling
                cal_Start = case_when(DOY == "Start_DOY" ~ paste0(substring(as.Date(value, origin = "2025-01-01"), first = 9), 
                                                                  "/", substring(as.Date(value, origin = "2025-01-01"), 
                                                                                 first = 6, last = 7)),
                                      TRUE ~ NA),
                cal_End = case_when(DOY == "End_DOY" ~ paste0(substring(as.Date(value, origin = "2025-01-01"), first = 9), 
                                                              "/", substring(as.Date(value, origin = "2025-01-01"), 
                                                                             first = 6, last = 7)),
                                    TRUE ~ NA)) %>% 
  ggplot2::ggplot(ggplot2::aes(x = value, y = fct_variable, colour = fct_variable, shape = fct_variable)) +
  # add invisible points for creating a common legend with other panels of this figure
  ggplot2::geom_point(size = 0, alpha = 0)  +
  ggplot2::geom_line(linewidth = 2) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::labs(x = "Day of year", y = "", colour = "", shape = "") +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::xlim(c(79, 265)) +
  ggplot2::scale_colour_manual(values = col_pal) +
  ggplot2::scale_shape_manual(values = shape_pal) +
  ggplot2::geom_vline(xintercept = 172, linetype = "dotted") +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1,
                                                 override.aes = list(linewidth = 1)),
                  shape = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2))) +
  # add labels for calendar dates to the data points
  ggplot2::geom_text(ggplot2::aes(label = cal_Start), hjust = 0.7, vjust = 2, show.legend = FALSE, size = 1.5) +
  ggplot2::geom_text(ggplot2::aes(label = cal_End), hjust = 0.3, vjust = 2, show.legend = FALSE, size = 1.5) 


## AIC patterns of sliding window runs per climate variable ####

# plot the AIC patterns of each climate variable from the fourth iteration
summary_4thIt <- summary_allModels %>%
  dplyr::filter(stringr::str_detect(variable, pattern = "4It")) %>% 
  dplyr::mutate(variable = dplyr::case_when(variable == "z_meanTempT0_4It" ~ "meanTempT0",
                                            variable == "z_maxTempT1_4It" ~ "maxTempT1",
                                            variable == "z_maxTempT2_4It" ~ "maxTempT2",
                                            variable == "z_PrecT0_4It" ~ "PrecT0",
                                            variable == "z_PrecT1_4It" ~ "PrecT1",
                                            variable == "z_PrecT2_4It" ~ "PrecT2"),
                # order variables
                fct_variable = factor(variable,
                                      levels = c("meanTempT0", "maxTempT1", "maxTempT2", "PrecT0", "PrecT1", "PrecT2"))) 

## Figure S3: AIC patterns of fourth iteration ####

## plot temperature and precipitation separately, otherwise the the difference in AIC is so big that patterns in temperature become invisible

ggpubr::ggarrange(
  
  # temperature
  summary_4thIt %>%  
    dplyr::filter(variable %in% c("meanTempT0", "maxTempT1", "maxTempT2")) %>% 
    ggplot2::ggplot(ggplot2::aes(x = Start_DOY, y = AIC, colour = End_DOY)) +
    ggplot2::geom_point(size = 1, alpha = 0.8) +
    ggplot2::theme_minimal(base_size = 17) +
    ggplot2::guides(colour = ggplot2::guide_colourbar(position = "bottom",
                                                      theme = theme(legend.key.width  = unit(10, "lines")))) +
    ggplot2::scale_colour_viridis_c(option = "D", end = 0.9) +
    ggplot2::labs(x = "Start DOY of window", colour = "End DOY of window") +
    ggplot2::facet_wrap(~ fct_variable) +
    ggplot2::scale_x_continuous(breaks = c(seq(80, 270, by = 50))) +
    ggplot2::scale_y_continuous(breaks = c(seq(11900, 12150, by = 50))),
  
  # precipitation
  summary_4thIt %>%  
    dplyr::filter(variable %in% c("PrecT0", "PrecT1", "PrecT2")) %>% 
    ggplot2::ggplot(ggplot2::aes(x = Start_DOY, y = AIC, colour = End_DOY)) +
    ggplot2::geom_point(size = 1, alpha = 0.8) +
    ggplot2::theme_minimal(base_size = 17) +
    ggplot2::guides(colour = ggplot2::guide_colourbar(position = "bottom",
                                                      theme = theme(legend.key.width  = unit(10, "lines")))) +
    ggplot2::scale_colour_viridis_c(option = "D", end = 0.9) +
    ggplot2::labs(x = "Start DOY of window", colour = "End DOY of window") +
    ggplot2::facet_wrap(~ fct_variable) +
    ggplot2::scale_x_continuous(breaks = c(seq(80, 270, by = 50))) +
    ggplot2::scale_y_continuous(breaks = c(seq(11900, 12500, by = 100))),
  
  nrow = 2, ncol = 1, common.legend = TRUE, legend = "bottom"
)

ggplot2::ggsave(plot = ggplot2::last_plot(),
                file =  here::here("plots", "Figure_S3.png"),
                units = "cm", dpi = 600, width = 40, height = 25)


# Temporal changes in climatic drivers in respective windows ----------------------

# get real values for climate variables in best windows
d_tempInWin <- purrr::map2(.x = c("z_maxTempT1_4It", "z_maxTempT2_4It", "z_meanTempT0_4It", "z_PrecT1_4It", "z_PrecT0_4It", "z_PrecT2_4It"),
                           .y = c("meanMaxWinTemp", "meanMaxWinTemp", "meanMeanWinTemp", "sumWinPrec", "sumWinPrec", "sumWinPrec"),
                           .f = ~{
                             
                             
                             start <- sens_wind %>% 
                               dplyr::filter(variable == .x) %>% 
                               pull(Start_DOY)
                             
                             end <- sens_wind %>% 
                               dplyr::filter(variable == .x) %>% 
                               pull(End_DOY)
                             
                             climate <- all_windows %>% 
                               dplyr::arrange(year) %>% 
                               dplyr::filter(Start_DOY == start,
                                             End_DOY == end) %>% 
                               dplyr::select(all_of(.y))
                             
                             variable_name <- sapply(strsplit(.x, "_"), function(x) x[2])
                             colnames(climate)[1] <- variable_name
                             
                             return(climate)
                             
                           }) %>% dplyr::bind_cols() %>% 
  dplyr::bind_cols(tibble::tibble(year = c(1974:2025),
                                  ord_year = c(1:length(year))))


# re-calculate z-scores from real climate data to use all years (1974 to 2025) in models & predictions
z_climate <- d_tempInWin %>% 
  dplyr::mutate(z_tempT1 = scale(maxTempT1, center = TRUE, scale = TRUE)[,1],
                z_tempT0 = scale(meanTempT0, center = TRUE, scale = TRUE)[,1],
                z_tempT2 = scale(maxTempT2, center = TRUE, scale = TRUE)[,1],
                z_precT0 = scale(PrecT0, center = TRUE, scale = TRUE)[,1],
                z_precT1 = scale(PrecT1, center = TRUE, scale = TRUE)[,1],
                z_precT2 = scale(PrecT2, center = TRUE, scale = TRUE)[,1])


# Model temporal trends of each climate variable with linear model #####

# temperature in T1 
lm_z_tempT1 <- lm(z_tempT1 ~ ord_year, data = z_climate)

summary(lm_z_tempT1)
# slightly significant

# temperature in T2 
lm_z_tempT2 <- lm(z_tempT2 ~ ord_year, data = z_climate)

summary(lm_z_tempT2)
# non significant

# temperature in T0 
lm_z_tempT0 <- lm(z_tempT0 ~ ord_year, data = z_climate)

summary(lm_z_tempT0)
# significant

# precipitation in T1 
lm_z_precT1 <- lm(z_precT0 ~ ord_year, data = z_climate)

summary(lm_z_precT1)
# not significant

# precipitation in T2 
lm_z_precT2 <- lm(z_precT1 ~ ord_year, data = z_climate)

summary(lm_z_precT2)
# not significant

# precipitation in T0 
lm_z_precT0 <- lm(z_precT2 ~ ord_year, data = z_climate)

summary(lm_z_precT0)
# not significant


# add predictions for plotting of significant terms
z_climate$pred.tempT1 <- predict(lm_z_tempT1, newdata = z_climate)
z_climate$pred.tempT0 <- predict(lm_z_tempT0, newdata = z_climate)

# pivot data for easier plotting
df_plot_clim <- d_tempInWin %>% 
  dplyr::left_join(z_climate %>%  
                     dplyr::select("year", "ord_year", "pred.tempT1", "pred.tempT0") %>% 
                     dplyr::distinct(), 
                   by = "year") %>% 
  dplyr::mutate(pred.realtempT1 = mean(maxTempT1) + (pred.tempT1 * sd(maxTempT1)),
                pred.realtempT0 = mean(meanTempT0) + (pred.tempT0 * sd(meanTempT0))) %>% 
  tidyr::pivot_longer(cols = c("meanTempT0", "maxTempT1"), names_to = "temperature_variable")



# Figure 2: panel b (trends in temperature) #####
panel_b_temp <- ggplot2::ggplot(df_plot_clim, 
                                ggplot2::aes(x = year, y = value, colour = temperature_variable, shape = temperature_variable)) +
  ggplot2::geom_point(size = 1.25) +
  ggplot2::geom_line(ggplot2::aes(y = pred.realtempT1), linewidth = 0.75, colour = col_pal %>% pluck("maxTempT1")) +
  ggplot2::geom_line(ggplot2::aes(y = pred.realtempT0), linewidth = 0.75, colour = col_pal %>% pluck("meanTempT0")) +
  ggplot2::geom_point(data = d_tempInWin %>% 
                        tidyr::pivot_longer(cols = "maxTempT2", names_to = "temperature_variable"), size = 1.25) +
  ggplot2::geom_smooth(data = d_tempInWin %>% 
                         tidyr::pivot_longer(cols = "maxTempT2", names_to = "temperature_variable"), method = "lm", linetype = 8, se = FALSE, linewidth = 0.75) +
  ggplot2::labs(x = "Year", y = "Mean daily temperature [°C]", colour = "", shape = "") +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::scale_colour_manual(values = col_pal) + 
  ggplot2::scale_shape_manual(values = shape_pal) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1)) +
  ggplot2::scale_y_continuous(breaks = c(seq(0, 40, by = 5))) +
  ggplot2::scale_x_continuous(breaks = c(seq(1975, 2025, by = 10)))


# Figure 2: panel c (trends in precipitation) #####
panel_c_temp <- ggplot2::ggplot(d_tempInWin %>% 
                                  tidyr::pivot_longer(cols = c("PrecT0", "PrecT1", "PrecT2"), names_to = "prec_variable"), 
                                ggplot2::aes(x = year, y = value, colour = prec_variable, shape = prec_variable)) +
  ggplot2::geom_point(size = 1.25) +
  ggplot2::labs(x = "Year", y = "Daily precipitation sum [mm]", colour = "", shape = "") +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::scale_y_continuous(breaks = c(seq(0, 350, by = 50))) + 
  ggplot2::scale_x_continuous(breaks = c(seq(1975, 2025, by = 10))) +
  ggplot2::scale_colour_manual(values = col_pal) + 
  ggplot2::scale_shape_manual(values = shape_pal) +
  ggplot2::geom_smooth(method = "lm", linetype = 8, se = FALSE, linewidth = 0.75) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1))


# Make Figure 3 --------------------------------------------------------

# tweaking the aesthetics of the panels for nicer arrangement

## extract legend for climate variables
leg_climate <- cowplot::get_legend(p_windows + ggplot2::theme(legend.position = "bottom"))

## remove legend from panel a
p_windows <- p_windows + ggplot2::theme(legend.position = "none")

## convert legend in plotable object for ggarrange
leg_climate_plot <- ggpubr::as_ggplot(leg_climate)

## arrange panels and legend
ggpubr::ggarrange(
  
  ggpubr::ggarrange(
    
    # panel a
    ggpubr::ggarrange(p_windows, labels = "a", font.label = list(size = 9)),
    
    ggpubr::ggarrange(
      
      panel_b_temp, panel_c_temp,
      
      legend = "none", ncol = 1, nrow = 2, labels = c("b", "c"), font.label = list(size = 9)),
    widths = c(1.5, 1)),
  
  # add common legend as "plot"
  leg_climate_plot,
  
  ncol = 1, nrow = 2, heights = c(1, 0.12))


ggplot2::ggsave(plot = ggplot2::last_plot(), 
                file = here::here("plots", "Figure_3.png"), 
                units = "mm", height = 129, width = 180, dpi = 900)

# Analyse drivers of seed production -------------------------------------------

# make the final model: with all climate variables set to values in their best windows
m_final <- glmmTMB::glmmTMB(TotalNuts ~
                              z_temp_meanMean_4It_T0 + z_temp_meanMax_4It_T1 + z_temp_meanMax_4It_T2 +
                              z_prec_4It_T0 +  z_prec_4It_T1 + z_prec_4It_T2 + 
                              TotalNuts_T1 + (1|TreeID),
                            data = di_win_4It,
                            family = nbinom2(link = "log"),
                            ziformula = ~ .)

summary(m_final)

# model diagnostics
diagnose(m_final)
plot(DHARMa::simulateResiduals(m_final))

# correct p-values for the number of windows tested & see if they are still significant 
bind_rows(tibble::enframe(coef(summary(m_final))$cond [, 4], name =  "variable", value = "p"),
          tibble::enframe(coef(summary(m_final))$zi [, 4], name = "variable", value = "p")) %>% 
  mutate(corr_signlevel = 0.05 / 15343, # Option 1: correct significance level for number of tests
         sign = p < corr_signlevel, # check whether p-values are below corrected significance level
         cor_pvalue = p.adjust(p, n = 15343, method = "bonferroni"), # Option 2: Bonferroni correction of p-values
         check_cor_pvalue = cor_pvalue < 0.001) # check whether correct p-values are still highly significant


## compare final model with base model (model with windows defined in literature)

m_lit <- di %>% 
  dplyr::mutate(TreeID = as.factor(TreeID)) %>% 
  dplyr::left_join(di_win_4It %>%  dplyr::select(WinterYear, TotalNuts_T1, TreeID), by = c("WinterYear", "TreeID")) %>% 
  glmmTMB::glmmTMB(TotalNuts ~
                   z_temp_meanGrow + z_maxTempSummer_T1 + z_maxTempSummer_T2 +
                   z_prec_Grow +  z_prec_Summer_T1 + z_prec_Summer_T2 + 
                   TotalNuts_T1 + (1|TreeID),
                 data = .,
                 family = nbinom2(link = "log"),
                 ziformula = ~ .)

summary(m_lit)
### prec_Grow and prec_T1 are no longer significant in both parts of the model, maxTemp_T2 no longer in zero-inflated part of the model


# compare both models based on their AIC
AIC(m_final) - AIC(m_lit)

# compare the variance explained by the the literature- and the final model 
## look at the marginal R^2 as that includes random and fixed effects; 
## component = "all" includes conditional and zero-inflated part of the model
performance::r2(m_final, component = "all")
performance::r2(m_lit, component = "all")

## Check model fit by comparing observed with predicted values ####

# get model data frame, as that ommited some rows 
fit_final <- m_final[["frame"]] %>% 
  dplyr::left_join(di_win_4It, 
                   by = c("TreeID", "TotalNuts", "z_temp_meanMean_4It_T0", "z_temp_meanMax_4It_T1", 
                          "z_temp_meanMax_4It_T2", "z_prec_4It_T0", "z_prec_4It_T1", "z_prec_4It_T2")) 

# predict number nuts on final model under observed climate
fit_final$pred.nuts_realClim <- predict(m_final, type = "response")

## Figure 4 - panel a: Observed and fitted values ####
Fig_4a <- ggplot2::ggplot(fit_final, ggplot2::aes(x = TotalNuts, y = pred.nuts_realClim, colour = WinterYear)) +
  ggplot2::geom_point(size = 1.5, alpha = 0.5) +
  ggplot2::xlim(0, 2600) +
  ggplot2::ylim(0, 2600) +
  ggplot2::scale_colour_viridis_c(option = "H") +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  ggplot2::labs(y = "Predicted number nuts [per m²]", x = "Observed number nuts [per m²]", colour = "Year") +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(legend.position = c(0.15, 0.85),
                 legend.title = element_text(size = 9), 
                 legend.text = element_text(size = 7),
                 legend.key.height = unit(0.45, 'cm'),
                 legend.key.width = unit(0.4, 'cm')) 

# get correlation coefficient of observed and fitted values
cor(fit_final$TotalNuts, fit_final$pred.nuts_realClim, method = "spearman")
# coef = 0.861


# calculate population means 
pop_Total <- fit_final %>%  
  dplyr::summarise(pop_TotalNuts = mean(TotalNuts, na.rm = TRUE),
                   pop_pred.nuts_realClim = mean(pred.nuts_realClim, na.rm = TRUE),
                   .by = WinterYear)

# get correlation coefficient of observed and fitted values
cor(pop_Total$pop_TotalNuts, pop_Total$pop_pred.nuts_realClim)
# coef = 0.9060949


## Figure 4 - panel b: Temporal patterns of model predictions and observed values
Fig_4b <- ggplot2::ggplot(fit_final, ggplot2::aes(x = WinterYear)) +
  ggplot2::geom_point(ggplot2::aes(y = pred.nuts_realClim, group = TreeID, colour = "Model predictions"), alpha = 0.2, size = 0.5) +
  ggplot2::geom_point(ggplot2::aes(y = TotalNuts, group = TreeID, colour = "Observed"), alpha = 0.2, size = 0.5) +
  ggplot2::geom_line(ggplot2::aes(y = pred.nuts_realClim, group = TreeID, colour = "Model predictions"), linewidth = 0.25, alpha = 0.1) +
  ggplot2::geom_line(ggplot2::aes(y = TotalNuts, group = TreeID, colour = "Observed"), linewidth = 0.25, alpha = 0.1) +
  ggplot2::labs(x = "Year", y = "Total number beechnuts [per m²]", colour = "") +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::geom_line(data = pop_Total, ggplot2::aes(y = pop_pred.nuts_realClim, colour = "Model predictions"), linewidth = 0.75) +
  ggplot2::geom_line(data = pop_Total, ggplot2::aes(y = pop_TotalNuts, colour = "Observed"), linewidth = 0.75) +
  ggplot2::scale_x_continuous(breaks = c(seq(1975, 2025, by = 10))) +
  ggplot2::scale_colour_manual(values = c("Model predictions" = col_pred.real_clim, "Observed" = col_per_tree)) +
  ggplot2::theme(legend.position = c(0.85, 0.95),
                 legend.title = element_text(size = 9), 
                 legend.text = element_text(size = 7))


## Figure 4: arrange panels to make final figure ####
ggpubr::ggarrange(Fig_4a, Fig_4b, 
                  ncol = 2, nrow = 1, labels = "auto", widths = c(0.7, 1), font.label = list(size = 9))

ggplot2::ggsave(plot = ggplot2::last_plot(), 
                file = here::here("plots", "Figure_4.png"), 
                units = "mm", width = 180, height = 110, dpi = 900)


# Test for unexplained variation in the model -----------------------------


## Test for remaining year effect by fitting year back into the model ####

m_final_year <- glmmTMB::glmmTMB(TotalNuts ~
                                   z_temp_meanMean_4It_T0 +z_temp_meanMax_4It_T1 + z_temp_meanMax_4It_T2 +
                                   z_prec_4It_T0 +  z_prec_4It_T1 + z_prec_4It_T2 + ord_year +
                                   TotalNuts_T1 + (1|TreeID),
                                 data = di_win_4It,
                                 family = glmmTMB::nbinom2(link = "log"),
                                 ziformula = ~ .)

summary(m_final_year)

# model diagnostics
diagnose(m_final_year)
plot(DHARMa::simulateResiduals(m_final_year))

# correct p-values for the number of windows tested & see if they are still highly significant (<0.001)
bind_rows(tibble::enframe(coef(summary(m_final_year))$cond [, 4], name =  "variable", value = "p"),
          tibble::enframe(coef(summary(m_final_year))$zi [, 4], name = "variable", value = "p")) %>% 
  mutate(corr_signlevel = 0.05 / 15343,
         sign = p < corr_signlevel,
         cor_pvalue = p.adjust(p, n = 15343, method = "bonferroni"), # Option 2: Bonferroni correction of p-values
         check_cor_pvalue = cor_pvalue < 0.001)


# extract data frame on which the model was fit 
## (exclusion of some values by the model cause issues in the model prediction in the next step when using the input file to predict)
fit_final_year <- m_final_year[["frame"]] 

# make model predictions based on all observed values (fitted in the model) 
fit_final_year$pred.nuts_realClim <- predict(m_final_year, type = "response")

# calculate correlation coefficient between observed and fitted
cor(fit_final_year$TotalNuts, fit_final_year$pred.nuts_realClim)
# coef =  0.7684516

# visually compare observed and predicted 
ggplot2::ggplot(fit_final_year, ggplot2::aes(x = ord_year, y = pred.nuts_realClim, group = TreeID)) +
  ggplot2::geom_line() +
  ggplot2::geom_line(ggplot2::aes(y = TotalNuts), alpha = 0.5, colour = "steelblue")


# Supplementary analysis --------------------------------------------------


## Supporting Information D: Calculate resource reserves ####
## Following analysis in: Kelly et al., 2025

# get ordinal year (as a subsitute for cummulative years, because it better accounts for the missing year in the data)
ord_year <- di_win_4It %>%
  dplyr::select(WinterYear, ord_year) %>% 
  dplyr::distinct(ord_year, .keep_all = TRUE)

# calculate the cumulative total nuts produced over time
di <- di %>% 
  # make sure values per tree are arranged by year to calculalte cumsum correctly
  dplyr::arrange(TreeID, WinterYear) %>% 
  dplyr::mutate(TreeID = as.factor(TreeID),
                cum.totalNuts = cumsum(TotalNuts), .by = TreeID) %>% 
  dplyr::left_join(ord_year, by = "WinterYear") 

# look at the cumulative number of nuts over time
ggplot2::ggplot(di, ggplot2::aes(x = as.numeric(WinterYear), y = cum.totalNuts, group = TreeID, color = TreeID)) +
  ggplot2::geom_line() +
  ggplot2::theme(legend.position = "none")

# fit linear mixed model to calculate the resource use
lmer_resources <- lme4::lmer(cum.totalNuts ~ -1 + ord_year + (-1 + ord_year|TreeID),
                             data = di)

# extract the negative of residuals as resource value ( to account for the inverse relationship of stored resources and the residuals)
di$resources <- -resid(lmer_resources)

# add resources from the previous year to the seed production of this year
di <- di %>% 
  dplyr::mutate(resources = lag(resources, n = 1), .by = TreeID)

## look at how resource reserves have changed over time
ggplot2::ggplot(di, ggplot2::aes(x = as.numeric(WinterYear), y = resources, group = TreeID, color = TreeID)) +
  ggplot2::geom_line() +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::theme(legend.position = "none")

## look how reproductive effort changes with resource reserves 
ggplot2::ggplot(di, ggplot2::aes(x = resources, y = TotalNuts, colour = TreeID)) +
  ggplot2::geom_point(alpha = 0.7) +
  ggplot2::theme(legend.position = "none")

# add calculated resources to the input data for the final model
df_final_resource <- di_win_4It %>% 
  dplyr::left_join(di %>% 
                     dplyr::select("WinterYear", "TreeID", "resources"), 
                   by = c("WinterYear" = "WinterYear", "TreeID" = "TreeID"))


#  add resource reserves to the final model as interaction effect with tempT1
m_final_resource <- glmmTMB::glmmTMB(TotalNuts ~
                                       z_temp_meanMean_4It_T0 + z_temp_meanMax_4It_T1 * resources + z_temp_meanMax_4It_T2 +
                                       z_prec_4It_T0  +  z_prec_4It_T1 + z_prec_4It_T2 + 
                                       TotalNuts_T1 + (1|TreeID),
                                     data = df_final_resource,
                                     family = glmmTMB::nbinom2(link = "log"),
                                     # for zero-inflation part, we don't expect an interaction but rather an additive effect of resources
                                     ziformula = ~ z_temp_meanMean_4It_T0 + z_temp_meanMax_4It_T1 + resources + z_temp_meanMax_4It_T2 +
                                       z_prec_4It_T0  +  z_prec_4It_T1 + z_prec_4It_T2 + TotalNuts_T1 + (1|TreeID))

summary(m_final_resource)


# drop interaction effect and refit model with additive effect (in conditional and zero-inflated part)
m_final_resource2 <- glmmTMB::glmmTMB(TotalNuts ~
                                        z_temp_meanMean_4It_T0 + z_temp_meanMax_4It_T1 + z_temp_meanMax_4It_T2 +
                                        z_prec_4It_T0  +  z_prec_4It_T1 + z_prec_4It_T2 + 
                                        resources + TotalNuts_T1 + (1|TreeID),
                                      data = df_final_resource,
                                      family = glmmTMB::nbinom2(link = "log"),
                                      ziformula = ~ .)

summary(m_final_resource2)

# fit interaction of resource and tempT1 only (excluding all other variables)

m_final_resource3 <- glmmTMB::glmmTMB(TotalNuts ~
                                        z_temp_meanMax_4It_T1 * resources + z_temp_meanMax_4It_T2 +
                                        TotalNuts_T1 + (1|TreeID),
                                      data = df_final_resource,
                                      family = glmmTMB::nbinom2(link = "log"),
                                      ziformula = ~ z_temp_meanMax_4It_T1 + resources + z_temp_meanMax_4It_T2 +
                                        TotalNuts_T1 + (1|TreeID))


summary(m_final_resource3)

## Supporting information C.4.: Check for colinearity of the climate explanatory variables ####

# reduce climate variables of final model to only one observation per year 
## (because di_win_4It duplicates each climate variable per year for each tree)

df_cor <- di_win_4It %>% 
  dplyr::summarise(tempT1 = mean(z_temp_meanMax_4It_T1),
                   tempT2 = mean(z_temp_meanMax_4It_T2),
                   tempT0 = mean(z_temp_meanMean_4It_T0),
                   precT1 = mean(z_prec_4It_T1),
                   precT2 = mean(z_prec_4It_T2),
                   precT0 = mean(z_prec_4It_T0),
                   .by = WinterYear) %>%
  # drop non-climate variable(s)
  dplyr::select(!WinterYear)

# calculate correlation matrix (Pearson correlation)
cor_matrix <- cor(df_cor)

# Figure S5: Visualize colinearity matrix ####

# as we use the corrplot package to do so, which is based on base R, we cannot use ggplot
png(here::here("plots", "Figure_S5.png"), 
    width = 15, height = 10, units = "cm", res = 600, bg = "white")

corrplot::corrplot(cor_matrix, method = "color", type = "lower", 
                   diag = FALSE, tl.col = "black", addCoef.col = "black")

dev.off()
