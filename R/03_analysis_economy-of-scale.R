# Analysis: Economies of scale ####

# Author:       Cherine C. Jantzen
# Created:      2025-02-07
# Last updated: 2026-02-02

################################
## Content: This script looks at the selective benefits of masting by analysing temporal trends in the ratio of predated and pollinated nuts, and by looking at 
## three hypothesised economies of scale: starvation, satiation and pollination efficiency. 
################################

# load libraries 
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(ggplot2)
library(ggpubr)
library(here)
library(ggeffects)
library(cowplot)
library(grid) 

# get colour palette for plots
source(here::here("R", "colour_palette.R"))

# load data & floor counts from years that correct for missing plot positions
di <- read.csv(here::here("data", "di_input-data.csv"))

# create folder for plots if not existing
if (!dir.exists(here::here("plots"))) {dir.create(here::here("plots"))}

# Prepare data ---------------------------------------------------------

# calculate population sum of total nuts per year corrected by number of trees measured
pop_Total <- di %>% 
  dplyr::summarise(sum_TotalNuts = sum(TotalNuts, na.rm = TRUE),
                   n_trees = dplyr::n(),
                   .by = WinterYear)

# calculate within-site synchrony
df_CVp <- di %>%
  dplyr::summarise(sd = sd(TotalNuts, na.rm = TRUE),
                   mean = mean(TotalNuts, na.rm = TRUE),
                   CVp = sd/mean,
                   .by = WinterYear) %>% 
  # for full zero years, it results in NaN values, but as we know there is no variation between trees then, we set it to 0
  dplyr::mutate(CVp = dplyr::case_when(CVp == "NaN" ~ 0,
                                       TRUE  ~ CVp))

# calculate economies of scale measures and create input data frame for models
eos <- di %>% 
  dplyr::mutate(year_min1 = WinterYear - 1) %>%
  # add number TotalNuts of previous year to data
  dplyr::left_join(di %>% 
                     dplyr::select("WinterYear", "TotalNuts_T1" = "TotalNuts", "TreeID"), 
                   by = c("year_min1" = "WinterYear", "TreeID" = "TreeID")) %>% 
  # add population means to data
  dplyr::left_join(pop_Total, by = "WinterYear") %>% 
  # add within-year synchrony to data
  dplyr::left_join(df_CVp %>%  
                     dplyr::select("WinterYear", "CVp"), 
                   by = "WinterYear") %>% 
  dplyr::mutate(prop_pred = TotalPredated / TotalNuts, # proportion predated nuts
                starvation = TotalNuts / TotalNuts_T1, # ratio nuts between T and T1
                ln_starvation = log(starvation + 1), 
                conspec_TotalNuts = (sum_TotalNuts -  TotalNuts)/(n_trees - 1), # total nuts of conspecifics
                pollinated = TotalWhole + TotalPredated, # successfully pollinated nuts (= TotalNuts - TotalEmpty)
                prop_poll = pollinated / TotalNuts, # proportion of pollinated nuts
                TreeID = as.factor(TreeID),
                ord_year = WinterYear - (min(WinterYear) - 1),
                factor_year = as.factor(ord_year)) %>% 
  dplyr::filter(TotalNuts != 0) # filter out all records with zero total nuts as the proportions are then undefined


# 1. Satiation & Starvation effect -----------------------------------------------------

# The satiation effect describes the functional response of the predator, i.e., the more nuts there are in the population 
## the lower the proportion of predated nuts, as the predators become satiated and cannot consume all available nuts.
## The starvation effects describes the numerical response of the predators, that is depended on the change in available nuts between years.
## A high ratio between this year's seed production and last years seed production causes a lower proportion of predated seeds, 
## as high predator populations starve through low food availability in the previous year. High starvation = little nuts in the previous year, many nuts this year


# filter out Inf and NA values for starvation to run models
eos_starv <- eos %>%
  dplyr::filter(ln_starvation != Inf, !is.na(ln_starvation))

## Model: binomial model with temporal auto-correlation

## Try full model as in literature
m_pred <- glmmTMB::glmmTMB(cbind(TotalPredated, TotalNuts - TotalPredated) ~ ord_year : TotalNuts + ord_year : ln_starvation +
                             TotalNuts + ord_year + ln_starvation +
                             I(ln_starvation^2) + I(ord_year^2) + I(TotalNuts^2) +
                             (1 | TreeID) +  ar1(factor_year + 0|TreeID),
                           data = eos_starv, 
                           family = binomial(link = "logit"))

summary(m_pred)


## Drop non-significant term of I(ord_year^2)
m_pred2 <- glmmTMB::glmmTMB(cbind(TotalPredated, TotalNuts - TotalPredated) ~ ord_year : TotalNuts + ord_year : ln_starvation +
                              TotalNuts + ord_year + ln_starvation +
                              I(ln_starvation^2) +  I(TotalNuts^2) +
                              (1 | TreeID) +  ar1(factor_year + 0|TreeID),
                            data = eos_starv,
                            family = binomial(link = "logit"))

summary(m_pred2)
# non-significant quadratic term for TotalNuts is kept in the model as it was fitted on a priori assumptions

# Model diagnostics
glmmTMB::diagnose(m_pred2)
plot(DHARMa::simulateResiduals(m_pred2))
DHARMa::testOutliers(DHARMa::simulateResiduals(m_pred2), type = "bootstrap")
DHARMa::testDispersion(DHARMa::simulateResiduals(m_pred2))


## Predict on final model for both effects #####

# select (ordinal) years for which to predict
years_to_predict <- c(4, 17, 32, 48)

# predict individually for each selected year and only on the range of observed TotalNuts of that year
sat_pred <- purrr::map(.x = years_to_predict,
                       .f = ~{
                         
                         # get the range of values for TotalNuts for year .x
                         nut_values <- eos_starv %>% 
                           filter(ord_year == .x) %>% 
                           pull(TotalNuts)
                         
                         df <- tibble::tibble(TotalNuts = nut_values, 
                                              ord_year = rep(.x, length(TotalNuts))) 
                         
                         # predict on model only for the range of observed TotalNuts in .x
                         predictions <- ggeffects::predict_response(m_pred2,
                                                                    terms = df,
                                                                    type = "fixed") %>% 
                           as.data.frame() 
                         
                         return(predictions)
                         
                       }
) %>% dplyr::bind_rows() 

# reassign calendar years for plotting
pred_m_pred2 <- sat_pred %>% 
  dplyr::mutate(Year = dplyr::case_when(group == 4 ~ 1979, 
                                        group == 17 ~ 1992,
                                        group == 32 ~ 2007,
                                        group == 48 ~ 2023))  

## Figure 2 - panel c: Satiation effect ####
plot_sat <- ggplot2::ggplot(data = pred_m_pred2, ggplot2::aes(x = x, y = predicted)) +
  ggplot2::geom_point(data = eos, ggplot2::aes(x = TotalNuts, y = prop_pred), alpha = 0.1, size = 2) +
  ggplot2::geom_line(ggplot2::aes(colour = as.factor(Year)), linewidth = 1.5) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high, 
                                    colour = as.factor(Year), fill = as.factor(Year)), alpha = 0.1) +
  ggplot2::labs(x = "Total number of beechnuts", 
                y = "Proportion of predated nuts",
                title = "Satiation effect",
                colour = "Year", fill = "Year") +
  ggplot2::scale_colour_manual(values = col_years) + 
  ggplot2::scale_fill_manual(values = col_years) +
  ggplot2::theme_classic(base_size = 17) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 17, hjust = 0.5))


# predict for year effect on proportion predated
pred_pred_year <- ggeffects::predict_response(m_pred2,
                                              terms = c("ord_year [all]"))


# re-assign calendar years
pred_pred_year <- as.data.frame(pred_pred_year) %>% 
  dplyr::left_join(eos %>% 
                     dplyr::distinct(WinterYear, .keep_all = TRUE) %>% 
                     dplyr::select("WinterYear", "ord_year"), 
                   by = c("x" = "ord_year"))

## Figure 2 - panel a: Predation ratio over time ####
plot_starv_year <- ggplot2::ggplot(data = pred_pred_year, ggplot2::aes(x = WinterYear, y = predicted)) +
  ggplot2::geom_point(data = eos, ggplot2::aes(x = WinterYear, y = prop_pred), alpha = 0.1, size = 2) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  ggplot2::geom_line(linewidth = 1.5) +
  ggplot2::scale_x_continuous(breaks = c(seq(1975, 2025, by = 5))) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(x = "Year", 
                y = "Proportion of predated nuts") +
  ggplot2::theme_classic(base_size = 17)


## Starvation effect ----------------------------------------------------

# predict individually for each selected year and only on the range of observed Totalnuts of that year
starv_pred <- purrr::map(.x = years_to_predict,
                         .f = ~{
                           
                           # get only range of observed values for ln_starvation for year .x
                           starv_values <- eos_starv %>% 
                             filter(ord_year == .x) %>% 
                             pull(ln_starvation)
                           
                           df <- tibble::tibble(ln_starvation = starv_values, 
                                                ord_year = rep(.x, length(ln_starvation))) 
                           
                           # predict only on the range of observed values for year .x
                           predictions <- ggeffects::predict_response(m_pred2,
                                                                      terms = df,
                                                                      type = "fixed") %>% 
                             as.data.frame() 
                           
                           return(predictions)
                           
                         }
) %>% dplyr::bind_rows() 


# reassign calendar years for plotting
pred_starv_pred <- starv_pred %>% 
  dplyr::mutate(Year =  dplyr::case_when(group == 4 ~ 1979, 
                                         group == 17 ~ 1992,
                                         group == 32 ~ 2007,
                                         group == 48 ~ 2023))  

# Figure 2 - panel d: Starvation effect ####
plot_starv <- ggplot2::ggplot(data = pred_starv_pred, ggplot2::aes(x = x, y = predicted)) +
  ggplot2::geom_point(data = eos_starv, ggplot2::aes(x = ln_starvation, y = prop_pred), size = 2, alpha = 0.1) +
  ggplot2::geom_line(ggplot2::aes(colour = as.factor(Year)), linewidth = 1.5) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high, 
                                    colour = as.factor(Year), fill = as.factor(Year)), alpha = 0.1) +
  ggplot2::labs(x = "Ln Seed production ratio (T/T-1)", 
                y = "Proportion of predated nuts", 
                title = "Starvation effect",
                colour = "Year", fill = "Year") +
  ggplot2::scale_colour_manual(values = col_years) + 
  ggplot2::scale_fill_manual(values = col_years) +
  ggplot2::theme_classic(base_size = 17) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 17, hjust = 0.5))



# 2. Pollination efficiency -----------------------------------------------

## The pollination efficiency depends on the synchrony of flowering events in the population as the 
## probability of successful pollination increases the more pollen there is. 

## Model of the literature
m_poll <- glmmTMB::glmmTMB(cbind(pollinated, TotalNuts - pollinated) ~ conspec_TotalNuts * CVp  + ord_year +
                             I(conspec_TotalNuts ^ 2) + I(ord_year^2) +
                             (1 | TreeID) +  ar1(factor_year + 0|TreeID),
                           data = eos,
                           family = binomial(link = "logit"))

summary(m_poll)

## drop non-significant quadratic year term
m_poll2 <- glmmTMB::glmmTMB(cbind(pollinated, TotalNuts - pollinated) ~ conspec_TotalNuts * CVp  + ord_year + 
                              I(conspec_TotalNuts ^ 2) +
                              (1 | TreeID) +  ar1(factor_year + 0|TreeID),
                            data = eos,
                            family = binomial(link = "logit"))

summary(m_poll2)

# model diagnostics
glmmTMB::diagnose(m_poll2)
plot(DHARMa::simulateResiduals(m_poll2))

# look for correlation between CVp and nuts of conspecifics
ggplot2::ggplot(eos, ggplot2::aes(x = CVp, y = conspec_TotalNuts)) +
  ggplot2::geom_point()


## Predict on final model ####

# get model predictions on the data scale for plotting
pred_poll <- ggeffects::predict_response(m_poll2,
                                         terms = c("conspec_TotalNuts [n=30]", "CVp [0.5, 2, 3.5, 5]"))

## Figure 2 - panel e: Pollination efficiency ####
plot_poll <- ggplot2::ggplot(data = pred_poll, ggplot2::aes(x = x, y = predicted)) +
  ggplot2::geom_point(data = eos, ggplot2::aes(x = conspec_TotalNuts, y = prop_poll), size = 2, alpha = 0.1) +
  ggplot2::geom_line(linewidth = 1.5, ggplot2::aes(colour = group)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high, colour = group, fill = group), alpha = 0.2) +
  ggplot2::labs(x = "Number nuts of conspecifics", 
                y = "Proportion of pollinated nuts", 
                title = "Pollination efficiency",
                colour = "CVp", fill = "CVp") +
  ggplot2::scale_colour_manual(values = col_CVp) +
  ggplot2::scale_fill_manual(values = col_CVp) +
  ggplot2::theme_classic(base_size = 17) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 17, hjust = 0.5))


# predict for the year effect
pred_poll_time <- ggeffects::predict_response(m_poll2,
                                              terms = c("ord_year [1:49]"))

# re-assign calendar years
pred_poll_time <- as.data.frame(pred_poll_time) %>% 
  dplyr::left_join(eos %>% 
                     dplyr::distinct(WinterYear, .keep_all = TRUE) %>% 
                     dplyr::select("WinterYear", "ord_year"), 
                   by = c("x" = "ord_year"))

## Figure 2 - panel b: Porportion pollinated over time ####
plot_poll_year <- ggplot2::ggplot(data = pred_poll_time, ggplot2::aes(x = WinterYear, y = predicted)) +
  ggplot2::geom_point(data = eos, ggplot2::aes(x = WinterYear, y = prop_poll), alpha = 0.1, size = 2) +
  ggplot2::geom_line(linewidth = 2) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  ggplot2::labs(x = "Year", 
                y = "Proportion of pollinated nuts") +
  ggplot2::scale_x_continuous(breaks = c(seq(1975, 2025, by = 5))) +
  ggplot2::theme_classic(base_size = 17)


# Make final figure 2 -------------------------------------------------------

# aligning the panels with different legends requires some tweaking

# Give all plots the same outer margins to avoid small offsets
common_margin <- ggplot2::theme(plot.margin = grid::unit(c(4, 4, 4, 4), "pt"))
plot_sat <- plot_sat + common_margin
plot_starv <- plot_starv + common_margin
plot_poll <- plot_poll + common_margin

# Extract legends (set legend.position for extraction)
leg_pred_grob <- cowplot::get_legend(plot_sat + ggplot2::theme(legend.position = "bottom"))  
leg_poll__grob  <- cowplot::get_legend(plot_poll + ggplot2::theme(legend.position = "bottom"))  

# Remove legends from the plots themselves
plot_sat_noleg <- plot_sat + ggplot2::theme(legend.position = "none")
plot_starv_noleg <- plot_starv + ggplot2::theme(legend.position = "none")
plot_poll_noleg <- plot_poll + ggplot2::theme(legend.position = "none")

# Convert legend grobs into ggplot objects so ggarrange can place them
leg_pred_plot <- ggpubr::as_ggplot(leg_pred_grob)
leg_poll_plot  <- ggpubr::as_ggplot(leg_poll__grob)

# make lower half of Figure by arranging EOS panels
EOS_plots <- ggpubr::ggarrange(plot_sat_noleg, plot_starv_noleg, plot_poll_noleg,
                               ncol = 3, nrow = 1, labels = c("c", "d", "e"))  

# Put the two extracted legends side-by-side in their own row
legend_row <- ggpubr::ggarrange(leg_pred_plot, leg_poll_plot, ncol = 2, widths = c(1, 0.4))

# Make Figure 2
ggpubr::ggarrange(
  
  ggpubr::ggarrange(plot_starv_year, plot_poll_year, labels = c("a", "b")),
  
  EOS_plots, legend_row, 
  
  ncol = 1, nrow = 3, heights = c(0.8, 1, 0.12)
)
 

# save Figure 2 as PNG
ggplot2::ggsave(plot = ggplot2::last_plot(), file = here::here("plots", "Figure_2.png"), 
               units = "cm", height = 25, width = 40, dpi = 600)
