## Colour palette plots ####

# individual trees
col_per_tree <- "#008171"

# population mean (over all trees)
col_pop_mean <- "black"

# temporal trend in seed production split into conditional & zero-inflated model
nut_trend_split <- "#6A0136"

# predicted seed production under observed climate
col_pred.real_clim<- "#A06CD5"  

# colours climate variables
col_pal <- c("PrecT2" = "#e69f00", 
             "PrecT1" =  "#009e73", 
             "PrecT0" = "#cc79a7", 
             "maxTempT2" = "#56b4e9",
             "maxTempT1" = "#d55e00", 
             "meanTempT0" = "#0072b2")

# shapes climate variables
shape_pal <- c("PrecT2" = 15, 
             "PrecT1" = 16, 
             "PrecT0" = 17, 
             "maxTempT2" = 15,
             "maxTempT1" = 16, 
             "meanTempT0" = 17)

# years for which we predict EOS effects 
col_years <- c("1979" = "#0094ff",
               "1992" = "#34AB53", 
               "2007" = "#e8780f", 
               "2023" = "#54428E") 

# levels of CVp (= inverse of synchrony) for which we predict pollination ratios
col_CVp <- c("0.5" =  "#184e77",
             "2" =  "#168aad",
             "3.5" =  "#52b69a",
             "5" =  "#99d98c")
