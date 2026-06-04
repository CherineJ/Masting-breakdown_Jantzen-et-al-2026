# README

This repository contains all code and data used for the analyses in Jantzen et al. (2026) (DOI: [10.1002/ece3.73809](https://doi.org/10.1002/ece3.73809)). 
Please be aware that the data in this repository is a derivative, as used for analyses, of the full, openly accessible dataset on beechnut production in the Netherlands that can be found on DataverseNL under the DOI: [10.34894/TQY74M](https://doi.org/10.34894/TQY74M).


## How to use

The renv lock file allows to reproduce the packages adn their versions as used for the original analyses. Scripts should be run within this R project, as they refer to relative paths within the project to read and write data and plots (by using the here package).

## Structure

All R scripts are stored in the folder "R". It consists of six main scripts used for different components of the analysis and three helper scripts, containing colour set up and functions that are used within the other scripts (see below for more details). The folder "plots" contains all figures of the main text and the supplements (marked through the letter "S"). The "data" folder contains all input data files, as well as data files generated throughout the analyses. This allows to run scripts independent from each other.

### Content of each script

-   *01-data-preparation.R*: This script prepares the climate data for analysis, calculating climate values for selected periods of the year and combining this with the seed production data, ultimately creating the final input file (di_input-data.csv) used for all further analyses.

-   *02-analysis-masting-metrics.R*: This script analyses the temporal change of seed production and the masting characteristics (synchrony and inter-annual variation).

-   *03-analysis_economy-of-scale.R*: This script looks at the fitness benefits of masting, analysing changes in pollination and predation and the economies of scale (predator satiation, starvation and pollination efficiency).

-   *04-analysis_sens-windows.R*: This script finds the windows of highest sensitivity of seed production to different climate variables.

-   *05-analysis_climate-drivers-seed-production.R*: This script processes the results of script 04, analyses changes in the identified weather cues and analyses the effects of the weather cues and resource reserves on seed production.

-   *06-supplementary_sens-window-Journé-et-al-2024.R*: This script contains supplementary analysis to compare the results of the our approach of identifying weather cues with the approach reported in Journé et al. (2024).

-   *colour_palette.R*: contains the colour palette used throughout scripts for plotting

-   *function_calc-mean-climate-in-window.R*: This script contains a function that is used to more easily calculate the mean or sum of a climate variable in a defined period of the year.

-   *functions_sensitivity-window-analyses.R*: This script contains three functions used within script 04 to facilitate the calculation and visualisation of the sensitivity windows.
