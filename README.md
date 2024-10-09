# Code and data for the statistical analysis in "Rapid temporal adaptation structures tolerance to toxic cyanobacteria in a natural population of the water flea Daphnia."

This repository contains the code and data for the manuscript "Rapid temporal adaptation structures tolerance to toxic cyanobacteria in a natural population of the water flea Daphnia."

The R-file [`Main_script.R`](./Main_script.R) constitutes the principal script for reading, processing, visualizing, and analyzing the experimental data. The csv-file [`summary_data.csv`](./summary_data.csv) contains the survival data for each of the time intervals summarized per experimental unit and is read in the main script. The Stan-files [`Binomial_GLMM.stan`](./Binomial_GLMM.stan) and [`Interval_censored_survival.stan`](./Interval_censored_survival.stan) contain the probabilistic models used for the Bayesian analyses in the main script.

The R-file [`Seasonal_dynamics_visualization.R`](./Seasonal_dynamics_visualization.R) constitutes a supplementary script for reading and visualizing observational data on temporal dynamics of phytoplankton and zooplankton in Langerodevijver from February to November 2019, based on the data in the csv-file [`Seasonal_dynamics_data.csv`](./Seasonal_dynamics_data.csv).
