Before running the bayesian_model.R script, replace "~/Desktop/folder-location" with the folder location of bayesian_model.R

Introduction: The bayesian_model.R script runs the predictive model of resistant infection counts using different metagenomic covariates. It assesses the model fits using different information using leave-one-out cross-validation and visual inspection. 


In particular, the script

1) loads clinical and metagenomic information
2) prepares the covariates (imputation and normalisation)
3) runs the the hierarchical model fits using rstan and the Stan models located in the stan folder
4) performs the model comparison
5) visualises the model fit analogous to Figure 5 in the paper  

This script was tested with R version 3.6.1 (2019-07-05), and requires the following libraries: 
-tidyverse
-rstan
-bayesplot
-loo
-bayestestR