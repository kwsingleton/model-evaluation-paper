# Replicating Published Survival Prediction Models for Local Applications

This repository accompanies the paper "Replicating Published Survival Prediction Models for Local Applications" and provides the data and R code to repeat the reported results and figures.

### R Scripts

`model_replication_analysis.R`

Applys four published Cox survival models to local data by: 
 - Computing linear predictors and survival probabilies
 - Generating ROC curves with associated AUC values
 - Plotting calibration curves

`data_preprocessing.R`

Loads data files and handles required data manipulations to prepare for analysis. Included automatically when running primary model_replication_analysis script.

### Data

`Helseth_CoxData.csv, Michaelsen_CoxData.csv, Gutman_CoxData.csv, Kumar_CoxData.csv`

Local data from UCLA are split into four files with relevant features for the four comparison published models .

`Gutman2013_SurvivalData.csv`

Subset of relevant features from shared data by [Gutman et.al. 2013](https://doi.org/10.1148/radiol.13120118). Used to determine Gutman model feature means and median survival time.

### R Packages Used
 - `rms` - https://cran.r-project.org/web/packages/rms/
 - `survivalROC` - https://cran.r-project.org/web/packages/survivalROC/
