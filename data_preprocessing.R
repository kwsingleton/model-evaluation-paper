################################################################################################
# Data Preprocessing
# 
# Imports, converts, and adjusts CSV files of local UCLA patient data and shared data from
# Gutman et.al. 2013 (https://doi.org/10.1148/radiol.13120118).
#
# Written By: Kyle W. Singleton
################################################################################################

wd <- getwd()
data_dir <- paste0(wd, "/data/")
setwd(data_dir)

## Import pre-processed UCLA data
file_helseth <- "Helseth_CoxData.csv"
file_michaelsen <- "Michaelsen_CoxData.csv"
file_gutman <- "Gutman_CoxData.csv"
file_kumar <- "Kumar_CoxData.csv"

ucla_helseth <- read.csv(file=file_helseth, head=TRUE, sep=",") #125 cases
ucla_michaelsen <- read.csv(file=file_michaelsen, head=TRUE, sep=",") #125 cases
ucla_gutman <- read.csv(file=file_gutman, head=TRUE, sep=",") #125 cases
ucla_kumar <- read.csv(file=file_kumar, head=TRUE, sep=",") #125 cases

## Correct value difference between UCLA and Gutman size measurements (cm vs mm)
ucla_gutman$F29 <- ucla_gutman$F29 * 10

## Import pre-processed Gutman 2013 data (SUpplied by Gutman, based on TCGA Dataset)
file_gutman_original <- "Gutman2013_SurvivalData.csv"
gutman_data <- read.csv(file=file_gutman_original, head=TRUE, sep=",") #68 cases

# Recode Censorship Values: Death = 1, Alive = 0
ucla_helseth <- cbind(ucla_helseth[1], OutcomeDead = model.matrix(~0+Outcome, ucla_helseth)[,2], ucla_helseth[-1])
ucla_michaelsen <- cbind(ucla_michaelsen[1], OutcomeDead = model.matrix(~0+Outcome, ucla_michaelsen)[,2], ucla_michaelsen[-1])
ucla_gutman <- cbind(ucla_gutman[1], OutcomeDead = model.matrix(~0+Outcome, ucla_gutman)[,2], ucla_gutman[-1])
ucla_kumar <- cbind(ucla_kumar[1], OutcomeDead = model.matrix(~0+Outcome, ucla_kumar)[,2], ucla_kumar[-1])

gutman_data <- cbind(gutman_data[1], OutcomeDead = model.matrix(~0+VITALSTATUS, gutman_data)[,2], gutman_data[-1])

######### Cleanup
rm(data_dir)
rm(file_helseth)
rm(file_michaelsen)
rm(file_gutman)
rm(file_kumar)
rm(file_gutman_original)
setwd(wd)