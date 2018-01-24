################################################################################################
# Model Evaluation Analysis
# 
# Code applying four published Survival models to local data from UCLA.
#
# Written By: Kyle W. Singleton
################################################################################################

### Step 0
### Preload libraries (Installs missing libraries if required)
list.of.packages <- c("survival","rms", "survivalROC", "pROC")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(rms)
library(survivalROC)
library(pROC)

### Step 1
### Load pre-processed UCLA and shared Gutman data
source(paste0(getwd(),'/data_preprocessing.R'))

### Step 2
### Gather important published values

# Published Hazard Coefficients
coeffs <- list(helseth=log(c(1.44, 1.02, 2.13, 2.72, 2.31))) ## Sex, Age, ECOG, Surgery, TumorLocation
coeffs$michaelsen <- log(c(1.22, 2.06, 2.06, 1.31)) ## ECOG, ECOG 2v0, Corticosteroids, Age by 10
coeffs$gutman <- log(c(0.972, 1.016, 7.745)) ## KPS, Major Axis Length, Proportion CET
coeffs$kumar <- log(c(2.34, 1.52, 2.03, 0.44)) ## Tumor Site, Tumor Location, Radiation Dose, Chemo Given

# Published Feature means
# Note: Gutman values not summarized in paper, derived from shared data in gutman_data dataframe
means <- list(helseth=c(0.411, 63.7, 0.188, 0.089, 0.147))
means$michaelsen <- c(0.2933333, 0.08444444, 0.733, 6)
means$gutman <- c(mean(gutman_data$KPS), mean(gutman_data$LesionSizeMajor), mean(gutman_data$PropEnhance))
means$kumar <- c(0.044872, 0.564103, 0.224359, 0.2852564)

# Local data (UCLA) median survival time - 506 Days for all local data
ucla_median_survival = median(ucla_gutman$SurvivalDays)

# Survival probabilities at t=median survival
# Determined using graph digitizer from published Kaplan-Meier curves
s_t <- list(ucla_median=c(0.257206, 0.406393, 0.367255, 0.159468))

# Calculate values for average patient (i.e., all values set to feature means)
correction <- list(helseth = sum(coeffs$helseth * means$helseth))
correction$michaelsen <- sum(coeffs$michaelsen * means$michaelsen)
correction$gutman <- sum(coeffs$gutman * means$gutman)
correction$kumar <- sum(coeffs$kumar * means$kumar)

### Step 3
### Calculate local data (UCLA) linear predictor and survival probability values

# Helseth
temp <- ucla_helseth[,4:8]
temp <- data.frame(mapply(`*`,temp,coeffs$helseth))
ucla_helseth$total <- apply(temp, 1, sum)
ucla_helseth$total_corr <- ucla_helseth$total - correction$helseth
ucla_helseth$prob_ucla <- (s_t$ucla_median[1]^exp(ucla_helseth$total_corr))

# Michaelsen
temp <- ucla_michaelsen[,4:7]
temp <- data.frame(mapply(`*`,temp,coeffs$michaelsen))
ucla_michaelsen$total <- apply(temp, 1, sum)
ucla_michaelsen$total_corr <- ucla_michaelsen$total - correction$michaelsen
ucla_michaelsen$prob_ucla <- (s_t$ucla_median[2]^exp(ucla_michaelsen$total_corr))

# Gutman
temp <- ucla_gutman[,4:6]
temp <- data.frame(mapply(`*`,temp,coeffs$gutman))
ucla_gutman$total <- apply(temp, 1, sum)
ucla_gutman$total_corr <- ucla_gutman$total - correction$gutman
ucla_gutman$prob_ucla <- (s_t$ucla_median[3]^exp(ucla_gutman$total_corr))

# Kumar
temp <- ucla_kumar[,4:7]
temp <- data.frame(mapply(`*`,temp,coeffs$kumar))
ucla_kumar$total <- apply(temp, 1, sum)
ucla_kumar$total_corr <- ucla_kumar$total - correction$kumar
ucla_kumar$prob_ucla <- (s_t$ucla_median[4]^exp(ucla_kumar$total_corr))

### Step 4
### Calculate ROC curves and AUC values

# ROC Curves
hROC <- survivalROC(Stime=ucla_helseth$SurvivalDays, status=ucla_helseth$OutcomeDead, 
                    marker = -ucla_helseth$prob_ucla, predict.time = ucla_median_survival, method="KM")
mROC <- survivalROC(Stime=ucla_michaelsen$SurvivalDays, status=ucla_michaelsen$OutcomeDead, 
                    marker = -ucla_michaelsen$prob_ucla, predict.time = ucla_median_survival, method="KM")
gROC <- survivalROC(Stime=ucla_gutman$SurvivalDays, status=ucla_gutman$OutcomeDead, 
                    marker = -ucla_gutman$prob_ucla, predict.time = ucla_median_survival, method="KM")
kROC <- survivalROC(Stime=ucla_kumar$SurvivalDays, status=ucla_kumar$OutcomeDead, 
                    marker = -ucla_kumar$prob_ucla, predict.time = ucla_median_survival, method="KM")

plot(hROC$FP, hROC$TP, type="l", col="red", lty=3, xlab="", ylab="")
lines(mROC$FP, mROC$TP, type="l", col="black", lty=1, xlab="", ylab="")
lines(gROC$FP, gROC$TP, type="l", col="green", lty=4, xlab="", ylab="")
lines(kROC$FP, kROC$TP, type="l", col="orange", lty=5, xlab="", ylab="")
abline(a=0, b= 1, lty=2, col="gray")
title(main = "ROC Curve Comparisons", xlab="False Positive Rate", ylab="True Positive Rate")
legend("topleft", legend=c("Helseth", "Michaelsen", "Gutman", "Kumar"), col=c("red","black","green","orange"), lty=c(3,1,4,5), cex=0.8)

# Collect AUC values and print
auc_survROC <- list()
auc_survROC$helseth <- hROC$AUC
auc_survROC$michaelsen <- mROC$AUC
auc_survROC$gutman <- gROC$AUC
auc_survROC$kumar <- kROC$AUC
unlist(auc_survROC)

### Step 5
### Calibration Curves
par(mfrow=c(2,2))
helseth_cal <- groupkm(ucla_helseth$prob_ucla, Surv(ucla_helseth$SurvivalDays, ucla_helseth$OutcomeDead),
        g=5, u=506, pl=TRUE, xlim=c(0,1), ylim=c(0,1), cex.subtitle = FALSE,
        xlab="Observed Probability\nN = 125, 25 patients per group", 
        ylab="Predicted Probability")
lines(x = c(0,1), y = c(0,1), lty=5, col="gray")
title("Helseth")

michaelsen_cal <- groupkm(ucla_michaelsen$prob_ucla, Surv(ucla_michaelsen$SurvivalDays, ucla_michaelsen$OutcomeDead),
        g=5, u=506, pl=TRUE, xlim=c(0,1), ylim=c(0,1), cex.subtitle = FALSE, 
        xlab="Observed Probability\nN = 125, 25 patients per group", 
        ylab="Predicted Probability")
lines(x = c(0,1), y = c(0,1), lty=5, col="gray")
title("Michaelsen")

gutman_cal <- groupkm(ucla_gutman$prob_ucla, Surv(ucla_gutman$SurvivalDays, ucla_gutman$OutcomeDead),
        g=5, u=506, pl=TRUE, xlim=c(0,1), ylim=c(0,1), cex.subtitle = FALSE,
        xlab="Observed Probability\nN = 125, 25 patients per group", 
        ylab="Predicted Probability")
lines(x = c(0,1), y = c(0,1), lty=5, col="gray")
title("Gutman")

kumar_cal <- groupkm(ucla_kumar$prob_ucla, Surv(ucla_kumar$SurvivalDays, ucla_kumar$OutcomeDead),
        g=20, u=506, pl=TRUE, xlim=c(0,1), ylim=c(0,1), cex.subtitle = FALSE,
        xlab="Observed Probability\nN = 125, 25 patients per group", 
        ylab="Predicted Probability")
lines(x = c(0,1), y = c(0,1), lty=5, col="gray")
title("Kumar")
par(mfrow=c(1,1))

### Step 6
### Membership Model and Linear Predictor comparisons from framework 
### by Debray et.al. (https://doi.org/10.1016/j.jclinepi.2014.06.018)

### Approach 1: Distinguishing between individuals of validation and development sets
# Copy data harmonize, and combine two datasets
gutman_debray = gutman_data[2:6]
colnames(gutman_debray) <- c("OutcomeDead", "SurvivalDays", "KPS", "F29", "F5")
ucla_debray = ucla_gutman[2:6]
gutman_debray$Sample = 1
ucla_debray$Sample = 0
combined_debray = rbind(gutman_debray, ucla_debray)

# Create Membership Model with Logistic Regression
membership_model <- glm(Sample ~ ., family=binomial(link='logit'), data=combined_debray)
exp(cbind("Odds ratio" = coef(membership_model), confint.default(membership_model, level = 0.95)))

# Get AUC for prediction of sample set
combined_debray$lp <- predict(membership_model, combined_debray, type="response")

# Compute standard ROC
roc_curve <- roc(combined_debray$Sample ~ combined_debray$lp, ci=TRUE)
plot(roc_curve, print.auc=TRUE, main="Membership Model ROC Curve")
auc(roc_curve)

### Approach 2: Comparing the predicted risks between development and validation samples

# Generate linear predictors for Original Gutman cases
temp <- gutman_data[,4:6]
temp <- data.frame(mapply(`*`,temp,coeffs$gutman))
gutman_data$total <- apply(temp, 1, sum)
gutman_data$total_corr <- gutman_data$total - correction$gutman
gutman_data$prob_ucla <- (s_t$ucla_median[1]^exp(gutman_data$total_corr))

# Compare Difference of LP Means
u_dev = mean(gutman_data$total)
u_val = mean(ucla_gutman$total)
message("Difference in Means")
u_val - u_dev

# Compare Ratio of LP Standard Deviations
sd_dev = sd(gutman_data$total)
sd_val = sd(ucla_gutman$total)
message("Ratio of Deviations")
sd_val/sd_dev
