################################################################################################
# Appendix 1: Sensitivity Analysis of Calibrations
# (Requires prior run of model_evaluation_analysis.R)
# 
# Code evaluating calibration plots for 1, 2, 3, and 4 year timepoints.
# Primary plots with groupkm. val.surv used to get fit estimates using maxdim=2.
#
# Written By: Kyle W. Singleton
################################################################################################

### Test Calibration at multiple values of time
# - Requires baseline Survival, S(t), from KM at each timepoint
# - Nearest S(t) value extracted from WebPlotDigitizer output

# Load WebPlotDigitizer Data
wd <- getwd()
data_dir <- paste0(wd, "/data/KM_Survival")
setwd(data_dir)

file_helseth_km <- "helseth_extracted_survival_5Xby5Y.csv"
file_michaelsen_km <- "michaelsen_extracted_survival_5Xby5Y.csv"
file_gutman_km <- "gutman_reestimated_extracted_survival_5Xby5Y.csv"
file_kumar_km <- "kumar_extracted_survival_5Xby5Y.csv"

ucla_helseth_km <- read.csv(file=file_helseth_km, head=TRUE, sep=",")
ucla_michaelsen_km <- read.csv(file=file_michaelsen_km, head=TRUE, sep=",")
ucla_gutman_km <- read.csv(file=file_gutman_km, head=TRUE, sep=",")
ucla_kumar_km <- read.csv(file=file_kumar_km, head=TRUE, sep=",")

setwd(wd)

# Find nearest timepoint and extract S(t)
months_of_interest = c(12, 24, 36, 48)
days_of_interest = c(365, 365*2, 365*3, 365*4)

extracted_S_T <- data.frame(Days = days_of_interest)

# Reported KM's in Months
for (i in seq(1,length(months_of_interest))){
  month = months_of_interest[i]
  extracted_S_T$Helseth[i] = ucla_helseth_km$surv[which(abs(ucla_helseth_km$time-month)==min(abs(ucla_helseth_km$time-month)))]
  extracted_S_T$Michaelsen[i] = ucla_michaelsen_km$surv[which(abs(ucla_michaelsen_km$time-month)==min(abs(ucla_michaelsen_km$time-month)))]
  extracted_S_T$Kumar[i] = ucla_kumar_km$surv[which(abs(ucla_kumar_km$time-month)==min(abs(ucla_kumar_km$time-month)))]
}

# Re-estimated Gutman in Days
for (i in seq(1,length(months_of_interest))){
  days = days_of_interest[i]
  extracted_S_T$Gutman[i] = ucla_gutman_km$surv[which(abs(ucla_gutman_km$time-days)==min(abs(ucla_gutman_km$time-days)))]
}

# Calculate Survival Probabilities for each timepoint and plot 
for(i in seq(1, length(extracted_S_T$Days))){
  u_Days = extracted_S_T$Days[i]
  plot_valsurv = FALSE #Set to TRUE to view val.surv plots with groupkm
  
  ucla_helseth$prob_temp <- (extracted_S_T$Helseth[i]^exp(ucla_helseth$total_corr))
  ucla_michaelsen$prob_temp <- (extracted_S_T$Michaelsen[i]^exp(ucla_michaelsen$total_corr))
  ucla_gutman$prob_temp <- (extracted_S_T$Gutman[i]^exp(ucla_gutman$total_corr))
  ucla_kumar$prob_temp <- (extracted_S_T$Kumar[i]^exp(ucla_kumar$total_corr))
  
  set_max = 2 # Fixed to 2 (except for Kumar - set to 3) to get simplest HARE fit
  w_h_temp <- val.surv(est.surv = ucla_helseth$prob_temp, S=Surv(ucla_helseth$SurvivalDays, ucla_helseth$OutcomeDead), u=u_Days, fun=function(p)log(-log(p)), maxdim=set_max)
  w_m_temp <- val.surv(est.surv = ucla_michaelsen$prob_temp, S=Surv(ucla_michaelsen$SurvivalDays, ucla_michaelsen$OutcomeDead), u=u_Days, fun=function(p)log(-log(p)), maxdim=set_max)
  w_g_temp <- val.surv(est.surv = ucla_gutman$prob_temp, S=Surv(ucla_gutman$SurvivalDays, ucla_gutman$OutcomeDead), u=u_Days, fun=function(p)log(-log(p)), maxdim=set_max)
  w_k_temp <- val.surv(est.surv = ucla_kumar$prob_temp, S=Surv(ucla_kumar$SurvivalDays, ucla_kumar$OutcomeDead), u=u_Days, fun=function(p)log(-log(p)), maxdim=set_max+1)
  
  par(mfrow=c(2,2))
  helseth_cal <- groupkm(ucla_helseth$prob_temp, Surv(ucla_helseth$SurvivalDays, ucla_helseth$OutcomeDead),
                         g=5, u=u_Days, pl=TRUE, xlim=c(0,1), ylim=c(0,1), cex.subtitle = FALSE,
                         xlab="Observed Probability\nN = 125, 25 patients per group", 
                         ylab="Predicted Probability")
  lines(x = c(0,1), y = c(0,1), lty=5, col="gray")
  title(paste0("Helseth", " (", u_Days, " Days)"))
  if (plot_valsurv){
    par(new="TRUE")
    plot(w_h_temp, xlab="", ylab="", lim=c(0,1))    
  }
  
  michaelsen_cal <- groupkm(ucla_michaelsen$prob_temp, Surv(ucla_michaelsen$SurvivalDays, ucla_michaelsen$OutcomeDead),
                            g=5, u=u_Days, pl=TRUE, xlim=c(0,1), ylim=c(0,1), cex.subtitle = FALSE, 
                            xlab="Observed Probability\nN = 125, 25 patients per group", 
                            ylab="Predicted Probability")
  lines(x = c(0,1), y = c(0,1), lty=5, col="gray")
  title(paste0("Michaelsen", " (", u_Days, " Days)"))
  if (plot_valsurv){
    par(new="TRUE")
    plot(w_m_temp, xlab="", ylab="", lim=c(0,1))
  }
  
  gutman_cal <- groupkm(ucla_gutman$prob_temp, Surv(ucla_gutman$SurvivalDays, ucla_gutman$OutcomeDead),
                        g=5, u=u_Days, pl=TRUE, xlim=c(0,1), ylim=c(0,1), cex.subtitle = FALSE,
                        xlab="Observed Probability\nN = 125, 25 patients per group", 
                        ylab="Predicted Probability")
  lines(x = c(0,1), y = c(0,1), lty=5, col="gray")
  title(paste0("Gutman", " (", u_Days, " Days)"))
  if (plot_valsurv){
    par(new="TRUE")
    plot(w_g_temp, xlab="", ylab="", lim=c(0,1))
  }
  
  kumar_cal <- groupkm(ucla_kumar$prob_temp, Surv(ucla_kumar$SurvivalDays, ucla_kumar$OutcomeDead),
                       g=20, u=u_Days, pl=TRUE, xlim=c(0,1), ylim=c(0,1), cex.subtitle = FALSE,
                       xlab="Observed Probability\nN = 125, 25 patients per group", 
                       ylab="Predicted Probability")
  lines(x = c(0,1), y = c(0,1), lty=5, col="gray")
  title(paste0("Kumar", " (", u_Days, " Days)"))
  if (plot_valsurv){
    par(new="TRUE")
    plot(w_k_temp, xlab="", ylab="", lim=c(0,1))
  }
  par(mfrow=c(1,1))
  
  # Print HARE fit coefficients
  message(paste0("Evaluation: ", u_Days, " Days"))
  message("Helseth")
  print(w_h_temp$harefit$fcts)
  message("Michaelsen")
  print(w_m_temp$harefit$fcts)
  message("Gutman")
  print(w_g_temp$harefit$fcts)
  message("Kumar")
  print(w_k_temp$harefit$fcts)
  
  # Include Median calibration from primary analysis
  if (i == 1){
    u_Days = 506 # UCLA Median
    
    w_h <- val.surv(est.surv = ucla_helseth$prob_ucla, S=Surv(ucla_helseth$SurvivalDays, ucla_helseth$OutcomeDead), u=u_Days, fun=function(p)log(-log(p)), maxdim=set_max)
    w_m <- val.surv(est.surv = ucla_michaelsen$prob_ucla, S=Surv(ucla_michaelsen$SurvivalDays, ucla_michaelsen$OutcomeDead), u=u_Days, fun=function(p)log(-log(p)), maxdim=set_max)
    w_g <- val.surv(est.surv = ucla_gutman$prob_ucla, S=Surv(ucla_gutman$SurvivalDays, ucla_gutman$OutcomeDead), u=u_Days, fun=function(p)log(-log(p)), maxdim=set_max)
    w_k <- val.surv(est.surv = ucla_kumar$prob_ucla, S=Surv(ucla_kumar$SurvivalDays, ucla_kumar$OutcomeDead), u=u_Days, fun=function(p)log(-log(p)), maxdim=set_max+1)
    
    par(mfrow=c(2,2))
    helseth_cal <- groupkm(ucla_helseth$prob_ucla, Surv(ucla_helseth$SurvivalDays, ucla_helseth$OutcomeDead),
                           g=5, u=u_Days, pl=TRUE, xlim=c(0,1), ylim=c(0,1), cex.subtitle = FALSE,
                           xlab="Observed Probability\nN = 125, 25 patients per group", 
                           ylab="Predicted Probability")
    lines(x = c(0,1), y = c(0,1), lty=5, col="gray")
    title(paste0("Helseth", " (", u_Days, " Days)"))
    if (plot_valsurv){
      par(new="TRUE")
      plot(w_h, xlab="", ylab="", lim=c(0,1))
    }
    
    
    michaelsen_cal <- groupkm(ucla_michaelsen$prob_ucla, Surv(ucla_michaelsen$SurvivalDays, ucla_michaelsen$OutcomeDead),
                              g=5, u=u_Days, pl=TRUE, xlim=c(0,1), ylim=c(0,1), cex.subtitle = FALSE, 
                              xlab="Observed Probability\nN = 125, 25 patients per group", 
                              ylab="Predicted Probability")
    lines(x = c(0,1), y = c(0,1), lty=5, col="gray")
    title(paste0("Michaelsen", " (", u_Days, " Days)"))
    if (plot_valsurv){
      par(new="TRUE")
      plot(w_m, xlab="", ylab="", lim=c(0,1))
    }
    
    gutman_cal <- groupkm(ucla_gutman$prob_ucla, Surv(ucla_gutman$SurvivalDays, ucla_gutman$OutcomeDead),
                          g=5, u=u_Days, pl=TRUE, xlim=c(0,1), ylim=c(0,1), cex.subtitle = FALSE,
                          xlab="Observed Probability\nN = 125, 25 patients per group", 
                          ylab="Predicted Probability")
    lines(x = c(0,1), y = c(0,1), lty=5, col="gray")
    title(paste0("Gutman", " (", u_Days, " Days)"))
    if (plot_valsurv){
      par(new="TRUE")
      plot(w_g, xlab="", ylab="", lim=c(0,1))
    }
    
    kumar_cal <- groupkm(ucla_kumar$prob_ucla, Surv(ucla_kumar$SurvivalDays, ucla_kumar$OutcomeDead),
                         g=20, u=u_Days, pl=TRUE, xlim=c(0,1), ylim=c(0,1), cex.subtitle = FALSE,
                         xlab="Observed Probability\nN = 125, 25 patients per group", 
                         ylab="Predicted Probability")
    lines(x = c(0,1), y = c(0,1), lty=5, col="gray")
    title(paste0("Kumar", " (", u_Days, " Days)"))
    if (plot_valsurv){
      par(new="TRUE")
      plot(w_k, xlab="", ylab="", lim=c(0,1))
    }
    par(mfrow=c(1,1))
    
    # Print HARE fit coefficients
    message(paste0("Evaluation: ", u_Days, " Days"))
    message("Helseth")
    print(w_h$harefit$fcts)
    message("Michaelsen")
    print(w_m$harefit$fcts)
    message("Gutman")
    print(w_g$harefit$fcts)
    message("Kumar")
    print(w_k$harefit$fcts)
  }
  
}

# Clean up data frames
ucla_helseth <- subset(ucla_helseth, select = -c(prob_temp))
ucla_michaelsen <- subset(ucla_michaelsen, select = -c(prob_temp))
ucla_gutman <- subset(ucla_gutman, select = -c(prob_temp))
ucla_kumar <- subset(ucla_kumar, select = -c(prob_temp))