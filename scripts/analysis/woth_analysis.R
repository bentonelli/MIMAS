### Script to look at paramter estimates for WOTH over time ####

library(dplyr)
library(ggplot2)
#Read in prior ranges
M1_Pr <- lapply(Sys.glob("data/training_runs/WOTH_M1/3_28_22/Step_1/WOTH*.csv"), read.csv)
M1_Pr <- do.call(bind_rows,M1_Pr[1:1000])

M2_Pr <- lapply(Sys.glob("data/training_runs/WOTH_M2/3_30_22/Step_1/WOTH*.csv"), read.csv)
M2_Pr <- do.call(bind_rows,M2_Pr[1:1000])

M3_Pr <- lapply(Sys.glob("data/training_runs/WOTH_M3/3_31_22/Step_1/WOTH*.csv"), read.csv)
M3_Pr <- do.call(bind_rows,M3_Pr[1:1000])


#Read in fit parameter sets

WOTH_M1 <- readRDS("data/param_sets/WOTH_M1/3_28_22/WOTH_M1_S2.rds")
WOTH_M2 <- readRDS("data/param_sets/WOTH_M2/3_30_22/WOTH_M2_S2.rds")
WOTH_M3 <- readRDS("data/param_sets/WOTH_M3/3_31_22/WOTH_M3_S2.rds")

pdf(file = "~/Documents/Coding/R/M-IBM/data/WOTH_analysis_4_23_22.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches

rec <- c()
for (each_param in 2:26){
  #Eliminate parameters that are constant (not drawn from priors)
  if(!each_param %in% c(10,11,20,21,22)){
    param <- colnames(WOTH_M1)[each_param]
    
    x_min <- min(M3_Pr[,each_param])
    x_max <- max(M3_Pr[,each_param])
    dens_temp <- density(WOTH_M1[,each_param])
    y_min <- min(dens_temp$y)
    y_max <- max(dens_temp$y)*1.25
    bw_in <- dens_temp$bw
    plot(density(WOTH_M1[,each_param]),col=alpha("dodgerblue4",.8),lwd=3,
         main=param,ylim=c(y_min,y_max),xlim=c(x_min,x_max))
    points(density(WOTH_M2[,each_param],bw=bw_in),col=alpha("firebrick3",.8),lwd=3,type="l")
    points(density(WOTH_M3[,each_param],bw=bw_in),col=alpha("black",.8),lwd=3,type="l")
    
    points(density(M1_Pr[,each_param],bw=bw_in),col=alpha("dodgerblue4",.8),lwd=2,type="l",lty=2)
    points(density(M2_Pr[,each_param],bw=bw_in),col=alpha("firebrick3",.8),lwd=2,type="l",lty=2)
    points(density(M3_Pr[,each_param],bw=bw_in),col=alpha("black",.8),lwd=2,type="l",lty=2)
    
    #Get mean estimate
    med1 <- mean(WOTH_M1[,each_param])
    med2 <- mean(WOTH_M2[,each_param])
    med3 <- mean(WOTH_M3[,each_param])
    
    #Get sd estimate
    sd1 <- sd(WOTH_M1[,each_param])
    sd2 <- sd(WOTH_M2[,each_param])
    sd3 <- sd(WOTH_M3[,each_param])
    
    #Get quantiles
    quantile_m1 <- quantile(WOTH_M1[,each_param],probs = c(.025,.975))
    quantile_m2 <- quantile(WOTH_M2[,each_param],probs = c(.025,.975))
    quantile_m3 <- quantile(WOTH_M3[,each_param],probs = c(.025,.975))
    
    #Get 95% overlap
    WOTH_M1_CI <- WOTH_M1[which(WOTH_M1[,each_param] > quantile_m1[1] & WOTH_M1[,each_param] < quantile_m1[2]),each_param]
    WOTH_M2_CI <- WOTH_M2[which(WOTH_M2[,each_param] > quantile_m2[1] & WOTH_M2[,each_param] < quantile_m2[2]),each_param]
    overlap_m3_m2 <- sum(WOTH_M2_CI > quantile_m3[1] & WOTH_M2_CI < quantile_m3[2])/length(WOTH_M2_CI)
    overlap_m2_m1 <- sum(WOTH_M1_CI > quantile_m2[1] & WOTH_M1_CI < quantile_m2[2])/length(WOTH_M1_CI)
    
    quant_add <- round(c(med1,med2,med3,sd1,sd2,sd3,overlap_m3_m2,overlap_m2_m1),4)
    p_add <- c(param,quant_add)
    rec <- rbind(rec,p_add) 
  }
}
rec <- as.data.frame(rec)
colnames(rec) <- c("param","M1_mean_est","M2_mean_est","M3_mean_est",
                   "M1_sd","M2_sd","M3_sd",
                   "M3_to_M2_overlap","M2_to_M1_overlap")
dev.off()