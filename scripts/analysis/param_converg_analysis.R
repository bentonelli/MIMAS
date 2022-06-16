# Analysis of step 2 model run
# The goal here is to look at how the model fitting process has led to convergence
# of each parameter. To do this, both the steps need to be analyzed here

library(dplyr)
library(ggplot2)

#Set number of simulations in each step
total_sims <- 100000

#Read in first fit parameter set
step1_param_set <- lapply(Sys.glob("data/training_runs/WOTH_M1/3_28_22/Step_1/WOTH*.csv"), read.csv)
step1_compiled_df <- do.call(bind_rows,step1_param_set[1:1000])

# Read in second fit information into single dataframe, then process and analyze
step2_param_set <- lapply(Sys.glob("data/training_runs/WOTH_M1/3_28_22/Step_2/WOTH*.csv"), read.csv)
step2_compiled_df <- do.call(bind_rows,step2_param_set[1:1000])

#Get best of the second step runs
perc_lim_plus <- quantile(step2_compiled_df$sum_err,probs = c(.01)) *1.03
one_perc_plus_sims <- step2_compiled_df[which(step2_compiled_df$sum_err<=perc_lim_plus),]

pdf(file = "~/Documents/Coding/R/M-IBM/data/output/WOTH_M1/3_28_22/WOTH_S2_Convergence.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches

param_convergence_rec <- c()
for (n in 2:26){
  set_range <- density(step1_compiled_df[,n])
  xx_max <- max(set_range$x)
  xx_min <- min(set_range$x)
  set_bw <- density(one_perc_plus_sims[,n])
  set_bw <- set_bw$bw
  
  #Get 95% range estimate for each stage of the fitting process
  og_range <- abs(quantile(step1_compiled_df[,n],probs=c(.025)) - quantile(step1_compiled_df[,n],probs=c(.975)))
  first_range <- abs(quantile(step2_compiled_df[,n],probs=c(.025)) - quantile(step2_compiled_df[,n],probs=c(.975)))
  second_range <- abs(quantile(one_perc_plus_sims[,n],probs=c(.025)) - quantile(one_perc_plus_sims[,n],probs=c(.975)))
  
  #Get the total convergence from the first to last step
  tot_range_convergence <- 1- (second_range/og_range)
  
  #Get median estimate of each parameter in each step, with 95% CI for final
  median_est_prior <- median(step1_compiled_df[,n])
  median_est_s1 <- median(step2_compiled_df[,n])

  median_est_s2 <- median(one_perc_plus_sims[,n])
  lower_est <- quantile(one_perc_plus_sims[,n],probs=c(.025)) 
  upper_est <- quantile(one_perc_plus_sims[,n],probs=c(.975)) 
  
  add_v <- c(colnames(one_perc_plus_sims)[n],og_range,first_range,second_range,
             median_est_prior,median_est_s1,median_est_s2,lower_est,upper_est,
             tot_range_convergence)
  param_convergence_rec <- rbind(param_convergence_rec,add_v)
  
  plot(density(one_perc_plus_sims[,n]),col="dodgerblue4",lwd=3,xlim=c(xx_min,xx_max),main=colnames(one_perc_plus_sims)[n])
  points(density(step2_compiled_df[,n],bw=set_bw),col="firebrick3",lwd=3,type="l")
  points(density(step1_compiled_df[,n],bw=set_bw),col="black",lwd=3,type="l")
}
dev.off()
colnames(param_convergence_rec) <- c("param_name","prior range","S1 range","S2 range",
                                          "prior med. est","S1 med. est.","S2 med. est","2.5%",
                                          "97.5%","% Converg.")
param_convergence_rec <- as.data.frame(param_convergence_rec)
View(param_convergence_rec)



