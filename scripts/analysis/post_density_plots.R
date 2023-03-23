# Look at estimates of traits across species
library(dplyr)
library(HDInterval)
library(RColorBrewer)
library(ggplot2)

source(file = "scripts/sim_functions/sim_funct.R")

run_date <- "8_21_22"

flc_all <- c("BRSP","BUOR","CCSP","HOWA","OROR","RNSA","TOWA","VATH","WOTH_M1","YBSA")
par(mfrow=c(2,2))
par(mar=c(4,4,4,4))
par(pty="s")

plot(NULL, xlim=c(0,1000),ylim=c(0,.005),ylab="Posterior Density",yaxt="n",xlab="Km/Flight")
for (all_spec in flc_all){
  four_letter_code <- all_spec
  #Read in all S1 runs
  final_sim_set <- readRDS(paste("data/output/",four_letter_code,"_",run_date,"/",four_letter_code,"_best.rds",sep=""))
  
  d_2 <- derived_params(final_sim_set,path = FALSE)
  points(density(d_2$fall_speed_adj),type="l",lwd=2,col=alpha("tan4",.5))
  #points(density(d_2$spring_speed_adj),type="l",lwd=4,col=alpha("forestgreen",.6))
}

plot(NULL, xlim=c(220,300),ylim=c(0,.15),ylab="Posterior Density",yaxt="n",xlab="Start Date")
for (all_spec in flc_all){
  four_letter_code <- all_spec
  #Read in all S1 runs
  final_sim_set <- readRDS(paste("data/output/",four_letter_code,"_",run_date,"/",four_letter_code,"_best.rds",sep=""))
  
  d_2 <- derived_params(final_sim_set,path = FALSE)
  points(density(final_sim_set$start_date_u_f),type="l",lwd=2,col=alpha("tan4",.5))
  #points(density(d_2$spring_speed_adj),type="l",lwd=4,col=alpha("forestgreen",.6))
}

plot(NULL, xlim=c(0,1000),ylim=c(0,.005),ylab="Posterior Density",yaxt="n",xlab="Km/Flight")
for (all_spec in flc_all){
  four_letter_code <- all_spec
  #Read in all S1 runs
  final_sim_set <- readRDS(paste("data/output/",four_letter_code,"_",run_date,"/",four_letter_code,"_best.rds",sep=""))
  
  d_2 <- derived_params(final_sim_set,path = FALSE)
  points(density(d_2$spring_speed_adj),type="l",lwd=2,col=alpha("forestgreen",.5))
  #points(density(d_2$spring_speed_adj),type="l",lwd=4,col=alpha("forestgreen",.6))
}

plot(NULL, xlim=c(60,130),ylim=c(0,.15),ylab="Posterior Density",yaxt="n",xlab="Start Date")
for (all_spec in flc_all){
  four_letter_code <- all_spec
  #Read in all S1 runs
  final_sim_set <- readRDS(paste("data/output/",four_letter_code,"_",run_date,"/",four_letter_code,"_best.rds",sep=""))
  
  d_2 <- derived_params(final_sim_set,path = FALSE)
  points(density(final_sim_set$start_date_u_s),type="l",lwd=2,col=alpha("forestgreen",.5))
  #points(density(d_2$spring_speed_adj),type="l",lwd=4,col=alpha("forestgreen",.6))
}

