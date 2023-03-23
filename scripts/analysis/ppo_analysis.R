# Get PPO between initial parameter set and final run
library(dplyr)
library(HDInterval)
library(RColorBrewer)
library(ggplot2)

source(file = "scripts/sim_functions/sim_funct.R")
#Test
run_date <- "9_20_22"

flc_all <- c("BRSP","BUOR","CCSP","HOWA","OROR","RNSA","TOWA","VATH",
             "WOTH_M1","YBSA")

all_spec_rec <- c()
for (all_spec in flc_all){
  four_letter_code <- all_spec
  #Read in all S1 runs
  training_files <- lapply(Sys.glob(paste("data/training_runs/",four_letter_code,"_",run_date,"/S1/*.csv",sep="")), read.csv)
  compiled_df <- do.call(bind_rows,training_files)
  
  final_sim_set <- readRDS(paste("data/output/",four_letter_code,"_",run_date,"/",four_letter_code,"_best.rds",sep=""))
  
  ppo_pr <- postPriorOverlap(final_sim_set$start_date_u_s,compiled_df$start_date_u_s)
  
  d_1 <- derived_params(compiled_df,path = FALSE)
  d_2 <- derived_params(final_sim_set,path = FALSE)
  
  ppo_start_date_u_s <- postPriorOverlap(final_sim_set$start_date_u_s,compiled_df$start_date_u_s)
  ppo_start_date_sd_s <- postPriorOverlap(final_sim_set$start_date_sd_s,compiled_df$start_date_sd_s)
  ppo_start_date_u_f <- postPriorOverlap(final_sim_set$start_date_u_f,compiled_df$start_date_u_f)
  ppo_start_date_sd_f <- postPriorOverlap(final_sim_set$start_date_sd_f,compiled_df$start_date_sd_f)
  
  ppo_bear_err_mean_s <- postPriorOverlap(final_sim_set$bear_err_mean_s,compiled_df$bear_err_mean_s)
  ppo_bear_err_sd_s <- postPriorOverlap(final_sim_set$bear_err_sd_s,compiled_df$bear_err_sd_s)
  
  ppo_bear_err_mean_f <- postPriorOverlap(final_sim_set$bear_err_mean_f,compiled_df$bear_err_mean_f)
  ppo_bear_err_sd_f <- postPriorOverlap(final_sim_set$bear_err_sd_f,compiled_df$bear_err_sd_f)
  
  ppo_max_energy_s <- postPriorOverlap(final_sim_set$max_energy_s,compiled_df$max_energy_s)
  ppo_max_energy_f <- postPriorOverlap(final_sim_set$max_energy_f,compiled_df$max_energy_f)
  
  ppo_migr_timing_lat_s <- postPriorOverlap(final_sim_set$migr_timing_lat_s,compiled_df$migr_timing_lat_s)
  ppo_migr_timing_lat_f <- postPriorOverlap(final_sim_set$migr_timing_lat_f,compiled_df$migr_timing_lat_f)
  
  ppo_migr_con <- postPriorOverlap(final_sim_set$mig_con,compiled_df$mig_con)
  
  ppo_spring_speed_adj <- postPriorOverlap(d_2$spring_speed_adj,d_1$spring_speed_adj)
  ppo_fall_speed_adj <- postPriorOverlap(d_2$fall_speed_adj,d_1$fall_speed_adj)
  ppo_rec_rate_s <- postPriorOverlap(d_2$rec_rate_s,d_1$rec_rate_s)
  ppo_rec_rate_f <- postPriorOverlap(d_2$rec_rate_f,d_1$rec_rate_f)

  sp_add <- c(ppo_start_date_u_s,ppo_start_date_sd_s,ppo_start_date_u_f,ppo_start_date_sd_f,
              ppo_bear_err_mean_s,ppo_bear_err_sd_s,ppo_bear_err_mean_f,ppo_bear_err_sd_f,
              ppo_max_energy_s,ppo_max_energy_f,ppo_migr_timing_lat_s,ppo_migr_timing_lat_f,
              ppo_migr_con,ppo_spring_speed_adj,ppo_fall_speed_adj,ppo_rec_rate_s,ppo_rec_rate_f)
  all_spec_rec <- rbind(all_spec_rec,sp_add)
}
dev.off()
all_spec_rec <- as.data.frame(all_spec_rec)
colnames(all_spec_rec) <- c("ppo_start_date_u_s","ppo_start_date_sd_s","ppo_start_date_u_f","ppo_start_date_sd_f",
                            "ppo_bear_err_mean_s","ppo_bear_err_sd_s","ppo_bear_err_mean_f","ppo_bear_err_sd_f",
                            "ppo_max_energy_s","ppo_max_energy_f","ppo_migr_timing_lat_s","ppo_migr_timing_lat_f",
                            "ppo_migr_con","ppo_spring_speed_adj","ppo_fall_speed_adj","ppo_rec_rate_s","ppo_rec_rate_f")
all_spec_rec$spec <- flc_all
saveRDS(all_spec_rec,"ppo_analysis_2_24.rds")

all_spec_rec$spec_num <- 1:nrow(all_spec_rec)

ax_lab <- c("Start date, mean - Spring","Start date, s.d. - Spring",
             "Start date, mean Fall","Start date, s.d. - Fall",
  "Bear. error, mean - Spring","Bear. error, s.d. - Spring","Bear. error, mean Fall","Bear. error, s.d. - Fall",
  "Max. energy - Spring","Max. energy - Fall","Migr. timing by lat. - Spring","Migr. timing by lat. - Fall",
  "Migratory connectivity","Flight distance, mean - Spring","Flight distance, mean - Fall",
  "Recovery rate, Spring","Recovery rate, Fall")

cr <- c(brewer.pal(12,"Set3"),"black","forestgreen","firebrick3","dodgerblue4","orchid")

par(mar=(c(6,20,2,2)))
plot(x=NULL,y=NULL,xlim=c(0,1),ylim=c(1,17),yaxt="n",ylab="",xlab="PPO",cex.lab=1.65,cex.axis=1.65)
par(las=1)
axis(2,at=1:17,labels=rev(ax_lab),cex.lab=1.25,cex.axis=1.25)
for (nn in rev(1:17)){
  print(nn)
  points(all_spec_rec[,nn],18-rep(nn,length(all_spec_rec[,nn])),pch=19,col=alpha("dodgerblue3",.6),cex=2.25)
  points(mean(all_spec_rec[,nn]),18-nn,pch="|",col=alpha("firebrick3",.8),cex=2)
}

#Add letters
points(all_spec_rec[8,4],14,pch=19,col=alpha("sienna2",1),cex=3)
text(all_spec_rec[8,4],14,"B",cex=1.25) #Start date, SD, Fall VATH

points(all_spec_rec[4,13],5,pch=19,col=alpha("sienna2",1),cex=3)
text(all_spec_rec[4,13],5,"C",cex=1.25) # Migratory connectivity, Hooded Warbler

points(all_spec_rec[3,15],3,pch=19,col=alpha("sienna2",1),cex=3)
text(all_spec_rec[3,15],3,"D",cex=1.25) #Flight distance, CCSP 

# Plot 1 by 3 square plots of examples
par(mfrow=c(3,1))
par(mar=(c(5,2,2,2)))
par(pty="s")
four_letter_code <- "VATH"
run_date <- "9_20_22"
#Read in all S1 runs
training_files <- lapply(Sys.glob(paste("data/training_runs/",four_letter_code,"_",run_date,"/S1/*.csv",sep="")), read.csv)
compiled_df <- do.call(bind_rows,training_files)

final_sim_set <- readRDS(paste("data/output/",four_letter_code,"_",run_date,"/",four_letter_code,"_best.rds",sep=""))

ppo_pr <- postPriorOverlap(final_sim_set$start_date_sd_f,compiled_df$start_date_sd_f,xlab="Start date, s.d., Spring")

four_letter_code <- "HOWA"
run_date <- "9_20_22"
#Read in all S1 runs
training_files <- lapply(Sys.glob(paste("data/training_runs/",four_letter_code,"_",run_date,"/S1/*.csv",sep="")), read.csv)
compiled_df <- do.call(bind_rows,training_files)

final_sim_set <- readRDS(paste("data/output/",four_letter_code,"_",run_date,"/",four_letter_code,"_best.rds",sep=""))

ppo_pr <- postPriorOverlap(final_sim_set$mig_con,compiled_df$mig_con,xlab="Migratory connectivity")

four_letter_code <- "CCSP"
run_date <- "9_20_22"
#Read in all S1 runs
training_files <- lapply(Sys.glob(paste("data/training_runs/",four_letter_code,"_",run_date,"/S1/*.csv",sep="")), read.csv)
compiled_df <- do.call(bind_rows,training_files)

final_sim_set <- readRDS(paste("data/output/",four_letter_code,"_",run_date,"/",four_letter_code,"_best.rds",sep=""))
d_1 <- derived_params(compiled_df,path = FALSE)
d_2 <- derived_params(final_sim_set,path = FALSE)

ppo_pr <- postPriorOverlap(d_2$fall_speed_adj,d_1$fall_speed_adj,xlab="Flight distance, mean - Fall")



