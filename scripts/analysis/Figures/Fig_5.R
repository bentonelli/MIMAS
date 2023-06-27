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
  
  final_sim_set <- readRDS(paste("data/output/Spec_IBM_output/",four_letter_code,"_",run_date,"/",four_letter_code,"_best.rds",sep=""))
  
  ppo_pr <- postPriorOverlap(final_sim_set$start_date_u_s,compiled_df$start_date_u_s)
  
  d_1 <- derived_params(compiled_df,path = FALSE)
  d_2 <- derived_params(final_sim_set,path = FALSE)
  
  # 1-6 = timing
  ppo_start_date_u_s <- postPriorOverlap(final_sim_set$start_date_u_s,compiled_df$start_date_u_s) #1
  ppo_start_date_u_f <- postPriorOverlap(final_sim_set$start_date_u_f,compiled_df$start_date_u_f) #1
  
  ppo_start_date_sd_s <- postPriorOverlap(final_sim_set$start_date_sd_s,compiled_df$start_date_sd_s) #2
  ppo_start_date_sd_f <- postPriorOverlap(final_sim_set$start_date_sd_f,compiled_df$start_date_sd_f) #2
  
  ppo_migr_timing_lat_s <- postPriorOverlap(final_sim_set$migr_timing_lat_s,compiled_df$migr_timing_lat_s) #3
  ppo_migr_timing_lat_f <- postPriorOverlap(final_sim_set$migr_timing_lat_f,compiled_df$migr_timing_lat_f) #3
  
  # 7-8 = speed
  ppo_spring_speed_adj <- postPriorOverlap(d_2$spring_speed_adj,d_1$spring_speed_adj) #6
  ppo_fall_speed_adj <- postPriorOverlap(d_2$fall_speed_adj,d_1$fall_speed_adj) #6
  
  # 9-12 = orientation
  ppo_bear_err_mean_s <- postPriorOverlap(final_sim_set$bear_err_mean_s,compiled_df$bear_err_mean_s) #7
  po_bear_err_mean_f <- postPriorOverlap(final_sim_set$bear_err_mean_f,compiled_df$bear_err_mean_f) #7
  
  ppo_bear_err_sd_s <- postPriorOverlap(final_sim_set$bear_err_sd_s,compiled_df$bear_err_sd_s) #8
  ppo_bear_err_sd_f <- postPriorOverlap(final_sim_set$bear_err_sd_f,compiled_df$bear_err_sd_f) #8
  
  #13 MC
  ppo_migr_con <- postPriorOverlap(final_sim_set$mig_con,compiled_df$mig_con) #9
  
  # 14-17 Energetics
  ppo_max_energy_s <- postPriorOverlap(final_sim_set$max_energy_s,compiled_df$max_energy_s) #10
  ppo_max_energy_f <- postPriorOverlap(final_sim_set$max_energy_f,compiled_df$max_energy_f) #10
  
  ppo_rec_rate_s <- postPriorOverlap(d_2$rec_rate_s,d_1$rec_rate_s) #12
  ppo_rec_rate_f <- postPriorOverlap(d_2$rec_rate_f,d_1$rec_rate_f) #12
  
  sp_add <- c(ppo_start_date_u_s,ppo_start_date_u_f,
              ppo_start_date_sd_s,ppo_start_date_sd_f,
              ppo_migr_timing_lat_s,ppo_migr_timing_lat_f,
              ppo_spring_speed_adj,ppo_fall_speed_adj,
              ppo_bear_err_mean_s,po_bear_err_mean_f,
              ppo_bear_err_sd_s,ppo_bear_err_sd_f,
              ppo_migr_con,
              ppo_max_energy_s,ppo_max_energy_f,
              ppo_rec_rate_s,ppo_rec_rate_f)
  all_spec_rec <- rbind(all_spec_rec,sp_add)
}
dev.off()
all_spec_rec <- as.data.frame(all_spec_rec)
colnames(all_spec_rec) <- c("ppo_start_date_u_s","ppo_start_date_u_f",
                            "ppo_start_date_sd_s","ppo_start_date_sd_f",
                            "ppo_migr_timing_lat_s","ppo_migr_timing_lat_f",
                            "ppo_spring_speed_adj","ppo_fall_speed_adj",
                            "ppo_bear_err_mean_s","po_bear_err_mean_f",
                            "ppo_bear_err_sd_s","ppo_bear_err_sd_f",
                            "ppo_migr_con",
                            "ppo_max_energy_s","ppo_max_energy_f",
                            "ppo_rec_rate_s",'ppo_rec_rate_f')
all_spec_rec$spec <- flc_all
saveRDS(all_spec_rec,"ppo_analysis_2_24.rds")

all_spec_rec$spec_num <- 1:nrow(all_spec_rec)

ax_lab <- as.character(c(1,1,2,2,3,3,6,6,7,7,8,8,9,10,10,12,12))

spring_col <- "forestgreen"
fall_col <- "tan3"

spring_shp <- 19
fall_shp <- 18
both_shp <- 17

fall_spring_cols <- c(spring_col,fall_col,
                      spring_col,fall_col,
                      spring_col,fall_col,
                      spring_col,fall_col,
                      spring_col,fall_col,
                      spring_col,fall_col,
                      "dodgerblue4",
                      spring_col,fall_col,
                      spring_col,fall_col)

fall_spring_shape <- c(spring_shp,fall_shp,
                       spring_shp,fall_shp,
                       spring_shp,fall_shp,
                       spring_shp,fall_shp,
                       spring_shp,fall_shp,
                       spring_shp,fall_shp,
                       both_shp,
                       spring_shp,fall_shp,
                       spring_shp,fall_shp)

cr <- c(brewer.pal(12,"Set3"),"black","forestgreen","firebrick3","dodgerblue4","orchid")

par(mar=(c(6,4,2,2)))
plot(x=NULL,y=NULL,xlim=c(0,1),ylim=c(1,17),yaxt="n",ylab="",xlab="Prior-posterior overlap",cex.lab=1.65,cex.axis=1.65)
#Guidelines
abline(h=1:17,lty=2,lwd=2,col="grey80")
par(las=1)
axis(2,at=1:17,labels=rev(ax_lab),cex.lab=1.65,cex.axis=1.65)
for (nn in rev(1:17)){
  print(nn)
  points(all_spec_rec[,nn],18-rep(nn,length(all_spec_rec[,nn])),pch=fall_spring_shape[nn],col=alpha(fall_spring_cols[nn],.7),cex=2.5)
  points(mean(all_spec_rec[,nn]),18-nn,pch="|",col=alpha("firebrick3",.8),cex=2)
}

#Add letters
points(all_spec_rec[8,2],16,pch=fall_shp,col=alpha("grey75",1),cex=3)
text(all_spec_rec[8,2]+.0015,16,"B",cex=1.15) #Start date, mean, fall VATH

points(all_spec_rec[1,8],10,pch=fall_shp,col=alpha("grey75",1),cex=3)
text(all_spec_rec[1,8]+.0015,10,"C",cex=1.15) #Flight distance, BRSP 

points(all_spec_rec[4,13],5,pch=both_shp,col=alpha("grey75",1),cex=3)
text(all_spec_rec[4,13]+.0015,5,"D",cex=1.15) # Migratory connectivity, Hooded Warbler

# Plot 1 by 3 square plots of examples
par(mfrow=c(3,1))
par(mar=(c(5,2,0,2)))
par(pty="s")
four_letter_code <- "VATH"
run_date <- "9_20_22"
#Read in all S1 runs
training_files <- lapply(Sys.glob(paste("data/training_runs/",four_letter_code,"_",run_date,"/S1/*.csv",sep="")), read.csv)
compiled_df <- do.call(bind_rows,training_files)

final_sim_set <- readRDS(paste("data/output/",four_letter_code,"_",run_date,"/",four_letter_code,"_best.rds",sep=""))


ppo_pr <- postPriorOverlap(final_sim_set$start_date_u_f,compiled_df$start_date_u_f,
                           xlab=expression(paste("Departure date, ",mu," - Fall")),
                           mainColor="dodgerblue3", priorColor="firebrick3", overlapColor="purple4",
                           cex.lab=2, cex=1.75, cex.axis=1.5)

four_letter_code <- "BRSP"
run_date <- "9_20_22"
#Read in all S1 runs
training_files <- lapply(Sys.glob(paste("data/training_runs/",four_letter_code,"_",run_date,"/S1/*.csv",sep="")), read.csv)
compiled_df <- do.call(bind_rows,training_files)

final_sim_set <- readRDS(paste("data/output/",four_letter_code,"_",run_date,"/",four_letter_code,"_best.rds",sep=""))
d_1 <- derived_params(compiled_df,path = FALSE)
d_2 <- derived_params(final_sim_set,path = FALSE)

ppo_pr <- postPriorOverlap(d_2$fall_speed_adj,d_1$fall_speed_adj,
                           xlab=expression(paste("Flight distance, ",mu," - Fall")),
                           mainColor="dodgerblue3", priorColor="firebrick3", overlapColor="purple4",
                           cex.lab=2, cex=1.75, cex.axis=1.5)

four_letter_code <- "HOWA"
run_date <- "9_20_22"
#Read in all S1 runs
training_files <- lapply(Sys.glob(paste("data/training_runs/",four_letter_code,"_",run_date,"/S1/*.csv",sep="")), read.csv)
compiled_df <- do.call(bind_rows,training_files)

final_sim_set <- readRDS(paste("data/output/",four_letter_code,"_",run_date,"/",four_letter_code,"_best.rds",sep=""))

ppo_pr <- postPriorOverlap(final_sim_set$mig_con,compiled_df$mig_con,
                           xlab="Migratory connectivity",
                           mainColor="dodgerblue3", priorColor="firebrick3", overlapColor="purple4",
                           cex.lab=2, cex=1.75, cex.axis=1.5)


