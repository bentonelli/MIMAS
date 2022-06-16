### Script to look at paramter estimates for multiple species ####
library(dplyr)
source(file = "scripts/sim_functions/sim_funct.R")
#Read in fit parameter sets

param_convergence_files <- c("data/output/BRSP/4_8_22/param_convergence.csv",
                             "data/output/BUOR/4_6_22/param_convergence.csv",
                             "data/output/CCSP/3_7_22/param_convergence.csv",
                             "data/output/HOWA/3_12_22/param_convergence.csv",
                             "data/output/OROR/4_7_22/param_convergence.csv",
                             "data/output/RNSA/4_5_22/param_convergence.csv",
                             "data/output/TOWA/4_4_22/param_convergence.csv",
                             "data/output/VATH/3_11_22/param_convergence.csv",
                             "data/output/WOTH_M1/3_28_22/param_convergence.csv",
                             "data/output/YBSA/4_10_22/param_convergence.csv")

S2_files<- c("data/param_sets/BRSP/4_8_22/BRSP_S2.rds",
                             "data/param_sets/BUOR/4_6_22/BUOR_S2.rds",
                             "data/param_sets/CCSP/3_7_22/CCSP_S2.rds",
                             "data/param_sets/HOWA/3_12_22/HOWA_S2.rds",
                             "data/param_sets/OROR/4_7_22/OROR_S2.rds",
                             "data/param_sets/RNSA/4_5_22/RNSA_S2.rds",
                             "data/param_sets/TOWA/4_4_22/TOWA_S2.rds",
                             "data/param_sets/VATH/3_11_22/VATH_S2.rds",
                             "data/param_sets/WOTH_M1/3_28_22/WOTH_M1_S2.rds",
                             "data/param_sets/YBSA/4_10_22/YBSA_S2.rds")

spec_list <- c("Brewer's Sparrow","Bullock's Oriole","Clay-colored Sparrow",
               "Hooded Warbler","Orchard Oriole","Red-naped Sapsucker",
               "Townsend's Warbler","Varied Thrush","Wood Thrush","Yellow-bellied Sapsucker")

taxon <- c("Sparrow","Oriole","Sparrow","Warbler","Oriole","Woodpecker","Warbler","Thrush","Thrush","Woodpecker")
e_c_w <- c("West","West","Central",
           "East","East","West",
           "West","West","East","East")

spec_weight <- c(11,36,12,11,20,50,9,79,45,48)
all_spec_rec <- c()
flight_speeds <- c()
for (n in 1:length(param_convergence_files)){
  
  spec_name <- spec_list[n]
  spec_in <- readRDS(S2_files[n])
  
  dr_params <- derived_params(S2_files[n])
  
  mean_fs_s <- mean(dr_params[[1]])
  fs_s_q <- quantile(dr_params[[1]],probs = c(.045,.945))
  
  mean_fs_f <- mean(dr_params[[2]])
  fs_f_q <- quantile(dr_params[[2]],probs = c(.045,.945))
  
  rec_rate_s <- mean(dr_params[[3]])
  rec_rate_f <- mean(dr_params[[4]])
  
  flight_speeds <- rbind(flight_speeds,c(spec_list[n],mean_fs_s,fs_s_q,mean_fs_f,fs_f_q))
  
  #Get breeding data
  #breeding_file <- read.csv(paste("data/species_maps/breeding/",gsub(" ","_",spec_name),"_breeding_all.csv",sep=""))
  
  #Get random locations - as vector
  #rand_loc <- sample(nrow(breeding_file),num_pulls,prob=breeding_file$perc_occupancy,replace=TRUE)
  #Save lat-lng
  #rand_lat_long <- breeding_file[rand_loc,1:2]
  #rand_lat_long$long_trans <- rand_lat_long$longitude
  #rand_lat_long$lat_trans <- rand_lat_long$latitude
  
  #rand_lat_long$long_trans <- ebird_lat_long_tran(rand_lat_long$longitude,rand_lat_long$lat_trans)[,1]
  #rand_lat_long$lat_trans <- ebird_lat_long_tran(rand_lat_long$longitude,rand_lat_long$lat_trans)[,2]
  
  #save_points <- as.matrix(rand_lat_long[,3:4])
  
  #avg_lat <- mean(save_points[,2])
  #sd_lat <- sd(save_points[,2])
  
  s2_mean <- mean(spec_in$start_date_u_s)
  s2_range <- quantile(spec_in$start_date_u_s,probs = c(.045,.945))
  
  s3_mean <- mean(spec_in$bear_err_mean_s)
  s3_range <- quantile(spec_in$bear_err_mean_s,probs = c(.045,.945))
  
  s4_mean <- mean(spec_in$bear_err_mean_f)
  s4_range <- quantile(spec_in$bear_err_mean_f,probs = c(.045,.945))
  
  s5_mean <- mean(spec_in$bear_err_sd_s)
  s5_range <- quantile(spec_in$bear_err_sd_s,probs = c(.045,.945))
  
  s6_mean <- mean(spec_in$bear_err_sd_f)
  s6_range <- quantile(spec_in$bear_err_sd_f,probs = c(.045,.945))

  add <- c(spec_name,s2_mean,s2_range,s3_mean,s3_range,s4_mean,s4_range,s5_mean,s5_range,s6_mean,s6_range)
  all_spec_rec <- rbind(all_spec_rec,add)
}
flight_speeds <- as.data.frame(flight_speeds)
colnames(flight_speeds) <- c("Species","Mean Spring Flight Dist.","Spring 2.5%","Spring 94.5%",
                             "Mean Fall Flight Dist.","Fall 2.5%","Fall 94.5%")
all_spec_rec <- as.data.frame(all_spec_rec)

colnames(all_spec_rec) <- c("Species","Mean_SD_spring","SD_spring 4.5%","SD_spring 94.5%",
                            "Mean_bearing_spring","bearing_s 4.5%","bearing_s 94.5%",
                            "Mean_bearing_fall","bearing_f 4.5%","bearing_f 94.5%",
                            "Mean_bearing_sd_spring","bearing_s_sd 4.5%","bearing_s_sd 94.5%",
                            "Mean_bearing_sd_fall","bearing_f_sd 4.5%","bearing_f_sd 94.5%")

all_spec_rec_f <- merge(all_spec_rec,flight_speeds,by.x="Species")
all_spec_rec_f$taxon <- taxon
all_spec_rec_f$loc <- e_c_w
all_spec_rec_f <- arrange(all_spec_rec_f,taxon)


par(mar = c(6, 15, 2, 2))
plot(NULL,xlim=c(70,125),ylim=c(0.75,10.25),ylab="",xlab="Julian Day",yaxt="n")
axis(2, at = 1:10,
     labels = rev(all_spec_rec_f$Species),
     las=1,
     main = "Horizontal")
points(all_spec_rec_f$Mean_SD_spring,10:1,pch=19,cex=1.5)
arrows(as.numeric(all_spec_rec_f$`SD_spring 4.5%`),10:1,as.numeric(all_spec_rec_f$`SD_spring 94.5%`),10:1,length = 0,lwd=2)

par(mfrow=c(1,3))

par(mar = c(6, 12, 1, 1))
plot(NULL,xlim=c(0,750),ylim=c(0.75,10.25),ylab="",xlab="Flight Dist.",yaxt="n")
axis(2, at = 1:10,
     labels = rev(all_spec_rec_f$Species),
     las=1,
     main = "Horizontal")
abline(h=c(seq(1.5,9.5,by=1)),lty=2)
points(all_spec_rec_f$`Mean Spring Flight Dist.`,10:1-.25,pch=19,cex=1.5,col="forestgreen")
points(all_spec_rec_f$`Mean Fall Flight Dist.`,10:1+.25,pch=19,cex=1.5,col="orange3")
arrows(as.numeric(all_spec_rec_f$`Spring 2.5%`),10:1-.25,as.numeric(all_spec_rec_f$`Spring 94.5%`),10:1-.25,length = 0,lwd=2,col="forestgreen")
arrows(as.numeric(all_spec_rec_f$`Fall 2.5%`),10:1+.25,as.numeric(all_spec_rec_f$`Fall 94.5%`),10:1+.25,length = 0,lwd=2,col="orange3")

par(mar = c(6, 12, 1, 1))
plot(NULL,xlim=c(0,80),ylim=c(0.75,10.25),ylab="",xlab="Bearing Dispersion",yaxt="n")
#axis(2, at = 1:10,
#     labels = rev(all_spec_rec_f$Species),
#     las=1,
#     main = "Horizontal")
abline(h=c(seq(1.5,9.5,by=1)),lty=2)
points(all_spec_rec_f$Mean_bearing_sd_spring,10:1-.25,pch=19,cex=1.5,col="forestgreen")
points(all_spec_rec_f$Mean_bearing_sd_fall,10:1+.25,pch=19,cex=1.5,col="orange3")
arrows(as.numeric(all_spec_rec_f$`bearing_s_sd 4.5%`),10:1-.25,as.numeric(all_spec_rec_f$`bearing_s_sd 94.5%`),10:1-.25,length = 0,lwd=2,col="forestgreen")
arrows(as.numeric(all_spec_rec_f$`bearing_f_sd 4.5%`),10:1+.25,as.numeric(all_spec_rec_f$`bearing_f_sd 94.5%`),10:1+.25,length = 0,lwd=2,col="orange3")

par(mar = c(6, 12, 1, 1))
plot(NULL,xlim=c(-20,20),ylim=c(0.75,10.25),ylab="",xlab="Deviation, Shortest Route",yaxt="n")
#axis(2, at = 1:10,
#     labels = rev(all_spec_rec_f$Species),
#     las=1,
#     main = "Horizontal")
abline(h=c(seq(1.5,9.5,by=1)),lty=2)
abline(v=0,lty=3)
points(all_spec_rec_f$Mean_bearing_spring,10:1-.25,pch=19,cex=1.5,col="forestgreen")
points(all_spec_rec_f$Mean_bearing_fall,10:1+.25,pch=19,cex=1.5,col="orange3")
arrows(as.numeric(all_spec_rec_f$`bearing_s 4.5%`),10:1-.25,as.numeric(all_spec_rec_f$`bearing_s 94.5%`),10:1-.25,length = 0,lwd=2,col="forestgreen")
arrows(as.numeric(all_spec_rec_f$`bearing_f 4.5%`),10:1+.25,as.numeric(all_spec_rec_f$`bearing_f 94.5%`),10:1+.25,length = 0,lwd=2,col="orange3")
