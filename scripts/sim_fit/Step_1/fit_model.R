# Script to train model by running multiple simulations and assessing error
library(data.table)
library(ggmap)
library(ebirdst)
library(sf)
library(dplyr)
library(lubridate)
library(rgeos)
library(cowplot)
library(readr)
library(maps)
library(spData)
library(dggridR)
library(geosphere)
library(truncnorm)
library(sigmoid)
library(pryr)
sf::sf_use_s2(FALSE)

#Number of simulations to run
num_simulations <- 10

#Number of birds
num_pulls <- 1000

#Species name
species_target <- "Varied Thrush"

#Define spacing of error grid
spacing <- 200

#Define starting season
season <- 1

#Set the Julian date based on the middle f the breeding season, as defined by eBird
start_date <- (ebirdst_runs[which(ebirdst_runs$common_name==species_target),7:8])
start_date <- as.numeric(start_date[1] + (start_date[2]-start_date[1])/2)
start_date <- yday(as.Date(start_date,origin = "1970-01-01"))
acceptable_week_dates <- seq(1,365,by=7)
start_date <- acceptable_week_dates[which(abs(acceptable_week_dates - start_date) == min(abs(acceptable_week_dates - start_date)))]
print(start_date)

#This is the Julian date of the start date for simulations
start_date <- 162

#Get functions
source(file = "scripts/sim_functions/sim_funct.R")

#Import breeding, non-breeding maps
breeding_file <- read.csv(paste("data/species_maps/breeding/",gsub(" ","_",species_target),"_breeding_all.csv",sep=""))
nonbreeding_file <- read.csv(paste("data/species_maps/non_breeding/",gsub(" ","_",species_target),"_non_breeding_all.csv",sep=""))

#Read in weekly error maps 
weekly_err <- readRDS(paste("data/species_maps/error_maps/",gsub(" ","_",species_target),"_s",spacing,".rds",sep=""))

#Read in weekly modeled areas 
weekly_modeled <- readRDS(paste("data/species_maps/modeled_area_maps/",gsub(" ","_",species_target),"_s",spacing,".rds",sep=""))

#Set up err record
err_record <- c()
error_run_list <- c()
for (each_sim in 1:num_simulations){
  print(each_sim)
  
  #Parameter sets are drawn from a script in the data folder
  source(file = "data/param_sets/VATH/3_11_22/VATH_S0.R")
  
  #Initialize model
  initial_dfs <- initialize_MIBM(num_pulls,breeding_file,nonbreeding_file,speed_mean_s,
                                 speed_sd_s,speed_mean_f,speed_sd_f,start_date_u_s,
                                 start_date_sd_s,start_date_u_f,start_date_sd_f,
                                 max_mig_time_s,max_mig_time_f,bear_err_mean_s,
                                 bear_err_sd_s,bear_err_mean_f,bear_err_sd_f,
                                 max_energy_s,max_energy_f,recovery_rate_s,
                                 recovery_rate_f,season,start_date,goal_radius,
                                 mig_con,mig_con_type,migr_timing_lat_s,migr_timing_lat_f)
  #Run model, get error
  error_run <- model_w_error(initial_dfs,start_date,spacing,weekly_err,weekly_modeled)
  #sum(error_run[[1]])
  #mean(as.numeric(error_run[[2]]))
  #par(mfrow=c(2,1))
  #plot(as.numeric(error_run[[2]]),type="l",ylim=c(.9,1))
  #plot(error_run[[1]],type="l",ylim=c(0,15),ylab="",xlab="")
  #Add the error to a record with all of the simulation variables
  err_add <- c(speed_mean_s,speed_sd_s,speed_mean_f,speed_sd_f,start_date_u_s,
               start_date_sd_s,start_date_u_f,start_date_sd_f,
               max_mig_time_s,max_mig_time_f,bear_err_mean_s,
               bear_err_sd_s,bear_err_mean_f,bear_err_sd_f,
               max_energy_s,max_energy_f,recovery_rate_s,
               recovery_rate_f,season,start_date,goal_radius,
               mig_con,mig_con_type,migr_timing_lat_s,migr_timing_lat_f,
               sum(error_run[[1]]),mean(as.numeric(error_run[[2]])))
  
  #Save the weekly error for each simulation 
  error_run_list <- rbind(error_run_list,as.numeric(error_run[[1]]))
  
  #Save the simulation variables and error
  err_record <- rbind(err_record,err_add)

}
#Give the names of variables
err_record <- as.data.frame(err_record)
colnames(err_record) <- c("speed_mean_s","speed_sd_s","speed_mean_f","speed_sd_f","start_date_u_s",
                          "start_date_sd_s","start_date_u_f","start_date_sd_f",
                          "max_mig_time_s","max_mig_time_f","bear_err_mean_s",
                          "bear_err_sd_s","bear_err_mean_f","bear_err_sd_f",
                          "max_energy_s","max_energy_f","recovery_rate_s",
                          "recovery_rate_f","season","start_date","goal_radius",
                          "mig_con","mig_con_type","migr_timing_lat_s","migr_timing_lat_f",
                          "sum_err","in_range_err")

#Combine parameter sets with weekly error sets, write to CSV
err_record <- cbind(err_record,error_run_list)
#fn <- paste("VATH_8_19_",ceiling(runif(1,1,99999)),".csv",sep="")
#write.csv(err_record,fn)
plot(NULL,xlim=c(0,52),ylim=c(0,40))
for (n in 1:nrow(err_record)){
  points(as.numeric(err_record[n,28:79]),type="l",col="blue")
}

median(err_record$sum_err)
sd(err_record$sum_err)

median(err_record$in_range_err)
sd(err_record$in_range_err)



