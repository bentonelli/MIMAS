#Fit parameters from recombined posteriors

Sys.setenv(TZ='America/Los_Angeles')
setwd("~/R/M-IBM/")
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

#Load in posterior of step 1
posterior_params <- readRDS("data/param_sets/YBSA/4_10_22/YBSA_S1.rds")

# Determine names of output files
file_name_start <- "YBSA_S2_4_10_22_"

#Number of simulations to run
num_simulations <- 100

#Number of birds
num_pulls <- 5000

#Species name
species_target <- "Yellow-bellied Sapsucker"

#Define spacing of error grid
spacing <- 200

### Draw parameter sets ####
# The first step is to re-draw parameters from the posteriors, taking into account
# that the parameters come from two different, and more or less independent, 
# fall and spring models. To do this, independent parameter estimates can be drawn,
# and used to re-run models that now represent a greatly reduced parameter space.

top_f_models <- posterior_params[[1]]
top_s_models <- posterior_params[[2]]
top_f_s_models <- posterior_params[[3]]

fall_param_inx <- c(4,5,8,9,11,14,15,17,19,26)
spring_param_inx <- c(2,3,6,7,10,12,13,16,18,25)
shared_inx <- c(20,21,22,23,24)

n_param_sets <- num_simulations
param_sets <- as.data.frame(matrix(NA, ncol = length(c(fall_param_inx,spring_param_inx,shared_inx)),
                                   nrow = n_param_sets))
for (n in 1:n_param_sets){
  
  #Grab a random set of fall parameters
  fll_rand <- (sample(1:nrow(top_f_models),length(fall_param_inx)))
  fall_set <- diag(as.matrix(top_f_models[fll_rand,fall_param_inx]))
  #Grab a random set of spring parameters
  spr_rand <- (sample(1:nrow(top_s_models),length(spring_param_inx)))
  spring_set <- diag(as.matrix(top_s_models[spr_rand,spring_param_inx]))
  
  #grab a random shared set - coin flip between the two param sets chosen above
  all_rand <- (sample(1:nrow(top_f_s_models),length(shared_inx)))
  shared_set <- diag(as.matrix(top_f_s_models[all_rand,shared_inx]))
  
  new_param_set <- c(fall_set,spring_set,shared_set)
  param_sets[n,] <- new_param_set
}
colnames(param_sets) <- c(colnames(top_f_models)[fall_param_inx],
                          colnames(top_f_models)[spring_param_inx],
                          colnames(top_f_s_models)[shared_inx])

#This is the Julian date of the start date for simulations
start_date <- param_sets$start_date[1]

#Define starting season
season <- param_sets$season[1]

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
  
  #Parameter sets are drawn from the read .rds file
  for (each_var in 1:ncol(param_sets)){
    assign(colnames(param_sets)[each_var],param_sets[each_sim,each_var])
  }
  
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
fn <- paste(file_name_start,ceiling(runif(1,1,9999999)),".csv",sep="")
write.csv(err_record,fn)



