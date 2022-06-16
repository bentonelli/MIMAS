#This script runs a single M-IBM simulation

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
sf::sf_use_s2(FALSE)

setwd("~/Documents/Coding/R/M-IBM")

#Get functions
source(file = "scripts/sim_functions/sim_funct.R")

#Parameter sets should be drawn from a script in the data folder
source(file = "data/param_sets/Varied_Thrush_MLM_9_20_21.R")

#Define species
species_target <- "Varied Thrush"

#Set number of birds to include in simulation
num_pulls <- 10000

#Define spacing of error maps
spacing <- 200

#Import breeding, non-breeding maps
breeding_file <- read.csv(paste("data/species_maps/breeding/",gsub(" ","_",species_target),"_breeding_all.csv",sep=""))
nonbreeding_file <- read.csv(paste("data/species_maps/non_breeding/",gsub(" ","_",species_target),"_non_breeding_all.csv",sep=""))

#Read in weekly error maps - not explicitly needed here in the simple run, but good to know how to do
weekly_err <- readRDS(paste("data/species_maps/error_maps/",gsub(" ","_",species_target),"_s",spacing,".rds",sep=""))

#Read in weekly modeled areas - also not explicitly needed
weekly_modeled <- readRDS(paste("data/species_maps/modeled_area_maps//",gsub(" ","_",species_target),"_s",spacing,".rds",sep=""))

#Define starting season
season <- 1

#This is the Julian date for May 28th.
start_date <- 148

#Initialize model
initial_dfs <- initialize_MIBM(num_pulls,breeding_file,nonbreeding_file,speed_mean_s,
                               speed_sd_s,speed_mean_f,speed_sd_f,start_date_u_s,
                               start_date_sd_s,start_date_u_f,start_date_sd_f,
                               max_mig_time_s,max_mig_time_f,bear_err_mean_s,
                               bear_err_sd_s,bear_err_mean_f,bear_err_sd_f,
                               max_energy_s,max_energy_f,recovery_rate_s,
                               recovery_rate_f,season,start_date,goal_radius,
                               mig_con,mig_con_type,migr_timing_lat_s,migr_timing_lat_f)
static_df <- initial_dfs[[1]]
upd_df <- initial_dfs[[2]]
locs_rec <- as.data.frame(matrix(NA,ncol = 4,nrow = (365*num_pulls)))
colnames(locs_rec) <- c("ID","Lon","Lat","Day")
for (days in 1:365){
  print(days)
  upd_df <- run_day(upd_df,static_df)
  locs <- save_locations(upd_df)
  loc_ind_s <- (days - 1)*num_pulls + 1 
  loc_ind_e <- (days - 1)*num_pulls + num_pulls 
  locs_rec[loc_ind_s:loc_ind_e,] <- locs
}


