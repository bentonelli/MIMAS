# Load necessary packages (make sure to install, if you don't have these already)
library(ebirdst)
library(lubridate)
library(geosphere)
library(dplyr)
library(truncnorm)
library(sf)
library(spData)
#Read in functions (you may have just done this in the last code block)
source(file = "scripts/sim_functions/sim_funct.R")

#Define species
species_target <- "Townsend's Warbler"

# Read in parameter set of trained simulations. Make sure to check out which species
# have trained models. Also note these use species 4-letter codes.
param_set <- readRDS("data/output/Spec_IBM_output/TOWA_9_20_22/TOWA_best.rds")

# Check out the structure of this data frame: each row is an accepted simulation
# Each row is a parameter value.
head(param_set[,1:10])

# We are going to pull one of these simulations/parameter sets at random to use,
# and we are going to save these to the environment
rand_sim <- sample(1:nrow(param_set),1)
for (each_var in 1:ncol(param_set)){
  assign(colnames(param_set)[each_var],param_set[rand_sim,each_var])
}

#Set number of birds to include in simulation, 1000 birds is a reasonable starting point
num_pulls <- 1000

#Set number of days to run the simulation for.
num_days <- 365

#Import breeding, non-breeding maps
breeding_file <- read.csv(paste("data/species_maps/breeding/",gsub(" ","_",species_target),"_breeding_all.csv",sep=""))
nonbreeding_file <- read.csv(paste("data/species_maps/non_breeding/",gsub(" ","_",species_target),"_non_breeding_all.csv",sep=""))

#Define starting season, leave this at 1 (breeding season)
season <- 1

#We will set the start date as the midpoint of the breeding season.
start_date <- (ebirdst_runs[which(ebirdst_runs$common_name==species_target),7:8])
start_date <- as.numeric(start_date[1] + (start_date[2]-start_date[1])/2)
start_date <- yday(as.Date(start_date,origin = "1970-01-01"))
acceptable_week_dates <- seq(1,365,by=7)
start_date <- acceptable_week_dates[which(abs(acceptable_week_dates - start_date) == min(abs(acceptable_week_dates - start_date)))]

#Initialize model using parameters pulled above
initial_dfs <- initialize_MIBM(num_pulls,breeding_file,nonbreeding_file,speed_mean_s,
                               speed_sd_s,speed_mean_f,speed_sd_f,start_date_u_s,
                               start_date_sd_s,start_date_u_f,start_date_sd_f,
                               max_mig_time_s,max_mig_time_f,bear_err_mean_s,
                               bear_err_sd_s,bear_err_mean_f,bear_err_sd_f,
                               max_energy_s,max_energy_f,recovery_rate_s,
                               recovery_rate_f,season,start_date,goal_radius,
                               mig_con,migr_timing_lat_s,migr_timing_lat_f)

# The IBMs rely on two dataframes - one that stays the same (Static) and the other
# that updates based on the location, status of individual birds
static_df <- initial_dfs[[1]]
upd_df <- initial_dfs[[2]]

#Set up a dataframe to record the locations of birds on each day
locs_rec <- as.data.frame(matrix(NA,ncol = 5,nrow = (num_days*num_pulls)))
colnames(locs_rec) <- c("ID","Lon","Lat","Day","Season")
for (days in 1:num_days){
  upd_df <- run_day(upd_df,static_df)
  locs <- save_locations(upd_df,save_season = TRUE)
  loc_ind_s <- (days - 1)*num_pulls + 1 
  loc_ind_e <- (days - 1)*num_pulls + num_pulls 
  locs_rec[loc_ind_s:loc_ind_e,] <- locs
}
head(locs_rec)