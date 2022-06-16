#This script runs a single M-IBM simulation and creates visualizations of the results.

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
library(gganimate)
library(PBSmapping)
library(data.table)
library(gifski)
sf::sf_use_s2(FALSE)

setwd("~/Documents/Coding/R/M-IBM")

#Get functions
source(file = "scripts/sim_functions/sim_funct.R")

# Parameter sets should be drawn from a script in the data folder, or from a 
# fit parameter set
#Script
#source(file = "data/param_sets/Clay_colored_Sparrow_BG_3_1_22.R")
#Fit Data
param_set <- readRDS("data/param_sets/WOTH_M1/3_28_22/WOTH_M1_S2.rds")
for (each_var in 1:ncol(param_set)){
  #Random
  assign(colnames(param_set)[each_var],param_set[sample(1:nrow(param_set),1),each_var])
  #Get best
  #assign(colnames(param_set)[each_var],param_set[which(param_set$sum_err == min(param_set$sum_err)),each_var])
}

#Define species
species_target <- "Wood Thrush"

#Set number of birds to include in simulation
num_pulls <- 2000

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

start_date <- (ebirdst_runs[which(ebirdst_runs$common_name==species_target),8:9])
start_date <- as.numeric(start_date[1] + (start_date[2]-start_date[1])/2)
start_date <- yday(as.Date(start_date,origin = "1970-01-01"))
acceptable_week_dates <- seq(1,365,by=7)
start_date <- acceptable_week_dates[which(abs(acceptable_week_dates - start_date) == min(abs(acceptable_week_dates - start_date)))]
print(start_date)

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
  summary(upd_df)
  locs <- save_locations(upd_df)
  loc_ind_s <- (days - 1)*num_pulls + 1 
  loc_ind_e <- (days - 1)*num_pulls + num_pulls 
  locs_rec[loc_ind_s:loc_ind_e,] <- locs
}


#Get map
myMap <- get_stamenmap(bbox = c(left = -165,
                                bottom = 30,
                                right = -100,
                                top = 72),
                       maptype = "terrain-background",
                       color="bw",
                       crop = TRUE,
                       zoom = 2,
                       force = TRUE,
) 

countries <- map_data("world")
states <- map_data("state")
countries <- filter(countries,region %in% c("Canada","USA","Mexico","Guatemala","Bahamas",
                                            "Belize","Honduras","El Salvador","Nicaragua",
                                            "Costa Rica","Panama","Cuba","Jamaica","Haiti",
                                            "Dominican Republic","Puerto Rico","Colombia",
                                            "Venezuela","Brazil","Guyana","Suriname","French Guiana"))
#p <- ggplot() + coord_map("ortho", xlim=c(-145,-80),ylim=c(20,75)) +
#  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="lightgrey", color="darkgrey")
#p

p <- ggplot() + coord_map("mollweide",xlim=c(-140,-65),ylim=c(15,70)) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  
  theme(panel.background = element_rect(fill = alpha("lightskyblue2",.5)))
p
#Plot with points
#Location of birds on a single day
#loc_day <- filter(locs_rec,Day == 250)
#mappoints <- p + 
#   geom_point(data = loc_day,aes(x=Lon,y=Lat),color=alpha("black",.8),size=.5,pch=19)
#plot(mappoints)

#Animate

head(locs_rec)

days_to_inc <- seq(from=1,to=365, by=3)

days_to_inc

loc_day <- filter(locs_rec, Day %in% days_to_inc)

head(loc_day)

pp <- p +
  geom_point(data = loc_day,aes(x=Lon,y=Lat,group=ID),color=alpha("tan3",.6),size=1.25,pch=19)
                 

anim <- pp + 
  transition_states(Day,transition_length = 1,state_length = 1) #+ ggtitle('Julian Day: {closest_state}')
                    #subtitle = 'Frame {frame} of {nframes}')

#This will display the plot in the viewer tab of rstudio

#the fps should put the gif in the ~13 second range, similar to eBird S&T, but
# for exact matching a third-party software will need to be used.
animate(anim, nframes = 2*length(days_to_inc),fps=20,height=800,width=800)


#Plot every location for every bird
#p <- ggplot() + coord_map("mollweide",xlim=c(-165,-110),ylim=c(30,75)) +
#  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="lightgrey", color="darkgrey")
#p

