#Phenology

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
library(ggplot2)
library(rnaturalearth)
library(ggpubr)
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
num_pulls <- 5000

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
locs_rec <- as.data.frame(matrix(NA,ncol = 5,nrow = (365*num_pulls)))
colnames(locs_rec) <- c("ID","Lon","Lat","Day","season")
for (days in 1:365){
  print(days)
  upd_df <- run_day(upd_df,static_df)
  summary(upd_df)
  locs <- save_locations(upd_df,save_season = TRUE)
  loc_ind_s <- (days - 1)*num_pulls + 1 
  loc_ind_e <- (days - 1)*num_pulls + num_pulls 
  locs_rec[loc_ind_s:loc_ind_e,] <- locs
}

#Find birds in Yucatan, only during springtime - 50: 200
yuc_records <- filter(locs_rec, Lon > -92 & Lon < -85 & Lat < 24 & Lat > 19,Day < 200 & Day > 50)
nrow(yuc_records)

#Save box for later
boxes<-data.frame(maxlat = 24,minlat = 19,maxlong = -85,minlong = -92, id="1")
boxes<-transform(boxes, laby=(maxlat +minlat )/2, labx=(maxlong+minlong )/2)

#Simulate yes or no - probability based on nothing!
yuc_records$tick_attach <- rbinom(nrow(yuc_records),size = 1,prob = .50)
sum(yuc_records$tick_attach)

#Simulate detach date
yuc_records$tick_detach_date <- (yuc_records$Day + rpois(nrow(yuc_records),1)) * yuc_records$tick_attach

detach_ticks <- filter(yuc_records,tick_detach_date>0)
detach_ticks$Lon_drop_off <- NA
detach_ticks$Lat_drop_off <- NA

tick_record <- c()
for (each_day in unique(detach_ticks$tick_detach_date)){
  
  start_locations <- detach_ticks %>%
    filter(tick_detach_date == each_day) %>%
    select(ID,Lon,Lat)
  colnames(start_locations) <- c("ID","Lon_start","Lat_start")
  id_birds <- as.numeric(unlist(start_locations$ID))
  
  bird_w_tick_num <- as.data.frame(table(id_birds))
  colnames(bird_w_tick_num) <- c("ID","tick_num")
  
  bird_locations <- filter(locs_rec,Day == each_day & ID %in% id_birds)
  
  tick_locations <- merge(bird_locations,bird_w_tick_num,by="ID")
  tick_locations <- merge(tick_locations,start_locations,by="ID")
  
  diff_loc <- tick_locations$Lon - tick_locations$Lon_start
  
  tick_locations <- select(tick_locations,Lon,Lat,Day,season,tick_num)
  tick_locations <- tick_locations[which(diff_loc!=0),]
  tick_record <- rbind(tick_record,tick_locations)
}

countries <- map_data("world")
states <- map_data("state")
countries <- filter(countries,region %in% c("Canada","USA","Mexico","Guatemala","Bahamas",
                                            "Belize","Honduras","El Salvador","Nicaragua",
                                            "Costa Rica","Panama","Cuba","Jamaica","Haiti",
                                            "Dominican Republic","Puerto Rico","Colombia",
                                            "Venezuela","Brazil","Guyana","Suriname","French Guiana"))
p <- ggplot() + coord_map("mollweide",xlim=c(-110,-75),ylim=c(10,45)) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="lightgrey", color="darkgrey") +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill="lightgrey", color="darkgrey") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank()) +
  theme(panel.background = element_rect(fill = alpha("lightskyblue2",.5)))
#Plot with points
mappoints <- p + 
  geom_point(data = tick_record,aes(x=Lon,y=Lat),color=alpha("black",.8),size=.5,pch=19)
plot(mappoints)


cell_size <- dgconstruct(spacing=200,metric=TRUE) 
tick_record$cell <- dgGEO_to_SEQNUM(cell_size, tick_record$Lon, tick_record$Lat)$seqnum
tick_record_grp <- tick_record %>%
  group_by(cell) %>%
  summarise(total_tick_count = sum(tick_num))


grid <- dgcellstogrid(cell_size, tick_record_grp$cell, frame=TRUE)
grid$cell <- as.numeric(grid$cell)
grid <- left_join(tick_record_grp,grid,by="cell")

worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf')
world_cropped <- st_crop(worldmap, xmin = -179, xmax = -20,
                         ymin = -20, ymax = 85)
p<- ggplot() + geom_sf(data = world_cropped) +
  coord_sf(xlim = c(-115, -70), ylim = c(12, 50), expand = FALSE) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="lightgrey", color="darkgrey") +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill="lightgrey", color="darkgrey") +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=total_tick_count), alpha=0.4)    +
  geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradient(low="blue", high="red")+ theme_bw() + 
  geom_rect(data=boxes, aes(xmin=minlong , xmax=maxlong, ymin=minlat, ymax=maxlat ), color="red", fill="transparent") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank()) +
  theme(panel.background = element_rect(fill = alpha("lightskyblue2",.5)))
p
