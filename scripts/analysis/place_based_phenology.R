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
  #assign(colnames(param_set)[each_var],param_set[sample(1:nrow(param_set),1),each_var])
  #Get best
  assign(colnames(param_set)[each_var],param_set[which(param_set$sum_err == min(param_set$sum_err)),each_var])
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
#Enter city coordianates
#Atlanta
city_name <- "Atlanta"
p2 <- c(-84.4,33.75)
#New York
#city_name <- "New York"
#p2 <- c(-74,40.75)
#Chicago
#city_name <- "Chicago"
#p2 <- c(-87.7,41.84)

locs_rec_interest <- filter(locs_rec,Lon < p2[1]+2 & Lon >p2[1]-2 & Lat > p2[2]-2 & Lat < p2[2]+2)
to_check <- as.matrix(locs_rec_interest[,2:3])
dist_list <- distHaversine(to_check,p2)/1000
locs_rec_interest <- locs_rec_interest[which(dist_list <= 100),]
summary(locs_rec_interest)

birds_by_season <- locs_rec_interest %>%
  group_by(Day) %>%
  summarise(breeding_count = sum(season==1),
            fall_migrant_count = sum(season==2),
            spring_migrant_count = sum(season==4))
birds_by_season$perc_migrant <- 1 - (birds_by_season$breeding_count/
  (birds_by_season$breeding_count+birds_by_season$fall_migrant_count + birds_by_season$spring_migrant_count))
birds_by_season$total_birds <- birds_by_season$breeding_count+birds_by_season$fall_migrant_count+birds_by_season$spring_migrant_count


#Smooth, buffer
library(zoo)
day_merge <- as.data.frame(1:365)
colnames(day_merge) <- "Day"
birds_by_season <- merge(day_merge,birds_by_season,by="Day",all=TRUE)

birds_by_season$breeding_count[which(is.na(birds_by_season$breeding_count))] <- 0
birds_by_season$fall_migrant_count[which(is.na(birds_by_season$fall_migrant_count))] <- 0
birds_by_season$spring_migrant_count[which(is.na(birds_by_season$spring_migrant_count))] <- 0


birds_by_season$breeding_count_roll <- rollmean(birds_by_season$breeding_count,5,fill = NA)
birds_by_season$spring_migrant_count_roll <- rollmean(birds_by_season$spring_migrant_count,5,fill = NA)
birds_by_season$fall_migrant_count_roll <- rollmean(birds_by_season$fall_migrant_count,5,fill = NA)

par(pty="s")
par(las=1)
#plot(birds_by_season$Day,birds_by_season$perc_migrant)
plot(birds_by_season$Day,birds_by_season$breeding_count_roll,type="l",
     col=alpha("orchid",.7),xlim=c(80,170),
     ylim=c(0,max(birds_by_season$total_birds,na.rm = TRUE)*.9),
     cex.lab=1.4,cex.axis=1.2,cex.main=1.4,
     lwd=10,ylab="Number of Birds",xlab="Julian Day",main=city_name)
points(birds_by_season$Day,birds_by_season$spring_migrant_count_roll,type = "l",col=alpha("forestgreen",.9),lwd=10,lty=1)
legend("topright",legend=c("Migrants","Local Breeders"),lty=1,lwd=10,col=c("forestgreen","orchid"))
#points(birds_by_season$Day,birds_by_season$fall_migrant_count_roll,type = "l",col=alpha("orange3",.9),lwd=6,lty=2)

#Plot map
#Save boxes
boxes1<-data.frame(maxlat = 35.75,minlat = 31.75,maxlong = -82.4,minlong = -86.4, id="1")
boxes1<-transform(boxes1, laby=(maxlat +minlat )/2, labx=(maxlong+minlong )/2)

boxes2<-data.frame(maxlat = 42.75,minlat = 38.75,maxlong = -72,minlong = -76, id="1")
boxes2<-transform(boxes2, laby=(maxlat +minlat )/2, labx=(maxlong+minlong )/2)

boxes3<-data.frame(maxlat = 43.84,minlat = 39.84,maxlong = -85.7,minlong = -89.5, id="1")
boxes3<-transform(boxes3, laby=(maxlat +minlat )/2, labx=(maxlong+minlong )/2)


worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf')
world_cropped <- st_crop(worldmap, xmin = -179, xmax = -20,
                         ymin = -20, ymax = 85)
p<- ggplot() + geom_sf(data = world_cropped) +
  coord_sf(xlim = c(-105, -65), ylim = c(20, 50), expand = FALSE) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="lightgrey", color="darkgrey") +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill="lightgrey", color="darkgrey") +
  scale_fill_gradient(low="blue", high="red")+ theme_bw() + 
  geom_rect(data=boxes1, aes(xmin=minlong , xmax=maxlong, ymin=minlat, ymax=maxlat ), color="black", fill=alpha("grey99",.8)) + 
  geom_rect(data=boxes2, aes(xmin=minlong , xmax=maxlong, ymin=minlat, ymax=maxlat ), color="black", fill=alpha("grey99",.8)) + 
  geom_rect(data=boxes3, aes(xmin=minlong , xmax=maxlong, ymin=minlat, ymax=maxlat ), color="black", fill=alpha("grey99",.8)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank()) +
  theme(panel.background = element_rect(fill = alpha("lightskyblue2",.5)))
p

