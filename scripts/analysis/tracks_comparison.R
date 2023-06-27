# Script to compare tracking data to MIMAS-estimated migratory paths.

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
library(oce)
sf::sf_use_s2(FALSE)

setwd("~/Documents/Coding/R/MIMAS/")

#Get functions
source(file = "scripts/sim_functions/sim_funct.R")

#Define species
species_target <- "Wood Thrush"

#Set number of birds to include in simulation
num_pulls <- 200

#Define spacing of error maps
spacing <- 200

#Define starting season
season <- 1

# Read in geolocator tracks, convert locations to lat longs
# To recreate, these need to be downloaded from: https://esajournals.onlinelibrary.wiley.com/doi/pdf/10.1002/ecs2.3421
emp_tracks <- read_csv("data/GPS_tracks/Stanley_2021_GPSLocations.csv")
emp_tracks_ll <- utm2lonlat(easting=emp_tracks$utmEasting,
           northing=emp_tracks$utmNorthing,
           zone=as.numeric(gsub("N","",emp_tracks$utmZone)),
           hemisphere = "N")
emp_tracks$latitude <- emp_tracks_ll$latitude
emp_tracks$longitude <- emp_tracks_ll$longitude
emp_tracks_ll <- NULL

#plot(emp_tracks$longitude,emp_tracks$latitude)

#Get a look at empirically tracked breeding locations
cell_size <- dgconstruct(spacing=50,metric=TRUE) 

breed_emp_tracks <- filter(emp_tracks,stage=="breeding") %>% dplyr::select(tagNo,latitude,longitude)

breed_emp_tracks$cell <- dgGEO_to_SEQNUM(cell_size,breed_emp_tracks$longitude,breed_emp_tracks$latitude)$seqnum

table(breed_emp_tracks$tagNo,breed_emp_tracks$cell)

#Get the cell locations of each of the breed
cell_breeding_locs <- unique(breed_emp_tracks$cell)
cell_breeding_locs <- cell_breeding_locs[order(cell_breeding_locs)]
col_polys <- c("forestgreen","orchid4","dodgerblue3","firebrick3","orange2")

grid <- dgcellstogrid(cell_size, cell_breeding_locs, frame=TRUE)
#grid$cell <- as.numeric(grid$cell)

countries <- map_data("world")
countries <- countries %>% filter(long < -32)

p <- ggplot() + coord_map("mollweide",xlim=c(-100,-70),ylim=c(10,48)) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="grey70", color="grey50") +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group),col=col_polys[as.factor(grid$cell)],fill="transparent", alpha=1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  
  theme(panel.background = element_rect(fill = alpha("white",.5)))
p

#Read in breeding, wintering locations
#Import breeding, non-breeding maps
breeding_file <- read.csv(paste("data/species_maps/breeding/",gsub(" ","_",species_target),"_breeding_all.csv",sep=""))
nonbreeding_file <- read.csv(paste("data/species_maps/non_breeding/",gsub(" ","_",species_target),"_non_breeding_all.csv",sep=""))

#Here, we can "hack" the breeding file to only include cells where tracking data is available
breeding_file$cell <- dgGEO_to_SEQNUM(cell_size,breeding_file$longitude,breeding_file$latitude)$seqnum

#now remove cells
breeding_file <- breeding_file %>% filter(cell %in% cell_breeding_locs)

#Weight cells draws equally
cell_count <- breeding_file %>% group_by(cell) %>% summarise(count=n())
cell_count$prob <- 1/cell_count$count

breeding_file <- merge(breeding_file,cell_count,by="cell")
breeding_file$perc_occupancy <- breeding_file$prob
breeding_file <- dplyr::select(breeding_file,-c("prob","cell"))

tracked_locs_all <- c()
for (nn in 1:100){
  print(nn)
  # Parameter sets should be drawn from a script in the data folder, or from a 
  # fit parameter set
  #Script
  #source(file = "data/param_sets/Clay_colored_Sparrow_BG_3_1_22.R")
  #Fit Data
  param_set <- readRDS("data/param_sets/WOTH_M1/9_20_22/WOTH_M1_S2.rds")
  for (each_var in 1:ncol(param_set)){
    #Random
    assign(colnames(param_set)[each_var],param_set[sample(1:nrow(param_set),1),each_var])
    #Get best
    #assign(colnames(param_set)[each_var],param_set[which(param_set$sum_err == min(param_set$sum_err)),each_var])
  }
  
  start_date <- (ebirdst_runs[which(ebirdst_runs$common_name==species_target),7:8])
  start_date <- as.numeric(start_date[1] + (start_date[2]-start_date[1])/2)
  start_date <- yday(as.Date(start_date,origin = "1970-01-01"))
  acceptable_week_dates <- seq(1,365,by=7)
  start_date <- acceptable_week_dates[which(abs(acceptable_week_dates - start_date) == min(abs(acceptable_week_dates - start_date)))]
  #print(start_date)
  
  #Initialize model
  initial_dfs <- initialize_MIBM(num_pulls,breeding_file,nonbreeding_file,speed_mean_s,
                                 speed_sd_s,speed_mean_f,speed_sd_f,start_date_u_s,
                                 start_date_sd_s,start_date_u_f,start_date_sd_f,
                                 max_mig_time_s,max_mig_time_f,bear_err_mean_s,
                                 bear_err_sd_s,bear_err_mean_f,bear_err_sd_f,
                                 max_energy_s,max_energy_f,recovery_rate_s,
                                 recovery_rate_f,season,start_date,goal_radius,
                                 mig_con,migr_timing_lat_s,migr_timing_lat_f)
  static_df <- initial_dfs[[1]]
  upd_df <- initial_dfs[[2]]
  
  locs_rec <- as.data.frame(matrix(NA,ncol = 5,nrow = (365*num_pulls)))
  colnames(locs_rec) <- c("ID","Lon","Lat","Day","season")
  for (days in 1:365){
    #print(days)
    upd_df <- run_day(upd_df,static_df)
    summary(upd_df)
    locs <- save_locations(upd_df,save_season = TRUE)
    loc_ind_s <- (days - 1)*num_pulls + 1 
    loc_ind_e <- (days - 1)*num_pulls + num_pulls 
    locs_rec[loc_ind_s:loc_ind_e,] <- locs
  }
  
  #Find birds that breed in cells where geolocators were deployed
  s_df_tracks <- static_df
  s_df_tracks$cell <- dgGEO_to_SEQNUM(cell_size,s_df_tracks$breeding_lon,s_df_tracks$breeding_lat)$seqnum
  
  s_df_tracks <- s_df_tracks %>% filter(cell %in% cell_breeding_locs) %>% group_by(cell)
  
  table(s_df_tracks$cell)
  
  tracked_sim_id <- s_df_tracks$ID
  
  nrow(s_df_tracks)
  
  #Reduce to a sample of 10 birds from each cell (50 total)
  tracked_locs <- locs_rec %>% filter(ID %in% tracked_sim_id & season == 2)
  
  tracked_locs <- merge(tracked_locs,dplyr::select(s_df_tracks,c(ID,cell)),by="ID")
  id_keep_list <- tracked_locs %>% dplyr::select(ID,cell) %>% unique() %>% group_by(cell) %>% reframe(id_samp = ID[1:10])
  
  tracked_locs <- filter(tracked_locs,ID %in% id_keep_list$id_samp)
  cell_breeding_locs_col <- data.frame(cell = cell_breeding_locs,col=col_polys)
  
  tracked_locs_day <- merge(tracked_locs,cell_breeding_locs_col,by="cell")  %>% unique()
  tracked_locs <- merge(tracked_locs,cell_breeding_locs_col,by="cell") %>% select(-c(Day)) %>% unique()
  
  
  #p <- p + geom_point(data=tracked_locs,aes(x=Lon,y=Lat),col=tracked_locs$col,size=.25,alpha=.1)
  #p2
  tracked_locs_all <- rbind(tracked_locs_all,tracked_locs)
}

#saveRDS(tracked_locs_all,"tracked_locs_all_WOTH.rds")
tracked_locs_all <- readRDS("tracked_locs_all_WOTH.rds")
cell_breeding_locs_col <- data.frame(cell = cell_breeding_locs,col=col_polys)
tag_col <- merge(breed_emp_tracks,cell_breeding_locs_col,by="cell") %>% dplyr::select(c(tagNo,col)) %>% unique()

#Now plot empirically-derived migratory paths - Fall
fall_emp_tracks <- emp_tracks %>% filter(stage %in% c("autumnMigration")) %>% group_by(tagNo) %>% merge(tag_col)
breed_start <- breed_emp_tracks %>% dplyr::select(longitude,latitude) %>% summarise(longitude=round(longitude,0),latitude=round(latitude,0)) %>% unique()

pp1 <- ggplot() + coord_map("mollweide",xlim=c(-100,-70),ylim=c(10,48)) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="grey70", color="grey50") +
  geom_point(data=breed_start,aes(x=longitude,y=latitude),pch=21,stroke=2,
             col=col_polys[as.factor(c(3,5,2,4,1))],fill="black",
             size=15,alpha=.9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  theme(panel.background = element_rect(fill = alpha("white",.5)))
pp2 <- pp1 + 
  geom_point(data=tracked_locs_all,aes(x=Lon,y=Lat),col=tracked_locs_all$col,size=1,alpha=.15) +
  geom_point(data=fall_emp_tracks,aes(longitude, latitude,group=tagNo),size=7.5,col="black",pch=25,fill=fall_emp_tracks$col,alpha=.8) 
pp2

library(ks)
# Estimate coverage (% geolocator points within range)
par(mfrow=c(3,2))
par(mar=c(2.5,0,2.5,0))
par(pty="s")
for (each_col in unique(tracked_locs_all$col)){
  point_mat <- tracked_locs_all %>% filter(col==each_col) %>% dplyr::select(c("Lon","Lat")) %>% as.matrix()
  geolocator_points <- fall_emp_tracks %>% filter(col==each_col) %>% dplyr::select(c("longitude","latitude")) %>% as.matrix()
  
  z <- kde2d(point_mat[,1], point_mat[,2], n=100,h = 1)
  confidencebound <- quantile(z$z, probs= 0.90)
  
  plot(point_mat, xlab="", ylab="", pch=19, cex=.1,col = alpha(each_col,.3),xlim=c(-105,-70),ylim=c(10,45),cex.axis=1.65)
  contour(z, levels = confidencebound, col="black", add = TRUE,labels = "",method="simple",lwd=2)
  points(geolocator_points,pch=25,col="black",cex=1.75,bg="black") 
  
  #smoothScatter(point_mat,colramp=colorRampPalette(c("white", each_col)),xlim=c(-110,-75),ylim=c(10,45),nbin=100)
  #points(geolocator_points,pch=19,col="black",cex=.5) 
}

# 87 of 94 points lie within the 90% KDE
