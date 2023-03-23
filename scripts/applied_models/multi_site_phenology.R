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

setwd("~/Documents/Coding/R/MIMAS")

#Get functions
source(file = "scripts/sim_functions/sim_funct.R")

# Parameter sets should be drawn from a script in the data folder, or from a 
# fit parameter set
#Script
#source(file = "data/param_sets/Clay_colored_Sparrow_BG_3_1_22.R")
#Fit Data
cumulative_breeding_rec_lst <- list()
cumulative_passage_rec_lst <- list()

multi_sim_breeder_nb_locs <- list()
multi_sim_migrants_b_locs <- list()

multi_sim_all_nb <- list()
multi_sim_all_b <- list()


birds_by_season_rec <- list()

cell_size <- dgconstruct(spacing=300,metric=TRUE) 

#Enter city coordinates
city_name1 <- "Mexico City"
p1 <- c(-99.13,19.43)
targ_cell1 <- dgGEO_to_SEQNUM(cell_size, p1[1], p1[2])$seqnum

city_name2 <- "Los Angeles"
p2 <- c(-118.24,34.07)
targ_cell2 <- dgGEO_to_SEQNUM(cell_size, p2[1], p2[2])$seqnum

city_name3 <- "Seattle"
p3 <- c(-122.33,47.61)
targ_cell3 <- dgGEO_to_SEQNUM(cell_size, p3[1], p3[2])$seqnum

all_cell_nums <- c(targ_cell1,targ_cell2,targ_cell3)

p1_spring_passage <- c()
p2_spring_passage <- c()
p3_spring_passage <- c()

p1_fall_passage <- c()
p2_fall_passage <- c()
p3_fall_passage <- c()

#Define species
species_target <- "Townsend's Warbler"
param_set <- readRDS("data/param_sets/TOWA/9_20_22/TOWA_S2.rds")

#Import breeding, non-breeding maps
breeding_file <- read.csv(paste("data/species_maps/breeding/",gsub(" ","_",species_target),"_breeding_all.csv",sep=""))
nonbreeding_file <- read.csv(paste("data/species_maps/non_breeding/",gsub(" ","_",species_target),"_non_breeding_all.csv",sep=""))

#Start date
start_date <- (ebirdst_runs[which(ebirdst_runs$common_name==species_target),7:8])
start_date <- as.numeric(start_date[1] + (start_date[2]-start_date[1])/2)
start_date <- yday(as.Date(start_date,origin = "1970-01-01"))
acceptable_week_dates <- seq(1,365,by=7)
start_date <- acceptable_week_dates[which(abs(acceptable_week_dates - start_date) == min(abs(acceptable_week_dates - start_date)))]

for (each_sim in 1:100){
  print(paste("working on simulation number:",each_sim))
  
  sample_num <- sample(1:nrow(param_set),1)
  for (each_var in 1:ncol(param_set)){
    #Random
    assign(colnames(param_set)[each_var],param_set[sample_num,each_var])
    #Get best
    #assign(colnames(param_set)[each_var],param_set[which(param_set$sum_err == min(param_set$sum_err)),each_var])
  }
  
  #Set number of birds to include in simulation
  num_pulls <- 2000
  
  #Define spacing of error maps
  spacing <- 200
  
  #Define starting season
  season <- 1
  
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
    upd_df <- run_day(upd_df,static_df)
    summary(upd_df)
    locs <- save_locations(upd_df,save_season = TRUE)
    loc_ind_s <- (days - 1)*num_pulls + 1 
    loc_ind_e <- (days - 1)*num_pulls + num_pulls 
    locs_rec[loc_ind_s:loc_ind_e,] <- locs
  }
  
  # Get spring passage dates
  spring_migrant_records <- locs_rec %>% 
    filter(season == 4)
  spring_migrant_records$cell <- dgGEO_to_SEQNUM(cell_size, spring_migrant_records$Lon, spring_migrant_records$Lat)$seqnum
  
  spring_migrant_records <- spring_migrant_records %>% 
    filter(cell %in% all_cell_nums)
  
  targ_cell1_days<- as.numeric(unlist(spring_migrant_records %>% filter(cell == targ_cell1) %>% select(Day)))
  targ_cell2_days<- as.numeric(unlist(spring_migrant_records %>% filter(cell == targ_cell2) %>% select(Day)))
  targ_cell3_days<- as.numeric(unlist(spring_migrant_records %>% filter(cell == targ_cell3) %>% select(Day)))
    
  p1_spring_passage <-c(p1_spring_passage,targ_cell1_days)
  p2_spring_passage <-c(p2_spring_passage,targ_cell2_days)
  p3_spring_passage <-c(p3_spring_passage,targ_cell3_days)
  
  # Get fall passage dates
  fall_migrant_records <- locs_rec %>% 
    filter(season == 2)
  fall_migrant_records$cell <- dgGEO_to_SEQNUM(cell_size, fall_migrant_records$Lon, fall_migrant_records$Lat)$seqnum
  
  fall_migrant_records <- fall_migrant_records %>% 
    filter(cell %in% all_cell_nums)
  
  targ_cell1_days<- as.numeric(unlist(fall_migrant_records %>% filter(cell == targ_cell1) %>% select(Day)))
  targ_cell2_days<- as.numeric(unlist(fall_migrant_records %>% filter(cell == targ_cell2) %>% select(Day)))
  targ_cell3_days<- as.numeric(unlist(fall_migrant_records %>% filter(cell == targ_cell3) %>% select(Day)))
  
  p1_fall_passage <-c(p1_fall_passage,targ_cell1_days)
  p2_fall_passage <-c(p2_fall_passage,targ_cell2_days)
  p3_fall_passage <-c(p3_fall_passage,targ_cell3_days)
}

spring_passage <- list(p1_spring_passage,p2_spring_passage,p3_spring_passage)
saveRDS(spring_passage,"spring_passage_TOWA_2_27_23.rds")
hist(p1_spring_passage)
hist(p2_spring_passage)
hist(p3_spring_passage)

fall_passage <- list(p1_fall_passage,p2_fall_passage,p3_fall_passage)
saveRDS(fall_passage,"fall_passage_TOWA_2_27_23.rds")
hist(p1_fall_passage)
hist(p2_fall_passage)
hist(p3_fall_passage)

### Read in file ###
spring_passage <- readRDS("spring_passage_TOWA_2_27_23.rds")
p1_spring_passage <- spring_passage[[1]]
p2_spring_passage <- spring_passage[[2]]
p3_spring_passage <- spring_passage[[3]]

fall_passage <- readRDS("fall_passage_TOWA_2_27_23.rds")
p1_fall_passage <- fall_passage[[1]]
p2_fall_passage <- fall_passage[[2]]
p3_fall_passage <- fall_passage[[3]]
### Make plots showing spring passage rates ####

targ_cells_all <- data.frame(cell=c(targ_cell1,targ_cell2,targ_cell3))
grid <- dgcellstogrid(cell_size, targ_cells_all$cell, frame=TRUE)
cellcenters <- as.data.frame(dgSEQNUM_to_GEO(cell_size,targ_cells_all$cell))

cellcenters$name <- c("A","B","C")
#Get non-breeding points
#Set number of birds to include in simulation
num_pulls <- 10000

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

countries <- map_data("world")
states <- map_data("state")
countries <- filter(countries,region %in% c("Canada","USA","Mexico","Guatemala","Bahamas",
                                            "Belize","Honduras","El Salvador","Nicaragua",
                                            "Costa Rica","Panama","Cuba","Jamaica","Haiti",
                                            "Dominican Republic","Puerto Rico","Colombia",
                                            "Venezuela","Brazil","Guyana","Suriname","French Guiana",
                                            "Antigua","Barbados","Martinique","Montserrat",
                                            "Anguilla","Aruba","Curacao","Dominica","Cayman Islands",
                                            "Turks and Caicos Islands","Bermuda","Saint Lucia","Grenada","Virgin Islands"))
country_fill_color <- "grey70"
country_border_color <- "grey50"
p <- ggplot() + coord_map("mollweide",xlim=c(-125,-85),ylim=c(8,52)) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=country_fill_color, color=country_border_color) +
  #geom_polygon(data=states, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
  geom_point(data=static_df,aes(x=nonbreeding_lon,y=nonbreeding_lat),alpha=.15,size=.5,color="dodgerblue3",pch=19) +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group), alpha=0.3,col="orchid4",size=1.25)    +
  #geom_text(aes(x=cellcenters$lon_deg,y=cellcenters$lon_deg,label=cellcenters$name),stat="unique",alpha=.8,col="black",size=1.5) +
  geom_text(aes(x = cellcenters$lon_deg[1]+3.5, y = cellcenters$lat_deg[1]+.5,label = "1"),stat = "unique",size=8,color="black") +
  geom_text(aes(x = cellcenters$lon_deg[2]-3, y = cellcenters$lat_deg[2],label = "2"),stat = "unique",size=8,color="black") +
  geom_text(aes(x = cellcenters$lon_deg[3]-4, y = cellcenters$lat_deg[3],label = "3"),stat = "unique",size=8,color="black") +
  #geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  #geom_polygon(data=grid3,      aes(x=long, y=lat, group=group),alpha=1,color="black")    +
  #geom_path   (data=grid3,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  #scale_fill_gradient2(low="dodgerblue4", mid="grey",high="firebrick3",limits=c(-1,max(grid$diff))) +
  #labs(fill = "Prob.") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  theme(panel.background = element_rect(fill = alpha("white",.5)))
p

par(mar=c(5,5,5,5))
par(las=1)
par(pty="s")
size_labels <- 1.5
size_axis <- 1.5
plot(density(p1_spring_passage,bw=2),xlim=c(70,170),ylim=c(0,.045),
     xlab="Day of Year",ylab="",main="",
     lwd=5,col=alpha("forestgreen",.95),cex.lab=size_labels,cex.axis = size_axis)
points(density(p2_spring_passage,bw=2),type="l",lty=2,lwd=5,col=alpha("forestgreen",.95))
points(density(p3_spring_passage,bw=2),type="l",lty=3,lwd=5,col=alpha("forestgreen",.95))
title(ylab="Relative Migration Intensity", line=4.25, cex.lab=size_labels)

text(median(p1_spring_passage),.042,"1",cex=2.5)
text(median(p2_spring_passage),.04,"2",cex=2.5)
text(median(p3_spring_passage),.038,"3",cex=2.5)

plot(density(p1_fall_passage,bw=2),xlim=c(210,310),ylim=c(0,.031),
     xlab="Day of Year",ylab="",main="",
     lwd=5,col=alpha("tan3",.95),cex.lab=size_labels,cex.axis=size_axis)
points(density(p2_fall_passage,bw=2),type="l",lty=2,lwd=5,col=alpha("tan3",.95))
points(density(p3_fall_passage,bw=2),type="l",lty=3,lwd=5,col=alpha("tan3",.95))
title(ylab="Relative Migration Intensity", line=5, cex.lab=size_labels)
text(median(p1_fall_passage),.026,"1",cex=2.5)
text(median(p2_fall_passage),.027,"2",cex=2.5)
text(median(p3_fall_passage),.028,"3",cex=2.5)

