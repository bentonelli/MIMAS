#Model to visualize error spatially and temporally
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

### Run a single simulation ####

#Number of simulations to run
num_simulations <- 1

#Number of birds
num_pulls <- 5000

#Species name
species_target <- "Wood Thrush"

#Define spacing of error grid
spacing <- 200

#Define starting season
season <- 1

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
  param_set <- readRDS("data/fit_param_sets/WOTH_M1/Step_2/WOTH_M1_3_29_22.rds")
  #For best model
  samp_ind <- which(param_set$sum_err == min(param_set$sum_err))
  #samp_ind <- sample(1:nrow(param_set),1)
  for (each_var in 1:ncol(param_set)){
    assign(colnames(param_set)[each_var],param_set[samp_ind,each_var])
  }
  
  
  #source(file = "data/param_sets/Clay_colored_Sparrow_3_2_22.R")
  
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

### Create map of error ####
library(ggplot2)
library(sf)
library(rnaturalearth)
library(ggpubr)
#Basic temporal

plot(NULL,xlim=c(0,52),ylim=c(0,40))
for (n in 1:nrow(err_record)){
  points(as.numeric(err_record[n,28:79]),type="l",col="blue")
}

### Basic spatial, for single week ####
week_err <- as.data.frame(error_run[[3]][1])
cell_size <- dgconstruct(spacing=200,metric=TRUE) 
grid <- dgcellstogrid(cell_size, week_err$cell, frame=TRUE)
grid$cell <- as.numeric(grid$cell)
grid <- left_join(week_err,grid,by="cell")

worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf')
world_cropped <- st_crop(worldmap, xmin = -179, xmax = -20,
                          ymin = -20, ymax = 85)
p<- ggplot() + geom_sf(data = world_cropped) +
  coord_sf(xlim = c(-130, -55), ylim = c(10, 60), expand = FALSE) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=perc_pop), alpha=0.4)    +
  geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradient(low="blue", high="red")+ theme_bw()
p

### For multiple weeks ####

#First get the max error across all weeks to get a value for the color gradient
max_err <- 0
for (m in 1:52){
  week_err <- as.data.frame(error_run[[3]][m])
  max_week <- max(abs(week_err$count - week_err$perc_pop))
  total_week_err <- sum(abs(week_err$count - week_err$perc_pop))*num_pulls
  if (max_week > max_err){
    max_err <- max_week
  }
}
week_cumulative_err <- as.data.frame(seq(1,52,by=1))
week_cumulative_err$err <- as.numeric(error_run[[1]])
colnames(week_cumulative_err)[1] <- "week"

pdf(file = "spatial_err_analysis.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
for (n in 1:52){
  week_err <- as.data.frame(error_run[[3]][n])
  cell_size <- dgconstruct(spacing=200,metric=TRUE) 
  grid <- dgcellstogrid(cell_size, week_err$cell, frame=TRUE)
  grid$cell <- as.numeric(grid$cell) 
  grid <- left_join(week_err,grid,by="cell")
  
  #Get error for each cell, corrected 
  
  worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                           returnclass = 'sf')
  world_cropped <- st_crop(worldmap, xmin = -179, xmax = -20,
                           ymin = -20, ymax = 85)
  p<- ggplot() + geom_sf(data = world_cropped) +
    coord_sf(xlim = c(-130, -55), ylim = c(10, 60), expand = FALSE) +
    geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
    geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=count-perc_pop), alpha=0.4)    +
    geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
    #ggtitle(round(sum(abs(week_err$count - week_err$perc_pop)),4)) +
    scale_fill_gradient2(low="dodgerblue4",mid = "white", high="firebrick4",midpoint=0,limits=c(-max_err,max_err))+ theme_bw()
  week_cumulative_err$cols_group <- as.factor(c(rep(1,n),rep(2,52-n)))
  time_plot <- ggplot(data=week_cumulative_err, aes(x=week, y=err,color=cols_group,size=3)) +
    xlim(0,52) +
    ylim(0,20) +
    geom_line()+
    geom_vline(xintercept = n, linetype="solid", 
               color = "red", size=3) +
    theme_bw() + theme(panel.grid.major = element_blank(), legend.position = "none",panel.grid.minor = element_blank()) +
    ylab("Error") + xlab("Week") + scale_color_manual(values=c("red","black"))
  #time_plot
  print(ggarrange(p, time_plot,heights = c(2, 0.7),
            ncol = 1, nrow = 2))
}
dev.off()


