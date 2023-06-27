# Script to look at areas, times where eBird ST and MIMAS consistently disagree


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
library(viridis)
sf::sf_use_s2(FALSE)

### Run a single simulation ####

#Number of simulations to run
num_simulations <- 100

#Number of birds
num_pulls <- 2000

#Species name
species_target <- "Wood Thrush"

#Define spacing of error grid
spacing <- 200

#Define starting season
season <- 1

#This is the Julian date of the start date for simulations
start_date <- (ebirdst_runs[which(ebirdst_runs$common_name==species_target),7:8])
start_date <- as.numeric(start_date[1] + (start_date[2]-start_date[1])/2)
start_date <- yday(as.Date(start_date,origin = "1970-01-01"))
acceptable_week_dates <- seq(1,365,by=7)
start_date <- acceptable_week_dates[which(abs(acceptable_week_dates - start_date) == min(abs(acceptable_week_dates - start_date)))]
#print(start_date)

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
week_err_all <- c()
for (each_sim in 1:num_simulations){
  print(each_sim)
  
  #Parameter sets are drawn from a script in the data folder
  param_set <- readRDS("data/param_sets/WOTH_M1/9_20_22/WOTH_M1_S2.rds")
  #For best model
  #samp_ind <- which(param_set$sum_err == min(param_set$sum_err))
  samp_ind <- sample(1:nrow(param_set),1)
  for (each_var in 1:ncol(param_set)){
    assign(colnames(param_set)[each_var],param_set[sample(1:nrow(param_set),1),each_var])
  }
  
  #Initialize model
  initial_dfs <- initialize_MIBM(num_pulls,breeding_file,nonbreeding_file,speed_mean_s,
                                 speed_sd_s,speed_mean_f,speed_sd_f,start_date_u_s,
                                 start_date_sd_s,start_date_u_f,start_date_sd_f,
                                 max_mig_time_s,max_mig_time_f,bear_err_mean_s,
                                 bear_err_sd_s,bear_err_mean_f,bear_err_sd_f,
                                 max_energy_s,max_energy_f,recovery_rate_s,
                                 recovery_rate_f,season,start_date,goal_radius,
                                 mig_con,migr_timing_lat_s,migr_timing_lat_f)
  #Run model, get error
  error_run <- model_w_error(initial_dfs,start_date,spacing,weekly_err,weekly_modeled)
  
  #For each week, get the difference in MIMAS-estimate abundance and eBird ST-estimate abundance
  err_by_cell <- error_run[[3]]
  
  for (each_week in 1:52){
    week_err_in <- err_by_cell[[each_week]]
    week_err_in$sim_num <- each_sim
    week_err_in$week <- each_week
    week_err_in$model_diff <- week_err_in$count - week_err_in$perc_pop
    week_err_in <- week_err_in %>% select(cell,week,model_diff,sim_num)
    week_err_all <- rbind(week_err_all,week_err_in)
  }
}

#saveRDS(week_err_all,"week_err_all_WOTH.rds")
week_err_all <- readRDS("week_err_all_WOTH.rds")
avg_err_cell <- week_err_all %>% 
  group_by(cell) %>% 
  summarise(mean_diff = mean(model_diff))

#avg_err_cell <- avg_err_cell %>% filter(abs(mean_diff*100) > .02)

cell_size <- dgconstruct(spacing=150,metric=TRUE) 
grid <- dgcellstogrid(cell_size, avg_err_cell$cell, frame=TRUE)
grid$cell <- as.numeric(grid$cell)
grid <- left_join(avg_err_cell,grid,by="cell")

legend_position <- c(0.85, 0.6)
legend_title_size <- 15
legend_ax_size <- 15
legend_size <- 1

countries <- map_data("world")
countries <- countries %>% filter(long < -32)
p <-  ggplot() + coord_map("mollweide",xlim=c(-105,-65),ylim=c(5,50)) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="grey70", color="grey50") +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=mean_diff), alpha=.8)    +
  geom_path (data=grid,      aes(x=long, y=lat, group=group), alpha=0.2, color="white") + 
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="transparent", color="grey50") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(fill = "Error") +
  theme(panel.background = element_rect(fill = alpha("white",.5))) +
  scale_fill_gradient2(low="firebrick4",mid="white",high="dodgerblue4",
                       limits=c(min(avg_err_cell$mean_diff),max(avg_err_cell$mean_diff)),
                       labels = scales::label_percent()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        #panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.position = legend_position,
        legend.title = element_text(size=legend_title_size),
        legend.text = element_text(size=legend_ax_size),
        legend.key.size = unit(legend_size,"cm"),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"))
  #scale_fill_viridis(option="turbo",direction=-1,begin=.15,end=.95,
  #                   limits=c(min(avg_err_cell$mean_diff*100),max(avg_err_cell$mean_diff*100)))
p

