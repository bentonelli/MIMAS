#Run a simulation of the IBM, then use those records to create a SIR model "on top"
# of the existing model
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

setwd("~/Documents/Coding/R/MIMAS/")

#Get functions
source(file = "scripts/sim_functions/sim_funct.R")

#Define species
species_target <- "Wood Thrush"

#Set number of birds to include in simulation
num_pulls <- 5000

#Define spacing of error maps
spacing <- 200

#Set number of days
num_days <- 365

# Align function
align_legend <- function(p, hjust = 0.5)
{
  # extract legend
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]
  
  # extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  
  # there can be multiple guides within one legend box  
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]
    
    # add extra column for spacing
    # guides$width[5] is the extra spacing from the end of the legend text
    # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
    # both sides, we get an aligned legend
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1-hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2
    
    # reconstruct guides and write back
    legend$grobs[[gi]] <- guides
  }
  
  # reconstruct legend and write back
  g$grobs[[legend_index]] <- legend
  g
}

#Import breeding, non-breeding maps
breeding_file <- read.csv(paste("data/species_maps/breeding/",gsub(" ","_",species_target),"_breeding_all.csv",sep=""))
nonbreeding_file <- read.csv(paste("data/species_maps/non_breeding/",gsub(" ","_",species_target),"_non_breeding_all.csv",sep=""))

#Read in weekly error maps - not explicitly needed here in the simple run, but good to know how to do
weekly_err <- readRDS(paste("data/species_maps/error_maps/",gsub(" ","_",species_target),"_s",spacing,".rds",sep=""))

#Read in weekly modeled areas - also not explicitly needed
weekly_modeled <- readRDS(paste("data/species_maps/modeled_area_maps//",gsub(" ","_",species_target),"_s",spacing,".rds",sep=""))

#Define starting season
season <- 1

start_date <- (ebirdst_runs[which(ebirdst_runs$common_name==species_target),7:8])
start_date <- as.numeric(start_date[1] + (start_date[2]-start_date[1])/2)
start_date <- yday(as.Date(start_date,origin = "1970-01-01"))
acceptable_week_dates <- seq(1,365,by=7)
start_date <- acceptable_week_dates[which(abs(acceptable_week_dates - start_date) == min(abs(acceptable_week_dates - start_date)))]
print(start_date)

tick_record_list <- list()
for (nn in 1:50){
  print(paste("working on simulation number:",nn))
  param_set <- readRDS("data/param_sets/WOTH_M1/9_20_22/WOTH_M1_S2.rds")
  sample_int <- sample(1:nrow(param_set),1)
  for (each_var in 1:ncol(param_set)){
    #Random
    assign(colnames(param_set)[each_var],param_set[sample_int,each_var])
    #Get best
    #assign(colnames(param_set)[each_var],param_set[which(param_set$sum_err == min(param_set$sum_err)),each_var])
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
  static_df <- initial_dfs[[1]]
  upd_df <- initial_dfs[[2]]
  
  locs_rec <- as.data.frame(matrix(NA,ncol = 5,nrow = (num_days*num_pulls)))
  colnames(locs_rec) <- c("ID","Lon","Lat","Day","season")
  for (days in 1:num_days){
    #print(days)
    upd_df <- run_day(upd_df,static_df)
    summary(upd_df)
    locs <- save_locations(upd_df,save_season = TRUE)
    loc_ind_s <- (days - 1)*num_pulls + 1 
    loc_ind_e <- (days - 1)*num_pulls + num_pulls 
    locs_rec[loc_ind_s:loc_ind_e,] <- locs
  }
  
  #Find birds in Yucatan, only during springtime - 50: 200
  yuc_records <- filter(locs_rec, Lon > -92 & Lon < -86 & Lat < 23 & Lat > 19,Day < 150 & Day > 50)
  nrow(yuc_records)
  
  #Save box for later
  boxes<-data.frame(maxlat = 23,minlat = 19,maxlong = -86,minlong = -92, id="1")
  boxes<-transform(boxes, laby=(maxlat +minlat )/2, labx=(maxlong+minlong )/2)
  
  #Create a attachment rate based on date, here increasing linearly over the course
  # of the spring
  
  yuc_records$tick_attach <- rbinom(nrow(yuc_records),size = 1,prob = (yuc_records$Day*.0005))
  table(yuc_records$Day,yuc_records$tick_attach)
  sum(yuc_records$tick_attach)
  
  #Simulate detach date
  yuc_records$tick_detach_date <- (yuc_records$Day + rpois(nrow(yuc_records),2)) * yuc_records$tick_attach
  
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
    
    tick_locations <- select(tick_locations,c("Lon","Lat","Day","season","tick_num"))
    #get rid of tick records that don't leave original location
    tick_locations <- tick_locations[which(diff_loc!=0),]
    tick_locations <- filter(tick_locations,Lat > 23 | Lat < 19 | Lon > -86 | Lon < -92)
    
    #Get
    tick_record <- rbind(tick_record,tick_locations)
  }
  tick_record_list[[nn]] <- tick_record
}

#saveRDS(tick_record_list,"tick_record_list_WOTH.rds")

### Dispersal Plots #### 

# Looking at two different species here, HOWA and WOTH the goal is to make
# comparisons between the two
tick_record_list <- readRDS("tick_record_list_WOTH.rds")
total_sim_birds <- 5000
pop_size <- 12000000
total_ticks <- c()
tick_record_all <- c()
for (mm in 1:50){
  tick_record <- tick_record_list[[mm]]
  tick_record_all <- rbind(tick_record_all,tick_record)
  total_ticks <- c(total_ticks,nrow(tick_record))
}
total_ticks_pop <- (total_ticks * (pop_size/total_sim_birds))
plot(density(total_ticks_pop/1000),main="",lwd=8,col="tan4",ylim=c(0,.0055),
     xlab="",ylab="",yaxt='n',cex.axis=1.25, cex.lab=1.75)


tick_record_list2 <- readRDS("tick_record_list_HOWA.rds")
total_sim_birds <- 5000
pop_size <- 5400000
total_ticks2 <- c()
tick_record_all2 <- c()
for (mm in 1:50){
  tick_record2 <- tick_record_list2[[mm]]
  tick_record_all2 <- rbind(tick_record_all2,tick_record2)
  total_ticks2 <- c(total_ticks2,nrow(tick_record2))
}
total_ticks_pop2 <- (total_ticks * (pop_size/total_sim_birds))
points(density(total_ticks_pop2/1000),type="l",main="",lwd=8,col="goldenrod2")
legend("topright",legend=c("Wood Thrush","Hooded Warbler"),col=c("tan4","goldenrod2"),lwd=8,cex=2)

### Species 1 map ####
tick_record_list <- readRDS("tick_record_list_WOTH.rds")
total_sim_birds <- 250000
pop_size <- 12000000
total_ticks <- c()
tick_record_all <- c()
for (mm in 1:50){
  tick_record <- tick_record_list[[mm]]
  tick_record_all <- rbind(tick_record_all,tick_record)
  total_ticks <- c(total_ticks,nrow(tick_record))
}


cell_size <- dgconstruct(spacing=150,metric=TRUE) 
tick_record_all$cell <- dgGEO_to_SEQNUM(cell_size, tick_record_all$Lon, tick_record_all$Lat)$seqnum
tick_record_grp <- tick_record_all %>%
  group_by(cell) %>%
  summarise(total_tick_count = sum(tick_num))

tick_record_grp$total_tick_count <- tick_record_grp$total_tick_count * (pop_size/total_sim_birds)

tick_record_grp <- filter(tick_record_grp,total_tick_count > 100)

grid <- dgcellstogrid(cell_size, tick_record_grp$cell, frame=TRUE)
grid$cell <- as.numeric(grid$cell)
grid <- left_join(tick_record_grp,grid,by="cell")
grid$adj_tick <- grid$total_tick_count


legend_position <- c(.875,.6)
legend_title_size <- 22
legend_ax_size <- 18
legend_size <- 1.75
country_fill_color <- "grey70"
country_border_color <- "grey50"
par(mar=c(0,0,0,0))

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

boxes<-data.frame(maxlat = 23,minlat = 19,maxlong = -86,minlong = -92, id="1")
boxes<-transform(boxes, laby=(maxlat +minlat )/2, labx=(maxlong+minlong )/2)



p <-  ggplot() + coord_map("mollweide",xlim=c(-105,-65),ylim=c(5,50)) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=country_fill_color, color=country_border_color) +
  #geom_polygon(data=states, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=adj_tick), alpha=.8)    +
  scale_fill_gradient(low="grey", high="firebrick3",limits=c(0,12000)) + 
  geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.7, color="white") +
  labs(fill = "# Hitchhikers") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.position = legend_position,
        legend.title = element_text(size=legend_title_size),
        legend.text = element_text(size=legend_ax_size),
        legend.key.size = unit(legend_size,"cm"),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  theme(panel.background = element_rect(fill = alpha("white",.5)))
p <- ggdraw(align_legend(p))
p

### Species 2 map ####
tick_record_list <- readRDS("tick_record_list_HOWA.rds")
total_sim_birds <- 250000
pop_size <- 5400000
total_ticks <- c()
tick_record_all <- c()
for (mm in 1:50){
  tick_record <- tick_record_list[[mm]]
  tick_record_all <- rbind(tick_record_all,tick_record)
  total_ticks <- c(total_ticks,nrow(tick_record))
}


cell_size <- dgconstruct(spacing=150,metric=TRUE) 
tick_record_all$cell <- dgGEO_to_SEQNUM(cell_size, tick_record_all$Lon, tick_record_all$Lat)$seqnum
tick_record_grp <- tick_record_all %>%
  group_by(cell) %>%
  summarise(total_tick_count = sum(tick_num))

tick_record_grp$total_tick_count <- tick_record_grp$total_tick_count * (pop_size/total_sim_birds)

tick_record_grp <- filter(tick_record_grp,total_tick_count > 100)

grid <- dgcellstogrid(cell_size, tick_record_grp$cell, frame=TRUE)
grid$cell <- as.numeric(grid$cell)
grid <- left_join(tick_record_grp,grid,by="cell")
grid$adj_tick <- grid$total_tick_count

legend_position <- c(.875,.6)
legend_title_size <- 22
legend_ax_size <- 18
legend_size <- 1.75
country_fill_color <- "grey70"
country_border_color <- "grey50"
par(mar=c(0,0,0,0))

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

boxes<-data.frame(maxlat = 23,minlat = 19,maxlong = -86,minlong = -92, id="1")
boxes<-transform(boxes, laby=(maxlat +minlat )/2, labx=(maxlong+minlong )/2)



p2 <-  ggplot() + coord_map("mollweide",xlim=c(-105,-65),ylim=c(5,50)) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=country_fill_color, color=country_border_color) +
  #geom_polygon(data=states, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=adj_tick), alpha=.8)    +
  scale_fill_gradient(low="grey", high="firebrick3",limits=c(0,12000)) + 
  geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.7, color="white") +
  labs(fill = "# Hitchhikers") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.position = legend_position,
        legend.title = element_text(size=legend_title_size),
        legend.text = element_text(size=legend_ax_size),
        legend.key.size = unit(legend_size,"cm"),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  theme(panel.background = element_rect(fill = alpha("white",.5)))
p2 <- ggdraw(align_legend(p2))
p2

comb_plots <- ggpubr::ggarrange(p, p2, 
                                ncol = 2, nrow = 1)
comb_plots

