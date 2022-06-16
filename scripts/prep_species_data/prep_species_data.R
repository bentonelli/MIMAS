# IN DEVELOPMENT

# This script takes in a species name, and using ebird data creates:
# 1) Breeding, and non-breeding probability maps
# 2) Weekly probability maps
# 3) Weekly modeled area maps

# Saving each to folders in the data file

# For 2+3:
# This script reads abundance confidence intervals for a species to create weekly 
# expected maps in the dggridR format. These error maps can then be used for
# calculations for simulation error assessment.

# An important note here - specifying dggridR spacing here is necessary and important
# Here I am setting the spacing to 200 square km hexagons as a default. The
# spacing is relevant because it determines at what scale the error is calculated.
# A high spacing value will capture large-scale movements, but will not capture
# diffferences at a finer scale. 

library(ebirdst)
library(raster)
library(sf)
library(lubridate)
library(viridis)
library(dggridR)
library(dplyr)
library(readr)
sf::sf_use_s2(FALSE)

species_name <- "Yellow-bellied Sapsucker"
spacing <- 200
#Assign number of draws for each week of breeding, non-breeding season
num_draws <- 1000000
num_draws_b_nb <- 1000000

#Download data from ebird service
#set_ebirdst_access_key("fgjd8iu8kuar",overwrite = TRUE)
sp_path <- ebirdst_download(species = species_name,force=TRUE)
abunds <- load_raster("abundance", path = sp_path)
abunds_lower <- load_raster("abundance_lower", path = sp_path)
abunds_upper <- load_raster("abundance_upper", path = sp_path)

breed_dates <- week(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==species_name),8:9]))
nbreed_dates <- week(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==species_name),12:13]))

breed_weeks <- seq(from=breed_dates[1],to=breed_dates[2],by=1)

#Non-breeding seasons go over the new year normally, so make the adjustment to get the right weeks here
if (nbreed_dates[1]>nbreed_dates[2]){
  nbreed_weeks <- c(seq(from=nbreed_dates[1],to=52,by=1),seq(from=1,to=nbreed_dates[2],by=1))
} else {
  nbreed_weeks <- seq(from=nbreed_dates[1],to=nbreed_dates[2],by=1)
}

breed_points <- c()
nbreed_points <- c()

#Create list 
err_list <- list()
modeled_cell_list <- c()

#For each week
for (each_week in 1:52){
  print(each_week)
  #Calculate SD of ebird data from 10%-90% confidence intervals
  print("Calculating CI")
  conf_sd <- ((abunds_upper[[each_week]] - abunds_lower[[each_week]])/2)/1.282
  
  #Project, convert to df
  conf_proj <- projectRaster(conf_sd, crs = "+init=epsg:4326", method = "ngb")
  conf_points <- rasterToPoints(conf_proj)
  conf_points <- as.data.frame(conf_points)
  
  #Wipe memory
  conf_sd <- conf_proj <- c()
  
  #Name columns
  colnames(conf_points)[1:3] <- c("Lon","Lat","abundance_sd")
  
  #Do the same for the abundance estimates
  print("Trimming")
  abund_est <- trim(abunds[[each_week]])
  print("Done trimming")
  abunds_proj <- projectRaster(abund_est, crs = "+init=epsg:4326", method = "ngb")
  abunds_points <- rasterToPoints(abunds_proj)
  abunds_points <- as.data.frame(abunds_points)
  colnames(abunds_points)[1:3] <- c("Lon","Lat","abundance_umean")
  
  #Filter to non-zero cells
  abunds_points_trim <- filter(abunds_points,abundance_umean > 0)
  conf_points_trim <- filter(conf_points,abundance_sd > 0)
  
  #Construct grid based on spacing specified at the top
  dggs <- dgconstruct(spacing=spacing, metric=TRUE, resround='down')
  
  # Snap the eBird cells to their dggridR hexagon for both the confidence intervals,
  # as well as the SD estimates. Note here that this is not perfect,
  # as the eBird rectangles won't fit perfectly into those shapes!
  
  # For testing - get the modeled range in dggridR format. The goal here is just to get
  # a unique list of cells that have the abundance modeled. The lit can then be used
  modeled_cells <- dgGEO_to_SEQNUM(dggs,conf_points$Lon,conf_points$Lat)$seqnum
  #Get unique cells
  modeled_cells <- unique(modeled_cells)
  
  conf_points_trim$cell <- dgGEO_to_SEQNUM(dggs,conf_points_trim$Lon,conf_points_trim$Lat)$seqnum
  
  #Group by cell, get count mean of estimate, 
  conf_grp_by_cell <- conf_points_trim %>%
    group_by(cell) %>%
    summarise(cell_mean = mean(abundance_sd),n = n())
  
  abunds_points_trim$cell <- dgGEO_to_SEQNUM(dggs,abunds_points_trim$Lon,abunds_points_trim$Lat)$seqnum
  
  #Group by cell, get count mean of estimate, 
  abund_grp_by_cell <- abunds_points_trim %>%
    group_by(cell) %>%
    summarise(cell_mean = mean(abundance_umean),n = n())
  
  #Get percent of whole population estimated in each cell
  
  #The relative population is equal to number of cells multiplied by the mean abundance
  abund_grp_by_cell$rel_pop <- abund_grp_by_cell$cell_mean*abund_grp_by_cell$n
  
  #Then, take the value of each cell and calculate the percentage of the population in that cell
  abund_grp_by_cell$perc_pop <- abund_grp_by_cell$rel_pop/sum(abund_grp_by_cell$rel_pop)
  
  #Get error in the same terms (percent of population)
  
  # Relative error is equal to the 
  conf_grp_by_cell$rel_err <- conf_grp_by_cell$cell_mean*conf_grp_by_cell$n
  conf_grp_by_cell$err_adj <- conf_grp_by_cell$rel_err/sum(abund_grp_by_cell$rel_pop)
  
  #Merge percent population and adjusted err
  abund_w_error <- merge(abund_grp_by_cell,conf_grp_by_cell,by="cell")
  abund_w_error <- select(abund_w_error,cell,perc_pop,err_adj)
  
  #Instead of writing all weeks to a dataframe, save as two lists
  err_list[[each_week]] <- abund_w_error
  modeled_cell_list[[each_week]] <- modeled_cells
  
  
  ### Add abundance to breeding, non-breeding lists if week matches up ####
  if (each_week %in% c(nbreed_weeks,breed_weeks)){
    rel_abd_vlow <- abunds_points
    
    rel_abd_vlow$rel_abund <- rel_abd_vlow$abundance/(sum(rel_abd_vlow$abundance))
    
    # Draw points for sample. These points are not populated equally among the 2.96km squared area
    # as that minor adjustment is irrelevent for the error measurement.
    point_sample_index <- sample(1:nrow(rel_abd_vlow),num_draws,replace=TRUE,prob=rel_abd_vlow$rel_abund)
    point_sample <- rel_abd_vlow[point_sample_index,1:2]
    
    #Draw for the breeding, non-breeding maps, with
    
    #For weeks in breeding season, add to list
    if (each_week %in% breed_weeks){
      breed_points <- rbind(breed_points,rel_abd_vlow[point_sample_index,])
    }
    
    if (each_week %in% nbreed_weeks){
      nbreed_points <- rbind(nbreed_points,rel_abd_vlow[point_sample_index,])
    }
  }
}

#For weekly error, area modeled
saveRDS(err_list,paste("data/species_maps/error_maps/",gsub(" ","_",species_name),"_","s",spacing,".rds",sep=""))
saveRDS(modeled_cell_list,paste("data/species_maps/modeled_area_maps/",gsub(" ","_",species_name),"_","s",spacing,".rds",sep=""))

# Convert breeding points to relative abundance - note here that the breeding maps
# will lose some level of  
breed_group <- breed_points %>% 
  group_by(Lon,Lat) %>% 
  summarise(abundance = n())
breed_group <- as.data.frame(breed_group)
breed_group$perc_occupancy <- breed_group$abundance/sum(breed_group$abundance)
colnames(breed_group)[1:4] <- c("longitude","latitude","abundance","perc_occupancy")

#Do the same for non-breeding
nbreed_group <- nbreed_points %>% 
  group_by(Lon,Lat) %>% 
  summarise(abundance = n())
nbreed_group <- as.data.frame(nbreed_group)
nbreed_group$perc_occupancy <- nbreed_group$abundance/sum(breed_group$abundance)
colnames(nbreed_group)[1:4] <- c("longitude","latitude","abundance","perc_occupancy")

#Write both to csv format
write_csv(breed_group,paste("data/species_maps/breeding/",gsub(" ", "_", species_name),"_breeding_all",sep="",".csv"))
write_csv(nbreed_group,paste("data/species_maps/non_breeding/",gsub(" ", "_", species_name),"_non_breeding_all",sep="",".csv"))


