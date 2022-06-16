#IN DEVELOPMENT

# This script reads abundance confidence intervals for a species to create weekly 
# expected maps in the dggridR format. These error maps can then be used for
# calculations for simulation error assessment.

# An important note here - specifying dggridR spacing here is necessary and important
# Here I am setting the spacing to 200 square km hexagons as a default. The
# spacing is relevant because it determines at what scale the error is calculated.
# A high spacing value will capture large-scale movements, but will not capture
# high density flyways. 

# The end goal here is a weekly error dggridR-usable dataframe that describes expected 
# distributions of birds across their range.

library(ebirdst)
library(raster)
library(sf)
library(lubridate)
library(viridis)
library(dggridR)
library(dplyr)

species_name <- "Wood Thrush"
spacing <- 200

#Download data from ebird service
sp_path <- ebirdst_download(species = species_name)
abunds <- load_raster("abundance", path = sp_path)
abunds_lower <- load_raster("abundance_lower", path = sp_path)
abunds_upper <- load_raster("abundance_upper", path = sp_path)

#Create list 
err_list <- list()
modeled_cell_list <- c()

#For each week
for (each_week in 1:52){
  print(each_week)
  #Calculate SD of ebird data from 10%-90% confidence intervals
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
  abund_est <- trim(abunds[[each_week]])
  abunds_proj <- projectRaster(abund_est, crs = "+init=epsg:4326", method = "ngb")
  abunds_points <- rasterToPoints(abunds_proj)
  abunds_points <- as.data.frame(abunds_points)
  colnames(abunds_points)[1:3] <- c("Lon","Lat","abundance_umean")
  
  #Wipe memory
  abund_est <- abunds_proj <- c()
  
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
  
  # Relative error 
  conf_grp_by_cell$rel_err <- conf_grp_by_cell$cell_mean*conf_grp_by_cell$n
  conf_grp_by_cell$err_adj <- conf_grp_by_cell$rel_err/sum(abund_grp_by_cell$rel_pop)
  
  #Merge percent population and adjusted err
  abund_w_error <- merge(abund_grp_by_cell,conf_grp_by_cell,by="cell")
  abund_w_error <- select(abund_w_error,cell,perc_pop,err_adj)
  
  #Instead of writing all weeks to a dataframe, save as two lists
  err_list[[each_week]] <- abund_w_error
  modeled_cell_list[[each_week]] <- modeled_cells
}

saveRDS(err_list,paste("data/species_maps/error_maps/",gsub(" ","_",species_name),"_","s",spacing,".rds",sep=""))
saveRDS(modeled_cell_list,paste("data/species_maps/modeled_area_maps/",gsub(" ","_",species_name),"_","s",spacing,".rds",sep=""))

# week_in <- Varied_Thrush_s200[[11]]
# #Plotting for testing
# grid  <- dgcellstogrid(dggs,week_in$cell,frame=TRUE,wrapcells=TRUE)
# grid  <- merge(grid,week_in,by.y="cell")
# countries <- map_data("world")
# p<- ggplot() +
#   geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
#   geom_polygon(data=grid,      aes(x=long, y=lat, group=group,fill=perc_pop), alpha=0.4)    +
#   geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
#   scale_fill_gradient(low="blue", high="red") + xlim(-170,-100) + ylim(20,75)
# 
# p
