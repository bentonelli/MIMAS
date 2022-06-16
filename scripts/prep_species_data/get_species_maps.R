#IN DEVELOPMENT

#This script takes in data for a single bird, and breaks down the ebird data to
#a number of point draws, then writes a csv with those points. To be used in conjunction
# with model data to map accuracy.

library(ebirdst)
library(raster)
library(sf)
library(lubridate)
library(readr)
library(dplyr)

species_name <- "Wood Thrush"

#Assign number of draws for weekly maps
num_draws <- 100000

#Assign number of draws for each week of breeding, non-breeding season
num_draws_b_nb <- 1000000

#Download data
sp_path <- ebirdst_download(species = species_name)
abunds <- load_raster("abundance", path = sp_path)

breed_dates <- week(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==species_name),5:6]))
nbreed_dates <- week(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==species_name),7:8]))

breed_weeks <- seq(from=breed_dates[1],to=breed_dates[2],by=1)

#Non-breeding seasons go over the new year normally, so make the adjustment to get the right weeks here
if (nbreed_dates[1]>nbreed_dates[2]){
  nbreed_weeks <- c(seq(from=nbreed_dates[1],to=52,by=1),seq(from=1,to=nbreed_dates[2],by=1))
} else {
  nbreed_weeks <- seq(from=nbreed_dates[1],to=nbreed_dates[2],by=1)
}

breed_points <- c()
nbreed_points <- c()

#For each week, 
for (each_week in 1:52){
  print(each_week)
  #Split to week
  wk_abunds <- trim(abunds[[each_week]])
  
  #Convert to csvs
  abunds_vlow <- projectRaster(wk_abunds, crs = "+init=epsg:4326",method = "ngb")
  rel_abd_vlow <- as.data.frame(rasterToPoints(abunds_vlow))
  colnames(rel_abd_vlow)[1:3] <- c("Lon","Lat","abundance")
  rel_abd_vlow$rel_abund <- rel_abd_vlow$abundance/(sum(rel_abd_vlow$abundance))
  
  # Draw 100k points for sample. These points are not populated equally among the 2.96km squared area
  # as that minor adjustment is irrelevent for the error measurement.
  point_sample_index <- sample(1:nrow(rel_abd_vlow),num_draws,replace=TRUE,prob=rel_abd_vlow$rel_abund)
  point_sample <- rel_abd_vlow[point_sample_index,1:2]
  
  #Save to csv
  write_csv(point_sample,paste("data/species_maps/weekly/",gsub(" ", "_", species_name),
                               "/",gsub(" ", "_", species_name),"_Week_",each_week,sep="",".csv"))
  
  #Draw for the breeding, non-breeding maps, with
  
  #For weeks in breeding season, add to list
  if (each_week %in% breed_weeks){
    breed_points <- rbind(breed_points,rel_abd_vlow[point_sample_index,])
  }
  
  if (each_week %in% nbreed_weeks){
    nbreed_points <- rbind(nbreed_points,rel_abd_vlow[point_sample_index,])
  }
  
}

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


