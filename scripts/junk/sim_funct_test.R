### This script contains all the relevant functions for the first jab at the 
### MIBM program. All the scripts here will need to be read into memory via 
### the script that defines parameters.

# Initialization function - create initial dataframes for simulations
initialize_MIBM <- function(num_pulls,breeding_file,nonbreeding_file,speed_mean_s,
                            speed_sd_s,speed_mean_f,speed_sd_f,start_date_u_s,
                            start_date_sd_s,start_date_u_f,start_date_sd_f,
                            max_mig_time_s,max_mig_time_f,bear_err_mean_s,
                            bear_err_sd_s,bear_err_mean_f,bear_err_sd_f,
                            max_energy_s,max_energy_f,recovery_rate_s,
                            recovery_rate_f,season,start_date,goal_radius,
                            migr_con=0,mig_con_type=0,migr_timing_lat_s=0,
                            migr_timing_lat_f=0){
  
  ### Create Static DF - This dataframe contains information for the simulations
  # that does not change. 
  
  #For each of these sets, draw number of random points, adjust for the kernel
  #size, then add together. The breeding points will represent start points while
  #the non-breeding points represent the goal points. This is done randomly here,
  #but there is no doubt that some species with migratory connectivity will need
  #a non-random assorting here.
  
  #Record for points
  matched_points <- c()
  
  #Read in the abundance lists from the ebird data and construct breeding/non-
  #breeding points
  
  #Create counter
  for (abundance_lists in list(breeding_file,nonbreeding_file)){
    
    #Get random locations - as vector
    rand_loc <- sample(nrow(abundance_lists),num_pulls,prob=abundance_lists$perc_occupancy,replace=TRUE)
    #Save lat-lng
    rand_lat_long <- abundance_lists[rand_loc,1:2]
    rand_lat_long$long_trans <- rand_lat_long$longitude
    rand_lat_long$lat_trans <- rand_lat_long$latitude
    
    rand_lat_long$long_trans <- ebird_lat_long_tran(rand_lat_long$longitude,rand_lat_long$lat_trans)[,1]
    rand_lat_long$lat_trans <- ebird_lat_long_tran(rand_lat_long$longitude,rand_lat_long$lat_trans)[,2]
    #For each of these points, transform to fit grid size of ebird (see function)
    # for (n in 1:nrow(rand_lat_long)){
    #   rand_lat_long$long_trans[n] <- ebird_lat_long_tran(rand_lat_long$longitude[n],rand_lat_long$lat_trans[n])[1]
    #   rand_lat_long$lat_trans[n] <- ebird_lat_long_tran(rand_lat_long$longitude[n],rand_lat_long$lat_trans[n])[2]
    # }
    #save points as matrices
    save_points <- as.matrix(rand_lat_long[,3:4])
    
    matched_points <-cbind(matched_points,save_points)
    
  }
  #Save as dataframe
  matched_points <- as.data.frame(matched_points)
  
  # A note on the above -> this does not take into account that some birds may get 
  # locations that are "over water" according to the over_water function due to 
  # minor mismatches in the eBird package. For now I ignore this minor issue, in
  # order to maintain birds with locations on small islands and coastlines.
  
  #Add bird IDs. These IDs will follow birds and allow for additional functions
  # to be inteegrated in later on.
  matched_points <- data.frame(ID = 1:num_pulls, matched_points)
  
  #Update column names
  colnames(matched_points)[2] <- "breeding_lon"
  colnames(matched_points)[3] <- "breeding_lat"
  colnames(matched_points)[4] <- "nonbreeding_lon"
  colnames(matched_points)[5] <- "nonbreeding_lat"
  
  # This is where migratory connectivity comes in. Here, migratory connectivity
  # is modeled to be a function of latitude, that is, birds breeding father north
  # are more likely to breed further south (leapfrog migration, type 1). Or alternatively,
  # those breeding farther north are more likely to winter farther north (type 2).
  # However, this  section could be adapted in a number of ways to take into account different
  # types of connectivity. 
  
  if (mig_con_type == 1){
    nb_ranks <- 1:num_pulls
    probs <- (nb_ranks/num_pulls)*mig_con
    
    #Add noise
    new_ind <- rnorm(num_pulls,0,1-mig_con)
    new_ind <- rank(probs+new_ind)
    
    #Rank whole df breeding latitudes
    matched_points <- matched_points[order(matched_points$breeding_lat,decreasing=TRUE),]
    
    #Then order non-breeding latitudes,longitudes
    nb_loc <- matched_points[,4:5]
    nb_loc <- nb_loc[order(nb_loc$nonbreeding_lat,decreasing=FALSE),]
    #Now match based on new_ind
    nb_loc <- cbind(nb_loc$nonbreeding_lon[new_ind],nb_loc$nonbreeding_lat[new_ind])
    matched_points[,4:5] <- nb_loc
    
    #Revert sorting to ID
    matched_points <- matched_points[order(matched_points$ID,decreasing=FALSE),]
  } else if (mig_con_type == 2) {
    
    nb_ranks <- 1:num_pulls
    probs <- (nb_ranks/num_pulls)*mig_con
    
    #Add noise
    new_ind <- rnorm(num_pulls,0,1-mig_con)
    new_ind <- rank(probs+new_ind)
    
    #Rank whole df breeding latitudes
    matched_points <- matched_points[order(matched_points$breeding_lat,decreasing=TRUE),]
    
    #Then order non-breeding latitudes,longitudes
    nb_loc <- matched_points[,4:5]
    nb_loc <- nb_loc[order(nb_loc$nonbreeding_lat,decreasing=TRUE),]
    #Now match based on new_ind
    nb_loc <- cbind(nb_loc$nonbreeding_lon[new_ind],nb_loc$nonbreeding_lat[new_ind])
    matched_points[,4:5] <- nb_loc
    
    #Revert sorting to ID
    matched_points <- matched_points[order(matched_points$ID,decreasing=FALSE),]
  }
  
  
  # Assign parameter values for all. Note that these are unchanging and do not
  # vary by individual, unless otherwise noted (see departure timing)
  
  #Migratory flight length variables
  matched_points$speed_mean_s <- speed_mean_s
  matched_points$speed_sd_s <- speed_sd_s
  matched_points$speed_mean_f <- speed_mean_f
  matched_points$speed_sd_f <- speed_sd_f
  
  #Orientation variables
  matched_points$bear_err_mean_s <- bear_err_mean_s
  matched_points$bear_err_sd_s <- bear_err_sd_s
  matched_points$bear_err_mean_f <- bear_err_mean_f
  matched_points$bear_err_sd_f <- bear_err_sd_f
  
  #Energetic Variables
  matched_points$max_energy_s <- max_energy_s
  matched_points$max_energy_f <- max_energy_f
  matched_points$recovery_rate_s <- recovery_rate_s
  matched_points$recovery_rate_f <- recovery_rate_f
  
  #Phenology variables - These vary by bird
  matched_points$start_date_s <- round(rnorm(num_pulls,start_date_u_s,start_date_sd_s))
  matched_points$start_date_f <- round(rnorm(num_pulls,start_date_u_f,start_date_sd_f))
  
  # Here, if migration start timing varies by latitude, with birds at higher breeding 
  # latitudes initiating flight earlier, we adjust reorder these start dates
  
  if (migr_timing_lat_f > 0 | migr_timing_lat_s > 0){
    ranks <- 1:num_pulls
    
    probs_s <- (ranks/num_pulls)*migr_timing_lat_s
    probs_f <- (ranks/num_pulls)*migr_timing_lat_f
    
    #Draw ranks for spring, fall
    new_ind_s <- rnorm(num_pulls,0,1-migr_timing_lat_s)
    new_ind_s <- rank(probs_s+new_ind_s)
    
    new_ind_f <- rnorm(num_pulls,0,1-migr_timing_lat_f)
    new_ind_f <- rank(probs_f+new_ind_f)
    
    #Rank whole df breeding latitudes
    matched_points <- matched_points[order(matched_points$breeding_lat,decreasing=TRUE),]
    
    #Then order non-breeding latitudes,longitudes
    sd_s <- sort(matched_points$start_date_s)[new_ind_s]
    sd_f <- sort(matched_points$start_date_f)[new_ind_f]
    
    matched_points$start_date_s <- sd_s
    matched_points$start_date_f <- sd_f
    
    #Revert sorting to ID
    matched_points <- matched_points[order(matched_points$ID,decreasing=FALSE),]
  }
  
  # Max migration timing (this is a stopgap for birds that get really lost, with
  # the assumption that those individuals eventually "give up" and stay put).
  # This variable is important from a modeleing point of view because it stops birds from
  # being stuck in a migratory compartment.
  matched_points$max_mig_time_s <- max_mig_time_s
  matched_points$max_mig_time_f <- max_mig_time_f
  
  #Goal radius defines the area where a bird "homes," no longer obeying the general
  #migratory program, and instead finding its breeeding/wintering location.
  matched_points$goal_radius <- goal_radius
  
  ### Create Updatable dataframe ####
  #This dataframe changes as simulations progress, and stores the current state
  #of the system. Of the most interest is the current location of indivudals in the system.
  
  #Build dataframe
  upd_df <- as.data.frame(1:num_pulls)
  
  #Assign IDs
  colnames(upd_df)[1] <- "ID"
  
  #Assign current location as starting locations
  upd_df$current_lon <- matched_points$breeding_lon
  upd_df$current_lat <- matched_points$breeding_lat
  
  #Assign start date of simulation
  upd_df$day <- start_date
  
  #Initialize season (1 = breeding, 2 = spring migration, 3= non-breeding, 4 = fall migration)
  upd_df$season <- season
  
  #Assign birds to have max energy
  upd_df$energy <- matched_points$max_energy_f
  
  #No birds in stopover
  upd_df$stopover_counter <- 0
  
  #Start migration counter
  upd_df$migr_counter <- 0
  
  #Convert both to dataframes for speed considerations
  matched_points <- as.data.frame(matched_points)
  upd_df <- as.data.frame(upd_df)
  
  return(list(matched_points,upd_df))
}

#Function to take a centroid point of Ebird data and spread to the entire sample area. 
#In this case the kernel size is 2.96km in both directions for eBird data
ebird_lat_long_tran <- function(longs,lats){
  #Radius of the earth
  r_earth <- 6378
  #Kernal size is 2.96 so points need to be blown up by 1.48 in either direction
  dx <- runif(length(longs), min=-1.48, max=1.48)
  dy <- runif(length(lats), min=-1.48, max=1.48)
  new_latitudes  <- lats  + (dy / r_earth) * (180 / pi);
  new_longitudes <- longs + (dx / r_earth) * (180 / pi) / cos(lats * pi/180);
  return(cbind(new_longitudes,new_latitudes))
}

#Function to pass back updated dataframe with updated seasons. For birds that are
# transitioning between seasonal compartments
start_migration <- function(updatable_df,static_df){
  
  #Get current day
  day_num <- updatable_df$day[1]
  
  #For current breeding birds - start date is equal to day, the season still shows breeding
  update_list <- which(static_df$start_date_f == day_num & updatable_df$season == 1)
  
  #Update season - fall migration
  updatable_df[update_list,"season"] <- 2
  #Change energy to max
  updatable_df[update_list,"energy"] <- static_df[update_list,"max_energy_f"] 
  #Reset stopover counter
  updatable_df[update_list,"stopover_counter"] <- 0 
  #Reset migration counter
  updatable_df[update_list,"migr_counter"] <- 0
  
  #For current non-breeding birds - start date is equal or less than day, the season still shows breeding
  update_list <- which(static_df$start_date_s == day_num & updatable_df$season == 3)
  #Update season - spring migration
  updatable_df[update_list,"season"] <- 4
  #Change energy to max
  updatable_df[update_list,"energy"] <- static_df[update_list,"max_energy_s"] 
  #Reset stopover counter
  updatable_df[update_list,"stopover_counter"] <- 0 
  #Reset migration counter
  updatable_df[update_list,"migr_counter"] <- 0
  
  return(updatable_df)
}

#Function to determine if bird has reached goal
at_goal <- function(current_locations,goal_locations,goal_radii){
  if (nrow(current_locations) > 0){
    
    #Get the current distance to goal in km
    distance_to_goal <- distRhumb(current_locations,goal_locations)/1000
    
    # If the bird is within the goal radius, it "tracks" this day to its goal point
    at_goal_b <- distance_to_goal < goal_radii
    return(at_goal_b)
  }
  #Currently there is no energy deduction here. It may be worth considering one
  #depending on the context.
  return(FALSE)
}

#Function that takes in matrix of migrants and determines which should initiate stopover
# Forced stopovers are ethe result of the inability of an indiviudal to find a flight
# that doesn't end overwater.
initiate_stopover <- function(current_migrants,force=FALSE){
  
  #Energy equation -- percent chance of stopover = 1 - (energy/max_energy)
  #This equation is author-estimated, could be changed depending on the context
  # Example, stopover is determined in part by habitat,weather etc.
  stopover_perc <- 1 - (current_migrants$energy/current_migrants$max_energy)
  
  #If stopover is forced, no need to pull probabilities
  if (force == FALSE){
    #Draw random number
    rn <- runif(nrow(current_migrants),0,1)
  } else {
    #Always greater than 1
    rn <- numeric(nrow(current_migrants))
  }
  
  #Get birds that stopover
  stp_indices <- which(rn < stopover_perc)
  
  #Bird initiates stopover, add counter to calculate the number of days until
  #full energy is reached (this is rounded down to avoid 'over-fattening'. 
  #Recovery rate should be a proportion
  
  #Again, this is a place where changes could be made to account for varying
  # recovery rates due to habitat, weather, etc.
  current_migrants[stp_indices,"stopover_counter"] <- 
    floor(((current_migrants[stp_indices,"max_energy"] - 
              current_migrants[stp_indices,"energy"]) / 
             current_migrants[stp_indices,"max_energy"]) / 
            current_migrants[stp_indices,"recovery_rate"])
  
  return(current_migrants)
  
}

get_speeds <- function(speeds_mean,speeds_sd,energy,distances_to_goal){
  
  #Calculate today's speed - correct for values below zero
  todays_speed <- rtruncnorm(length(speeds_mean),a=0,mean = speeds_mean,sd=speeds_sd)
  
  #If the speed is more than the birds energy, limit to max energy of the bird
  todays_speed[which(todays_speed > energy)] <- energy[which(todays_speed > energy)]
  
  #If the speed is more than the distance the bird needs to fly, correct to that distance.
  # This stops birds from continually overshooting their target
  todays_speed[which(todays_speed > distances_to_goal)] <- distances_to_goal[which(todays_speed > energy)]
  
  return(todays_speed)
}

### Function for overwater check. This operates with the assumption that birds
# should not be able to to land overwater (something to change if the model is 
#setup for waterbirds)
over_water <- function(locations){
  ow_time_s <- Sys.time()
  #Get locations
  locations <- as.data.frame(locations)
  
  # A failure of the Rhumb line system may throw an error normally, so here 
  # we change any NaNs (which indicate this failure) to the coordinates (-150,90 - the North Pole)
  # which is over water
  if (sum(is.na(locations))>0){
    locations[which(is.na(locations[,1])),] <- c(-150,90)
  }
  #Convert points to sf format
  pts <- st_as_sf(locations, coords=1:2, crs=4326)
  
  #This will return the values that are over water/land (TRUE,FALSE) and should be repeated
  #this function  has given me issues previously, and is a good place to start debugging :)
  water <- is.na(as.numeric(suppressMessages(st_intersects(pts, world))))
  ow_time_e <- Sys.time()
  ow_time_add <- (as.numeric(ow_time_e) - as.numeric(ow_time_s))
  ow_time <<- ow_time + ow_time_add
  return(water)
}

#Function to update season due to completed or failed migration (over limit)
update_season <- function(current_migrants) {
  
  #Only run if length > 0
  if (nrow(current_migrants) > 0){
    
    #Update energy to max (should  be redundant)
    current_migrants[,"energy"] <- current_migrants[,"max_energy"]
    
    #Reset migration counter,stopover counter (again, probably redundant)
    current_migrants[,"stopover_counter"] <- 0
    current_migrants[,"migr_counter"] <- 0
    
    #Change season to breeding, non/breeding
    current_migrants[,"season"] <- current_migrants[,"season"]%%4+1
  }
  return(current_migrants)
}

# This is one of the larger functions designed to move birds in migration
move_bird <- function(current_migrants_to_move){
  #Skip if there are no migrants
  if (nrow(current_migrants_to_move) > 0){
    
    #Get current locations
    current_locations <- cbind(current_migrants_to_move$current_lon,current_migrants_to_move$current_lat)
    
    #Get goal locations
    goal_locations <- cbind(current_migrants_to_move$goal_lon,current_migrants_to_move$goal_lat)
    
    #Pull ideal bearing
    ideal_bearings <- bearingRhumb(current_locations,goal_locations)
    
    #Pull distance to goal
    distances_to_goal <- distRhumb(current_locations,goal_locations)
    
    #Get speeds
    speeds <- current_migrants_to_move$speed_mean
    speeds_sd <- current_migrants_to_move$speed_sd
    energy <- current_migrants_to_move$energy
    
    #Get speeds in km, convert to meters
    todays_speeds <- get_speeds(speeds,speeds_sd,energy,distances_to_goal)
    
    #Pull bearing w/ error
    todays_bearing_errors <- rnorm(nrow(current_migrants_to_move),
                                   mean = current_migrants_to_move$bear_err_u, 
                                   sd = current_migrants_to_move$bear_err_sd)
    
    # Now, each bird will attempt to fly. Some birds will end their flights over water. 
    # We assume no water landings for our landbirds, instead we check whether the flight ended
    # over water and redo the calculations for that flight. There are various possibilities
    # for how to parameterize this, but I have chosen to start with the assumption that a bird
    # would not initiate a flight that would end over water, and
    new_locations <- suppressWarnings(destPointRhumb(current_locations,(ideal_bearings+todays_bearing_errors),todays_speeds*1000))
    
    #Check which are overwater
    on_water <- over_water(new_locations)
    
    #Create df for the current locations, new locations, and the energy cost. Note that each bird gets a new record, and 
    # because some birds will end over water (See below) their "new location" will need to be their current location
    flight_record <- as.data.frame(cbind(current_migrants_to_move$ID,
                                         current_migrants_to_move$current_lon,
                                         current_migrants_to_move$current_lat,
                                         new_locations[,1],new_locations[,2],todays_speeds,on_water))
    
    colnames(flight_record)[1:7] <- c("ID","current_lon","current_lat","new_lon","new_lat","dist","water")
    
    #Give flight record a stopover index 
    flight_record$stopover_counter <- 0
    
    # For those that are overwater, run 10 attempts at alternative flights. Because significant water
    # barriers may mean that birds need to either fly much further than they would normally, or "stage"
    # before the flights, their behavior in these cases changes slightly. The normal distribution for the speed 
    # and bearing error of the flight now has an adjusted mean, increasing each iteration. 
    # In addition, birds that don't find a suitable landing (over land) after 10 iterations will 
    # be forced to take a stopover. If the birds cannot make a flight without landing
    # on land, it probably suggests that they don't have enough fuel in the tank, and need to stopover.
    
    #Get birds that are over water
    over_water_migrants <- current_migrants_to_move[on_water,]
    
    #Only run if some birds are over water
    if (nrow(over_water_migrants) > 0){
      #I am now setting up this section to do all the calculations at once, rather
      # than having it run in a loop. My hope is that this will greatly shorten
      # the runtime for this section
      
      #Get the multipliers for adjusted factors below
      water_mult <- rep(1.1^(1:10),nrow(over_water_migrants))
      
      #Replicate df 10 times for 10 tries
      over_water_migrants_rep <- do.call("rbind", replicate(10, over_water_migrants, simplify = FALSE))

      #Add sequence to dataframe for use later
      over_water_migrants_rep$test_num <- 1:nrow(over_water_migrants_rep)
      
      #Add water_mult to to df
      over_water_migrants_rep$water_mult <- water_mult
      
      #Pull ideal bearing
      ideal_bearings <- bearingRhumb(cbind(over_water_migrants_rep$current_lon,over_water_migrants_rep$current_lat),
                                     cbind(over_water_migrants_rep$goal_lon,over_water_migrants_rep$goal_lat))
      
      #Pull distance to goal
      distances_to_goal <- distRhumb(cbind(over_water_migrants_rep$current_lon,over_water_migrants_rep$current_lat),
                                     cbind(over_water_migrants_rep$goal_lon,over_water_migrants_rep$goal_lat))
      
      #Get speeds
      speeds <- over_water_migrants_rep$speed_mean * over_water_migrants_rep$water_mult
      speeds_sd <- over_water_migrants_rep$speed_sd * over_water_migrants_rep$water_mult
      energy <- over_water_migrants_rep$energy
      todays_speeds_ow <- get_speeds(speeds,speeds_sd,energy,distances_to_goal)
      
      #Pull bearing w/error
      todays_bearing_errors <- rnorm(nrow(over_water_migrants_rep),
                                     mean = over_water_migrants_rep$bear_err_u, 
                                     sd = (over_water_migrants_rep$bear_err_sd*over_water_migrants_rep$water_mult))
      
      #Get new locations
      new_locations_ow <- suppressWarnings(destPointRhumb(cbind(over_water_migrants_rep$current_lon,over_water_migrants_rep$current_lat),
                                         (ideal_bearings+todays_bearing_errors),todays_speeds_ow*1000))
      #Check which are overwater, save locations, speed
      over_water_migrants_rep$on_water_ow <- over_water(new_locations_ow)
      over_water_migrants_rep$new_lon <- new_locations_ow[,1]
      over_water_migrants_rep$new_lat <- new_locations_ow[,2]
      over_water_migrants_rep$dist <- todays_speeds_ow
      
      # Now the goal is to get the first draw that has a bird not overwater and 
      # use that as the new location of the bird
      
      #If at least one bird has one draw over land
      if (sum(over_water_migrants_rep$on_water_ow == FALSE) > 0){
        trim_list <- over_water_migrants_rep %>%
          filter(on_water_ow == FALSE) %>%
          group_by(ID) %>%
          summarise(min_ow = min(test_num))
        # Get birds that found an adequate location over land, and use that location
        # to update their current state
        over_water_migrants_ol <- over_water_migrants_rep[trim_list$min_ow,]
        #Now take the update flight_record frame, and use it to update the current_migrants_to_move df
        replace_index <- match(as.numeric(over_water_migrants_ol$ID),as.numeric(flight_record$ID))
        #Update location
        flight_record[replace_index,c("new_lon","new_lat")] <- over_water_migrants_ol[,c("new_lon","new_lat")]
        flight_record[replace_index,"dist"] <- over_water_migrants_ol[,"dist"]
        flight_record[replace_index,"water"] <- 0
        # See which birds didn't have any draws that were over land. These birds
        # will be forced to stopover
        force_stopover_birds <- unique(over_water_migrants_rep$ID[ ! over_water_migrants_rep$ID %in% trim_list$ID])
        
      } else {
        #All birds need to be go into stopover
        force_stopover_birds <- unique(over_water_migrants_rep$ID)
      }
      
      #Only force stopover if there are birds that need to be!
      if (length(force_stopover_birds)>0){
        #Force a stopover of these birds
        over_water_migrants_fs <- filter(over_water_migrants,ID %in% force_stopover_birds)
        over_water_migrants_fs <- initiate_stopover(over_water_migrants_fs,force=TRUE)
        #Now match to main df and update stopover counter, don't update location
        replace_index <- match(as.numeric(over_water_migrants_fs$ID),as.numeric(flight_record$ID))
        flight_record[replace_index,"stopover_counter"] <- over_water_migrants_fs[,"stopover_counter"]
        flight_record[replace_index,c("new_lon","new_lat")]  <- flight_record[replace_index,c("current_lon","current_lat")]
        flight_record[replace_index,"dist"] <- 0
        flight_record[replace_index,"water"] <- 1
      }
      
      replace_index <- match(as.numeric(flight_record$ID),as.numeric(current_migrants_to_move$ID))
      #Update location
      current_migrants_to_move[replace_index,c("current_lon","current_lat")] <- flight_record[,c("new_lon","new_lat")]
      #Update energy
      current_migrants_to_move[replace_index,c("energy")] <- current_migrants_to_move[replace_index,c("energy")] - flight_record[,"dist"]
      current_migrants_to_move[replace_index,c("stopover_counter")] <- flight_record[,c("stopover_counter")]
    }
    
    
  }
  
  #Return updated bird records
  return(current_migrants_to_move)
}

#Long function for running a day in the model 
run_day <- function(updatable_df,static_df){
  
  #Advance day
  updatable_df[,"day"] <- updatable_df[,"day"] + 1
  
  if (updatable_df[1,"day"] == 366){
    updatable_df[,"day"]  <- 1
  }
  
  #Update seasons for birds starting migration
  updatable_df <- start_migration(updatable_df,static_df)
  
  #Split migrating birds (seasons are 1:4 - breeding, fall migration, wintering, spring migration)
  current_migrants <- filter(updatable_df,season == 4 | season == 2)
  not_migrating_birds <- filter(updatable_df,season == 1 | season == 3)
  
  
  #Code for migrants - assumption here is that breeding/wintering birds stay put.
  #Only run if there are migrating birds
  if (nrow(current_migrants)>0){
    
    #Attach relevant info - migration seasons are coded different so the static info needs
    #to match up. 
    
    #SPEED. This may be one section of the code that could be streamlined
    
    #Get the IDS of birds in either migration season - ids should directly match to rows in static df
    spr_indices <- as.numeric(current_migrants[which(current_migrants$season == 4),"ID"])
    fall_indices <- as.numeric(current_migrants[which(current_migrants$season == 2),"ID"])
    
    #Get the static information
    static_info_spr <- static_df[spr_indices,c("ID","breeding_lon","breeding_lat","speed_mean_s",
                                               "speed_sd_s","bear_err_mean_s","bear_err_sd_s",
                                               "recovery_rate_s","max_mig_time_s","max_energy_s","goal_radius")]
    
    static_info_fall <- static_df[fall_indices,c("ID","nonbreeding_lon","nonbreeding_lat","speed_mean_f",
                                                 "speed_sd_f","bear_err_mean_f","bear_err_sd_f",
                                                 "recovery_rate_f","max_mig_time_f","max_energy_f","goal_radius")]
    colnames(static_info_spr) <- c("ID","goal_lon","goal_lat","speed_mean",
                                   "speed_sd","bear_err_u","bear_err_sd","recovery_rate",
                                   "migr_max_time","max_energy","goal_radius")
    colnames(static_info_fall) <- c("ID","goal_lon","goal_lat","speed_mean",
                                    "speed_sd","bear_err_u","bear_err_sd","recovery_rate",
                                    "migr_max_time","max_energy","goal_radius")
    
    static_info <- bind_rows(static_info_spr,static_info_fall)
    
    #Combine into current migrants df, with both spring and fall migrants
    current_migrants <- merge(current_migrants,static_info,by="ID")
    
    #Update migration counter
    current_migrants$migr_counter <- current_migrants$migr_counter + 1
    
    #If the bird has exceeded max migration time, it passes into the next season
    over_time_limit <- current_migrants$migr_counter > current_migrants$migr_max_time
    
    #Split the birds that have timed out from the rest
    birds_over_time_limit <- current_migrants[over_time_limit,]
    
    #Update season for these birds
    birds_over_time_limit <- update_season(birds_over_time_limit)
    
    #Save these records for the end of this function, for now just cut them out
    current_migrants <- current_migrants[!over_time_limit,]
    
    # If the bird has passed within the goal radius, it "tracks" this day to its goal point
    # Migration is considered over
    current_locations <- cbind(current_migrants$current_lon,current_migrants$current_lat)
    goal_locations <- cbind(current_migrants$goal_lon,current_migrants$goal_lat)
    goal_radii <- current_migrants$goal_radius

    at_goal_binary <- at_goal(current_locations,goal_locations,goal_radii)

    #New location = goal location, migration over
    birds_arrived <- current_migrants[at_goal_binary,]
    
    #Update season
    birds_arrived <- update_season(birds_arrived)
    
    #Update location - bird "homes" to goal location
    birds_arrived[,c("current_lon","current_lat")] <- birds_arrived[,c("goal_lon","goal_lat")]
    
    #Save these records for the end of this function, for now just cut them out
    current_migrants <- current_migrants[!at_goal_binary,]
    
    #Update stopover birds
    
    #Energy - current energy + recovery rate * max energy
    current_migrants[which(current_migrants$stopover_counter>0),"energy"] <- 
      current_migrants[which(current_migrants$stopover_counter>0),"energy"] + 
      current_migrants[which(current_migrants$stopover_counter>0),"recovery_rate"] * 
      current_migrants[which(current_migrants$stopover_counter>0),"max_energy"]
    
    #Counter = counter - 1
    current_migrants[which(current_migrants$stopover_counter>0),"stopover_counter"] <- 
      current_migrants[which(current_migrants$stopover_counter>0),"stopover_counter"] - 1
    
    #For all migrants that are not in stopover, check to see if they should start stopover
    #Here, there is a small chance that birds that have just ended a stopover may start a
    #one day stopover, but this will be very uncommon
    #This initiates stopover for some birds.
    
    current_migrants <- initiate_stopover(current_migrants)

    #Filter birds that are not in stopover
    current_migrants_to_move <- filter(current_migrants,stopover_counter == 0)
    current_migrants_stopover <- filter(current_migrants,stopover_counter != 0)
    
    #Move birds, update location/energy, force stopover for overwater birds
    current_migrants_to_move <- move_bird(current_migrants_to_move)
    
    #Combine records for birds that stopped over timed out, arrived at destination, and are still moving
    #SPEED - this will be slow. Might be worth pre-assigning dataframe
    current_migrants_to_move <- bind_rows(current_migrants_to_move,current_migrants_stopover,birds_arrived,birds_over_time_limit)
    
    #Truncate current_migrants_to_move df to updatable format, sort by ID return
    current_migrants_to_move <- current_migrants_to_move[,c("ID","current_lon","current_lat","day","season","energy","stopover_counter","migr_counter")]
    
    #Add non-migrating birds
    current_migrants_to_move <- bind_rows(current_migrants_to_move,not_migrating_birds)
    
    #Sort
    updatable_df <- current_migrants_to_move[order(current_migrants_to_move$ID),]
    
    return(updatable_df)
  } else {
    
    return(updatable_df)
  }

}

#This function runs the model and just spits out the error for each week
# This is the function useful for training models
model_w_error <- function(initial_dfs,start_date,grid_size,err_grid,modeled_area){
  
  #Get dataframes in
  static_info <- as.data.frame(initial_dfs[1])
  updatable_info <- as.data.frame(initial_dfs[2])
  
  #Create dataframe that contains all information about the simulation, 
  # this resets each week
  # SPEED - this can probably be improved with a dataframe of preassigned length
  all_days_rec <- updatable_info
  
  #Construct grids for analysis - this creates a grid system where the centers
  #of the grids are "spacing" distance apart. 
  dggs <- dgconstruct(spacing=grid_size, metric=TRUE, resround='down')
  
  #Get the start date of the model input as the week
  week_starts <- seq(1,365,by=7)
  week_count = which(week_starts==start_date)
  
  #Set up a counter for the days in the week
  day_in_week_count <- 0
  
  #Create a record for the weekly error
  week_error_record <- rep(NA,52)
  bc_err_list <- list()
  #Run multiple days of the model
  #For each day in a year (not calendar year, but a year from the start date supplied)
  for (n in 1:365){
    
    #Run model, save results
    updatable_info <- run_day(updatable_info,static_info)
    all_days_rec <- bind_rows(all_days_rec,updatable_info)
    
    #For the sake of the training set being run on weeks, when the model reaches day 365, 
    #the output of the model is ignored and no validation/day in week is updated
    if (updatable_info[1,4]!=365){
      #Update day counter
      day_in_week_count <- day_in_week_count + 1
      if (day_in_week_count == 7){
        
        #Calculate model accuracy for each week
        
        #Read in relevant ebird derived data
        birdcounts_ebird <- weekly_err[[week_count]]
        modeled_area_w <- modeled_area[[week_count]] 
        
        #Do the same for the modeled data
        all_days_rec$cell <- dgGEO_to_SEQNUM(dggs,all_days_rec$current_lon,all_days_rec$current_lat)$seqnum
        
        birdcounts_mibm <- all_days_rec %>% 
          group_by(cell) %>% 
          summarise(count=n(),.groups="drop")
        
        #Only keep records within the modeled range
        birdcounts_mibm <- as.data.frame(filter(birdcounts_mibm,cell %in% modeled_area_w))
        
        #Convert to % of population within modeled range
        birdcounts_mibm$count <- birdcounts_mibm$count/(sum(birdcounts_mibm$count))
        
        #Merge, calculate error
        birdcounts_comb <- merge(birdcounts_ebird,birdcounts_mibm,by="cell",all=TRUE)
        
        #replace NA values with zero
        birdcounts_comb[is.na(birdcounts_comb)] <- 0
        
        #Get error - sd from mean
        birdcounts_comb$comb_e <- abs(birdcounts_comb$perc_pop-birdcounts_comb$count)/birdcounts_comb$err_adj
        
        #Use sigmoid to minimize runaway error
        birdcounts_comb$comb_e <- sigmoid(birdcounts_comb$comb_e/2)
        
        #Correct for number of cells with birds in them
        birdcounts_comb$comb_e <- birdcounts_comb$comb_e/nrow(birdcounts_comb)
        
        bc_err_list[[week_count]] <- birdcounts_comb 
        
        week_error <- sum(birdcounts_comb$comb_e)
        week_error_record[week_count] <- week_error
        
        #Update week
        week_count <- week_count + 1
        if (week_count > 52){
          week_count <- 1
        }
        #Reset day counter
        day_in_week_count <- 0
        #Reset updatable info
        all_days_rec <- data.frame(data=NULL)
      }
    }
  }
  return(list(week_error_record,bc_err_list))
}

#This function takes in the updated daily dataframe, and returns a dataframe with just the
# locations of the birds, day and their IDs
save_locations <- function(updatable_df){
  locations_IDs <- updatable_df %>% select(ID,current_lon,current_lat,day)
  return(locations_IDs)
}

