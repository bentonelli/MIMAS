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
  
  # This is where migratory connectivity comes in. Migratory connectivity is modeled
  # here as a function of the distance between breeding and non-breeding points.
  # Ultimately, two types of migr. connectivity are modeled - chain migration (-1) and leapfrog (1)
  # Random assortment exists at the center of "perfect" leapfrog and chain migration
  
  #Get distances in km
  ll_b <- matched_points[,2:3]
  ll_nb <- matched_points[,4:5]
    
  distance_matrix <- distm(ll_b,ll_nb)/1000
  
  #Get average distance of all breeding points to all other non-breeding points
  dist_avg_points <- rowMeans(distance_matrix)
  
  if(mig_con < 0){
    #Rank where the values with largest values get highest priority - chain
    rnk_type <- rank(-dist_avg_points,ties.method = "random")
  } else {
    #Rank where the values with largest values get lowest priority - leapfrog
    rnk_type <- rank(dist_avg_points,ties.method = "random")
  }
  
  #Convert ranks, add noise
  rank_adj <- (rnk_type/length(dist_avg_points))*abs(mig_con) + rnorm(length(dist_avg_points),0,1-abs(mig_con))
  
  #Set final priority ranks
  fin_ranks <- rank(rank_adj)
  
  #For each point, get best available match
  match_list <- c()
  for (nn in 1:length(dist_avg_points)){
    ind_by_rank <- which(fin_ranks == nn)
    #Get costs for all points
    min_match <- min(distance_matrix[ind_by_rank,],na.rm=TRUE)
    min_ind <- which(distance_matrix[ind_by_rank,]==min_match)
    
    ml_add <- c(ind_by_rank,min_ind)
    match_list <- rbind(match_list,ml_add)
    #NA out col
    distance_matrix[,min_ind] <- NA
  }
  #save memory
  distance_matrix <- NULL
  
  sort_br_p <- ll_b[match_list[,1],]
  sort_nb_p <- ll_nb[match_list[,2],]
  
  #Rescramble points in order to make sure theres no corr with ID and location
  rand_sort_ind <- sample(1:num_pulls)
  matched_points[,2:3] <- sort_br_p[rand_sort_ind,]
  matched_points[,4:5] <- sort_nb_p[rand_sort_ind,]
  
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
  
  if (migr_timing_lat_f != 0 | migr_timing_lat_s != 0){
    ranks <- 1:num_pulls
    
    probs_s <- (ranks/num_pulls)*abs(migr_timing_lat_s)
    probs_f <- (ranks/num_pulls)*abs(migr_timing_lat_f)
    
    #Draw ranks for spring, fall
    new_ind_s <- rnorm(num_pulls,0,1-abs(migr_timing_lat_s))
    
    if (migr_timing_lat_s > 0){
      new_ind_s <- rank(probs_s+new_ind_s)
    } else {
      new_ind_s <- rank(-(probs_s+new_ind_s))
    }
    
    new_ind_f <- rnorm(num_pulls,0,1-abs(migr_timing_lat_f))
    if (migr_timing_lat_f > 0){
      new_ind_f <- rank(probs_f+new_ind_f)
    } else {
      new_ind_f <- rank(-(probs_f+new_ind_f))
      
    }
    
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
  
  # Here, we can force start dates in the fall to not be before the start date of the
  # simulation. This forces those birds that would be really early migrants to instead 
  # just very early migrants! Another option to fix this is to start the simulation earlier.
  matched_points$start_date_f[which(matched_points$start_date_f < (start_date+1))] <- start_date + 1
  
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
  todays_speed[which(todays_speed > distances_to_goal)] <- 
    distances_to_goal[which(todays_speed > distances_to_goal)]
  return(todays_speed)
}

### Function for overwater check. This operates with the assumption that birds
# should not be able to to land overwater (something to change if the model is 
#setup for waterbirds)
over_water <- function(locations){
  #Get locations
  locations <- as.data.frame(locations)
  
  # A failure of the Rhumb line system may throw an error normally, so here 
  # we change any NaNs (which indicate this failure) to the coordinates (-150,90 - the North Pole)
  # which is over water.
  if (sum(is.na(locations[,1]))>0){
    na_ind <- which(is.na(locations$lon))
    locations[na_ind,1] <- -150
    locations[na_ind,2] <- 90
  }
  #Convert points to sf format
  pts <- suppressMessages(st_as_sf(locations, coords=1:2, crs=4326))
  
  #This will return the values that are over water/land (TRUE,FALSE) and should be repeated
  #this function  has given me issues previously, and is a good place to start debugging :)
  water <- is.na(as.numeric(suppressMessages(st_intersects(pts, world))))
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
    
    #Pull distance to goal, convert to km
    distances_to_goal <- distRhumb(current_locations,goal_locations)/1000
    
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
                                     cbind(over_water_migrants_rep$goal_lon,over_water_migrants_rep$goal_lat))/1000
      
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
    }
    replace_index <- match(as.numeric(flight_record$ID),as.numeric(current_migrants_to_move$ID))
    #Update location
    current_migrants_to_move[replace_index,c("current_lon","current_lat")] <- flight_record[,c("new_lon","new_lat")]
    #Update energy
    current_migrants_to_move[replace_index,c("energy")] <- current_migrants_to_move[replace_index,c("energy")] - flight_record[,"dist"]
    current_migrants_to_move[replace_index,c("stopover_counter")] <- flight_record[,c("stopover_counter")]
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
# This is the function is built for training models
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
  dggs <- dgconstruct(spacing=grid_size, metric=TRUE, resround='down',show_info = FALSE)
  
  #Get the start date of the model input as the week
  week_starts <- seq(1,365,by=7)
  week_count = which(week_starts==start_date)
  
  #Set up a counter for the days in the week
  day_in_week_count <- 0
  
  #Create a record for the weekly error
  week_error_record <- rep(NA,52)
  week_within_range <- rep(NA,52)
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
        
        #Calculate model error for each week
        
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
        birdcounts_comb$comb_e <- sigmoid(abs(birdcounts_comb$perc_pop-birdcounts_comb$count)/birdcounts_comb$err_adj)
        
        #Correct for number of cells with birds in them
        birdcounts_comb$comb_e <- birdcounts_comb$comb_e/nrow(birdcounts_comb)
        
        bc_err_list[[week_count]] <- birdcounts_comb 
        
        week_error <- sum(abs(birdcounts_comb$perc_pop - birdcounts_comb$count))/sum(birdcounts_comb$err_adj)
        week_error_record[week_count] <- week_error
        
        week_within_range[week_count] <- birdcounts_comb %>% 
          filter(perc_pop >0) %>%
          summarise(count_in_range <- sum(count))
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
  return(list(week_error_record,week_within_range,bc_err_list))
}

#This function takes in the updated daily dataframe, and returns a dataframe with just the
# locations of the birds, day and their IDs
save_locations <- function(updatable_df,save_season=FALSE){
  if (save_season){
    locations_IDs <- updatable_df %>% select(ID,current_lon,current_lat,day,season)
  } else {
    locations_IDs <- updatable_df %>% select(ID,current_lon,current_lat,day)
  }
  return(locations_IDs)
}

# This function reads in the simulation output of fit models, and returns the top
# models based on the sum_err value of each model

# sim_output_loc should be directory location of simulation output
# perc_plus should be percent of additional models accepted (expressed as .03 = 3%,etc)
# total_sims should be the total number of simulations to consider

best_models <- function(sim_output_loc, perc_plus,total_sims){
  #Read in files
  training_files <- lapply(Sys.glob(sim_output_loc), read.csv)
  #Add to df
  compiled_df <- do.call(bind_rows,training_files)
  #Select only the first XX number of simulation
  compiled_df <- compiled_df[1:total_sims,]
  #Get best models
  perc_lim_plus <- quantile(compiled_df$sum_err,probs = c(.01)) * (1+perc_plus)
  one_perc_plus_sims <- compiled_df[which(compiled_df$sum_err<=perc_lim_plus),]
  return(one_perc_plus_sims)
}

# This analysis function creates a pdf and csv file demonstrating the convergence 
# (or lack therof) of each parameter.
param_converge <- function(spec_abr,run_date,output_loc=paste("data/output/",spec_abr,"_",run_date,sep=""),
                           num_sims_in_file=100,total_sims = 100000,buffer = .03){
  library(dplyr)
  library(ggplot2)
  library(readr)
  
  s1_loc <- paste("data/training_runs/",spec_abr,"_",run_date,"/S1/*.csv",sep="")
  s2_loc <- paste("data/training_runs/",spec_abr,"_",run_date,"/S2/*.csv",sep="")
  #Read in first fit parameter set
  step1_param_set <- lapply(Sys.glob(s1_loc), read.csv)
  step1_compiled_df <- do.call(bind_rows,step1_param_set[1:(total_sims/num_sims_in_file)])
  step1_compiled_df <- step1_compiled_df[complete.cases(step1_compiled_df),]
  
  # Read in second fit information into single dataframe, then process and analyze
  step2_param_set <- lapply(Sys.glob(s2_loc), read.csv)
  step2_compiled_df <- do.call(bind_rows,step2_param_set[1:(total_sims/num_sims_in_file)])
  step2_compiled_df <- step2_compiled_df[complete.cases(step2_compiled_df),]
  #Get best of the second step runs
  perc_lim_plus <- quantile(step2_compiled_df$sum_err,probs = c(.01)) *(1+buffer)
  one_perc_plus_sims <- step2_compiled_df[which(step2_compiled_df$sum_err<=perc_lim_plus),]
  
  pdf(file = paste(output_loc,"/",spec_abr,"_convergence.pdf",sep=""),   # The directory you want to save the file in
      width = 4, # The width of the plot in inches
      height = 4) # The height of the plot in inches
  
  param_convergence_rec <- c()
  for (n in 2:26){
    set_range <- density(step1_compiled_df[,n],na.rm = TRUE)
    xx_max <- max(set_range$x)
    xx_min <- min(set_range$x)
    set_bw <- density(one_perc_plus_sims[,n],na.rm = TRUE)
    set_bw <- set_bw$bw
    
    #Get 95% range estimate for each stage of the fitting process
    og_range <- abs(quantile(step1_compiled_df[,n],probs=c(.025)) - quantile(step1_compiled_df[,n],probs=c(.975)))
    first_range <- abs(quantile(step2_compiled_df[,n],probs=c(.025)) - quantile(step2_compiled_df[,n],probs=c(.975)))
    second_range <- abs(quantile(one_perc_plus_sims[,n],probs=c(.025)) - quantile(one_perc_plus_sims[,n],probs=c(.975)))
    
    #Get the total convergence from the first to last step
    tot_range_convergence <- 1- (second_range/og_range)
    
    #Get median estimate of each parameter in each step, with 95% CI for final
    median_est_prior <- median(step1_compiled_df[,n])
    median_est_s1 <- median(step2_compiled_df[,n])
    
    median_est_s2 <- median(one_perc_plus_sims[,n])
    lower_est <- quantile(one_perc_plus_sims[,n],probs=c(.025)) 
    upper_est <- quantile(one_perc_plus_sims[,n],probs=c(.975)) 
    
    add_v <- c(colnames(one_perc_plus_sims)[n],og_range,first_range,second_range,
               median_est_prior,median_est_s1,median_est_s2,lower_est,upper_est,
               tot_range_convergence)
    param_convergence_rec <- rbind(param_convergence_rec,add_v)
    
    plot(density(one_perc_plus_sims[,n]),col="dodgerblue4",lwd=3,xlim=c(xx_min,xx_max),main=colnames(one_perc_plus_sims)[n])
    points(density(step2_compiled_df[,n],bw=set_bw),col="firebrick3",lwd=3,type="l")
    points(density(step1_compiled_df[,n],bw=set_bw),col="black",lwd=3,type="l")
  }
  dev.off()
  colnames(param_convergence_rec) <- c("param_name","prior range","S1 range","S2 range",
                                       "prior med. est","S1 med. est.","S2 med. est","2.5%",
                                       "97.5%","% Converg.")
  param_convergence_rec <- as.data.frame(param_convergence_rec)
  write_csv(param_convergence_rec,paste(paste("data/output/",spec_abr,"_",run_date,"/param_convergence.csv",sep="")))
}

#Code that extracts parameter sets associated with the best simulation runs
best_sims <- function(spec_abr,run_date,output_loc=paste("data/output/",spec_abr,"_",run_date,sep=""),
                      output_loc2 = paste("data/param_sets/",spec_abr,"/",run_date,sep=""),
                      num_sims_in_file=100,total_sims = 100000,buffer = .03){
  library(dplyr)
  print(paste(output_loc,"/",spec_abr,"_best.rds",sep=""))
  print(paste(output_loc2,"/",spec_abr,"_S2.rds",sep=""))
  
  #If output folder doesn't exist, create one
  if(!dir.exists(paste("data/output/",spec_abr,sep=""))){
    dir.create(paste("data/output/",spec_abr,sep=""))
  }
  
  #If output folder doesn't exist, create one
  if(!dir.exists(paste("data/output/",spec_abr,"_",run_date,sep=""))){
    dir.create(paste("data/output/",spec_abr,"_",run_date,sep=""))
  }
  
  #Location of saved simulation files
  sims_location <- paste("data/training_runs/",spec_abr,"_",run_date,"/S2/*.csv",sep="")
  
  #Define buffer for error - here set to 3%
  sim_best_models <- best_models(sims_location,buffer,total_sims)
  
  #Save to RDS in two locations
  saveRDS(sim_best_models,paste(output_loc,"/",spec_abr,"_best.rds",sep=""))
  saveRDS(sim_best_models,paste(output_loc2,"/",spec_abr,"_S2.rds",sep=""))
  
}

#Code to create parameter correlation plot
param_corr <- function(spec_abr,run_date,output_loc=paste("data/output/",spec_abr,"_",run_date,sep="")){
  # Script to compare parameter correlations across best fit models
  library(corrplot)
  
  #If output folder doesn't exist, create one
  if(!dir.exists(paste("data/output/",spec_abr,"_",run_date,sep=""))){
    dir.create(paste("data/output/",spec_abr,"_",run_date,sep=""))
  }
  
  #Read in parameter set
  fit_params <- readRDS(paste("data/output/",spec_abr,"_",run_date,"/",spec_abr,"_best.rds",sep=""))
  
  plot(fit_params$speed_mean_s,fit_params$start_date_u_s,cex=.1)
  wanted_corr_params <- as.matrix(fit_params[,2:26])
  cor_data <- round(cor(wanted_corr_params),2)
  
  #Remove NA rows
  cor_data <- cor_data[-as.numeric(which(rowSums(cor_data,na.rm=TRUE)==TRUE)),-as.numeric(which(rowSums(cor_data,na.rm=TRUE)==TRUE))]
  
  #Plot correlations. Note here that some parameters might be not be continuous,
  #so interpret with caution!
  corrplot(cor_data, method = 'color',order = 'AOE')
  
  pdf(file = paste(output_loc,"/",spec_abr,"_correlation.pdf",sep=""),   # The directory you want to save the file in
      width = 8, # The width of the plot in inches
      height = 8) # The height of the plot in inches
  corrplot(cor_data, method = 'color',order = 'AOE')
  dev.off()
  
}

#'This script takes in the model simulation runs from the initial parameter space
#'exploration, performs rejection-ABC for each season - fall and spring - 
#'and then resamples from the posteriors to get new parameter sets. These new 
#'parameter sets can then be used to perform another cycle of rejection ABC to 
#'arrive at better estimates for a given species. The output of this script
#'is a list of top models for fall, spring and both combined
recombine <- function(spec_abr,run_date,output_loc=paste("data/param_sets/",spec_abr,"/",run_date,sep="")){
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  # Read in training information into single dataframe, then process and analyze
  n_files <- 1000
  n_in_each <- 100
  
  if(!dir.exists(output_loc)){
    dir.create(output_loc)
  }
  
  #Set location to save file
  output_loc_adj <- paste(output_loc,"/",spec_abr,"_S1.rds",sep = "")
  
  #Set location of simulation files
  training_files <- lapply(Sys.glob(paste("data/training_runs/",spec_abr,"_",run_date,"/S1/*.csv",sep="")), read.csv)
  
  compiled_df <- do.call(bind_rows,training_files)
  
  #If there was an issue with a simulation causing and NA, remove those here
  
  compiled_df <- compiled_df[complete.cases(compiled_df),]
  
  # Do to some simulations failing, the number of actual simulations likely exceeds
  # the number wanted/required. Here, this code truncates to the first ### of simulations
  compiled_df <- compiled_df[1:(n_files*n_in_each),]
  
  #Get the start date of the model run -- the -1 accounts for indexing
  m_start_date <- as_date(compiled_df$start_date[1]-1)
  #This should be a day that starts a week
  m_start_week <- week(m_start_date)
  
  all_weeks <- seq(1,52,by=1)
  fall_weeks <- (seq(m_start_week,m_start_week+25,by=1)-1)%%52 + 1
  spring_weeks <- all_weeks[!all_weeks %in% fall_weeks]
  #Calculate error terms for spring, fall
  #Split and sum error for 2 sections of the year, weeks 1:26, and 27:52
  err_all <-  compiled_df[1:nrow(compiled_df),29:80]
  err_s <- rowSums(err_all[,spring_weeks])
  err_f <- rowSums(err_all[,fall_weeks])
  
  #add to df, cut weekly error
  compiled_df <- compiled_df[,1:28]
  compiled_df$err_spring <- err_s
  compiled_df$err_fall <- err_f
  
  # Get 1% best models from both + models with error within 5% of the 'worst' model
  # in that bunch. Adding back more models is necessary in order to
  top_f_models <- filter(compiled_df,quantile(compiled_df$err_fall,.01)*1.05 > err_fall)
  
  top_s_models <- filter(compiled_df,quantile(compiled_df$err_spring,.01)*1.05 > err_spring)
  
  # Run the recombination step - take fall and spring model parameter sets, and combine
  # them to create a list of simulation parameter sets. Because there are a large number 
  # of potential combinations, this process can be done randomly
  fall_param_inx <- c(4,5,8,9,11,14,15,17,19,26)
  spring_param_inx <- c(2,3,6,7,10,12,13,16,18,25)
  shared_inx <- c(20,21,22,23,24)
  colnames(top_f_models)[fall_param_inx]
  colnames(top_s_models)[spring_param_inx]
  colnames(top_f_models)[shared_inx]
  
  #Combine best models from spring and fall for use later
  top_f_s_models <- rbind(top_f_models,top_s_models)
  three_batch_out <- list()
  three_batch_out[[1]] <- top_f_models
  three_batch_out[[2]] <- top_s_models
  three_batch_out[[3]] <- top_f_s_models
  
  saveRDS(three_batch_out,output_loc_adj)
}

derived_params <- function(param_file,path = TRUE){
  #Script to get derived parameters - mean speed km/day and recovery rate, km/day
  # 8-26-22: Adding combined migratory connectivity and strength 
  library(truncnorm)
  
  if (path){
    #Read in parameter set
    param_set <- readRDS(param_file)
  } else {
    param_set <- param_file
  }
  
  ### Spring ####
  speed_theor_weighted <- c()
  prob_theor <- c()
  for (n in 0:9){
    prob <- rep(1,length(param_set$max_energy_s))
    for (m in 0:n){
      prob <- prob * (param_set$max_energy_s-m*param_set$speed_mean_s)/param_set$max_energy_s
      prob[prob<0] <- 0
    }
    prob_theor <- cbind(prob_theor,prob)
    mean_speed <- etruncnorm(a=0,b=(param_set$max_energy_s-n*param_set$speed_mean_s),
                             mean = param_set$speed_mean_s,sd = param_set$speed_sd_s)
    mean_speed[mean_speed<0] <- 0
    speed_theor_weighted <- cbind(speed_theor_weighted,mean_speed)
  }
  spring_speed_adj <- (rowSums(speed_theor_weighted*prob_theor)/rowSums(prob_theor))
  
  ### Fall ####
  speed_theor_weighted <- c()
  prob_theor <- c()
  for (n in 0:9){
    prob <- rep(1,length(param_set$max_energy_f))
    for (m in 0:n){
      prob <- prob * (param_set$max_energy_f-m*param_set$speed_mean_f)/param_set$max_energy_f
      prob[prob<0] <- 0
    }
    prob_theor <- cbind(prob_theor,prob)
    mean_speed <- etruncnorm(a=0,b=(param_set$max_energy_f-n*param_set$speed_mean_f),
                             mean = param_set$speed_mean_f,sd = param_set$speed_sd_f)
    mean_speed[mean_speed<0] <- 0
    speed_theor_weighted <- cbind(speed_theor_weighted,mean_speed)
  }
  fall_speed_adj <- (rowSums(speed_theor_weighted*prob_theor)/rowSums(prob_theor))
  
  #Recovery rate
  rec_rate_f <- (param_set$max_energy_f*param_set$recovery_rate_f)
  rec_rate_s <- (param_set$max_energy_s*param_set$recovery_rate_s)
  
  # Migratory connectivity type x strength 
  # First converts type to -1,1, then multiplies to get continuous gradient
  mig_con_str <- param_set$mig_con * ((param_set$mig_con_type - 1.5)*2)
  
  
  r_list <- list(spring_speed_adj,fall_speed_adj,rec_rate_s,rec_rate_f,mig_con_str)
  names(r_list) <- c("spring_speed_adj","fall_speed_adj",
                     "rec_rate_s","rec_rate_f","mig_con_str")
  return(r_list)
  
}

#The function below is adapted from the function of the same name in the currently
#broken package, BEST
postPriorOverlap <-
  function( paramSampleVec, prior, ..., yaxt="n", ylab="",
            xlab="Parameter", main="", cex.lab=1.5, cex=1.4,
            xlim=range(paramSampleVec), breaks=NULL,
            mainColor="dodgerblue4", priorColor="grey", overlapColor="tan4") {
    
    # Does a posterior histogram for a single parameter, adds the prior,
    #   displays and calculates the overlap.
    # Returns the overlap.
    
    oldpar <- par(xpd=NA) ; on.exit(par(oldpar))
    
    # get breaks: a sensible number over the hdi; cover the full range (and no more);
    #   equal spacing.
    if (is.null(breaks)) {
      nbreaks <- ceiling(diff(range(paramSampleVec)) / as.numeric(diff(hdi(paramSampleVec))/18))
      breaks <- seq(from=min(paramSampleVec), to=max(paramSampleVec), length.out=nbreaks)
    }
    # plot posterior histogram.
    histinfo <- hist(paramSampleVec, xlab=xlab, yaxt=yaxt, ylab=ylab,
                     freq=FALSE, border='white', col=mainColor,
                     xlim=xlim, main=main, cex=cex, cex.lab=cex.lab,
                     breaks=breaks)
    
    if (is.numeric(prior))  {
      # plot the prior if it's numeric
      priorInfo <- hist(prior, breaks=c(-Inf, breaks, Inf), add=TRUE,
                        freq=FALSE, col=priorColor, border='white')$density[2:length(breaks)]
    } else if (is.function(prior)) {
      if(class(try(prior(0.5, ...), TRUE)) == "try-error")
        stop(paste("Incorrect arguments for the density function", substitute(prior)))
      priorInfo <- prior(histinfo$mids, ...)
    }
    # get (and plot) the overlap
    minHt <- pmin(priorInfo, histinfo$density)
    rect(breaks[-length(breaks)], rep(0, length(breaks)-1), breaks[-1], minHt, col=overlapColor,
         border='white')
    overlap <- sum(minHt * diff(histinfo$breaks))
    # Add curve if prior is a function
    if (is.function(prior))
      lines(histinfo$mids, priorInfo, lwd=2, col=priorColor)
    # Add text
    text(mean(breaks), 0, paste0("overlap = ", round(overlap*100), "%"), pos=3, cex=cex)
    
    return(overlap)
  }
