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
library(zoo)
sf::sf_use_s2(FALSE)

setwd("~/Documents/Coding/R/MIMAS/")
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
#Set seed
set.seed(188)

#Get functions
source(file = "scripts/sim_functions/sim_funct.R")

# Parameter sets should be drawn from a script in the data folder, or from a 
# fit parameter set
param_set <- readRDS("data/param_sets/CCSP/9_20_22/CCSP_S2.rds")
for (each_var in 1:ncol(param_set)){
  #Get best
  assign(colnames(param_set)[each_var],param_set[which(param_set$sum_err == min(param_set$sum_err)),each_var])
}

#Define species
species_target <- "Clay-colored Sparrow"

#Set number of birds to include in simulation
num_pulls <- 2000

#Define spacing of error maps
spacing <- 200

#Set number of days
num_days <- 365*2

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

locs_rec <- as.data.frame(matrix(NA,ncol = 4,nrow = (num_days*num_pulls)))
colnames(locs_rec) <- c("ID","Lon","Lat","Day")
for (days in 1:num_days){
  print(days)
  upd_df <- run_day(upd_df,static_df)
  summary(upd_df)
  locs <- save_locations(upd_df)
  loc_ind_s <- (days - 1)*num_pulls + 1 
  loc_ind_e <- (days - 1)*num_pulls + num_pulls 
  locs_rec[loc_ind_s:loc_ind_e,] <- locs
}

### Simulate pathogen spread continuous version ####
set.seed(190)
N <- num_pulls
inf_time <- 30
rec_time <- 90
disease_radius <- 100
lambda <- 1/inf_time
rnaught <- 60
beta <- (rnaught*lambda)/N
start_inf_num <- 50
start_num_rec <- 200

locs_rec_p <- locs_rec

#Add column to denote simulation time
locs_rec_p$sim_days <- rep(1:num_days, each=num_pulls)

#Add infection status to new column
locs_rec_p$inf <- "S"

#Add infection timer to new column
locs_rec_p$timer <- NA

#Infect starting birds
locs_rec_p$inf[1:start_inf_num] <- "I"
locs_rec_p$timer[1:start_inf_num] <- sample(1:inf_time,start_inf_num,replace=TRUE)

#Recover birds
locs_rec_p$inf[start_inf_num+1:start_num_rec+1+start_num_rec] <- "R"
locs_rec_p$timer[start_inf_num+1:start_num_rec+1+start_num_rec] <- sample(1:rec_time,start_num_rec,replace=TRUE)

#Get first day
start_day <- locs_rec_p$Day[1]
end_day <- locs_rec_p$Day[nrow(locs_rec_p)]

day_list <- unique(locs_rec_p$sim_days)

#Get breeding, non-breeding locations
breed_loc <- select(static_df,c(ID,breeding_lon,breeding_lat))
nbreed_loc <- select(static_df,c(ID,nonbreeding_lon,nonbreeding_lat))

# For each day in simulation
index_count <- 0
migr_tracker <- c()
lat_tracker <- c()

for (each_day in day_list){
  print(each_day)
  
  #Get where next day data exists
  next_day_indices <- which(locs_rec_p$sim_days==(each_day+1))
  
  #Extract df of birds for a given day
  day_data <- filter(locs_rec_p, sim_days==each_day)
  
  #Get number of birds migrating by comparing to breed, non breeding points
  current_lons <- select(day_data,c(ID,Lon))
  current_lats <- select(day_data,c(ID,Lat))
  breed_num <- sum(current_lons$Lon %in% breed_loc$breeding_lon & current_lats$Lat %in% breed_loc$breeding_lat)
  nbreed_num <- sum(current_lons$Lon %in% nbreed_loc$nonbreeding_lon & current_lats$Lat %in% nbreed_loc$nonbreeding_lat)
  
  migr_num <- num_pulls - (breed_num + nbreed_num)
  migr_tracker <- c(migr_tracker,migr_num)
  
  #Get latitude of birds
  lat_tracker <- c(lat_tracker, mean(current_lats$Lat))
  
  if (sum(day_data$inf == "I") > 0 & sum(day_data$inf=="S") > 0){
    inf_day_locs <- filter(day_data,inf=="I") %>% select(Lon,Lat)
    sus_day_locs <- filter(day_data,inf=="S") %>% select(Lon,Lat)
    sus_day_ids <- filter(day_data,inf=="S") %>% select(ID)
    
    dist_i_s <- distm(inf_day_locs,sus_day_locs)/1000
    sus_day_ids$dis_risk <- colSums(dist_i_s<disease_radius)
    #sus_day_ids$dis_risk[which(sus_day_ids$dis_risk > 1)] <- 1
    
    sus_day_ids$dis_y_n <- rbinom(n = nrow(sus_day_ids),size = sus_day_ids$dis_risk,prob = beta)
    sus_day_ids$dis_y_n[which(sus_day_ids$dis_y_n>1)] <- 1

    s_to_i <- sus_day_ids$ID[which(sus_day_ids$dis_y_n == 1)]
    
    day_data$inf[which(day_data$ID %in% s_to_i)] <- "I"
    day_data$timer[which(day_data$ID %in% s_to_i)] <- inf_time
  }
  
  day_data$timer <- day_data$timer - 1
  
  #Move to recovered, set new timer
  i_to_r <- which(day_data$timer <= 0 & day_data$inf =="I")
  day_data$timer[i_to_r] <- rec_time+1
  day_data$inf[i_to_r] <- "R"
  
  r_to_s <- which(day_data$timer <= 0 & day_data$inf =="R")
  day_data$inf[r_to_s] <- "S"
  day_data$timer[r_to_s] <- NA
  
  #Update status for next day
  locs_rec_p[next_day_indices,c("inf","timer")] <- day_data[,c("inf","timer")]
}

sir_df <- locs_rec_p %>% group_by(sim_days) %>% summarise(sus_num = sum(inf=="S"),
                                                          inf_num = sum(inf=="I"),
                                                          rec_num = sum(inf=="R"))
par(mfrow=c(2,1))
par(las=1)
par(mar=(c(5,4,4,4)))
plot(sir_df$sim_days,100*sir_df$sus_num/2000,type="l",col=alpha("forestgreen",.9),
     lwd=5,ylim=c(0,100),ylab="",xlab="",xaxt="n",cex.axis=1.75)
points(sir_df$sim_days,100*sir_df$rec_num/2000,type="l",col=alpha("dodgerblue4",.9),lwd=5)
points(sir_df$sim_days,100*sir_df$inf_num/2000,type="l",col=alpha("firebrick4",.9),lwd=5)
legend("topright",inset=.02,c("Susceptible","Infected","Recovered"),col=c("forestgreen","firebrick4","dodgerblue4"),lty=1,lwd=5,cex=1)
#Rate of change, infections
perc_change_inf <- (sir_df$inf_num[2:length(sir_df$inf_num)] - sir_df$inf_num[1:(length(sir_df$inf_num)-1)])/num_pulls
perc_change_inf_roll <- rollmean(perc_change_inf,20,na.pad=TRUE,align="center")
plot(sir_df$sim_days,100*migr_tracker/2000,type="l",col=alpha("black",.8),lwd=5,ylab="",xlab="",xaxt="n",cex.axis=1.65,ylim=c(0,100))
abline(h=50,lty=2,col="grey80",lwd=3)
par(new=TRUE)
plot(sir_df$sim_days[2:length(sir_df$inf_num)],perc_change_inf_roll*100,col=alpha("firebrick4",.8),lwd=5,type="l",
     axes = FALSE, xlab = "", ylab = "",cex.axis=1.65)
axis(side = 4, at = pretty(range(perc_change_inf*100)),cex.axis=1.65)      # Add second axis
mtext("", side = 4, line = 3)             # Add second axis label
axis(1,                         # Define x-axis manually
     at = seq(30,700,by=60),
     labels = c("Aug","Oct","Dec","Feb","Apr","Jun","Aug","Oct","Dec","Feb","Apr","Jun"),las=2,cex.axis=1.65)
legend("topright",inset=.02,c("% Migrating","% Change, Inf."),col=c("black","firebrick4"),lty=1,lwd=5,cex=1)

# Get sites of infections

inf_sites <- filter(locs_rec_p,inf=="I" & timer == inf_time-1)

#Randomly sample 10% of infections
inf_sites <- inf_sites[sample(1:nrow(inf_sites),nrow(inf_sites)/10),]

countries <- map_data("world")
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
p <- ggplot() + coord_map("mollweide",xlim=c(-125,-60),ylim=c(18,60)) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=country_fill_color, color=country_border_color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(data=breed_loc, aes(x=breeding_lon,y=breeding_lat),colour="firebrick3",alpha=.15,size=4) + 
  geom_point(data=nbreed_loc, aes(x=nonbreeding_lon,y=nonbreeding_lat),colour="dodgerblue3",alpha=.15,size=4) + 
  geom_point(data=inf_sites,aes(x=Lon,y=Lat),colour="black",alpha=.7,shape=18,size=6.5) +
  labs(colour = c("black","firebrick3","dodgerblue3"),shape=c(19,19,8)) + 
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.background = element_rect(fill = alpha("white",.5)))
p

### Simulate pathogen spread Grid Version ####

#Create SIR function #adapted from "Primer of Ecology using R"
density_SIR <- function(yy, beta,lambda, N) {
  S <- yy[1]
  I <- yy[2]
  R <- yy[3]
  
  dS.dt <- -beta*I*S
  dI.dt <- beta*I*S - lambda*I
  dR.dt <- lambda*I
  
  prop_ind_inf <- -dS.dt/S
  
  #If probability of infection, set to 1
  if (prop_ind_inf>1){
    prop_ind_inf <- 1
  }
  #Return probability of infection
  return(prop_ind_inf)
}

N <- num_pulls
inf_time <- 14
rec_time <- 90
lambda <- 1/inf_time
rnaught <- 20
beta <- (rnaught*lambda)/N
start_inf_num <- 50

locs_rec_p <- locs_rec

#Add column to denote simulation time
locs_rec_p$sim_days <- rep(1:num_days, each=num_pulls)

#Add infection status to new column
locs_rec_p$inf <- "S"

#Add infection timer to new column
locs_rec_p$timer <- NA

#Infect starting birds
locs_rec_p$inf[1:start_inf_num] <- "I"
locs_rec_p$timer[1:start_inf_num] <- inf_time

#Get first day
start_day <- locs_rec_p$Day[1]
end_day <- locs_rec_p$Day[nrow(locs_rec_p)]

day_list <- unique(locs_rec_p$sim_days)

#Create grid to evaluate SIR model, big grids here
sir_grid <- dgconstruct(spacing=300, metric=TRUE, resround='down',show_info = FALSE)
#Get where they are in grid system
locs_rec_p$cell <- dgGEO_to_SEQNUM(sir_grid,locs_rec_p$Lon,locs_rec_p$Lat)$seqnum

# For each day in simulation
index_count <- 0
for (each_day in day_list){
  print(each_day)
  
  #Get where next day data exists
  next_day_indices <- which(locs_rec_p$sim_days==(each_day+1))
  
  #Extract df of birds for a given day
  day_data <- filter(locs_rec_p, sim_days==each_day)
  
  #day_data$timer <- as.numeric(day_data$timer)
  
  #For each cell, calculate individual infection probability
  inf_by_cl <- day_data %>% group_by(cell) %>% summarise(sus_cl = sum(inf=="S"),inf_cl = sum(inf=="I"),rec_cl = sum(inf=="R"))
  inf_by_cl <- as.data.frame(inf_by_cl)
  #For each cell
  for (nn in 1:nrow(inf_by_cl)){
    
    #If there are infected AND sus. birds in the cell, otherwise move on
    if (inf_by_cl[nn,3] != 0 & inf_by_cl[nn,2] != 0 ){
      
      cl_target <- inf_by_cl[nn,1]
      inf_by_cl[nn,2:4]
      #Query probability of infection
      prob_in_inf <- density_SIR(inf_by_cl[nn,2:4],beta,lambda,N)
      
      #Get location of S birds
      s_cl_locs <- which(day_data$cell==cl_target & day_data$inf == "S")
      new_status <- sample(c("S","I"),length(s_cl_locs),replace = TRUE,prob = c(1-prob_in_inf,prob_in_inf))
      
      #Update status to infected, add counter
      day_data$inf[s_cl_locs] <- new_status
      day_data$timer[s_cl_locs] <- as.integer(as.logical(new_status== "I"))*inf_time
    }
  }
  
  day_data$timer <- day_data$timer -1
  
  #Move to recovered, set new timer
  i_to_r <- which(day_data$timer <= 0 & day_data$inf =="I")
  day_data$timer[i_to_r] <- rec_time+1
  day_data$inf[i_to_r] <- "R"
  
  r_to_s <- which(day_data$timer <= 0 & day_data$inf =="R")
  day_data$inf[r_to_s] <- "S"
  day_data$timer[r_to_s] <- NA
  
  #Update status for next day
  locs_rec_p[next_day_indices,c("inf","timer")] <- day_data[,c("inf","timer")]
}

sir_df <- locs_rec_p %>% group_by(sim_days) %>% summarise(sus_num = sum(inf=="S"),
                                                          inf_num = sum(inf=="I"),
                                                          rec_num = sum(inf=="R"))
plot(sir_df$sim_days,sir_df$sus_num,type="l",col="forestgreen",lwd=2,ylim=c(0,5000))
points(sir_df$sim_days,sir_df$inf_num,type="l",col="firebrick4",lwd=2)
points(sir_df$sim_days,sir_df$rec_num,type="l",col="dodgerblue4",lwd=2)

plot(sir_df$sim_days,sir_df$inf_num,type="l",col="firebrick4",lwd=2)

fall_mig_peak <- c(80,365+80,365*2+80,365*3+80,365*4+80)
spring_mig_peak <- fall_mig_peak - 180
abline(v=fall_mig_peak)
abline(v=spring_mig_peak)
table(locs_rec_p$sim_days,locs_rec_p$inf,type=2)


cl_sir <- locs_rec_p %>% group_by(sim_days,cell) %>% summarise(ss = sum(inf=="S"),
                                                     ii = sum(inf=="I"),
                                                     rr = sum(inf=="R"))

### Plotting and animation #### 

countries <- map_data("world")
countries <- filter(countries,region %in% c("Canada","USA","Mexico","Guatemala","Bahamas",
                                            "Belize","Honduras","El Salvador","Nicaragua",
                                            "Costa Rica","Panama","Cuba","Jamaica","Haiti",
                                            "Dominican Republic","Puerto Rico","Colombia",
                                            "Venezuela","Brazil","Guyana","Suriname","French Guiana"))
#p <- ggplot() + coord_map("ortho", xlim=c(-145,-80),ylim=c(20,75)) +
#  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="lightgrey", color="darkgrey")
#p
pdf("sir_test_plts.pdf")
for (day_target in seq(1,max(cl_sir$sim_days),by=14)){
  cl_sir_day <- filter(cl_sir,sim_days==day_target)
  grid  <- dgcellstogrid(sir_grid,cl_sir_day$cell,frame=TRUE,wrapcells=TRUE)
  
  cl_sir_day_dgg <- merge(grid,cl_sir_day,by="cell")
  p <- ggplot() + coord_map("mollweide",xlim=c(-125,-70),ylim=c(15,65)) +
    geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="lightgrey", color="darkgrey") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_point(data=breed_loc, aes(x=Lon,y=Lat)) + 
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),axis.text.y = element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank())  +
    scale_fill_gradient2(low="dodgerblue4", mid="grey",high="firebrick4",midpoint=0,limits=c(0,1))+ theme_bw()
  print(p)
}
dev.off()

#Plot with points
#Location of birds on a single day
loc_day <- filter(locs_rec_p,Day == 300)
mappoints <- p + 
  geom_point(data = loc_day,aes(x=Lon,y=Lat),color=alpha(loc_day$col,.7),size=2,pch=19)
plot(mappoints)

#Add column to adjust days
locs_rec_p$adj_dates <- floor((as.numeric(rownames(locs_rec_p))-1)/2000)+1

#Animate
days_to_inc <- seq(from=1,to=365, by=3)
loc_day <- filter(locs_rec_p, Day %in% days_to_inc)
pp <- p +
  geom_point(data = loc_day,aes(x=Lon,y=Lat,group=ID),color=alpha(loc_day$col,.7),size=2,pch=19)
#plot(p)                  

anim <- pp + 
  transition_states(adj_dates,transition_length = 1,state_length = 1) # + labs(title = 'Year: {frame_time}') #+ 
#ggtitle('Julian Day: {closest_state}',
#subtitle = 'Frame {frame} of {nframes}')

#This will display the plot in the viewer tab of rstudio

#the fps should put the gif in the ~13 second range, similar to eBird S&T, but
# for exact matching a third-party software will need to be used.
animate(anim, nframes = 2*length(days_to_inc),fps=20,height=800,width=800)



