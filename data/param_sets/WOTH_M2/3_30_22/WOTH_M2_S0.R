#Parameter set based on all data up until 2012

#Define species name
species_target <- "Wood Thrush"

# Defines variables for flight length
#Mean daily flight speed, spring
speed_mean_s <- rtruncnorm(1,a=0,mean=370,sd=200)

#SD of daily flight speed, spring
speed_sd_s <- rtruncnorm(1,a=0,mean=300,sd=300)

#Mean daily flight speed, fall
speed_mean_f <- rtruncnorm(1,a=0,mean=220,sd=200)

#SD of daily flight speed, fall
speed_sd_f <- rtruncnorm(1,a=0,mean=300,sd=300)

#Migration start date programming  - as Julian day

#Spring start date (constrained to Julian days). This varies among individuals
# so all birds don't initiate migration on same date
start_date_u_s <- rtruncnorm(1,a=1,b=365,mean=120,sd=30)

#Spring migration dates standard deviation
start_date_sd_s <- rtruncnorm(1,a=0,mean=30,sd=30)

#Fall migration dates 
start_date_u_f <- rtruncnorm(1,a=1,b=365,mean=290,sd=30)

start_date_sd_f <- rtruncnorm(1,a=0,mean=30,sd=30)

#Orientation parameters - means should be 0
bear_err_mean_s <- rnorm(1,0,10)
bear_err_sd_s <- rtruncnorm(1,a=0,mean=30,sd=20)
bear_err_mean_f <- rnorm(1,0,10)
bear_err_sd_f <- rtruncnorm(1,a=0,mean=30,sd=20)

#Energetics info
max_energy_s <- rtruncnorm(1,a=0,mean=2000,sd=1000)
max_energy_f <- rtruncnorm(1,a=0,mean=2000,sd=1000)
recovery_rate_s <- rtruncnorm(1,a=0,b=1,mean=.15,sd=.1)
recovery_rate_f <- rtruncnorm(1,a=0,b=1,mean=.15,sd=.1)

#Migration time max
max_mig_time_s <- 120
max_mig_time_f <- 120

#Goal radius - how far from "home" does a bird need before it "snaps" to it's location
goal_radius <- 100

#Migratory connectivity - strength
mig_con <- rtruncnorm(1,a=0,b=1,mean=.8,sd=.4)

# Migratory connectivity type
# Type - 1 is leapfrog (breeding lat is neg. correlated with non-breeding lat),
# Type 2 is chain, the opposite, (breeding lat is pos. correlated with non-breeding lat)
mig_con_type <- ceiling(runif(1,0,2))

# Connection between breeding latitude and migration timing. 1 means a strong
# correlation between early departures and higher breeding latitudes.
migr_timing_lat_s <- runif(1,0,1)
migr_timing_lat_f <- runif(1,0,1)

