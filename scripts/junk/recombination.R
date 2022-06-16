# Script to comapre two parameter sets for a single species - 
# This script is designed to explore the relationships between different posteriors
# under different prior conditions.

library(dplyr)
library(ggplot2)
# Read in training information into single dataframe, then process and analyze
n_files <- 1000
n_in_each <- 100

training_files <- lapply(Sys.glob("data/training_runs/CCSP_3_1_22/CCSP*.csv"), read.csv)

compiled_df <- do.call(bind_rows,training_files)

# Do to some simulations failing, the number of actual simulations likely exceeds
# the number wanted/required. Here, this code truncates to the first simulations
compiled_df <- compiled_df[1:(n_files*n_in_each),]

#Calculate error terms for spring, fall
#Split and sum error for 2 sections of the year, weeks 1:26, and 27:52
err_all <-  compiled_df[1:nrow(compiled_df),29:80]
err_1_26 <- rowSums(err_all[,1:26])
err_27_52 <- rowSums(err_all[,27:52])

#add to df, cut weekly error
compiled_df <- compiled_df[,1:28]
compiled_df$err_spring <- err_1_26
compiled_df$err_fall <- err_27_52

# At this point, the model is doing a better job under better-parameterized conditions, 
# which makes sense.

#The next step is to analyze the best fit models under both training parameter sets.
# The goal is to what degree the density estimates overlap with one another.

#Get 100 best models from both
#top100_f <- filter(compiled_df,quantile(compiled_df$err_fall,.01)*1.05 > err_fall)
top100_f <- filter(compiled_df,quantile(compiled_df$err_fall,.01) > err_fall)

top100_s <- filter(compiled_df,quantile(compiled_df$err_spring,.01) > err_spring)

# Run the recombination step - take fall and spring model parameter sets, and combine
# them to create a list of simulation parameter sets. Because there are a large number 
# of potential combinations, this process can be done randomly
fall_param_inx <- c(4,5,8,9,11,14,15,17,19,26)
spring_param_inx <- c(2,3,6,7,10,12,13,16,18,25)
shared_inx <- c(21,22,23,24)
colnames(top100_f)[fall_param_inx]
colnames(top100_s)[spring_param_inx]
colnames(top100_s)[shared_inx]

# Decide how many simulation sets are required, pull random fall and spring
# parameter sets. This could be rewritten to equally sample all of the parameter
# sets, but if more than 1000 are required, this is likely the best way to do this 
n_param_sets <- 2000
param_sets <- as.data.frame(matrix(NA, ncol = 24, nrow = n_param_sets))
for (n in 1:n_param_sets){
  #Grab a random set of fall parameters
  fll_rand <- (sample(1:nrow(top100_f),1))
  fall_set <- top100_f[fll_rand,fall_param_inx]
  #Grab a random set of spring parameters
  spr_rand <- (sample(1:nrow(top100_s),1))
  spring_set <- top100_s[spr_rand,spring_param_inx]
  #grab a random shared set - coin flip between the two param sets chosen above
  shr_rand <- sample(1:2,1)
  if(shr_rand == 1){
    shared_set <- top100_s[spr_rand,shared_inx]
  } else {
    shared_set <- top100_f[spr_rand,shared_inx]
  }
  new_param_set <- as.data.frame(c(fall_set,spring_set,shared_set))
  param_sets[n,] <- new_param_set
}
colnames(param_sets) <- colnames(as.data.frame(c(fall_set,spring_set,shared_set)))
saveRDS(param_sets,"data/fit_param_sets/VATH_10_29.rds")
# Look at parameter correlations
#pairs(param_sets[c(1,2,4,6,7,8,9)])
#pairs(param_sets[c(10,11,12,13,15,16,17)])
