#'This script takes in the model simulation runs from the initial parameter space
#'exploration, performs rejection-ABC for each season - fall and spring - 
#'and then resamples from the posteriors to get new parameter sets. These new 
#'parameter sets can then be used to perform another cycle of rejection ABC to 
#'arrive at better estimates for a given species. The output of this script
#'is a list of top models for fall, spring and both combined

library(dplyr)
library(ggplot2)
library(lubridate)
# Read in training information into single dataframe, then process and analyze
n_files <- 1000
n_in_each <- 100

#Set location to save file
output_loc <- "data/param_sets/BUOR/8_21_22/BUOR_S1.rds"

#Set location of simulation files
training_files <- lapply(Sys.glob("data/training_runs/BRSP_8_21_22/S1/BRSP*.csv"), read.csv)

compiled_df <- do.call(bind_rows,training_files)

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
# in that bunch. Adding back more models is neccesary in order to
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

saveRDS(three_batch_out,output_loc)
