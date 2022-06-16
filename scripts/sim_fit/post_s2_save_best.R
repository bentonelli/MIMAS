# Simple script to save best models

source(file = "scripts/sim_functions/sim_funct.R")

#Location of saved simulation files
sims_location <- "~/Documents/Coding/R/M-IBM/data/training_runs/Step_2/Wood_Thrush_M2/WOTH*.csv"

#Define buffer for error - here set to 3%
sim_best_models <- best_models(sims_location,.03,100000)

#Save to RDS
saveRDS(sim_best_models,"~/Documents/Coding/R/M-IBM/data/fit_param_sets/WOTH_M2/Step_2/WOTH_M2_3_30_22.rds")
