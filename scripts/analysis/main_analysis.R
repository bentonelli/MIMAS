#Script to create output files for a training run. All output goes to designated 
#output folder (created automatically) unless otherwise specified
spec_abr <- "VATH"
run_date <- "9_20_22"

#Get functions
source(file = "scripts/sim_functions/sim_funct.R")

#Recombine step (Step 1)
#recombine(spec_abr,run_date)

#Save best simulations
best_sims(spec_abr,run_date)

#See how much parameters converged
param_converge(spec_abr,run_date)

#See how correlated parameters are
param_corr(spec_abr,run_date)
