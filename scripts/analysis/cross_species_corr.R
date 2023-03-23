

spec_list <- c("BRSP","BUOR","CCSP","HOWA","OROR","RNSA","TOWA","VATH","WOTH_M1","YBSA")
run_date <- "9_20_22"

#Get functions
source(file = "scripts/sim_functions/sim_funct.R")

count <- 0
corr_sp_list <- list()
for (each_spec in spec_list){
  count <- count + 1
  #See how correlated parameters are
  corr_sp <- param_corr(each_spec,run_date,return_matrix=TRUE)
  corr_sp_list[[count]] <- corr_sp
}
totals_mat <- corr_sp_list[[1]]
full_corr_vec <- c()
for (nn in 2:10){
  totals_mat <- totals_mat + corr_sp_list[[nn]]
  temp_mat <- corr_sp_list[[nn]]
  diag(temp_mat) <- NA
  full_corr_vec <- c(full_corr_vec,unlist(as.numeric(temp_mat)))
}
summary(full_corr_vec^2)

corrplot(totals_mat/10, method = 'color',order = 'AOE')
diag(totals_mat) <- NA
summary(as.numeric(totals_mat/10),na.rm=TRUE)


### Wood Thrush ####
corr_sp_list[[9]]
