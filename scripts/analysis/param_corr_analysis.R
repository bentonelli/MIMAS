# Script to compare parameter correlations across best fit models
library(corrplot)

#Read in parameter set
fit_params <- readRDS("data/fit_param_sets/WOTH_M1/Step_2/WOTH_M1_3_29_22.rds")

plot(fit_params$speed_mean_s,fit_params$start_date_u_s,cex=.1)
wanted_corr_params <- as.matrix(fit_params[,2:26])
cor_data <- round(cor(wanted_corr_params),2)

#Remove NA rows
cor_data <- cor_data[-as.numeric(which(rowSums(cor_data,na.rm=TRUE)==TRUE)),-as.numeric(which(rowSums(cor_data,na.rm=TRUE)==TRUE))]

#Plot correlations. Note here that some parameters might be not be continuous,
#so interpret with caution!
corrplot(cor_data, method = 'color',order = 'AOE')

#Get summary stats, ploot histogram
summary(cor_data[which(cor_data!=1)])
hist(cor_data[which(cor_data!=1)],xlim=c(-1,1),breaks=5)
