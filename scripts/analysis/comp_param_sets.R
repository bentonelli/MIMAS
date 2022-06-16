# Script to comapre two parameter sets for a single species - 
# This script is designed to explore the relationships between different posteriors
# under different prior conditions.

library(dplyr)
library(ggplot2)
# Read in training information into single dataframe, then process and analyze
n_files <- 1000
n_in_each <- 100


training_files_set1 <- lapply(Sys.glob("data/training_runs/Step_1/Varied_Thrush/VATH_S1_3_11_22/VATH_*.csv"), read.csv)
training_files_set2 <- lapply(Sys.glob("data/training_runs/VATH_10_27_WI/VATH_WI*.csv"), read.csv)

compiled_df_set1 <- do.call(bind_rows,training_files_set1)
compiled_df_set2 <- do.call(bind_rows,training_files_set2)

# Do to some simulations failing, the number of actual simulations likely exceeds
# the number wanted/required. Here, this code truncates to the first simulations
compiled_df_set1 <- compiled_df_set1[1:(n_files*n_in_each),]
compiled_df_set2 <- compiled_df_set2[1:(n_files*n_in_each),]

#Calculate error terms for spring, fall
#Split and sum error for 2 sections of the year, weeks 1:26, and 27:52
err_all <-  compiled_df_set1[1:nrow(compiled_df_set1),29:80]
err_1_26 <- rowSums(err_all[,1:26])
err_27_52 <- rowSums(err_all[,27:52])

#add to df, cut weekly error
compiled_df_set1 <- compiled_df_set1[,1:28]
compiled_df_set1$err_spring <- err_1_26
compiled_df_set1$err_fall <- err_27_52

#Repeat
err_all <-  compiled_df_set2[1:nrow(compiled_df_set2),29:80]
err_1_26 <- rowSums(err_all[,1:26])
err_27_52 <- rowSums(err_all[,27:52])

#add to df, cut weekly error
compiled_df_set2 <- compiled_df_set2[,1:28]
compiled_df_set2$err_spring <- err_1_26
compiled_df_set2$err_fall <- err_27_52

#Take a quick look at the distribution of the error terms
plot(density(compiled_df_set2$sum_err,bw=10),col="blue",main="Full Model Error",ylab="Density",xlab="Sum Error",xlim=c(150,1000))
points(density(compiled_df_set1$sum_err,bw=10),col="red",type="l")
legend("topright",legend=c("VWI","WI"),col=c("red","blue"),lty=1)

#Next, fall model error
plot(density(compiled_df_set2$err_fall,bw=5),col="blue",main="Fall Model Error",ylab="Density",xlab="Sum Error",xlim=c(50,600))
points(density(compiled_df_set1$err_fall,bw=5),col="red",type="l")
legend("topright",legend=c("VWI","WI"),col=c("red","blue"),lty=1)

#It should be that the model that is better parameterized has a lower error in general
plot(density(compiled_df_set2$err_spring,bw=5),col="blue",main="Spring Model Error",ylab="Density",xlab="Sum Error",xlim=c(50,600))
points(density(compiled_df_set1$err_spring,bw=5),col="red",type="l")
legend("topright",legend=c("VWI","WI"),col=c("red","blue"),lty=1)
# At this point, the model is doing a better job under better-parameterizd conditions, 
# which makes sense.

#The next step is to analyze the best fit models under both training parameter sets.
# The goal is to what degree the density estimates overlap with one another.

#Get 100 best models from both
top100_s1_f <- filter(compiled_df_set1,quantile(compiled_df_set1$err_fall,.01) > err_fall)
top100_s2_f <- filter(compiled_df_set2,quantile(compiled_df_set2$err_fall,.01) > err_fall)

plot(density(top100_s2_f$err_fall),col="blue",main="Fall Model Error - 1%",ylab="Density",xlab="Sum Error",xlim=c(105,150))
points(density(top100_s1_f$err_fall),type="l",col="red")
legend("topright",legend=c("VWI","WI"),col=c("red","blue"),lty=1)
top100_s1_s <- filter(compiled_df_set1,(quantile(compiled_df_set1$err_spring,.01)*1.10) > err_spring)
top100_s2_s <- filter(compiled_df_set2,(quantile(compiled_df_set2$err_spring,.01)*1.10) > err_spring)


plot(density(top100_s2_s$err_spring),col="blue",main="Spring Model Error - 1%",ylab="Density",xlab="Sum Error",xlim=c(90,170))
points(density(top100_s1_s$err_spring),type="l",col="red")
legend("topright",legend=c("VWI","WI"),col=c("red","blue"),lty=1)

#Test with a single parameter, fall migration speed
plot(density(top100_s2_f$start_date_u_f,bw=2),col="blue",xlim=c(220,320),main="Fall Start Date")
points(density(top100_s1_f$start_date_u_f,bw=2),col="red",type="l")

#Show the priors
points(density(compiled_df_set2$start_date_u_f,bw=2),col="blue",type="l",lty=2)
points(density(compiled_df_set1$start_date_u_f,bw=2),col="red",type="l",lty=2)
legend("topright",legend=c("WD Post.","WD Prior","VWD Post.","VWD Prior"),lty=c(1,2,1,2),col=c("blue","blue","red","red"))

#Do the above for each parameter for fall 
fall_param_inx <- c(4,5,8,9,15,17,19,23,24,26)
spring_param_inx <- c(2,3,6,7,13,16,18,23,24,25)
pdf(file="param_set_comparison.pdf")
for (n in fall_param_inx){

  plot(density(top100_s2_f[,n]),col="blue",main=colnames(top100_s1_f[n]))
  points(density(top100_s1_f[,n]),col="red",type="l")
  
  #Show the priors
  points(density(compiled_df_set2[,n]),col="blue",type="l",lty=2)
  points(density(compiled_df_set1[,n]),col="red",type="l",lty=2)
  legend("topright",legend=c("WI Post.","WI Prior","VWI Post.","VWI Prior"),lty=c(1,2,1,2),col=c("blue","blue","red","red"))
}

#Do the above for each parameter for spring models
for (n in spring_param_inx){

  plot(density(top100_s2_s[,n]),col="blue",main=colnames(top100_s1_s[n]))
  points(density(top100_s1_s[,n]),col="red",type="l")
  
  #Show the priors
  points(density(compiled_df_set2[,n]),col="blue",type="l",lty=2)
  points(density(compiled_df_set1[,n]),col="red",type="l",lty=2)
  legend("topright",legend=c("WI Post.","WI Prior","VWI Post.","VWI Prior"),lty=c(1,2,1,2),col=c("blue","blue","red","red"))
}
dev.off()






