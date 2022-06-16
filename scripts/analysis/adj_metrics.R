#Script to get derived parameters - mean speed km/day and recovery rate, km/day
library(truncnorm)

#Read in final parameter set
param_set <- readRDS("data/param_sets/WOTH_M1/3_28_22/WOTH_M1_S2.rds")
head(param_set)

hist(param_set$speed_mean_s)
mean_speed_s <- etruncnorm(a=0,b=param_set$max_energy_s,mean = param_set$speed_mean_s,sd = param_set$speed_sd_s)
hist(mean_speed_s)
summary(mean_speed_s)

hist(param_set$speed_mean_f)
mean_speed_f <- etruncnorm(a=0,b=param_set$max_energy_f,mean = param_set$speed_mean_f,sd = param_set$speed_sd_f)
hist(mean_speed_f)
summary(mean_speed_f)

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
spring_speeds <- rowSums(speed_theor_weighted*prob_theor)/rowSums(prob_theor)
#summary(rowSums(speed_theor_weighted*prob_theor)/rowSums(prob_theor))

spring_speed_adj <- summary(rowSums(speed_theor_weighted*prob_theor)/rowSums(prob_theor))

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
#fall_speeds <- rowSums(speed_theor_weighted*prob_theor)/rowSums(prob_theor)
 
fall_speed_adj <- summary(rowSums(speed_theor_weighted*prob_theor)/rowSums(prob_theor))

#Recovery rate
rec_rate_f <- (param_set$max_energy_f*param_set$recovery_rate_f)
rec_rate_s <- (param_set$max_energy_s*param_set$recovery_rate_s)

summary(rec_rate_f)
summary(rec_rate_s)
