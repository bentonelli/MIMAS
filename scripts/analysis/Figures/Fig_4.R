# Script to demonstrate ABC process
library(truncnorm)
library(ggplot2)

set.seed(123)
num_trials <- 50
alpha_trials <- .5
density_thick <- 20
best_guess <- 200
prior_range <- rtruncnorm(10000,a=0,mean=200,sd=75)
prior_range <- rtruncnorm(10000,a=0,mean=200,sd=75)

par(mar=c(6,6,6,6))
par(pty="s")
par(las=1)
#hist(prior_range)
#jpeg(filename=paste("prior.jpg",sep=""),width = 375,height = 375)
plot(NULL,xlim=c(0,475),ylim=c(0,.01),ylab="",yaxt='n',xlab="Flight Distance, km",cex.lab=1.5,cex.axis = 1.5)# + 
title(ylab="Density", line=.5,cex.lab=1.5) 
#abline(v=best_guess,lty=2,lwd=2) 
points(density(prior_range,density_thick),type="l",lwd=4)
legend("topright",legend=c("Prior"),lty=1,col="black",lwd=4,cex=1.5)
p1 <- recordPlot()

rec_keep <- c()
col_list <- c()
keep_sims <- c()
for (n in 1:num_trials){
  r_sim <- sample(1:10000,1)
  rand_chance <- runif(1)
  if (rand_chance < .1){
    col_new <- "forestgreen"
    keep_sims <- c(keep_sims,prior_range[r_sim])
  } else if (rand_chance > .7) {
    col_new <- "orchid4"
  } else if (prior_range[r_sim] > 110 & prior_range[r_sim] < 210){
    col_new <- "forestgreen"
    keep_sims <- c(keep_sims,prior_range[r_sim])
  } else {
    col_new <- "orchid4"
  }
  rec_keep <- c(rec_keep,prior_range[r_sim])
  col_list <- c(col_list,col_new)
  
  #jpeg(filename=paste("plot_",n,".jpg",sep=""),width = 400,height = 400)
  if (n == num_trials){
    p2 <- plot(NULL,xlim=c(0,475),ylim=c(0,.01),ylab="",yaxt='n',xlab="Flight Distance, km",cex.lab=1.5,cex.axis = 1.5)
    title(ylab="Density",line=.5, cex.lab=1.5) 
    #abline(v=best_guess,lty=2,lwd=2)
    points(density(prior_range,density_thick),type="l",lwd=4,ylab="") 
    abline(v=rec_keep,col=alpha(col_list,.6),lwd=2,lty=1) 
    legend("topright",legend=c("Prior","Accept","Reject"),lty=c(1,1,1),
           col=c("black","forestgreen","orchid4"),lwd=c(4,2,2),cex=1.5)
    p2 <- recordPlot()
  }
  
  #dev.off()
  #Sys.sleep(1)
}
#jpeg(filename=paste("best_sims.jpg",sep=""),width = 400,height = 400)

#plot(NULL,xlim=c(0,475),ylim=c(0,.01),ylab="",yaxt='n',xlab="Flight Distance, km",cex.lab=1.5,cex.axis = 1.5)
#title(ylab="Likelihood", line=0.5, cex.lab=1.5)
#points(density(prior_range),type="l",lwd=5,ylab="")
#abline(v=keep_sims,col=alpha("forestgreen",.5),lwd=4,lty=2)

#dev.off()

#jpeg(filename=paste("poster.jpg",sep=""),width = 400,height = 400)

plot(NULL,xlim=c(0,475),ylim=c(0,.01),ylab="",yaxt='n',xlab="Flight Distance, km",cex.lab=1.5,cex.axis = 1.5)
title(ylab="Density", line=.5,cex.lab=1.5)
points(density(prior_range,density_thick),type="l",lwd=4,ylab="")
points(density(keep_sims,density_thick),type="l",lwd=4,ylab="",col="forestgreen")
legend("topright",legend=c("Prior","Posterior"),lty=c(1,1),
       col=c("black","forestgreen"),lwd=c(4,4),cex=1.5)
#abline(v=keep_sims,col=alpha("forestgreen",.5),lwd=2,lty=1)
p3 <- recordPlot()
#dev.off()


