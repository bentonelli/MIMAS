library(dplyr)
library(ggplot2)
library(lubridate)
# Read in training information into single dataframe, then process and analyze
total_sims <- 100000
training_files <- lapply(Sys.glob("data/training_runs/Step_1/Varied_Thrush/VATH_S1_3_11_22/VATH*.csv"), read.csv)
#training_files[[1]]

compiled_df <- do.call(bind_rows,training_files)
compiled_df <- compiled_df[1:total_sims,]

#Get start date
m_start_week <- week(as_date(compiled_df$start_date[1]-1))
all_weeks <- seq(1,52,by=1)
fall_weeks <- (seq(m_start_week,m_start_week+25,by=1)-1)%%52 + 1
spring_weeks <- all_weeks[!all_weeks %in% fall_weeks]

#Pull first 100 and plot weekly error
plot(NULL,xlim=c(1,52),ylim=c(0,40))
for (n in 1:100){
  err_in <- as.numeric(compiled_df[n,29:80])
  points(1:52,err_in,type="l",col=alpha("blue",.2),lwd=.2)
}

#Split and sum error for 2 sections of the year, based on simulation start date
err_all <-  compiled_df[1:nrow(compiled_df),29:80]
err_s <- rowSums(err_all[,spring_weeks])
err_f <- rowSums(err_all[,fall_weeks])

#Plot against each other
plot(err_s,err_f,cex=.1,col=alpha("dodgerblue4",.2))

#add to df, cut weekly error
compiled_df <- compiled_df[,1:28]
compiled_df$err_spring <- err_s
compiled_df$err_fall <- err_f

# Get 1% best models from both + models with error within 5% of the 'worst' model
# in that bunch. Adding back more models is neccesary in order to
top_f_models <- filter(compiled_df,quantile(compiled_df$err_fall,.01)*1.05 > err_fall)

top_s_models <- filter(compiled_df,quantile(compiled_df$err_spring,.01)*1.05 > err_spring)

top_f_s_models <- rbind(top_f_models,top_s_models)

# Run the recombination step - take fall and spring model parameter sets, and combine
# them to create a list of simulation parameter sets. Because there are a large number 
# of potential combinations, this process can be done randomly
fall_param_inx <- c(4,5,8,9,11,14,15,17,19,26)
spring_param_inx <- c(2,3,6,7,10,12,13,16,18,25)
shared_inx <- c(20,21,22,23,24)
colnames(top_f_models)[fall_param_inx]
colnames(top_s_models)[spring_param_inx]
colnames(top_f_models)[shared_inx]

pdf(file = "VATH_S1_Convergence_3_11_22.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
for (n in 2:26){
  set_range <- density(compiled_df[,n])
  xx_max <- max(set_range$x)
  xx_min <- min(set_range$x)
  if (n %in% fall_param_inx){
    set_bw <- density(top_f_models[,n])
    set_bw <- set_bw$bw
    plot(density(top_f_models[,n]),col="dodgerblue4",lwd=3,xlim=c(xx_min,xx_max),main=colnames(top_f_models)[n])
    points(density(compiled_df[,n],bw=set_bw),col="firebrick3",lwd=3,type="l")
  } else if (n %in% spring_param_inx) {
    set_bw <- density(top_s_models[,n])
    set_bw <- set_bw$bw
    plot(density(top_s_models[,n]),col="dodgerblue4",lwd=3,xlim=c(xx_min,xx_max),main=colnames(top_s_models)[n])
    points(density(compiled_df[,n],bw=set_bw),col="firebrick3",lwd=3,type="l")
  } else if (n %in% shared_inx){
    set_bw <- density(top_f_s_models[,n])
    set_bw <- set_bw$bw
    plot(density(top_f_s_models[,n]),col="dodgerblue4",lwd=3,xlim=c(xx_min,xx_max),main=colnames(top_f_s_models)[n])
    points(density(compiled_df[,n],bw=set_bw),col="firebrick3",lwd=3,type="l")
  }

}
dev.off()


