rm(list = ls())

setwd("~/Dropbox/PINLab/Projects/R_Code/TPMNE/")
#load current results
load("Current_Itteration_Results.Rdata")

#1. take a look at intercepts over metropolis iterations
Metropolis_Results <- Return_List[[1]]

#function to plot intercept over time
plot_intercepts <- function(intercept){
    intercepts <- rep(0,1000)
    for(i in 1:1000){
        intercepts[i] <- Metropolis_Results[[1000+i]][intercept]
    }
    plot(intercepts, main = paste("Topic:",intercept),ylab= "Value")
}
    
#generate pdf of intercepts   
pdf(file='~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/Intercepts.pdf', width = 25, height = 15)
par(mfrow= c(5,10))
sapply(1:50,plot_intercepts)
dev.off()

par(mfrow= c(1,1))


#2. take a look at betas over time
#function to plot intercept over time
plot_betas <- function(topic){
    betas <- matrix(0,ncol = 4, nrow = 1000)
    for(i in 1:1000){
        betas[i,] <- Metropolis_Results[[2000+i]][topic,]
    }
    matplot(betas, main = paste("Topic:",topic),ylab= "Value")
}

#generate pdf of intercepts   
pdf(file='~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/Betas.pdf', width = 25, height = 15)
par(mfrow= c(5,10))
sapply(1:50,plot_betas)
dev.off()

#black = 
par(mfrow= c(1,1))




