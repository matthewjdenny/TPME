rm(list = ls())

setwd("~/Dropbox/PINLab/Projects/R_Code/TPMNE/")
#load current results
#load("Current_Itteration_Results.Rdata")
load("Current_Itteration_Results_McDowell_3-4-14.Rdata")
Itterations <- 1000
Actors <- 17
Latent_Spaces <- 2
Topics <- 50

Metropolis_Results <- Return_List[[1]]
sum(unlist(Metropolis_Results[3001:4000]))

#1.=========== take a look at metropolis iterations ===========#
#function to plot intercept over time
plot_intercepts <- function(intercept){
    intercepts <- rep(0,Itterations)
    for(i in 1:Itterations){
        intercepts[i] <- Metropolis_Results[[Itterations+i]][intercept]
    }
    plot(intercepts, main = paste("Topic:",intercept),ylab= "Value")
}
    
#generate pdf of intercepts   
pdf(file='~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/Intercepts.pdf', width = 25, height = 12)
par(mfrow= c(5,10))
sapply(1:Topics,plot_intercepts)
dev.off()

#function to plot betas over time
plot_betas <- function(topic){
    betas <- matrix(0,ncol = 4, nrow = Itterations)
    for(i in 1:Itterations){
        betas[i,] <- Metropolis_Results[[2*Itterations+i]][topic,]
    }
    matplot(betas, main = paste("Topic:",topic),ylab= "Value")
}

#generate pdf of intercepts   
pdf(file='~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/Betas.pdf', width = 25, height = 12)
par(mfrow= c(5,10))
sapply(1:Topics,plot_betas)
dev.off()

plot_LS_Positions <- function(topic){
    LSP <- matrix(0,nrow=Itterations, ncol= Latent_Spaces)
    for(i in 1:Itterations){
        LSP[i,] <- Metropolis_Results[[i]][,topic,2]
    }
    corel <- cor(LSP[,1], LSP[,2])
    matplot(LSP, main = paste("Topic:",topic, "\n LS Correlations:", round(corel,3)),ylab= "Value")
}
#generate pdf of intercepts   
pdf(file='~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/LS_Positions.pdf', width = 25, height = 12)
par(mfrow= c(5,10))
sapply(1:Topics,plot_LS_Positions)
dev.off()


par(mfrow= c(1,1))

