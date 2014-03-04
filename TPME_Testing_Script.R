rm(list = ls())


#take a look at intercepts
intercepts <- rep(0,1000)
Metropolis_Results <- Return_List[[1]]

for(i in 1:1000){
    intercepts[i] <- Metropolis_Results[[2000+i]][20]
}
plot(intercepts)
    
    