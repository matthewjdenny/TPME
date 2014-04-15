#this script reads in hand codings of topics based on whther they fall into one of 5 categories. Asortativity and properties of email networks are then compared across areas, particularly with regards to congestion.

#0 = other, 1 = broadcast, 2 = general work, 3= social, 4= coordination

#preliminaries
rm(list = ls())
setwd("~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/")
library(foreign)
library(ggplot2)
library(gridExtra)
codings <- read.csv("Email_Area_Coding_MJD_3-25-14.csv")

#load in topic-word-beta data 
load("Columbus_10M_3-13-14Topic_Assortativity_Data.Rdata")
topic_data_col <- cbind(codings[1:50,4],topic_data)
load("New_Hannover_10M_3-14-14Topic_Assortativity_Data.Rdata")
topic_data_nhc <- cbind(codings[1:50,5],topic_data)
load("Transylvania_50M_3-13-14Topic_Assortativity_Data.Rdata")
topic_data_tran <- cbind(codings[1:50,3],topic_data)
load("McDowell_50M_3-13-14Topic_Assortativity_Data.Rdata")
topic_data_mcd <- cbind(codings[,2],topic_data)

data_list <- list(topic_data_col, topic_data_nhc, topic_data_tran, topic_data_mcd)
counties <- c("Columbus", "New Hannover", "Transylvania", "McDowell")
#
i = 4
for(i in 1:4){
    #for coordination topics
    temp <- data_list[[i]]
    coord <- temp[which(temp[,1] == 4),]
    aves_sd_coord <- matrix(0,nrow =4,ncol = 2)
    for(j in 1:4){
        aves_sd_coord[j,1] <- mean(coord[,j+1])
        aves_sd_coord[j,2] <- mean(coord[,j+5])    
    }
    aves_sd_all <- matrix(0,nrow =4,ncol = 2)
    for(j in 1:4){
        aves_sd_all[j,1] <- mean(temp[,j+1])
        aves_sd_all[j,2] <- mean(temp[,j+5])    
    }
    aves_sd <- rbind(aves_sd_coord,aves_sd_all)
    aves_sd <- as.data.frame(aves_sd)
    aves_sd <- cbind(c("M-M", "M-F","F-M", "F-F","M-M (All)", "M-F (All)","F-M (All)", "F-F (All)"),aves_sd)
    aves_sd <- as.data.frame(aves_sd)
    names(aves_sd) <- c("Tie","Parameter_Estimate","SE")
    plot <- ggplot(aves_sd[1:8,], aes(x=Tie, y=Parameter_Estimate))+geom_errorbar(aes(ymin=Parameter_Estimate-SE, ymax=Parameter_Estimate+SE),col = c("black","red", "green", "blue","grey","pink", "lightgreen", "lightblue"), width=.1, size = 1.3)+geom_point(fill= c("black","red", "green", "blue","grey","pink", "lightgreen", "lightblue"),colour="black",pch=c(21,22,23,24,21,22,23,24), size = 4)+ labs(title = paste(counties[i],"County Asortativity Comparison for Coordination vs. All Topics") )
#     
#     #now for broadcast topics
#     temp <- data_list[[i]]
#     coord <- temp[which(temp[,1] == 2),]
#     aves_sd_coord <- matrix(0,nrow =4,ncol = 2)
#     for(j in 1:4){
#         aves_sd_coord[j,1] <- mean(coord[,j+1])
#         aves_sd_coord[j,2] <- mean(coord[,j+5])    
#     }
#     aves_sd_all <- matrix(0,nrow =4,ncol = 2)
#     for(j in 1:4){
#         aves_sd_all[j,1] <- mean(temp[,j+1])
#         aves_sd_all[j,2] <- mean(temp[,j+5])    
#     }
#     aves_sd <- rbind(aves_sd_coord,aves_sd_all)
#     aves_sd <- as.data.frame(aves_sd)
#     aves_sd <- cbind(c("M-M", "M-F","F-M", "F-F","M-M (All)", "M-F (All)","F-M (All)", "F-F (All)"),aves_sd)
#     aves_sd <- as.data.frame(aves_sd)
#     names(aves_sd) <- c("Tie","Parameter_Estimate","SE")
#     plot2 <- ggplot(aves_sd[1:8,], aes(x=Tie, y=Parameter_Estimate))+geom_errorbar(aes(ymin=Parameter_Estimate-SE, ymax=Parameter_Estimate+SE),col = c("black","red", "green", "blue","grey","pink", "lightgreen", "lightblue"), width=.1, size = 1.3)+geom_point(fill= c("black","red", "green", "blue","grey","pink", "lightgreen", "lightblue"),colour="black",pch=c(21,22,23,24,21,22,23,24), size = 4)+ labs(title = paste(counties[i],"County Asortativity Comparison for Broadcast vs. All Topics") )
    
    
    #transform with plogis
    temp <- data_list[[i]]
    coord <- temp[which(temp[,1] == 4),]
    
    aves_sd_coord <- matrix(0,nrow =4,ncol = 2)
    temp2 <- rep(0,length(coord[,1]))
    for(j in 1:4){
        for(m in 1:length(temp2)){
            temp2[m] <- plogis((coord[m,2]+coord[m,j+4] -coord[m,4]) )
        }
        aves_sd_coord[j,1] <- mean(temp2)
        aves_sd_coord[j,2] <- sd(temp2)    
    }
    
    coord2 <- temp[which(temp[,1] == 2),]
    aves_sd_coord2 <- matrix(0,nrow =4,ncol = 2)
    for(j in 1:4){
        temp3 <- plogis(coord2[,2]+coord2[,j+4]-coord2[,4])
        aves_sd_coord2[j,1] <- mean(temp3)
        aves_sd_coord2[j,2] <- sd(temp3)    
    }
    
    
    aves_sd_all <- matrix(0,nrow =4,ncol = 2)
    for(j in 1:4){
        temp4 <- plogis(temp[,2]+temp[,j+4]-temp[,4])
        aves_sd_all[j,1] <- mean(temp4)
        aves_sd_all[j,2] <- sd(temp4)    
    }
    
    aves_sd <- rbind(aves_sd_coord,aves_sd_coord2,aves_sd_all)
    aves_sd <- as.data.frame(aves_sd)
    aves_sd <- cbind(c("C:M-M", "C:M-F","C:F-M", "C:F-F ","B:M-M", "B:M-F","B:F-M", "B:F-F", "A:M-M", "A:M-F","A:F-M", "A:F-F"),aves_sd)
    aves_sd <- as.data.frame(aves_sd)
    names(aves_sd) <- c("Tie","Parameter_Estimate","SE")
    plot <- ggplot(aves_sd[1:12,], aes(x=Tie, y=Parameter_Estimate))+geom_errorbar(aes(ymin=Parameter_Estimate-SE, ymax=Parameter_Estimate+SE),col = c("black","red", "green", "blue","black","red", "green", "blue","grey","pink", "lightgreen", "lightblue"), width=.1, size = 1.3)+geom_point(fill= c("black","red", "green", "blue","black","red", "green", "blue","grey","pink", "lightgreen", "lightblue"),colour="black",pch=c(21,22,23,24,21,22,23,24,21,22,23,24), size = 4)+ labs(title = paste(counties[i],"County Asortativity Comparison for Coordination vs. All Topics") )

    plot

    
    grid.arrange( plot,plot2, ncol=2)
    
}







