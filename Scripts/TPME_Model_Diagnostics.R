
Generate_Model_Diagnsotics <- function(input_folder_path = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", input_file = "Test",LS_Actor = 2, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/",Thin_Itterations = 20, vocab = vocabulary,county_name = "McDowell_County",skip_first = 1){
    #load current results
    #load("Current_Itteration_Results.Rdata")
    #for testing load("~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/Transylvania_Sample_10M_2011_3-13-14.Rdata")
    
    #plots to generate
    library(coda)
    #library(statnet)
    library(gregmisc)
    library(ggplot2)
    print("Loading Data...")
    load(paste(input_folder_path,input_file,".Rdata", sep = ""))
    
    print("Extracting Reduced Data")
    #extract current metropolis results
    Metropolis_Results <- Return_List[[1]]
    Topic_Model_Results <- Return_List[[2]]
    #free up memory
    rm(Return_List)
    
    
    
    
    
    
    skip_first= skip_first+1
    #remove the first skip_first itterations of each sublist and recombine
    Itterations <- length(Metropolis_Results)/7
    temp <- append(Metropolis_Results[skip_first:Itterations],
                Metropolis_Results[(Itterations+skip_first):(2*Itterations)])
    temp <- append(temp,Metropolis_Results[(2*Itterations+skip_first):(3*Itterations)])
    temp <- append(temp,Metropolis_Results[(3*Itterations+skip_first):(4*Itterations)])
    temp <- append(temp,Metropolis_Results[(4*Itterations+skip_first):(5*Itterations)])
    temp <- append(temp,Metropolis_Results[(5*Itterations+skip_first):(6*Itterations)])
    temp <- append(temp,Metropolis_Results[(6*Itterations+skip_first):(7*Itterations)])
    Metropolis_Results <- temp
    
    #thin out the data by taking every Thin_Itterations itteration for the metropolis step
    Metropolis_Results <- Metropolis_Results[seq(1, length(Metropolis_Results),Thin_Itterations)]

    #get model information and extract data
    Itterations <- length(Metropolis_Results)/7
    Latent_Spaces <- length(Metropolis_Results[[1]][,1,1])
    Topics <- length(Metropolis_Results[[1]][1,,1])
    Actors <- length(Metropolis_Results[[1]][1,1,])
    Token_Topic_Assignments <- Topic_Model_Results[[1]]
    Topic_Present_Edge_Counts <- Topic_Model_Results[[2]]
    Topic_Absent_Edge_Counts <- Topic_Model_Results[[3]]
    Word_Type_Topic_Counts <- Topic_Model_Results[[4]]
    Edge_Topic_Assignments <- Topic_Model_Results[[5]]
    
    
    #get the total number of tokens assigned to each topic
    Topic_Token_Totals <- apply(Word_Type_Topic_Counts,2,sum)
    #get the total number of present edges assigned to each topic
    Email_Assignments <- apply(Topic_Present_Edge_Counts,3,sum)
    #number of words
    num_words <- length(vocab)
    
    
    #this list will be used to store lists of three vectors: word indicies in descending order of likelihood, word probability in topic, actual word. 
    Topic_Top_Words <- list()
    Top_Ten_Words <- rep("",Topics)
    Top_Four_Words <- rep("",Topics)
    
    
    for(i in 1:Topics){
        indicies <- order(Word_Type_Topic_Counts[,i],decreasing = TRUE)
        probabilities <- Word_Type_Topic_Counts[indicies,i]/Topic_Token_Totals[i]
        #reduce to only words with non-zero probability
        indicies <- indicies[which(probabilities > 0)]
        probabilities <- probabilities[which(probabilities > 0)]
        top_words <- vocab[indicies,]
        Topic_List <- list(top_words,probabilities,indicies)
        Topic_Top_Words <- append(Topic_Top_Words,list(Topic_List))
        #add top words
        words <- ""
        words4 <- ""
        for(j in 1:10){
            words <- paste(words,top_words[j],", ",sep = "")
        }
        for(j in 1:4){
            words4 <- paste(words4,top_words[j],", ",sep = "")
        }
        Top_Ten_Words[i] <- words
        Top_Four_Words[i] <- words4
    }
    
    
    
    
    
    
    
    
    #display current accept rate
    start <- (3*Itterations) + 1
    end <- 4*Itterations
    print(paste("Current_Accept_Rate:", sum(unlist(Metropolis_Results[start:end]))/Itterations))
         
    #======= plotting function definitions =======#
    #function to plot intercept over time
    plot_intercepts <- function(intercept){
        intercepts <- rep(0,Itterations)
        for(i in 1:Itterations){
            intercepts[i] <- Metropolis_Results[[Itterations+i]][intercept]
        }
        plot(intercepts, main = paste("Topic:",intercept,"Geweke:",round(geweke.diag(intercepts)$z,2)),ylab= "Value",pch = 20)
    }
    #function to plot betas over time
    plot_betas <- function(topic){
        betas <- matrix(0,ncol = 4, nrow = Itterations)
        for(i in 1:Itterations){
            betas[i,] <- Metropolis_Results[[2*Itterations+i]][topic,]
        }
        matplot(betas, main = paste("Topic:",topic,"Geweke","\n MM:",round(geweke.diag(betas[,1])$z,2) , "MF:",round(geweke.diag(betas[,2])$z,2) ,"\n FM:",round(geweke.diag(betas[,3])$z,2) ,"FF:",round(geweke.diag(betas[,4])$z,2)),ylab= "Value",pch = 20)
    }
    #function to plot latent spaces for a given actor
    plot_LS_Positions <- function(topic){
        LSP <- matrix(0,nrow=Itterations, ncol= Latent_Spaces)
        for(i in 1:Itterations){
            LSP[i,] <- Metropolis_Results[[i]][,topic,LS_Actor]
        }
        corel <- cor(LSP[,1], LSP[,2])
        matplot(LSP, main = paste("Topic:",topic, "\n LS Correlations:", round(corel,3),"\n Geweke - LS1:",round(geweke.diag(LSP[,1])$z,2),"LS2:",round(geweke.diag(LSP[,2])$z,2)),ylab= "Value",pch = 20)
    }
    plot_Topic_Network <- function(topic){
        slice <- Topic_Present_Edge_Counts[,,topic]
        coordinates <- cbind(Metropolis_Results[[Itterations]][1,topic,],Metropolis_Results[[Itterations]][2,topic,])
        coordinates <- as.data.frame(coordinates)
        plot(coordinates,main = paste("Topic:",topic,"Number of Edges:",sum(slice),"\n", Top_Four_Words[topic]),pch = 20, col = "red")
        #add in lines between actors
        for(i in 1:Actors){
            for(j in 1:Actors){
                if(slice[i,j] >0){
                    lines(c(coordinates[i,1],coordinates[j,1]) , c(coordinates[i,2],coordinates[j,2]), col = "black")
                }
            }
        }
        
    }
    
    #make plots will all top words for one per page output
    plot_Topic_Network_Full <- function(topic){
        slice <- Topic_Present_Edge_Counts[,,topic]
        coordinates <- cbind(Metropolis_Results[[Itterations]][1,topic,],Metropolis_Results[[Itterations]][2,topic,])
        coordinates <- as.data.frame(coordinates)
        plot(coordinates,main = paste("Topic:",topic,"Number of Edges:",sum(slice),"\n", Top_Ten_Words[topic]),pch = 20, col = "red")
        #add in lines between actors
        for(i in 1:Actors){
            for(j in 1:Actors){
                if(slice[i,j] >0){
                    lines(c(coordinates[i,1],coordinates[j,1]) , c(coordinates[i,2],coordinates[j,2]), col = "black")
                }
            }
        }
        
    }
    
    #function to plot log likelihoods over time
    plot_likelihoods <- function(Itterations){
        likelihoods <- matrix(0,ncol = 2, nrow = Itterations)
        for(i in 1:Itterations){
            likelihoods[i,1] <- Metropolis_Results[[4*Itterations+i]]
            likelihoods[i,2] <- Metropolis_Results[[5*Itterations+i]]
        }
        matplot(likelihoods, main = "Log Likelihoods of Current and Proposed Positions over Time",ylab= "Value",pch = 20)
    }
    
    #function to plot log ratios vs log uniform draws over time
    plot_ratio_lud <- function(Itterations){
        #likelihoods <- matrix(0,ncol = 2, nrow = Itterations)
        accept <- rep(0,Itterations)
        likelihoods <- rep(0,Itterations)
        lud <- rep(0,Itterations)
        colors <- rep("",Itterations)
        for(i in 1:Itterations){
            likelihoods[i] <- Metropolis_Results[[4*Itterations+i]] - Metropolis_Results[[5*Itterations+i]] #- Metropolis_Results[[6*Itterations+i]]
            accept[i] <- Metropolis_Results[[3*Itterations+i]]
            lud[i] <- Metropolis_Results[[6*Itterations+i]]
            if(Metropolis_Results[[3*Itterations+i]] == 1){
                col <- "blue"
            }else{
                col <- "red"
            }
            colors[i] <- col
            
            accepted <- accept*likelihoods
            accepted <- accepted[which(accept == 1)]
            
        }
        #fit <- lm(accepted~1:10000)
        print(paste("Average Accepted Prob:",mean(accepted)))
        len <- length(accepted)
        print(paste("Average Accepted Prob last 10 percent:",mean(accepted[(len - len/10):len])))
        print(paste("Average Log Ratio:",mean(likelihoods)))
        print(paste("Average Log Uniform Draw:",mean(lud)))
        scatter.smooth(x = 1:length(accepted), y = accepted, main = "Log Likelihood Ratio of Accepted Proposals",ylab= "Value",pch = 20)
        plot(likelihoods, main = "Red represents rejected proposals, blue represents accepted proposals",ylab= "Value",pch = 20, col = colors)
    }
    
    #function to plot mean beta values over time
    plot_beta_estimates <- function(topic){
        betas <- matrix(0,ncol = 4, nrow = Itterations)
        for(i in 1:Itterations){
            betas[i,] <- Metropolis_Results[[2*Itterations+i]][topic,]
        }
        mean_se <- as.data.frame(matrix(0,nrow=4, ncol = 2))
        mean_se <- cbind(c("M-M", "M-F","F-M", "F-F"),mean_se)
        names(mean_se) <- c("Tie","Parameter_Estimate","SE")
        for(j in 1:length(mean_se[,1])){
            mean_se[j,2] <- mean(betas[,j])
            mean_se[j,3] <- sd(betas[,j]) 
        }
        plot <- ggplot(mean_se, aes(x=Tie, y=Parameter_Estimate))+geom_errorbar(aes(ymin=Parameter_Estimate-SE, ymax=Parameter_Estimate+SE),col = c("black","red", "green", "blue"), width=.1, size = 1.3)+geom_point(fill= c("black","red", "green", "blue"),colour="black",pch=c(21,22,23,24), size = 4)+ labs(title = paste("Topic:",topic))
        return(list(plot))
    }
    # ================================================ #
              
    # ======= Generate Plots and save as PDFs =========#      
    make_plots <- function(Width, Height, parrow,parcol){
        #generate pdf of intercepts
        print("Plotting Intercepts...")
        pdf(file=paste(out_directory,county_name,"_Intercepts.pdf",sep = ""), width = Width, height = Height)
        par(mfrow= c(parrow,parcol))
        sapply(1:Topics,plot_intercepts)
        dev.off()
        
        #generate pdf of intercepts 
        print("Plotting Betas...")
        pdf(file=paste(out_directory,county_name,"_Betas.pdf",sep = ""), width = Width, height = Height)
        par(mfrow= c(parrow,parcol))
        sapply(1:Topics,plot_betas)
        dev.off()
        
        #generate pdf of intercepts 
        print("Plotting Latent Space Positions...")
        pdf(file=paste(out_directory,county_name,"_LS_Positions.pdf",sep = ""), width = Width, height = Height)
        par(mfrow= c(parrow,parcol))
        sapply(1:Topics,plot_LS_Positions)
        dev.off()
        
        #generate network plots
        print("Plotting Networks...")
        pdf(file=paste(out_directory,county_name,"_Network_Plots.pdf",sep = ""), width = Width, height = Height)
        par(mfrow= c(parrow,parcol))
        sapply(1:Topics,plot_Topic_Network)
        dev.off()
        
        
        #generate network plots one per page
        print("Plotting Networks...")
        pdf(file=paste(out_directory,county_name,"_Network_Plots_Full_Page.pdf",sep = ""), width = Width, height = Height)
        par(mfrow= c(1,1))
        sapply(1:Topics,plot_Topic_Network_Full)
        dev.off()
        
        
        #generate network plots
        print("Plotting Beta Estimates...")
        pdf(file=paste(out_directory,county_name,"_Beta_Estimates.pdf",sep = ""), width = Width, height = Height)
        par(mfrow= c(parrow,parcol))
        plots <- sapply(1:Topics,plot_beta_estimates)
        #for(j in 1:length(plots)){
        #    print(plots[[j]])
        #}
        multiplot(plotlist = plots, cols = 10)
        dev.off()
    }
    
    #stolen from the R cookbook
    multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
        require(grid)
        
        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)
        
        numPlots = length(plots)
        
        # If layout is NULL, then use 'cols' to determine layout
        if (is.null(layout)) {
            # Make the panel
            # ncol: Number of columns of plots
            # nrow: Number of rows needed, calculated from # of cols
            layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                             ncol = cols, nrow = ceiling(numPlots/cols), byrow= T)
        }
        
        if (numPlots==1) {
            print(plots[[1]])
            
        } else {
            # Set up the page
            grid.newpage()
            pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
            
            # Make each plot, in the correct location
            for (i in 1:numPlots) {
                # Get the i,j matrix positions of the regions that contain this subplot
                matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                
                print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                layout.pos.col = matchidx$col))
            }
        }
    }
     
    #now actually generate the plots based on the number of topics     
    if(Topics <= 25){
        make_plots(12,12,5,5)
    }else if(Topics <= 50){
        make_plots(25,12,5,10)
    }else if(Topics <= 100){
        make_plots(25,25,10,10)
    }else{
        make_plots(50,25,10,20)   
    }
         
    par(mfrow= c(1,1))
    
    
    # ======= Now generate top words and topic model diagnostics ====== #
    
    print("Generating Topic Model Output")
    
    
    
    #Establish Ordering 
    ordering <- order(Email_Assignments, decreasing = FALSE)
    Email_Assignments <- Email_Assignments[ordering]
    Top_Ten_Words <- Top_Ten_Words[ordering]
    
    
    
    pdf(paste(out_directory,county_name,"_Top_Words.pdf", sep = ""),height=12,width=10,pointsize=7)
    par(mfrow= c(1,1))
    bp <- barplot2(Email_Assignments, beside = TRUE, horiz = TRUE, col = "lightblue1", border= "lightblue1", ylab = "Topic", xlab = "Number of Edges Assigned to Topic",main = paste("Top Words by Number of Words Assigned to Topic for ",county_name,sep = "")) 
    text(0,bp,Top_Ten_Words,cex=1,pos=4)# label on the bars themselves 
    dev.off()
    
    pdf(paste(out_directory,county_name,"_Likelihoods.pdf", sep = ""),height=8,width=12,pointsize=7)
    par(mfrow= c(1,1))
    plot_likelihoods(Itterations) 
    dev.off()
    
    pdf(paste(out_directory,county_name,"_Log_Ratio_LUD.pdf", sep = ""),height=8,width=12,pointsize=7)
    par(mfrow= c(2,1))
    plot_ratio_lud(Itterations) 
    dev.off()
    
    
    #check to see if likelihood ever goes down
#     last <- -1000000000
#     for(i in 1:Itterations){
#         cur <- Metropolis_Results[[5*Itterations+i]]
#         if(cur < last){
#             print(paste("Likelihood decreased at itteration:",i,"Current:",cur,"Last:",last))
#         }
#         last <- cur
#     }
     
}#end of function definition



