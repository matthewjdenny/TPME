
Generate_Model_Diagnsotics <- function(input_folder_path = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", input_file = "Test",LS_Actor = 2, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/",Thin_Itterations = 20, vocab = vocabulary,county_name = "McDowell_County"){
    #load current results
    #load("Current_Itteration_Results.Rdata")
    library(statnet)
    library(gregmisc)
    print("Loading Data...")
    load(paste(input_folder_path,input_file,".Rdata", sep = ""))
    
    print("Extracting Reduced Data")
    #extract current metropolis results
    Metropolis_Results <- Return_List[[1]]
    Topic_Model_Results <- Return_List[[2]]
    #free up memory
    rm(Return_List)
    
    #thin out the data by taking every Thin_Itterations itteration for the metropolis step
    Metropolis_Results <- Metropolis_Results[seq(1, length(Metropolis_Results),Thin_Itterations)]
    
    #get model information and extract data
    Itterations <- length(Metropolis_Results)/6
    Latent_Spaces <- length(Metropolis_Results[[1]][,1,1])
    Topics <- length(Metropolis_Results[[1]][1,,1])
    Actors <- length(Metropolis_Results[[1]][1,1,])
    Token_Topic_Assignments <- Topic_Model_Results[[1]]
    Topic_Present_Edge_Counts <- Topic_Model_Results[[2]]
    Topic_Absent_Edge_Counts <- Topic_Model_Results[[3]]
    Word_Type_Topic_Counts <- Topic_Model_Results[[4]]
    Edge_Topic_Assignments <- Topic_Model_Results[[5]]
    
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
        plot(intercepts, main = paste("Topic:",intercept),ylab= "Value",pch = 20)
    }
    #function to plot betas over time
    plot_betas <- function(topic){
        betas <- matrix(0,ncol = 4, nrow = Itterations)
        for(i in 1:Itterations){
            betas[i,] <- Metropolis_Results[[2*Itterations+i]][topic,]
        }
        matplot(betas, main = paste("Topic:",topic),ylab= "Value",pch = 20)
    }
    #function to plot latent spaces for a given actor
    plot_LS_Positions <- function(topic){
        LSP <- matrix(0,nrow=Itterations, ncol= Latent_Spaces)
        for(i in 1:Itterations){
            LSP[i,] <- Metropolis_Results[[i]][,topic,LS_Actor]
        }
        corel <- cor(LSP[,1], LSP[,2])
        matplot(LSP, main = paste("Topic:",topic, "\n LS Correlations:", round(corel,3)),ylab= "Value",pch = 20)
    }
    plot_Topic_Network <- function(topic){
        slice <- Topic_Present_Edge_Counts[,,topic]
        net <- as.network(slice)
        plot(net, main = paste("Topic:",topic,"Number of Edges:",sum(slice)),pch = 20)
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
    #get the total number of tokens assigned to each topic
    Topic_Token_Totals <- apply(Word_Type_Topic_Counts,2,sum)
    #get the total number of present edges assigned to each topic
    Email_Assignments <- apply(Topic_Present_Edge_Counts,3,sum)
    #number of words
    num_words <- length(vocab)
    
    
    #this list will be used to store lists of three vectors: word indicies in descending order of likelihood, word probability in topic, actual word. 
    Topic_Top_Words <- list()
    Top_Ten_Words <- rep("",Topics)
    
    
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
        for(j in 1:10){
            words <- paste(words,top_words[j],", ",sep = "")
        }
        Top_Ten_Words[i] <- words
    }
    
    
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
    
    #check to see if likelihood ever goes down
    last <- -1000000000
    for(i in 1:Itterations){
        cur <- Metropolis_Results[[5*Itterations+i]]
        if(cur < last){
            print(paste("Likelihood decreased at itteration:",i,"Current:",cur,"Last:",last))
        }
        last <- cur
    }
    
    
    
     
}#end of function definition



