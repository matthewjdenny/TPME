
Generate_Model_Diagnsotics <- function(output_file = "Test",Actors = 17,Itterations = 1000, Latent_Spaces = 2, Topics = 50,LS_Actor = 2, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/",Thin_Itterations = 1){
    #load current results
    #load("Current_Itteration_Results.Rdata")
    print("Loading Data...")
    load(paste("./Output/",output_file,".Rdata", sep = ""))
    
    print("Extracting Reduced Data")
    #extract current metropolis results
    Metropolis_Results <- Return_List[[1]]
    Topic_Model_Results <- Return_List[[2]]
    #free up memory
    rm(Return_List)
    
    #thin out the data by taking every Thin_Itterations itteration for the metropolis step
    Metropolis_Results <- Metropolis_Results[seq(1, length(Metropolis_Results),Thin_Itterations)]
         
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
        plot(intercepts, main = paste("Topic:",intercept),ylab= "Value")
    }
    #function to plot betas over time
    plot_betas <- function(topic){
        betas <- matrix(0,ncol = 4, nrow = Itterations)
        for(i in 1:Itterations){
            betas[i,] <- Metropolis_Results[[2*Itterations+i]][topic,]
        }
        matplot(betas, main = paste("Topic:",topic),ylab= "Value")
    }
    #function to plot latent spaces for a given actor
    plot_LS_Positions <- function(topic){
        LSP <- matrix(0,nrow=Itterations, ncol= Latent_Spaces)
        for(i in 1:Itterations){
            LSP[i,] <- Metropolis_Results[[i]][,topic,LS_Actor]
        }
        corel <- cor(LSP[,1], LSP[,2])
        matplot(LSP, main = paste("Topic:",topic, "\n LS Correlations:", round(corel,3)),ylab= "Value")
    }
    # ================================================ #
              
    # ======= Generate Plots and save as PDFs =========#      
    make_plots <- function(Width, Height, parrow,parcol){
        #generate pdf of intercepts
        
        pdf(file=paste(out_directory,"Intercepts.pdf",sep = ""), width = Width, height = Height)
        par(mfrow= c(parrow,parcol))
        sapply(1:Topics,plot_intercepts)
        dev.off()
        
        #generate pdf of intercepts   
        pdf(file=paste(out_directory,"Betas.pdf",sep = ""), width = Width, height = Height)
        par(mfrow= c(parrow,parcol))
        sapply(1:Topics,plot_betas)
        dev.off()
        
        #generate pdf of intercepts   
        pdf(file=paste(out_directory,"LS_Positions.pdf",sep = ""), width = Width, height = Height)
        par(mfrow= c(parrow,parcol))
        sapply(1:Topics,plot_LS_Positions)
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
    
    
    
    
     
}#end of function definition



