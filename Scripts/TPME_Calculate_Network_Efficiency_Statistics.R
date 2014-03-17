Calculate_Network_Efficiency_Statistics <- function(input_file = "Current_Itteration_McDowell_2011_3-7-14", output_file = "McDowell_2011_3-7-14", data_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", method = "Latora-Machiori",Thin_Itterations = 1,skip_first = 6200){
    
    #rm(list = ls())
    #for testing load("~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/Transylvania_Sample_10M_2011_3-13-14.Rdata")
    #Thin_Itterations = 1
    #skip_first = 6200
    
    #load in data
    load(paste(data_directory,input_file,".Rdata", sep = ""))
    
    print("Loading Data")
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
    
    
    Itterations <- length(Metropolis_Results)/7
    
    #get model information and extract data
    Latent_Dimensions <- length(Metropolis_Results[[1]][,1,1])
    Number_Of_Topics <- length(Metropolis_Results[[1]][1,,1])
    Number_Of_Authors <- length(Metropolis_Results[[1]][1,1,])
    Token_Topic_Assignments <- Topic_Model_Results[[1]]
    Topic_Present_Edge_Counts <- Topic_Model_Results[[2]]
    Topic_Absent_Edge_Counts <- Topic_Model_Results[[3]]
    Word_Type_Topic_Counts <- Topic_Model_Results[[4]]
    Edge_Topic_Assignments <- Topic_Model_Results[[5]]
    Latent_Space_Positions <- Metropolis_Results[[1000]]
    Latent_Space_Intercepts <- Metropolis_Results[[2*1000]]
    Betas <- Metropolis_Results[[3*1000]]
    
    require(RcppArmadillo)
    set.seed(1234)
    
    
    #calculate the efficiency of static communication networks using the method presented by: Latora, V., & Marchiori, M. (2001). Efficient Behavior of Small-World Networks. Physical Review Letters. Retrieved from http://prl.aps.org/abstract/PRL/v87/i19/e198701.
    if(method == "Latora-Machiori"){
        
        get_beta_estimates <- function(topic){
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
            
            #now calculate whether the parameter estimate was signifincatn and whether it was significcant positive or negative
            ret <- rep(0,4)
            for(k in 1:length(mean_se[,1])){
                if(mean_se[k,2] >0){
                    if((mean_se[k,2] - 1.645*mean_se[k,3]) > 0 ){
                        ret[k] <- 1
                    }else{
                        ret[k] <- 0
                    }
                }else{
                    if(mean_se[k,2] + 1.645*mean_se[k,3] < 0){
                        ret[k] <- -1
                    }else{
                        ret[k] <- 0
                    }
                }
                
            }
            return(mean_se[,2])
        }
        
        #takes a weighted adjacency matrix as its argument
        calculate_global_efficiency <- function(net){
            #load igraph
            library(igraph)
            #create a network object
            network <- graph.adjacency(net , mode = "directed", weighted = TRUE)
            #calculate the shortest path lengths between nodes in the network
            d_ij <- shortest.paths(network,mode = "all", weights = network$weights)
            num_nodes <- length(net[1,])
            E_G <- 0
            for(i in 1: num_nodes){
                for(j in 1:num_nodes){
                    if(i != j){
                        E_G <- E_G + 1/d_ij[i,j]
                    }
                }
            }
            #normalize
            E_G <- E_G/(num_nodes*(num_nodes -1))
            if(sum(net) > 0){
                E_G <- E_G/sum(net)
            }
            return(E_G)
        }
        
        #calculate efficiency and beta statistics
        efficiencies <- rep(0,Number_Of_Topics)
        beta_averages <- matrix(0, ncol = 4,nrow = Number_Of_Topics)
        for(i in 1:Number_Of_Topics){
            topic <- i
            beta_averages[topic,] <- get_beta_estimates(topic)
            efficiencies[topic] <- calculate_global_efficiency(Topic_Present_Edge_Counts[,,topic])
        }
        
        
        
        print("Plotting Beta-efficiency scatter plots...")
        pdf(file=paste(data_directory,output_file,"_Efficiency_Plots.pdf",sep = ""), width = 16, height = 16)
        par(mfrow= c(2,2))
        scatter.smooth(x= beta_averages[,1], y= efficiencies, pch =20, col = "red", main = paste("Male-Male Mixing Parameter Estimates vs. Tie Weighted Network Efficiency \n", "Pearson Correlation Test p-value:",round(cor.test(beta_averages[,1],efficiencies)$p.value,3)), ylab = "Per-Tie Efficiency", xlab = "Topic Specific Parameter Estimate")
        scatter.smooth(x= beta_averages[,2], y= efficiencies,pch =20, col = "red", main = paste("Male-Female Mixing Parameter Estimates vs. Tie Weighted Network Efficiency \n", "Pearson Correlation Test p-value:",round(cor.test(beta_averages[,2],efficiencies)$p.value,3)) , ylab = "Per-Tie Efficiency", xlab = "Topic Specific Parameter Estimate")
        scatter.smooth(x= beta_averages[,3], y= efficiencies,pch =20, col = "red", main = paste("Female-Male Mixing Parameter Estimates vs. Tie Weighted Network Efficiency \n", "Pearson Correlation Test p-value:",round(cor.test(beta_averages[,3],efficiencies)$p.value,3)), ylab = "Per-Tie Efficiency", xlab = "Topic Specific Parameter Estimate")
        scatter.smooth(x= beta_averages[,4], y= efficiencies,pch =20, col = "red", main = paste("Female-Female Mixing Parameter Estimates vs. Tie Weighted Network Efficiency \n", "Pearson Correlation Test p-value:",round(cor.test(beta_averages[,4],efficiencies)$p.value,3)), ylab = "Per-Tie Efficiency", xlab = "Topic Specific Parameter Estimate")
        dev.off()
    }
    
    
    
    
}#end of function definition
    
    
    
       
        