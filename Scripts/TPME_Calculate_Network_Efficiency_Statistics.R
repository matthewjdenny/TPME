Calculate_Network_Efficiency_Statistics <- function(input_file = "Current_Itteration_McDowell_2011_3-7-14", output_file = "Network_Efficiency_McDowell_2011_3-7-14", data_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", method = "Latora-Machiori"){
    
    #load in data
    load(paste(data_directory,input_file,".Rdata", sep = ""))
    
    print("Loading Data")
    #extract current metropolis results
    Metropolis_Results <- Return_List[[1]]
    Topic_Model_Results <- Return_List[[2]]
    #free up memory
    rm(Return_List)
    
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
    if(method = "Latora-Machiori"){
        
        #takes a weighted adjacency matrix as its argument
        
        
        
        calculate_global_efficiency <- function(net){
            #load igraph
            library(igraph)
            #create a network object
            network <- graph.adjacency(net , mode = "undirected", weighted = TRUE)
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
            
            return(E_G)
        }
        
        
        calculate_global_efficiency(net)
        
        
    }
    
    
    
    
}
    
    
    
       
        