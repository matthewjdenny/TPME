
library(Rcpp)

cppFunction('
            double Log_Porbability_Of_Edge_Cpp(int number_of_latent_dimensions, int actual_edge, double topic_intercept ,NumericVector author_position ,NumericVector recipient_position){

            double distance = 0;
            
            for(int k = 0; k < number_of_latent_dimensions; ++k){
                distance += pow((author_position[k] - recipient_position[k]),2);
            }

            double eta = topic_intercept - pow(distance,.5);
            double log_prob = 0;
            if(eta > 0){
                if(actual_edge == 1){
                    log_prob = eta -log(1 + exp(eta));
                }
                else{
                    log_prob = -log(1 + exp(eta));
                }
            }
            else{
                if(actual_edge == 1){
                    log_prob = -log(1 + exp(-eta));
                }
                else{
                    log_prob = -eta -log(1 + exp(-eta));
                }
            }

            return log_prob;
            }
            ')




Log_Probability_Of_Edge <- function(topic,author,recipient,edge_present,for_new_intercept = 0,proposal = 0){
    
    
    
    #if we know its a new value, then 
    if(proposal == 1){
        proposed_author_position <- Proposed_Edge_Information[topic,author,recipient][[1]]
        proposed_recipient_position <- Proposed_Edge_Information[topic,author,recipient][[2]]
        
        #just calculate porbability of edge using the cached values for 
        log_prob_of_edge <- Log_Porbability_Of_Edge_Cpp(Current_Edge_Information[topic,author,recipient][[6]],edge_present, Latent_Space_Intercepts[topic],proposed_author_position,proposed_recipient_position)
    }else{
        #get the current information on latent positions and edge log likelihood -- we will need this regardless
        stored_author_position <- Current_Edge_Information[topic,author,recipient][[1]]
        stored_recipient_position <- Current_Edge_Information[topic,author,recipient][[2]]
        stored_intercept <- Current_Edge_Information[topic,author,recipient][[3]]
        
        #if for_new_intercept is true, then we know we have to calculate the probability regardless of the latent postions so we jsut jump right to that. If not, then we need to check
        
        if(for_new_intercept == 1){
            #just calculate porbability of edge using the cached values for 
            log_prob_of_edge <- Log_Porbability_Of_Edge_Cpp(Current_Edge_Information[topic,author,recipient][[6]],edge_present, Latent_Space_Intercepts[topic],stored_author_position,stored_recipient_position)
        }else{
            #get author and recipeint current latent positons
            proposed_author_LS_position <- Latent_Space_Positions[,topic,author]
            proposed_recipient_LS_position <- Latent_Space_Positions[,topic,recipient]
            
            #check them against the latent postions associtated with the edge likelihood stored in the Current_Edge_Log_Probability array. If they are the same, then just return that value, otherwise, go ahead and calculate the value. 
            if(identical(stored_author_position,proposed_author_LS_position) & identical(stored_recipient_position,proposed_recipient_LS_position) & stored_intercept == Latent_Space_Intercepts[topic]){
                #if everything is the same jsut check to see if the y value is the same and if it is not, then calculate its reciprocal
                if(edge_present == Current_Edge_Information[topic,author,recipient][[4]]){
                    log_prob_of_edge <- Current_Edge_Information[topic,author,recipient][[5]]
                }else{
                    log_prob_of_edge <- 1- exp(Current_Edge_Information[topic,author,recipient][[5]])
                }
                
            }else{ #if the stored values are not the same as the input values then we have to calculate
                log_prob_of_edge <- Log_Porbability_Of_Edge_Cpp(Current_Edge_Information[topic,author,recipient][[6]],edge_present, Latent_Space_Intercepts[topic],proposed_author_LS_position,proposed_recipient_LS_position)
            }
            
            
        }#end else statement if we are not just checking for a new intercept where we know we have to check

    
    }#end of else statement for calculating for a proposed value
    
  return(log_prob_of_edge)
}
