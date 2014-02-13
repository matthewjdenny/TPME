#implement Gibbs Sampling in Rcpp for an 18-20x speedup over the straight R script.

#alpha, the current token index, the maximum toke index reached already, current token number, the vector of token probabilities in each topic, the topic assignments already.


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




Log_Probability_Of_Edge <- function(topic,author,recipient,for_new_intercept = F){
    
    #get the current information on latent positions and edge log likelihood -- we will need this regardless
    stored_author_position <- Proposed_Edge_Log_Probability[topic,author,recipient,,1]
    stored_author_position <- stored_author_position[-length(stored_author_position)]
    
    stored_recipient_position <- Proposed_Edge_Log_Probability[topic,author,recipient,,2]
    stored_recipient_position <- stored_recipient_position[-length(stored_recipient_position)]
    
    stored_intercept <- tail(Proposed_Edge_Log_Probability[topic,author,recipient,-1,2],1)
    
    #if for_new_intercept is true, then we know we have to calculate the probability regardless of the latent postions so we jsut jump right to that. If not, then we need to check
    
    if(for_new_intercept){
        #just calculate porbability of edge
        log_prob_of_edge <- Log_Porbability_Of_Edge_Cpp(length(stored_author_position),1, 1.234,c(1.7,.36,-.89),c(2,-1,-.5))
        
    }
    
    #get author and recipeint current latent positons
    
    
    #check them against the latent postions associtated with the edge likelihood stored in the Current_Edge_Log_Probability array. If they are the same, then just return that value, otherwise, go ahead and calculate the value.  
    
    Latent_Space_Positions
    
    
    
    Current_Edge_Log_Probability
    
    Log_Porbability_Of_Edge_Cpp(2,1, 1.234,c(1.7,.36,-.89),c(2,-1,-.5))
    
    test <- array(0,c(2,2,2,2))
    tes
    
}
