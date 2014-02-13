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




Log_Probability_Of_Edge <- function(topic,author,recipient){
    
    
    
    Log_Porbability_Of_Edge_Cpp(2,1, 1.234,c(1.7,.36,-.89),c(2,-1,-.5))
    
    
    
}
