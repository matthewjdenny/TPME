#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector SAMPLE_SINGLE_TOKEN_TOPIC_ASSIGNMENT_CPP(
    int number_of_tokens, 
    int number_of_topics,
    int number_of_actors, 
    int author, 
    int document, 
    double beta, 
    NumericVector alpha_m, 
    NumericMatrix latent_positions1, 
    NumericMatrix latent_positions2, 
    NumericVector intercepts, 
    int number_of_latent_dimensions, 
    NumericVector edge_topic_assignments, 
    int token_topic_assignments,
    NumericVector topic_token_sums,
    NumericVector token_type_topic_counts,
    int number_of_word_types,
    NumericVector observed_edges
    ){
    //this function is special for docuemnts with only one token because having matricies
    //of one row crashes c++ so we need a special function jsut for this case. 
    //only implemented with two dimensions. This function makes one call to R to do multinomial draws.
    //This appears to have almost zero overhead in R so it should not have a noticeable effect on runtime.
    Function log_multinomial_draw("log_multinomial_draw");
    
    //sample new topic assignment for each token
    for(int w = 0; w < number_of_tokens; ++w){
        //account for R indexing starting at one and c++ starting at 0
        int token = w + 1;
        
        //initialize the probability vector to be fed into the sampler
        NumericVector token_topic_distribution(number_of_topics);
        for(int t = 0; t < number_of_topics; ++t){
            int topic = t + 1;
            
            //get author latent position for this topic
            NumericVector author_position(number_of_latent_dimensions);
            author_position[0] = latent_positions1(t,(author-1));
            author_position[1] = latent_positions2(t,(author-1));
            //get topic specific intercept
            double topic_intercept = intercepts[t];
            //this calculates the addition to the probability that the token was sampled from the topic by
            //adding together edge likelihoods associated with that topic in the document
            double additional_edge_probability = 0;
            //go through and check for an edge associated with each other actor
            for(int a = 0; a < number_of_actors; ++a){
                int actor = a + 1;
                //do not check for edges from the author to the author
                if(author != actor){
                    //get the topic assignemnt of the current author-actor edge
                    int edge_assignment = edge_topic_assignments[a];
                    //if the edge is assigned in the current topic, then we get its likelihood
                    if(edge_assignment == topic){
                        
                        //get whether the edge was actually observed so this can be factored into the probability calculation
                        int actual_edge = observed_edges[a];
                        
                        //get recipient latent position for this topic
                        NumericVector recipient_position(number_of_latent_dimensions);
                        recipient_position[0] = latent_positions1(t,a);
                        recipient_position[1] = latent_positions2(t,a);
                        
                        //intitalize the distance
                        double distance = 0;
                        //calculate the distance
                        for(int k = 0; k < number_of_latent_dimensions; ++k){
                            distance += pow((author_position[k] - recipient_position[k]),2);
                        }
                        //subtract the square root of the distance from the intercept
                        double eta = topic_intercept - pow(distance,.5);
                        double log_prob = 0;
                        if(eta > 0){
                            if(actual_edge == 1){
                                log_prob = eta -log(1 + exp(eta));
                            }
                            else{
                                log_prob = 0 -log(1 + exp(eta));
                            }
                        }
                        else{
                            if(actual_edge == 1){
                                log_prob = 0 -log(1 + exp(-eta));
                            }
                            else{
                                log_prob = 0 -eta -log(1 + exp(-eta));
                            }
                        }
                        additional_edge_probability += log_prob;
                    }   
                } 
            }
            
            
            //number of times the topic appears in the document minus the only token will always be zero
            int ntd = 0;
            //number of times this token has been assigned to the topic across all documents
            int wttac = token_type_topic_counts[t];
            //the number of tokens assigned to the topic across all docuemnts
            int ntt = topic_token_sums[t];
            //decrement if the current token accounts for one of these
            if(topic == token_topic_assignments){
                ntt -= 1;
                wttac -=1;
            }
            //now we calculate the first and second terms in the likelihood of of the token being from the current topic
            double first_term = ntd + (alpha_m[t]/number_of_topics);
            double second_term = (wttac + (beta/number_of_word_types))/(ntt + beta);
            token_topic_distribution[t] = log(first_term) + log(second_term) + additional_edge_probability;
        }
        //sample the new token topic assignment
        token_topic_assignments = as<double>(log_multinomial_draw(token_topic_distribution));        
    }
            
    return token_topic_assignments;
}


//d= 1
//SAMPLE_TOKEN_TOPIC_ASSIGNMENTS_CPP(Token_Topic_Assignments[[d]][[1]],Number_Of_Topics,Number_Of_Authors,Document_Authors[d],d,Beta,Alpha_Base_Measure_Vector)



