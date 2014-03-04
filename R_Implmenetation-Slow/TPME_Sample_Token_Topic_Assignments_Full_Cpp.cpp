#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector SAMPLE_TOKEN_TOPIC_ASSIGNMENTS_CPP(
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
    NumericVector token_topic_assignments,
    NumericVector topic_token_sums,
    NumericMatrix token_type_topic_counts,
    int number_of_word_types,
    NumericVector observed_edges
    ){
        
    //only implemented with two dimensions
    //, NumericVector cur_token_topic_assignments, NumericVector edge_values, NumericVector edge_topic_assignments
    Function get_token_topic_assignment("get_token_topic_assignment");
    Function get_observed_edge_value("get_observed_edge_value");
    Function log_multinomial_draw("log_multinomial_draw");
    Function get_edge_topic_assignment("get_edge_topic_assignment");
    Function get_sum_token_topic_assignments("get_sum_token_topic_assignments");
    Function get_number_of_unique_words("get_number_of_unique_words");
    Function get_word_type_topic_assignemnt_count("get_word_type_topic_assignemnt_count");
    Function get_number_of_tokens_assigned_to_topic("get_number_of_tokens_assigned_to_topic");
     
    //int number_of_word_types = as<int>(get_number_of_unique_words());
    
    //NumericVector token_topic_assignments(number_of_tokens);
    //NumericVector token_topic_assignments = cur_token_topic_assignments;
    for(int w = 0; w < number_of_tokens; ++w){
        int token = w + 1;
        
        
        NumericVector token_topic_distribution(number_of_topics);
        for(int t = 0; t < number_of_topics; ++t){
            int topic = t + 1;
            
            NumericVector author_position(number_of_latent_dimensions);
            author_position[0] = latent_positions1(t,(author-1));
            author_position[1] = latent_positions2(t,(author-1));
            double topic_intercept = intercepts[t];
            //this calculates the addition to the probability that the token was sampled from the topic by
            //adding together edge likelihoods associated with that topic in the document
            double additional_edge_probability = 0;
            for(int a = 0; a < number_of_actors; ++a){
                int actor = a + 1;
                
                if(author != actor){
                    
                    int actual_edge = observed_edges[a];
                    
                    int edge_assignment = edge_topic_assignments[a]; 
                    //int edge_assignment = as<int>(get_edge_topic_assignment(document, actor));
                    if(edge_assignment == topic){
                        
                        //change one element in the array/cube
                        NumericVector recipient_position(number_of_latent_dimensions);
                        recipient_position[0] = latent_positions1(t,a);
                        recipient_position[1] = latent_positions2(t,a);
                        
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
            
            
            //now we calculate the first and second terms in the likelihood of of the token being from the current topic
            int ntd = 0;
            for(int b = 0; b < number_of_tokens; ++b){
                if(b != w){
                    if(token_topic_assignments[b] == topic){
                        ntd +=1;
                    }
                }
            }
            //int ntd = as<int>(get_sum_token_topic_assignments(document,token,topic));
            int wttac = token_type_topic_counts(w,t);
            //int wttac = as<int>(get_word_type_topic_assignemnt_count(document,token,topic));
            
            int ntt = topic_token_sums[t];
            if(topic == token_topic_assignments[w]){
                ntt -= 1;
                wttac -=1;
            }
            //int ntt = as<int>(get_number_of_tokens_assigned_to_topic(document,token,topic));
            //double first_term = ntd + (alpha_m[t]/number_of_topics);
            double first_term = ntd + alpha_m[t];
            double second_term = (wttac + (beta/number_of_word_types))/(ntt + beta);
            
            
            token_topic_distribution[t] = log(first_term) + log(second_term) + additional_edge_probability;
        }
        int old_topic = token_topic_assignments[w];
        token_topic_assignments[w] = as<double>(log_multinomial_draw(token_topic_distribution));        
        int new_topic = token_topic_assignments[w];
        
        //now we need to update all of the internal counts for the same words as the current token
        if(old_topic != new_topic){
            topic_token_sums[old_topic] -=1;
            topic_token_sums[new_topic] += 1; 
            //now for all tokens that are the same
        }
        
    }
            
    return token_topic_assignments;
}


//d= 1
//SAMPLE_TOKEN_TOPIC_ASSIGNMENTS_CPP(Token_Topic_Assignments[[d]][[1]],Number_Of_Topics,Number_Of_Authors,Document_Authors[d],d,Beta,Alpha_Base_Measure_Vector)



