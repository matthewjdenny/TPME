#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Tests the array functionality of RcppArmadillo

// [[Rcpp::export]]
List test(NumericVector myArray, NumericVector probs){
    IntegerVector arrayDims = myArray.attr("dim");
    double x = Rf_rnorm(0.0,1.0);
    Function log_multinomial_draw("log_multinomial_draw");
    List out(2);
  arma::cube cubeArray(myArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
  //change one element in the array/cube
  cubeArray(0,0,0) = 518;  
  
  out[0] = cubeArray;
  out[1] = log_multinomial_draw(probs);
 
  return(out); 

}
//test(array(1:1000,c(10,10,10)), rnorm(100,5,1))

// [[Rcpp::export]]
List test2(
    int number_of_actors, 
    int number_of_topics,
    NumericVector tpec,
    NumericVector taec,
    NumericVector clp, 
    NumericVector current_intercepts,
    int number_of_latent_dimensions,
    NumericMatrix betas,
    int number_of_betas,
    NumericVector indicator_array,
    int number_of_itterations,
    double proposal_variance,
    int number_of_documents,
    NumericMatrix observed_edges,
    List token_topic_assignment_list,
    List token_word_type_list,
    NumericVector document_authors,
double beta,
NumericVector alpha_m,
NumericMatrix edge_topic_assignments,
NumericMatrix token_type_topic_counts,
NumericVector topic_token_sums,
int number_of_word_types,
NumericVector words_in_document
    ){
  
  
  //This function handles all of the token topic assingment sampling as well as the edge topic assignment sampling
    //still need to implement check for zero tokens in document
    Function log_multinomial_draw("log_multinomial_draw");
    
    IntegerVector arrayDims1 = tpec.attr("dim");
    arma::cube topic_present_edge_counts(tpec.begin(), arrayDims1[0], arrayDims1[1], arrayDims1[2], false);
    
    IntegerVector arrayDims2 = taec.attr("dim");
    arma::cube topic_absent_edge_counts(taec.begin(), arrayDims2[0], arrayDims2[1], arrayDims2[2], false);
    
    IntegerVector arrayDims3 = clp.attr("dim");
    arma::cube current_latent_positions(clp.begin(), arrayDims3[0], arrayDims3[1], arrayDims3[2], false);
    
    IntegerVector arrayDims5 = indicator_array.attr("dim");
    arma::cube beta_indicator_array(indicator_array.begin(), arrayDims5[0], arrayDims5[1], arrayDims5[2], false);

    double beta_val = 0;
    NumericVector current_author_position(number_of_latent_dimensions);
    NumericVector recipient_position(number_of_latent_dimensions);
    NumericVector current_topic_betas(number_of_betas);
  
    for(int i = 0; i < number_of_itterations; ++i){
  
    // ========== token topic assignment step =============
        for(int d = 0; d < number_of_documents; ++d){
            //set all document specific parameters 
            int document_author = document_authors[d] - 1;
            NumericVector token_topic_assignments1 = token_topic_assignment_list[d];
            int number_of_tokens = words_in_document[d];

            for(int w = 0; w < number_of_tokens; ++w){
                int token = w + 1;
                NumericVector token_topic_distribution(number_of_topics);
                for(int t = 0; t < number_of_topics; ++t){
                    int topic = t + 1;
                    
                    current_author_position[0] = current_latent_positions(0,t,document_author);
                    current_author_position[1] = current_latent_positions(1,t,document_author);
                    double current_topic_intercept = current_intercepts[t];
                    NumericVector current_topic_betas = betas.row(t);
                    //this calculates the addition to the probability that the token was sampled from the topic by
                    //adding together edge likelihoods associated with that topic in the document
                    double additional_edge_probability = 0;
                    for(int a = 0; a < number_of_actors; ++a){
                        
                        if(document_author != a){
                            
                            int actual_edge = observed_edges(d,a);
                            
                            int edge_assignment = edge_topic_assignments(d,a); 
                            //int edge_assignment = as<int>(get_edge_topic_assignment(document, actor));
                            if(edge_assignment == topic){
                                
                                //get current recipient position
                                recipient_position[0] = current_latent_positions(0,t,a);
                                recipient_position[1] = current_latent_positions(1,t,a);
                                
                                //initialize distance
                                double distance = 0;
                                //calculate distance
                                for(int k = 0; k < number_of_latent_dimensions; ++k){
                                    distance += pow((current_author_position[k] - recipient_position[k]),2);
                                }
                        
                                beta_val = 0;
                                for(int c = 0; c < number_of_betas; ++c){
                                    beta_val += current_topic_betas[c]*beta_indicator_array(document_author,a,c);
                                }
                                //calculate linear predictor
                                double eta = current_topic_intercept - pow(distance,.5) + beta_val;
                        
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
                            if(token_topic_assignments1[b] == topic){
                                ntd +=1;
                            }
                        }
                    }
                    int wttac = token_type_topic_counts(w,t);
                    int ntt = topic_token_sums[t];
                    if(topic == token_topic_assignments1[w]){
                        ntt -= 1;
                        wttac -=1;
                    }
                    //int ntt = as<int>(get_number_of_tokens_assigned_to_topic(document,token,topic));
                    //double first_term = ntd + (alpha_m[t]/number_of_topics);
                    double first_term = ntd + alpha_m[t];
                    double second_term = (wttac + (beta/number_of_word_types))/(ntt + beta);

                    //token_topic_distribution[t] = log(first_term) + log(second_term) + additional_edge_probability;
                    token_topic_distribution[t] = wttac;
                }

                int old_topic = token_topic_assignments1[w];
                //token_topic_assignments1[w] 
                int tester = as<int>(log_multinomial_draw(token_topic_distribution));        
                int new_topic = token_topic_assignments1[w];
                //now we need to update all of the internal counts for the same words as the current token
                if(old_topic != new_topic){
                    topic_token_sums[(old_topic-1)] -=1;
                    topic_token_sums[(new_topic-1)] += 1; 
                    //now for all tokens that are the same
                    token_type_topic_counts(w,(old_topic-1)) -=1;
                    token_type_topic_counts(w,(new_topic-1)) +=1;
                }
            }//end of loop over tokens 
        //now assign the new token topic assignments to the appropriate place in the list.
        token_topic_assignment_list[d] = token_topic_assignments1; 
        }//end of token topic assignemnt sampler document loop
    }
    List to_return(3);
    to_return[1] = words_in_document;
    return to_return;
    }  
  
    
    
    
    