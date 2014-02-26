#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
List Topic_Assignment_Step_CPP(
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
    NumericMatrix word_type_topic_counts
    ){
    
    
    //This function handles all of the token topic assingment sampling as well as the edge topic assignment sampling
    
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
    
    //loop over the number of metropolis itterations (default 1000)
    for(int i = 0; i < number_of_itterations; ++i){
        
        
        // ========== token topic assignment step =============
        for(int d = 0; d < number_of_documents; ++d){
            
        }
        
        
        
        // ============= edge topic assignment step ==============
        for(int d = 0; d < number_of_documents; ++d){
            
            //take care of document specific assignment tasks
            int document_author = document_authors[d] - 1;
            NumericVector edge_topic_assignments(number_of_actors);
            NumericVector token_topic_assignments = token_topic_assignment_list[d];
            int number_of_tokens = token_topic_assignments.length();
            
            
            
            for(int a = 0; a < number_of_actors; ++a){
                if(document_author != a){
                    int recipient = a + 1;
                    
                    NumericVector edge_log_probabilities(number_of_tokens);
                    for(int w = 0; w < number_of_tokens; ++w){
                        int token = w + 1;
                        int topic_assignment = token_topic_assignments[w] -1;
                        //int topic_assignment = as<int>(get_token_topic_assignment(document,token));
                        int actual_edge = observed_edges(d,a);
                        //get current topic intercept
                        double current_topic_intercept = current_intercepts[topic_assignment];
                        //get current topic betas
                        NumericVector current_topic_betas = betas.row(topic_assignment);
    
                        //get current author position
                        current_author_position[0] = current_latent_positions(0,topic_assignment,document_author);
                        current_author_position[1] = current_latent_positions(1,topic_assignment,document_author);
                        //get current recipient position
                        recipient_position[0] = current_latent_positions(0,topic_assignment,a);
                        recipient_position[1] = current_latent_positions(1,topic_assignment,a);
                        
                        //initialize distance
                        double distance = 0;
                        //calculate distance
                        for(int k = 0; k < number_of_latent_dimensions; ++k){
                            distance += pow((current_author_position[k] - recipient_position[k]),2);
                        }
                        
                        beta_val = 0;
                        for(int c = 0; c < number_of_betas; ++c){
                            beta_val += current_topic_betas[c]*beta_indicator_array(a,b,c);
                        }
                        //calculate linear predictor
                        double eta = current_topic_intercept - pow(distance,.5) + beta_val;
                        
                        double log_prob = 0;
                        if (eta != 0){
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
                        }
                        
                        edge_log_probabilities[w] = log_prob;
                    }
                    int sampled_token = as<int>(log_multinomial_draw(edge_log_probabilities));
                    edge_topic_assignments(d,a) = token_topic_assignments[sampled_token-1];
                    
                    //do some updating 
                }
                
            }//end of loop over edges for for current document
        }// end of loop over docuemnts for edge-topic assignemnt step
        
        
    }//end of big number of itterations loop for entire step.    

    return to_return;
}

//Metropolis_Step_CPP(30,100,array(1:90000,c(30,30,100)),array(1:90000,c(30,30,100)),array(runif(6000),c(2,100,30)),c(1:100),2,array(runif(6000),c(2,100,30)),c(1:100))
            


