#include <RcppArmadillo.h>
#include <random>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
List Cluster_Integrated_Sampler(
    int number_of_outer_itterations,
    int number_of_Gibbs_itterations,
    int number_of_MH_itterations,
    int number_of_actors, 
    int number_of_topics,
    int number_of_clusters,
    int number_of_latent_dimensions,
    int number_of_documents,
    double proposal_variance,
    NumericVector topic_cluster_assignments,
    NumericVector tpec,
    NumericVector taec,
    NumericVector clp, 
    NumericVector plp
    NumericVector current_intercepts,
    NumericMatrix betas,
    int number_of_betas,
    NumericVector indicator_array,
    NumericMatrix observed_edges,
    List token_topic_assignment_list,
    List token_word_type_list,
    NumericVector document_authors,
    double beta,
    NumericVector alpha_m,
    NumericMatrix edge_topic_assignments,
    NumericMatrix token_type_topic_counts,
    NumericVector topic_token_sums,
    int number_of_word_types
    ){
        
    //This function handles all of the token topic assingment sampling as well as the edge topic assignment sampling
    //still need to implement check for zero tokens in document
    Function log_multinomial_draw("log_multinomial_draw");
    
    //set seed and designate RNG
    srand((unsigned)time(NULL));
    std::default_random_engine generator;
    
    //read in topic present edge counts array [num actors x num actors x topics]
    IntegerVector arrayDims1 = tpec.attr("dim");
    arma::cube topic_present_edge_counts(tpec.begin(), arrayDims1[0], arrayDims1[1], arrayDims1[2], false);
    
    //read in topic absent edge counts array [num actors x num actors x topics]
    IntegerVector arrayDims2 = taec.attr("dim");
    arma::cube topic_absent_edge_counts(taec.begin(), arrayDims2[0], arrayDims2[1], arrayDims2[2], false);
    
    //read in latent positions array [num dimensions x topics x actors]
    IntegerVector arrayDims3 = clp.attr("dim");
    arma::cube current_latent_positions(clp.begin(), arrayDims3[0], arrayDims3[1], arrayDims3[2], false);
    
    //read in beta indicator array (0,1) [number of topics x number of actors x number of betas]
    IntegerVector arrayDims5 = indicator_array.attr("dim"); 
    arma::cube beta_indicator_array(indicator_array.begin(), arrayDims5[0], arrayDims5[1], arrayDims5[2], false);

    //initialize variables that will be used across iterations
    double beta_val = 0;
    NumericVector current_author_position(number_of_latent_dimensions);
    NumericVector recipient_position(number_of_latent_dimensions);
    NumericVector current_topic_betas(number_of_betas);
    
    //loop over the number of metropolis itterations (default 1000)
    for(int i = 0; i < number_of_itterations; ++i){
        
        
        // ========== token topic assignment step =============
        for(int d = 0; d < number_of_documents; ++d){
            //set all document specific parameters 
            int document_author = document_authors[d] - 1;
            NumericVector token_topic_assignments1 = token_topic_assignment_list[d];
            int number_of_tokens = token_topic_assignments1.length();
            NumericVector token_word_types = token_word_type_list[d];
            
            //loop over tokens
            for(int w = 0; w < number_of_tokens; ++w){
                
                //initialize token topic distribution vector
                NumericVector token_topic_distribution(number_of_topics);
                
                //loop over topics
                for(int t = 0; t < number_of_topics; ++t){
                    int topic = t + 1;
                    
                    //get current author latent positions (could be updated to run in a loop)
                    current_author_position[0] = current_latent_positions(0,t,document_author);
                    current_author_position[1] = current_latent_positions(1,t,document_author);
                    
                    //get current topic betas and intercepts
                    double current_topic_intercept = current_intercepts[t];
                    NumericVector current_topic_betas = betas.row(t);
                    //this calculates the addition to the probability that the token was sampled from the topic by
                    //adding together edge likelihoods associated with that topic in the document
                    double additional_edge_probability = 0;
                    for(int a = 0; a < number_of_actors; ++a){
                        
                        //no self loops
                        if(document_author != a){
                            
                            //get whether or not there was an observed edge
                            int actual_edge = observed_edges(d,a);
                            //get the edge's current topic assignment
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
                                //pull out the correct beta value 
                                for(int c = 0; c < number_of_betas; ++c){
                                    beta_val += current_topic_betas[c]*beta_indicator_array(document_author,a,c);
                                }
                                //calculate linear predictor
                                double eta = current_topic_intercept - pow(distance,.5) + beta_val;
                        
                                double log_prob = 0;
                                if(eta < 0){
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
                    //calculate the number of times a token in the current document has been assigned to the current topic
                    int ntd = 0;
                    for(int b = 0; b < number_of_tokens; ++b){
                        if(b != w){
                            if(token_topic_assignments1[b] == topic){
                                ntd +=1;
                            }
                        }
                    }
                    
                    int current_word = token_word_types[w] -1;
                    //get the number of times this word type has been assigned in the current topic
                    int wttac = token_type_topic_counts(current_word,t);
                    //int ntd = as<int>(get_sum_token_topic_assignments(document,token,topic));
                    //int wttac = as<int>(get_word_type_topic_assignemnt_count(document,token,topic));
                    
                    //number of tokens assigned to the topic
                    int ntt = topic_token_sums[t];
                    //subtract one from these counts if the current token was assigned to this topic
                    if(topic == token_topic_assignments1[w]){
                        ntt -= 1;
                        wttac -=1;
                    }
                    
                    //helps deal with wierd initializations by no allowing negative counts
                    if(wttac < 0){
                        wttac = 0;
                    }
                    //int ntt = as<int>(get_number_of_tokens_assigned_to_topic(document,token,topic));
                    double first_term = ntd + alpha_m[t];
                    double second_term = (wttac + (beta/number_of_word_types))/(ntt + beta);
                    
                    //combine terms and add to conditional posterior distribution
                    token_topic_distribution[t] = log(first_term) + log(second_term) + additional_edge_probability;
                }
                //save old topic asignemnt
                int old_topic = token_topic_assignments1[w];
                
                //generate the un-logged distribution to sample from
                NumericVector topic_distribution(number_of_topics);
                for(int x = 0; x < number_of_topics; ++x){
                    topic_distribution[x] = exp(token_topic_distribution[x]);
                }
                //use the above to initialize a discrete distribution
                std::discrete_distribution<int> distribution (topic_distribution.begin(),topic_distribution.end());
                //take an integer topic assignment draw from the distribtuion (needed to work with some c++ compilers)
                int temp = distribution(generator) +1;
                //cast the integer draw as a double so it can be assigned to token topic assignemnt vector
                token_topic_assignments1[w] = double (temp); 
                //token_topic_assignments1[w] = as<double>(log_multinomial_draw(token_topic_distribution));        
                int new_topic = token_topic_assignments1[w];
                
                //now we need to update all of the internal counts for the same words as the current token if we sampled a new assignment
                if(old_topic != new_topic){
                    topic_token_sums[(old_topic-1)] -=1;
                    topic_token_sums[(new_topic-1)] += 1;
                    
                    int current_word_type = token_word_types[w] -1;
                    //now for all tokens that are the same
                    token_type_topic_counts(current_word_type,(old_topic-1)) -=1;
                    token_type_topic_counts(current_word_type,(new_topic-1)) +=1;
                }
                
            }//end of loop over tokens 
        
        //now assign the new token topic assignments to the appropriate place in the list.
        token_topic_assignment_list[d] = token_topic_assignments1; 
        }//end of token topic assignemnt sampler document loop
        
        
        
        // ============= edge topic assignment step ==============
        for(int d = 0; d < number_of_documents; ++d){
            
            //take care of document specific assignment tasks
            int document_author = document_authors[d] - 1;
            NumericVector token_topic_assignments2 = token_topic_assignment_list[d];
            int number_of_tokens2 = token_topic_assignments2.length();
            
            
            
            for(int a = 0; a < number_of_actors; ++a){
                if(document_author != a){
                    //int recipient = a + 1;
                    int actual_edge = observed_edges(d,a);
                    
                    NumericVector edge_log_probabilities(number_of_tokens2);
                    for(int w = 0; w < number_of_tokens2; ++w){
                        //int token = w + 1;
                        int topic_assignment = token_topic_assignments2[w] -1;
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
                            beta_val += current_topic_betas[c]*beta_indicator_array(document_author,a,c);
                        }
                        //calculate linear predictor
                        double eta = current_topic_intercept - pow(distance,.5) + beta_val;
                        
                        double log_prob = 0;
                        if (eta != 0){
                            if(eta < 0){
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
                    
                    NumericVector edge_probabilities(number_of_tokens2);
                    for(int x = 0; x < number_of_tokens2; ++x){
                        edge_probabilities[x] = exp(edge_log_probabilities[x]);
                    }
                    std::discrete_distribution<int> distribution2 (edge_probabilities.begin(),edge_probabilities.end());
                    int sampled_token = distribution2(generator);
                    //sampled_token = as<int>(log_multinomial_draw(edge_log_probabilities));
                    edge_topic_assignments(d,a) = token_topic_assignments2[sampled_token];
                    
                    //if we are on the last itteration do some updating
                    if(i == (number_of_itterations -1)){
                        if(actual_edge == 1){
                            topic_present_edge_counts(document_author,a,(token_topic_assignments2[sampled_token] -1)) +=1;
                        }
                        else{
                            topic_absent_edge_counts(document_author,a,(token_topic_assignments2[sampled_token] -1)) +=1;
                        } 
                    }//end update if statement
                    
                }//end recipeints loop
                
            }//end of loop over edges for for current document
        }// end of loop over docuemnts for edge-topic assignemnt step
        
        
    }//end of big number of itterations loop for entire step. 
    
    //return something
    List to_return(6);
    to_return[0] = token_topic_assignment_list;
    to_return[1] = topic_present_edge_counts;
    to_return[2] = topic_absent_edge_counts;
    to_return[3] = token_type_topic_counts;
    to_return[4] = edge_topic_assignments;
    to_return[5] = number_of_documents;
    return to_return;
    
    
    
    ///////////// LS MODEL BEGIN ///////////////
    
    srand((unsigned)time(NULL));
    std::default_random_engine generator;
    
    IntegerVector arrayDims1 = tpec.attr("dim");
    arma::cube topic_present_edge_counts(tpec.begin(), arrayDims1[0], arrayDims1[1], arrayDims1[2], false);
    
    IntegerVector arrayDims2 = taec.attr("dim");
    arma::cube topic_absent_edge_counts(taec.begin(), arrayDims2[0], arrayDims2[1], arrayDims2[2], false);
    
    IntegerVector arrayDims3 = clp.attr("dim");
    arma::cube current_latent_positions(clp.begin(), arrayDims3[0], arrayDims3[1], arrayDims3[2], false);
    
    
    //arma::cube proposed_latent_positions(plp.begin(), arrayDims4[0], arrayDims4[1], arrayDims4[2], false);
    //Create an array to hold new latent positions
    //arma::cube proposed_latent_positions = current_latent_positions;
    
    //NumericVector proposed_intercepts = current_intercepts;
    
    //NumericMatrix proposed_betas = betas;
    
    IntegerVector arrayDims5 = indicator_array.attr("dim");
    arma::cube beta_indicator_array(indicator_array.begin(), arrayDims5[0], arrayDims5[1], arrayDims5[2], false);
    
    
    
    //this is what we return -- it must contain intercepts, betas, latent positions and whether accepted proposal for all iiterations.
    int list_length = (7*number_of_metropolis_itterations);
    List to_return(list_length);
    
    
    
    //loop over the number of metropolis itterations (default 1000)
    for(int i = 0; i < number_of_metropolis_itterations; ++i){
        
        double beta_val = 0;
        NumericVector current_author_position(number_of_latent_dimensions);
        NumericVector proposed_author_position(number_of_latent_dimensions);
        NumericVector recipient_position(number_of_latent_dimensions);
        double sum_log_probability_of_current_positions = 0;
        double sum_log_probability_of_proposed_positions = 0;
        NumericVector proposed_intercepts(number_of_topics);
        IntegerVector arrayDims4 = plp.attr("dim");
        arma::cube proposed_latent_positions(plp.begin(), arrayDims4[0], arrayDims4[1], arrayDims4[2], false);
        NumericMatrix proposed_betas(number_of_topics,number_of_betas);
        
        //calculate proposed intercepts,latent positions, betas  double x = Rf_rnorm(mean,st. dev);
        int t = current_cluster;
        //for intercepts
        std::normal_distribution<double> distribution1(current_intercepts[t],proposal_variance);
        proposed_intercepts[t] = distribution1(generator);
        //proposed_intercepts[t] = Rf_rnorm(current_intercepts[t],proposal_variance);
        //for latent positions
        for(int a = 0; a < number_of_actors; ++a){
            for(int l = 0; l < number_of_latent_dimensions; ++l){
                std::normal_distribution<double> distribution2(current_latent_positions(l,t,a),proposal_variance);
                proposed_latent_positions(l,t,a) = distribution2(generator);
                //proposed_latent_positions(l,t,a) = Rf_rnorm(current_latent_positions(l,t,a),proposal_variance);
            }
        }
        //for betas
        for(int b = 0; b < number_of_betas; ++b){
            std::normal_distribution<double> distribution3(betas(t,b),proposal_variance);
            proposed_betas(t,b) = distribution3(generator);
            //proposed_betas(t,b) = Rf_rnorm(betas(t,b),proposal_variance);
        }
        
        
        
        //main loop
        //for(int t = 0; t < 0; ++t){
        int t = current_cluster;
        //get current topic intercept
        double current_topic_intercept = current_intercepts[t];
        double proposed_topic_intercept = proposed_intercepts[t];
        //get current topic betas
        NumericVector current_topic_betas = betas.row(t);
        NumericVector proposed_topic_betas = proposed_betas.row(t);
        
        //for(int a = 0; a < 0; ++a){
        for(int a = 0; a < number_of_actors; ++a){
            
            //get current position
            current_author_position[0] = current_latent_positions(0,t,a);
            current_author_position[1] = current_latent_positions(1,t,a);
            proposed_author_position[0] = proposed_latent_positions(0,t,a);
            proposed_author_position[1] = proposed_latent_positions(1,t,a);
            
            for(int b = 0; b < number_of_actors; ++b){
                if(b!= a){
                    
                    
                    //get number of actual edge present and absent for this topic sender reciever combo
                    int num_actual_edge = topic_present_edge_counts(a,b,t) + topic_present_edge_counts(b,a,t);
                    int num_non_edge = topic_absent_edge_counts(a,b,t) + topic_absent_edge_counts(b,a,t) ;
                    
                    // ===== calculate likelihood for current position =====//
                        //get current recipient position
                    recipient_position[0] = current_latent_positions(0,t,b);
                    recipient_position[1] = current_latent_positions(1,t,b);
                    
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
                    
                    //calculate likelihoods for both
                    double log_prob_edge = 0;
                    double log_prob_no_edge = 0;
                    if (eta != 0){
                        if(eta < 0){
                            log_prob_edge = eta -log(1 + exp(eta));
                            log_prob_no_edge = 0 -log(1 + exp(eta));
                        }
                        else{
                            log_prob_edge = 0 -log(1 + exp(-eta));
                            log_prob_no_edge = 0 -eta -log(1 + exp(-eta));
                        }
                    }
                    
                    //multiply and add to sum
                    sum_log_probability_of_current_positions += num_actual_edge*log_prob_edge;
                    sum_log_probability_of_current_positions += num_non_edge*log_prob_no_edge;
                    
                    
                    // ======== Now calculate for new positions ==========//
                        //get current recipient position
                    recipient_position[0] = proposed_latent_positions(0,t,b);
                    recipient_position[1] = proposed_latent_positions(1,t,b);
                    
                    //initialize distance
                    distance = 0;
                    //calculate distance
                    for(int k = 0; k < number_of_latent_dimensions; ++k){
                        distance += pow((proposed_author_position[k] - recipient_position[k]),2);
                    }
                    
                    beta_val = 0;
                    for(int c = 0; c < number_of_betas; ++c){
                        beta_val += proposed_topic_betas[c]*beta_indicator_array(a,b,c);
                    }
                    //calculate linear predictor
                    eta = proposed_topic_intercept - pow(distance,.5) + beta_val;
                    //calculate likelihoods for both
                    log_prob_edge = 0;
                    log_prob_no_edge = 0;
                    if (eta != 0){
                        if(eta < 0){
                            log_prob_edge = eta -log(1 + exp(eta));
                            log_prob_no_edge = 0 -log(1 + exp(eta));
                        }
                        else{
                            log_prob_edge = 0 -log(1 + exp(-eta));
                            log_prob_no_edge = 0 -eta -log(1 + exp(-eta));
                        }
                    }
                    
                    //multiply and add to sum
                    sum_log_probability_of_proposed_positions += num_actual_edge*log_prob_edge;
                    sum_log_probability_of_proposed_positions += num_non_edge*log_prob_no_edge;
                }
                
                
            }
        }
        
        
        
        //now calculate log ratio between two
        double log_ratio = sum_log_probability_of_proposed_positions - sum_log_probability_of_current_positions;
        
        double rand_num=((double)rand()/(double)RAND_MAX);
        double lud = log(rand_num);
        
        //take the log of a uniform draw on 0 to  1
        //double lud = as<double>(log_uniform_draw());
        
        if(log_ratio < lud){
            //if the log ratio is smaller then reject the new positions
            to_return[i] = current_latent_positions; 
            to_return[number_of_metropolis_itterations+i] = current_intercepts;
            to_return[2*number_of_metropolis_itterations+i] = betas;
            to_return[3*number_of_metropolis_itterations+i] = 0;
            to_return[4*number_of_metropolis_itterations+i] = sum_log_probability_of_proposed_positions;
            to_return[5*number_of_metropolis_itterations+i] = sum_log_probability_of_current_positions;
            to_return[6*number_of_metropolis_itterations+i] = lud;
            
        }
        else{
            //accept the new positions 
            to_return[i] = proposed_latent_positions; 
            to_return[number_of_metropolis_itterations+i] = proposed_intercepts;
            to_return[2*number_of_metropolis_itterations+i] = proposed_betas;
            to_return[3*number_of_metropolis_itterations+i] = 1;
            to_return[4*number_of_metropolis_itterations+i] = sum_log_probability_of_proposed_positions;
            to_return[5*number_of_metropolis_itterations+i] = sum_log_probability_of_current_positions;
            to_return[6*number_of_metropolis_itterations+i] = lud;
            //update current data structures with proposed positions
            
            
            current_latent_positions = proposed_latent_positions;
            current_intercepts = proposed_intercepts;
            betas = proposed_betas;
        }
    }

    return to_return;
}

//Metropolis_Step_CPP(30,100,array(1:90000,c(30,30,100)),array(1:90000,c(30,30,100)),array(runif(6000),c(2,100,30)),c(1:100),2,array(runif(6000),c(2,100,30)),c(1:100))
            

