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
    NumericVector plp,
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
    //Function log_multinomial_draw("log_multinomial_draw");
    
    Function report("Report_Probs");
    
    int list_length = (10 + number_of_outer_itterations*number_of_Gibbs_itterations + 3*number_of_MH_itterations);
    List to_return(list_length);
    int Gibbs_Counter = 1;
    int MH_Counter = 1;
    
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

    //outer loop over the number of itterations (default 1000)
    for(int n = 0; n < number_of_outer_itterations; ++n){
        
        report(topic_cluster_assignments);
    
        // ===================================================== //
        // ===================== LDA Step ====================== //
        // ===================================================== //
    
        for(int i = 0; i < number_of_Gibbs_itterations; ++i){
            
            
            //initialize variables that will be used across iterations
            double beta_val = 0;
            NumericVector current_author_position(number_of_latent_dimensions);
            NumericVector recipient_position(number_of_latent_dimensions);
            NumericVector current_topic_betas(number_of_betas);
            
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
                        
                        //get the current topic's cluster assignment
                        double cluster = topic_cluster_assignments[t] - 1;
                        
                        //get current author latent positions (could be updated to run in a loop)
                        current_author_position[0] = current_latent_positions(0,cluster,document_author);
                        current_author_position[1] = current_latent_positions(1,cluster,document_author);
                        
                        //get current topic betas and intercepts
                        double current_cluster_intercept = current_intercepts[cluster];
                        NumericVector current_cluster_betas = betas.row(cluster);
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
                                    recipient_position[0] = current_latent_positions(0,cluster,a);
                                    recipient_position[1] = current_latent_positions(1,cluster,a);
                                    
                                    //initialize distance
                                    double distance = 0;
                                    //calculate distance
                                    for(int k = 0; k < number_of_latent_dimensions; ++k){
                                        distance += pow((current_author_position[k] - recipient_position[k]),2);
                                    }
                            
                                    beta_val = 0;
                                    //pull out the correct beta value 
                                    for(int c = 0; c < number_of_betas; ++c){
                                        beta_val += current_cluster_betas[c]*beta_indicator_array(document_author,a,c);
                                    }
                                    //calculate linear predictor
                                    double eta = current_cluster_intercept - pow(distance,.5) + beta_val;
                            
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
            
            
            
            // ============= edge topic assignment step ============== // 
            
            
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
                            
                            //get the current topic's cluster assignment
                            double cluster = topic_cluster_assignments[topic_assignment] -1;
                        
                            //get current topic intercept
                            double current_cluster_intercept = current_intercepts[cluster];
                            //get current topic betas
                            NumericVector current_cluster_betas = betas.row(cluster);
        
                            //get current author position
                            current_author_position[0] = current_latent_positions(0,cluster,document_author);
                            current_author_position[1] = current_latent_positions(1,cluster,document_author);
                            //get current recipient position
                            recipient_position[0] = current_latent_positions(0,cluster,a);
                            recipient_position[1] = current_latent_positions(1,cluster,a);
                            
                            //initialize distance
                            double distance = 0;
                            //calculate distance
                            for(int k = 0; k < number_of_latent_dimensions; ++k){
                                distance += pow((current_author_position[k] - recipient_position[k]),2);
                            }
                            
                            beta_val = 0;
                            for(int c = 0; c < number_of_betas; ++c){
                                beta_val += current_cluster_betas[c]*beta_indicator_array(document_author,a,c);
                            }
                            //calculate linear predictor
                            double eta = current_cluster_intercept - pow(distance,.5) + beta_val;
                            
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
                        
                        // do some updating
                        if(actual_edge == 1){
                            topic_present_edge_counts(document_author,a,(token_topic_assignments2[sampled_token] -1)) +=1;
                        }
                        else{
                            topic_absent_edge_counts(document_author,a,(token_topic_assignments2[sampled_token] -1)) +=1;
                        } 
                    }//end recipeints loop
                }//end of loop over edges for for current document
            }// end of loop over docuemnts for edge-topic assignemnt step
            
            
            // ===================== Hard Cluster Topics ====================== //
            
            
            for(int t = 0; t < number_of_topics; ++t){
                
                NumericVector topic_cluster_distribution(number_of_clusters);
                
                for(int k = 0; k < number_of_clusters; ++k){
                    
                    NumericVector current_author_position(number_of_latent_dimensions);
                    NumericVector recipient_position(number_of_latent_dimensions);
                
                    //get current cluster intercept
                    double current_cluster_intercept = current_intercepts[k];
                    //get current cluster betas
                    NumericVector current_cluster_betas = betas.row(k);
                    double sum_log_probability_of_current_positions = 0;
                    
                    //for(int a = 0; a < 0; ++a){
                    for(int a = 0; a < number_of_actors; ++a){
                        
                        //get current position
                        current_author_position[0] = current_latent_positions(0,k,a);
                        current_author_position[1] = current_latent_positions(1,k,a);
                        
                        for(int b = 0; b < number_of_actors; ++b){
                            if(b!= a){
                                
                                
                                //get number of actual edge present and absent for this topic sender reciever combo
                                int num_actual_edge = topic_present_edge_counts(a,b,t) + topic_present_edge_counts(b,a,t);
                                int num_non_edge = topic_absent_edge_counts(a,b,t) + topic_absent_edge_counts(b,a,t) ;
                                
                                // ===== calculate likelihood for current position =====//
                                    //get current recipient position
                                recipient_position[0] = current_latent_positions(0,k,b);
                                recipient_position[1] = current_latent_positions(1,k,b);
                                
                                //initialize distance
                                double distance = 0;
                                //calculate distance
                                for(int j = 0; j < number_of_latent_dimensions; ++j){
                                    distance += pow((current_author_position[j] - recipient_position[j]),2);
                                }
                                
                                double beta_val = 0;
                                for(int c = 0; c < number_of_betas; ++c){
                                    beta_val += current_cluster_betas[c]*beta_indicator_array(a,b,c);
                                }
                                //calculate linear predictor
                                double eta = current_cluster_intercept - pow(distance,.5) + beta_val;
                                
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
                                 
                            }
                        }
                    }
                    //assign this probability to the cluster probability vector
                    topic_cluster_distribution[k] = sum_log_probability_of_current_positions;
                    
                } // end loop over clusters
                report(topic_cluster_distribution);
                NumericVector cluster_probabilities(number_of_clusters);
                for(int x = 0; x < number_of_clusters; ++x){
                    cluster_probabilities[x] = exp(topic_cluster_distribution[x]);
                }
                
                std::discrete_distribution<int> distribution_clusters (cluster_probabilities.begin(),cluster_probabilities.end());
                int cluster_assignment = distribution_clusters(generator) + 1;
                //report(cluster_probabilities);
                //set the new cluster assignment for the topic 
                topic_cluster_assignments[t] = double(cluster_assignment);
            }//end loop over topics
            
            to_return[10+Gibbs_Counter] = topic_cluster_assignments;
            Gibbs_Counter += 1;
        }//end of big number of itterations loop for entire step. 
        
        //return something
        
        
        
        // ===================================================== //
        // ===================== LSM Step ====================== //
        // ===================================================== //
        
        //loop over the number of metropolis itterations (default 1000)
        for(int i = 0; i < number_of_MH_itterations; ++i){
            
            double beta_val = 0;
            NumericVector current_author_position(number_of_latent_dimensions);
            NumericVector proposed_author_position(number_of_latent_dimensions);
            NumericVector recipient_position(number_of_latent_dimensions);
            NumericVector proposed_intercepts(number_of_clusters);
            IntegerVector arrayDims4 = plp.attr("dim");
            arma::cube proposed_latent_positions(plp.begin(), arrayDims4[0], arrayDims4[1], arrayDims4[2], false);
            NumericMatrix proposed_betas(number_of_clusters,number_of_betas);
            
            //calculate proposed intercepts,latent positions, betas  double x = Rf_rnorm(mean,st. dev);
            for(int k = 0; k < number_of_clusters; ++k){
                //for intercepts
                std::normal_distribution<double> distribution1(current_intercepts[k],proposal_variance);
                proposed_intercepts[k] = distribution1(generator);
                //for latent positions
                for(int a = 0; a < number_of_actors; ++a){
                    for(int l = 0; l < number_of_latent_dimensions; ++l){
                        std::normal_distribution<double> distribution2(current_latent_positions(l,k,a),proposal_variance);
                        proposed_latent_positions(l,k,a) = distribution2(generator);
                    }
                }
                //for betas
                for(int b = 0; b < number_of_betas; ++b){
                    std::normal_distribution<double> distribution3(betas(k,b),proposal_variance);
                    proposed_betas(k,b) = distribution3(generator);
                }
            }//end of loop over generating new potenttial LS positions
            
            
            //main loop
            for(int k = 0; k < number_of_clusters; ++k){
                double sum_log_probability_of_current_positions = 0;
                double sum_log_probability_of_proposed_positions = 0;
                //get current topic intercept
                double current_cluster_intercept = current_intercepts[k];
                double proposed_cluster_intercept = proposed_intercepts[k];
                //get current topic betas
                NumericVector current_cluster_betas = betas.row(k);
                NumericVector proposed_cluster_betas = proposed_betas.row(k);
                
                //for(int a = 0; a < 0; ++a){
                for(int a = 0; a < number_of_actors; ++a){
                    
                    //get current position
                    current_author_position[0] = current_latent_positions(0,k,a);
                    current_author_position[1] = current_latent_positions(1,k,a);
                    proposed_author_position[0] = proposed_latent_positions(0,k,a);
                    proposed_author_position[1] = proposed_latent_positions(1,k,a);
                    
                    for(int b = 0; b < number_of_actors; ++b){
                        if(b!= a){
                            
                            int num_actual_edge = 0;
                            int num_non_edge = 0;
                            //get number of actual edge present and absent for this cluster sender reciever combo
                            for(int t = 0; t < number_of_topics; ++t){
                                double topic_cluster = topic_cluster_assignments[t] -1;
                                if(topic_cluster == k){
                                    num_actual_edge += topic_present_edge_counts(a,b,t) + topic_present_edge_counts(b,a,t);
                                    num_non_edge += topic_absent_edge_counts(a,b,t) + topic_absent_edge_counts(b,a,t) ;
                                }
                            }
                            
                            
                            // ===== calculate likelihood for current position =====//
                                //get current recipient position
                            recipient_position[0] = current_latent_positions(0,k,b);
                            recipient_position[1] = current_latent_positions(1,k,b);
                            
                            //initialize distance
                            double distance = 0;
                            //calculate distance
                            for(int j = 0; j < number_of_latent_dimensions; ++j){
                                distance += pow((current_author_position[j] - recipient_position[j]),2);
                            }
                            
                            beta_val = 0;
                            for(int c = 0; c < number_of_betas; ++c){
                                beta_val += current_cluster_betas[c]*beta_indicator_array(a,b,c);
                            }
                            //calculate linear predictor
                            double eta = current_cluster_intercept - pow(distance,.5) + beta_val;
                            
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
                            recipient_position[0] = proposed_latent_positions(0,k,b);
                            recipient_position[1] = proposed_latent_positions(1,k,b);
                            
                            //initialize distance
                            distance = 0;
                            //calculate distance
                            for(int j = 0; j < number_of_latent_dimensions; ++j){
                                distance += pow((proposed_author_position[j] - recipient_position[j]),2);
                            }
                            
                            beta_val = 0;
                            for(int c = 0; c < number_of_betas; ++c){
                                beta_val += proposed_cluster_betas[c]*beta_indicator_array(a,b,c);
                            }
                            //calculate linear predictor
                            eta = proposed_cluster_intercept - pow(distance,.5) + beta_val;
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
                    // DO NOTHING
                }
                else{
                    //if we accept UPDATE
                    //for intercepts
                    current_intercepts[k] = proposed_intercepts[k];
                    //for latent positions
                    for(int a = 0; a < number_of_actors; ++a){
                        for(int l = 0; l < number_of_latent_dimensions; ++l){
                            current_latent_positions(l,k,a) = proposed_latent_positions(l,k,a);
                        }
                    }
                    //for betas
                    for(int b = 0; b < number_of_betas; ++b){
                        betas(k,b) = proposed_betas(k,b) ;
                    }
                }
            }//end of loop over clusters for which we are doing separate proposals and acceptances
            
            //if we are in the last outter iteration save everything
            if(n == (number_of_outer_itterations -1)){
                to_return[10 + number_of_outer_itterations*number_of_Gibbs_itterations + MH_Counter] = current_intercepts;
                to_return[10 + number_of_outer_itterations*number_of_Gibbs_itterations + number_of_MH_itterations + MH_Counter] = current_latent_positions;
                to_return[10 + number_of_outer_itterations*number_of_Gibbs_itterations + 2*number_of_MH_itterations + MH_Counter] = betas;
                MH_Counter += 1;
            }
            
        }// end of MH loop

    }// end of main outer itteration loop 
    
    
    to_return[0] = token_topic_assignment_list;
    to_return[1] = topic_present_edge_counts;
    to_return[2] = topic_absent_edge_counts;
    to_return[3] = token_type_topic_counts;
    to_return[4] = edge_topic_assignments;
    to_return[5] = number_of_documents;
    to_return[6] = number_of_outer_itterations;
    to_return[7] = number_of_Gibbs_itterations;
    to_return[8] = number_of_MH_itterations;
    to_return[9] = number_of_clusters;
    return to_return;
}


            

