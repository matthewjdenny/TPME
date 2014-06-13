#include <RcppArmadillo.h>
#include <random>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
List Cluster_Integrated_Final_MH_Sampler(
    int number_of_MH_itterations,
    int number_of_actors, 
    int number_of_topics,
    int number_of_latent_dimensions,
    arma::vec proposal_variance,
    arma::vec topic_cluster_assignments,
    NumericVector tpec,
    NumericVector taec,
    NumericVector clp,
    NumericVector plp, 
    arma::vec current_intercepts,
    arma::mat betas,
    int number_of_betas,
    NumericVector indicator_array,
    arma::mat edge_topic_assignments,
    int burnin,
    int number_of_clusters,
    int current_cluster,
    int store_every_x_rounds,
    double post_burnin_variance_multiplier
    ){
    
    Function report("Report_Probs");
    //Function report2("Report_2");
    
    //count up to a storage round
    int storage_counter = 0;

    //one less than the number of unique objects I want to store at the beginning of the list
    int list_length = 8;
    List to_return(list_length);
    
    int MH_Counter = 0;
    arma::vec accept(number_of_MH_itterations/store_every_x_rounds);
    arma::vec Proposed_MH_Likelihoods(number_of_MH_itterations/store_every_x_rounds);
    arma::vec Current_MH_Likelihoods(number_of_MH_itterations/store_every_x_rounds);
    arma::vec store_intercepts(number_of_MH_itterations/store_every_x_rounds);
    arma::mat store_betas(number_of_betas,number_of_MH_itterations/store_every_x_rounds);
    arma::cube store_latent_positions(number_of_latent_dimensions,number_of_actors,number_of_MH_itterations/store_every_x_rounds);
    
    //set seed and designate RNG
    srand(1234);
    std::default_random_engine generator;
    
    
    //read in topic present edge counts array [num actors x num actors x topics]
    IntegerVector arrayDims1 = tpec.attr("dim");
    arma::cube topic_present_edge_counts(tpec.begin(), arrayDims1[0], arrayDims1[1], arrayDims1[2], false);
    
    //read in topic absent edge counts array [num actors x num actors x topics]
    IntegerVector arrayDims2 = taec.attr("dim");
    arma::cube topic_absent_edge_counts(taec.begin(), arrayDims2[0], arrayDims2[1], arrayDims2[2], false);
    
    //read in latent positions array [num dimensions x clusters x actors]
    IntegerVector arrayDims3 = clp.attr("dim");
    arma::cube current_latent_positions(clp.begin(), arrayDims3[0], arrayDims3[1], arrayDims3[2], false);
    
    //read in beta indicator array (0,1) [number of topics x number of actors x number of betas]
    IntegerVector arrayDims5 = indicator_array.attr("dim"); 
    arma::cube beta_indicator_array(indicator_array.begin(), arrayDims5[0], arrayDims5[1], arrayDims5[2], false);
    
    //double beta_val = 0;
    arma::vec current_author_position(number_of_latent_dimensions);
    arma::vec proposed_author_position(number_of_latent_dimensions);
    arma::vec recipient_position(number_of_latent_dimensions);
    //NumericVector proposed_intercepts;
    IntegerVector arrayDims4 = plp.attr("dim");
    arma::cube proposed_latent_positions(plp.begin(), arrayDims4[0], arrayDims4[1], arrayDims4[2], false);
    
    arma::mat proposed_betas(number_of_clusters,number_of_betas);    
    double standard_deviation = 2;
    double dist_center = 0;
    arma::vec current_cluster_betas(number_of_betas);
    arma::vec proposed_cluster_betas(number_of_betas);
    //arma::vec bets(number_of_betas);
    //arma::mat lat_pos(number_of_latent_dimensions,number_of_actors);
    
    
    //arma::vec proposed_int(number_of_clusters);
    
    
    //this sets k for the entire time
    int k = current_cluster - 1;
    
    //create a double since we are only needing one for the durration
    double proposed_intercept = 0;
    double current_intercept = current_intercepts[k];
    
    //we are holding out the MM homophily parameter in this case as our baseline
    for(int k = 0; k < number_of_clusters; ++k){
        betas(k,0) = 0;
        proposed_betas(k,0) = 0;
    }
    
    for(int i = 0; i < number_of_MH_itterations; ++i){

            if(i == burnin){
                double* tempvar;
                tempvar = new double;
                *tempvar = double(proposal_variance[k]);
                proposal_variance[k] = post_burnin_variance_multiplier*(*tempvar);
                delete tempvar;
            }
            //report(i);
            //calculate proposed intercepts,latent positions, betas  double x = Rf_rnorm(mean,st. dev);
            std::normal_distribution<double> distribution1(current_intercept,proposal_variance[k]);
            proposed_intercept = double(distribution1(generator));
            
            for(int a = 0; a < number_of_actors; ++a){
                for(int l = 0; l < number_of_latent_dimensions; ++l){
                    std::normal_distribution<double> distribution2(current_latent_positions(l,k,a),proposal_variance[k]);
                    proposed_latent_positions(l,k,a) = distribution2(generator);
                }
            }
            
            for(int b = 1; b < number_of_betas; ++b){
                std::normal_distribution<double> distribution3(betas(k,b),proposal_variance[k]);
                proposed_betas(k,b) = distribution3(generator);
            }
            
            //main loop
            double lsm_sum_log_probability_of_current_positions = 0 ;
            double lsm_sum_log_probability_of_proposed_positions = 0 ;
            
            for(int b = 1; b < number_of_betas; ++b){
                current_cluster_betas(b) = betas(k,b);
                proposed_cluster_betas(b) = proposed_betas(k,b);
            }
            
            
            double lsm_prior_current_positions =0;
            double lsm_prior_proposed_positions = 0;
            
            lsm_prior_current_positions += log(( 1 / ( standard_deviation * sqrt(2*M_PI) ) ) * exp( -0.5 * pow( (current_intercept-dist_center)/standard_deviation, 2.0 ) ));
            lsm_prior_proposed_positions += log(( 1 / ( standard_deviation * sqrt(2*M_PI) ) ) * exp( -0.5 * pow( (proposed_intercept-dist_center)/standard_deviation, 2.0 ) ));
            for(int c = 0; c < number_of_betas; ++c){
                lsm_prior_current_positions += log(( 1 / ( standard_deviation * sqrt(2*M_PI) ) ) * exp( -0.5 * pow( (current_cluster_betas[c]-dist_center)/standard_deviation, 2.0 ) ));
                lsm_prior_proposed_positions += log(( 1 / ( standard_deviation * sqrt(2*M_PI) ) ) * exp( -0.5 * pow( (proposed_cluster_betas[c]-dist_center)/standard_deviation, 2.0 ) ));
            }
            
            for(int a = 0; a < number_of_actors; ++a){
                
                //get current position
                current_author_position[0] = current_latent_positions(0,k,a);
                current_author_position[1] = current_latent_positions(1,k,a);
                proposed_author_position[0] = proposed_latent_positions(0,k,a);
                proposed_author_position[1] = proposed_latent_positions(1,k,a);
                //add of for each latent position
                for(int c = 0; c < number_of_latent_dimensions; ++c){
                    lsm_prior_current_positions += log(( 1 / ( standard_deviation * sqrt(2*M_PI) ) ) * exp( -0.5 * pow( (current_author_position[c]-dist_center)/standard_deviation, 2.0 ) ));
                    lsm_prior_proposed_positions += log(( 1 / ( standard_deviation * sqrt(2*M_PI) ) ) * exp( -0.5 * pow( (proposed_author_position[c]-dist_center)/standard_deviation, 2.0 ) ));
                }
            
                for(int b = 0; b < number_of_actors; ++b){
                    if(b!= a){
                        int num_actual_edge =0;
                        int num_non_edge =0;
                        //get number of actual edge present and absent for this cluster sender reciever combo
                        for(int t = 0; t < number_of_topics; ++t){
                            double topic_cluster = topic_cluster_assignments[t] -1;
                            if(topic_cluster == k){
                                num_actual_edge += topic_present_edge_counts(a,b,t); //+ topic_present_edge_counts(b,a,t);
                                num_non_edge += topic_absent_edge_counts(a,b,t); //+ topic_absent_edge_counts(b,a,t) ;
                            }
                        }
                        

                        // ===== calculate likelihood for current position =====//
                            //get current recipient position
                        recipient_position[0] = current_latent_positions(0,k,b);
                        recipient_position[1] = current_latent_positions(1,k,b);

                        //initialize distance
                        double distance =0;


                        //calculate distance
                        for(int j = 0; j < number_of_latent_dimensions; ++j){
                            distance += pow((current_author_position[j] - recipient_position[j]),2);
                        }

                        double beta_val =0;


                        for(int c = 0; c < number_of_betas; ++c){
                            beta_val += current_cluster_betas[c]*beta_indicator_array(a,b,c);
                        }
                        

                        //calculate linear predictor
                        double eta =0;

                        eta = current_intercept - pow(distance,.5) + beta_val;


                        //calculate likelihoods for both
                        double log_prob_edge =0;
                        double log_prob_no_edge =0;
                        
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
                        lsm_sum_log_probability_of_current_positions += num_actual_edge*log_prob_edge;
                        
                        lsm_sum_log_probability_of_current_positions += num_non_edge*log_prob_no_edge;
                        
                        
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
                        
                        beta_val= 0;
                        for(int c = 0; c < number_of_betas; ++c){
                            beta_val += proposed_cluster_betas[c]*beta_indicator_array(a,b,c);
                        }
                        
                        //calculate linear predictor

                        eta = proposed_intercept - pow(distance,.5) + beta_val;

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
                        lsm_sum_log_probability_of_proposed_positions += num_actual_edge*log_prob_edge;
                        lsm_sum_log_probability_of_proposed_positions += num_non_edge*log_prob_no_edge;
                        
                    }   
                }
            }// end of actor a loop
            
            //now calculate log ratio between two
            lsm_sum_log_probability_of_proposed_positions += lsm_prior_proposed_positions;
            
            lsm_sum_log_probability_of_current_positions += lsm_prior_current_positions;
            double accepted = 0;
            
            double log_ratio =0;
            log_ratio= lsm_sum_log_probability_of_proposed_positions - lsm_sum_log_probability_of_current_positions;
            double rand_num=((double)rand()/(double)RAND_MAX);
            double lud =0;
            lud= log(rand_num);
            if(log_ratio == 0){
                log_ratio = -100;
            }
            if(log_ratio < lud){
                //if the log ratio is smaller then reject the new positions
                // DO NOTHING
                accepted = 0;
            }
            else{
                //report(12345);
                accepted = 1;
                current_intercept = proposed_intercept;
                //for latent positions
                for(int a = 0; a < number_of_actors; ++a){
                    for(int l = 0; l < number_of_latent_dimensions; ++l){
                        current_latent_positions(l,k,a) = proposed_latent_positions(l,k,a);
                    }
                }
                //for betas
                for(int b = 0; b < number_of_betas; ++b){
                    betas(k,b) = proposed_betas(k,b);
                }
            }//end of update 
            
            storage_counter +=1;
            if(store_every_x_rounds == storage_counter){
                report(MH_Counter);
                //arma::vec bets(number_of_betas);
                //arma::mat lat_pos(number_of_latent_dimensions,number_of_actors);
                storage_counter =0;
                double ints =0;
                ints = current_intercept;
                for(int b = 0; b < number_of_betas; ++b){
                    store_betas(b,MH_Counter) = betas(k,b);
                    //store_betas(k,b,MH_Counter) = 1;
                }
                
                for(int a = 0; a < number_of_actors; ++a){
                    for(int l = 0; l < number_of_latent_dimensions; ++l){
                        store_latent_positions(l,a,MH_Counter) = current_latent_positions(l,k,a);
                        //store_latent_positions(l,a,MH_Counter) = 1;
                    }
                }
                //store everything
                accept[MH_Counter] = accepted;
                Proposed_MH_Likelihoods[MH_Counter] = lsm_sum_log_probability_of_proposed_positions;
                Current_MH_Likelihoods[MH_Counter] = lsm_sum_log_probability_of_current_positions;
                store_intercepts[MH_Counter] = ints;
                //store_intercepts[MH_Counter] = 1;
                MH_Counter += 1;
            }//end of save 
            
    }//end of outer i loop
   
    to_return[0] = current_cluster;
    to_return[1] = number_of_MH_itterations;
    to_return[2] = Proposed_MH_Likelihoods;
    to_return[3] = Current_MH_Likelihoods;
    to_return[4] = accept;
    to_return[5] = store_intercepts;
    to_return[6] = store_latent_positions;
    to_return[7] = store_betas;

 return to_return;
}

            

