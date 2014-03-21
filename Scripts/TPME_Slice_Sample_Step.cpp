#include <RcppArmadillo.h>
#include <random>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
List Slice_Sample_Step_CPP(
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
    double step_size,
    NumericVector plp
    ){
 
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
    int list_length = (7*number_of_itterations);
    List to_return(list_length);
    
    
    
    //loop over the number of metropolis itterations (default 1000)
    for(int i = 0; i < number_of_itterations; ++i){
        
        NumericVector current_author_position(number_of_latent_dimensions);
        NumericVector proposed_author_position(number_of_latent_dimensions);
        NumericVector recipient_position(number_of_latent_dimensions);
        double sum_log_probability_of_current_positions = 0;
        double sum_log_probability_of_proposed_positions = 0;
        NumericVector proposed_intercepts(number_of_topics);
        IntegerVector arrayDims4 = plp.attr("dim");
        arma::cube proposed_latent_positions(plp.begin(), arrayDims4[0], arrayDims4[1], arrayDims4[2], false);
        NumericMatrix proposed_betas(number_of_topics,number_of_betas);
        
        
        // =================== Current Position Probability ==================== //
        //calculate the likelihood of the current position
        for(int t = 0; t < number_of_topics; ++t){
            //get current topic intercept
            double current_topic_intercept = current_intercepts[t];
            //get current topic betas
            NumericVector current_topic_betas = betas.row(t);
            for(int a = 0; a < number_of_actors; ++a){
                //get current position
                for(int k = 0; k < number_of_latent_dimensions; ++k){
                    current_author_position[k] = current_latent_positions(k,t,a);
                }
                for(int b = 0; b < number_of_actors; ++b){
                    if(b!= a){
                        //get number of actual edge present and absent for this topic sender reciever combo
                        int num_actual_edge = topic_present_edge_counts(a,b,t) + topic_present_edge_counts(b,a,t);
                        int num_non_edge = topic_absent_edge_counts(a,b,t) + topic_absent_edge_counts(b,a,t) ;
                        
                        //get recipient positions
                        for(int k = 0; k < number_of_latent_dimensions; ++k){
                            recipient_position[k] = current_latent_positions(k,t,b);
                        }
                        
                        //initialize distance
                        double distance = 0;
                        //calculate distance
                        for(int k = 0; k < number_of_latent_dimensions; ++k){
                            distance += pow((current_author_position[k] - recipient_position[k]),2);
                        }
                        
                        double beta_val = 0;
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
                    } 
                }
            }
        }
        
        //take a log uniform draw and add it to the probability of the current sample to get a floor on probabilites of new slice samples we can accept
        double rand_num=((double)rand()/(double)RAND_MAX);
        double lud = log(rand_num);
        double slice_probability_floor = sum_log_probability_of_current_positions + lud;
        
        
        //now go ahead and define hte slice for which we will be unifromly sampling new values from 
        //we need a left and right value
        NumericVector left_intercepts(number_of_topics);
        arma::cube left_latent_positions(plp.begin(), arrayDims4[0], arrayDims4[1], arrayDims4[2], false);
        NumericMatrix left_betas(number_of_topics,number_of_betas);
        NumericVector right_intercepts(number_of_topics);
        arma::cube right_latent_positions(plp.begin(), arrayDims4[0], arrayDims4[1], arrayDims4[2], false);
        NumericMatrix right_betas(number_of_topics,number_of_betas);
        
        //get the left and right bounds on the slice for intercepts,latent positions, betas  
        for(int t = 0; t < number_of_topics; ++t){
            //for intercepts
            rand_num=((double)rand()/(double)RAND_MAX);
            left_intercepts[t] = current_intercepts[t] - rand_num*step_size;
            right_intercepts[t] = left_intercepts[t] + step_size;
            for(int a = 0; a < number_of_actors; ++a){
                for(int l = 0; l < number_of_latent_dimensions; ++l){
                    rand_num=((double)rand()/(double)RAND_MAX);
                    left_latent_positions(l,t,a) = current_latent_positions(l,t,a) - rand_num*step_size;
                    right_latent_positions(l,t,a) = left_latent_positions(l,t,a) + step_size;
                }
            }
            //for betas
            for(int b = 0; b < number_of_betas; ++b){
                rand_num=((double)rand()/(double)RAND_MAX);
                left_betas(t,b) = betas(t,b) - rand_num*step_size;
                right_betas(t,b) = left_betas(t,b) + step_size;
            }
        }
        
        int in_slice = 0; 
        while(in_slice < 1){
            
            //get new values for the slice for intercepts,latent positions, betas  
            for(int t = 0; t < number_of_topics; ++t){
                //for intercepts
                rand_num=((double)rand()/(double)RAND_MAX);
                proposed_intercepts[t] =  left_intercepts[t] + rand_num*(right_intercepts[t] -left_intercepts[t]);
                for(int a = 0; a < number_of_actors; ++a){
                    for(int l = 0; l < number_of_latent_dimensions; ++l){
                        rand_num=((double)rand()/(double)RAND_MAX);
                        proposed_latent_positions(l,t,a) =  left_latent_positions(l,t,a) + rand_num*(right_latent_positions(l,t,a) - left_latent_positions(l,t,a));
                    }
                }
                //for betas
                for(int b = 0; b < number_of_betas; ++b){
                    rand_num=((double)rand()/(double)RAND_MAX);
                    proposed_betas(t,b) = left_betas(t,b) + rand_num*(right_betas(t,b) - left_betas(t,b));
                }
            }    
            
            
            
            
            
            
            // =================== Slice Sample Probability ==================== //
            //calculate the likelihood of the current draw from the slice
            for(int t = 0; t < number_of_topics; ++t){
                //get current topic intercept
                double proposed_topic_intercept = proposed_intercepts[t];
                //get current topic betas
                NumericVector proposed_topic_betas = proposed_betas.row(t);
                for(int a = 0; a < number_of_actors; ++a){
                    //get get the proposed positions
                    for(int k = 0; k < number_of_latent_dimensions; ++k){
                        proposed_author_position[k] = proposed_latent_positions(k,t,a);
                    }
                
                    for(int b = 0; b < number_of_actors; ++b){
                        if(b!= a){
                        
                            //get number of actual edge present and absent for this topic sender reciever combo
                            int num_actual_edge = topic_present_edge_counts(a,b,t) + topic_present_edge_counts(b,a,t);
                            int num_non_edge = topic_absent_edge_counts(a,b,t) + topic_absent_edge_counts(b,a,t) ;
                        
                            //get recipient positions
                            for(int k = 0; k < number_of_latent_dimensions; ++k){
                                recipient_position[k] = current_latent_positions(k,t,b);
                            }
                        
                            //initialize distance
                            double distance = 0;
                            //calculate distance
                            for(int k = 0; k < number_of_latent_dimensions; ++k){
                                distance += pow((proposed_author_position[k] - recipient_position[k]),2);
                            }
                        
                            double beta_val = 0;
                            for(int c = 0; c < number_of_betas; ++c){
                                beta_val += proposed_topic_betas[c]*beta_indicator_array(a,b,c);
                            }
                            //calculate linear predictor
                            double eta = proposed_topic_intercept - pow(distance,.5) + beta_val;
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
                            sum_log_probability_of_proposed_positions += num_actual_edge*log_prob_edge;
                            sum_log_probability_of_proposed_positions += num_non_edge*log_prob_no_edge;
                        }
                    }
                }
            }
            
            
            if(sum_log_probability_of_proposed_positions > slice_probability_floor){
                in_slice = 1;
            }
            else{
                //if the positions we tried were outside of the slice, set them as the new boundary
                //get the left and right bounds on the slice for intercepts,latent positions, betas  
                for(int t = 0; t < number_of_topics; ++t){
                    //for intercepts
                    if(proposed_intercepts[t] < current_intercepts[t]){
                        left_intercepts[t] = proposed_intercepts[t];
                    }
                    else{
                        right_intercepts[t] = proposed_intercepts[t];
                    }
                    
                    for(int a = 0; a < number_of_actors; ++a){
                        for(int l = 0; l < number_of_latent_dimensions; ++l){
                            if(proposed_latent_positions(l,t,a) < current_latent_positions(l,t,a)){
                                left_latent_positions(l,t,a)  = proposed_latent_positions(l,t,a);   
                            }
                            else{
                                right_latent_positions(l,t,a) = proposed_latent_positions(l,t,a); 
                            }
                        }
                    }
                    //for betas
                    for(int b = 0; b < number_of_betas; ++b){
                        if(proposed_betas(t,b) < betas(t,b)){
                            left_betas(t,b) = proposed_betas(t,b);
                        }
                        else{
                            right_betas(t,b) = proposed_betas(t,b);
                        }
                    }
                }
                
            }
        }// end while loop over checking to see if hte current positions are in the slice
        
        
        //update all values
        to_return[i] = proposed_latent_positions; 
        to_return[number_of_itterations+i] = proposed_intercepts;
        to_return[2*number_of_itterations+i] = proposed_betas;
        to_return[3*number_of_itterations+i] = 1;
        to_return[4*number_of_itterations+i] = sum_log_probability_of_proposed_positions;
        to_return[5*number_of_itterations+i] = sum_log_probability_of_current_positions;
        to_return[6*number_of_itterations+i] = lud;
        //update current data structures with proposed positions
        current_latent_positions = proposed_latent_positions;
        current_intercepts = proposed_intercepts;
        betas = proposed_betas;
        
        
    }

    return to_return;
}

//Metropolis_Step_CPP(30,100,array(1:90000,c(30,30,100)),array(1:90000,c(30,30,100)),array(runif(6000),c(2,100,30)),c(1:100),2,array(runif(6000),c(2,100,30)),c(1:100))
            


