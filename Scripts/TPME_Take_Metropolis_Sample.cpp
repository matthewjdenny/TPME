#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

    

// [[Rcpp::export]]
List Metropolis_Sample_CPP(
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
    int number_of_metropolis_itterations,
    double proposal_variance,
    NumericVector plp,
    int sample_interval
    ){
        
    Function log_uniform_draw("log_uniform_draw");
    
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
    int divided_itterations = number_of_metropolis_itterations/sample_interval;
    int list_length = (7*divided_itterations);
    List to_return(list_length);
    
    int take_draw = 0;
    int sample_number = 0;
    
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
        for(int t = 0; t < number_of_topics; ++t){
            //for intercepts
            proposed_intercepts[t] = Rf_rnorm(current_intercepts[t],proposal_variance);
            //for latent positions
            for(int a = 0; a < number_of_actors; ++a){
                for(int l = 0; l < number_of_latent_dimensions; ++l){
                    proposed_latent_positions(l,t,a) = Rf_rnorm(current_latent_positions(l,t,a),proposal_variance);
                }
            }
            //for betas
            for(int b = 0; b < number_of_betas; ++b){
                proposed_betas(t,b) = Rf_rnorm(betas(t,b),proposal_variance);
            }
        }
        
        
        //main loop
        //for(int t = 0; t < 0; ++t){
        for(int t = 0; t < number_of_topics; ++t){
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
        }
        
        
        //now calculate log ratio between two
        double log_ratio = sum_log_probability_of_proposed_positions - sum_log_probability_of_current_positions;
        
        //take the log of a uniform draw on 0 to  1
        double lud = as<double>(log_uniform_draw());
        lud -= 0.1;
        
        take_draw += 1;
        if(log_ratio < (0-1.4) ){
            //if the log ratio is smaller then reject the new positions
            if(take_draw == sample_interval){
            to_return[sample_number] = current_latent_positions; 
            to_return[divided_itterations+sample_number] = current_intercepts;
            to_return[2*divided_itterations+sample_number] = betas;
            to_return[3*divided_itterations+sample_number] = 0;
            to_return[4*divided_itterations+sample_number] = sum_log_probability_of_proposed_positions;
            to_return[5*divided_itterations+sample_number] = sum_log_probability_of_current_positions;
            to_return[6*divided_itterations+sample_number] = lud;
            take_draw = 0;
            sample_number += 1;
            }
            
            
        }
        else{
            //accept the new positions 
            if(take_draw == sample_interval){
            to_return[sample_number] = proposed_latent_positions; 
            to_return[divided_itterations+sample_number] = proposed_intercepts;
            to_return[2*divided_itterations+sample_number] = proposed_betas;
            to_return[3*divided_itterations+sample_number] = 1;
            to_return[4*divided_itterations+sample_number] = sum_log_probability_of_proposed_positions;
            to_return[5*divided_itterations+sample_number] =sum_log_probability_of_current_positions;
            to_return[6*divided_itterations+sample_number] = lud;
            take_draw = 0;
            sample_number += 1;
            }
            //update current data structures with proposed positions
            
            
            current_latent_positions = proposed_latent_positions;
            current_intercepts = proposed_intercepts;
            betas = proposed_betas;
        }
    }

    return to_return;
}

//Metropolis_Step_CPP(30,100,array(1:90000,c(30,30,100)),array(1:90000,c(30,30,100)),array(runif(6000),c(2,100,30)),c(1:100),2,array(runif(6000),c(2,100,30)),c(1:100))
            

