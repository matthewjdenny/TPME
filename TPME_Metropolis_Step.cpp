#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
List Metropolis_Step_CPP(
    int number_of_actors, 
    int number_of_topics,
    NumericVector tpec,
    NumericVector taec,
    NumericVector clp, 
    NumericVector current_intercepts,
    int number_of_latent_dimensions,
    NumericVector plp, 
    NumericVector proposed_intercepts
    ){
        
    Function log_uniform_draw("log_uniform_draw");
    
    IntegerVector arrayDims1 = tpec.attr("dim");
    arma::cube topic_present_edge_counts(tpec.begin(), arrayDims1[0], arrayDims1[1], arrayDims1[2], false);
    
    IntegerVector arrayDims2 = taec.attr("dim");
    arma::cube topic_absent_edge_counts(taec.begin(), arrayDims2[0], arrayDims2[1], arrayDims2[2], false);
    
    IntegerVector arrayDims3 = clp.attr("dim");
    arma::cube current_latent_positions(clp.begin(), arrayDims3[0], arrayDims3[1], arrayDims3[2], false);
    
    IntegerVector arrayDims4 = plp.attr("dim");
    arma::cube proposed_latent_positions(plp.begin(), arrayDims4[0], arrayDims4[1], arrayDims4[2], false);
    
    double sum_log_probability_of_current_positions = 0;
    double sum_log_probability_of_proposed_positions = 0;
    NumericVector current_author_position(number_of_latent_dimensions);
    NumericVector proposed_author_position(number_of_latent_dimensions);
    NumericVector recipient_position(number_of_latent_dimensions);
    
    
    for(int t = 0; t < number_of_topics; ++t){
		//get current topic intercept
		double current_topic_intercept = current_intercepts[t];
        double proposed_topic_intercept = proposed_intercepts[t];

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
                    //calculate linear predictor
                    double eta = current_topic_intercept - pow(distance,.5);
                    
                    //calculate likelihoods for both
                    double log_prob_edge = 0;
                    double log_prob_no_edge = 0;
                    if (eta != 0){
                        if(eta > 0){
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
                    //calculate linear predictor
                    eta = proposed_topic_intercept - pow(distance,.5);
                    //calculate likelihoods for both
                    log_prob_edge = 0;
                    log_prob_no_edge = 0;
                    if (eta != 0){
                        if(eta > 0){
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
    log_ratio = sum_log_probability_of_proposed_positions - sum_log_probability_of_current_positions;
    
    //take the log of a uniform draw on 0 to  1
    double lud = log_uniform_draw();
    
    if(log_ratio < lud){
        //if the log ratio is smaller then reject the new positions
        
        
    }
    else{
       //accept the new positions 
        
        
    }
    
            
    List to_return(4);
    to_return[0] = log_ratio;
    to_return[1] = sum_log_probability_of_current_positions; 
    to_return[2] = sum_log_probability_of_proposed_positions; 
    return to_return;
}

//Metropolis_Step_CPP(10,10,array(1:1000,c(10,10,10)),array(1:1000,c(10,10,10)),array(runif(200),c(2,10,10)),c(1:10),2,array(runif(200),c(2,10,10)),c(1:10))
            

