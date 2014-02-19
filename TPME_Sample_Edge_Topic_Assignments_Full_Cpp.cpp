#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericVector SAMPLE_EDGE_TOPIC_ASSIGNMENTS_CPP(
    int number_of_actors, 
    int author, 
    int number_of_tokens,
    int document,
    NumericVector token_topic_assignments,
    NumericVector observed_edges,
    NumericMatrix latent_positions1, 
    NumericMatrix latent_positions2, 
    NumericVector intercepts, 
    int number_of_latent_dimensions
    ){
    
    //Function get_token_topic_assignment("get_token_topic_assignment");
    //Function get_observed_edge_value("get_observed_edge_value");
    Function log_multinomial_draw("log_multinomial_draw");
    //Function Log_Probability_Of_Edge("Log_Probability_Of_Edge");
    int document_author = author - 1;
    
    NumericVector edge_topic_assignments(number_of_actors);
     
    
    for(int a = 0; a < number_of_actors; ++a){
        if(document_author != a){
            int recipient = a + 1;
            
            NumericVector edge_log_probabilities(number_of_tokens);
            for(int w = 0; w < number_of_tokens; ++w){
                int token = w + 1;
                int topic_assignment = token_topic_assignments[w] -1;
                //int topic_assignment = as<int>(get_token_topic_assignment(document,token));
                int actual_edge = observed_edges[a];
                //int actual_edge = as<int>(get_observed_edge_value(document,recipient));
                //edge_log_probabilities[w] = as<double>(Log_Probability_Of_Edge(topic_assignment,author,recipient,actual_edge,0,0));
                //change one element in the array/cube
                NumericVector author_position(number_of_latent_dimensions);
                author_position[0] = latent_positions1(topic_assignment,document_author);
                author_position[1] = latent_positions2(topic_assignment,document_author);
                double topic_intercept = intercepts[topic_assignment];
                NumericVector recipient_position(number_of_latent_dimensions);
                recipient_position[0] = latent_positions1(topic_assignment,a);
                recipient_position[1] = latent_positions2(topic_assignment,a);
                        
                double distance = 0;
            
                for(int k = 0; k < number_of_latent_dimensions; ++k){
                    distance += pow((author_position[k] - recipient_position[k]),2);
                }

                double eta = topic_intercept - pow(distance,.5);
                
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
            //edge_topic_assignments[a] = as<int>(get_token_topic_assignment(document,sampled_token));
            edge_topic_assignments[a] = token_topic_assignments[sampled_token-1];
        }
        
    }
    return edge_topic_assignments;
}


            
