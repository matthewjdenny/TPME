#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericVector SAMPLE_EDGE_TOPIC_ASSIGNMENTS_CPP(int number_of_actors, int author, int number_of_tokens ,int document){
    
    Function get_token_topic_assignment("get_token_topic_assignment");
    Function get_observed_edge_value("get_observed_edge_value");
    Function log_multinomial_draw("log_multinomial_draw");
    Function Log_Probability_Of_Edge("Log_Probability_Of_Edge");
    int document_author = author - 1;
    
    NumericVector edge_token_assignments(number_of_actors);
     
    
    for(int a = 0; a < (number_of_actors -1); ++a){
        if(document_author != a){
            int recipient = a + 1;
            NumericVector edge_log_probabilities(number_of_tokens);
            for(int w = 0; w < (number_of_tokens -1); ++w){
                int token = w + 1;
                int topic_assignment = as<int>(get_token_topic_assignment(document,token));
                int actual_edge = as<int>(get_observed_edge_value(document,recipient));
                edge_log_probabilities[w] = as<double>(Log_Probability_Of_Edge(topic_assignment,author,recipient,actual_edge,0,0));  
            }
            edge_token_assignments[a] = as<double>(log_multinomial_draw(edge_log_probabilities));
        }
        
    }
    return edge_token_assignments;
}


            
