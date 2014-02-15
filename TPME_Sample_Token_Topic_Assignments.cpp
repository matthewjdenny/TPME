#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector SAMPLE_TOKEN_TOPIC_ASSIGNMENTS_CPP(int number_of_tokens, int number_of_topics,int number_of_actors, int author, int document){
    
    Function get_token_topic_assignment("get_token_topic_assignment");
    Function get_observed_edge_value("get_observed_edge_value");
    Function log_multinomial_draw("log_multinomial_draw");
    Function Log_Probability_Of_Edge("Log_Probability_Of_Edge");
    Function get_edge_topic_assignment("get_edge_topic_assignment");
    
    int document_author = author - 1;
    
    NumericVector token_topic_assignments(number_of_tokens);
    for(int w = 0; w < (number_of_tokens -1); ++w){
        int token = w + 1;
        
        NumericVector token_topic_distribution(number_of_topics);
        for(int t = 0; t < (number_of_topics -1); ++t){
            int topic = t + 1;
            
            //this calculates the addition to the probability that the token was sampled from the topic by
            //adding together edge likelihoods associated with that topic in the document
            double additional_edge_probability = 0;
            for(int a = 0; a < (number_of_actors -1); ++a){
                int actor = a + 1;
                
                if(author != actor){
                    int actual_edge = as<int>(get_observed_edge_value(document,actor));
                    int edge_assignment = as<int>(get_edge_topic_assignment(document, actor));
                    if(edge_assignment == topic){
                        additional_edge_probability += as<double>(Log_Probability_Of_Edge(topic_assignment,author,actor,actual_edge,0,0));
                    }   
                } 
            }
            
            token_topic_distribution[t] = additional_edge_probability + 
            //need to come up with a new caching mechanism to keep track of token topic assignemnts
            //and number of tokens assigned to the topic across all documents. Wrap the call in an 
            //R function.
        }
        
                
    }
            
    return token_topic_assignments;
}






int topic_assignment = as<int>(get_token_topic_assignment(document,token));
int actual_edge = as<int>(get_observed_edge_value(document,recipient));
edge_log_probabilities[w] = as<double>(Log_Probability_Of_Edge(topic_assignment,author,recipient,actual_edge,0,0));  
            
