#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector SAMPLE_TOKEN_TOPIC_ASSIGNMENTS_CPP(int number_of_tokens, int number_of_topics,int number_of_actors, int author, int document, double beta, NumericVector alpha_m){
    
    Function get_token_topic_assignment("get_token_topic_assignment");
    Function get_observed_edge_value("get_observed_edge_value");
    Function log_multinomial_draw("log_multinomial_draw");
    Function Log_Probability_Of_Edge("Log_Probability_Of_Edge");
    Function get_edge_topic_assignment("get_edge_topic_assignment");
    Function get_sum_token_topic_assignments("get_sum_token_topic_assignments");
    Function get_number_of_unique_words("get_number_of_unique_words");
    Function get_word_type_topic_assignemnt_count("get_word_type_topic_assignemnt_count");
    Function get_number_of_tokens_assigned_to_topic("get_number_of_tokens_assigned_to_topic");
     
    int number_of_word_types = as<int>(get_number_of_unique_words());
    
    NumericVector token_topic_assignments(number_of_tokens);
    for(int w = 0; w < number_of_tokens; ++w){
        int token = w + 1;
        
        NumericVector token_topic_distribution(number_of_topics);
        for(int t = 0; t < number_of_topics; ++t){
            int topic = t + 1;
            
            //this calculates the addition to the probability that the token was sampled from the topic by
            //adding together edge likelihoods associated with that topic in the document
            double additional_edge_probability = 0;
            for(int a = 0; a < number_of_actors; ++a){
                int actor = a + 1;
                
                if(author != actor){
                    int actual_edge = as<int>(get_observed_edge_value(document,actor));
                    int edge_assignment = as<int>(get_edge_topic_assignment(document, actor));
                    if(edge_assignment == topic){
                        additional_edge_probability += as<double>(Log_Probability_Of_Edge(topic,author,actor,actual_edge,0,0));
                    }   
                } 
            }
            
            //now we calculate the first and second terms in the likelihood of of the token being from the current topic
            int ntd = as<int>(get_sum_token_topic_assignments(document,token,topic));
            int wttac = as<int>(get_word_type_topic_assignemnt_count(document,token,topic));
            int ntt = as<int>(get_number_of_tokens_assigned_to_topic(document,token,topic));
            double first_term = ntd + (alpha_m[t]/number_of_topics);
            double second_term = (wttac + (beta/number_of_word_types))/(ntt + beta);
            
            
            token_topic_distribution[t] = log(first_term) + log(second_term) + additional_edge_probability;
        }
        token_topic_assignments[w] = as<double>(log_multinomial_draw(token_topic_distribution));        
    }
            
    return token_topic_assignments;
}


//d= 1
//SAMPLE_TOKEN_TOPIC_ASSIGNMENTS_CPP(Token_Topic_Assignments[[d]][[1]],Number_Of_Topics,Number_Of_Authors,Document_Authors[d],d,Beta,Alpha_Base_Measure_Vector)



