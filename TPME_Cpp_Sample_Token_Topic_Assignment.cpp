#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericVector Gibbs_Sampler_Cpp(double alpha,int index,int max,int token_number,NumericVector word_topic_probabilities,NumericVector topic_assignments, int tokens_in_doc){

    int num_topics = word_topic_probabilities.size() ;
    int maximum = max -1;
    int ind = index - 1;
    //generate conditional posterior distribution to sample topic assignment from
    NumericVector conditional_posterior(num_topics);

    for(int i = 0; i < num_topics; ++i) {
        int N = 0;
    
    if(maximum > 0){
        for(int k = 0; k < maximum; ++k) {
            int tmp = i+1;
            if(topic_assignments[k] == tmp){
                if(ind != k){
                    N += 1;
                }
            }
        }
    }
    double denom = tokens_in_doc  - 1 + (num_topics*alpha);
    double numerator = N + alpha;
    double multiply = numerator/denom;
    conditional_posterior[i] = word_topic_probabilities[i]*multiply;

    }
    
    return conditional_posterior;
  }
