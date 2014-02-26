#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Tests the array functionality of RcppArmadillo

// [[Rcpp::export]]
List test(NumericVector myArray){
    IntegerVector arrayDims = myArray.attr("dim");
    double x = Rf_rnorm(0.0,1.0);
 
    List out(2);
  arma::cube cubeArray(myArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
 
  //change one element in the array/cube
  cubeArray(0,0,0) = 518;  
  
  out[0] = cubeArray;
  out[1] = x;
 
  return(out); 

}

// [[Rcpp::export]]
List test2(
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
    double proposal_variance,
    int number_of_documents,
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
  
    List to_return(1);
    to_return[0] = 2;
    return to_return;
    }  
  
    
    
    
    