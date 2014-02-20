#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Tests the array functionality of RcppArmadillo

// [[Rcpp::export]]
List test(NumericVector myArray){
    IntegerVector arrayDims = myArray.attr("dim");
 
    List out(2);
  arma::cube cubeArray(myArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
 
  //change one element in the array/cube
  cubeArray(0,0,0) = 518;  
  
  out[0] = cubeArray;
  out[1] = 4;
 
  return(out); 

}
  