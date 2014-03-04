#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double TEST_CPP( 
int author, 
int t, 
NumericMatrix latent_positions1, 
NumericMatrix latent_positions2,
NumericVector edge_values,
NumericVector edge_topic_assignments){
    NumericVector author_position(2);
    author_position[0] = latent_positions1(t,(author-1));
    author_position[1] = latent_positions2(t,(author-1));
    double topic_intercept = 10;
    double additional_edge_probability = 0;
    
    for(int a = 0; a < 10; ++a){
                int actor = a + 1;
                
                if(author != actor){
                    int actual_edge = edge_values[a];
                    int edge_assignment = edge_topic_assignments[a];
                    if(edge_assignment == 3){
                        
                        //change one element in the array/cube
                        NumericVector recipient_position(2);
                        recipient_position[0] = latent_positions1(t,a);
                        recipient_position[1] = latent_positions2(t,a);
                        
                        double distance = 0;
            
                        for(int k = 0; k < 2; ++k){
                            distance += pow((author_position[k] - recipient_position[k]),2);
                        }

                        double eta = topic_intercept - pow(distance,.5);
                        double log_prob = 0;
                        if(eta > 0){
                            //we only have to deal with actual edges
                            log_prob = eta -log(1 + exp(eta));
                        }
                        else{
                            log_prob = -log(1 + exp(-eta));
                        }
                        additional_edge_probability += log_prob;
                    }   
                } 
            }
    
    return additional_edge_probability;      
}


//mat = matrix(1:100,10,10)
// TEST_CPP(1,1,mat,mat,rep(0,10),rep(1,10))


