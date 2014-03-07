rm(list= ls())

load("~/Dropbox/PINLab/Projects/R_Code/TPMNE/Output/Current_Itteration_McDowell_2011_3-5-14.Rdata")


print("Extracting Reduced Data")
#extract current metropolis results
Metropolis_Results <- Return_List[[1]]
Topic_Model_Results <- Return_List[[2]]
#free up memory
rm(Return_List)

load("./Data/McDowell_2011_Data.Rdata")
Author_Attributes= author_attributes
Number_of_Betas <- 4


log_multinomial_draw <- function(probability_vector){
    return(which(rmultinom(1,1,exp(probability_vector)) == 1))
}

log_uniform_draw <- function(){
    return(log(runif(1, min=0, max=1)))
}

#get model information and extract data
Latent_Dimensions <- length(Metropolis_Results[[1]][,1,1])
Number_Of_Topics <- length(Metropolis_Results[[1]][1,,1])
Number_Of_Authors <- length(Metropolis_Results[[1]][1,1,])
Token_Topic_Assignments <- Topic_Model_Results[[1]]
Topic_Present_Edge_Counts <- Topic_Model_Results[[2]]
Topic_Absent_Edge_Counts <- Topic_Model_Results[[3]]
Word_Type_Topic_Counts <- Topic_Model_Results[[4]]
Edge_Topic_Assignments <- Topic_Model_Results[[5]]
Latent_Space_Positions <- Metropolis_Results[[1000]]
Latent_Space_Intercepts <- Metropolis_Results[[2*1000]]
Betas <- Metropolis_Results[[3*1000]]

require(RcppArmadillo)
set.seed(1234)
Rcpp::sourceCpp("./Scripts/TPME_Take_Metropolis_Sample.cpp")



#generate indicator array
Beta_Indicator_Array <- array(0,c(Number_Of_Authors,Number_Of_Authors,Number_of_Betas))
for(j in 1:Number_Of_Authors){
    for(k in 1:Number_Of_Authors){
        if(Author_Attributes$Gender[j] == "M" & Author_Attributes$Gender[k] == "M"){
            Beta_Indicator_Array[j,k,1] = 1   
        }
        if(Author_Attributes$Gender[j] == "M" & Author_Attributes$Gender[k] == "F"){
            Beta_Indicator_Array[j,k,2] = 1   
        }
        if(Author_Attributes$Gender[j] == "F" & Author_Attributes$Gender[k] == "M"){
            Beta_Indicator_Array[j,k,3] = 1   
        }
        if(Author_Attributes$Gender[j] == "F" & Author_Attributes$Gender[k] == "F"){
            Beta_Indicator_Array[j,k,4] = 1   
        }
    }
}


Results <- Metropolis_Sample_CPP(
    Number_Of_Authors, 
    Number_Of_Topics,
    Topic_Present_Edge_Counts,
    Topic_Absent_Edge_Counts,
    Latent_Space_Positions,
    Latent_Space_Intercepts,
    Latent_Dimensions,
    Betas,
    Number_of_Betas,
    Beta_Indicator_Array,
    10000,
    0.01,
    array(0,c(Latent_Dimensions,Number_Of_Topics,Number_Of_Authors)),
    100
)


