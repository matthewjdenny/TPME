Run_Sample_Step <- function(input_file = "Current_Itteration_McDowell_2011_3-7-14",data_source = "McDowell_2011_Data", output_file = "Sample_McDowell_2011_3-7-14", itterations = 1200000, proposal_variance = 0.01,sample_every = 100, data_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/",sample_step_burnin = 200000,post_burin_variance = 0.01, system_OS = "Linux", sampler = c("Slice","Metropolis", "Block"), slice_sample_step_size = 1,post_burnin_step_size = 0.001, within_block_iterations = 1){
    #load the data
    load(paste(data_directory,input_file,".Rdata", sep = ""))
    
    #if we are running linux then we need to add the appropriate c flags to use c++2011
    if(system_OS == "Linux"){
        PKG_CPPFLAGS = "-std=c++11"
        Sys.setenv(PKG_CPPFLAGS = PKG_CPPFLAGS)
    }
    
    
    
    Rcpp::sourceCpp("./Scripts/TPME_Take_Metropolis_Sample.cpp")
    Rcpp::sourceCpp("./Scripts/TPME_Take_Slice_Sample.cpp")
    Rcpp::sourceCpp('Scripts/TPME_Take_Block_Slice_Sample.cpp')
    
    
    
    print("Loading Data")
    #extract current metropolis results
    Sample_Results <- Return_List[[1]]
    Topic_Model_Results <- Return_List[[2]]
    #free up memory
    rm(Return_List)
    
    load(paste("./Data/",data_source,".Rdata", sep = ""))
    Author_Attributes= author_attributes
    Number_of_Betas <- 4
    
    
    log_multinomial_draw <- function(probability_vector){
        return(which(rmultinom(1,1,exp(probability_vector)) == 1))
    }
    
    log_uniform_draw <- function(){
        return(log(runif(1, min=0, max=1)))
    }
    
    #get model information and extract data
    Latent_Dimensions <- length(Sample_Results[[1]][,1,1])
    Number_Of_Topics <- length(Sample_Results[[1]][1,,1])
    Number_Of_Authors <- length(Sample_Results[[1]][1,1,])
    Token_Topic_Assignments <- Topic_Model_Results[[1]]
    Topic_Present_Edge_Counts <- Topic_Model_Results[[2]]
    Topic_Absent_Edge_Counts <- Topic_Model_Results[[3]]
    Word_Type_Topic_Counts <- Topic_Model_Results[[4]]
    Edge_Topic_Assignments <- Topic_Model_Results[[5]]
    Latent_Space_Positions <- Sample_Results[[1000]]
    Latent_Space_Intercepts <- Sample_Results[[2*1000]]
    Betas <- Sample_Results[[3*1000]]
    
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
    
    print("Running Model")
    if(sampler == "Slice"){
        Results <- Slice_Sample_CPP(
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
            itterations,
            slice_sample_step_size,
            array(0,c(Latent_Dimensions,Number_Of_Topics,Number_Of_Authors)),
            array(0,c(Latent_Dimensions,Number_Of_Topics,Number_Of_Authors)),
            array(0,c(Latent_Dimensions,Number_Of_Topics,Number_Of_Authors)),
            sample_every,
            sample_step_burnin,
            post_burnin_step_size
        )
    }
    
    if(sampler == "Metropolis"){
        #run final metropolis step with more itterations 
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
            itterations,
            Proposal_Variance,
            array(0,c(Latent_Dimensions,Number_Of_Topics,Number_Of_Authors)),
            ample_every,
            sample_step_burnin,
            post_burin_variance
        )
    }
    
    if(sampler == "Block"){
        Results <- Block_Slice_Sample_CPP(
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
            itterations,
            slice_sample_step_size,
            array(0,c(Latent_Dimensions,Number_Of_Topics,Number_Of_Authors)),
            array(0,c(Latent_Dimensions,Number_Of_Topics,Number_Of_Authors)),
            array(0,c(Latent_Dimensions,Number_Of_Topics,Number_Of_Authors)),
            sample_every,
            within_block_iterations,
            sample_step_burnin,
            post_burnin_step_size
        )
        
        
    }
    
    #unlist(Results[60001:72000])
    print("Saving Results")
    Return_List <- list(Results,Topic_Model_Results)
    
    save(Return_List, file=paste("~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/",output_file,".Rdata",sep = ""))

}

    
