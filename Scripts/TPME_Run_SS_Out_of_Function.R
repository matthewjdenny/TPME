
#use this to open a shell in CentOS 6 that can be used to run the scripts
#scl enable devtoolset-1.1 bash
#source("~/Dropbox/PINLab/Projects/R_Code/TPMNE/Scripts/TPME_Run_SS_Out_of_Function.R")

setwd("~/Dropbox/PINLab/Projects/R_Code/TPMNE")

#requires RcppArmadillo and all dependenices be installed
#install.packages("RcppArmadillo")

# 2. Load user defined functions
source("./Scripts/TPME_Run_Analysis.R")
source("./Scripts/TPME_Run_Cluster_Integrated_Analysis.R")
source("./Scripts/TPME_Model_Diagnostics.R")
#source("./Scripts/TPME_Take_Sample.R")
source("./Scripts/TPME_Cluster_Integrated_Take_Final_MH_Sample.R")
source("./Scripts/TPME_Calculate_Network_Efficiency_Statistics.R")

Report_Probs <- function(current){
    print(current)
    #print(paste("Current Save Iteration:",current) )
}

Report_2 <- function(current){
    print(sum(as.vector(apply(current,2,sum))))
    #print(paste("Current Save Iteration:",current) )
}





    input_file = "Model_Output_New_Hannover_2011_4-30-14"
    data_source = "New_Hannover_2011_Data"
    output_file = "New_Hannover_2011_5-20-14"
    itterations = 10200000
    sample_every = 2000
    data_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/"
    sample_step_burnin = 200000
    post_burin_variance_multiplier = 0.5
    system_OS = "Linux"
    max_cpus = 4

#load the data
require(Rcpp)
require(RcppArmadillo)
load(paste(data_directory,input_file,".Rdata", sep = ""))
#if we are running linux then we need to add the appropriate c flags to use c++2011
if(system_OS == "Linux"){
    PKG_CPPFLAGS = "-std=c++11 -g -og"
    Sys.setenv(PKG_CPPFLAGS = PKG_CPPFLAGS)
}


set.seed(1234)
#Rcpp::sourceCpp("./Scripts/TPME_Cluster_Integrated_Final_MH_Step.cpp")
Rcpp::sourceCpp("./Scripts/Test3.cpp")



print("Loading Data")
#extract current metropolis results
first_return <- 13
Topic_Model_Results <- Return_List[1:5]
Model_Parameters <- Return_List[6:first_return]
Cluster_Topic_Assignments <- Return_List[(first_return+1):(first_return+Model_Parameters[[2]])]
Last_Cluster_Topic_Assignments <- unlist(Cluster_Topic_Assignments[Model_Parameters[[2]]])
Metropolis_Results <- Return_List[(first_return+1+Model_Parameters[[2]]):length(Return_List)]
#free up memory
#rm(Return_List)

load(paste("./Data/",data_source,".Rdata", sep = ""))
Author_Attributes= author_attributes
Number_of_Betas <- 4


#remove the first skip_first itterations of each sublist and recombine
Itterations <- Model_Parameters[[4]]


#get model information and extract data
Latent_Spaces <- length(Metropolis_Results[[(Itterations + 1)]][,1,1])
Clusters <- length(Metropolis_Results[[1]])
Topics <- length(Cluster_Topic_Assignments[[1]])
Actors <- length(Metropolis_Results[[(Itterations + 1)]][1,1,])
Token_Topic_Assignments <- Topic_Model_Results[[1]]
Topic_Present_Edge_Counts <- Topic_Model_Results[[2]]
Topic_Absent_Edge_Counts <- Topic_Model_Results[[3]]
Word_Type_Topic_Counts <- Topic_Model_Results[[4]]
Edge_Topic_Assignments <- Topic_Model_Results[[5]]
#Proposal_Variances <- Model_Parameters[[6]][Model_Parameters[[2]],]
Cluster_Topic_Assigns <- Cluster_Topic_Assignments[[Model_Parameters[[2]]]]
LSPs <- Metropolis_Results[[2*Itterations]]
Intercepts <- Metropolis_Results[[Itterations]]
Betas <- Metropolis_Results[[3*Itterations]]

Proposal_Variances <- rep(0.2,Clusters)

#get the total number of tokens assigned to each topic
Topic_Token_Totals <- apply(Word_Type_Topic_Counts,2,sum)
#get the total number of present edges assigned to each topic
Email_Assignments <- apply(Topic_Present_Edge_Counts,3,sum)
clp <- array(0,c(Latent_Spaces,Clusters,Actors))



#generate indicator array
Beta_Indicator_Array <- array(0,c(Actors,Actors,Number_of_Betas))
for(j in 1:Actors){
    for(k in 1:Actors){
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
#     Single_Cluster_MH <- function(Cur_Cluster){
#         
#         Report_Probs <- function(current){
#             print(current)
#             #print(paste("Current Save Iteration:",current) )
#         }
#         
#         Report_2 <- function(current){
#             print(sum(as.vector(apply(current,2,sum))))
#             #print(paste("Current Save Iteration:",current) )
#         }
#         require(Rcpp)
#         require(RcppArmadillo)
#         set.seed(1234)
#         Rcpp::sourceCpp("/Users/matthewjdenny/Dropbox/PINLab/Projects/R_Code/TPMNE/Scripts/Test.cpp")
#         
#         Result <- Cluster_Integrated_Final_MH_Sampler(
#             itterations,
#             as.integer(Actors), 
#             as.integer(Topics),
#             as.integer(Latent_Spaces),
#             Proposal_Variances,
#             Cluster_Topic_Assigns,
#             Topic_Present_Edge_Counts,
#             Topic_Absent_Edge_Counts,
#             LSPs,
#             clp,
#             Intercepts,
#             Betas,
#             Number_of_Betas,
#             Beta_Indicator_Array,
#             as.matrix(Edge_Topic_Assignments),
#             sample_step_burnin,
#             as.integer(Clusters),
#             Cur_Cluster,
#             sample_every,
#             post_burin_variance_multiplier)
#         return(Result)
#     }

# if(max_cpus < Clusters){
#     numcpus <- max_cpus
# }else{
#     numcpus <- as.integer(Clusters)
# }
# 
#     require(doMC)
#     require(foreach)
#     registerDoMC(numcpus) #set number of cores
# 
#     Result_List <- foreach(i=1:Clusters) %dopar% {
#         temp <- Cluster_Integrated_Final_MH_Sampler(
#             itterations,
#             as.integer(Actors), 
#             as.integer(Topics),
#             as.integer(Latent_Spaces),
#             Proposal_Variances,
#             Cluster_Topic_Assigns,
#             Topic_Present_Edge_Counts,
#             Topic_Absent_Edge_Counts,
#             LSPs,
#             clp,
#             Intercepts,
#             Betas,
#             Number_of_Betas,
#             Beta_Indicator_Array,
#             as.matrix(Edge_Topic_Assignments),
#             sample_step_burnin,
#             as.integer(Clusters),
#             i,
#             sample_every,
#             post_burin_variance_multiplier)
#         temp <- list(temp)
#     }


#     system.time(
Result_List = list()
for(i in 1:8){
    temp <- Cluster_Integrated_Final_MH_Sampler(
        itterations,
        as.integer(Actors), 
        as.integer(Topics),
        as.integer(Latent_Spaces),
        Proposal_Variances,
        Cluster_Topic_Assigns,
        Topic_Present_Edge_Counts,
        Topic_Absent_Edge_Counts,
        LSPs,
        clp,
        Intercepts,
        Betas,
        Number_of_Betas,
        Beta_Indicator_Array,
        as.matrix(Edge_Topic_Assignments),
        sample_step_burnin,
        as.integer(Clusters),
        i,
        sample_every,
        post_burin_variance_multiplier)
    Result_List <- append(Result_List,list(temp))
    
}
#     )

#     library(snowfall) #load the snowfall package
#     # set up number of cpus
#     if(max_cpus < Clusters){
#         numcpus <- max_cpus
#     }else{
#         numcpus <- as.integer(Clusters)
#     }
#      #set the number of cpus to use
#     sfInit(parallel=TRUE, cpus=numcpus ) #initialize the snowfall cluster
#     if(sfParallel()){
#         cat( "Running in parallel mode on", sfCpus(), "nodes.\n" )
#     }else{
#         cat( "Running in sequential mode.\n" )
#     }
#     #export all packages currently loaded in your R session
#     for (i in 1:length(.packages())){
#         eval(call("sfLibrary", (.packages()[i]), character.only=TRUE))
#     }
#     #export a list of R data objects that your function will need
#     sfExport(list = ls())
#     
#     indexes <- 1:as.integer(Clusters)
#     #apply a function across the cluster
#     result <- sfClusterApplyLB(indexes,Single_Cluster_MH)
#     #stop the cluster (this needs to be done explicitly)
#     sfStop()


#store all data from the main model so it can be used again -- we jsut need to update the number of metropolis iterations
#Result_List <- Results_List
Ret_Lis <- Return_List[1:(first_return+Model_Parameters[[2]])]

#order of the composite list we must construct
# ints;
# lat_pos;
# bets;
# cluster_accepted;
# Cur_Proposed_MH_Likelihoods;
# Cur_Current_MH_Likelihoods;
#Metropolis_Results <- Metropolis_Results[seq(1, length(Metropolis_Results),Thin_Itterations)]

num_clusters <- Clusters
Bulk_Data_List <- list()
for(i in 1:num_clusters){
    Bulk_Data_List <- append(Bulk_Data_List,list(Result_List[[i]][1:5]))
    Result_List[[i]] <- Result_List[[i]][6:length(Result_List[[i]])]
}

Met_Res <- vector("list",3*length(Result_List[[1]][[1]]) )


for(i in 1:(length(Met_Res)/3)){
    vect <- rep(0,num_clusters)
    for(j in 1:num_clusters){
        vect[j] <- Result_List[[j]][[1]][i]
    }
    print(vect)
    print(i)
    Met_Res[i] <- list(vect)
}

for(i in ((length(Met_Res)/3)+1):(2*(length(Met_Res)/3))){
    arr <- array(0,c(Latent_Spaces,num_clusters,Actors))
    tmp <-  i -(length(Met_Res)/3)
    for(j in 1:num_clusters){
        arr[,j,] <- Result_List[[j]][[2]][,,tmp]
    }
    print(i)
    Met_Res[i] <- list(arr)
}

for(i in ((2*length(Met_Res)/3)+1):(3*(length(Met_Res)/3))){
    mat <- matrix(0,ncol = Number_of_Betas, nrow = num_clusters)
    tmp <-  i -2*(length(Met_Res)/3)
    for(j in 1:num_clusters){
        mat[j,] <- Result_List[[j]][[3]][,tmp]
    }
    print(i)
    Met_Res[i] <- list(mat)
}


#     for(j in 1:num_clusters){
#         for(i in 3:5){
#             Bulk_Data_List[[j]][[i]] <- Bulk_Data_List[[j]][[i]][seq(1, length(Bulk_Data_List[[j]][[i]]),100)]
#         } 
#     }
#now deal with accept rates, mh likelihoods and mh proposed likelihoods
Met_Res2 <- vector("list",3*length(Result_List[[1]][[1]]) )

for(i in 1:(length(Met_Res)/3)){
    vect <- rep(0,num_clusters)
    for(j in 1:num_clusters){
        vect[j] <- Bulk_Data_List[[j]][[5]][i]
    }
    print(i)
    print(vect)
    Met_Res2[i] <- list(vect)
}
for(i in ((length(Met_Res)/3)+1):(2*(length(Met_Res)/3))){
    vect <- rep(0,num_clusters)
    it <- i - (length(Met_Res)/3)
    for(j in 1:num_clusters){
        vect[j] <- Bulk_Data_List[[j]][[3]][it]
    }
    print(vect)
    print(i)
    Met_Res2[i] <- list(vect)
}
for(i in ((2*length(Met_Res)/3)+1):(3*(length(Met_Res)/3))){
    vect <- rep(0,num_clusters)
    it <- i - 2*(length(Met_Res)/3)
    for(j in 1:num_clusters){
        vect[j] <- Bulk_Data_List[[j]][[4]][it]
    }
    print(vect)
    print(i)
    Met_Res2[i] <- list(vect)
}

#get the right number of metropolis iterations
Ret_Lis[[9]] <- (length(Met_Res)/3)

REt_Lis <- append(Ret_Lis,Met_Res)
REt_Lis <- append(REt_Lis,Met_Res2)

Return_List <- REt_Lis
print("saving")
save(list = ls(), file=paste("~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/",output_file,".Rdata",sep = ""))