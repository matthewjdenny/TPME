# Main class to run TPME Model

#use this to open a shell in CentOS 6 that can be used to run the scripts
#scl enable devtoolset-1.1 bash
#source("~/Dropbox/PINLab/Projects/R_Code/TPMNE/TPME_Main.R")

# 1. Preliminaries
rm(list=ls())
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


# 3. Load data: vocab file, document word matrix, document edge matrix and actor covariates

# choose a dataset to work with:
#load("./Data/McDowell_2011_Data.Rdata")
#load("./Data/New_Hannover_2011_Data.Rdata")
#load("./Data/Transylvania_2011_Data.Rdata")
load("./Data/Columbus_2011_Data.Rdata")

# 4. Run analysis for 50,000 itterations by setting equal to 50


Results <- Run_Cluster_Integrated_Analysis(Number_Of_Iterations = 2000, Base_Alpha =1, Base_Beta = 0.01, Number_Of_Topics = 50, Topic_Step_Itterations = 1, Sample_Step_Itterations = 1000, output_file = "Columbus_7-15-14",Proposal_Variance = 0.5, seed = 1234, system_OS = "Linux", Number_of_Clusters = 10,Itterations_Before_Cluster_Assingment_Updates = 2, Adaptive_Metropolis_Target_Accept_Rate = 0.3, slice_sample_alpha_step_size = 1,TPME_Mode = F)



#Result <- Run_Cluster_Integrated_Analysis(Number_Of_Iterations = 100, Base_Alpha =1, Base_Beta = 0.01, Number_Of_Topics = 50, Topic_Step_Itterations = 1, Sample_Step_Itterations = 100, output_file = "Test2",Proposal_Variance = 0.2, system_OS = "Mac", Number_of_Clusters = 8,Itterations_Before_Cluster_Assingment_Updates = 0,TPME_Mode = F )



# 5. Model Diagnostic plots 

#Generate_Model_Diagnsotics(input_file = "New_Hannover_2011_5-20-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "New_Hannover_2011_5-21-14", Thin_Itterations = 1,skip_first = 600,Cluster_Integrated = T)






#5 run additional sample steps:

#Run_Sample_Step(input_file = "Model_Output_New_Hannover_2011_4-30-14",data_source = "New_Hannover_2011_Data", output_file = "New_Hannover_2011_5-17-14", data_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/",sample_step_burnin = 200000,itterations = 10200000,sample_every = 2000,post_burin_variance_multiplier = 0.5, system_OS = "Linux", max_cpus = 20)


# 5. Output and Analyze results
#Calculate_Network_Efficiency_Statistics(input_file = "Transylvania_Sample_10M_2011_3-13-14", output_file = "Transylvania_10M_2011_3-13-14", data_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", method = "Latora-Machiori",Thin_Itterations = 1,skip_first = 6200)

#Calculate_Network_Efficiency_Statistics(input_file = "McDowell_Sample_10M_2011_3-13-14", output_file = "McDowell_10M_2011_3-13-14", data_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", method = "Latora-Machiori",Thin_Itterations = 1,skip_first = 6200)

#Calculate_Network_Efficiency_Statistics(input_file = "Sample_Step_New_Hannover_2011_3-13-14", output_file = "New_Hannover_1M_2011_3-13-14", data_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", method = "Latora-Machiori",Thin_Itterations = 1,skip_first = 6200)







