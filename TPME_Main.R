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
source("./Scripts/TPME_Model_Diagnostics.R")
source("./Scripts/TPME_Take_Sample.R")
source("./Scripts/TPME_Calculate_Network_Efficiency_Statistics.R")

Report_Probs <- function(current){
    print(paste("Current Save Iteration:",current) )
}


# 3. Load data: vocab file, document word matrix, document edge matrix and actor covariates

# choose a dataset to work with:
#load("./Data/McDowell_2011_Data.Rdata")
#load("./Data/New_Hannover_2011_Data.Rdata")
#load("./Data/Transylvania_2011_Data.Rdata")
#load("./Data/Columbus_2011_Data.Rdata")

# 4. Run analysis for 50,000 itterations by setting equal to 50
#Model_Accept_Rate <- Run_Analysis(Number_Of_Iterations = 50,Run_Sample_Step = T,output_file = "Columbus_2011_3-13-14",Base_Alpha =.1, Base_Beta = 0.01, Number_Of_Topics = 50,Proposal_Variance_Vector = c(.5,.1),post_burin_variance = 0.05, system_OS = "Linux", sampler = "Slice", slice_sample_step_size = 1)

#Model_Accept_Rate <- Run_Analysis(Number_Of_Iterations = 50,Run_Sample_Step = T,output_file = "McDowell_2011_3-13-14",Base_Alpha =1, Base_Beta = 0.01, Number_Of_Topics = 50,Proposal_Variance_Vector = c(.5,.1),post_burin_variance = 0.05,Metropolis_Step_Itterations = 1000, system_OS = "Linux", sampler = "Slice", slice_sample_step_size = 1)

#Model_Accept_Rate <- Run_Analysis(Number_Of_Iterations = 50,Run_Sample_Step = T,output_file = "Transylvania_2011_3-13-14",Base_Alpha =.1, Base_Beta = 0.01, Number_Of_Topics = 50,Proposal_Variance_Vector = c(.5,.1),post_burin_variance = 0.05, system_OS = "Linux", sampler = "Slice", slice_sample_step_size = 1)

#Model_Accept_Rate <- Run_Analysis(Number_Of_Iterations = 50,Run_Sample_Step = T,output_file = "New_Hannover_2011_3-14-14",Base_Alpha =.1, Base_Beta = 0.01, Number_Of_Topics = 100,Proposal_Variance_Vector = c(.5,.1,.01),post_burin_variance = 0.005,Sample_Step_Itterations = 10200000,Sample_Every = 1000, system_OS = "Linux", sampler = "Slice", slice_sample_step_size = 1)

# ======== for benchmarking ========= #
#print(system.time(Model_Accept_Rate <- Run_Analysis(Number_Of_Iterations = 1,Run_Sample_Step = F,output_file = "McDowell_10K_2011_3-10-14",Base_Alpha =1, Base_Beta = 0.01, Number_Of_Topics = 50,Proposal_Variance_Vector = c(.5,.1),post_burin_variance = 0.05,Metropolis_Step_Itterations = 1000,Topic_Step_Itterations = 1000, system_OS = "Mac", sampler = "Slice", slice_sample_step_size = 1)))



# 5. Model Diagnostic plots 

#Generate_Model_Diagnsotics(input_file = "Sample_Step_McDowell_2011_3-13-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "McDowell_County_Sample_10_3-13-14", Thin_Itterations = 1,skip_first = 2000)

#Generate_Model_Diagnsotics(input_file = "Sample_Step_Transylvania_2011_3-13-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "Transylvania_County_Sample_10_3-13-14", Thin_Itterations = 1,skip_first = 10000)

#Generate_Model_Diagnsotics(input_file = "Sample_Step_Columbus_2011_3-13-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "Columbus_County_Sample_10_3-13-14", Thin_Itterations = 1,skip_first = 2000)

#Generate_Model_Diagnsotics(input_file = "Columbus_Sample_10M_2011_3-13-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "Columbus_County_Sample_10M_3-13-14", Thin_Itterations = 1,skip_first = 5200)

#Generate_Model_Diagnsotics(input_file = "Columbus_Sample_10M_2011_3-13-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "Columbus_10M_3-13-14", Thin_Itterations = 1,skip_first = 7200)

#Generate_Model_Diagnsotics(input_file = "Sample_Step_New_Hannover_2011_3-14-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "New_Hannover_10M_3-14-14", Thin_Itterations = 1,skip_first = 8200)

#Generate_Model_Diagnsotics(input_file = "Transylvania_Slice_2M_2011_3-13-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "Transylvania_2M_slice_3-13-14", Thin_Itterations = 1,skip_first = 8200)


#Generate_Model_Diagnsotics(input_file = "Transylvania_Sample_50M_2011_3-13-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "Transylvania_50M_3-13-14", Thin_Itterations = 1,skip_first = 40200)


#Generate_Model_Diagnsotics(input_file = "McDowell_Sample_50M_2011_3-13-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "McDowell_50M_3-13-14", Thin_Itterations = 1,skip_first = 7000)

#Generate_Model_Diagnsotics(input_file = "McDowell_Slice_2M_2011_3-13-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "McDowell_2M_slice_3-13-14", Thin_Itterations = 1,skip_first = 6200)


#5 run additional sample steps:

#Run_Sample_Step(input_file = "Current_Itteration_McDowell_2011_3-13-14",data_source = "McDowell_2011_Data", output_file = "McDowell_Slice_2M_2011_3-13-14", itterations = 2050000, proposal_variance = 0.25,sample_every = 200,sample_step_burnin = 50000,post_burin_variance = 0.2, system_OS = "Linux", sampler = "Block", slice_sample_step_size = 1,post_burnin_step_size = 0.5)

#Run_Sample_Step(input_file = "Current_Itteration_Columbus_2011_3-13-14",data_source = "Columbus_2011_Data", output_file = "Columbus_Slice_2M_2011_3-13-14", itterations = 2050000, proposal_variance = 0.1,sample_every = 200,sample_step_burnin = 50000,post_burin_variance = 0.02, system_OS = "Linux", sampler = "Block", slice_sample_step_size = 1,post_burnin_step_size = 0.5)

#Run_Sample_Step(input_file = "Current_Itteration_Transylvania_2011_3-13-14",data_source = "Transylvania_2011_Data", output_file = "Transylvania_Slice_2M_2011_3-13-14", itterations = 2050000, proposal_variance = 0.1,sample_every = 200,sample_step_burnin = 50000,post_burin_variance = 0.05, system_OS = "Linux", sampler = "Block", slice_sample_step_size = 1,post_burnin_step_size = 0.5)

#Run_Sample_Step(input_file = "Current_Itteration_New_Hannover_2011_3-14-14",data_source = "New_Hannover_2011_Data", output_file = "New_Hannover_Slice_2M_2011_3-14-14", itterations = 2050000, proposal_variance = 0.02,sample_every = 200,sample_step_burnin = 50000,post_burin_variance = 0.01, system_OS = "Linux", sampler = "Block", slice_sample_step_size = 1,post_burnin_step_size = 0.5)





# 5. Output and Analyze results
#Calculate_Network_Efficiency_Statistics(input_file = "Transylvania_Sample_10M_2011_3-13-14", output_file = "Transylvania_10M_2011_3-13-14", data_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", method = "Latora-Machiori",Thin_Itterations = 1,skip_first = 6200)

#Calculate_Network_Efficiency_Statistics(input_file = "McDowell_Sample_10M_2011_3-13-14", output_file = "McDowell_10M_2011_3-13-14", data_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", method = "Latora-Machiori",Thin_Itterations = 1,skip_first = 6200)

#Calculate_Network_Efficiency_Statistics(input_file = "Sample_Step_New_Hannover_2011_3-13-14", output_file = "New_Hannover_1M_2011_3-13-14", data_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", method = "Latora-Machiori",Thin_Itterations = 1,skip_first = 6200)







