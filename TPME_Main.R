# Main class to run TPME Model

#source("~/Dropbox/PINLab/Projects/R_Code/TPMNE/TPME_Main.R")

# 1. Preliminaries
rm(list=ls())
setwd("~/Dropbox/PINLab/Projects/R_Code/TPMNE")

#requires RcppArmadillo and all dependenices be installed
#install.packages("RcppArmadillo")

# 2. Load user defined functions
source("./Scripts/TPME_Run_Analysis.R")
source("./Scripts/TPME_Model_Diagnostics.R")

# 3. Load data: vocab file, document word matrix, document edge matrix and actor covariates

# choose a dataset to work with:
#load("./Data/McDowell_2011_Data.Rdata")
#load("./Data/New_Hannover_2011_Data.Rdata")
#load("./Data/Transylvania_2011_Data.Rdata")
load("./Data/Columbus_2011_Data.Rdata")

# 4. Run analysis for 50,000 itterations by setting equal to 50
Model_Accept_Rate <- Run_Analysis(Number_Of_Iterations = 50,Run_Sample_Step = T,output_file = "Columbus_2011_3-7-14",Base_Alpha =.1, Base_Beta = 0.01, Number_Of_Topics = 50,Proposal_Variance_Vector = c(.5,.1),post_burin_variance = 0.01)

#Model_Accept_Rate <- Run_Analysis(Number_Of_Iterations = 50,Run_Sample_Step = T,output_file = "McDowell_2011_3-7-14",Base_Alpha =1, Base_Beta = 0.01, Number_Of_Topics = 50,Proposal_Variance_Vector = c(.5,.1),post_burin_variance = 0.001)

#Model_Accept_Rate <- Run_Analysis(Number_Of_Iterations = 50,Run_Sample_Step = T,output_file = "Transylvania_2011_3-7-14",Base_Alpha =.1, Base_Beta = 0.01, Number_Of_Topics = 50,Proposal_Variance_Vector = c(.5,.1,.05),post_burin_variance = 0.005)

#Model_Accept_Rate <- Run_Analysis(Number_Of_Iterations = 50,Run_Sample_Step = T,output_file = "New_Hannover_2011_3-7-14",Base_Alpha =.1, Base_Beta = 0.01, Number_Of_Topics = 100,Proposal_Variance_Vector = c(.5,.1,.05),post_burin_variance = 0.005)

# 5. Model Diagnostic plots 
#Generate_Model_Diagnsotics(input_folder_path = "./Output/",input_file = "Sample_Step_Columbus_2011_3-5-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "Columbus_County")

#Generate_Model_Diagnsotics(input_folder_path = "./Output/",input_file = "Sample_Step_McDowell_2011_3-5-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "McDowell_County")

#Generate_Model_Diagnsotics(input_file = "Current_Itteration_McDowell_2011_3-6-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "McDowell_County_3-6-14", Thin_Itterations = 1)

#Generate_Model_Diagnsotics(input_file = "Current_Itteration_New_Hannover_2011_3-6-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "New_Hannover_County_3-6-14", Thin_Itterations = 1)

#Generate_Model_Diagnsotics(input_file = "Sample_Step_McDowell_2011_3-6-14",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", vocab = vocabulary,county_name = "McDowell_County_Sample_Step_3-6-14", Thin_Itterations = 1)



# 5. Output and Analyze results

