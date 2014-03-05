# Main class to run TPME Model



# 1. Preliminaries
rm(list=ls())
setwd("~/Dropbox/PINLab/Projects/R_Code/TPMNE")

#requires RcppArmadillo and all dependenices be installed
install.packages("RcppArmadillo")

# 2. Load user defined functions
source("./Scripts/TPME_Run_Analysis.R")

# 3. Load data: vocab file, document word matrix, document edge matrix and actor covariates

# choose a dataset to work with:
# load("./Data/McDowell_2011_Data.Rdata")
# load("./Data/New_Hannover_2011_Data.Rdata")
# load("./Data/Transylvania_2011_Data.Rdata")
# load("./Data/Columbus_2011_Data.Rdata")

# 4. Run analysis for 50,000 itterations by setting equal to 50
Model_Accept_Rate <- Run_Analysis(Number_Of_Iterations = 50,Run_Sample_Step = T,output_file = "Results_McDowell_3-4-14")
    

# 5. Output and Analyze results

