# Main class to run TPME Model
rm(list=ls())

# 1. Load packages
library(stringr)
library(statnet)
library(Rcpp)

# 2. Load user defined functions
source("TPME_Run_Analysis.R")

# 3. Load data: vocab file, document word matrix and document edge matrix


# 4. Run analysis
Result <- Run_Analysis( <- function(Number_Of_Iterations = 1000, Base_Alpha =1, Base_Beta = 0.01, Number_Of_Topics = 50, Author_Attributes, Document_Edge_Matrix ,Document_Word_Matrix, Vocabulary, Latent_Dimensions)
    

# 5. Output and Analyze results

