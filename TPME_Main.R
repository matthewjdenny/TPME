# Main class to run TPME Model
rm(list=ls())

# 1. Load packages
library(stringr)
library(statnet)
library(Rcpp)

# 2. Load user defined functions
source("TPME_Run_Analysis.R")

# 3. Load data: vocab file, document word matrix and document edge matrix

#======== McDowell - 2011 =========#
author_attributes =  read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/mcdowell/Email_Name_Dept_Gender_2011.csv", header= T, stringsAsFactors = F)
#making a judgement call that Robbin silvers in a woman
author_attributes[14,4] <- F
document_edge_matrix = read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/mcdowell/edge-matrix.csv", header= F, stringsAsFactors = F)
document_edge_matrix = document_edge_matrix[,-1]
document_edge_matrix[,1] <- document_edge_matrix[,1] + 1 #make sure that authors are indexed starting at 1
document_word_matrix = read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/mcdowell/word-matrix.csv", header= F, stringsAsFactors = F)
document_word_matrix = document_word_matrix[,-1]
vocabulary = read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/mcdowell/vocab.txt", header= F, stringsAsFactors = F)
#===================================#

# 4. Run analysis
Result <- Run_Analysis(Number_Of_Iterations = 50)
    

# 5. Output and Analyze results

