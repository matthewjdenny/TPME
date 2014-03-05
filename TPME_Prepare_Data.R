#prepare data for analysis with model
rm(list = ls())
setwd("~/Dropbox/PINLab/Projects/R_Code/TPMNE")

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
save(list = ls(),file = "./Data/McDowell_2011_Data.Rdata")
#===================================#

#======== New Hannover - 2011 =========#
author_attributes =  read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/nhc/Email_Name_Dept_Gender_2011.csv", header= T, stringsAsFactors = F)
document_edge_matrix = read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/nhc/edge-matrix.csv", header= F, stringsAsFactors = F)
document_edge_matrix = document_edge_matrix[,-1]
document_edge_matrix[,1] <- document_edge_matrix[,1] + 1 #make sure that authors are indexed starting at 1
document_word_matrix = read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/nhc/word-matrix.csv", header= F, stringsAsFactors = F)
document_word_matrix = document_word_matrix[,-1]
vocabulary = read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/nhc/vocab.txt", header= F, stringsAsFactors = F)

remove <- which(apply(document_word_matrix,1,sum) == 0)
document_edge_matrix <- document_edge_matrix[-remove,]
document_word_matrix <- document_word_matrix[-remove,]
save(list = ls(),file = "New_Hannover_2011_Data.Rdata")
#===================================#

#======== Transylvania - 2011 =========#
author_attributes =  read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/Transylvania/Email_Name_Dept_Gender_2011.csv", header= T, stringsAsFactors = F)
document_edge_matrix = read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/Transylvania/edge-matrix.csv", header= F, stringsAsFactors = F)
document_edge_matrix = document_edge_matrix[,-1]
document_edge_matrix[,1] <- document_edge_matrix[,1] + 1 #make sure that authors are indexed starting at 1
document_word_matrix = read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/Transylvania/word-matrix.csv", header= F, stringsAsFactors = F)
document_word_matrix = document_word_matrix[,-1]
vocabulary = read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/Transylvania/vocab.txt", header= F, stringsAsFactors = F)

remove <- which(apply(document_word_matrix,1,sum) == 0)
document_edge_matrix <- document_edge_matrix[-remove,]
document_word_matrix <- document_word_matrix[-remove,]
save(list = ls(),file = "Transylvania_2011_Data.Rdata")
#===================================#

#======== Columbus - 2011 =========#
author_attributes =  read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/columbus/Email_Name_Dept_Gender_2011.csv", header= T, stringsAsFactors = F)
document_edge_matrix = read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/columbus/edge-matrix.csv", header= F, stringsAsFactors = F)
document_edge_matrix = document_edge_matrix[,-1]
document_edge_matrix[,1] <- document_edge_matrix[,1] + 1 #make sure that authors are indexed starting at 1
document_word_matrix = read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/columbus/word-matrix.csv", header= F, stringsAsFactors = F)
document_word_matrix = document_word_matrix[,-1]
vocabulary = read.csv("~/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/columbus/vocab.txt", header= F, stringsAsFactors = F)

remove <- which(apply(document_word_matrix,1,sum) == 0)
document_edge_matrix <- document_edge_matrix[-remove,]
document_word_matrix <- document_word_matrix[-remove,]
save(list = ls(),file = "Columbus_2011_Data.Rdata")
#===================================#