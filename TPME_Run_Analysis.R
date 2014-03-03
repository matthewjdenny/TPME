# Defines function which runs main analysis
# gloabl variables have first letter of every word capitalized
# local variables are not capitalized
# user defined functions are in all caps

#======== TESTING =========#
rm(list = ls())

#Number_Of_Iterations = 1
#Base_Alpha =1
#Base_Beta = 0.01
#Number_Of_Topics = 50
author_attributes =  read.csv("/Users/matthewjdenny/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/mcdowell/Email_Name_Dept_Gender_2011.csv", header= T, stringsAsFactors = F)
#making a judgement call that Robbin silvers in a woman
author_attributes[14,4] <- F

document_edge_matrix = read.csv("/Users/matthewjdenny/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/mcdowell/edge-matrix.csv", header= F, stringsAsFactors = F)
document_edge_matrix = Document_Edge_Matrix[,-1]
document_edge_matrix[,1] <- Document_Edge_Matrix[,1] + 1 #make sure that authors are indexed starting at 1
document_word_matrix = read.csv("/Users/matthewjdenny/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/mcdowell/word-matrix.csv", header= F, stringsAsFactors = F)
document_word_matrix = Document_Word_Matrix[,-1]
vocabulary = read.csv("/Users/matthewjdenny/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/mcdowell/vocab.txt", header= F, stringsAsFactors = F)
#Latent_Dimensions = 2

#==========================#
}

Run_Analysis <- function(Number_Of_Iterations = 1, Base_Alpha =1, Base_Beta = 0.01, Number_Of_Topics = 50, Author_Attributes= author_attributes, Document_Edge_Matrix = document_edge_matrix ,Document_Word_Matrix = document_word_matrix, Vocabulary = vocabulary, Latent_Dimensions = 2){
    
    #================ set working driectory and source all functions ====================#
    setwd("~/Dropbox/PINLab/Projects/R_Code/TPMNE")
    require(Rcpp)
    #Rcpp::sourceCpp("TPME_Sample_Token_Topic_Assignments.cpp")
    #Rcpp::sourceCpp("TPME_Sample_Token_Topic_Assignments_Full_Cpp.cpp")
    #Rcpp::sourceCpp("TPME_Sample_Single_Token_Topic_Assignment_Full_Cpp.cpp")
    #Rcpp::sourceCpp("TPME_Sample_Edge_Topic_Assignments_Full_Cpp.cpp")
    #source("TPME_Sample_Author_Topic_Latent_Space.R")
    #source("TPME_Get_Probability_of_Edge.R")
    Rcpp::sourceCpp("TPME_Metropolis_Step.cpp")
    Rcpp::sourceCpp("TPME_Topic_assignment_Step.cpp")
    source("TPME_R_Get_Wrapper_Functions.R")
    
    #================= Initialize all variables, latent spaces edge assingments and topic assignments ==============#
    
    Latent_Space_Intercepts <- rep(10, Number_Of_Topics) #this is set to 10 becasue it can only get smaller
    
    Number_Of_Documents <- length(Document_Word_Matrix[,1]) # the number of documents is equal to the number of rows 
                                  
    Metropolis_Hastings_Control_Parameter <- 0 #this is used to shrink the proposal variace of the metropolis hastings portion of the algorithm 
    
    Number_Of_Authors <- length(Author_Attributes[,1]) 
    
    Number_Of_Words <- length(Vocabulary[,1]) #the number of unique words in the corpus
    
    Beta <- Base_Beta*Number_Of_Words 
    
    #we define alpha to be a vector so that it can accomodate an asymmetric base measure in the future
    Alpha_Base_Measure_Vector <- rep(Base_Alpha/Number_Of_Topics,Number_Of_Topics)
                                     
    Document_Authors <- Document_Edge_Matrix[,1] #make a vector of document authors
    
    Document_Edge_Matrix <- Document_Edge_Matrix[,-1] # remove the authors from the docuemnt edge matrix
    
    #token topic assignemnts are stores in a list of vectors data structure
    Token_Topic_Assignments <- list()
    for(d in 1:Number_Of_Documents){
        #allocate a vector of zeros equal to the number of tokens in the document
        #cur_token_assignments <- rep(0,sum(Document_Word_Matrix[d,])) #assign a zero vector of topic assignments if we then wanted to do a more complicated sampling proceedure 
        cur_token_assignments <- sample(1:Number_Of_Topics,sum(Document_Word_Matrix[d,]),replace= T) #samples from a discrete uniform distribution
        Token_Topic_Assignments <- append(Token_Topic_Assignments,list(cur_token_assignments))
    }
    
    #sample latent space postions for each actor for each topic from a uniform distribution. This will be a list of matricies data structure with rows in each matrix being the latent dimensions and columns being each topic 
    Latent_Space_Positions <- array(0,c(Latent_Dimensions,Number_Of_Topics,Number_Of_Authors))
    for(a in 1:Number_Of_Authors){ 
        Latent_Space_Positions[,,a] <- matrix(0,nrow = Latent_Dimensions, ncol = Number_Of_Topics)
        for(s in 1:Latent_Dimensions){
            Latent_Space_Positions[s,,a] <- runif(Number_Of_Topics, min = -1, max = 1) #samples from a continuous uniform distribution on (-1,1)
        }
        
    }
    
    #store information about current edge log likelihoods. This is a number of topics to number of authors to number of recipients list of lists containing a list of edge information including author latent coordinates, recipient coordinates, intercept, edge value (whether it was set to 1 or 0 (also known as y)),edge log likelihood anmd number of latent dimensions for convenience. 
    test <-   list(rep(0,Latent_Dimensions),rep(0,Latent_Dimensions),10,0,0,Latent_Dimensions)       
    test2 <- list()
    for(i in 1:Number_Of_Authors){
       test2 = append(test2,list(test))
    }
    test3 <- list()
    for(i in 1:Number_Of_Authors){
       test3 = append(test3,list(test2))
    }
    test4 <- list()
    for(i in 1:Number_Of_Topics){
       test4 = append(test4,list(test3))
    }
    Current_Edge_Information <- test4
    Proposed_Edge_Information <- Current_Edge_Information
    #initialize edge topic assignments. this is a matrix that indexes documents by rows and the first column is the sender number and then there is one column for ever possible sender after that with zeros indicating the message was not sent to them and 1 indicating that it was sent to them. 
    Edge_Topic_Assignments <- Document_Edge_Matrix #jsut copying it so we get the right dimensions
    #now go in and replace all ones with a sampled edge topic assignment
    for(d in 1:Number_Of_Documents){
        for(a in 1:Number_Of_Authors){
                Edge_Topic_Assignments[d,] <- sample(1:Number_Of_Topics,Number_Of_Authors,replace= TRUE) #give it a topic edge assignment value
        }
    }
   
   #initialize a datastructure to keep a number of topics by number of unique words matrix 
   Word_Type_Topic_Counts <- matrix(0,nrow = Number_Of_Words, ncol = Number_Of_Topics)
   #this list keeps a vector for each document of word types for each token in the same order that topic assignemnts are kept
   Token_Word_Types <- list()
   for(d in 1:Number_Of_Documents){
       #create a vector of word types for each token
       word_indexes <- which(Document_Word_Matrix[d,] > 0)
       word_counts <- as.numeric(Document_Word_Matrix[d,word_indexes])
       already <- F
       for(i in 1:length(word_indexes)){
           if(!already){
               already <- T
               word_types <- rep(word_indexes[i],word_counts[i])
           }else{
               word_types <- c(word_types,rep(word_indexes[i],word_counts[i]))
           }
       }
       Token_Word_Types <- append(Token_Word_Types,list(word_types))
       #now get the token topic assignemnts for this document
       current_doc_assignments <- Token_Topic_Assignments[[d]]
       #now go through and increment based in intial draws
       for(i in 1:length(current_doc_assignments)){
           Word_Type_Topic_Counts[word_types[i],current_doc_assignments[i]] <- Word_Type_Topic_Counts[word_types[i],current_doc_assignments[i]] + 1
       }
    }
    
   #initialize betas 
    Number_of_Betas <- 4
    Betas <- matrix(runif(Number_Of_Topics*Number_of_Betas),nrow =Number_Of_Topics,ncol = Number_of_Betas)
   
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
    
   
   
    
    #==================== MAIN LOOP OVER NUMER OF ITTERATIONS ====================#                              
    for(i in 1:Number_Of_Iterations){
        
        #1. Set proposal variance for current itteration for metropolis hastings step
        Metropolis_Hastings_Control_Parameter <- Metropolis_Hastings_Control_Parameter + 1
        if(Metropolis_Hastings_Control_Parameter < 6){
            Proposal_Variance <- (3/Metropolis_Hastings_Control_Parameter) # this shrinks down the proposal variance to 1 as we reach the 100th itteration
        }
        
        #2. Sample token topic assignments
        
        Topic_Assignment_Results <- Topic_Assignment_Step_CPP(
            Number_Of_Authors, 
            Number_Of_Topics,
            array(0,c(Number_Of_Authors,Number_Of_Authors,Number_Of_Topics)),
            array(0,c(Number_Of_Authors,Number_Of_Authors,Number_Of_Topics)),
            Latent_Space_Positions, 
            Latent_Space_Intercepts,
            Latent_Dimensions,
            Betas,
            Number_of_Betas,
            Beta_Indicator_Array,
            1000,
            Proposal_Variance,
            Number_Of_Documents,
            as.matrix(Document_Edge_Matrix),
            Token_Topic_Assignments,
            Token_Word_Types,
            Document_Authors,
            Beta,
            Alpha_Base_Measure_Vector,
            as.matrix(Edge_Topic_Assignments),
            Word_Type_Topic_Counts,
            apply(Word_Type_Topic_Counts,2,sum),
            Number_Of_Words
            )
        
        
        #Assign Results
        Token_Topic_Assignments <- Topic_Assignment_Results[[1]]
        Topic_Present_Edge_Counts <- Topic_Assignment_Results[[2]]
        Topic_Absent_Edge_Counts <- Topic_Assignment_Results[[3]]
        #Test: (sum(Topic_Present_Edge_Counts) + sum(Topic_Absent_Edge_Counts)) == (Number_Of_Authors -1)* Number_Of_Documents
        Word_Type_Topic_Counts <- Topic_Assignment_Results[[4]]
        #Test: sum(Word_Type_Topic_Counts) == sum(unlist(lapply(Token_Topic_Assignments,length)))
        Edge_Topic_Assignments <- Topic_Assignment_Results[[5]]
        #sum(Edge_Topic_Assignments)
        
        
        #3. Perform MEtropolis step jointly for intercepts, latent positions and betas (not currently implemented)
        #calculate new latent positions, intercepts and betas
        
        Metropolis_Results <- Metropolis_Step_CPP(
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
            100000,
            Proposal_Variance
            )
        
        #assign metropolis results
        Latent_Space_Positions <- Metropolis_Results[[1000]]
        Latent_Space_Intercepts <- Metropolis_Results[[2000]]
        Betas <- Metropolis_Results[[3000]]
        #testing
        #Metropolis_Results[2001:2030]
        #sum(unlist(Metropolis_Results[30001:40000]))
        #intercepts <- rep(0,10000)
        #for(i in 1:10000){
        #    intercepts[i] <- Metropolis_Results[[i]][20]
        #}
        #plot(intercepts)
        
        
    }#end of main loop over number of itterations
   
   

    #get things ready to return a model object with all of the relevant info
   Return_List <- list(list(Metropolis_Results),list(Topic_Assignment_Results))
   
    
    #save everything

    return(Return_List)
} # End of Run_Analysis definition

# test