# Defines function which runs main analysis
# gloabl variables have first letter of every word capitalized
# local variables are not capitalized
# user defined functions are in all caps

#======== TESTING =========#
rm(list = ls())

Number_Of_Iterations = 1
Base_Alpha =1
Base_Beta = 0.01
Number_Of_Topics = 50
Author_Attributes = matrix(1:17,ncol =2,nrow =17)
Document_Edge_Matrix = read.csv("/Users/matthewjdenny/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/mcdowell/edge-matrix.csv", header= F, stringsAsFactors = F)
Document_Edge_Matrix = Document_Edge_Matrix[,-1]
Document_Edge_Matrix[,1] <- Document_Edge_Matrix[,1] + 1 #make sure that authors are indexed starting at 1
Document_Word_Matrix = read.csv("/Users/matthewjdenny/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/mcdowell/Train_Word_Matrix.csv", header= F, stringsAsFactors = F)
Document_Word_Matrix = Document_Word_Matrix[,-1]
Vocabulary = read.csv("/Users/matthewjdenny/Dropbox/PINLab/Projects/Denny_Working_Directory/Remove_Names_2011/mcdowell/vocab.txt", header= F, stringsAsFactors = F)
Latent_Dimensions = 2

#==========================#


Run_Analysis <- function(Number_Of_Iterations = 1000, Base_Alpha =1, Base_Beta = 0.01, Number_Of_Topics = 50, Author_Attributes, Document_Edge_Matrix ,Document_Word_Matrix, Vocabulary, Latent_Dimensions){
    
    #================ set working driectory and source all functions ====================#
    setwd("~/Dropbox/PINLab/Projects/R_Code/TPMNE")
    require(Rcpp)
    require(doMC)
    require(foreach)
    registerDoMC(8)
    #Rcpp::sourceCpp("TPME_Sample_Token_Topic_Assignments.cpp")
    Rcpp::sourceCpp("TPME_Sample_Token_Topic_Assignments_Full_Cpp.cpp")
    Rcpp::sourceCpp("TPME_Sample_Single_Token_Topic_Assignment_Full_Cpp.cpp")
    Rcpp::sourceCpp("TPME_Sample_Edge_Topic_Assignments.cpp")
    source("TPME_Sample_Author_Topic_Latent_Space.R")
    #Rcpp::sourceCpp("TPME_Sample_Latent_Space_Intercepts.cpp")
    source("TPME_Get_Probability_of_Edge.R")
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
   
   
   
   
    
    #==================== MAIN LOOP OVER NUMER OF ITTERATIONS ====================#                              
    for(i in 1:Number_Of_Iterations){
        
        #1. Set proposal variance for current itteration for metropolis hastings step
        Metropolis_Hastings_Control_Parameter <- Metropolis_Hastings_Control_Parameter + 1
        if(Metropolis_Hastings_Control_Parameter < 101){
            Proposal_Variance <- (100/Metropolis_Hastings_Control_Parameter) # this shrinks down the proposal variance to 1 as we reach the 100th itteration
        }
        
        #2. Sample token topic assignments
        for(d in 1:Number_Of_Documents){
            if(sum(Document_Word_Matrix[d,]) > 1){ #if there is atleast 2 tokens in the document
                #Token_Topic_Assignments[[d]] <- 
                SAMPLE_TOKEN_TOPIC_ASSIGNMENTS_CPP(length(Token_Topic_Assignments[[d]]),Number_Of_Topics,Number_Of_Authors,Document_Authors[d],d,Beta,Alpha_Base_Measure_Vector,Latent_Space_Positions[1,,],Latent_Space_Positions[2,,],Latent_Space_Intercepts,Latent_Dimensions, as.numeric(Edge_Topic_Assignments[d,]),Token_Topic_Assignments[[d]],apply(Word_Type_Topic_Counts,2,sum),Word_Type_Topic_Counts[Token_Word_Types[[d]],],Number_Of_Words,as.numeric(Document_Edge_Matrix[d,]))
                
                #update data structures
            }else{ 
                if(sum(Document_Word_Matrix[d,]) > 0){ #if there is one token in the document
                #print(paste("there was one token in document",d)) 
                SAMPLE_SINGLE_TOKEN_TOPIC_ASSIGNMENT_CPP(length(Token_Topic_Assignments[[d]]),Number_Of_Topics,Number_Of_Authors,Document_Authors[d],d,Beta,Alpha_Base_Measure_Vector,Latent_Space_Positions[1,,],Latent_Space_Positions[2,,],Latent_Space_Intercepts,Latent_Dimensions, as.numeric(Edge_Topic_Assignments[d,]),Token_Topic_Assignments[[d]],apply(Word_Type_Topic_Counts,2,sum),Word_Type_Topic_Counts[Token_Word_Types[[d]],],Number_Of_Words,as.numeric(Document_Edge_Matrix[d,]))
                #update data structures
                
                }else{
                    #assign the docuemnt a dummy word and dummy topic
                    #print(paste("there were no tokens in document",d)) 
                }   
            } 
        }
        
        
        #3. Sample document edge assingments 
        foreach(d=1:Number_Of_Documents) %dopar% {
            Edge_Topic_Assignments[d,] <- SAMPLE_EDGE_TOPIC_ASSIGNMENTS_CPP(Number_Of_Authors,Document_Authors[d],length(Token_Topic_Assignments[[d]]),d)
            
        }
        
        #4. Sample latent positions for each actor and topic
        for(t in 1:Number_Of_Topics){
            for(a in 1:Number_Of_Authors){ #Author_Attributes is a matrix with one row per author and one column per attribute
                Latent_Space_Positions[,t,a] <- SAMPLE_NEW_LATENT_SPACE_POSITION_FOR_CURRENT_ACTOR_AND_TOPIC(Edge_Topic_Assignments,Latent_Space_Positions[,t,a],Latent_Space_Intercepts,Latent_Dimensions,a,Document_Edge_Matrix)
                
            }
        }
        
        #5. Sample new intercept for current topic and perform a metropolis step
        foreach(t=1:Number_Of_Topics) %dopar% {
            Latent_Space_Intercepts[t] <- SAMPLE_NEW_TOPIC_INTERCEPT_CPP(t,Latent_Dimensions,Number_Of_Authors,Proposal_Variance)
        }
        
        
        
        
    }#end of main loop over number of itterations                     
    
    
    #get things ready to return a model object with all of the relevant info
    
    #save everything

    
} # End of Run_Analysis definition

# test