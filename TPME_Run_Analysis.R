# Defines function which runs main analysis
# gloabl variables have first letter of every word capitalized
# local variables are not capitalized
# user defined functions are in all caps

Run_Analysis <- function(Number_Of_Iterations = 1000, Base_Alpha =1, Base_Beta = 0.01, Number_Of_Topics = 50, Author_Attributes, Document_Edge_Matrix ,Document_Word_Matrix, Vocabulary, Latent_Dimensions){
    
    #================ set working driectory and source all functions ====================#
    setwd("~/Dropbox/PINLab/Projects/R_Code/TPMNE")
    source("TPME_Sample_Token_Topic_Assignments.R")
    source("TPME_Sample_Edge_Topic_Assignments.R")
    source("TPME_Sample_Author_Topic_Latent_Space.R")
    source("TPME_Sample_Topic_Latent_Space_Intercept.R")
    source("TPME_Get_Probability_of_Edge.R")
    
    #================= Initialize all variables, latent spaces edge assingments and topic assignments ==============#
    
    Latent_Space_Intercepts <- rep(10, Number_Of_Topics) #this is set to 10 becasue it can only get smaller
    
    Number_Of_Documents <- length(Document_Word_Matrix[,1]) # the number of documents is equal to the number of rows 
                                  
    Metropolis_Hastings_Control_Parameter <- 0 #this is used to shrink the proposal variace of the metropolis hastings portion of the algorithm 
    
    Number_Of_Authors <- length(Author_Attributes[,1]) 
    
    Number_Of_Words <- length(Vocabulary) #the number of unique words in the corpus
    
    Beta <- Base_Beta*Number_Of_Words 
    
    #we define alpha to be a vector so that it can accomodate an asymmetric base measure in the future
    Alpha_Base_Measure_Vector <- rep((Base_Alpha/Number_Of_Topics,Number_Of_Topics)
                                     
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
        latent_space_positions[,,a] <- matrix(0,nrow = Latent_Dimensions, ncol = Number_Of_Topics)
        for(s in 1:Latent_Dimensions){
            latent_space_positions[s,,a] <- runif(Number_Of_Topics, min = -1, max = 1) #samples from a continuous uniform distribution on (-1,1)
        }
        #now add to the list object
        Latent_Space_Positions <- append(Latent_Space_Positions,list(latent_space_positions))
    }
    
    #store information about current edge log likelihoods. This is an author by author by number of topics array of lists of edge information including author latent coordinates, recipient coordinates, intercept, edge value (whether it was set to 1 or 0 (also known as y)),edge log likelihood anmd number of latent dimensions for convenience. 
    Current_Edge_Information <- array(list(rep(0,Latent_Dimensions),rep(0,Latent_Dimensions),10,0,0,Latent_Dimensions),c(Number_Of_Topics,Number_Of_Authors,Number_Of_Authors)
    
    #store information about proposed edge log likelihoods. This is an author by author by number of topics array of lists of edge information including author latent coordinates, recipient coordinates, intercept, edge value (whether it was set to 1 or 0 (also known as y)), edge log likelihood anmd number of latent dimensions for convenience.
    Proposed_Edge_Information <- array(list(rep(0,Latent_Dimensions),rep(0,Latent_Dimensions),10,0,0,Latent_Dimensions),c(Number_Of_Topics,Number_Of_Authors,Number_Of_Authors)
    
            
    
    #initialize edge topic assignments. this is a matrix that indexes documents by rows and the first column is the sender number and then there is one column for ever possible sender after that with zeros indicating the message was not sent to them and 1 indicating that it was sent to them. 
    Edge_Topic_Assignments <- Document_Edge_Matrix #jsut assing it so we get the right dimensions
    #now go in and replace all ones with a sampled edge topic assignment
    for(d in 1:Number_Of_Documents){
        for(a in 1:Number_Of_Authors){
                Edge_Topic_Assignments[d,] <- sample(1:Number_Of_Topics,Number_Of_Authors,replace= TRUE) #give it a topic edge assignment value
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
            if(sum(Document_Word_Matrix[d,]) > 0){ #if there is atleast one token in the document
                Token_Topic_Assignments[[d]][[1]] <- SAMPLE_TOKEN_TOPIC_ASSIGNMENTS(Token_Topic_Assignments[[d]][[1]],Number_Of_Topics,Number_Of_Authors,Alpha_Base_Measure_Vector,Beta,Number_Of_Words,Edge_Topic_Assignments[d,],Latent_Space_Positions,Latent_Space_Intercepts,Latent_Dimensions,Document_Authors[d],Document_Edge_Matrix)
            }else{ #assign the docuemnt a dummy word and dummy topic
                
            }
            
            
        }
        
        #3. Sample document edge assingments 
        for(d in 1:Number_Of_Documents){
            Edge_Topic_Assignments[d,] <- SAMPLE_DOCUMENT_EDGE_ASSIGNMENTS(Edge_Topic_Assignments[d,],Latent_Space_Positions,Latent_Space_Intercepts,Latent_Dimensions,Document_Authors[d],Document_Edge_Matrix,sum(Document_Word_Matrix[d,]))
            
        }
        
        #4. Sample latent positions for each actor and topic
        for(t in 1:Number_Of_Topics){
            for(a in 1:Number_Of_Authors){ #Author_Attributes is a matrix with one row per author and one column per attribute
                Latent_Space_Positions[,t,a] <- SAMPLE_NEW_LATENT_SPACE_POSITION_FOR_CURRENT_ACTOR_AND_TOPIC(Edge_Topic_Assignments,Latent_Space_Positions[,t,a],Latent_Space_Intercepts,Latent_Dimensions,a,Document_Edge_Matrix)
                
            }
        }
        
        #5. Sample new intercept for current topic
        for(t in 1:Number_Of_Topics){
            Latent_Space_Intercepts <- SAMPLE_NEW_INTERCEPT_FOR_CURRENT_TOPIC(Edge_Topic_Assignments,Latent_Space_Positions,Latent_Space_Intercepts,Latent_Dimensions,Document_Edge_Matrix,t)
        }
        
        
        
        
    }#end of main loop over number of itterations                     
    
    
    #get things ready to return a model object with all of the relevant info
    
    #save everything

    
} # End of Run_Analysis definition

# test