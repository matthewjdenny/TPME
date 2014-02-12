# Defines function which runs main analysis
# gloabl variables have first letter of every word capitalized
# local variables are not capitalized
# user defined functions are in all caps

Run_Analysis <- function(Number_Of_Iterations = 1000, Base_Alpha =1, Base_Beta = 0.01, Number_Of_Topics = 50, Author_Attributes, Document_Edge_Matrix ,Document_Word_Matrix, Vocabulary, Latent_Dimensions){
    
    #================ set working driectory and source all functions ====================#
    setwd("~/Dropbox/PINLab/Projects/R_Code/TPMNE")
    source("TPME_Sample_Token_Topic_Assignments.R")
    
    #================= Initialize all variables, latent spaces edge assingments and topic assignments ==============#
    
    Latent_Space_Intercepts <- rep(10, Number_Of_Topics) #this is set to 10 becasue it can only get smaller
    
    Number_Of_Documents <- length(Document_Word_Matrix[,1]) # the number of documents is equal to the number of rows 
                                  
    Metropolis_Hastings_Control_Parameter <- 0 #this is used to shrink the proposal variace of the metropolis hastings portion of the algorithm 
    
    Number_Of_Authors <- length(Author_Attributes[,1]) 
    
    Beta <- Base_Beta*length(Vocabulary)
    
    #we define alpha to be a vector so that it can accomodate an asymmetric base measure in the future
    Alpha_Base_Measure_Vector <- rep((Base_Alpha/Number_Of_Topics,Number_Of_Topics)
    
    #token topic assignemnts are stores in a list of vectors data structure
    Token_Topic_Assignments <- list()
    for(d in 1:Number_Of_Documents){
        #allocate a vector of zeros equal to the number of tokens in the document
        #cur_token_assignments <- rep(0,sum(Document_Word_Matrix[d,])) #assign a zero vector of topic assignments if we then wanted to do a more complicated sampling proceedure 
        cur_token_assignments <- sample(1:Number_Of_Topics,sum(Document_Word_Matrix[d,]),replace= T) #samples from a discrete uniform distribution
        Token_Topic_Assignments <- append(Token_Topic_Assignments,list(cur_token_assignments))
    }
    
    #sample latent space postions for each actor for each topic from a uniform distribution. This will be a list of matricies data structure with rows in each matrix being the latent dimensions and columns being each topic 
    Latent_Space_Positions <- list()
    for(a in 1:Number_Of_Authors){ 
        latent_space_positions <- matrix(0,nrow = Latent_Dimensions, ncol = Number_Of_Topics)
        for(s in 1:Latent_Dimensions){
            latent_space_positions[s,] <- runif(Number_Of_Topics, min = -1, max = 1) #samples from a continuous uniform distribution on (-1,1)
        }
        #now add to the list object
        Latent_Space_Positions <- append(Latent_Space_Positions,list(latent_space_positions))
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
                Token_Topic_Assignments[[d]][[1]] <- SAMPLE_TOKEN_TOPIC_ASSIGNMENTS(Token_Topic_Assignments[[d]][[1]],Number_Of_Topics,Number_Of_Authors,Alpha_Base_Measure_Vector, Beta)
            }else{ #assign the docuemnt a dummy word and dummy topic
                
            }
            
            
        }
        
        #3. Sample document edge assingments 
        for(d in 1:Number_Of_Documents){
            
            
        }
        
        #4. Sample latent positions for each actor and topic
        for(t in 1:Number_Of_Topics){
            for(a in 1:Number_Of_Authors){ #Author_Attributes is a matrix with one row per author and one column per attribute
                
                
            }
        }
        
        #5. Sample new intercept for current topic
        for(t in 1:Number_Of_Topics){
            
        }
        
        
        
        
    }#end of main loop over number of itterations                     
    
    
    #get things ready to return a model object with all of the relevant info
    
    #save everything

    
} # End of Run_Analysis definition

# test