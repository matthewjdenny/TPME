# Defines function which runs main analysis

Run_Analysis <- function(Number_Of_Iterations = 1000, Alpha =1, Beta = 0.01, Number_Of_Topics = 50, Author_Attributes, Document_Edge_Matrix ,Document_Word_Matrix, Vocabulary, Latent_Dimensions, Save_Interval){
    
    #==== Initialize all variables, latent spaces edge assingments and topic assignments ====#
    
    Latent_Space_Intercepts <- rep(10, Number_Of_Topics) #this is set to 10 becasue it can only get smaller
    NUMBER_OF_DOCUMENTS <- length(Document_Word_Matrix[,1] # the number of documents is equal to the number of rows 
    Metropolis_Hastings_Control_Parameter <- 0 #this is used to shrink the proposal variace of the metropolis hastings portion of the algorithm 
    
    
                                  
    #===== MAIN LOOP OVER NUMER OF ITTERATIONS =====#                              
    for(i in 1:Number_Of_Iterations){
        
        #1. Set proposal variance for current itteration for metropolis hastings step
        Metropolis_Hastings_Control_Parameter <- Metropolis_Hastings_Control_Parameter + 1
        if(Metropolis_Hastings_Control_Parameter < 101){
            Proposal_Variance <- (100/Metropolis_Hastings_Control_Parameter) # this shrinks down the proposal variance to 1 as we reach the 100th itteration
        }
        
        
        #2. Sample token topic assignments
        for(d in 1:NUMBER_OF_DOCUMENTS){
            if(sum(Document_Word_Matrix[d,]) > 0){ #if there is atleast one token in the document
                
            }else{ #assign the docuemnt a dummy word and dummy topic
                
            }
            
            
        }
        
        #3. Sample document edge assingments 
        for(d in 1:NUMBER_OF_DOCUMENTS){
            
            
        }
        
        #4. Sample latent positions for each actor and topic
        for(t in 1:Number_Of_Topics){
            for(a in length(Author_Attributes[,1])){ #Author_Attributes is a matrix with one row per author and one column per attribute
                
                
            }
        }
        
        #5. Sample new intercept for current topic
        for(t in 1:Number_Of_Topics){
            
        }
        
        #6. If we are at a save interval, save everything to an .RData object
        
        
    }#end of main loop over number of itterations                     
    
    
    #get things ready to return a model object with all of the relevant info
    
    #save everything

    
} # End of Run_Analysis definition

# test