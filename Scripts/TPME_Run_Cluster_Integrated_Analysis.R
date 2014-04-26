# Defines function which runs main analysis
# gloabl variables have first letter of every word capitalized
# local variables are not capitalized
# user defined functions are in all caps


Run_Cluster_Integrated_Analysis <- function(Number_Of_Iterations = 1000, Base_Alpha =1, Base_Beta = 0.01, Number_Of_Topics = 50, Author_Attributes= author_attributes, Document_Edge_Matrix = document_edge_matrix ,Document_Word_Matrix = document_word_matrix, Vocabulary = vocabulary, Latent_Dimensions = 2, Topic_Step_Itterations = 1, Sample_Step_Itterations = 10, output_file = "Test",Proposal_Variance = 0.5, seed = 1234, output_folder_path = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", system_OS = "Linux", Number_of_Clusters = 10,Itterations_Before_Cluster_Assingment_Updates = 5, Adaptive_Metropolis_Target_Accept_Rate = 0.3){
    
    #================ set working driectory and source all functions ====================#
    require(Rcpp)
    require(RcppArmadillo)
    set.seed(seed)
    
    #if we are running linux then we need to add the appropriate c flags to use c++2011
    if(system_OS == "Linux"){
        #PKG_CPPFLAGS = "-std=c++11"
        PKG_CPPFLAGS = "-std=c++0x"
        Sys.setenv(PKG_CPPFLAGS = PKG_CPPFLAGS)
    }
    Rcpp::sourceCpp("./Scripts/TPME_Cluster_Integrated_Sampler.cpp")
    print("Source Files Loaded...")
    
    #================= Initialize all variables, latent spaces edge assingments and topic assignments ==============#
    
    Latent_Space_Intercepts <- rep(10, Number_of_Clusters) #this is set to 10 becasue it can only get smaller
    
    temp <- Proposal_Variance #create vector of proposal variances
    Proposal_Variance <- rep(temp, Number_of_Clusters)
    
    Number_Of_Documents <- length(Document_Word_Matrix[,1]) # the number of documents is equal to the number of rows 
    
    Metropolis_Hastings_Control_Parameter <- 0 #this is used to shrink the proposal variace of the metropolis hastings portion of the algorithm 
    
    Number_Of_Authors <- length(Author_Attributes[,1]) 
    
    Number_Of_Words <- length(Vocabulary[,1]) #the number of unique words in the corpus
    
    Beta <- Base_Beta*Number_Of_Words 
    
    #we define alpha to be a vector so that it can accomodate an asymmetric base measure in the future
    Alpha_Base_Measure_Vector <- rep(Base_Alpha/Number_Of_Topics,Number_Of_Topics)
    
    Document_Authors <- Document_Edge_Matrix[,1] #make a vector of document authors
    
    Document_Edge_Matrix <- Document_Edge_Matrix[,-1] # remove the authors from the docuemnt edge matrix
    print("Initializing Topic Assignments...")
    #token topic assignemnts are stores in a list of vectors data structure
    Token_Topic_Assignments <- list()
    for(d in 1:Number_Of_Documents){
        #allocate a vector of zeros equal to the number of tokens in the document
        #cur_token_assignments <- rep(0,sum(Document_Word_Matrix[d,])) #assign a zero vector of topic assignments if we then wanted to do a more complicated sampling proceedure 
        cur_token_assignments <- sample(1:Number_Of_Topics,sum(Document_Word_Matrix[d,]),replace= T) #samples from a discrete uniform distribution
        Token_Topic_Assignments <- append(Token_Topic_Assignments,list(cur_token_assignments))
    }
    print("Initiailizing Latent Space Positions...")
    #sample latent space postions for each actor for each topic from a uniform distribution. This will be a list of matricies data structure with rows in each matrix being the latent dimensions and columns being each topic 
    Latent_Space_Positions <- array(0,c(Latent_Dimensions,Number_of_Clusters,Number_Of_Authors))
    for(a in 1:Number_Of_Authors){ 
        Latent_Space_Positions[,,a] <- matrix(0,nrow = Latent_Dimensions, ncol = Number_of_Clusters)
        for(s in 1:Latent_Dimensions){
            Latent_Space_Positions[s,,a] <- runif(Number_of_Clusters, min = -1, max = 1) #samples from a continuous uniform distribution on (-1,1)
        }  
    }
    
    #now make them all the same
    for(c in 1:Number_of_Clusters){
        Latent_Space_Positions[,c,] <- Latent_Space_Positions[,1,]
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
    print("Initializing Edge Topic Assignments...")
    Edge_Topic_Assignments <- Document_Edge_Matrix #jsut copying it so we get the right dimensions
    #now go in and replace all ones with a sampled edge topic assignment
    for(d in 1:Number_Of_Documents){
        for(a in 1:Number_Of_Authors){
            Edge_Topic_Assignments[d,] <- sample(1:Number_Of_Topics,Number_Of_Authors,replace= TRUE) #give it a topic edge assignment value
        }
    }
    
    #initialize a datastructure to keep a number of topics by number of unique words matrix 
    print("Initializing Word Type Topic Counts...")
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
    print("Initializing Betas...")
    #initialize betas 
    Number_of_Betas <- 4
    Betas <- matrix(runif(Number_of_Clusters*Number_of_Betas),nrow =Number_of_Clusters,ncol = Number_of_Betas)
    
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
    
    
    #Assign topics to clusters
    Topic_Cluster_Assignments <- rep(0,Number_Of_Topics)
    for(k in 1:Number_Of_Topics){
        Topic_Cluster_Assignments[k] <- round(runif(1, min = 1, max = Number_of_Clusters),0)
    }
    
    #==================== MAIN Function ====================#                             
        
        Return_List <- Cluster_Integrated_Sampler(
            Number_Of_Iterations,
            Topic_Step_Itterations,
            Sample_Step_Itterations,
            Number_Of_Authors, 
            Number_Of_Topics,
            Number_of_Clusters,
            Latent_Dimensions,
            Number_Of_Documents,
            Proposal_Variance,
            Topic_Cluster_Assignments,
            array(0,c(Number_Of_Authors,Number_Of_Authors,Number_Of_Topics)),
            array(0,c(Number_Of_Authors,Number_Of_Authors,Number_Of_Topics)),
            Latent_Space_Positions, 
            array(0,c(Latent_Dimensions,Number_Of_Topics,Number_Of_Authors)),
            Latent_Space_Intercepts,
            Betas,
            Number_of_Betas,
            Beta_Indicator_Array,
            as.matrix(Document_Edge_Matrix),
            Token_Topic_Assignments,
            Token_Word_Types,
            Document_Authors,
            Beta,
            Alpha_Base_Measure_Vector,
            as.matrix(Edge_Topic_Assignments),
            Word_Type_Topic_Counts,
            apply(Word_Type_Topic_Counts,2,sum),
            Number_Of_Words,
            Itterations_Before_Cluster_Assingment_Updates,
            Adaptive_Metropolis_Target_Accept_Rate
        )

        #get things ready to return a model object with all of the relevant info 
        save(Return_List, file = paste(output_folder_path,"Model_Output_",output_file,".Rdata",sep = ""))
    
    
    return(Return_List)
} # End of Run_Analysis definition




