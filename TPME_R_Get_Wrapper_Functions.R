#Function that wraps taking a log multinomial draw from a given distribution

log_multinomial_draw <- function(probability_vector){
    #print(probability_vector)
    edge_selected <-which(rmultinom(1,1,exp(probability_vector)) == 1)
    #print(edge_selected)
    return(edge_selected)
}

#log_multinomial_draw(rep(-.25,10))

get_observed_edge_value <- function(document,recipient){
    return(Document_Edge_Matrix[document,recipient])
}


#get_observed_edge_value(1,2)

get_token_topic_assignment <- function(document,token){
    return(Token_Topic_Assignments[[document]][token])
}

#get_token_topic_assignment(1,1)

get_sum_token_topic_assignments <- function(document,token,topic){
    return(length(which(Token_Topic_Assignments[[document]][-token] == topic)))
}

#get_sum_token_topic_assignments(1,1,1)

get_edge_topic_assignment <- function(document,recipient){
    return(Edge_Topic_Assignments[document,recipient])
}

get_number_of_topics <- function(){
    return(Number_Of_Topics)
}

#get_number_of_topics()

get_number_of_unique_words <- function(){
    return(Number_Of_Words)
}

#get_number_of_unique_words()

get_word_type_topic_assignemnt_count <- function(document,token,topic){
    temp <- Word_Type_Topic_Counts[Token_Word_Types[[document]][token],topic]
    if(Token_Topic_Assignments[[document]][token] == topic){
        temp <- temp -1
    }
    return(temp)
}

#get_word_type_topic_assignemnt_count(1,1,1)

get_number_of_tokens_assigned_to_topic <- function(document,token,topic){
    temp <- sum(Word_Type_Topic_Counts[,topic])
    if(Token_Topic_Assignments[[document]][token] == topic){
        temp <- temp -1
    }
    return(temp)
}

#get_number_of_tokens_assigned_to_topic(1,1,1)

get_new_intercept_sample <- function(topic,variance){
    return()
}

