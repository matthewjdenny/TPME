#Function that wraps taking a log multinomial draw from a given distribution

log_multinomial_draw <- function(probability_vector){
    return(which(rmultinom(1,1,exp(probability_vector)) == 1))
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
    return(length(which(Token_Topic_Assignments[[document]] == topic)))
}

#get_sum_token_topic_assignments(1,1,1)

get_edge_topic_assignment <- function(document,recipient){
    return(Edge_Topic_Assignments[document,recipient])
}
