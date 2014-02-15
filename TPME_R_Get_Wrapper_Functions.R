#Function that wraps taking a log multinomial draw from a given distribution

log_multinomial_draw <- function(probability_vector){
    return(which(rmultinom(1,1,exp(probability_vector)) == 1))
}

get_observed_edge_value <- function(document,recipient){
    return(Document_Edge_Matrix[document,recipient])
}

get_token_topic_assignment <- function(document,token){
    return(Token_Topic_Assignments[[document][[token]]])
}

get_edge_topic_assignment <- function(document,recipient){
    return(Edge_Topic_Assignments[document,recipient])
}
