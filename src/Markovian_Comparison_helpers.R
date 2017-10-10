#
# De 5,2,0 a list(M=5,l=2,f=0)
#
Markovian_Comparison_helpers.extractHyperParamsFromRanking <- function(rank, num = length(rank)) {
  bestHyperparams <- (rank %>% strsplit(split = ","))[1:num] %>% 
    lapply(as.numeric) %>% 
    lapply(function (x) setNames(x,c("M","l","f"))) %>% 
    lapply(as.list)
           
 bestHyperparams
}

#
# Sca la posicion de unos hiper parametros dentro del ranking
#
Markovian_Comparison_helpers.rankPosition <- function(rank,M,l,f) {
  rank_ <- Markovian_Comparison_helpers.extractHyperParamsFromRanking(rank)
  rank_df <- matrix(unlist(rank_),ncol=3,byrow = T) %>% as.data.frame %>% setNames(c("M","l","f"))
  
  rankPosition <- which(rank_df$M == M & rank_df$l == l & rank_df$f == f)
  ifelse(length(rankPosition != 0),rankPosition,NA)
}


##
# Function to convert the grid search data frame to a rank aggregation format
# Input   gs: A grid search object for DCMM, each column represets a performance measure
##
Markovian_Comparison_helpers.getRanksAndWeights <- function (gs, measures = c("ll","BIC","NRS","SLL","AGP","MPGP","CHPHT","CHPVT")) {
  
  #Error control
  stopifnot(all(measures %in% names(gs)))
  
  aux <- lapply(X = measures, FUN = function (meas) {
    theHighertheBetter <- meas %in% c("ll","SLL","AGP", "MPGP","CHPHT","CHPVT")
    
    sortedGs <- gs %>% arrange_(meas)
    if (theHighertheBetter) sortedGs <- sortedGs %>% slice(nrow(sortedGs):1)
    sortedGs <- sortedGs %>% unite(col = hyperParams, M, orderHC, orderVC, sep=",")
    #sortedGs$hyperParams <- paste("$\\mu(",sortedGs$hyperParams,")$")
    
    list(
      ranks = sortedGs$hyperParams,
      weights = sortedGs[[meas]]
    )
  }) %>% setNames(measures)
  
  result <- list(
    ranks = bind_rows(lapply(aux,"[[","ranks")) %>% t,
    weights = bind_rows(lapply(aux,"[[","weights")) %>% t
  )
  result$weights_normalized <- result$weights %>% apply(1, function (row) (row - min(row))/(max(row) - min(row))) %>% t
  
  result
}

##
# Coefficient of variation (Cv), also called Relative Standard Deviation (RSD)
##
Markovian_Comparison_helpers.coefficient_of_variation <- function (x) {
  return(
    sd(x)/abs(mean(x))
  )
}

Markovian_Comparison_helpers.quartile_coefficient_of_dispersion <- function (x) {
  q1 <- quantile(x,probs = 0.25)
  q3 <- quantile(x,probs = 0.75)
  return(
    (q3-q1)/(q3+q1)
  )
}