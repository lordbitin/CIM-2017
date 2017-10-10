#INPUT: List with timestamps of each observation (MILLISECONDS)
#OUTPUT
hmmHelper.timeResolution <- function(timingData) {
  min(diff(timingData))
}

# INPUT
#   - observationDf: Dataframe with multiple (and timestamped) observation sequences
#   - idColname: Name of the column containing the identification of each observation sequence
#   - timingColname: Timestamp column name
# OUTPUT
#   - Time resolution : Minimum difference between two observation symbols in the training data (numeric)
hmmHelper.timeResolutionForMultipleSequences <- function (observationDf,idColname,timingColname) {
  trss <- observationDf %>%
    group_by_(idColname) %>%
    summarise_each(funs(hmmHelper.timeResolution),matches(timingColname))
  
  min(trss[timingColname])
}

# INPUT 
#   - model : class "hmm" or "hsmm" from mhsmm package
#   - state : Number of the state to treat
#
# OUPUT : List of ranked observation symbols
# NOTE: Only valid for categorical HMMs or HSMMs
hmmHelper.getObservableRanking <- function (model,state, format = "vector") {
  
  symbolRanking <- NULL
  
  switch(format,
         vector={
           sortedProbs <- sort.int(model$model$parms.emission$prob[[state]], 
                                   decreasing=TRUE, 
                                   index.return=TRUE)
           symbolRanking <- lapply(sortedProbs$ix, function (obsSymbolIndex) obsSymbolIndex-1)
           },
         matrix={
           nSymbols <- length(model$model$parms.emission$prob[[1]])
           symbolRanking <- llply(1:nSymbols,function (x) rep(0,nSymbols)) #Initialize
           sortedProbs <- sort.int(model$model$parms.emission$prob[[state]], 
                                   decreasing=TRUE, 
                                   index.return=TRUE)
           
           for (i in 1:length(sortedProbs$ix)) {
             symbolRanking[[i]][sortedProbs$ix[i]] <- 1
           }           
         })
  
  return(symbolRanking)
}

#
hmmHelper.randomInitialVector <- function (J) {
  aux <- runif(J)
  aux/sum(aux)
}

# INPUT: J (number of states)
# OUTPUT : random transition matrix for an HMM
hmmHelper.randomTransitionMatrix.hmm <- function (J) {
  if (J<2)
    stop("J must be >= 2")
  
  auxL <- lapply(1:J,function (i) hmmHelper.randomInitialVector(J))
  A0 <- do.call(rbind,auxL)
  A0
}

# INPUT: J (number of states)
# OUTPUT : random transition matrix for an HSMM
hmmHelper.randomTransitionMatrix.hsmm <- function (J) {
  
  if (J<2)
    stop("J must be >= 2")
  
  m <- matrix(rep(0,J*J),nrow=J)

  for (i in 1:J) {
    randRow <- hmmHelper.randomInitialVector(J-1)
    
    if (i>1)
      m[i,1:(i-1)] <- randRow[1:(i-1)]
    
    if (i<J)
      m[i,(i+1):J] <- randRow[i:(J-1)]
  }
  
  m
}

##
# METHOD: hmmHelper.emissionMatrix.seq
##
# Input:  stsList: Object of class stsList from package TraMineR. 
#         J: Integer. Number of states for the HMM
#         
# Output:
#         Matrix containing the emission distributions of each state based on the frequencies of 
#         each alphabet symbol in stsList
#
hmmHelper.emissionMatrix.seq <- function (stsList,J) {
  result <- laply(helper.splitIntoEqualChunks(seq(1,max(seqlength(stsList))),J),
                                    function (subseqsIndex) {
                                      seqstatf(stsList[,subseqsIndex])[,2] + 0.1
                                    })
  result <- result/rowSums(result)  
  result
}