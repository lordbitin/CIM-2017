#
# Validation metric (predictive capabilities metric) of an HSMM given an specific test dataset
#
validation.hsmm <- function (object, testData, predictionLengths = 2:min(testData$N), predictionDepth = 1:object$J, method = validation.RMSE.hsmm,...) {
  resultsDf <- data.frame(
    llply(predictionDepth, 
          function (pd) {
            
            rmses <- laply(predictionLengths, function (PL) {
              sequenceValidations <- laply(mhsmmHelper.getDataIndexes(testData), function (testSequenceIndexes) {
                method(object,
                                     testData$x[testSequenceIndexes],
                                     predictionLength = PL,
                                     predictionDepth = pd,
                                    ...)
              })
              
              #print(sequenceValidations)
              mean(sequenceValidations)
            })
            
            rmses
          })
  )
  colnames(resultsDf) <- predictionDepth
  rownames(resultsDf) <- predictionLengths

  resultsDf
}

#
# PREDICTION SCORE
#
validation.score.simple.hsmm <- function (object, testSequence, predictionLength = 10, predictionDepth = 1, C=1, observationWeights = rep(1/object$J,object$J), ...) {
  slidedTestSequence <- helper.slice(x = testSequence, n = predictionLength)
  slidedTestSequenceCumLengths <- cumsum(laply(slidedTestSequence,length))
  
  #
  # 1. Prediction Generation
  #
  currentMostLikelyPrediction <- NULL
  statePredictions <- numeric(length(testSequence))
  for (i in 1:length(slidedTestSequence)) {
    
    currentOwnPrediction <- .stateSequencePrediction.simple.hsmm(object,
                                                          lastPrediction = currentMostLikelyPrediction,
                                                          length = length(slidedTestSequence[[i]]))
    
    #Update the current most likely prediction via the viterbi algorithm, and get only the information 
    # about the current slide
    currentMostLikelyPrediction <- predict.hsmm(object = object,
                                                newdata = mhsmmHelper.formatMultinomialData(testSequence[1:slidedTestSequenceCumLengths[i]],
                                                                                            slidedTestSequenceCumLengths[i]),
                                                method = "viterbi")
    
    currentMostLikelyPrediction$s <- currentMostLikelyPrediction$s[(slidedTestSequenceCumLengths[i]-length(slidedTestSequence[[i]])+1):slidedTestSequenceCumLengths[i]]
    
    #Save own predictions
    statePredictions[(slidedTestSequenceCumLengths[i]-length(slidedTestSequence[[i]])+1):slidedTestSequenceCumLengths[i]] <- currentOwnPrediction$s
  }
  
  modelEmissionMatrix <- do.call(rbind,object$model$parms.emission$prob)
  
  #Once we have the state predictions for the whole test sequence, compute the score for each observation
  observationScores <- laply(which(!is.na(testSequence)), function (obsTimeStep) {
    currentObservedValue <- testSequence[obsTimeStep] #Value from 0 to V-1
    durationDeviation <- round(durationDeviation.hsmm(object,statePredictions[obsTimeStep]))
    
    #Get the observation probs vector in the range [t - sd, t + sd] that maximizes the probability of the currentObservedValu
    scopeObservationProbs <- laply((max(1,obsTimeStep-durationDeviation):min(length(testSequence),obsTimeStep+durationDeviation)),
                                   function (t) {
                                     return(modelEmissionMatrix[statePredictions[t],])  
                                   })
    #     print(currentObservedValue)
    #     print(durationDeviation)
    #     print(scopeObservationProbs)
    
    that <- which.max(apply(scopeObservationProbs,1, function (oProb) oProb[currentObservedValue+1]))
    # print(paste("that:",that))
    
    score <- exp(-1*C*((max(scopeObservationProbs[that,]) - scopeObservationProbs[that,currentObservedValue+1])))
    # print(paste("score",score))
    score
  })
  
  aux <- sapply(testSequence[!is.na(testSequence)], function (obs) {observationWeights[obs+1]})
  #print(matrix(c(testSequence[!is.na(testSequence)],observationScores,aux),ncol=3))
  #print(observationWeights)
  
  return(weighted.mean(observationScores,laply(testSequence[!is.na(testSequence)], function (obs) observationWeights[obs+1])))
}

#
# PREDICTION SCORE
#
validation.predictionScore.hsmm <- function (object, testSequence, predictionLength = 10, predictionDepth = object$J, C=1, observationWeights = rep(1/object$J,object$J), ...) {
  slidedTestSequence <- helper.slice(x = testSequence, n = predictionLength)
  slidedTestSequenceCumLengths <- cumsum(laply(slidedTestSequence,length))
  
  #
  # 1. Prediction Generation
  #
  currentMostLikelyPrediction <- NULL
  statePredictions <- matrix(nrow=length(testSequence),ncol=object$J)
  for (i in 1:length(slidedTestSequence)) {
    
    currentOwnPrediction <- .stateSequencePrediction.hsmm(object,
                                                          lastPrediction = currentMostLikelyPrediction,
                                                          length = length(slidedTestSequence[[i]]),
                                                          depth = predictionDepth)
    
    #Update the current most likely prediction via the viterbi algorithm, and get only the information 
    # about the current slide
#     currentMostLikelyPrediction <- predict.hsmm(object = object,
#                                                 newdata = mhsmmHelper.formatMultinomialData(testSequence[1:slidedTestSequenceCumLengths[i]],
#                                                                                             slidedTestSequenceCumLengths[i]),
#                                                 method = "smoothed")
    currentMostLikelyPrediction <- predict.hsmm(object = object,
                                            newdata = mhsmmHelper.formatMultinomialData(testSequence[1:slidedTestSequenceCumLengths[i]],
                                                                                        slidedTestSequenceCumLengths[i]),
                                            method = "viterbi")
    
    currentMostLikelyPrediction$s <- currentMostLikelyPrediction$s[(slidedTestSequenceCumLengths[i]-length(slidedTestSequence[[i]])+1):slidedTestSequenceCumLengths[i]]
#     currentMostLikelyPrediction$s <- tail(currentMostLikelyPrediction$s,n=1)
    currentMostLikelyPrediction$p <- object$model$transition[tail(currentMostLikelyPrediction$s,n=1),]
    #currentMostLikelyPrediction$p <- currentMostLikelyPrediction$p[(slidedTestSequenceCumLengths[i]-length(slidedTestSequence[[i]])+1):slidedTestSequenceCumLengths[i],,drop=FALSE]
#     currentMostLikelyPrediction$p <- t(apply(currentMostLikelyPrediction$p,
#                                              1,
#                                              function (row) .statesUncertainty.depthPrun(row,depth = predictionDepth)))
    
    #print(currentMostLikelyPrediction$s)
    #print(currentMostLikelyPrediction$p)
    #Save own predictions
    statePredictions[(slidedTestSequenceCumLengths[i]-length(slidedTestSequence[[i]])+1):slidedTestSequenceCumLengths[i],] <- currentOwnPrediction$p
  }
  
  modelEmissionMatrix <- do.call(rbind,object$model$parms.emission$prob)
  
  #Once we have the state predictions for the whole test sequence, compute the score for each observation
  observationScores <- laply(which(!is.na(testSequence)), function (obsTimeStep) {
    currentObservedValue <- testSequence[obsTimeStep] #Value from 0 to V-1
    durationDeviation <- .durationDeviation.hsmm(object,statePredictions[obsTimeStep,])
    #durationDeviation <- 1 #DEBUG
    #Get the observation probs vector in the range [t - sd, t + sd] that maximizes the probability of the currentObservedValu
    scopeObservationProbs <- laply((max(1,obsTimeStep-durationDeviation):min(length(testSequence),obsTimeStep+durationDeviation)),
                                 function (that) {
                                  return(t(modelEmissionMatrix) %*% statePredictions[obsTimeStep,])  
                                 })
#     print(currentObservedValue)
#     print(durationDeviation)
#     print(scopeObservationProbs)

    that <- which.max(apply(scopeObservationProbs,1, function (oProb) oProb[currentObservedValue+1]))
    # print(paste("that:",that))
    
    score <- exp(-1*C*((max(scopeObservationProbs[that,]) - scopeObservationProbs[that,currentObservedValue+1])))
    # print(paste("score",score))
    score
  })
  
  aux <- sapply(testSequence[!is.na(testSequence)], function (obs) {observationWeights[obs+1]})
  #print(matrix(c(testSequence[!is.na(testSequence)],observationScores,aux),ncol=3))
  #print(observationWeights)
  return(weighted.mean(observationScores,laply(testSequence[!is.na(testSequence)], function (obs) observationWeights[obs+1])))
}

#
# Validates the predictive capabilities of a model by a simple and direct comparison between the predicted states
# of a model and the ones returned by the viterbi algorithm. Test sequence may contain missing data in some time steps (NAs)
#
validation.RMSE.hsmm <- function (object, testSequence, predictionLength = 10, predictionDepth = 1,...) {

  slidedTestSequence <- helper.slice(x = testSequence, n = predictionLength)
  slidedTestSequenceCumLengths <- cumsum(laply(slidedTestSequence,length))
  
  #The sliding strategy avoid the predicition metric to be determinned by the beginning of the model
  #The idea is that we emit sets of [[predictionLength]]-length predictions with the current state of the model given a part (slide) 
  # of the test sequence.
  currentMostLikelyPrediction <- NULL
  errorAcc <- 0
  for (i in 1:length(slidedTestSequence)) {
    
    currentOwnPrediction <- .stateSequencePrediction.hsmm(object,
                                                       lastPrediction = currentMostLikelyPrediction,
                                                       length = length(slidedTestSequence[[i]]),
                                                       depth = predictionDepth)
    
    #Update the current most likely prediction via the viterbi algorithm, and get only the information 
    # about the current slide
    currentMostLikelyPrediction <- predict.hsmm(object = object,
                                                newdata = mhsmmHelper.formatMultinomialData(testSequence[1:slidedTestSequenceCumLengths[i]],
                                                                                            slidedTestSequenceCumLengths[i]),
                                                method = "smoothed")
    
    currentMostLikelyPrediction$s <- currentMostLikelyPrediction$s[(slidedTestSequenceCumLengths[i]-length(slidedTestSequence[[i]])+1):slidedTestSequenceCumLengths[i]]
    currentMostLikelyPrediction$p <- currentMostLikelyPrediction$p[(slidedTestSequenceCumLengths[i]-length(slidedTestSequence[[i]])+1):slidedTestSequenceCumLengths[i],,drop=FALSE]
    currentMostLikelyPrediction$p <- t(apply(currentMostLikelyPrediction$p,
                                             1,
                                             function (row) .statesUncertainty.depthPrun(row,depth = predictionDepth)))
    #print(currentOwnPrediction$p)
    
    # Compare our own state prediction with the most likely state prediction, for the current slided sequence
    NonNAEventIndexes <- which(!is.na(slidedTestSequence[[i]]))
    
    currentSlideErrors <- laply(NonNAEventIndexes, function (eventIndex) {
      b <- sapply(object$model$parms.emission$prob,"[[",(slidedTestSequence[[i]][eventIndex])+1)
      pPredicted <- sum(b*currentOwnPrediction$p[eventIndex,])
      pMostLikely <- sum(b*currentMostLikelyPrediction$p[eventIndex,])
      
      #(pMostLikely - pPredicted)^2
      abs(pMostLikely - pPredicted)
    })
    
    errorAcc <- errorAcc + sum(currentSlideErrors)
  }
  
  #Return the Root Mean Square Error over all the non NA events in the test sequence
  return((errorAcc/length(na.exclude(testSequence))))
}

##################################################
#
# STATE SEQUENCE PREDICTION (SIMPLE)
#
##################################################
.stateSequencePrediction.simple.hsmm <- function (object, lastPrediction, length) {
  # Create a first model-based prediction in the case that no last prediction is passed
  if (is.null(lastPrediction)) {
    lastPrediction$s <- which.max(object$model$init)
  }
  
  ret <- list()
  ret$s <- numeric(length)
  
  currentStates <- list()
  currentStates$s <- lastPrediction$s[length(lastPrediction$s)]
  
  predictionAssignmentsCount <- 0
  firstIterationFlag <- TRUE
  while (predictionAssignmentsCount < length) {
    currentExpectedDuration <- round(expectedDuration.hsmm(object,j=currentStates$s))
    
    #Restar tiempo actual de prediccion en la primera iteracion
    if (firstIterationFlag == TRUE) {
      currentExpectedDuration = currentExpectedDuration - tail(rle(lastPrediction$s)$lengths,1)
      firstIterationFlag <- FALSE
    }
    
    currentPredictionLength <- min(length-predictionAssignmentsCount,currentExpectedDuration)
    
    if (currentPredictionLength > 0) {
      ret$s[(predictionAssignmentsCount+1):(predictionAssignmentsCount+currentPredictionLength)] <- rep(currentStates$s,currentPredictionLength)
      
      #Update assignments count
      predictionAssignmentsCount = predictionAssignmentsCount + currentPredictionLength
    }
    
    #Update current state information based on the model transition matrix
    currentStates$s <- which.max(object$model$transition[currentStates$s,])
  }
  
  ret
  
}

##################################################
#
# STATE SEQUENCE PREDICTION (PROBABILISTIC)
#
##################################################
.stateSequencePrediction.hsmm <- function (object,lastPrediction,length,depth=object$J) {
  
  # Create a first model-based prediction in the case that no last prediction is passed
  if (is.null(lastPrediction)) {
    lastPrediction <- list()
#       lastPrediction$p <- matrix(.statesUncertainty.depthPrun(object$model$init,depth),
#                                  nrow = 1)
      lastPrediction$p <- .statesUncertainty.depthPrun(object$model$init,depth)
      lastPrediction$s <- c(which.max(lastPrediction$p))
  }
  
  ret <- list()
  ret$s <- numeric(length)
  ret$p <- matrix(nrow=length,ncol=length(lastPrediction$p))
  
  currentStates <- list()
    currentStates$s <- lastPrediction$s[length(lastPrediction$s)]
    #currentStates$p <- lastPrediction$p[nrow(lastPrediction$p),]
    #urrentStates$s <- lastPrediction$s
    currentStates$p <- lastPrediction$p
    
  predictionAssignmentsCount <- 0
  firstIterationFlag <- TRUE
  while (predictionAssignmentsCount < length) {
    #currentExpectedDuration <- .expectedDuration.hsmm(object,currentStates$p,depth=ncol(lastPrediction$p))
    currentExpectedDuration <- .expectedDuration.hsmm(object,currentStates$p)
    
    #Restar tiempo actual de prediccion en la primera iteracion TODO
    if (firstIterationFlag == TRUE) {
      currentExpectedDuration = currentExpectedDuration - tail(rle(lastPrediction$s)$lengths,1)
      firstIterationFlag <- FALSE
    }
    
    currentPredictionLength <- min(length-predictionAssignmentsCount,currentExpectedDuration)
    
    if (currentPredictionLength > 0) {
      ret$s[(predictionAssignmentsCount+1):(predictionAssignmentsCount+currentPredictionLength)] <- rep(currentStates$s,currentPredictionLength)
      ret$p[(predictionAssignmentsCount+1):(predictionAssignmentsCount+currentPredictionLength),] <- rep(currentStates$p,each=currentPredictionLength)
      
      #Update assignments count
      predictionAssignmentsCount = predictionAssignmentsCount + currentPredictionLength
    }
    
    #Update current state information based on the model transition matrix
    currentStates <- .nextStates(object,currentStates,depth = depth)
  }
  
  ret
}

#
# Weighted Expected Duration of a state, given the last states uncertainty vector
#
# .expectedDuration.hsmm <- function (object,statesUncertainty,depth=1) {
#   
#   if (depth > length(statesUncertainty))
#     stop("Prediction depth cannot be higher than the number of states in the model")
#   
#   # The expected duration of the current prediction is given by a weighted mean among the mean duration of each of the
#   # last states in the last predicion. The weight for each state is given by its prior prediction probability
#   sortedLastStates <- sort(statesUncertainty, index.return = TRUE)
#   
#   ret <- 0
#   for (i in 1:depth) {
#     ret <- ret + (sortedLastStates$x[i]/sum(sortedLastStates$x[1:depth]))*expectedDuration.hsmm(object,sortedLastStates$ix[i])
#   }
#   
#   round(ret)
# }

# The expected duration of the current prediction is given by a weighted mean among the mean duration of each of the
# last states in the last predicion. The weight for each state is given by its prior prediction probability
.expectedDuration.hsmm <- function (object,statesUncertainty) {
  round(
    weighted.mean(
      sapply(1:length(statesUncertainty),function (j) expectedDuration.hsmm(object,j)),
      statesUncertainty
      )
    )
}

#Weighted mean of the standard deviations of the duration distributions
.durationDeviation.hsmm <- function (object, statesUncertainty) {
  round(
    weighted.mean(
      sapply(1:length(statesUncertainty),function (j) durationDeviation.hsmm(object,j)),
      statesUncertainty
    )
  )
}


#
# Get the next possible states, along with the uncertainty of being in each of them, 
# given the current states uncertainty
#
.nextStates <- function (object,currentStates, depth) {
  ret <- list()
    ret$p <- c(t(object$model$transition) %*% currentStates$p)
    ret$p <- .statesUncertainty.depthPrun(ret$p,depth)
    ret$s <- which.max(ret$p)
    
  ret
}

#
# Prun a states uncertainty vector to get only the [depth] most probable states, normalizing
# their probabilities. All states out of the prun get p(s) = 0
#
.statesUncertainty.depthPrun <-function(statesUncertainty,depth) {
  ranking <- sort(statesUncertainty,decreasing = TRUE,index.return=TRUE)
  statesUncertainty[-ranking$ix[1:depth]] <- 0
  statesUncertainty <- statesUncertainty/sum(statesUncertainty)
  
  statesUncertainty
}

#
# Expected scorefor the predictions of a model, given some observation weights (STATE-ORIENTED)
#
validation.score.simple.expected.hsmm <- function (object,obsWeights, C = 1) {
  weighted.mean(
    laply(1:object$J, function (j) {
      weighted.mean(max(object$model$parms.emission$prob[[j]])-object$model$parms.emission$prob[[j]],obsWeights)
    }),
    laply(1:object$J, function (j) expectedDuration.hsmm(object,j))
  )
}

#
# Expected scorefor the predictions of a model, given some observation weights (OBSERVATION-ORIENTED)
#
validation.score.simple.expected.hsmm.2 <- function (object,obsWeights,C=1) {
  modelEmissionMatrix <- do.call(rbind,object$model$parms.emission$prob)
  weighted.mean(
    laply(1:length(obsWeights), function (o) {
      oProbs <- modelEmissionMatrix[,o]
      exp(-1*(max(oProbs) - mean(oProbs)))
    }),
    obsWeights
  )
}

#
# Tercer intento para acabar el EVO
#
validation.score.simple.expected.hsmm.3 <- function (object, obsWeights, C=1) {
  modelEmissionMatrix <- do.call(rbind,object$model$parms.emission$prob)
  weighted.mean(
    laply(1:length(obsWeights), function (o) {
      oMaxAvg <- mean(apply(modelEmissionMatrix,1,max))
      oProbs <- modelEmissionMatrix[,o]
      exp(-1*(oMaxAvg - mean(oProbs)))
    }),
    obsWeights
  )
}

#
# Expected scorefor the predictions of a model, given some observation weights
#
validation.score.simple.sd.hsmm <- function (object,obsWeights, C = 1) {
  wt.sd(
    laply(1:object$J, function (j) {
      weighted.mean(max(object$model$parms.emission$prob[[j]])-object$model$parms.emission$prob[[j]],obsWeights)
    }),
    laply(1:object$J, function (j) expectedDuration.hsmm(object,j))
  )
}

validation.score.simple.sd.hsmm.2 <- function (object,obsWeights,C = 1) {
  modelEmissionMatrix <- do.call(rbind,object$model$parms.emission$prob)
  wt.sd(
    laply(1:length(obsWeights), function (o) {
      oProbs <- modelEmissionMatrix[,o]
      exp(-1*(max(oProbs) - mean(oProbs)))
    }),
    obsWeights
  )
}

validation.score.simple.sd.hsmm.3 <- function (object,obsWeights,C = 1) {
  modelEmissionMatrix <- do.call(rbind,object$model$parms.emission$prob)
  wt.sd(
    laply(1:length(obsWeights), function (o) {
      oMaxAvg <- mean(apply(modelEmissionMatrix,1,max))
      oProbs <- modelEmissionMatrix[,o]
      exp(-1*(oMaxAvg - mean(oProbs)))
    }),
    obsWeights
  )
}

validation.score.simple.upperBound <- function (object,obsWeights,C=1) {
  modelEmissionMatrix <- do.call(rbind,object$model$parms.emission$prob)
  
  weighted.mean(
    laply(1:length(obsWeights), function (o) {
      #oMaxAvgSd <- mean(apply(modelEmissionMatrix,1,max)) + sd(apply(modelEmissionMatrix,1,max))
      oMaxAvgSd <- mean(apply(modelEmissionMatrix,1,max))
      oProbsAvgSd <- mean(modelEmissionMatrix[,o]) + sd(modelEmissionMatrix[,o])
      exp(-1*(oMaxAvgSd - oProbsAvgSd))
    }),
    obsWeights
  )
}

validation.score.simple.lowerBound <- function (object,obsWeights,C=1) {
  modelEmissionMatrix <- do.call(rbind,object$model$parms.emission$prob)
  weighted.mean(
    laply(1:length(obsWeights), function (o) {
#       oMaxAvgSd <- mean(apply(modelEmissionMatrix,1,max)) - sd(apply(modelEmissionMatrix,1,max))
      oMaxAvgSd <- mean(apply(modelEmissionMatrix,1,max))
      oProbsAvgSd <- mean(modelEmissionMatrix[,o]) - sd(modelEmissionMatrix[,o])
      exp(-1*(oMaxAvgSd - oProbsAvgSd))
    }),
    obsWeights
  )
}

caca <- function (obsProbs) {
  generalObsWeights <- (1/generalObsProbs)/sum(1/generalObsProbs)
  
}
