#
# Apply a funtion over every sequence of a "march.Dataset" object
#
# Arguments   y : An object of class "march.Dataset"
#             f : A function to apply. The function has an argument "s"
#
march.helper.apply_dataset <- function (y,f) {
  sapply(X = 1:y@N, FUN = function (i) f(march.dataset.h.extractSequence(y,i)))
}

#
# Compute the metric "Number of Rare States" (NRS) for a given object of class "march.dcmm"
#
# Input   d : An object of class "march.dcmm"
#         State_Visit_Frequency_Threshold : Sets the minimum frequency to consider a state as "rare"
#
march.helper.NRS <- function (d, State_Visit_Frequency_Threshold = 0.05) {
  mpp <- march.dcmm.viterbi(d = d, y = d@y)
  
  #Prepare STS object
  mpp_sts <- lapply(mpp,function (x) paste(x,collapse="-"))
  mpp_sts <- TraMineR::seqdef(as.data.frame(mpp_sts) %>% t)
  freqs <- apply(X = TraMineR::seqstatd(seqdata = mpp_sts)$Frequencies,
                 MARGIN = 1,
                 FUN = mean)
  return(freqs[freqs < State_Visit_Frequency_Threshold] %>% length)
}


#
# Predict the next observation of sequence "s" using the model "d"
#
# Input   d : An object of class "march.dcmm"
#         s : A test sequence (object of class "march.sequence")
#
#         This prediction method has been designed inspired by the work of Huang, 2009.
#         "Developing an Abstraction Layer for the Visualization of HSMM-Based Predictive Decision Support" (page 64)
#
#         TODO: Parallel computing here!?
#
march.helper.predictNextObservation <- function (d,s,return_probs=F,normalize_probs = T) {
  predictedLogProbs <- sapply(X = 1:length(d@y@dictionary),FUN = function (obs) {
    newSeq <- new(Class = "march.Sequence",
                  y = c(s@y,obs),
                  weight = s@weight,
                  N = as.integer(s@N+1))
    march.helper.sequenceLogLikelihood(d = d, s = newSeq)
  })
  names(predictedLogProbs) <- d@y@dictionary
  
  #Tomorrow. Normalize probs
  if (return_probs == F) return(which.max(predictedLogProbs))
  else {
    if (normalize_probs == F) return(predictedLogProbs)
    else {
      convertedProbs <- exp(predictedLogProbs + floor(abs(max(predictedLogProbs))))
      return(convertedProbs/sum(convertedProbs))
    }
  }
}

#
# Compute the Prediction Score of a given model for a given test sequence at time t.
# The prediction generation method is coded in the function "march.helper.predictNextObservation"
#
# Arguments   d : An object of class "march.dcmm"
#         s : A test sequence(s) (object of class "march.sequence")
#         timeSteps : Time steps in which the prediction score will be computed (numeric vector)
#
# Value  numeric vector (same length as t) of prediction scores, ranging between 0 (worst) and 1 (best)
#
# Details   The score calculation has been inspired by the works of Rodríguez-Fernández et.al in 
#         "A method for building predictive HSMMs in interactive environments"
#
# TODO    Parallel!?
march.helper.predictionScore <- function (d,s,timeSteps,C = 1) {
  sapply(X = timeSteps, FUN = function (t) {
    # [1-t-1] subsequence
    s_sub <- new(Class = "march.Sequence",
                  y = s@y[1:(t-1)],
                  weight = s@weight,
                  N = as.integer(t-1))
    
    actualObs <- s@y[t]
    predictionProbs <- march.helper.predictNextObservation(d,s_sub,return_probs=T,normalize_probs = T)
    
    #Exponential function!
    exp(-1*C*(max(predictionProbs) - predictionProbs[[actualObs]]))
  })
}

#
# Compute the metric "Model Accuracy Score" (MAS) for a given object of class "march.dcmm"
# Online Accuracy Score
#
# Arguments   d : An object of class "march.dcmm"
#             s : A test sequence(s) (object of class "march.sequence")
#             t : numeric. Time step (WARNING: t must be greater or equal than n)
#             n : Aggregation level, i.e, width of the sliding window
#
# Notes   The general equation is taken from the work of Cummings et al., but the 
#         quality of the predictions is computed following the works of WCCI
#
march.helper.MAS <- function (d,s,t, n = 10, penalty_rate = 2) {
  stopifnot(t >= n)
  
  # Compute the qualities of all the prediciotns from t-n to t
  qualities <- march.helper.predictionScore(
    d = d,
    s = s,
    timeSteps = (t-n+1):t,
    C = penalty_rate
  )
  
  return(mean(qualities))
}


# 
#
# Input   d : An object of class "march.dcmm"
#         s : A test sequence(s) (object of class "march.sequence")
#
#
march.helper.MMAS <- function (d,s, penalty_rate = 2) {
  march.helper.MAS(
    d = d, 
    s = s, 
    t = s@N, 
    n = s@N - d@orderHC - d@orderVC,
    penalty_rate = penalty_rate
  )
}

##
# Compute the confusion matrix obtained after predicting each of the symbols of a test sequence
# using the given model. The prediction generation method is coded in the function 
# "march.helper.predictNextObservation"
#
# Arguments  d : An object of class "march.dcmm"
#            s : A test sequence(s) (object of class "march.sequence")
# 
##
march.helper.confusion_matrix <- function(d,s) {
  confusionMatrix <- matrix(0,nrow=length(d@y@dictionary),ncol=length(d@y@dictionary))
  rownames(confusionMatrix) <- d@y@dictionary
  colnames(confusionMatrix) <- d@y@dictionary
  
  for (t in (d@orderHC + d@orderVC):s@N) {
    s_sub <- new(Class = "march.Sequence",
                 y = s@y[1:(t-1)],
                 weight = s@weight,
                 N = as.integer(t-1))
    
    actualObs <- s@y[t]
    predictionObs <- march.helper.predictNextObservation(d,s_sub,return_probs=F)
    confusionMatrix[predictionObs,actualObs] <- confusionMatrix[predictionObs,actualObs] + 1
  }
  
  confusionMatrix
}

##
# Compute the measure "Min Precision", which applies a multiclass precision calculation
# over the confusion matrix obtained by predicting sequence "s" with model "d"
#
##
march.helper.minPrecision <- function (d,s) {
  # Primero creamos la matriz de confusion multiclase
  confusionMatrix <- march.helper.confusion_matrix(d,s)
  
  precision_multiclass <- sapply(
    X = 1:nrow(confusionMatrix),
    FUN = function (i) {
      ifelse(
        test = all(confusionMatrix[i,] == 0),
        yes = NA,
        no = confusionMatrix[i,i]/sum(confusionMatrix[i,])
      )
    }
  )
  print(precision_multiclass)
  
  return(min(precision_multiclass,na.rm = TRUE))
}

#Logistic function. The order of the chain indicates the steepness of the function
.transition_relevancy_score <- function (prob,order,transition_relevancy_threshold) {
  return(
    1 / (1 + exp(-1*order*(prob - transition_relevancy_threshold)))
  )
}

##
# Coefficient of High-Probability hidden transitions (CHPHT)
# Calcula el nivel de relevancia de las transiciones de un DCMM para la cadena oculta
#
# Input:  d : An object of class "march.dcmm"
#
##
march.helper.CHPHT <- function (d, transition_relevancy_threshold = 0.75) {
  
  #hidden chain
  averageHiddenRTC <- apply(
    X = d@A,
    MARGIN = 1,
    FUN = max
  ) %>% 
    map_dbl(.f = ~ .transition_relevancy_score(
      prob = ., 
      order = d@orderHC, 
      transition_relevancy_threshold = transition_relevancy_threshold)
      ) %>% 
    mean
  
  averageHiddenRTC
}


##
# Coefficient of High-Probability Visible Transitions (CHPVT)
# Calcula el nivel de relevancia de las transiciones de un DCMM para la cadena oculta
#
# Input:  d : An object of class "march.dcmm"
#
##
march.helper.CHPVT <- function (d, transition_relevancy_threshold = 0.75) {
  #media entre los estados
  RVTC_by_state <- sapply(1:d@M, function (m) {
    if (class(d@RB[m,,]) != "matrix")
      visibleTransitions <- matrix(d@RB[m,,],nrow=1)
    else
      visibleTransitions <- d@RB[m,,]
    
    averageVisibleRTC <- apply(
      X = visibleTransitions,
      MARGIN = 1,
      FUN = function (row) {
        ifelse(
          test = sum(row) == 0,
          yes = NA,
          no = max(row)
        )
      }
    ) %>% 
      map_dbl(.f = ~ .transition_relevancy_score(
        prob = ., 
        order = d@orderVC,
        transition_relevancy_threshold = transition_relevancy_threshold
        )) %>% 
      mean(na.rm = T)
    
    averageVisibleRTC
  })
  
  return(min(RVTC_by_state))
}

#
# Compute the probability that a given model "d" can produce the sequence "s"
#
# Input   d : An object of class "march.dcmm"
#         s : A test sequence(s) (object of class "march.sequence")
#
#
march.helper.sequenceLogLikelihood <- function (d, s) {
  stopifnot(s@N >= d@orderHC + d@orderVC) #Error control
    (march:::march.dcmm.forward(d = d, s = s))$LL
}

#
# Compute the mean probability that a given model "d" can produce a sequence of the dataset "y"
#
#
# Input   d : An object of class "march.dcmm"
#         y : An object of class "march.dataset"
#
#
march.helper.MSLL <- function (d,y) {
  map_dbl(.x = 1:y@N, 
          .f = ~ march.helper.sequenceLogLikelihood(d,march.dataset.h.extractSequence(y,.))) %>% 
    mean
}