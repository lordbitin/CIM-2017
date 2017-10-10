#Input : HMM or HSMM
#Output: Number of iterations neccesary to fit this object
mhsmmHelper.iterations <- function(object) {
  if (class(object)!= "hmm" && class(object) !="hsmm")
    stop("Wrong object passed as argument. Need an hmm or an hsmm")
  
  return(length(object$loglik))
}

# INPUT: HSMM data object
# OUTPUT: List of indexes of each observation sequence
mhsmmHelper.getDataIndexes <- function(hsmmData) {
  if (class(hsmmData)!="hsmm.data") 
    stop("hsmmData must have class hsmm.data")
  
  split(1:sum(hsmmData$N),rep(1:length(hsmmData$N),hsmmData$N))
}

#INPUT: hsmm.data object from the package mhsmm
#OUTPUT: hsmm.data object, random resample of the input data
mhsmmHelper.resampleHSMMData <- function(hsmmData) {
  if (class(hsmmData)!="hsmm.data") 
    stop("hsmmData must have class hsmm.data")
  
  randIndexes <- sample(mhsmmHelper.getDataIndexes(hsmmData))
  
  resample <- list()
  resample$x <- unlist(lapply(randIndexes, function (sequenceIndexes) {
      if (class(hsmmData$x) == "matrix")
        hsmmData$x[sequenceIndexes,]
      else
        hsmmData$x[sequenceIndexes]
    })
    )
  if (!is.null(hsmmData$timing))
    resample$timing <- unlist(lapply(randIndexes, function (sequenceIndexes) hsmmData$timing[sequenceIndexes]))
  
  resample$N <- sapply(as.numeric(names(randIndexes)), function (nIndex) hsmmData$N[[nIndex]])
  
  class(resample) <- "hsmm.data"
  names(resample$x) <- NULL
  names(resample$timing) <- NULL
  
  resample
}


mhsmmHelper.getMultinomialVector <- function(option,minOption,maxOption) {
  v <- rep(0,(maxOption-minOption+1))
  v[option-minOption+1] <- 1
  v
}

#Input : Data with the option number of each trial (0,4,3,2,5,6,6,6,8...)
#Output : hmmdata class to fit hmms and hsmms with the package mhsmm
mhsmmHelper.formatMatrixMultinomialData <- function (data,numberOfSequences,minOption, maxOption) {
  
  data_ <- t(sapply(data,function (trial) mhsmmHelper.getMultinomialVector(trial,minOption,maxOption)))
  train <- list(x=data_,N = numberOfSequences)
  class(train) <- "hsmm.data"
  
  train
}

#Input : Data with the option number of each trial (0,4,3,2,5,6,6,6,8...)
#Output : hmmdata class to fit hmms and hsmms with the package mhsmm
mhsmmHelper.formatMultinomialData <- function (data,numberOfSequences) {
  train <- list(x=data,N = numberOfSequences)
  class(train) <- "hsmm.data"
  train
}

# train : hsmm.data object
#index : index of the observation sequence to retrieve
mhsmmHelper.getObservationSequenceByIndex <- function (train,index,getTiming = FALSE) {
  
  ret <- list()
  dataIndexes <- split(1:sum(train$N),rep(1:length(train$N),train$N))
  
  ret$x <- train$x[dataIndexes[[index]]]
  
  if (getTiming == TRUE) {
    ret$timing <- train$timing[dataIndexes[[index]]]
  }
  
  ret
}

#Input = multinommial data (hsmm.data)
mhsmmHelper.initEmissionParms.matrixMultinomial.segment <- function (x,nStates, randomSegmentation = TRUE) {
  
  if (class(x)!="hsmm.data")
    stop("x is not a hsmm data object")

  probs <- list()
  
  #Split the data into J chunks, where J is the number of states
  data <- x$x
  nMultiOptions <- ncol(data)
  s<-split(data, (seq(nrow(data))-1) %/% (nrow(data)/nStates))
  
  for (i in 1:nStates) {
    segment <- matrix(s[[i]],ncol=nMultiOptions)
    probs[[i]] <- apply(segment,2,function (dataFromOption) sum(dataFromOption,na.rm = TRUE)/sum(!is.na(dataFromOption)))
  }
  
  probs
}

#INPUT - Interger 1-size Multinomial Vector (Categorical variable)
#OUTPUT - Estimation of the e
mhsmmHelper.initEmissionParms.multinomial.segment <- function (x,nStates, nObsSymbols,
                                                               randomSegmentation = TRUE) {

  if (randomSegmentation) {
    rand <- sample(length((x)))
    x <- x[rand]
  }
  
  probs <- list()

  #Split the data into J chunks, where J is the number of states
  s<-split(x, (seq(length(x))-1) %/% (length(x)/nStates))
  
  for (i in 1:nStates) {
    probs[[i]] <- sapply(0:(nObsSymbols-1),function (obsSym) sum(s[[i]]==obsSym, na.rm = TRUE)/sum(!is.na(s[[i]])))
  }
  
  probs
}

#INPUT
#   - "hsmm.data"
#OUTPUT
#   - list("hsmm.data","hsmm.data) with the training and test datasets
mhsmmHelper.partitionHSMMData <- function(originalData, testSize, random = FALSE) {
  if (class(originalData)!= "hsmm.data")
    stop("class(data) must be hsmm.data (See mhsmm package)")
  
  #Randomize the data sequences (if needed)
  data <- originalData
  if (random) {
    data <- mhsmmHelper.resampleHSMMData(originalData)
  }
  
  #TEST : Simple partition (The first sequences are test sequences)
  test <- list()
  
    #x
    if (class(data$x) == "matrix")
      test$x <- data$x[1:sum(data$N[1:testSize]),]
    else
      test$x <- data$x[1:sum(data$N[1:testSize])]
    
    #timing
    if (!is.null(data$timing))
      test$timing <- data$timing[1:sum(data$N[1:testSize])]
    
    #N
    test$N <- data$N[1:testSize]
    class(test) <- "hsmm.data"

  #TRAINING : Exclude the test sequences from the training dataset
  training <- list()
  
    #x
    if (class(data$x) == "matrix")
      training$x <- data$x[-1:-sum(data$N[1:testSize]),]
    else
      training$x <- data$x[-1:-sum(data$N[1:testSize])]
    
    #timing
    if (!is.null(data$timing))
      training$timing <- data$timing[-1:-sum(data$N[1:testSize])]
    
    #N
    training$N <- data$N[-1:-testSize]
    class(training) <- "hsmm.data"
    
  list(training=training,test=test)
}

##
# Input
#   - data: class "hsmm.data" (mhsmm)
#   - J: Number of states
#   - dictionarySize: If NULL, it is a matrix multinomial (Only neccesary for vector multinomial)
#   - times: Numebr of times the HMM will be fitted using different init values
#   - mstep: Maximization step for the Baum-Welch Algoritm
##
mhsmmHelper.trainHMM <- function (data,J,dictionarySize=NULL,crossValidIter=10,times=10,testSequences=1,...) {
  
  if (class(data)!="hsmm.data")
    stop("data is not a hsmm.data object")
  if (times <=0)
    stop("times must be greater than 0")
  
  #Cross validation (k-fold cross validation) - TODO: Efficiency?
  MASs <- list()
  for (i in 1:crossValidIter) {
    print("################################")
    print(paste("BEGINNING CROSS VALIDATION (",i," of ",crossValidIter,")"))
    
    partitionedDataset <- mhsmmHelper.partitionHSMMData(originalData = data,
                                                        testSize = testSequences,
                                                        random=TRUE)
    
    bestHMM <- NULL
    newHMM <- NULL
    set.seed(1)
    for (j in 1:times) {
      
      #Random Initial probabilities vector (Sum up 1)
      init0 <- hmmHelper.randomInitialVector(J)
      
      #Random transition matrix
      A0 <- hmmHelper.randomTransitionMatrix.hmm(J)
      
      #Random emission distribution params (Segment the data ordered randmoly)
      if (class(data$x) == "matrix") {
        B0 <- list(probs=mhsmmHelper.initEmissionParms.matrixMultinomial.segment(partitionedDataset$training$x,
                                                                                 J,
                                                                                 randomSegmentation = TRUE)
        )    
      }
      else {
        B0 <- list(probs=mhsmmHelper.initEmissionParms.multinomial.segment(partitionedDataset$training$x,
                                                                           J,
                                                                           dictionarySize,
                                                                           randomSegmentation = TRUE)
        )       
      }
      
      #Call mhsmm fit
      startModel <- hmmspec(init = init0, trans = A0, parms.emission = B0, dens.emission = dmultinom.hsmm)
      print(paste("---Fitting HMM with",J,"states, validation",i,"of",crossValidIter,"iteration",j,"of",times,"---"))
      newHMM <- hmmfit(x =partitionedDataset$training , start.val = startModel, mstep = mstep.multinom,...)
      print(paste("---[FINISHED] logLik =",logLik(newHMM),"---"))
      
      if (is.null(bestHMM)) {
        bestHMM <- newHMM
      } else {
        if (logLik(newHMM) > logLik(bestHMM))
          bestHMM <- newHMM
      }
    }
    
    print(paste("ENDING CROSS VALIDATION (",i,"of",crossValidIter,")"))
    print("################################")
    #MAS calculation
    MASs[[i]] <- mhsmmHelper.MAS(model = bestHMM, data = partitionedDataset$test)
  }
  
  #MAS Average
  bestHMM$MAS <- mean(unlist(MASs))
  
  #Return the best model trained (Last?)
  return(bestHMM)
}

#
# INPUT: List of hmms or hsmms
# OUTPUT: Dataframe with a summary of different scores (loglik, BIC, MAS)
mhsmmHelper.createScoreSummary <- function (modelList) {
  
  #We suppose that all models are class-equally, so we only check the first one
  if (class(modelList[[1]])!="hmm" && class(modelList[[1]])!="hsmm")
    stop("ModelList must contain models of class hmm or hsmm")
  
  if (class(modelList[[1]]) == "hmm") {
    #MHSMM Models
    nstateList <- laply(modelList,function (hmm) hmm$K)
  } else if (class(modelList[[1]]) == "hsmm") {
    nstateList <- laply(modelList,function (hsmm) hsmm$J)
  }
  
  scoreSummary <- data.frame(nstates=nstateList,
                       logLik=laply(modelList,logLik),
                       BIC=laply(modelList,BIC),
                       class=laply(modelList,class),
                       MAS=laply(modelList, function (hmm) hmm$MAS))
  
  return(scoreSummary)
}

##
# Input
#   - data: class "hsmm.data" (mhsmm) (MULTINOMIAL)
#   - J: Number of states
#   - times: Number of times the HMM will be fitted using different init values
#   - M : Maximum state time units
#   - mstep: Maximization step for the Baum-Welch Algoritm
##
mhsmmHelper.trainHSMM <- function (data,
                                   J,
                                   dictionarySize,
                                   crossValidIter=10,
                                   times=10,
                                   M,
                                   testSequences=1,
                                   startModel = NULL,
                                   randomSegmentation = TRUE,
                                   initializeSojournDistribution = sojournDistribution.nonparametrical.uniform.hsmm,
                                   ...) {
  
  if (class(data)!="hsmm.data")
    stop("data is not a hsmm.data object")
  if (times <=0)
    stop("times must be grater than 0")
  
  #Cross validation (k-fold cross validation) - TODO: Efficiency?
  optimalModels <- list()
  MASs <- list()
  for (i in 1:crossValidIter) {
    print("################################")
    print(paste("BEGINNING CROSS VALIDATION (",i," of ",crossValidIter,")"))
    partitionedDataset <- mhsmmHelper.partitionHSMMData(originalData = data,
                                                        testSize = testSequences,
                                                        random = randomSegmentation)
    
    bestHSMM <- NULL
    newHSMM <- NULL
    
    #If a startmodel is given (hsmmspecs), the best hsmm possible only needs to get trained from that model.
    #If a start model is not given, we must try to search the best initial model
    if (!is.null(startModel)) {
      print(paste("---Fitting HSMM with",J,"states, validation",i,"of",crossValidIter,"---"))
      bestHSMM <- hsmmfit(x = partitionedDataset$training, 
                         model = startModel, 
                         mstep = mstep.multinom, 
                         M = M,
                         ...)
      print(paste("---[FINISHED] logLik =",logLik(bestHSMM),"---"))
      
    } else {
      set.seed(1)
      for (j in 1:times) {
        
        #Random Initial probabilities vector (Sum up 1, DIAGONAL 0)
        init0 <- hmmHelper.randomInitialVector(J)
        
        # Random transition matrix
        A0 <- hmmHelper.randomTransitionMatrix.hsmm(J)
        
        #TODO: Random sojourn distribution (Non-parametric)
        sojournD <- initializeSojournDistribution(M,J)
        
        #Random emission distribution params (Segment the data ordered randmoly)
        if (class(data$x) == "matrix") {
          B0 <- list(probs=mhsmmHelper.initEmissionParms.matrixMultinomial.segment(partitionedDataset$training$x,
                                                                                   J,
                                                                                   randomSegmentation = TRUE)
          )    
        }
        else {
          B0 <- list(probs=mhsmmHelper.initEmissionParms.multinomial.segment(partitionedDataset$training$x,
                                                                             J,
                                                                             dictionarySize,
                                                                             randomSegmentation = TRUE)
          )       
        }
        
        #Call mhsmm fit
        startModel_ <- hsmmspec(init = init0, 
                               transition = A0, 
                               parms.emission = B0, 
                               sojourn = sojournD,
                               dens.emission = dmultinom.hsmm)
        
        print(paste("---Fitting HSMM with",J,"states, validation",i,"of",crossValidIter,"iteration",j,"of",times,"---"))
        print(startModel_)
        newHSMM <- hsmmfit(x = partitionedDataset$training, 
                           model = startModel_, 
                           mstep = mstep.multinom, 
                           M = M,
                           ...)
        print(paste("---[FINISHED] logLik =",logLik(newHSMM),"---"))
        
        if (is.null(bestHSMM)) {
          bestHSMM <- newHSMM
          bestHSMM$startSpecs <- startModel_
        } else {
          if (logLik(newHSMM) > logLik(bestHSMM))
            bestHSMM <- newHSMM
            bestHSMM$startSpecs <- startModel_
        }
      }
    }
    
    print(paste("ENDING CROSS VALIDATION (",i,"of",crossValidIter,")"))
    print("################################")
    
    #MAS calculation
    #MASs[[i]] <- mhsmmHelper.MAS(model = bestHSMM,
     #                            data = partitionedDataset$test)
    optimalModels[[i]] <- bestHSMM
  }
  
  #Return the optimal model with information about the cross validation error and prediction power
  ret <- optimalModels[[which.max(laply(optimalModels,logLik))]]
  #ret$MAS <- mean(unlist(MASs))
  
  #Return the best model trained
  return(ret)
}

#
# INPUT : 
#   - train: "hsmm.data" object from the package mhsmm  with the original observation sequences
#   - timing: Vector with the timing data for each observation symbol in train$x
#   - TSR (Time Slot Resolution) (Milliseconds)
# 
mhsmmHelper.timeHandling <- function(train,timing,timeHandlingStrategy,TSR = 100) {
  
  if (class(train)!="hsmm.data")
    stop("train must be a hsmm.data object (See mhsmm package)")
  
  if (length(timing)!= sum(train$N))
    stop("timing length must be equal to sum(train$N)")
  
  if (missing(timeHandlingStrategy))
    stop("No timeHandling strategy specified")
  
  # Convert the training data into a list (for both vector and matrix input data)
  data <- as.list(data.frame(t(train$x)))
  
  dataIndexes <- split(1:sum(train$N),rep(1:length(train$N),train$N))
  newData <- llply(1:length(dataIndexes), function (i) {
    
    newObservationSeq <-timeHandlingStrategy(data[dataIndexes[[i]]],
                                             timing[dataIndexes[[i]]],
                                             TSR) 
    newObservationSeq
  })
  newN <- sapply(newData,length)
  
  #Create the new hsmm.data object (distinguish between vectors and matrix data types)
  if (class(train$x)=="matrix")
    newData <- matrix(unlist(m),ncol=ncol(train$x),byrow=TRUE)    # Multi-variate
  else 
    newData <- unlist(newData) #Single variate
  
  newTrain <- list(x=newData,N = newN)
  class(newTrain) <- "hsmm.data"
  
  newTrain
}

#
# Input:
# data: hsmm.data object
#
mhsmmHelper.getObservationFrequencies <- function (data,grammarSize) {
  sapply(0:(grammarSize-1), function (obsSymbol) {
    length(data$x[data$x == obsSymbol])/length(!is.na(data$x))
  })
}


#Internal function to calculate quality subscores
# Symbol ranking is a list of symbols
mhsmmHelper.getQualitySubScore <- function(symbol,symbolRanking) {
  ratings <- c(100,70,30)

  ranking <- Position(function(x) identical(x, as.numeric(symbol)), symbolRanking)
  #We suppose that the ranking object contains information about ALL possible symbols
  if (is.na(ranking)) 
    stop("getQualitySubscore - Bad ranking assignation")
  
  if (ranking <= length(ratings))
    return(ratings[ranking])
  else
    return(0)
}

#
# Model Accuracy Score
#
# INPUT
#   - model : HMM or HSMM
#   - data: hsmmdata with multiple test observation sequences
#   - n: Accumulative parameter of the metric (length(data) must be greater than n)
#   - alpha: Quality-Timing balance (1 for HMMs)
#
# OUTPUT
#   - Model Accuracy Score
mhsmmHelper.MAS <- function(model, data, n=10, alpha = 1, aggregate = TRUE,format ="vector") {
  
  #Control error
  if (class(data)!="hsmm.data")
    stop("data must be a hsmm.data object, from mhsmm package")
  
  if (class(model)!="hmm" && class(model)!= "hsmm")
    stop("Model must be a HMM or a HSMM")
  
  prediction <- predict(model,data,method="viterbi")
  dataIndexes <- split(1:sum(data$N),rep(1:length(data$N),data$N))
  
  dataDf <- as.data.frame(data$x)
  MAS <- laply(dataIndexes, 
               function (sequenceIndexes) mhsmmHelper.singleMAS(model, 
                                                                dataDf %>% slice(sequenceIndexes),
                                                                prediction$s[sequenceIndexes],
                                                                format=format
                                                                ))
  if (aggregate)
    mean(MAS)
  else
    MAS
}

#
# Calculate the Model Accuracy Score of a single observation sequence, given a vector of state predictions
# for that sequence
#
mhsmmHelper.singleMAS <- function (model,data,mostLikelyStateSequence, n=10,alpha=1,format="vector") {
  
  dataDf <- as.data.frame(data)
  if (nrow(dataDf)!=length(mostLikelyStateSequence))
    stop("Observable data and state predictions must have the same length")
  
  symbolRankings <- llply(mostLikelyStateSequence, 
                          function (state) hmmHelper.getObservableRanking(model,state,format))
  
  dataDf$ranking <- symbolRankings
  qualitiesSubscores <- apply(dataDf, 1, function (x) mhsmmHelper.getQualitySubScore(symbol = t(x[-length(x)]),
                                                                                     symbolRanking = x$ranking))
  
  MAS  <- laply(n:nrow(dataDf), function (i) sum(alpha*qualitiesSubscores[(i-n+1):i])/n)
  mean(MAS)  
}