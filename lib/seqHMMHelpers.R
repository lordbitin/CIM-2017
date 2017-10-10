#
# Custom BIC?
#
seqHMMHelper.customBIC <- function (object, ...) {
  #ll <- if (isNamespaceLoaded("stats4")) 
  #  stats4::logLik
  #else logLik
  ll <- seqHMM:::logLik.hmm
  Nobs <- if (isNamespaceLoaded("stats4")) 
    stats4::nobs
  else nobs
  #print(ll)
  #print(nobs)
  if (!missing(...)) {
    lls <- lapply(list(object, ...), ll)
    vals <- sapply(lls, function(el) {
      no <- attr(el, "nobs")
      c(as.numeric(el), attr(el, "df"), if (is.null(no)) NA_integer_ else no)
    })
    val <- data.frame(df = vals[2L, ], ll = vals[1L, ], 
                      nobs = vals[3L, ])
    nos <- na.omit(val$nobs)
    if (length(nos) && any(nos != nos[1L])) 
      warning("models are not all fitted to the same number of observations")
    unknown <- is.na(val$nobs)
    if (any(unknown)) 
      val$nobs[unknown] <- sapply(list(object, ...)[unknown], 
                                  function(x) tryCatch(Nobs(x), error = function(e) NA_real_))
    val <- data.frame(df = val$df, BIC = -2 * val$ll + log(val$nobs) * 
                        val$df)
    row.names(val) <- as.character(match.call()[-1L])
    val
  }
  else {
    lls <- ll(object)
    nos <- attr(lls, "nobs")
    if (is.null(nos)) 
      nos <- tryCatch(Nobs(object), error = function(e) NA_real_)
    -2 * as.numeric(lls) + log(nos) * attr(lls, "df")
  }
}


#
# Direct computation of BIC from log likelihood, number of observations and degrees of freedom
#
seqHMMHelper.customBIC_direct <- function (ll,nobs,df) {
  -2 *ll + log(nobs) * df
}


#
# Compute the BIC for the best models resulting of a training
#
seqHMMHelper.best_opt_restart_BIC <- function (fit_model_results, 
                                               opt_threshold = length(fit_model_results$em_results$best_opt_restart)) {
  
  laply(fit_model_results$em_results$best_opt_restart[1:opt_threshold], 
        function (optLL) seqHMMHelper.customBIC_direct(
          optLL,
          attr(fit_model_results$model,"nobs"),
          attr(fit_model_results$model,"df")
          )
        )
}

#
# Train and select a HMM using the package seqHMM.
# seqData : sequence (TramineR)
# nstates: number of states (or a vector for multiple training) 
#
# Output: List of trained HMMs, one for each possible state
seqHMMHelper.trainHMMs <- function (seqData,nstates,control_em=list(),threads=1, verbose=FALSE, returnTrainingInfo = FALSE) {
  #Create a set of models with different number of states
  hmmCandidates <- llply(nstates, function (J) {
    mc_init <- hmmHelper.randomInitialVector(J)
    mc_trans <- hmmHelper.randomTransitionMatrix.hmm(J)
    
    #TODO: Las emisiones se estan inicializando de forma aleatoria...hacerlo mejor?
    if (!is.data.frame(seqData)) {
      #mc_emissions <- llply(seqData, function (channel) hmmHelper.emissionMatrix.seq(channel,J))    
      mc_emissions <-llply(seqData, function (channel) {
        simulate_emission_probs(n_states = J,
                                n_symbols = length(alphabet(channel)),
                                n_clusters = 1)
      })
    } else {
      #mc_emissions <- hmmHelper.emissionMatrix.seq(seqData,J)
      mc_emissions <- simulate_emission_probs(n_states = J,
                                              n_clusters = 1,
                                              n_symbols = length(alphabet(seqData)))
    }
    
    mc_init_mod <- build_hmm(
      observations = seqData,
      initial_probs = mc_init,
      transition_probs = mc_trans,
      emission_probs = mc_emissions,
      channel_names = names(seqData)
    )
    
    if (verbose)
      print(paste("Fitting a",
                  mc_init_mod$n_states,
                  "state model, from",
                  mc_init_mod$n_sequences,
                  "sequences with",
                  mc_init_mod$n_channels,
                  "channels and",
                  mc_init_mod$n_symbols,
                  "symbols"))
    mc_fit <- fit_model(
      mc_init_mod,
      em_step = TRUE,
      control_em = control_em,
      threads = threads
    )
    
    #ifelse(returnTrainingInfo==TRUE,mc_fit,mc_fit)
    if (returnTrainingInfo) return(mc_fit)
    else return(mc_fit$model)
  }, .progress=ifelse(verbose,"text","none"))
  
  if (length(hmmCandidates) == 1) hmmCandidates[[1]]
  else hmmCandidates
}


##
# Train Dual HMM Classifier
##
seqHMMHelper.trainDualHMMClassifier <- function (seqData, 
                                                dataLabels, 
                                                possibleStates, 
                                                control_em = list(), 
                                                threads = 1,
                                                verbose = FALSE) {
  
  hmms <- llply(unique(dataLabels), function (classLabel) {
    if (!is.data.frame(seqData)) {
      classData <- llply(seqData,function (channel) channel[dataLabels == classLabel,])
      notClassData <- llply(seqData,function (channel) channel[dataLabels != classLabel,])
    } else {
      classData <- seqData[dataLabels == classLabel,]
      notClassData <- seqData[dataLabels != classLabel,]
    }
    
    if (verbose) print(paste("Training HMMs for class",classLabel))
    classHMMs <- seqHMMHelper.trainHMMs(seqData = classData,
                                       nstates=possibleStates,
                                       control_em = control_em,
                                       threads = threads,
                                       verbose=verbose)
    
    if (verbose) print(paste("Training HMMs for NOT class",classLabel))
    notClassHMMs <- seqHMMHelper.trainHMMs(seqData = notClassData,
                                          nstates=possibleStates,
                                          control_em = control_em,
                                          threads = threads,
                                          verbose=verbose)
    
    list(YES=classHMMs,NO=notClassHMMs)
  }, .progress=ifelse(verbose,"text","none")
  )
  names(hmms) <- unique(dataLabels)
  
  # Model selection
  if (verbose) print("************* Selecting the best HMMs via BIC metric ******************")
  dualHMMClassifier <- llply(hmms,function (classModels) {
    
    #TODO: La seleccion de modelos no funciona si se entrenan modelos con un solo numero de estados
    classModels$YES <- classModels$YES[[which.min(laply(classModels$YES, function (model) {
      seqHMMHelper.customBIC(object=model)
    }))
    ]]
    classModels$NO <- classModels$NO[[which.min(laply(classModels$NO, function (model) {
      seqHMMHelper.customBIC(object=model)
    }))
    ]]
    
    # llply(classModels, function (models) {
    #   models[[which.min(laply(models, function (model) {
    #       seqHMMHelper.customBIC(object=model)
    #     }))
    #   ]]
    # })
    classModels
  })
  
  dualHMMClassifier
}



##
# LOOCV
##
seqHMMHelper.LOOCV <- function (seqData, dataLabels, possibleStates, control_em, threads, verbose) {
  result <- llply(names(dataLabels), function (simulation) {
    
    #Remove a simulation from the analysis, to use it as test
    if (!is.data.frame(seqData)) {
      testData <- llply(seqData, function (channel) channel[rownames(channel) %in% simulation,])
      trainingData <- llply(seqData, function (channel) channel[!(rownames(channel) %in% simulation),])
    } else {
      testData <- seqData[rownames(seqData) %in% simulation,]
      trainingData <- seqData[!(rownames(seqData) %in% simulation),]
    }
    
    #Estimate the HMM classifier accurately with this subdataset
    print(paste("###### Training Dual HMM Classifier ######"))
    dualHMMClassifier <- seqHMMHelper.trainDualHMMClassifier(seqData = trainingData,
                                                            dataLabels = dataLabels[names(dataLabels) != simulation],
                                                            possibleStates = NSTATES,
                                                            control_em = CONTROL_EM,
                                                            threads = THREADS,
                                                            verbose = VERBOSE)
    dualHMMClassifier
    
    list(
      model=dualHMMClassifier,
      trainingData=trainingData,
      testData=testData
    )
    
  }, .progress="text")
}

##
# State frequency distribution
# Input: Object of class "hmm" from the package seqHMM
##
seqHMMHelper.state_frequency_distribution <- function (hmm) {
  hp <- hidden_paths(model=hmm)
  freqs <- seqstatd(seqdata = hp)$Frequencies
  
  return(apply(freqs,1,mean))
}

##
# rare states of an hmm
# Input: hmm (seqHMM)
#        frequency_threshold : [0,1]
# Output: Named Vector (State -> frequency)
##
seqHMMHelper.rare_states <- function (hmm, frequency_threshold = 0.05) {
  aux <- seqHMMHelper.state_frequency_distribution(hmm)
  return(aux[aux < frequency_threshold])
}
