# Interactions Dataset
dwrHelper.getInteractionsDataset <- function () {
  interactions <- students.simulationSnapshots %>%
    filter(importedMissionPlan.id != 4) %>%
    filter(cause.id == SnapshotCauses$USER_INPUT) %>%
    dplyr::select(cause.params.inputId,simulation,elapsedRealTime) %>%
    filter(elapsedRealTime!=0) %>%
  group_by(simulation) %>%
  filter(n() > Constants$INTERACTIONS_THRESHOLD) %>%
    arrange(simulation,elapsedRealTime) %>%
    ungroup %>%
    rename(interactionId=cause.params.inputId)
  
  return(interactions)
}

#Get interaction snapshots of a given mission ID
dwrHelper.getMissionInteractionSnapshots <- function (missionID) {
  missionInteractions <- students.simulationSnapshots %>%
    filter(importedMissionPlan.id == missionID) %>%
    filter(cause.id == SnapshotCauses$USER_INPUT) %>%
    filter(elapsedRealTime != 0) %>%
    group_by(simulation) %>%
    filter(n() > Constants$INTERACTIONS_THRESHOLD) %>%
    arrange(simulation,elapsedRealTime) %>%
    ungroup
  
  return(missionInteractions)
}

dwrHelper.getInteractionSnaps <- function() {
  
}

#
# OUTPUT : Vector with timegeneral time betwen interactions in DWR
#
dwrHelper.getTimeBetweenInteractions <- function () {
    
    l <- dlply(dwrHelper.getInteractionsDataset(),
               .(simulation),
               function (x) {diff(x$elapsedRealTime)})
    
    return(unlist(l))
}

##
## INPUT : -
## OUTPUT: "hsmm.data" object, from the package mhsmm with the interactions from DWR prepared as training data
dwrHelper.getMHSMMTrainingData <- function(format = "vector") {
  
  #Prepare the dataframe of interactions
  interactions <- dwrHelper.getInteractionsDataset()
  
  #Test Multinomial HMM in package MHSMM (Model example with interactions Data)
  interactionNames <- Constants$INTERACTION_NAMES
  numberOfObservationSequences <- as.numeric(table(interactions$simulation))
  
  # Data in suitable format
  switch(format,
         vector={
           dwrTrain <- mhsmmHelper.formatMultinomialData(data = interactions$interactionId,
                                                         numberOfSequences = numberOfObservationSequences,
                                                         minOption = 0,
                                                         maxOption = 5)           
         },
         matrix={
           dwrTrain <- mhsmmHelper.formatMatrixMultinomialData(data = interactions$interactionId,
                                                         numberOfSequences = numberOfObservationSequences,
                                                         minOption = 0,
                                                         maxOption = 5)           
         },
         {
           stop("Bad format parameter")
         })
  
  dwrTrain$timing <- interactions$elapsedRealTime
  dwrTrain
}

# Fake MHSMM training data for fast developing
# INPUT
#   - K : Number of observation sequences to retrieve
dwrHelper.getFakeHSMMTrainingData <- function (K=3, format="vector") {
  completeDataset <- dwrHelper.getMHSMMTrainingData(format)
  
  if (K > length(completeDataset$N))
    stop("K must be lower than the number of observation sequences in the dataset (N)")
  
  # X
  if (class(completeDataset$x)=="matrix")
    completeDataset$x <- completeDataset$x[1:sum(completeDataset$N[1:K]),]
  else
    completeDataset$x <- completeDataset$x[1:sum(completeDataset$N[1:K])]
  
  #timing
  completeDataset$timing <- completeDataset$timing[1:sum(completeDataset$N[1:K])]
  
  #N
  completeDataset$N <- completeDataset$N[1:K]
  
  
  completeDataset
}