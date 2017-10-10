
##
# Too short durations
##
outlier.dwr.isSimulationTooShort <- function (simulationID, 
                                              durationThreshold = Constants$OUTLIER_REAL_TIME_THRESHOLD) {
  return(dataHelper.getRealDuration(simulationID = simulationID) < durationThreshold)
}

##
# Too short durations
##
outlier.dwr.fewInteractions <- function (simulationID, 
                                         interactionThreshold = Constants$INTERACTIONS_THRESHOLD) {
  
  return(nrow(dataHelper.getSimulationInteractions(simulationID = simulationID)) < interactionThreshold)
}

##
# 
##
outlier.dwr.noIncidents <- function (simulationID) {
  incidentSnaps <- data.simulationSnapshots %>% 
    filter(simulation==simulationID,
           cause.id == SnapshotCauses$INCIDENT_STARTED
    )
  
  return(nrow(incidentSnaps) == 0)
}



##
# The simulation is a tutorial. We detect tutorials cuase the got more than two instructions sent
##
outlier.dwr.isTutorial <- function(simulationID) {
  instructionsSent <- data.simulationSnapshots %>%
    filter(
      simulation == simulationID,
      cause.id == SnapshotCauses$INSTRUCTION_SENT
    )
  return(nrow(instructionsSent) > 2)
}
