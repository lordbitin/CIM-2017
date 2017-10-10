#
# MultiUAV Grammar for HSMM creation
# Input:  snapsDf: Dataframe of simulation snapshots
#         uavs: Vector with UAV description (Ids and mode "ALL")
#         TSR: time slot resolution (ms)
# Output:
#         Object of class "hsmm.data" containing the interactions contained in "snapsDf" properly parsed
#
parse.grammar.dwr.multiUAV <- function (snapsDf,uavs=c(1,2,3,"ALL"), TSR = 1000) {
  # Columns - Monitoring, waypoint handling, replanning
  GRAMMAR_OPERATIONS <- c("MONITORING","WAYPOINT_HANDLING","REPLANNING")
  GRAMMAR_OPERANDS <- uavs
  GRAMMAR_SIZE <- (length(GRAMMAR_OPERATIONS)*(length(uavs)-1)) + 1
  
  #Input data
  #missionInteractions <- dwrHelper.getMissionInteractionSnapshots(missionID)
  missionInteractions <- snapsDf %>% 
    filter(cause.id == SnapshotCauses$USER_INPUT,elapsedRealTime!=0) %>% 
    filter(!is.na(cause.params.inputId)) %>%
    group_by(simulation) %>%
    filter(n() > Constants$INTERACTIONS_THRESHOLD) %>%
    arrange(simulation,elapsedRealTime) %>%
    ungroup
  
  #Grammar logic
  currentOperation <- NULL
  aux <- dlply(missionInteractions,.(simulation),function (snapshots) {
    
    interactions <- numeric(nrow(snapshots))
    for (i in 1:nrow(snapshots)) {
      ss <- snapshots[i,]
      
      #MONITORING OPERATION
      if (ss$cause.params.inputId == UserInputs$SELECT_DRONE | 
          ss$cause.params.inputId == UserInputs$SET_DRONE_SPEED | 
          ss$cause.params.inputId == UserInputs$SET_SIMULATION_TIME_RATIO |
          ss$cause.params.inputId == UserInputs$SET_WAYPOINT_SELECTION) {
        
        currentOperation <- "MONITORING"
        
        if (ss$cause.params.inputId == UserInputs$SELECT_DRONE)
          currentSelectedUAV <- ss$cause.params.droneId
        else
          currentSelectedUAV <- dataHelper.getCurrentSelectedUAV(ss)
        
        if (is.null(currentSelectedUAV)) currentSelectedUAV <- "ALL"
      }
      
      #WAYPOINT HANDLING AND REPLANNING OPERATIONS      
      else if (ss$cause.params.inputId == UserInputs$CHANGE_DRONE_PATH) {
        currentSelectedUAV <- ss$cause.params.droneId
        
        # Type of interactions in change drone path
        currentControlMode <- dataHelper.getControlMode(ss)
        if (currentControlMode == ControlModes$MONITOR) 
          currentOperation <- "WAYPOINT_HANDLING"
        else
          currentOperation <- "REPLANNING"
        
      }
      
      else if (ss$cause.params.inputId == UserInputs$SET_CONTROL_MODE) {
        currentSelectedUAV <- dataHelper.getCurrentSelectedUAV(ss)
        if (ss$cause.params.controlMode.id == ControlModes$MONITOR) 
          currentOperation <- "MONITORING"
        else if (ss$cause.params.controlMode.id == ControlModes$ADD_WAYPOINT)
          currentOperation <- "WAYPOINT_HANDLING"
        else if (ss$cause.params.controlMode.id == ControlModes$MANUAL) 
          currentOperation <- "REPLANNING"
      }
      
# print(paste("Interaction:",ss$cause.params.inputId))
# print(paste(currentSelectedUAV,currentOperation))
# print(paste(currentSelectedUAV,which(GRAMMAR_OPERANDS == currentSelectedUAV)))
# print(paste(currentOperation,which(GRAMMAR_OPERATIONS == currentOperation)))

      interactions[i] <- (which(GRAMMAR_OPERANDS == currentSelectedUAV)-1)*length(GRAMMAR_OPERATIONS) + 
        (which(GRAMMAR_OPERATIONS == currentOperation)-1)
    }
    return(interactions)
  })
  
  names(aux) <- NULL
  interactionList <- unlist(aux)
  
  #Create the mhsmm object
  numberOfObservationSequences <- as.numeric(table(missionInteractions$simulation))
  data <- mhsmmHelper.formatMultinomialData(data = interactionList,
                                                numberOfSequences = numberOfObservationSequences)
  
  #time handling
  data <- mhsmmHelper.timeHandling(data,
                                   timing = missionInteractions$elapsedRealTime,
                                   timeHandlingStrategy = timeHandling.equidistantSlots.replicate.hsmm,
                                   TSR = TSR)
  
  #Replicate interactions dataset
  return(list(data=data,size=GRAMMAR_SIZE))
}

##
# NOTE: Basic HSMM Grammar (Replicate Symbol)
##
# Input:  snapsDf: Dataframe of simulation snapshots
#         TSR: time slot resolution (ms)
# Output:
#         Object of class "hsmm.data" containing the interactions contained in "snapsDf" properly parsed
#
parse.grammar.dwr.basic <- function (snapsDf, TSR=1000) {
  
  GRAMMAR_SIZE = 6
  
  interactions <- snapsDf %>% 
    filter(cause.id == SnapshotCauses$USER_INPUT,elapsedRealTime!=0) %>% 
    filter(!is.na(cause.params.inputId)) %>%
    group_by(simulation) %>%
    filter(n() > Constants$INTERACTIONS_THRESHOLD) %>%
    arrange(simulation,elapsedRealTime) %>%
    ungroup %>%
    dplyr::select(cause.params.inputId,simulation,elapsedRealTime)
  
  
  #TODO : Arreglar esto...
  aux <- laply(interactions$cause.params.inputId, function (x) if (x>3) x-1 else x)
  interactions$cause.params.inputId <- aux
  
  #Create the mhsmm object
  numberOfObservationSequences <- as.numeric(table(interactions$simulation))
  data <- mhsmmHelper.formatMultinomialData(data = interactions$cause.params.inputId,
                                                numberOfSequences = numberOfObservationSequences)
  
  #Time handling
  data <- mhsmmHelper.timeHandling(data,
                                   timing = interactions$elapsedRealTime,
                                   timeHandlingStrategy = timeHandling.equidistantSlots.na.hsmm,
                                   TSR = TSR)
  
  return(list(data=data,size=GRAMMAR_SIZE))
}

##
# METHOD: Sylence symbol grammar for HSMM using mhsmm package
##
# Input:  snapsDf: Dataframe of simulation snapshots
#         TSR: time slot resolution (ms)
# Output:
#         Object of class "hsmm.data" containing the interactions contained in "snapsDf" properly parsed
#
parse.grammar.dwr.silenceSymbols <- function (snapsDf, TSR = TSR) {
  
  GRAMMAR_SIZE = 6
  
  interactions <- snapsDf %>% 
    filter(cause.id == SnapshotCauses$USER_INPUT,elapsedRealTime!=0) %>% 
    filter(!is.na(cause.params.inputId)) %>%
    group_by(simulation) %>%
    filter(n() > Constants$INTERACTIONS_THRESHOLD) %>%
    arrange(simulation,elapsedRealTime) %>%
    ungroup %>%
    dplyr::select(cause.params.inputId,simulation,elapsedRealTime)
  
  #TODO : Arreglar esto...
  aux <- laply(interactions$cause.params.inputId, function (x) if (x>3) x-1 else x)
  interactions$cause.params.inputId <- aux
  
  #Create the mhsmm object
  numberOfObservationSequences <- as.numeric(table(interactions$simulation))
  data <- mhsmmHelper.formatMultinomialData(data = interactions$cause.params.inputId,
                                            numberOfSequences = numberOfObservationSequences)
  
  #Time handling
  data <- mhsmmHelper.timeHandling(data,
                                   timing = interactions$elapsedRealTime,
                                   timeHandlingStrategy = timeHandling.equidistantSlots.sylenceSymbols.hsmm,
                                   TSR = TSR)
  
  return(list(data=data,size=GRAMMAR_SIZE+1))
}

##
# METHOD: parse.grammar.dwr.basic.seq: Create sequence objects with the interactions contained in 
#         a snapshots dataframe from DWR
#         NOTA: ESTE METODO TIENE EN CUENTA TODOS LOS SNAPSHOTS, NO SOLO LOS QUE ESTAN EN SIMULATIONS
##
# Input:  snapsDf: Dataframe of simulation snapshots from DWR
#         TSR: time slot resolution (ms) to sample the resulting sequence
# Output:
#         An object of class stslist (from TraMineR package)
#
parse.grammar.dwr.lastInteraction.seq <- function (snapsDf, TSR, useRealTime=FALSE, compress=TRUE) {

  #Data Frame
  data <- .dwrGrammar.interactionSnaps(snapsDf)
  
  if (useRealTime) {
    data <- data %>% rename(time=elapsedRealTime)
  } else {
    data <- data %>% rename(time=elapsedSimulationTime)
  }
  data <- data %>% 
    dplyr::select(simulation,time,cause.id,cause.params.inputId) %>%
    mutate(index=.dwrGrammar.getSnapshotIndex(time,TSR))
  
  rawSeqData <- seqformat(data=data,
                    id="simulation",
                    from = "SPELL",
                    to="STS",
                    status="cause.params.inputId",
                    begin="index",
                    end="index",
                    process=FALSE,
                    tmin=1,
                    tmax=max(data$index),
                    compressed = FALSE)
  
  compressedSeqData <- apply(rawSeqData,1,function (row) paste(row,collapse="-"))
  
  #Replicar el ultimo simbolo
  replications <- sequenceHelper.fillGaps(compressedSeqData,mode="replicate")
  seq <- seqdef(data.frame(Sequence=unlist(replications)),
                right=NA,
                alphabet=c(0,1,2,4,5,6),
                states=c("SelDrone","SetSpeed","ChangeTime","ChangePath","W.Table","ControlMode"),
                labels=c("Select UAV",
                         "Change UAV Speed",
                         "Change Simulation Speed",
                         "Chage UAV Path",
                         "Modify waypoints table",
                         "Change Control Mode"),
                xtstep = 100
  )
  
  #TODO: Tratar indices repetidos!
  
  
  return(seq)
}

##
# METHOD: parse.grammar.dwr.lastMissionEvent.seq: Parse a dataset of simulation sjnapshots in DWR to create
#         a sequence (stslist object from package TraMineR) ofwith the information of the "last Mission Event" happened
#         in the mission. Mission events do not include user interactions
##
# Input:  snapsDf: Dataframe of simulation snapshots from DWR
#         TSR: time slot resolution (ms) to sample the resulting sequence
# Output:
#         An object of class stslist (from TraMineR package)
#
parse.grammar.dwr.lastMissionEvent.seq <- function (snapsDf, TSR, useRealTime=FALSE, compress=TRUE) {
  
  #Data Frame
  data <- .dwrGrammar.interactionSnaps(snapsDf)
  
  if (useRealTime) {
    data <- data %>% rename(time=elapsedRealTime)
  } else {
    data <- data %>% rename(time=elapsedSimulationTime)
  }
  data <- data %>% 
    dplyr::select(simulation,time,cause.id,cause.params.inputId) %>%
    mutate(index=.dwrGrammar.getSnapshotIndex(time,TSR))
  
  rawSeqData <- seqformat(data=data,
                          id="simulation",
                          from = "SPELL",
                          to="STS",
                          status="cause.params.inputId",
                          begin="index",
                          end="index",
                          process=FALSE,
                          tmin=1,
                          tmax=max(data$index),
                          compressed = FALSE)
  
  compressedSeqData <- apply(rawSeqData,1,function (row) paste(row,collapse="-"))
  
  #Replicar el ultimo simbolo
  replications <- sequenceHelper.fillGaps(compressedSeqData,mode="replicate")
  seq <- seqdef(data.frame(Sequence=unlist(replications)),
                right=NA,
                alphabet=c(0,1,2,4,5,6),
                states=c("SelDrone","SetSpeed","ChangeTime","ChangePath","W.Table","ControlMode"),
                labels=c("Select UAV",
                         "Change UAV Speed",
                         "Change Simulation Speed",
                         "Chage UAV Path",
                         "Modify waypoints table",
                         "Change Control Mode"),
                xtstep = 100
  )
  
  #TODO: Tratar indices repetidos!
  
  
  return(seq)
}


##
# METHOD: parse.grammar.dwr.interactionsTransition.seq
#
# Get the sequence of interactions WITHOUT REPLICATION
#
# Input:  data: Dataframe of simulation snapshots
#         
# Output:
#         Dataframe containing the interactions of every simulation in the input. Ordered by simulation and elapsedRealtime
#
parse.grammar.dwr.interactionsTransition.seq <- function (data) {
  interactionsSubjectData <- dlply(data,.(simulation), function (subject) {
    
    interactions <- subject %>% 
      filter(
        !is.na(cause.id),
        cause.id == SnapshotCauses$USER_INPUT,
        elapsedRealTime!=0,
        !is.na(cause.params.inputId)
      ) %>%
      filter(
        cause.params.inputId != UserInputs$ACCEPT_INSTRUCTION & 
          cause.params.inputId != UserInputs$ADD_WAYPOINT
      ) %>%
      arrange(elapsedRealTime) %>%
      dplyr::select(simulation,elapsedRealTime,cause.params.inputId)
    
    interactions$index <- 1:nrow(interactions)

    seqInteractions <- seqformat(data=interactions,
                                 id="simulation",
                                 from = "SPELL",
                                 to="STS",
                                 status="cause.params.inputId",
                                 begin="index",
                                 end="index",
                                 process=FALSE,
                                 tmin=1,
                                 tmax=max(interactions$index),
                                 compressed = FALSE)

    seqInteractions <- apply(seqInteractions,1,function (row) paste(row,collapse="-")) #Compress sequences

    list(interactions=seqInteractions)
  })
  
  # Merge the results with rbind.fill (adjust subject sequence sizes)
  interactionsTransition <- do.call(rbind.fill,
                             llply(interactionsSubjectData, 
                                   function (subject) as.data.frame(seqdecomp(subject$interactions))))
  
  row.names(interactionsTransition) <- names(interactionsSubjectData)
  
  
  seqdef(interactionsTransition,
                          #right=NA,
                          alphabet=c(0,1,2,4,5,6),
                          states=c("SD",
                                   "SS",
                                   "ST",
                                   "CP",
                                   "MT",
                                   "CM"),
                          labels=c("Select UAV",
                                   "Change UAV Speed",
                                   "Change Simulation Speed",
                                   "Change UAV Path",
                                   "Modify waypoints table",
                                   "Change Control Mode")
  )
}



##
# METHOD: parse.grammar.dwr.lastState.multichannel.seq
##
# Input:  data: Dataframe of simulation snapshots
#         TSR: Time Slot Resolution (ms)
#         userRealTime: Whether to use elapsedRealtime or elapsedSimulationTime in the time axis
#         missing_left: What to do with left missing values (See TramineR::seqdef, left argument)
#         missing_right: What to do with right missing values (See TramineR::seqdef, right argument)
#         
# Output:
#         Dataframe containing the interactions of every simulation in the input. Ordered by simulation and elapsedRealtime
#
parse.grammar.dwr.lastState.multichannel.seq <- function(data,TSR,useRealTime=FALSE, missing_left = NA, missing_right = "DEL") {

  if (useRealTime) {
    data <- data %>% rename(time=elapsedRealTime)
  } else {
    data <- data %>% rename(time=elapsedSimulationTime)
  }
  
  multiChannelSubjectData <- dlply(data,.(simulation), function (subject) {
    
    interactions <- subject %>% 
      filter(
              !is.na(cause.id),
              cause.id == SnapshotCauses$USER_INPUT,
              time!=0,
              !is.na(cause.params.inputId)
             ) %>%
      filter(
        cause.params.inputId == UserInputs$SELECT_DRONE |
        cause.params.inputId == UserInputs$SET_SIMULATION_TIME_RATIO |
        cause.params.inputId == UserInputs$CHANGE_DRONE_PATH |
        cause.params.inputId == UserInputs$SET_CONTROL_MODE |
        cause.params.inputId == UserInputs$SET_WAYPOINT_SELECTION | 
        cause.params.inputId == UserInputs$SET_DRONE_SPEED
      ) %>%
      arrange(time) %>%
      dplyr::select(simulation,time,cause.params.inputId) %>%
      mutate(index = .dwrGrammar.getSnapshotIndex(time,TSR = TSR))
    
    missionEvents <- subject %>%
      filter(
              #cause.id == SnapshotCauses$DRONE_REACH_WAYPOINT | 
              cause.id == SnapshotCauses$INCIDENT_STARTED |
             cause.id ==SnapshotCauses$INCIDENT_ENDED |
             cause.id ==SnapshotCauses$TARGET_DETECTED |
             cause.id ==SnapshotCauses$DRONE_DESTROYED |
             cause.id == SnapshotCauses$ACTION_STARTED
             #cause.id == SnapshotCauses$ACTION_ENDED
      ) %>% 
      arrange(time) %>%
      dplyr::select(simulation,time,cause.id) %>%
      mutate(index = .dwrGrammar.getSnapshotIndex(time,TSR = TSR))
    
    #Create the sequence objects for interactions and mission events (multichannel sequence)
    subjectMaxSnapshotIndex <- max(max(interactions$index),max(missionEvents$index))

    seqInteractions <- seqformat(data=interactions,
                                 id="simulation",
                                 from = "SPELL",
                                 to="STS",
                                 status="cause.params.inputId",
                                 begin="index",
                                 end="index",
                                 process=FALSE,
                                 tmin=1,
                                 tmax=subjectMaxSnapshotIndex,
                                 compressed = FALSE)
    seqInteractions <- apply(seqInteractions,1,function (row) paste(row,collapse="-")) #Compress sequences
    #Replicate last symbol
    seqInteractions <- sequenceHelper.fillGaps(seqInteractions,mode="replicate")
    
    seqMissionEvents <- seqformat(data=missionEvents,
                                  id="simulation",
                                  from = "SPELL",
                                  to="STS",
                                  status="cause.id",
                                  begin="index",
                                  end="index",
                                  process=FALSE,
                                  tmin=1,
                                  tmax=subjectMaxSnapshotIndex,
                                  compressed = FALSE)
    seqMissionEvents <- apply(seqMissionEvents,1,function (row) paste(row,collapse="-")) #Compress sequences  
    #Replicate last symbol
    seqMissionEvents <- sequenceHelper.fillGaps(seqMissionEvents,mode="replicate")

    list(lastInteraction=seqInteractions,
         lastMissionEvent=seqMissionEvents)
  })
  
  # Merge the results with rbind.fill (adjust subject sequence sizes)
  lastInteraction <- do.call(rbind.fill,llply(multiChannelSubjectData, function (subject) as.data.frame(seqdecomp(subject$lastInteraction))))
  lastMissionEvent <- do.call(rbind.fill,llply(multiChannelSubjectData, function (subject) as.data.frame(seqdecomp(subject$lastMissionEvent))))
  
  row.names(lastInteraction) <- names(multiChannelSubjectData)
  row.names(lastMissionEvent) <- names(multiChannelSubjectData)
  
  list(
    lastInteraction= seqdef(lastInteraction,
                            left = missing_left,
                            right=missing_right,
                             alphabet=c(0,1,2,4,5,6),
                              #alphabet=c(0,2,4,6),
                             states=c("SD",
                                      "SS",
                                      "ST",
                                      "CP",
                                      "MT",
                                      "CM"),
                             labels=c("Select UAV",
                                      "Set UAV Speed",
                                      "Change Simulation Speed",
                                      "Change UAV Path",
                                      "Modify waypoints table",
                                      "Change Control Mode"),
                             xtstep = 100
    ),
    lastMissionEvent= seqdef(lastMissionEvent,
                             left = missing_left,
                             right= missing_right,                              
                             alphabet=c(
                                          2,
                                          3,
                                          4,
                                          5,
                                          6
                                          #7
                                          ),
                              states=c(
                                        "IS",
                                        "IE",
                                        "TD",
                                        "DD",
                                        "AS"
                                        #"AE"Target 
                                        ),
                              labels=c(
                                       "Incident Started",
                                       "Incident Ended",
                                       "Target Detected",
                                       "Drone Destroyed",
                                       "Action Started"
                                       #"Action Ended"
                                        ),
                              xtstep = 100
    )
  )
}

##
# METHOD: parse.grammar.dwr.lastState.interactions.seq
##
# Input:  data: Dataframe of simulation snapshots
#         TSR: Time Slot Resolution (ms)
#         userRealTime: Whether to use elapsedRealtime or elapsedSimulationTime in the time axis
#         
# Output:
#         Dataframe containing the interactions of every simulation in the input. Ordered by simulation and elapsedRealtime
#
parse.grammar.dwr.lastState.interactions.seq <- function(data,TSR,useRealTime=FALSE) {
  
  if (useRealTime) {
    data <- data %>% rename(time=elapsedRealTime)
  } else {
    data <- data %>% rename(time=elapsedSimulationTime)
  }
  
  multiChannelSubjectData <- dlply(data,.(simulation), function (subject) {
    
    interactions <- subject %>% 
      filter(
        !is.na(cause.id),
        cause.id == SnapshotCauses$USER_INPUT,
        time!=0,
        !is.na(cause.params.inputId)
      ) %>%
      arrange(time) %>%
      dplyr::select(simulation,time,cause.params.inputId) %>%
      mutate(index = .dwrGrammar.getSnapshotIndex(time,TSR = TSR))
    
    #Create the sequence objects for interactions and mission events (multichannel sequence)
    subjectMaxSnapshotIndex <- max(interactions$index)
    
    seqInteractions <- seqformat(data=interactions,
                                 id="simulation",
                                 from = "SPELL",
                                 to="STS",
                                 status="cause.params.inputId",
                                 begin="index",
                                 end="index",
                                 process=FALSE,
                                 tmin=1,
                                 tmax=subjectMaxSnapshotIndex,
                                 compressed = FALSE)
    seqInteractions <- apply(seqInteractions,1,function (row) paste(row,collapse="-")) #Compress sequences
    #Replicate last symbol
    seqInteractions <- sequenceHelper.fillGaps(seqInteractions,mode="replicate")
    
    list(lastInteraction=seqInteractions)
  })
  
  # Merge the results with rbind.fill (adjust subject sequence sizes)
  lastInteraction <- do.call(rbind.fill,
                             llply(multiChannelSubjectData, 
                                   function (subject) as.data.frame(seqdecomp(subject$lastInteraction))))
  
  row.names(lastInteraction) <- names(multiChannelSubjectData)

  lastInteraction= seqdef(lastInteraction,
                          #right=NA,
                          alphabet=c(0,1,2,4,5,6),
                          states=c("SD",
                                   "SS",
                                   "ST",
                                   "CP",
                                   "MT",
                                   "CM"),
                          labels=c("Select UAV",
                                   "Change UAV Speed",
                                   "Change Simulation Speed",
                                   "Change UAV Path",
                                   "Modify waypoints table",
                                   "Change Control Mode"),
                          xtstep = 100
  )
  
  lastInteraction
}


##
# METHOD: interactionSnaps
##
# Input:  snapsDf: Dataframe of simulation snapshots
#         
# Output:
#         Dataframe containing the interactions of every simulation in the input. Ordered by simulation and elapsedRealtime
#
.dwrGrammar.interactionSnaps <- function (snapsDf) {
  interactions <- snapsDf %>% 
    filter(cause.id == SnapshotCauses$USER_INPUT,elapsedRealTime!=0) %>% 
    filter(!is.na(cause.params.inputId)) %>%
    group_by(simulation) %>%
    filter(n() > Constants$INTERACTIONS_THRESHOLD) %>%
    arrange(simulation,elapsedRealTime) %>%
    ungroup
}

#
# OUTPUT : Vector with timegeneral time betwen interactions in DWR
#
.dwrGrammar.getTimeBetweenInteractions <- function (interactionsDf) {
  
  l <- dlply(interactionsDf,
             .(simulation),
             function (x) {diff(x$elapsedRealTime)})
  
  return(unlist(l,use.names = FALSE))
}

##
# METHOD: getSnaphotIndex
##
# Input:  snapTime: Integer. Timestamp of a snapshot
#         TSR: Integer. Time Slot Resolution
#         startTime: Integer. Offset for the index creation
#         
# Output:
#         Positive Integer indicating the index of the snapshot
#
.dwrGrammar.getSnapshotIndex <- function(snapTime,TSR,startTime=0) {
  return(
    ceiling(snapTime/TSR)
  )
} 
