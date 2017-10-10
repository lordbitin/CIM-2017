###########################################
## Time handling strategies
###########################################


# Equidistant Time Slots - Replicate symbols
#Input : List of symbols! (each symbol can be a vector or a numeric value)
#Output: New list of symbols
timeHandling.equidistantSlots.replicate.hsmm <- function (obsList,timing,tsr) {
  
  #TODO: Ultima interaccion... Revisar bucle
  
  #print(obsSeqDf)
  
  newObsList <- list()
  equidistantSlots <- seq(head(timing,1),tail(timing,1), by=tsr)
  currentSymbol <- NULL
  
  i <- 1
  symbolOcurrenceQueue <- NULL
  for (t in 1:(length(equidistantSlots)-1)) {
    aux <- which(timing >=equidistantSlots[t] & timing < equidistantSlots[t+1])
    #print(aux)
    symbolOcurrenceQueue <- append(symbolOcurrenceQueue,aux)
    
    #The first iteration always enter here
    if (!is.na(symbolOcurrenceQueue[1])) {
      currentSymbol <- obsList[[symbolOcurrenceQueue[1]]]
    }
    
    newObsList[[i]] <- currentSymbol
    
    #Remove the first element of the ocurrence queue
    if (length(symbolOcurrenceQueue) > 1) {
      symbolOcurrenceQueue <- symbolOcurrenceQueue[2:length(symbolOcurrenceQueue)] 
    }
    else 
      symbolOcurrenceQueue <- NULL
    
    i <- i+1
  }
  
  #Last symbol??
  return(newObsList)
}

## Equidistant Time Slots - Silence symbols
## Input : List of symbols! (each symbol can be a vector or a numeric value)
## Output: New list of symbols
timeHandling.equidistantSlots.sylenceSymbols.hsmm <- function (obsList,timing,tsr) {
  
  SilenceSymbol <- length(Constants$INTERACTION_NAMES) #Starts at 0!!
  
  newObsList <- numeric(floor(tail(timing,1)/tsr)+1)
  newObsList[] <- NA
  
  for (i in 1:length(timing)) {
    slot <- floor(timing[[i]]/tsr) + 1
    while(slot <= length(newObsList) && !is.na(newObsList[[slot]])) slot <- slot+1 #Avoid symbol collision
    
    newObsList[[slot]] <- obsList[[i]] #TODO: Estamos obviando colisiones en el ultimo intervalo
  }
  
  #Fill the rest of the interactions with empty symbols
  newObsList[is.na(newObsList)] <- SilenceSymbol
  
  return(newObsList)
}

###################################################################################
## Equidistant Time Slots - NA
## Input : List of symbols! (each symbol can be a vector or a numeric value)
## Output: New list of symbols
###################################################################################
timeHandling.equidistantSlots.na.hsmm <- function (obsList,timing,tsr) {
  newObsList <- numeric(floor(tail(timing,1)/tsr)+1)
  newObsList[] <- NA
  
  for (i in 1:length(timing)) {
    slot <- floor(timing[[i]]/tsr) + 1
    while(slot <= length(newObsList) && !is.na(newObsList[[slot]])) slot <- slot+1 #Avoid symbol collision
    
    newObsList[[slot]] <- obsList[[i]] #TODO: Estamos obviando colisiones en el ultimo intervalo
  }
  
  return(newObsList)
  
}


##
# METHOD: sequenceHelper.fillGaps 
#         NOTES: Leftmost states stay as missing values
##
# Input:  seq: Vector of sequences with NA
#         mode: String. "replicate" (default) replicate the last symbol. "silence" adds a special silence symbol
#         silenceSymbol= Character. In case mode="silence", it defines the silence symbol
#
# Output: An object of class stslist (from TraMineR package)
#         
#
sequenceHelper.fillGaps <- function (seqs,mode="replicate",silenceSymbol='S') {
  #TODO Mode replicate, silence (Currently onluy replication works)
  aux <- sapply(seqs,function (subject) {
    result <- NULL
    lastValidState <- NA
    for (state in seqdecomp(subject)) {
      if (!is.na(state)) {
        lastValidState <- state
      }
      if (is.null(result)) {
        result <- lastValidState
      } else {
        result <- paste(result,lastValidState,sep="-")
      }
    }
    result
  })
  
  aux
}