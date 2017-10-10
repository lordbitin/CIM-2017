# 
# This script loads data logs from the simulation environment "Drone Watch & Rescue", contained in the files data.simulations
# and data.simulationSnapshots, and processes them in order to be used by the package "marcj". An object of 
# class "march.Dataset" is created, containing the sequence of operator interactions for each simulation in the log. 
# The dataset is split into a training/test set
#
# Requirements:
#
# Parameters:
#   - LOAD_DATA: T/F
#   - SAVE_RESULTS: T/F
#
# Output:
#   - Plot of the distribution of the sequential data.
#   - objects: marchData, marchData.training, marchData.test

library(ProjectTemplate)
load.project()

LOAD_DATA <- TRUE
SAVE_RESULTS <- FALSE

#Load data into "march" format
if (LOAD_DATA) {
  seqData <- parse.grammar.dwr.interactionsTransition.seq(
    data = data.simulationSnapshots %>% 
      filter(
        simulation %in% data.simulations$simulationID
      )
  )
}

marchData <- march.dataset.loadFromDataFrame(
  dataframe = seqData,
  MARGIN = 2,
  weights = NA,
  missingDataRep = '%'
)

# Split the data set into train/test sets
set.seed(1234)
trainingRows <- sample(1:nrow(marchData@yRaw), 0.75*nrow(marchData@yRaw))
marchData.training <- march.dataset.loadFromDataFrame(marchData@yRaw[trainingRows,],MARGIN = 2, weights = NA,missingDataRep = '%')
marchData.test <- march.dataset.loadFromDataFrame(marchData@yRaw[-trainingRows,],MARGIN = 2, weights = NA,missingDataRep = '%')

#Output
if (SAVE_RESULTS)
{
  save(marchData, file = paste0("output/Markovian_Comparison/marchData_NO_TSR_MIN_INTERACTIONS_15","(",Sys.Date(),")",".RData"))
  save(marchData.training, file = paste0("output/Markovian_Comparison/marchData_training_NO_TSR_MIN_INTERACTIONS_15","(",Sys.Date(),")",".RData"))
  save(marchData.test, file = paste0("output/Markovian_Comparison/marchData_test_NO_TSR_MIN_INTERACTIONS_15","(",Sys.Date(),")",".RData"))
}

#Plot
seqHMM::ssplot(seqData,
       type = "d",
       #title = "State distribution plots, K = 55",
       title = FALSE,
       xlab = "Time Steps",
       ylab = "Last Interaction",
       border=NA,
       withlegend = "bottom",
       plots = "obs",
       with.missing = TRUE,
       legend.prop = 0.33,
       xlab.pos = 1.5,
       ncol.legend = 1,
       cex.legend = 1.6,
       cex.lab = 1.5,
       cex.title = 1.5,
       cex.axis = 1.6,
       title.n = FALSE
)