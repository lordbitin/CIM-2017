# 
# This script performs a Rank Aggregation process over a grid of evaluated DCMMs in order to create a final ranking
# that achieve a compromise between the predictability and the interpretability of the model.
#
# Requirements:
#     The object "gs" must exist, as a result of executing the script "2-Markovian_Comparison-DCMM"
#
# Parameters:
#     SAVE_OUTPUT: T/F
#     VIEW_PLOT: T/F
#     SAVE_PLOT: T/F
#     PREDICTABILITY_MEASURES: character vector. Names of the measures that conform the family of "predictability" measures
#     INTERPRETABILITY_MEASURES: character vector. Names of the measures that conform the family of "predictability" measures
#     K1: Size of the medium agrregated rankings
#     K2: Size of the final ranking
#     RANK_AGGREGATION_METHOD: Method used to performed rank aggregation (See RankAggreg package)
#     GA_MAX_ITER: Maximum number of iterations of the GA (See RankAggreg package)
#     GA_CONVERGENCE_ITERATIONS: Converege iterations of the GA (See RankAggreg package)
#     GA_POPULATION_SIZE: Population size in the GA (See RankAggreg package)
#     GA_CROSS_OVER_PROBABILITY: (See RankAggreg package)
#     GA_MUTATION_PROBABILITY: (See RankAggreg package)
#     GA_RERUN_TIMES: Number of times that the GA will be run, in order to choose the best solution.
#     GA_VERBOSE: T/F
#     PREDICTABILITY_IMPORTANCE: Real [0,1]. The Interpretability importance is computed as 1 - PREDICTABILITY_IMPORTANCE
#
# Output:
#     object raggr of class "raggr" (See RankAggreg package), containing the final DCMM ranking
#     Plot of the probability transition matrix of the hidden chain of the best DCMM found
#
library(ProjectTemplate)
load.project()

#Metaparameters (Control the script)
if (!exists("EXTERNAL_CALL")) {
  VIEW_PLOT <- TRUE
  SAVE_PLOT <- FALSE
  SAVE_OUTPUT <- FALSE
  SET_PARAMETERS <- TRUE
}

# Script requirements - The grid search data frame must exists
stopifnot(exists("gs"))

#Parameters
if (SET_PARAMETERS) {
  #PREDICTABILITY_MEASURES <- c("BIC","SLL","AGP", "MPGP")
  PREDICTABILITY_MEASURES <- c("BIC","AGP", "MPGP")
  #INTERPRETABILITY_MEASURES <- c("BIC","NRS",,"CHPVT")
  INTERPRETABILITY_MEASURES <- c("BIC","CHPHT","CHPVT")
  K1 <- 30
  K2 <- 5
  RANK_AGGREGATION_METHOD <- "GA"
  GA_MAX_ITER <- 2000
  GA_CONVERGENCE_ITERATIONS <- 15
  GA_POPULATION_SIZE <- 100
  GA_CROSS_OVER_PROBABILITY <- 0.4
  GA_MUTATION_PROBABILITY <- 0.01
  GA_RERUN_TIMES <- 20
  GA_VERBOSE <- FALSE
  PREDICTABILITY_IMPORTANCE <- 0.65
  INTERPRETABILITY_IMPORTANCE <- 1 - PREDICTABILITY_IMPORTANCE
}

# LOAD FOREACH PACKAGE and prepare log capabilities
library(doParallel)
registerDoParallel(cores = 10)
library(foreach)

##
# Function to convert the grid search data frame to a rank aggregation format
# Input   gs: A grid search object for DCMM, each column represets a performance measure
##
.getRanksAndWeights <- function (gs, measures = c("ll","BIC","NRS","SLL","AGP","MPGP","CHPHT","CHPVT")) {
  
  #Error control
  stopifnot(all(measures %in% names(gs)))
  
  aux <- lapply(X = measures, FUN = function (meas) {
    theHighertheBetter <- meas %in% c("ll","SLL","AGP", "MPGP","CHPHT","CHPVT")
    
    sortedGs <- gs %>% arrange_(meas)
    if (theHighertheBetter) sortedGs <- sortedGs %>% slice(nrow(sortedGs):1)
    sortedGs <- sortedGs %>% unite(col = hyperParams, M, orderHC, orderVC, sep=",")
    #sortedGs$hyperParams <- paste("$\\mu(",sortedGs$hyperParams,")$")
    
    list(
      ranks = sortedGs$hyperParams,
      weights = sortedGs[[meas]]
    )
  }) %>% setNames(measures)
  
  result <- list(
    ranks = bind_rows(lapply(aux,"[[","ranks")) %>% t,
    weights = bind_rows(lapply(aux,"[[","weights")) %>% t
  )
  result$weights_normalized <- result$weights %>% apply(1, function (row) (row - min(row))/(max(row) - min(row))) %>% t

  result
}

#Create ordered rankings (predictability and interpretability)
ranks_predictability_measures <- .getRanksAndWeights(gs, measures = PREDICTABILITY_MEASURES)
ranks_interpretability_measures <- .getRanksAndWeights(gs, measures = INTERPRETABILITY_MEASURES)

##
# Predictability Rank Aggregation
# TODO: No se están haciendo las GA_RERUN_TIMES repeticiones del genetico, solo se hacen 2
##
predictability_raggr_pool <- llply(1:GA_RERUN_TIMES, function (i) {
  RankAggreg::RankAggreg(
    x = ranks_predictability_measures$ranks,
    k = K1,
    weights = ranks_predictability_measures$weights,
    method = RANK_AGGREGATION_METHOD,
    distance = "Spearman",
    importance=apply(ranks_predictability_measures$weights_normalized,1,var),
    maxIter = GA_MAX_ITER,
    convIn = GA_CONVERGENCE_ITERATIONS,
    popSize = GA_POPULATION_SIZE,
    CP = GA_CROSS_OVER_PROBABILITY,
    MP = GA_MUTATION_PROBABILITY,
    verbose = GA_VERBOSE
  )
},
.inform = FALSE,
.parallel = TRUE, 
.paropts = list(
  .packages = c("RankAggreg"),
  .export = c("ranks_predictability_measures")
))
predictability_raggr_rank <- predictability_raggr_pool[[which.min(lapply(predictability_raggr_pool, function (sol) sol$optimal.value))]]

##
# Interpretability Rank Aggregation
##
interpretability_raggr_pool <- llply(1:GA_RERUN_TIMES, function (i) {
  RankAggreg::RankAggreg(
    x = ranks_interpretability_measures$ranks,
    k = K1,
    weights = ranks_interpretability_measures$weights,
    method = RANK_AGGREGATION_METHOD,
    distance = "Spearman",
    importance=apply(ranks_interpretability_measures$weights_normalized,1,var),
    maxIter = GA_MAX_ITER,
    convIn = GA_CONVERGENCE_ITERATIONS,
    popSize = GA_POPULATION_SIZE,
    CP = GA_CROSS_OVER_PROBABILITY,
    MP = GA_MUTATION_PROBABILITY,
    verbose = GA_VERBOSE
  )
},
.inform = FALSE,
.parallel = TRUE, 
.paropts = list(
  .packages = c("RankAggreg"),
  .export = c("ranks_interpretability_measures")
))
interpretability_raggr_rank <- interpretability_raggr_pool[[which.min(lapply(interpretability_raggr_pool, function (sol) sol$optimal.value))]]


##
# Final Rank Aggregation
##
intersectionRanking <- rbind.fill(
  lapply(list(
    predictability_raggr_rank$top.list[predictability_raggr_rank$top.list %in% interpretability_raggr_rank$top.list],
    interpretability_raggr_rank$top.list[interpretability_raggr_rank$top.list %in% predictability_raggr_rank$top.list]
  ),
  function (x) {
    t(x) %>% as.data.frame
  })
)
rownames(intersectionRanking) <- c("Predictability", "Interpretability")

final_raggr_pool <- llply(.data = 1:GA_RERUN_TIMES,.fun = function (i) {
  RankAggreg::RankAggreg(
    x = as.matrix(intersectionRanking),
    k = K2,
    weights = NULL,
    method = RANK_AGGREGATION_METHOD,
    distance = "Spearman",
    importance= c(PREDICTABILITY_IMPORTANCE,INTERPRETABILITY_IMPORTANCE),
    maxIter = GA_MAX_ITER,
    convIn = GA_CONVERGENCE_ITERATIONS,
    popSize = GA_POPULATION_SIZE,
    CP = GA_CROSS_OVER_PROBABILITY,
    MP = GA_MUTATION_PROBABILITY,
    verbose = GA_VERBOSE
  )
},
.inform = TRUE,
.parallel = TRUE,
.paropts = list(
  .packages = c("RankAggreg"),
  .export = c("intersectionRanking"),
  .errorhandling = "stop"
))
raggr <- final_raggr_pool[[which.min(lapply(final_raggr_pool, function (sol) sol$optimal.value))]]

#Save output
if (SAVE_OUTPUT)
  saveRDS(gs, file = paste("output/raggr","(",Sys.Date(),")",".Rds"))


######################
## ANALYSIS
######################

#View & Plot
str(raggr$top.list)


#Plot the model! (Hidden & visible transitions)
bestHyperparams <- (raggr$top.list %>% strsplit(split = ","))[[1]] %>% 
  as.numeric %>% 
  setNames(c("M","l","f")) %>% 
  as.list
topModel <- gs %>% filter(M == bestHyperparams$M,orderHC == bestHyperparams$l, orderVC == bestHyperparams$f)
topModel <- topModel$fit[[1]]

# Step 1 - Visible transition matrix

# Step 2 - Hidden transition Matrix
if (VIEW_PLOT) {
  MmgraphR::march.Dcmm.trmatplot(
    d = topModel,
    seed = NULL,
    type = "hidden",
    hstate = 1,
    cspal = NULL,
    cpal = NULL,
    main = "Probability Transition Matrix",
    xlab = "Time",
    ylab = "States",
    ylim = NULL,
    xtlab = NULL,
    ytlab = NULL,
    pfilter = "tmax",
    shade.col = "white",
    num = 3,
    hide.col = "white",
    lorder = NULL,
    plot = TRUE,
    verbose = FALSE
  )
}

#End
if (exists("EXTERNAL_CALL")) rm(EXTERNAL_CALL)