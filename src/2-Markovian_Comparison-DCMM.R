# 
# This script trains and evaluates DCMMs using the march package. The configuration of the DCMM (M,l,f) is explored
# in a grid search fashion.
#
# Requirements:
#   - The objects "marchData", "marchData.training" and "marchData.test" must exist, as a result of executing the
#     script "1-Markovian_Comparison-InputData".
#
# Parameters:
#     SAVE_OUTPUT: T/F
#     MEMETIC_GENERATIONS: Integer. Number of generations (memetic algorithm)
#     MEMETIC_POP_SIZE: Integer. Population size (memetic algorithm)
#     ITER_BW: Integer. Maximum number of iterations in the Baum-Welch algorithm.
#     STOP_BW: Real. Tolerance in the Baum-Welch algorithm.
#     PENALTY_RATE. Real. Penalty rate for the measure AGP (Accuracy of Generated Predictions)
#     TRANSITION_RELEVANCY_THRESHOLD. Real. Used in the measure CHPHT (Coefficient of High Probability Hidden Transitions)
#     M_RANGE. Integer vector. Range of M (Number of states) in the grid search 
#     f_RANGE. Integer vector. Range of f (order of the visible chain) in the grid search
#     l_RANGE. Integer vector. Range of l (order of the hidden chain) in the grid search
#
# Output:
#     object gs: Dataframe containing the grid of DCMMs learnt along with the evaluation measures computed for all of them
#
library(ProjectTemplate)
load.project()
SAVE_OUTPUT <- FALSE


#Script conditions
stopifnot(exists("marchData"))
stopifnot(exists("marchData.training"))
stopifnot(exists("marchData.test"))

#Parameters - learning
MEMETIC_GENERATIONS <- 5
MEMETIC_POP_SIZE <- 5
ITER_BW <- 4
STOP_BW <- 1e-02

#Parameter - Evaluation
STATE_VISIT_FREQUENCY_THRESHOLD <- 0.05 #NRS
PENALTY_RATE <- 2 # AGP
TRANSITION_RELEVANCY_THRESHOLD <- 0.9 #CHPHT

M_RANGE <- 3:4
f_RANGE <- 0:2
l_RANGE <- 1:3

# Define a named list of parameter values
gs <- list(
  orderHC = f_RANGE,
  orderVC = l_RANGE,
  M = M_RANGE
) %>% purrr::cross_d() %>% 
  filter(orderHC <= M)

# LOAD FOREACH PACKAGE and prepare log capabilities
library(doParallel)
registerDoParallel(cores = 10)
library(foreach)

# Parallel Grid Search
print("Computing DCMM models...")
gs$fit <- foreach(i = 1:nrow(gs),
               .multicombine = TRUE,
               .inorder = TRUE,
               .packages = c("march"),
               .errorhandling = "stop",
               .export = c("gs","marchData.training"),
               .verbose = FALSE) %dopar% {
  #Log("Processing block %d of %d [DCMM(%d,%d,%d)]", i, nrow(gs),gs[[i,"M"]],gs[[i,"orderHC"]],gs[[i,"orderVC"]])
  set.seed(1234)
  march.dcmm.construct(
    y = marchData.training,
    orderHC = gs[[i,"orderHC"]],
    orderVC = gs[[i,"orderVC"]],
    M = gs[[i,"M"]],
    gen = MEMETIC_GENERATIONS,
    popSize = MEMETIC_POP_SIZE,
    maxOrder = max(gs$orderVC),
    seedModel = NULL,
    iterBw = ITER_BW,
    stopBw = 1e-03
  )
}

# Compute metrics for the model (algunas se computan en paralelo para ser mas eficiente)
print("Computing evaluation measures...Interpretability")
gs <- gs %>% mutate(
  ll = map_dbl(fit, ~ .@ll),
  param = map_int(fit, ~ march.summary(.)[2] %>% as.integer),
  BIC = map_dbl(fit,march.BIC),
  NRS = map_int(fit, ~ march.helper.NRS(d=.,State_Visit_Frequency_Threshold = STATE_VISIT_FREQUENCY_THRESHOLD)),
  CHPHT = map_dbl(fit,~march.helper.CHPHT(.,transition_relevancy_threshold = TRANSITION_RELEVANCY_THRESHOLD)),
  CHPVT = map_dbl(fit,~march.helper.CHPVT(.,transition_relevancy_threshold = TRANSITION_RELEVANCY_THRESHOLD))
  )

print("Computing evaluation measures...Predictability")
#Sequence (Log)-Likelihood
gs$SLL <- foreach(i = 1:nrow(gs), .combine = "c", .packages = c("march"), .errorhandling = "stop", .export = c("gs","marchData.test")) %dopar%  {
  mean(march.helper.apply_dataset(y=marchData.test, f = function (s) march.helper.sequenceLogLikelihood(gs$fit[[i]],s)))
}

#Accuracy of Generated Predictions
gs$AGP <- foreach(i = 1:nrow(gs), .combine = "c", .packages = c("march"), .errorhandling = "stop", .export = c("gs","marchData.test")) %dopar%  {
  mean(march.helper.apply_dataset(y=marchData.test, f = function (s) march.helper.MMAS(gs$fit[[i]],s,PENALTY_RATE)))
}

#Min precision
gs$MPGP <- foreach(i = 1:nrow(gs), .combine = "c", .packages = c("march"), .errorhandling = "stop", .export = c("gs","marchData.test")) %dopar%  {
  mean(march.helper.apply_dataset(y=marchData.test, f = function (s) march.helper.minPrecision(gs$fit[[i]],s)))
}

# Coefficient of High Probability Hidden Transitions
gs$CHPHT <- foreach(i = 1:nrow(gs), .combine = "c", .packages = c("march"), .errorhandling = "stop", .export = c("gs","marchData.test")) %dopar%  {
  march.helper.CHPHT(gs$fit[[i]],transition_relevancy_threshold = TRANSITION_RELEVANCY_THRESHOLD)
}

gs$CHPVT <- foreach(i = 1:nrow(gs), .combine = "c", .packages = c("march"), .errorhandling = "stop", .export = c("gs","marchData.test")) %dopar%  {
  march.helper.CHPVT(gs$fit[[i]],transition_relevancy_threshold = TRANSITION_RELEVANCY_THRESHOLD)
}

#Save data
if (SAVE_OUTPUT)
{
  save(gs, file = paste0("output/",
                         "M=[",min(M_RANGE),",",max(M_RANGE),"]",
                         "f=[",min(f_RANGE),",",max(f_RANGE),"]",
                         "l=[",min(l_RANGE),",",max(l_RANGE),"]",
                         "(",Sys.Date(),")",".RData"))
  
}

#View & Plot
View(gs)
