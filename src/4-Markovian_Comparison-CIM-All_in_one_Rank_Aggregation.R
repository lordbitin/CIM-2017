# 
# This script performs several rank aggregation process of a set of trained DCMMs. Each of the rank aggregation processes
# is executed for a different balance of the predictability/interpretability importance
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
#     Q: Size of the aggregated ranking
#     RANK_AGGREGATION_METHOD: Method used to performed rank aggregation (See RankAggreg package)
#     GA_MAX_ITER: Maximum number of iterations of the GA (See RankAggreg package)
#     GA_CONVERGENCE_ITERATIONS: Converege iterations of the GA (See RankAggreg package)
#     GA_POPULATION_SIZE: Population size in the GA (See RankAggreg package)
#     GA_CROSS_OVER_PROBABILITY: (See RankAggreg package)
#     GA_MUTATION_PROBABILITY: (See RankAggreg package)
#     GA_RERUN_TIMES: Number of times that the GA will be run, in order to choose the best solution.
#     GA_VERBOSE: T/F
#
# Output:
#     object m, containing a table of DCMM rankings for each evaluation measure
#     Heatmaps showing the best position of a DCMM in the top-10 ranking in terms of the three hyperparameters 
#     of a DCMM: $M$, $l$ and $f$. The $x$-axis shows how the results vary when we modify the importance of the 
#     predictability of the model, and thus, modifying importance of the interpretability too. 
#     White color is used when no model with such hyperparameters is found in the top-10 ranking. 
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

#Dependencies
base::source("./src/Markovian_Comparison_helpers.R")
#extrafont::font_import(pattern="[A/a]rial")

#Constraints
stopifnot(exists("gs"))

#Parameters
if (SET_PARAMETERS) {
  #PREDICTABILITY_MEASURES <- c("BIC","SLL","AGP", "MPGP")
  PREDICTABILITY_MEASURES <- c("BIC", "SLL","AGP", "MPGP")
  #INTERPRETABILITY_MEASURES <- c("BIC","NRS",,"CHPVT")
  INTERPRETABILITY_MEASURES <- c("BIC","CHPHT")
  Q <- 10
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

###################################
## SCRIPT CODE
###################################

# LOAD FOREACH PACKAGE FOR PARALLELIZATION
library(doParallel)
registerDoParallel(cores = 10)
library(foreach)

#Create weighted rankings for all the measures (predictability + interpretability)
ranks_all_measures <- Markovian_Comparison_helpers.getRanksAndWeights(
  gs, 
  measures = c(PREDICTABILITY_MEASURES,INTERPRETABILITY_MEASURES) %>% unique
  )

.getImportanceVector <- function (P_Importance, dispersion_measure) {
  importance_vector <- sapply(rownames(ranks_all_measures$weights_normalized),
    function (measure_name) {
      if (measure_name %in% PREDICTABILITY_MEASURES & 
          measure_name %in% INTERPRETABILITY_MEASURES) {
        ifelse(
          PREDICTABILITY_IMPORTANCE <= INTERPRETABILITY_IMPORTANCE,
          dispersion_measure(ranks_all_measures$weights_normalized[measure_name,])*P_Importance,
          dispersion_measure(ranks_all_measures$weights_normalized[measure_name,])*(1-P_Importance)
        )
      }
      else if (measure_name %in% PREDICTABILITY_MEASURES) {
        dispersion_measure(ranks_all_measures$weights_normalized[measure_name,])*P_Importance
      }
      else {
        dispersion_measure(ranks_all_measures$weights_normalized[measure_name,])*(1-P_Importance)
      }
    })
  
  importance_vector
}

evolution_pred_df <- tibble(P_Importance = seq(0,1,0.05))
evolution_pred_df$raggr <- llply(
  .data = 1:nrow(evolution_pred_df),
  .fun = function (i) {
    P_Importance <- evolution_pred_df[[i,"P_Importance"]]
    all_raggr_pool <- llply(1:GA_RERUN_TIMES, function (i) {
      RankAggreg::RankAggreg(
        x = ranks_all_measures$ranks,
        k = Q,
        weights = ranks_all_measures$weights,
        method = RANK_AGGREGATION_METHOD,
        distance = "Spearman",
        importance=.getImportanceVector(P_Importance = P_Importance, 
                                        dispersion_measure = Markovian_Comparison_helpers.coefficient_of_variation),
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
      .export = c("ranks_all_measures")
    ))
    all_raggr_rank <- all_raggr_pool[[which.min(lapply(all_raggr_pool, function (sol) sol$optimal.value))]]
    all_raggr_rank
  }
)

rank_position_df <- evolution_pred_df %>% mutate(
  topList = purrr::map(raggr,function (x) x$top.list)
) %>% 
  unnest(topList) %>% 
  separate(col = topList,into = c("M","l","f"),sep = ",", convert = T)
rank_position_df$rank_position <- 1:Q

#Plot rank position df
if (VIEW_PLOT) {
  y_axis <- base::expression(paste("Number of states (",italic("M"),")"))
  plotM <- ggplot(rank_position_df %>% group_by(P_Importance) %>% distinct(M,.keep_all = T)) + 
    geom_tile(mapping = aes(x=P_Importance,y=factor(M),fill=rank_position)) + 
    scale_fill_continuous(name = "Best\nRank\nPosition", breaks = c(1,5,10), limits = c(1,10)) +
    labs(x = "Predictability Importance", y = y_axis) + 
    theme(axis.title = element_text(size = 8), 
    legend.text = element_text(size = 8), 
    legend.title = element_text(size = 8), 
    panel.background = element_rect(fill = NA)) + 
    theme(axis.text = element_text(size = 8), 
    axis.text.y = element_text(vjust = 0.35))
  
  y_axis <- base::expression(paste("Hidden chain order (",italic("l"),")"))
  plotL <- ggplot(rank_position_df %>% group_by(P_Importance) %>% distinct(l,.keep_all = T)) + 
    geom_tile(mapping = aes(x=P_Importance,y=factor(l),fill=rank_position)) + 
    scale_fill_continuous(name = "Best\nRank\nPosition", breaks = c(1,5,10), limits = c(1,10)) +
    labs(x = "Predictability Importance", y = y_axis)+ 
    theme(axis.title = element_text(size = 8), 
          axis.text.y = element_text(vjust = 0), 
          legend.text = element_text(size = 8), 
          legend.title = element_text(size = 8), 
          panel.background = element_rect(fill = NA)) + 
    theme(axis.text = element_text(size = 8), 
          axis.text.y = element_text(vjust = 0.35))
  
  y_axis <- base::expression(paste("Visible chain order (",italic("f"),")"))
  plotF <-ggplot(rank_position_df %>% group_by(P_Importance) %>% distinct(f,.keep_all = T)) + 
    geom_tile(mapping = aes(x=P_Importance,y=factor(f),fill=rank_position)) + 
    scale_fill_continuous(name = "Best\nRank\nPosition", breaks = c(1,5,10), limits = c(1,10)) +
    labs(x = "Predictability Importance", y = y_axis)+ 
    theme(axis.title = element_text(size = 8), 
          axis.text.y = element_text(vjust = 0), 
          legend.text = element_text(size = 8), 
          legend.title = element_text(size = 8), 
          panel.background = element_rect(fill = NA)) + 
    theme(axis.text = element_text(size = 8), 
          axis.text.y = element_text(vjust = 0.35))
  
  legend <- cowplot::get_legend(plotM)
  cowplot::plot_grid(
    plotM + theme(legend.position = "none"),
    plotL + theme(legend.position = "none"),
    plotF + theme(legend.position = "none"),
    legend,
    rel_widths = c(1,1,1,0.3),
    nrow  =1)
}

#To print in latex
m <- paste0("$\\mu(",ranks_all_measures$ranks,")",
            "[", 
            ranks_all_measures$weights %>% round(digits=3),
            "]$")
m <- matrix(m, nrow = nrow(ranks_all_measures$ranks),
            ncol=ncol(ranks_all_measures$ranks))
rownames(m) <- rownames(ranks_all_measures$ranks)
m <- t(m)
m <- m[,c("SLL","AGP","MPGP","BIC","CHPHT")]

if (SAVE_OUTPUT)
{
  write.csv(x = m[1:5,], file="output/Markovian_Comparison/measureRankings.csv",row.names = TRUE)
}


#Output
if (SAVE_OUTPUT)
{
  save(rank_position_df, 
       file = paste0("output/Markovian_Comparison/rank_position_df",
                     "(",Sys.Date(),")",".RData"))
}
