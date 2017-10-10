load("output/Markovian_Comparison/marchData_NO_TSR_MIN_INTERACTIONS_15(2017-01-21).RData")
load("output/Markovian_Comparison/marchData_training_NO_TSR_MIN_INTERACTIONS_15(2017-01-21).RData")
load("output/Markovian_Comparison/marchData_test_NO_TSR_MIN_INTERACTIONS_15(2017-01-21).RData")

#Schedule model creation
start.time <- Sys.time()

M_RANGE_total <- 3:8
f_RANGE_total <- 1:3
l_RANGE_total <- 0:2

for (f in f_RANGE_total) {
  for (l in l_RANGE_total) {
    print(paste(
      "Starting DCMM learning with hyperparameters:",
      "M=[",min(M_RANGE),",",max(M_RANGE),"]",
      "f=[",min(f_RANGE),",",max(f_RANGE),"]",
      "l=[",min(l_RANGE),",",max(l_RANGE),"]"))
    
    M_RANGE <- M_RANGE_total
    f_RANGE <- f:f
    l_RANGE <- l:l
    
    source("src/Markovian_Comparison/2-Markovian_Comparison-DCMM.R")
    
    print(paste(
      "Finished DCMM learning with hyperparameters:",
      "M=[",min(M_RANGE),",",max(M_RANGE),"]",
      "f=[",min(f_RANGE),",",max(f_RANGE),"]",
      "l=[",min(l_RANGE),",",max(l_RANGE),"]"))
  }
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken