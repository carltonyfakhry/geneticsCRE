#library(devtools)
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
devtools::load_all()

#data <- '/home/bichkd/workspace/gcre-ci/data/test/testdata.gz'
data <- '/home/bichkd/workspace/gcre-ci/data/sim/SimData_LoF_100_LDLR-CDK19-APOE-LPL-LRP1_4.txt.gz'
num_cc <- 100

res2 <- GetBestPaths(data, nCases = num_cc, nControls = num_cc, method = 1,
  threshold_percent = 0.05, K = 8, pathLength = 4, iterations = 100, strataF = NA, nthreads = 4)

res2[,c("Paths", "Scores", "Pvalues")]

newpvals <- getDecoratedPvalues(dataset = data, BestPaths = res2,
  pathLength = 4, nCases = num_cc, nControls = num_cc,
  method = 1, n_permutations = 100, strataF = NA)

newpvals
