library(devtools)
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
load_all()

res2 <- GetBestPaths('/home/bichkd/workspace/gcre-ci/data/test/segfault-random.gz', nCases = 1537, nControls = 1537, method = 1,
  threshold_percent = 0.05, K = 4, pathLength = 4, iterations = 10, strataF = NA, nthreads = 0)

res2[,c("Paths", "Scores", "Pvalues")]
