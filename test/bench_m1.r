#options(echo = TRUE)

t <- 0
i <- 100
l <- 3
k <- 4

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0)
 file <- args[1]

if(length(args) > 1)
 t <- as.numeric(args[2])

if(length(args) > 2)
 i <- as.numeric(args[3])

if(length(args) > 3)
 l <- as.numeric(args[4])

if(length(args) > 4)
 k <- as.numeric(args[5])

if(length(args) > 5)
 c <- as.numeric(args[6])

library(geneticsCRE)

res <- GetBestPaths(file, nCases = c, nControls = c, method = 'enrich', threshold_percent = 1, K = k, pathLength = l, iterations = i, strataF = NA, nthreads = t)
res[,c("Paths", "Scores", "Pvalues")]

