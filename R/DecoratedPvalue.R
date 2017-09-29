# This function computes the decorated p-values for each path.
#
# @description This function computes the decorated p-values for each path.
#
# @usage getDecoratedPvalues(dataset, BestPaths, pathLength, nCases,
#                                     nControls, method, n_permutations = 100, strataF = NA)
#
# @param dataset Filename of the phenotype data. First column must correspond
#                to the gene symbols, and all other columns correspond to the patients.
#                First column must be named \code{GENE.Symbols} and must contain strings
#                corresponding to the symbols of the genes. All other columns must
#                contain 0-1 values depending if the patient has a variant in the
#                corresponding gene (a value of 1 if the variant exists and 0 otherwise).
# @param BestPaths A data frame returned from calling \link{ProcessPaths}.
# @param pathLength The maximum path length for which the top K paths
#                   will be computed. The default value is 5 (The maximum
#                   value for \emph{pathLength} is 5).
# @param nCases The number of cases (Must be greater than 1).
# @param nControls The number of controls (Must be greater than 1).
# @param method The method to be invoked. The value 1 calls the unsigned method
#               while the value 2 calls the signed method.
# @param threshold_percent All genes which occur in less than emph{threshold_percent}
#                          of patients are kept for processing.
# @param n_permutations The number of permutations to be used in computing the
#                   p-values. The default value is 100.
# @param strataF Filename for the strata file. \emph{strataF} must be a file
#                of two columns, first column named \emph{subjid }contains
#                the patient names similar to the dataset columns, and the
#                second column named \emph{stratum} must contain the integer
#                corresponding to the stratum.
#
# @return This function returns a list with the following items:
#         \item{SignedPaths}{The top K signed paths for each length.}
#         \item{Paths}{The top K paths for each length.}
#         \item{Subpaths1}{The first subpath.}d
#         \item{Subpaths1_cases}{The number of cases for Subpaths1.}
#         \item{Subpaths1_controls}{The number of controls for Subpaths1.}
#         \item{Subpaths2}{The second subpath.}
#         \item{Subpaths2_cases}{The number of cases for Subpaths2.}
#         \item{Subpaths2_controls}{The number of controls for Subpaths2.}
#         \item{Direction}{Direction of splitting along the path.}
#         \item{Lengths}{The lengths of each path.}
#         \item{Scores}{The scores of each path.}
#         \item{Pvalues}{The p-values of each path.}
#         \item{DecoratedPvalues}{The decorated p-values for the second subpath.}
#         \item{Total_Cases}{The total number of cases for each path.}
#         \item{Total_Controls}{The total number of controls for each path.}
#
getDecoratedPvalues <- function(dataset, BestPaths, pathLength, nCases, nControls,
                                method, threshold_percent, n_permutations, strataF,  ValueTable = NULL){

  # # Precompute values table
  # print("Computing Values Table...")
  if(is.null(ValueTable)){
    ValueTable <- getValuesTable(nCases, nControls)
  }

  # Preprocessing Phenotype dataset, all genes which occur in more
  # than threshold_percent of patients will be removed from the
  # dataset, also genes which do not occur in any patients will be removed
  dataset <- PreprocessTable(dataset, threshold_percent, nCases, nControls)
  genes <- dataset$genes
  data <- dataset$data
  genes_data <- data.frame(genes = genes, data, stringsAsFactors = F)
  patients <- dataset$patients
  rm(dataset)

  # Handling Strata File
  stratagroups <- NA
  stratanumbers <- NA
  if(length(strataF) == 1 & !is.na(strataF)){
    strata <- read.table(strataF, header = T,stringsAsFactors = F)
    if(names(strata)[1] != "subjid") stop("First column name in the strata file should be named subjid!")
    if(names(strata)[2] != "stratum") stop("Second column name in the strata file should be named stratum!")
    if(!all(strata[,1] %in% patients)) stop("All patient names in the strata file should be present in the column names of dataset!")
    if(!all(patients %in% strata[,1])) stop("All patient names in the dataset should be present in the strata file!")
    inds <- match(patients, strata[,1])
    strata <- strata[inds,] # strata patients now have the same order as the column names in the dataset
    stratanumbers <- as.integer(as.vector(strata[,2]))
    stratagroups <- as.vector(unique(stratanumbers))
  }

  print("Computing Decorated Pvalues...")

  DecoratedBestPaths <- data.frame(SignedPaths = character(0), Paths = character(0), Subpaths1 = character(0),
                                   Subpaths1_cases = integer(0), Subpaths_controls = integer(0),
                                   Subpaths2 = character(0), Subpaths1_cases = integer(0), Subpaths_controls = integer(0),
                                   Direction = integer(0), Lengths = integer(0), Scores = numeric(0),
                                   Pvalues = numeric(0), DecoratedPvalues = numeric(0),
                                   Total_Cases = integer(0), Total_Controls = integer(0))
  total_paths <- 1

  # Handle paths of length 1
  bestpaths <- BestPaths[BestPaths$Lengths == 1,]
  for(i in 1:nrow(bestpaths)){
    DecoratedBestPaths <- rbind(DecoratedBestPaths, rep(0,ncol(DecoratedBestPaths)))
    DecoratedBestPaths[total_paths,c(1,2,10,11,12,14,15)] = bestpaths[i,-ncol(bestpaths)]
    DecoratedBestPaths[total_paths,13] <- bestpaths$Pvalues[i]
    DecoratedBestPaths[total_paths,9] <- ""  # Direction
    DecoratedBestPaths[total_paths,8] <- 0
    DecoratedBestPaths[total_paths,7] <- 0
    DecoratedBestPaths[total_paths,6] <- ""
    DecoratedBestPaths[total_paths,5] <- bestpaths$Controls[i]
    DecoratedBestPaths[total_paths,4] <- bestpaths$Cases[i]
    DecoratedBestPaths[total_paths,3] <- bestpaths$Paths[i]
    total_paths <- total_paths + 1
  }

  # Handle paths of length 2 and greater
  for(path_len in 2:pathLength){

    bestpaths <- BestPaths[BestPaths$Lengths == path_len,]

    for(i in 1:nrow(bestpaths)){

      # Split a path which is a string of the form "Gene1 (-) -> Gene2 (+)" as an example
      split_path <- unlist(strsplit(bestpaths$SignedPaths[i], split = " -> "))

      # Extract the signs from the split path
      sign <- sapply(split_path, function(x) unlist(strsplit(x, split = " "))[2])
      sign <- ifelse(sign == "(+)", 1, -1)

      # Get the genes in the path
      genes <- sapply(split_path, function(x) unlist(strsplit(x, split = " "))[1])

      # Get the data for each gene in the path
      path_data_pos <- as.matrix(genes_data[match(genes, genes_data$genes),2:ncol(genes_data)],drop = FALSE)
      path_data_neg <- matrix(0, nrow = nrow(path_data_pos), ncol = ncol(path_data_pos))
      if(method == 2){
        path_data_neg[which(sign == -1),] = path_data_pos[which(sign == -1),]
        path_data_pos[which(sign == -1),] = 0
      }

      flipped <- ifelse(sign[1] == -1, T, F)
      names(flipped) <- NULL

      # Compute the p-values for the different splitted paths
      for(j in 1:(path_len-1)){
        DecoratedBestPaths <- rbind(DecoratedBestPaths, rep(0,ncol(DecoratedBestPaths)))

        DecoratedBestPaths[total_paths,c(1,2,10,11,12,14,15)] = bestpaths[i,-ncol(bestpaths)]
        subpath_pos1 <- colSums(path_data_pos[1:j,,drop=FALSE])
        subpath_pos1[subpath_pos1 != 0] <- 1
        subpath_neg1 <- colSums(path_data_neg[1:j,,drop=FALSE])
        subpath_neg1[subpath_neg1 != 0] <- 1
        subpath_pos2 <- path_data_pos[j+1,]
        subpath_neg2 <- path_data_neg[j+1,]
        lst <- computeDecoratedPvalue(subpath_pos1, subpath_neg1, subpath_pos2, subpath_neg2,
                                                                 path_len, nCases, nControls, method, n_permutations, ValueTable, stratagroups, stratanumbers, strataF, flipped)
        DecoratedBestPaths[total_paths,13] <- lst$decorated_pvalue
        DecoratedBestPaths[total_paths,9] <- "Forward"  # Direction
        DecoratedBestPaths[total_paths,8] <- lst$controls2
        DecoratedBestPaths[total_paths,7] <- lst$cases2
        DecoratedBestPaths[total_paths,6] <- paste(genes[j+1])
        DecoratedBestPaths[total_paths,5] <- lst$controls1
        DecoratedBestPaths[total_paths,4] <- lst$cases1
        DecoratedBestPaths[total_paths,3] <- paste(genes[1:j], collapse = " -> ")
        total_paths <- total_paths + 1
      }

      # Compute the p-values for the different splitted paths going in the opposite directions
      for(j in path_len:2){
        DecoratedBestPaths <- rbind(DecoratedBestPaths, rep(0,ncol(DecoratedBestPaths)))
        DecoratedBestPaths[total_paths,c(1,2,10,11,12,14,15)] = bestpaths[i,-ncol(bestpaths)]
        subpath_pos1 <- colSums(path_data_pos[path_len:j,,drop=FALSE])
        subpath_pos1[subpath_pos1 != 0] <- 1
        subpath_neg1 <- colSums(path_data_neg[path_len:j,,drop=FALSE])
        subpath_neg1[subpath_neg1 != 0] <- 1
        subpath_pos2 <- path_data_pos[j-1,]
        subpath_neg2 <- path_data_neg[j-1,]
        lst <- computeDecoratedPvalue(subpath_pos1, subpath_neg1, subpath_pos2, subpath_neg2,
                                                                 path_len, nCases, nControls, method, n_permutations, ValueTable, stratagroups, stratanumbers, strataF, flipped)
        DecoratedBestPaths[total_paths,13] <- lst$decorated_pvalue
        DecoratedBestPaths[total_paths,9] <- "Backward"  # Direction
        DecoratedBestPaths[total_paths,8] <- lst$controls2
        DecoratedBestPaths[total_paths,7] <- lst$cases2
        DecoratedBestPaths[total_paths,6] <- paste(genes[j-1])
        DecoratedBestPaths[total_paths,5] <- lst$controls1
        DecoratedBestPaths[total_paths,4] <- lst$cases1
        DecoratedBestPaths[total_paths,3] <- paste(genes[path_len:j], collapse = " -> ")
        total_paths <- total_paths + 1
      }

    }

  }

  names(DecoratedBestPaths) = c("SignedPaths", "Paths", "Subpaths1", "Subpaths1_Cases", "Subpaths1_Controls", "Subpaths2",
                                "Subpaths2_Cases", "Subpaths2_Controls", "Direction", "Lengths", "Scores",
                                "Pvalues", "DecoratedPvalues", "Cases", "Controls")
  return(DecoratedBestPaths)

}



# This function computes the corrected p-value of two given sub-paths
computeDecoratedPvalue <- function(subpath_pos1, subpath_neg1, subpath_pos2, subpath_neg2, path_len,
                                   nCases, nControls, method, n_permutations, ValueTable, stratagroups, stratanumbers, strataF, flipped){

  inds_pos1 <- which(subpath_pos1 != 0)
  inds_neg1 <- which(subpath_neg1 != 0)
  inds_pos2 <- which(subpath_pos2 != 0)
  inds_neg2 <- which(subpath_neg2 != 0)

  subpath_pos2[intersect(inds_pos1, inds_pos2)] <- 0
  subpath_neg2[intersect(inds_neg1, inds_neg2)] <- 0
  inds_pos2 <- setdiff(inds_pos2, inds_pos1)
  inds_neg2 <- setdiff(inds_neg2, inds_neg1)

  # compute the score of the path
  # cases1 <- sum(subpath_pos1[1:nCases]) + sum(subpath_neg1[(nCases+1):(nCases+nControls)])
  # controls1 <- sum(subpath_pos1[(nCases+1):(nCases+nControls)]) + sum(subpath_neg1[1:nCases])
  # cases2 <- sum(subpath_pos2[1:nCases]) + sum(subpath_neg2[(nCases+1):(nCases+nControls)])
  # controls2 <- sum(subpath_pos2[(nCases+1):(nCases+nControls)]) + sum(subpath_neg2[1:nCases])
  # score <- ifelse(flipped & method == 1, ValueTable[(controls1+controls2+1),(cases1+cases2+1)], ValueTable[(cases1+cases2+1),(controls1+controls2+1)])
  case_pos1 <- sum(subpath_pos1[1:nCases])
  case_neg1 <- sum(subpath_neg1[(nCases+1):(nCases+nControls)])
  control_pos1 <- sum(subpath_pos1[(nCases+1):(nCases+nControls)])
  control_neg1 <- sum(subpath_neg1[1:nCases])
  case_pos2 <- sum(subpath_pos2[1:nCases])
  case_neg2 <- sum(subpath_neg2[(nCases+1):(nCases+nControls)])
  control_pos2 <-sum(subpath_pos2[(nCases+1):(nCases+nControls)])
  control_neg2 <- sum(subpath_neg2[1:nCases])
  if(method == 1){
    score <- ValueTable[(case_pos1+case_pos2+case_neg1+case_neg2+1),(control_pos1+control_pos2+control_neg1+control_neg2+1)]
  }else if(method == 2){
    score <- ValueTable[(case_pos1+case_pos2+1),(control_pos1+control_pos2+1)] + ValueTable[(case_neg1+case_neg2+1),(control_neg1+control_neg2+1)]
  }
  scores <- c()

  toSample_pos <- setdiff(1:length(subpath_pos1), inds_pos1)
  toSample_neg <- setdiff(1:length(subpath_neg1), inds_neg1)

  ncol_data <- nCases + nControls
  stratagroup_numbers <- list()
  strata_inds_pos <- list()
  strata_inds_neg <- list()
  if(length(strataF) == 1 & !is.na(strataF)){
    inds_1 <- unique(c(inds_pos1, inds_neg1))
    for(i in stratagroups){
      stratagroup_numbers[[toString(i)]] <- setdiff(which(stratanumbers == i), inds_1)
      strata_inds_pos[[toString(i)]] <- inds_pos2[which(inds_pos2 %in% stratagroup_numbers[[toString(i)]])]
      strata_inds_neg[[toString(i)]] <- inds_neg2[which(inds_neg2 %in% stratagroup_numbers[[toString(i)]])]
    }
  }

  # Process the permutations
  for(i in 1:n_permutations){

    sample_inds_pos <- NA
    sample_inds_neg <- NA
    stratagroup_numbers2 <- stratagroup_numbers
    strata_inds_pos2 <- strata_inds_pos
    strata_inds_neg2 <- strata_inds_neg
    if(length(strataF) == 1 & !is.na(strataF)){
      sample_inds_pos <- c()
      sample_inds_neg <- c()
      for(i in stratagroups){
        stratagroup <- stratagroup_numbers2[[toString(i)]]
        strata_group_pos <- strata_inds_pos2[[toString(i)]]
        sample_group <- sample(stratagroup, length(strata_group_pos), replace = F)
        sample_inds_pos <- c(sample_inds_pos, sample_group)
        stratagroup <- setdiff(stratagroup, sample_group)
        strata_group_neg <- strata_inds_neg2[[toString(i)]]
        sample_inds_neg <- c(sample_inds_neg, sample(stratagroup, length(strata_group_neg), replace = F))
      }
    }else{
      sample_inds_pos <- sample(toSample_pos, length(inds_pos2), replace = F)
      sample_inds_neg <- sample(toSample_neg, length(inds_neg2), replace = F)
    }
    # cases_pos <- length(which(sample_inds_pos <= nCases))
    # cases_neg <- length(which(sample_inds_neg > nCases))
    # perm_cases <- cases_pos + cases_neg
    # perm_controls <- (length(sample_inds_pos) - cases_pos) + (length(sample_inds_neg) - cases_neg)
    # perm_score <- ifelse(flipped & method == 1, ValueTable[(controls1+perm_controls+1),(cases1+perm_cases+1)], ValueTable[(cases1+perm_cases+1),(controls1+perm_controls+1)])
    perm_case_pos <- length(which(sample_inds_pos <= nCases))
    perm_case_neg <- length(which(sample_inds_neg > nCases))
    perm_control_pos <- length(sample_inds_pos) - perm_case_pos
    perm_control_neg <- length(sample_inds_neg) - perm_case_neg
    if(method == 1){
      perm_score <- ValueTable[(case_pos1+perm_case_pos+case_neg1+perm_case_neg+1),(control_pos1+perm_control_pos+control_neg1+perm_control_neg+1)]
    }else if(method == 2){
      perm_score <- ValueTable[(case_pos1+perm_case_pos+1),(control_pos1+perm_control_pos+1)] + ValueTable[(case_neg1+perm_case_neg+1),(control_neg1+perm_control_neg+1)]
    }
    scores <- c(scores, perm_score)

  }

  pvalue <- length(which(scores >= score))/length(scores)

  lst <- list()
  # lst[["cases1"]] <- cases1
  # lst[["cases2"]] <- cases2
  # lst[["controls1"]] <- controls1
  # lst[["controls2"]] <- controls2
  lst[["cases1"]] <- case_pos1 + case_neg1
  lst[["cases2"]] <- case_pos2 + case_neg2
  lst[["controls1"]] <- control_pos1 + control_neg1
  lst[["controls2"]] <- control_pos2 + control_neg2
  lst[["decorated_pvalue"]] <- pvalue
  return(lst)

}
