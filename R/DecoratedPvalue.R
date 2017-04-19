#' This function computes the decorated p-values for each path.
#'
#' @description This function computes the decorated p-values for each path.
#'
#' @usage getDecoratedPvalues(dataset, BestPaths, pathLength, nCases,
#'                                     nControls, method, n_permutations = 100, strataF = NA)
#'
#' @param dataset Filename of the phenotype data. First column must correspond
#'                to the gene symbols; all other columns correspond to the patients.
#'                First column must be named symbols and must contrain strings corresponding
#'                to the symbols of the genes. All other columns either contain a 0 if the
#'                gene is not present in the patient, and 1 or 2 if the gene is present
#'                in the patient.
#' @param BestPaths A data frame returned from calling \link{ProcessPaths}.
#' @param pathLength The maximum path length for which the top K paths
#'                   will be computed. The default value is 5 (The maximum
#'                   value for \emph{pathLength} is 5).
#' @param nCases The number of cases (Must be greater than 1).
#' @param nControls The number of controls (Must be greater than 1).
#' @param method The method to be invoked. The value 1 calls the unsigned method
#'               while the value 2 calls the signed method.
#' @param n_permutations The number of permutations to be used in computing the
#'                   p-values. The default value is 100.
#' @param strataF Filename for the strata file. \emph{strataF} must be a file
#'                of two columns, first column named \emph{subjid }contains
#'                the patient names similar to the dataset columns, and the
#'                second column named \emph{stratum} must contain the integer
#'                corresponding to the stratum.
#'
#' @return This function returns a list with the following items:
#'         \item{SignedPaths}{The top K signed paths for each length.}
#'         \item{Paths}{The top K paths for each length.}
#'         \item{Subpaths1}{The first subpath.}
#'         \item{Subpaths2}{The second subpath.}
#'         \item{Direction}{Direction of splitting along the path.}
#'         \item{Lengths}{The lengths of each path.}
#'         \item{Scores}{The scores of each path.}
#'         \item{Pvalues}{The p-values of each path.}
#'         \item{DecoratedPvalues}{The decorated p-values for the second subpath.}
#'
#' @author Carl Tony Fakhry
#'
#' @references Franceschini, A (2013). STRING v9.1: protein-protein interaction networks, with increased coverage
#'             and integration. In:'Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094.
#'             Epub 2012 Nov 29'.
#'
#'@export
getDecoratedPvalues <- function(dataset, BestPaths, pathLength, nCases, nControls, method, n_permutations = 100, strataF = NA){

  # Precompute values table
  print("Computing Values Table...")
  ValueTable <- geneticsCRE:::getValuesTable(nCases, nControls)

  # Preprocessing Phenotype dataset, all genes which occur in more
  # than threshold_percent of patients will be removed from the
  # dataset, also genes which do not occur in any patients will be removed
  print("Preprocessing Phenotype dataset...")
  threshold_percent <- 1
  dataset <- geneticsCRE:::PreprocessTable(dataset, threshold_percent, nCases, nControls)
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
                                   Subpaths2 = character(0), Direction = integer(0), Lengths = integer(0), Scores = numeric(0),
                                   Pvalues = numeric(0), DecoratedPvalues = numeric(0))
  total_paths <- 1

  # Handle paths of length 1
  bestpaths <- BestPaths[BestPaths$Lengths == 1,]
  for(i in 1:nrow(bestpaths)){
    DecoratedBestPaths <- rbind(DecoratedBestPaths, rep(0,ncol(DecoratedBestPaths)))
    DecoratedBestPaths[total_paths,c(1,2,6,7,8)] = bestpaths[i,]
    DecoratedBestPaths[total_paths,9] <- bestpaths$Pvalues[i]
    DecoratedBestPaths[total_paths,5] <- ""
    DecoratedBestPaths[total_paths,3] <- bestpaths$Paths[i]
    DecoratedBestPaths[total_paths,4] <- ""
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
        DecoratedBestPaths[total_paths,c(1,2,6,7,8)] = bestpaths[i,]
        subpath_pos1 <- colSums(path_data_pos[1:j,,drop=FALSE])
        subpath_pos1[subpath_pos1 != 0] <- 1
        subpath_neg1 <- colSums(path_data_neg[1:j,,drop=FALSE])
        subpath_neg1[subpath_neg1 != 0] <- 1
        subpath_pos2 <- path_data_pos[j+1,]
        subpath_neg2 <- path_data_neg[j+1,]
        decorated_pvalue <- geneticsCRE:::computeDecoratedPvalue(subpath_pos1, subpath_neg1, subpath_pos2, subpath_neg2,
                                                                 path_len, nCases, nControls, method, n_permutations, ValueTable, stratagroups, stratanumbers, strataF, flipped)
        DecoratedBestPaths[total_paths,9] <- decorated_pvalue
        DecoratedBestPaths[total_paths,5] <- "Forward"
        DecoratedBestPaths[total_paths,3] <- paste(genes[1:j], collapse = " -> ")
        DecoratedBestPaths[total_paths,4] <- paste(genes[j+1])
        total_paths <- total_paths + 1
      }

      # Compute the p-values for the different splitted paths going in the opposite directions
      for(j in path_len:2){
        DecoratedBestPaths <- rbind(DecoratedBestPaths, rep(0,ncol(DecoratedBestPaths)))
        DecoratedBestPaths[total_paths,c(1,2,6,7,8)] = bestpaths[i,]
        subpath_pos1 <- colSums(path_data_pos[path_len:j,,drop=FALSE])
        subpath_pos1[subpath_pos1 != 0] <- 1
        subpath_neg1 <- colSums(path_data_neg[path_len:j,,drop=FALSE])
        subpath_neg1[subpath_neg1 != 0] <- 1
        subpath_pos2 <- path_data_pos[j-1,]
        subpath_neg2 <- path_data_neg[j-1,]
        decorated_pvalue <- geneticsCRE:::computeDecoratedPvalue(subpath_pos1, subpath_neg1, subpath_pos2, subpath_neg2,
                                                                 path_len, nCases, nControls, method, n_permutations, ValueTable, stratagroups, stratanumbers, strataF, flipped)
        DecoratedBestPaths[total_paths,9] <- decorated_pvalue
        DecoratedBestPaths[total_paths,5] <- "Backward"
        DecoratedBestPaths[total_paths,3] <- paste(genes[path_len:j], collapse = " -> ")
        DecoratedBestPaths[total_paths,4] <- paste(genes[j-1])
        total_paths <- total_paths + 1
      }

    }

  }

  names(DecoratedBestPaths) = c("SignedPaths", "Paths", "Subpaths1", "Subpaths2", "Direction", "Lengths", "Scores",
                                "Pvalues", "DecoratedPvalues")
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
  cases1 <- sum(subpath_pos1[1:nCases]) + sum(subpath_neg1[nCases:(nCases+nControls)])
  controls1 <- sum(subpath_pos1[nCases:(nCases+nControls)]) + sum(subpath_neg1[1:nCases])
  cases2 <- sum(subpath_pos2[1:nCases]) + sum(subpath_neg2[nCases:(nCases+nControls)])
  controls2 <- sum(subpath_pos2[nCases:(nCases+nControls)]) + sum(subpath_neg2[1:nCases])
  score <- ifelse(flipped & method == 1, ValueTable[(controls1+controls2+1),(cases1+cases2+1)], ValueTable[(cases1+cases2+1),(controls1+controls2+1)])
  scores <- c()

  toSample_pos <- setdiff(1:length(subpath_pos1), inds_pos1)
  toSample_neg <- setdiff(1:length(subpath_neg1), inds_neg1)

  ncol_data <- nCases + nControls
  stratagroup_numbers <- list()
  if(length(strataF) == 1 & !is.na(strataF)){
    inds_1 <- unique(c(inds_pos1, inds_neg1))
    for(i in stratagroups){
      stratagroup_numbers[[toString(i)]] <- setdiff(which(stratanumbers == i), inds_1)
    }
  }

  # Process the permutations
  for(i in 1:n_permutations){

    sample_inds_pos <- NA
    sample_inds_neg <- NA
    stratagroup_numbers2 <- stratagroup_numbers
    if(length(strataF) == 1 & !is.na(strataF)){
      sample_inds_pos <- unlist(lapply(stratagroup_numbers2, function(x){sample(x, length(inds_pos2), replace = F)}, inds_pos2))
      stratagroup_numbers2 <- lapply(stratagroup_numbers2, function(x){setdiff(x,inds_pos2)})
      sample_inds_neg <- unlist(lapply(stratagroup_numbers2, function(x){sample(x, length(inds_neg2), replace = F)}, inds_neg2))
    }else{
      sample_inds_pos <- sample(toSample_pos, length(inds_pos2), replace = F)
      sample_inds_neg <- sample(toSample_neg, length(inds_neg2), replace = F)
    }
    cases_pos <- length(which(sample_inds_pos <= nCases))
    cases_neg <- length(which(sample_inds_neg > nCases))
    perm_cases <- cases_pos + cases_neg
    perm_controls <- (length(sample_inds_pos) - cases_pos) + (length(sample_inds_neg) - cases_neg)
    perm_score <- ifelse(flipped & method == 1, ValueTable[(controls1+perm_controls+1),(cases1+perm_cases+1)], ValueTable[(cases1+perm_cases+1),(controls1+perm_controls+1)])
    scores <- c(scores, perm_score)

  }

  pvalue <- length(which(scores >= score))/length(scores)

  return(pvalue)

}
