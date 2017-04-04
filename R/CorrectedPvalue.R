#' This function corrects P-values of the top K paths of each length.
#'
#' @description This function corrects P-values of the top K paths of each length.
#'
#' @usage getDecoratedPvalues(dataset, BestPaths, pathLength, nCases,
#'                                     nControls, method, n_permutations)
#'
#' @param dataset Filename of the phenotype data. First column must correspond
#'                to the gene symbols; all other columns correspond to the patients.
#'                First column must be named symbols and must contrain strings corresponding
#'                to the symbols of the genes. All other columns either contain a 0 if the
#'                gene is not present in the patient, and 1 or 2 if the gene is present
#'                in the patient.
#' @param BestPaths
#' @param pathLength The maximum path length for which the top K paths
#'                   will be computed. The default value is 5 (The maximum
#'                   value for \emph{pathLength} is 5).
#' @param nCases The number of cases (Must be greater than 1).
#' @param nControls The number of controls (Must be greater than 1).
#' @param method The method to be invoked. The value 1 calls the unsigned method
#'               while the value 2 calls the signed method.
#' @param n_permutations The number of permutations to be used in computing the
#'                   p-values. The default value is 100.
#'
#' @return This function returns a list with the following items:
#'         \item{SignedPaths}{The top K signed paths for each length.}
#'         \item{UserKpaths}{The top K paths for each length.}
#'         \item{UserLengths}{The lengths of each path.}
#'         \item{UserScores}{The scores of each path.}
#'         \item{UserPvalues}{The p-values of each path.}
#'
#' @author Carl Tony Fakhry
#'
#' @references Franceschini, A (2013). STRING v9.1: protein-protein interaction networks, with increased coverage
#'             and integration. In:'Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094.
#'             Epub 2012 Nov 29'.
#'
#'@export
getDecoratedPvalues <- function(dataset, BestPaths, pathLength, nCases, nControls, method, n_permutations){

  # Precompute values table
  print("Computing Values Table...")
  ValueTable <- geneticsCRE:::getValuesTable(nCases, nControls)

  # Preprocessing Phenotype dataset, all genes which occur in more
  # than threshold_percent of patients will be removed from the
  # dataset, also genes which do not occur in any patients will be removed
  print("Preprocessing Phenotype dataset...")
  dataset <- geneticsCRE:::PreprocessTable(dataset, threshold_percent, nCases, nControls)
  genes <- dataset$genes
  data <- dataset$data
  genes_data <- data.frame(genes = genes, data, stringsAsFactors = F)

  CorrectedBestPaths <- data.frame(SignedPaths = character(0), Paths = character(0), Subpaths1 = character(0),
                           Subpaths2 = character(0), Lengths = integer(0), Scores = numeric(0),
                           Pvalues = numeric(0), CorrectedPvalues = numeric(0))
  total_paths <- 1

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
      path_data_pos <- as.matrix(genes_data[match(genes, genes_data$genes),3:ncol(genes_data)],drop = FALSE)
      path_data_neg <- matrix(0, nrow = nrow(path_data_pos), ncol = ncol(path_data_pos))
      if(method == 2){
        path_data_pos[which(sign == -1),] = 0
        path_data_neg[which(sign == 1),] = 0
      }

      # Compute the p-values for the different splitted paths
      for(j in 1:(path_len-1)){
        CorrectedBestPaths <- rbind(CorrectedBestPaths, rep(0,ncol(CorrectedBestPaths)))
        CorrectedBestPaths[total_paths,c(1,2,5,6,7)] = bestpaths[i,]
        subpath_pos1 <- colSums(path_data_pos[1:j,,drop=FALSE])
        subpath_pos1[subpath_pos1 != 0] <- 1
        subpath_neg1 <- colSums(path_data_neg[1:j,,drop=FALSE])
        subpath_neg1[subpath_neg1 != 0] <- 1
        subpath_pos2 <- path_data_pos[j+1,]
        subpath_neg2 <- path_data_neg[j+1,]
        pvalue <- geneticsCRE:::computeCorrectedPvalue(subpath_pos1, subpath_neg1, subpath_pos2, subpath_neg2,
                                                         path_len, nCases, nControls, method, n_permutations, ValueTable)
        CorrectedBestPaths[total_paths,8] <- pvalue
        CorrectedBestPaths[total_paths,3] <- paste(genes[1:j], collapse = " -> ")
        CorrectedBestPaths[total_paths,4] <- paste(genes[j+1])
        total_paths <- total_paths + 1
        j <- j + 1
      }

    }

  }

  names(CorrectedBestPaths) = c("SignedPaths", "Paths", "Subpath1", "Subpath2", "Lengths", "Scores",
                                "Pvalues", "CorrectedPvalues")
  return(CorrectedBestPaths)

}



# This function computes the corrected p-value of two given sub-paths
computeCorrectedPvalue <- function(subpath_pos1, subpath_neg1, subpath_pos2, subpath_neg2, path_len,
                                    nCases, nControls, method, n_permutations, ValueTable){

  inds_pos1 <- which(subpath_pos1 != 0)
  inds_neg1 <- which(subpath_neg1 != 0)
  inds_pos2 <- which(subpath_pos2 != 0)
  inds_neg2 <- which(subpath_neg2 != 0)

  inds_pos2 <- setdiff(inds_pos2, inds_pos1)
  subpath_pos2[inds_pos2] <- 0
  inds_neg2 <- setdiff(inds_neg2, inds_neg1)
  subpath_neg2[inds_neg2] <- 0

  cases1 <- sum(subpath_pos1[1:nCases]) + sum(subpath_neg1[nCases:(nCases+nControls)])
  controls1 <- sum(subpath_pos1[nCases:(nCases+nControls)]) + sum(subpath_neg1[1:nCases])
  cases2 <- sum(subpath_pos2[1:nCases]) + sum(subpath_neg2[nCases:(nCases+nControls)])
  controls2 <- sum(subpath_pos2[nCases:(nCases+nControls)]) + sum(subpath_neg2[1:nCases])
  score <- ValueTable[(cases1+cases2+1),(controls1+controls2+1)]
  scores <- c()

  for(i in 1:n_permutations){

    sample_inds_pos <- sample(inds_pos2, length(inds_pos2), replace = F)
    sample_inds_neg <- sample(inds_neg2, length(inds_neg2), replace = F)
    perm_path_pos <- subpath_pos2
    perm_path_neg <- subpath_neg2
    perm_path_pos[sample_inds_pos] <- subpath_pos2[inds_pos2]
    perm_path_neg[sample_inds_neg] <- subpath_neg2[inds_neg2]
    perm_cases <- sum(perm_path_pos[1:nCases]) + sum(perm_path_neg[nCases:(nCases+nControls)])
    perm_controls <- sum(perm_path_pos[nCases:(nCases+nControls)]) + sum(perm_path_neg[1:nCases])
    perm_score <- ValueTable[(cases1+perm_cases+1),(controls1+perm_controls+1)]
    scores <- c(scores, perm_score)

  }

  pvalue <- length(which(scores >= score))/length(scores)

  return(pvalue)

}
