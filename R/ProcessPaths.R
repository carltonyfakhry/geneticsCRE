#' This function runs the causal engine for a given phenotype dataset.
#'
#' @description This function runs the causal engine for a given phenotype dataset.
#'
#' @usage GetBestPaths(dataset, nCases, nControls, method = 1,
#'        threshold_percent = 0.05, K = 10, pathLength = 5, iterations = 100,
#'        strataF = NA)
#'
#' @param dataset Filename of the phenotype data. First column must correspond
#'                to the gene symbols; all other columns correspond to the patients.
#'                First column must be named symbols and must contrain strings corresponding
#'                to the symbols of the genes. All other columns either contain a 0 if the
#'                gene is not present in the patient, and 1 or 2 if the gene is present
#'                in the patient.
#' @param nCases The number of cases (Must be greater than 1).
#' @param nControls The number of controls (Must be greater than 1).
#' @param path The path where the package will write and read from files
#'              created by the package.
#' @param method The method to be invoked. The value 1 calls the unsigned method
#'               while the value 2 calls the signed method.
#' @param threshold_percent All genes which occur in less than emph{threshold_percent}
#'                          of patients are kept for processing. The default
#'                          value is 0.05.
#' @param K The top K paths to be returned of each path length. The default
#'          value is 10.
#' @param pathLength The maximum path length for which the top K paths
#'                   will be computed. The default value is 5 (The maximum
#'                   value for \emph{pathLength} is 5).
#' @param iterations The number of iterations to be used in computing the
#'                   p-values. The default value is 100 (The minimum acceptable
#'                   value is 1).
#' @param strataF Filename for the strata file. \emph{strataF} must be a file
#'                of two columns, first column named \emph{subjid }contains
#'                the patient names similar to the dataset columns, and the
#'                second column named \emph{stratum} must contain the integer
#'                corresponding to the stratum.
#' @param nthreads The number of threads to be used in the parallel region
#'                 of the code. Default value is \emph{nthreads = NA}, in
#'                 which case \emph{nthreads} will be set to the maximum number of
#'                 threads that can be created. If \emph{nthreads} is greater
#'                 than the maximum number of threads that can be created, then
#'                 \emph{nthreads} will default back to the maximum number of threads
#'                 that can be created.
#'
#' @return This function returns a list with the following items:
#'         \item{SignedPaths}{The top K signed paths for each length.}
#'         \item{Paths}{The top K paths for each length.}
#'         \item{Lengths}{The lengths of each path.}
#'         \item{Scores}{The scores of each path.}
#'         \item{Pvalues}{The p-values of each path.}
#'
#' @author Carl Tony Fakhry
#'
#' @references Franceschini, A (2013). STRING v9.1: protein-protein interaction networks, with increased coverage
#'             and integration. In:'Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094.
#'             Epub 2012 Nov 29'.
#'
#'@export
GetBestPaths <- function(dataset, nCases, nControls, method = 'method1', threshold_percent = 0.05, K = 10, pathLength = 5, iterations = 100, strataF = NA, nthreads = NA){

  # accept numeric method ids to match old interface --dmitri
  if(is.numeric(method))
    method <- paste0('method', toString(method))
  else if(grepl('^enrich', method))
    method <- 'method1'
  else if(grepl('^direct', method))
    method <- 'method2'

  # Check the input parameters
  geneticsCRE:::check_input(nCases, nControls, method, threshold_percent, K, pathLength, iterations, nthreads)
  if(is.na(nthreads)) nthreads <- -1
  # path <- normalizePath(path)
  # if(substr(path, nchar(path),1) != "/") path <- paste(path, "/", sep = "")

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
  genes_data <- data.frame(genes = genes, data)
  patients <- dataset$patients
  rm(dataset)

  # Read in Stringdb (Already parsed)
  print("Processing network and dataset...")
  KB <- getStringKB()
  Rels <- KB$Rels
  Ents <- KB$Ents
  rm(KB)

  # Keep gene symbols in Ents which are present in the dataset
  Ents <- Ents[Ents$symbol != "-1",]
  Ents <- Ents[!duplicated(Ents$symbol),]
  Ents <- Ents[Ents$symbol %in% genes_data$genes,]
  Ents <- Ents[order(Ents$uid),]

  # Keep gene symbols in data which are present in Ents
  genes_data <- genes_data[genes_data$genes %in% Ents$symbol,]

  # Remove edges in Stringdb which have genes not present
  # in Ents
  Rels <- Rels[,c(1,2,3)]
  names(Rels) <- c("srcuid", "trguid", "sign")
  Rels <- unique(Rels)
  Rels <- Rels[Rels$srcuid %in% Ents$uid,]
  Rels2 <- Rels   # Rels2 could contain target genes which are not present in Ents (and dataset)
  Rels <- Rels[Rels$trguid %in% Ents$uid,] # Rels only contains both source and target genes present in Ents (and dataset)
  Rels <- Rels[Rels$srcuid != Rels$trguid,]

  # genes2, data2 and geneuids2 have the same order as Ents2$symbol
  # Keep only genes which are still present in source nodes of Network
  leftuids2 <- unique(Rels2$srcuid)
  Ents2 <- Ents[Ents$uid %in% leftuids2,]
  Ents2 <- Ents2[order(Ents2$uid),]
  genes_data2 <- genes_data[match(Ents2$symbol, genes_data$genes),]
  data2 <- as.matrix(genes_data2[,2:ncol(genes_data2),drop=FALSE])
  genes_data2 <- cbind(uid = Ents2$uid, genes_data2)

  # genes, data and geneuids have the same order as Ents$symbol
  # Keep only genes which are still present in Network
  leftuids <- unique(c(Rels$srcuid, Rels$trguid))
  Ents <- Ents[Ents$uid %in% leftuids,]
  Ents <- Ents[order(Ents$uid),]
  genes_data <- genes_data[match(Ents$symbol, genes_data$genes),]
  data <- as.matrix(genes_data[,2:ncol(genes_data),drop=FALSE])
  genes_data <- cbind(uid = Ents$uid, genes_data)

  stratagroups <- NA
  stratanumbers <- NA

  # Handling Strata File
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

  # Precompute the permuted columns to be used in producing randomized paths for computing
  # the p-values of the paths
  CaseORControl <- matrix(0, nrow = 0, ncol = 0)
  if(iterations > 0) {
    RandIndicesMat <- geneticsCRE:::getRandIndicesMat(ncol(data), iterations, strataF, nCases, nControls, stratagroups, stratanumbers)
    CaseORControl <- geneticsCRE:::getCaseORControl(RandIndicesMat, nCases, nControls)
  }

  UserKpaths <- c()
  UserKsignpaths <- c()
  UserSymbPaths <- c()
  UserSymbSignPaths <- c()
  UserScores <- c()
  UserLengths <- c()
  UserPvalues <- c()

  # Create the relations data frame for the input data. This data frame will be used in the joining
  # of the paths of length 1 and 2.
  Rels_data <- data.frame(srcuid = Ents$uid)
  Rels_data2 <- data.frame(srcuid = Ents2$uid)

  # Create the relations data frame for paths of length 3
  Rels <- Rels[order(Rels$srcuid, Rels$trguid),]
  Rels2 <- Rels
  names(Rels2) <- c("trguid", "trguid2", "sign2")

  srcuids1 <- Rels_data$srcuid
  trguids1 <- Rels_data$srcuid
  joining_gene_sign1 <- rep(1,length(Ents$uid))
  uids_CountLoc1 <- geneticsCRE:::getUidsCountsLocations(srcuids1, Rels_data$srcuid)
  data_inds1 <- geneticsCRE:::matchDataIndices(genes_data, Rels_data$srcuid) - 1

  srcuids1_2 <-Rels_data2$srcuid
  trguids1_2 <- Rels_data2$srcuid
  joining_gene_sign1_2  <- rep(1, length(Ents2$uid))
  uids_CountLoc1_2 <- geneticsCRE:::getUidsCountsLocations(srcuids1_2, Rels_data2$srcuid)
  data_inds1_2 <- geneticsCRE:::matchDataIndices(genes_data2, Rels_data2$srcuid) - 1

  srcuids2 <- Rels_data$srcuid
  trguids2 <- Rels_data$srcuid
  joining_gene_sign2 <- Rels$sign
  uids_CountLoc2 <- geneticsCRE:::getUidsCountsLocations(srcuids2, Rels$srcuid)
  data_inds2 <- geneticsCRE:::matchDataIndices(genes_data, Rels$trguid) - 1

  srcuids3 <- Rels$srcuid
  trguids3 <- Rels$trguid
  joining_gene_sign3 <- Rels$sign
  uids_CountLoc3 <- geneticsCRE:::getUidsCountsLocations(trguids3, Rels$srcuid)
  data_inds3 <- geneticsCRE:::matchDataIndices(genes_data, Rels$trguid) - 1

  # Create the relations data frame for paths of length 3
  Rels3 <- geneticsCRE:::getRels3(Rels$srcuid, Rels$trguid, Rels$sign, uids_CountLoc3)

  # Get the sign of the third gene in a path of length 3 assuming the first
  # gene in the path has a positive sign
  third_gene_sign <- rep(1, nrow(Rels3))
  third_gene_sign[Rels3$sign == -1 & Rels3$sign2 == 1] <- -1
  third_gene_sign[Rels3$sign == 1 & Rels3$sign2 == -1] <- -1

  srcuids4 <- Rels3$srcuid
  trguids4 <- Rels3$trguid2
  joining_gene_sign4 <- third_gene_sign
  uids_CountLoc4 <- geneticsCRE:::getUidsCountsLocations(trguids4, Rels$srcuid)


  srcuids5 <- Rels3$srcuid
  trguids5 <- Rels3$trguid2
  joining_gene_sign5 <- third_gene_sign
  uids_CountLoc5 <- geneticsCRE:::getUidsCountsLocations(trguids5, Rels3$srcuid)


  lsts <- geneticsCRE:::ProcessPaths(srcuids1, trguids1, uids_CountLoc1, joining_gene_sign1,
                                     srcuids1_2, trguids1_2, uids_CountLoc1_2, joining_gene_sign1_2,
                                     srcuids2, trguids2, uids_CountLoc2, joining_gene_sign2,
                                     srcuids3, trguids3, uids_CountLoc3, joining_gene_sign3,
                                     srcuids4, trguids4, uids_CountLoc4, joining_gene_sign4,
                                     srcuids5, trguids5, uids_CountLoc5, joining_gene_sign5,
                                     data_inds1, data_inds1_2, data_inds2, data_inds3,
                                     data, data2, ValueTable, nCases, nControls, K,
                                     iterations, CaseORControl, method, pathLength, nthreads)

  for(path_length in 1:pathLength){

    lst <- NULL
    lst_paths <- NULL

    if(path_length == 1){
      lst <- lsts$lst1
      lst_paths <- geneticsCRE:::getPaths(lst, path_length, Rels_data2, Rels_data2)
    }else if(path_length == 2){
      lst <- lsts$lst2
      lst_paths <- geneticsCRE:::getPaths(lst, path_length, Rels_data, Rels)
    }else if(path_length == 3){
      lst <- lsts$lst3
      lst_paths <- geneticsCRE:::getPaths(lst, path_length, Rels, Rels)
    }else if(path_length == 4){
      lst <- lsts$lst4
      lst_paths <- geneticsCRE:::getPaths(lst, path_length, Rels3, Rels)
    }else{
      lst <- lsts$lst5
      lst_paths <- geneticsCRE:::getPaths(lst, path_length, Rels3, Rels3)
    }
    UserScores <- c(UserScores, lst$scores)
    UserKpaths <- c(UserKpaths, lst_paths$path)
    if(path_length == 1){
      UserSymbPaths <- c(UserSymbPaths, geneticsCRE:::Uid2Symbol(Ents2, lst_paths$path))
    }else{
      UserSymbPaths <- c(UserSymbPaths, geneticsCRE:::Uid2Symbol(Ents, lst_paths$path))
    }
    UserKsignpaths <- c(UserKsignpaths, lst$signpath)
    if(path_length == 1){
      UserSymbSignPaths <- c(UserSymbSignPaths, geneticsCRE:::SignUid2SignSymbol(Ents2, lst_paths$signpath))
    }else{
      UserSymbSignPaths <- c(UserSymbSignPaths, geneticsCRE:::SignUid2SignSymbol(Ents, lst_paths$signpath))
    }
    TestScores <- lst$TestScores

    Pvalues <- sapply(lst$scores, function(x, y) {length(which(TestScores >= x))/length(TestScores)}, y = TestScores)
    UserPvalues <- c(UserPvalues, Pvalues)
    UserLengths <- c(UserLengths, rep(path_length, length(Pvalues)))

  }

  BestPaths <- data.frame(UserSymbSignPaths, UserSymbPaths, UserLengths, UserScores, UserPvalues, stringsAsFactors = F)
  BestPaths <- BestPaths[order(BestPaths$UserPvalues),]
  names(BestPaths) <- c("SignedPaths", "Paths", "Lengths", "Scores", "Pvalues")
  return(BestPaths)

}
