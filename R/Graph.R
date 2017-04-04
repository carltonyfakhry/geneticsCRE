# # Parse String-db from scratch
# parseStringNetwork <- function(){
#
#   # Get the full file name containing the STRINGdb relations
#   ff <- system.file("extdata", "9606.protein.actions.v10.txt.gz", package="QuaternaryProd")
#   all_rels <- read_tsv(gzfile(ff), col_names = TRUE)
#
#   # Set new names for columns
#   names(all_rels) <- c("srcuid", "trguid", "mode", "action", "direction","score")
#   Rels <- all_rels[, c("srcuid", "trguid", "mode", "direction")]
#
#   # Get all rows with causal relations
#   Rels <- Rels[Rels$mode %in% c("activation", "inhibition","expression"),]
#
#   # Get causal relations where direction is not specified, and consider reversed
#   # direction of causality as a valid causal relation
#   Bidirectional <- Rels[Rels$direction == 0 , c("trguid", "srcuid", "mode", "direction")]
#   names(Bidirectional) <- c("srcuid", "trguid", "mode", "direction")
#   Rels <- unique(bind_rows(Rels, Bidirectional))
#   Rels$direction <- NULL
#
#   # Rename activation as increases, inhibition as decreases, expression
#   # as regulates
#   Rels$mode <- sub("activation", "increases", Rels$mode)
#   Rels$mode <- sub("inhibition", "decreases", Rels$mode)
#   Rels$mode <- sub("expression", "regulates", Rels$mode)
#   Rels <- unique(Rels)
#
#   # Duplicate relations that are of type regulate
#   inds <- which(Rels$mode == "regulates")
#   Rels$mode[inds] <- "increases"
#   Rels2 <- Rels[inds,]
#   Rels2$mode <- rep("decreases", nrow(Rels2))
#   Rels <- rbind(Rels, Rels2)
#   Rels <- unique(Rels)
#   Rels$mode[Rels$mode == "increases"] = 1
#   Rels$mode[Rels$mode == "decreases"] = -1
#   Rels$mode <- as.integer(Rels$mode)
#
#
#   # Create entities dataframe
#   map <- org.Hs.egENSEMBLPROT2EG
#   proteins <- unique(Rels$srcuid)
#   proteins <- gsub("9606.", "", proteins)
#   proteinids <- unlist(mget(proteins, map, ifnotfound=NA))
#   proteinids[is.na(proteinids)] <- "-1"
#
#   mRNAs <- unique(Rels$trguid)
#   mRNAs <- gsub("9606.", "", mRNAs)
#   mRNAids <- unlist(mget(mRNAs, map, ifnotfound=NA))
#   mRNAids[is.na(mRNAids)] <- "-1"
#
#   id <- c(as.integer(proteinids), as.integer(mRNAids))
#
#   type <- c(rep("Protein", length(proteinids)), rep("mRNA", length(mRNAids)))
#
#   ensembleid <- c(names(proteinids), names(mRNAids))
#   ensembleid <- paste("9606.", ensembleid, sep = "")
#
#   uid <- c(1:length(id))
#
#   map <- org.Hs.egSYMBOL
#   symbol <- unlist(mget(c(proteinids, mRNAids), map, ifnotfound=NA))
#   symbol[is.na(symbol)] <- "-1"
#
#   Ents <- data.frame(uid, id, symbol, ensembleid, type)
#
#   # Make sure entities in Rnts are in Rels
#   # and entities in Rels are in Ents
#   srcents <- Ents[Ents$type == "Protein",]
#   trgents <- Ents[Ents$type == "mRNA",]
#
#   Rels <- Rels[Rels$srcuid %in% srcents$ensembleid,]
#   Rels <- Rels[Rels$trguid %in% trgents$ensembleid,]
#
#   srcents <- srcents[srcents$ensembleid %in% Rels$srcuid,]
#   trgents <- trgents[trgents$ensembleid %in% Rels$trguid,]
#   Ents <- rbind(srcents, trgents)
#
#   # Replace ensemble names in Rels with uids
#   srcs <- unique(Rels$srcuid)
#   srcents <- Ents[Ents$type == "Protein",]
#   for(i in srcs){
#     ind1 <- which(Rels$srcuid == i)
#     ind2 <- which(srcents$ensembleid == i)
#     if(length(ind2) > 1 || length(ind2) == 0){
#       stop("Duplicates found in srcents")
#     }
#     Rels$srcuid[ind1] <- srcents$uid[ind2]
#   }
#
#   trgs <- unique(Rels$trguid)
#   trgents <- Ents[Ents$type == "mRNA",]
#   for(i in trgs){
#     ind1 <- which(Rels$trguid == i)
#     ind2 <- which(trgents$ensembleid == i)
#     if(length(ind2) > 1 || length(ind2) == 0){
#       stop("Duplicates found in srcents")
#     }
#     Rels$trguid[ind1] <- trgents$uid[ind2]
#   }
#
#   Rels$srcuid <- as.integer(Rels$srcuid)
#   Rels$trguid <- as.integer(Rels$trguid)
#   Rels$mode <- as.integer(Rels$mode)
#   Rels <- unique(Rels)
#   Rels <- Rels[complete.cases(Rels),]
#   Rels <- cbind(Rels, uid = c(1:nrow(Rels)))
#
#   KB = List(Rels = Rels, Ents = Ents)
#   return(KB)
# }



#' This function loads the string-db network.
#'
#' @description This functions loads the string-db network.
#'
#' @export
getStringKB <- function(){
  Entspath <- system.file("extdata", "Ents.dat", package="geneticsCRE")
  Relspath <- system.file("extdata", "Rels.dat", package="geneticsCRE")
  Ents <- read.table(Entspath, stringsAsFactors = F)
  Rels <- read.table(Relspath, stringsAsFactors = F)
  KB <- list()
  KB$Ents <- Ents
  KB$Rels <- Rels
  return(KB)
}
