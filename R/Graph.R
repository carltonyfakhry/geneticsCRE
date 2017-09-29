# This function loads the string-db network.
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
