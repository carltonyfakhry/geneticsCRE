\name{geneticsCRE-package}
\alias{geneticsCRE-package}
\docType{package}
\title{
geneticsCRE is an R package that performs pathway-based genome-wide association study (PGWAS) to
        identify statistically significant associations between variants on a gene regulatory pathways
        and a given phenotype
}
\description{
geneticsCRE is an R package that performs pathway-based genome-wide association study (PGWAS) to
             identify statistically significant associations between variants on a gene regulatory pathways
             and a given phenotype. Unlike Genome-wide association study (GWAS), that seeks to assign
             statistical significance to associations of variations in single genes to a  phenotype, PGWAS
             accumulates statistical power by examining rare variant along gene-gene interaction pathways.
             PGWAS  uses prior causal information a gene regulatory interactions to infer statistically
             significant associations between causal pathways and a the phenotype. Given phenotype data with
             case/control information, geneticsCRE computes PGWAS for all valid pathways as identified by the
             Homo Sapien STRINGdb causal network. Examples on usage are provided in the vignette.
}
\details{
\tabular{ll}{
Package: \tab geneticsCRE\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2016-06-01\cr
License: \tab GPL (>= 3)\cr
}
}
\author{
Carl Tony Fakhry and Kourosh Zarringhalam

Maintainer: Carl Tony Fakhry <cfakhry@cs.umb.edu>
}
\references{

Franceschini, A (2013). STRING v9.1: protein-protein interaction networks, with increased coverage and integration. In:'Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094. Epub 2012 Nov 29'.
}
