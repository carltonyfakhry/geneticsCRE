## geneticsCRE
An R package that performs pathway-based genome-wide association study (PGWAS) to identify statistically significant associations between variants on gene regulatory pathways and a given phenotype. Unlike Genome-wide association study (GWAS), that seeks to assign statistical significance to associations of variations in single genes to a phenotype, PGWAS accumulates statistical power by examining rare variants along gene-gene interaction pathways. PGWAS uses prior causal information from gene regulatory interactions to infer statistically significant associations between causal pathways and a phenotype. Given phenotype data with case/control information, geneticsCRE computes PGWAS for all valid pathways as identified by the Homo Sapien STRINGdb causal network. Examples on usage are provided in the vignette.

### Quickstart

geneticsCRE can be run from a Docker image:

    docker run -it -v </path/to/workspace>:/work gcre/gcre:latest

The image is about 1GB in size, and provides a pre-configured Rcpp development environment. On startup, it will install the R package from GitHub - using build settings appropriate for the hardware - and then drop into a shell.

* `/path/to/workspace` is the host directory that contains any necessary data and scripts.
* If running on macOS/OSX, **Docker for Mac** (https://www.docker.com/docker-mac) should work better than the legacy "Docker Toolbox".
* More details, build scripts, and Dockerfiles can be found here: https://github.com/dbichko/gcre-ci
* Docker Hub organization: https://hub.docker.com/r/gcre/

It's possible to bypass the automatic installation by specifying a Docker run command (most likely a shell); the install can then be run manually:

```
docker run -it -v </path/to/workspace>:/work gcre/gcre:latest bash
docker@09caf52a75ef:~$ ./install
```

### Installation
Before installing this package, make sure you have the latest versions of *Rstudio*, *R*, and the *devtools* package.

A recent version of Clang (https://clang.llvm.org/) is highly recommended. To configure the build environment, create **~/.R/Makevars**:
```properties
CXX1X=clang++
CXX=clang++
CC=clang

CXX1XFLAGS=-O3 -march=native
CXXFLAGS=$(CXX1XFLAGS)
```

You can then install the R package using the following command:
```r
devtools::install_github("carltonyfakhry/geneticsCRE", build_vignettes = TRUE, local = FALSE)
```

For an introduction to `geneticsCRE` and for an example on how to compute GWASPA over the publicly available network STRINGdb, please see the Vignette for this package using the following:

```r
browseVignettes("geneticsCRE")
```
