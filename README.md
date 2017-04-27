## geneticsCRE
geneticsCRE computes a statistic on causal graphs to identify rare and significant pathways.

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
devtools::install_github("carltonyfakhry/geneticsCRE", local = FALSE)
```
