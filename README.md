# geneticsCRE
geneticsCRE computes a statistic on causal graphs to identify rare and significant pathways.

## Installation
Before installing this package, make sure you have the latest versions of *Rstudio*, *R*, and the *devtools* package.

A recent version of Clang (https://clang.llvm.org/) is highly recommended. To configure the build environment, create **~/.R/Makevars**:
```javascript
CXX1X=clang++
CXX=clang++
CC=clang

CXX1XFLAGS=-O3 -march=native
CXXFLAGS=$(CXX1XFLAGS)
```

You can install this R pacakge using the following command:
```javascript
devtools::install_github("carltonyfakhry/geneticsCRE", local = FALSE)
```
