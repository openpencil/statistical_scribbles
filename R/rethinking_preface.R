#### Preface ####

#### Load necessary packages for rethinking ####
install.packages(c("coda", "mvtnorm", "devtools"))
devtools::install_github("rmcelreath/rethinking")

#### Check if the C++ compiler works ####
library(rstan)
fx <- inline::cxxfunction( signature(x = "integer", y = "numeric" ) , '
    return ScalarReal( INTEGER(x)[0] * REAL(y)[0] ) ;
' )
# should be 10
fx( 2L, 5 )

#### Configuration for RStan from stan-dev pages on github ####
## https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Mac-or-Linux
dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
# Setting the optimization level to 3 may prevent some other R packages
# from installing from source if they are only tested with the
# stock R configuration
cat("\nCXXFLAGS=-O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function",
    file = M, sep = "\n", append = TRUE)
cat("\nCC=clang", "CXX=clang++ -arch x86_64 -ftemplate-depth-256",
    file = M, sep = "\n", append = TRUE)
## Optimization flags for C++
cat(readLines(M), sep = "\n")
cat(M)

#### Install RStan as specified on stan-dev pages on github ####
Sys.setenv(MAKEFLAGS = "-j4")
# note: omit the 's' in 'https' if you cannot handle https downloads
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
