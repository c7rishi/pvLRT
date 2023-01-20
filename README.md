# R package `pvLRT`

 ![](http://cranlogs.r-pkg.org/badges/grand-total/pvLRT)



`pvLRT` is an R package that implements a suite of likelihood ratio test based methods to use in pharmacovigilance. It can handle adverse events data on several simultaneous drugs, with possibly zero-inflated report counts. The signals identified through these methods enjoy rigorous statistical gurantees, inlcuding strictly controlled type I error and false discovery rate.  Several testing and post-processing functions are implemented.

## Installation

The stable version of `pvLRT` can be installed from CRAN:

```
install.packages("pvLRT")
```

The development version is available from GitHub:

```
# if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("c7rishi/pvLRT")
```
