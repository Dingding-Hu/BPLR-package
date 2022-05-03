# R package BPLR
This package is developed to estimate the ROC curve and its summary statistics including AUC, Youden index and cutoff point based on a single biomarker. Default Bernstein polynomial approach deals with the case that the density ratio of x and y is monotone. Other approaches are also available. See more details below.
# Table of Contents
[Installation]

[Functions]

[Usage]

[References]

# Installation

Install this package from Github with

```r
# Install package "devtool" first if hasn't
#install.packages("devtools")
library(devtools)
devtools::install_github("Dingding-Hu/BPLR-package")
library(BPLR)
```


# Functions

This package contains the following functions:

- An approximation of the ROC curve and point estimate of the Youden index, cutoff point and AUC (`BPLR`) using different methods (by changing argument `method`)

# Usage

We provide two examples.

- Example 1: estimating the ROC curve and its summary statistics using the default "Bernstein" approach

```r
library(BPLR)
# generate observed data
x=rnorm(100,10,1)
y=rnorm(200,11.75,1)
# output the estimated ROC curve, Youden index, cutoff point and AUC in a list
BPLR(x,y)
```

- Example 2: estimating the ROC curve and its summary statistics using an empirical CDF based approach:

```r
BPLR(x,y,method="ECDF")
```




# References
arXiv:2205.00505, 
https://doi.org/10.48550/arXiv.2205.00505.


[Installation]: <https://github.com/Dingding-Hu/BPLR-package/blob/main/README.md#installation>
[Functions]: <https://github.com/Dingding-Hu/BPLR-package/blob/main/README.md#functions>
[Usage]: <https://github.com/Dingding-Hu/BPLR-package/blob/main/README.md#usage>
[References]: <https://github.com/Dingding-Hu/BPLR-package/blob/main/README.md#references>
