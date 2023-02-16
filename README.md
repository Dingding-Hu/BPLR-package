# R package BPLR
This package is developed for the proposed method in Hu et al.(2022+) to estimate the ROC curve and its summary statistics including AUC, Youden index and cutoff point based on a single biomarker. The default Bernstein polynomial approach deals with the case that the density ratio of two groups of biomarkers is monotone (diseased/non-diseased). The codes for other approaches in the simulation of Hu et al.(2022+) are also provided. Details are provided below.

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

- `BPLR`: An approximation of the ROC curve and point estimate of the Youden index, cutoff point and AUC using different methods (by changing the argument `method`)

# Usage

The function `BPLR` inputs the biomarkers from the diseased (`y`) and non-diseased group (`x`), and outputs a list consists of the ROC curve, AUC, Youden index and optimal cutoff point estimation. Here, we provide the code for the real data application in Hu et al.(2022+). The csv file of the data is provided in the R folder of this package.

- Example: estimating the ROC curve and its summary statistics using the default "BPLR" approach

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
