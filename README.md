# R package BPLR
This package is developed for the proposed Bernstein polynomial (BP) method in Hu et al. (2022) to estimate the ROC curve and its summary statistics including the AUC, Youden index, and optimal cutoff point based on a single biomarker under likelihood ratio ordering. The default BP approach deals with the case that the density ratio of two groups (diseased/non-diseased) of biomarkers is monotonically increasing. The code for other competitive methods in  Hu et al. (2022) are also included. 

# Table of Contents
[Installation]

[Functions]

[Usage]

[References]

# Installation

Install this package from Github with

```r
# Install package "devtool" first if hasn't
# install.packages("devtools")
library(devtools)
devtools::install_github("Dingding-Hu/BPLR-package")
library(BPLR)
```


# Functions

This package contains the following functions:

- `ROC`: Return the estimates of the ROC curve, AUC, Youden index, and optimal cutoff point using different methods (by changing the argument `method`)

# Usage

The function `ROC` inputs the biomarkers from the non-diseased group (`x`) and the diseased group (`y`) , and outputs a list consisting of the ROC curve, AUC, Youden index, and optimal cutoff point estimation. In the following, we provide the code for the real data application in Hu et al. (2022). The csv file of the DMD data is provided in the R folder of this package. To import the data, one needs to download the csv file and put it in the working directory.

- Example: estimating the ROC curve and its summary statistics using the proposed "BP" approach for the DMD data

```r
library(BPLR)
# read in the data from the csv file
data=read.csv("DMD_CK.csv")

# store biomarker CK
CK=data$CK

# separate diseased and non-diseased biomarkers
x=CK[data$Status==0]
y=CK[data$Status==1]

# output the estimated ROC curve, AUC, Youden index and optimal cutoff point in a list
ROC(x,y,method="BP")
```

- Other methods can be implemented as follows:

(1) The Box-Cox method in Bantis et al. (2019)
```r
ROC(x,y,method="Box-Cox")
```

(2) The ZL method in Zhou and Lin (2008)
```r
ROC(x,y,method="ZL")
```

(3) The LZL method in Lin et al. (2012)
```r
ROC(x,y,method="LZL")
```

(4) The ECDF method
```r
ROC(x,y,method="ECDF")
```

(5) The MNLE method in Dykstra et al. (1995)
```r
ROC(x,y,method="MNLE")
```

(6) The Kernel method in Bantis et al. (2019)
```r
ROC(x,y,method="Kernel")
```

(7) The MSLE method in Yu et al. (2017)
```r
ROC(x,y,method="MSLE")
```

# References

Bantis LE, Nakas CT, Reiser B (2019), "Construction of confidence intervals for the maximum of the Youden index and the
corresponding cutoff point of a continuous biomarker", Biometrical Journal, 61, 138-156.

Dykstra R, Kochar S, Robertson T (1995), "Inference for likelihood ratio ordering in the two-sample problem", Journal of the
American Statistical Association, 90, 1034-1040.

Hu D, Yuan M, Yu T, Li P (2022), "Statistical inference for the two-sample problem under
likelihood ratio ordering, with application to the ROC curve estimation", [arXiv:2205.00505](http://arxiv.org/abs/2205.00505).

Lin H, Zhou XH, Li G (2012), "A direct semiparametric receiver operating characteristic curve regression with unknown link and
baseline functions", Statistica Sinica, 22, 1427-1456.

Yu T, Li P, Qin J (2017), "Density estimation in the two-sample problem with likelihood ratio ordering", Biometrika, 104,
141-152.

Zhou XH, Lin H (2008), "Semi-parametric maximum likelihood estimates for ROC curves of continuous-scale tests", Statistics in
Medicine, 27, 5271-5290.




[Installation]: <https://github.com/Dingding-Hu/BPLR-package/blob/main/README.md#installation>
[Functions]: <https://github.com/Dingding-Hu/BPLR-package/blob/main/README.md#functions>
[Usage]: <https://github.com/Dingding-Hu/BPLR-package/blob/main/README.md#usage>
[References]: <https://github.com/Dingding-Hu/BPLR-package/blob/main/README.md#references>
