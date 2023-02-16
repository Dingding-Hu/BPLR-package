# R package BPLR
This package is developed for the proposed Bernstein polynomial (BP) method in Hu et al. (2022) to estimate the ROC curve and its summary statistics including the AUC, Youden index, and optimal cutoff point based on a single biomarker under likelihood ratio ordering. The default BP approach deals with the case that the density ratio of two groups (diseased/non-diseased) of biomarkers is monotonically increasing. The code for other competitive method in  Hu et al. (2022) are also included. Details are provided below.

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

The function `BPLR` inputs the biomarkers from the diseased (`y`) and non-diseased group (`x`), and outputs a list consists of the ROC curve, AUC, Youden index and optimal cutoff point estimation. Here, we provide the code for the real data application in Hu et al.(2022+). The csv file of the DMD data is provided in the R folder of this package.

- Example: estimating the ROC curve and its summary statistics using the default "BPLR" approach

```r
library(BPLR)
# read in the data from the csv file
data=read.csv("DMD_CK.csv")

# store biomarker CK
CK=data$CK

# separate diseased and non-diseased biomarkers
x=CK[data1$Status..0....normal..1...carries.==0]
y=CK[data1$Status..0....normal..1...carries.==1]

# output the estimated ROC curve, Youden index, cutoff point and AUC in a list
BPLR(x,y,method="BP")
```

- Other methods are also available:

(1) The Box-Cox method in Bantis et al.(2019)
```r
BPLR(x,y,method="Box-Cox")
```

(2) The ZL method in Zhou and Lin(2008)
```r
BPLR(x,y,method="ZL")
```

(3) The ECDF method in Zhou et al.(2011)
```r
BPLR(x,y,method="ECDF")
```

(4) The MNLE method in Dykstra et al.(1995)
```r
BPLR(x,y,method="MNLE")
```

(5) The LZL method in Lin et al.(2012)
```r
BPLR(x,y,method="LZL")
```

(6) The Kernel method in Bantis et al.(2019)
```r
BPLR(x,y,method="Kernel")
```

(6) The MSLE method in Yu et al.(2017)
```r
BPLR(x,y,method="MSLE")
```

# References
D. Hu, M. Yuan, T. Yu, and P. Li (2022+), Statistical inference for the two-sample problem under
likelihood ratio ordering, with application to the ROC curve estimation. 

Dykstra R, Kochar S, Robertson T. Inference for likelihood ratio ordering in the two-sample problem. Journal of the
American Statistical Association 1995; 90: 1034-1040.

Bantis LE, Nakas CT, Reiser B. Construction of confidence intervals for the maximum of the Youden index and the
corresponding cutoff point of a continuous biomarker. Biometrical Journal 2019; 61: 138-156.

Yu T, Li P, Qin J. Density estimation in the two-sample problem with likelihood ratio ordering. Biometrika 2017; 104:
141-152.

Zhou XH, Lin H. Semi-parametric maximum likelihood estimates for ROC curves of continuous-scale tests. Statistics in
Medicine 2008; 27: 5271-5290.

Lin H, Zhou XH, Li G. A direct semiparametric receiver operating characteristic curve regression with unknown link and
baseline functions. Statistica Sinica 2012; 22: 1427-1456.

Zhou X, McClish DK, Obuchowski NA. Statistical Methods in Diagnostic Medicine. New York: Wiley. second ed. 2011.

[Installation]: <https://github.com/Dingding-Hu/BPLR-package/blob/main/README.md#installation>
[Functions]: <https://github.com/Dingding-Hu/BPLR-package/blob/main/README.md#functions>
[Usage]: <https://github.com/Dingding-Hu/BPLR-package/blob/main/README.md#usage>
[References]: <https://github.com/Dingding-Hu/BPLR-package/blob/main/README.md#references>
