% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BPLR.R
\name{ROC}
\alias{ROC}
\title{ROC curve analysis}
\usage{
ROC(x, y, method = "BP", nss = 10^4)
}
\arguments{
\item{x}{vector contains the sample of biomarkers from the "healthy" group}

\item{y}{vector contains the sample of biomarkers from the "diseased" group}

\item{method}{the method to estimate ROC curve and its summary statistics. It can be "BP","Box-Cox","ZL","ECDF","LZL","MNLE","Kernel",or "MSLE".}

\item{nss}{control the number of point estimates of the ROC curve in the range [0,1], default is 10^4.}
}
\description{
To estimate the ROC curve and its summary statistics including the AUC, Youden index, and optimal cutoff point based on a single biomarker by eight differet methods.
The default is the Bernstein polynomial approach in Hu et al. (2022), which is developed under the likelihood ratio ordering assumption.
}
\details{
The detail for each method is provided in Hu et al. (2022) or https://github.com/Dingding-Hu/BPLR-package.
}
\examples{
x=rnorm(100,10,1)
y=rnrom(100,12,1)
ROC(x,y,method="BP")

}
