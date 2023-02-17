

#' Estimating the ROC curve
#'
#' To estimate ROC curve and its summary statistics including AUC, Youden index and cutoff point based on a single biomarker.
#' Default Bernstein polynomial approach deals with the case that the density ratio of x and y is monotone.
#' Other approaches are also available. See more details below.
#'
#' @param x vector contains the sample of biomarkers from the "healthy" group
#' @param y vector contains the sample of biomarkers from the "diseased" group
#' @param method the method to be used to estimate ROC curve and its summary statistics.
#'        The default method "BP" uses the Bernstein polynomials to approximate the log of density ratio incorporating the
#'        likelihood ratio ordering.
#'        In each case, the order of the Bernstein polynomials being used is selected by a BIC criteria.
#' @param nss control the number of point estimates of the ROC curve in the range [0,1], default is 10^4.
#'
#' @details The function
#'
#'
#' @examples x=rnorm(100,10,1)
#' y=rnrom(100,12,1)
#' ROC(x,y,method="Bernstein",nss=10^4)
#'
#' @import glmnet fdrtool rootSolve
#'
#' @export



ROC=function(x,y,method="BP",nss=10^4) {
  #check for arguments
  stopifnot("needs to use one of provided methods"= (method %in% c("BP","Box-Cox","ZL","ECDF","LZL","MNLE","Kernel","MSLE")))
  stopifnot("too many points chosen"= (nss<=10^6))

  if( method=="BP"){
    estbp(x,y,nss=nss)
  } else if( method=="ZL"){
    estzl(x,y,nss=nss)
  } else if( method=="Box-Cox"){
    estbc(x,y,nss=nss)
  } else if(method=="ECDF") {
    estemp(x,y,nss=nss)
  } else if(method=="LZL") {
    estsinica(x,y,nss=nss)
  } else if(method=="MNLE") {
    estlr(x,y,nss=nss)
  } else if(method=="Kernel") {
    estkern(x,y,nss=nss)
  } else if(method=="MSLE"){
    estsmlr(x,y,nss=nss)
  }
}


