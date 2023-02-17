#' ROC curve analysis
#'
#' To estimate the ROC curve and its summary statistics including the AUC, Youden index, and optimal cutoff point based on a single biomarker by eight differet methods.
#' The default is the Bernstein polynomial approach in Hu et al. (2022), which is developed under the likelihood ratio ordering assumption.
#'
#' @param x vector contains the sample of biomarkers from the "healthy" group
#' @param y vector contains the sample of biomarkers from the "diseased" group
#' @param method the method to estimate ROC curve and its summary statistics. It can be "BP","Box-Cox","ZL","ECDF","LZL","MNLE","Kernel",or "MSLE".
#' @param nss control the number of point estimates of the ROC curve in the range [0,1], default is 10^4.
#'
#' @details The detail for each method is provided in Hu et al. (2022) or 
#'
#'
#' @examples x=rnorm(100,10,1)
#' y=rnrom(100,12,1)
#' ROC(x,y,method="BP")
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
