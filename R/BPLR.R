

#' Estimating the ROC curve
#'
#' To estimate ROC curve and its summary statistics including AUC, Youden index and cutoff point based on a single biomarker.
#' Default Bernstein polynomial approach deals with the case that the density ratio of x and y is monotone.
#' Other approaches are also available. See more details below.
#'
#' @param x vector contains the sample of biomarkers from the "healthy" group
#' @param y vector contains the sample of biomarkers from the "diseased" group
#' @param method the method to be used to estimate ROC curve and its summary statistics.
#'        The default method "Bernstein" uses the Bernstein polynomials to approximate the log of density ratio incorporating the
#'        likelihood ratio ordering.
#'
#' @param type control the type of Bernstein polynomials used to approximate the log density ratio.
#'        "type=1" uses the Bernstein polynomial of x and y.
#'        "type=2" uses both Bernstein polynomial of x, y and log(x), log(y).
#'         In each case, the order of the Bernstein polynomials being used is selected by a BIC criteria.
#' @param nss control the number of point estimates of the ROC curve in the range [0,1], default is 10^4.
#'
#' @details The function
#'
#'
#' @examples x=rnorm(100,10,1)
#' y=rnrom(100,12,1)
#' BPLR(x,y,method="Bernstein",nss=10^4,type=1)
#'
#' @import glmnet fdrtool rootSolve
#'
#' @export



BPLR=function(x,y,method="Bernstein",nss=10^4,type=2) {
  #check for arguments
  stopifnot("needs to use one of provided methods"= (method %in% c("Bernstein","Box-Cox","ECDF","MLELR","Kernel","MSLELR")))
  stopifnot("too many points chosen"= (nss<=10^6))
  stopifnot("type should be 1 or 2 "= (type %in% 1:2))

  disroc=function(s,hatF0,hatF1)
  {
    n=length(hatF0)

    ind=min(c((1:n)[hatF0>=1-s],n))
    1-hatF1[ind]
  }

  estemp= function(x,y,nss=10^4)
  {
    hatF0 = as.stepfun(ecdf(x))
    hatF1 = as.stepfun(ecdf(y))

    tt = sort(unique(c(x,y)) )


    estF0=hatF0(tt)
    estF1=hatF1(tt)
    ss=((1:nss)-0.5)/nss

    ROC=as.numeric(sapply(ss,disroc,hatF0=estF0,hatF1=estF1))
    AUC=mean(ROC)
    n=length(tt)

    dif = hatF0(tt)-hatF1(tt)
    hatJ = max(dif)

    index = c(1:length(tt))[abs(dif-hatJ)<1/n^2 ]

    hatC = mean(tt[index])

    list(ROC=ROC,AUC=AUC,Youden=hatJ,cutoff=hatC)
  }


  pava=function(y,w)
  {
    ##Adapted from the isotone function in Package DDHFm##
    x=y
    wt=w
    nn=length(x)
    ip <- (1:nn)
    dx <- diff(x)
    nx <- length(x)
    while ((nx > 1) && (min(dx) < 0)) {
      jmax <- (1:nx)[c(dx <= 0, FALSE) & c(TRUE, dx > 0)]
      jmin <- (1:nx)[c(dx > 0, TRUE) & c(FALSE, dx <= 0)]
      for (jb in (1:length(jmax))) {
        ind <- (jmax[jb]:jmin[jb])
        wtn <- sum(wt[ind])
        x[jmax[jb]] <- sum(wt[ind] * x[ind])/wtn
        wt[jmax[jb]] <- wtn
        x[(jmax[jb] + 1):jmin[jb]] <- NA
      }
      ind <- !is.na(x)
      x <- x[ind]
      wt <- wt[ind]
      ip <- ip[ind]
      dx <- diff(x)
      nx <- length(x)
    }
    jj <- rep(0, nn)
    jj[ip] <- 1
    z <- x[cumsum(jj)]
    return(z)
  }


  datafreq=function(Tvalue,group)
  {
    tt=sort(unique(Tvalue))
    delta=tt
    tot=tt
    for(i in 1:length(tt))
    {
      delta[i]=sum((Tvalue==tt[i])&(group==1))
      tot[i]=sum(Tvalue==tt[i])
    }
    cbind(tt,delta,tot)
  }

  estlr=function(x,y,nss=10^4)
  {

    Tvalue=c(x,y)
    group=c(rep(0,length(x)),rep(1,length(y)) )

    out=datafreq(Tvalue,group)
    tt=out[,1]
    z=out[,2]/out[,3]
    w=out[,3]
    theta=pava(z,w)

    lam=length(y)/length(Tvalue)
    n=length(Tvalue)

    hatp1=w*theta/n/lam
    hatp0=w*(1-theta)/n/(1-lam)



    estF0=cumsum(hatp0)
    estF1=cumsum(hatp1)

    ss=((1:nss)-0.5)/nss
    ROC=as.numeric(sapply(ss,disroc,hatF0=estF0,hatF1=estF1))
    AUC=mean(ROC)
    hatF1=cumsum(hatp1)
    hatF0=cumsum(hatp0)

    dif = hatF0-hatF1
    hatJ = max(dif)

    index = c(1: length(tt) )[abs(dif-hatJ)<1/n^2 ]

    hatC = mean(tt[index])

    list(ROC=ROC,AUC=AUC,Youden=hatJ,cutoff=hatC)
  }


  kerncdf=function(nx,x,type="d")
  {
    ##Kernel estimation of pdf (type="d") or CDF (type="p") in Yuan et al's CJS paper##


    sx = sd(x)
    iqrx = IQR(x)
    n = length(x)
    h= 0.9*min(sx, (iqrx/1.34))*n^(-0.2)

    if(type=="d")
    {

      obj=function(tt)
      {
        mean(dnorm(tt,x,h))
      }

    }

    if(type=="p")
    {
      obj=function(tt)
      {
        mean(pnorm( (tt-x)/h ))
      }

    }

    out=sapply(nx,obj)

    as.numeric(out)
  }


  estkern= function(x,y,nss=10^4)
  {
    sx = sd(x)
    iqrx = IQR(x)
    nx = length(x)
    hx= 0.9*min(sx, (iqrx/1.34))*nx^(-0.2)

    sy = sd(y)
    iqry = IQR(y)
    ny = length(y)
    hy= 0.9*min(sy, (iqry/1.34))*ny^(-0.2)
    h=max(hx,hy)

    data=c(x,y)
    tt=seq(min(data)-4*h,max(data)+4*h,length=2*10^4)
    estF0=kerncdf(tt,x,type="p")
    estF1=kerncdf(tt,y,type="p")

    ss=((1:nss)-0.5)/nss
    ROC=as.numeric(sapply(ss,disroc,hatF0=estF0,hatF1=estF1))
    AUC=mean(ROC)
    obj=function(tt)
    {
      kerncdf(tt,x,type="d")- kerncdf(tt,y,type="d")
    }

    ts=uniroot.all(obj,c(min(data),max(data)))
    estF1 = kerncdf(ts, y,type="p")
    estF0=  kerncdf(ts, x,type="p")

    hatJ=max(estF0-estF1)
    hatC=mean(ts[estF0-estF1>=hatJ])

    list(ROC=ROC,AUC=AUC,Youden=hatJ,cutoff=hatC)
  }



  estsmlr=function(x,y,nss=10^4)
  {

    sx = sd(x)
    iqrx = IQR(x)
    nx = length(x)
    hx= 0.9*min(sx, (iqrx/1.34))*nx^(-0.2)

    sy = sd(y)
    iqry = IQR(y)
    ny = length(y)
    hy= 0.9*min(sy, (iqry/1.34))*ny^(-0.2)
    h=max(hx,hy)

    data=c(x,y)
    tt=seq(min(data)-4*h,max(data)+4*h,length=2*10^4)

    lam=length(y)/(length(x)+length(y))

    estF1 = kerncdf(tt, y,type="p")
    estF0=  kerncdf(tt, x,type="p")
    estW = (1-lam)*estF0 + lam*estF1

    indW=duplicated(estW)
    ut=tt[!indW]
    uW=estW[!indW]
    uF1=estF1[!indW]
    uF0=estF0[!indW]

    lcmfit = gcmlcm(uW, lam*uF1, type = "gcm")

    estF=function(nw)
    {
      xw=lcmfit$x.knots
      sp=lcmfit$slope.knots

      ind=max(c(1,(1:length(xw))[xw<nw]) )
      dxnw =nw-xw[ind]
      spnw=sp[ind]

      esty=lcmfit$y.knots[ind]+dxnw*spnw

      esty
    }

    outy=sapply(uW,estF)

    estF1=outy/lam
    estF0=(uW-outy)/(1-lam)
    ss=((1:nss)-0.5)/nss
    ROC=as.numeric(sapply(ss,disroc,hatF0=estF0,hatF1=estF1))
    AUC=mean(ROC)
    hatdif=abs(lcmfit$slope.knots-lam)

    len=length(hatdif)
    index=(1:len)[hatdif<=min(hatdif)]

    ind=index[1]
    hatW=lcmfit$x.knots[ind]
    hatF1=lcmfit$y.knots[ind]/lam
    hatJ=(hatW-hatF1)/(1-lam)
    indC=(1:nss)[estW==hatW]
    hatC=mean(tt[indC])

    list(ROC=ROC,AUC=AUC,Youden=hatJ,cutoff=hatC)

  }

  ###X-matrix based on K+1 Bernstein polynomials###
  BernPoly=function(x,K)
  {
    obj=function(sx)
    {
      dbinom(0:K,K,prob=sx)
    }
    out=t(sapply(x,obj))
    out
  }

  glmplus=function(group,xmat)
  {

    xmat=as.matrix(xmat)
    K=ncol(xmat)
    if(K==1)
    {
      out=suppressWarnings(glm(group~xmat,family="binomial"))
      if(out$coef[2]<0)
      {
        out=glm(group~1,family="binomial")
        hatp1=out$fitted.values
        hatp0=1-hatp1
        estb=as.numeric(c(out$coef,0))
        df=0
      }
      else
      {
        hatp1=out$fitted.values
        hatp0=1-hatp1
        estb=as.numeric(out$coef)
        df=1
      }
    }

    else
    {
      out=glmnet(xmat,group,family="binomial",lambda=0,lower.limits=0,standardize=F)
      estb=as.numeric(c(as.numeric(out$a0),as.numeric(out$beta) ))

      hatp1=predict(out,newx=xmat,type="response")
      hatp0=1-hatp1
      df=out$df
    }

    list(hatp0=hatp0,hatp1=hatp1,estb=estb,df=df)
  }


  ICBP=function(x,y,K,type=1)
  {
    n0=length(x)
    n1=length(y)
    tt=c(x,y)
    group=c(rep(0,n0),rep(1,n1))

    ind=order(tt)
    tt=tt[ind]
    group=group[ind]

    n=n0+n1
    lam=n1/n

    t01=( tt-min(tt) )/(max(tt)-min(tt))

    logt01=( log(tt)-log( min(tt) )  )/( log( max(tt)) - log( min(tt) ) )


    tmat=BernPoly(t01,K)
    logtmat=BernPoly(logt01,K)

    ###(K+1)*(K+1) lower-trainagle matrix###
    Bmat=1-upper.tri(matrix(1,nrow=K+1,ncol=K+1), diag = F)

    ###X-matrix %*% Bmat###

    trtmat=tmat%*%Bmat
    tlgmat=logtmat%*%Bmat

    if(type==1)
    {
      xmat=trtmat[,-1]
    }

    if(type==2)
    {
      xmat=cbind(trtmat[,-1],tlgmat[,-1])
    }


    out=glmplus(group,xmat)

    hatp1=out$hatp1
    hatp0=1-hatp1
    df=out$df

    lik=sum(group*log(hatp1+1e-30)+(1-group)*log(hatp0+1e-30))

    bic=-2*lik+log(n)*df

    c(lik,bic)
  }

  BPLR=function(x,y,type=1,nss=10^4)
  {
    n0=length(x)
    n1=length(y)
    tt=c(x,y)
    group=c(rep(0,n0),rep(1,n1))

    output=c()
    for(i in 1:15)
    {
      output=rbind(output,ICBP(x,y,i,type=type))
    }

    IC=min(output[,2])
    K=(1:15)[output[,2]<=IC][1]


    ind=order(tt)
    tt=tt[ind]
    group=group[ind]
    n=n0+n1
    lam=n1/n


    t01=( tt-min(tt) )/(max(tt)-min(tt))

    logt01=( log(tt)-log( min(tt) )  )/( log( max(tt)) - log( min(tt) ) )

    tmat=BernPoly(t01,K)
    logtmat=BernPoly(logt01,K)

    ###(K+1)*(K+1) lower-trainagle matrix###
    Bmat=1-upper.tri(matrix(1,nrow=K+1,ncol=K+1), diag = F)

    ###X-matrix %*% Bmat###

    trtmat=tmat%*%Bmat
    tlgmat=logtmat%*%Bmat

    if(type==1)
    {
      xmat=trtmat[,-1]
    }

    if(type==2)
    {
      xmat=cbind(trtmat[,-1],tlgmat[,-1])
    }


    out=glmplus(group,xmat)

    hatp1=out$hatp1
    hatp0=1-hatp1

    hatF=function(x,prob)
    {
      sum(prob[tt<=x])
    }

    hatp1=hatp1/n1
    hatp0=hatp0/n0
    ut=sort(unique(tt))
    estF0=sapply(ut,hatF,prob=hatp0)
    estF1=sapply(ut,hatF,prob=hatp1)
    ss=((1:nss)-0.5)/nss
    ROC=as.numeric(sapply(ss,disroc,hatF0=estF0,hatF1=estF1))
    AUC=mean(ROC)

    if(type==1)
    {
      cutobj=function(x)
      {

        estbeta=Bmat%*%(out$estb)

        sx=(x-min(tt))/(max(tt)-min(tt))

        part1=BernPoly(sx,K)%*%estbeta
        part1-log( lam/(1-lam) )
      }

    }


    if(type==2)
    {

      cutobj=function(x)
      {
        estb=as.numeric( (out$estb)[-1])
        estb1=estb[1:K]
        estb2=estb[-(1:K)]
        estb0=as.numeric(out$estb [1])

        nb1=c(estb0/2,estb1)
        nb2=c(estb0/2,estb2)
        estbeta1=Bmat%*%nb1
        estbeta2=Bmat%*%nb2

        sx=(x-min(tt))/(max(tt)-min(tt))
        logsx=(log(x)-log(min(tt)) )/(log(max(tt))- log(min(tt)) )

        part1=BernPoly(sx,K)%*%estbeta1
        part2=BernPoly(logsx,K)%*%estbeta2
        part1+part2-log( lam/(1-lam) )
      }

    }

    out=uniroot.all(cutobj,lower=min(tt),upper=max(tt))

    hatJ0=as.numeric( sapply(out,hatF,prob=hatp0) )
    hatJ1=as.numeric( sapply(out,hatF,prob=hatp1) )


    hatJ=max( hatJ0-hatJ1 )

    hatC=mean(out[ hatJ0-hatJ1 >=hatJ])


    list(ROC=ROC,AUC=AUC,Youden=hatJ,cutoff=hatC,BIC=IC)

  }

  estbp=function(x,y,type=1,nss=10^4)
  {
    out1=BPLR(x,y,type=1,nss=nss)
    out2=BPLR(x,y,type=2,nss=nss)
    if(type==1){
      out1
    } else if(type==2){
      out2
    }
  }





  pow=function(y,lam)
  {
    if(abs(lam)<1e-5)
    {
      tem=log(y)
    }
    else
    {
      tem=( (y)^lam-1 )/lam
    }
    return(tem)
  }

  estbc=function(x,y,nss=10^4)
  {

    obj=function(lam){
      tx=pow(x,lam)
      ty=pow(y,lam)
      mx = mean(tx)
      varx = mean((tx-mx)^2)
      my = mean(ty)
      vary = mean((ty-my)^2)
      -length(tx)/2*log(varx)-length(ty)/2*log(vary) + (lam-1)*(sum(log(x))+sum(log(y)))
    }

    model = optimize(obj,c(-5,5),maximum = T)
    lam = model$maximum

    tx=pow(x,lam)
    ty=pow(y,lam)

    mx = mean(tx)
    sx = sqrt(mean((tx-mx)^2))
    my = mean(ty)
    sy = sqrt(mean((ty-my)^2))
    ss=((1:nss)-0.5)/nss
    ROC=1-pnorm(qnorm(1-ss,mx,sx),my,sy)
    AUC=mean(ROC)
    difpdf=function(tt)
    {
      dnorm(tt,mx,sx)-dnorm(tt,my,sy)
    }

    low=min(c(tx,ty))
    upp=max(c(tx,ty))

    ts=uniroot.all(difpdf,c(low,upp))

    estF1 = pnorm(ts, my,sy)
    estF0=  pnorm(ts, mx,sx)

    hatJ=max(estF0-estF1)
    c=mean(ts[estF0-estF1>=hatJ])

    hatC = (c*lam+1)^(1/lam)

    list(ROC=ROC,AUC=AUC,Youden=hatJ,cutoff=hatC)
  }
  if( method=="Bernstein"){
    estbp(x,y,nss=nss,type=type)
  } else if( method=="Box-Cox"){
    estbc(x,y,nss=nss)
  } else if(method=="ECDF") {
    estemp(x,y,nss=nss)
  } else if(method=="MLELR") {
    estlr(x,y,nss=nss)
  } else if(method=="Kernel") {
    estkern(x,y,nss=nss)
  } else if(method=="MSLELR"){
    estsmlr(x,y,nss=nss)
  }
}


