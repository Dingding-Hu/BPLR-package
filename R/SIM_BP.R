# R codes for BP methods and competitive methods in the simulations


############################# Some functions that are used in some of the methods below


disroc=function(s,hatF0,hatF1)
{
  n=length(hatF0)

  ind=min(c((1:n)[hatF0>=1-s],n))
  1-hatF1[ind]
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

############################## Berstein polynomial methods (BPLR)

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

# for oreder selection
glmplus=function(group,xmat)
{

  xmat=as.matrix(xmat)
  out=glmnet(xmat,group,family="binomial",lambda=0,lower.limits=0,standardize =FALSE)
  estb=as.numeric(c(as.numeric(out$a0),as.numeric(out$beta) ))

  hatp1=predict(out,newx=xmat,type="response")
  hatp0=1-hatp1
  df=out$df

  list(hatp0=hatp0,hatp1=hatp1,estb=estb,df=df)
}


# for order selection
glmplus2=function(group,xmat)
{
  xmat=as.matrix(xmat)
  out=glmnet(xmat,group,family="binomial",alpha = 0, lambda = 1e-06,standardize =FALSE)
  estb=as.numeric(c(as.numeric(out$a0),as.numeric(out$beta) ))

  hatp1=predict(out,newx=xmat,type="response")
  hatp0=1-hatp1
  df=out$df

  list(hatp0=hatp0,hatp1=hatp1,estb=estb,df=df)
}

# calculate the BIC
ICBP=function(x,y,K)
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
  xmat=cbind(trtmat[,-1],tlgmat[,-1])

  out=glmplus2(group,xmat)

  hatp1=out$hatp1
  hatp0=1-hatp1
  df=out$df

  lik=sum(group*log(hatp1+1e-30)+(1-group)*log(hatp0+1e-30))

  bic=-2*lik+log(n)*df

  c(lik,bic)
}


# estimation of the ROC curve and its summary statistics
estbp=function(x,y,nss=10^4)
{
  n0=length(x)
  n1=length(y)
  tt=c(x,y)
  group=c(rep(0,n0),rep(1,n1))

  output=c()
  for(i in 1:5)
  {
    output=rbind(output,ICBP(x,y,i))
  }

  IC=min(output[,2])
  K=(1:5)[output[,2]<=IC][1]


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

  xmat=cbind(trtmat[,-1],tlgmat[,-1])

  out=glmplus(group,xmat)

  hatp1=out$hatp1
  hatp0=1-hatp1

  hatF=function(x,prob)
  {
    sum(prob[tt<=x])
  }

  hatp1=hatp1/sum(hatp1)
  hatp0=hatp0/sum(hatp0)
  ut=sort(unique(tt))
  estF0=as.numeric(sapply(ut,hatF,prob=hatp0))
  estF1=as.numeric(sapply(ut,hatF,prob=hatp1))
  ss=((1:nss)-0.5)/nss
  ROC=cbind(ss,as.numeric(sapply(ss,disroc,hatF0=estF0,hatF1=estF1)))
  pi=c(estF0[1],diff(estF0))
  qi=c(estF1[1],diff(estF1))
  AUC=sum((estF0-0.5*pi)*qi)

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

  outroot=uniroot.all(cutobj,lower=min(tt),upper=max(tt))
  hatJ0=as.numeric( sapply(outroot,hatF,prob=hatp0) )
  hatJ1=as.numeric( sapply(outroot,hatF,prob=hatp1) )


  hatJ=max( hatJ0-hatJ1 )

  hatC=mean(outroot[ (hatJ0-hatJ1) >=hatJ])


  list(ROC=ROC,AUC=AUC,Youden=hatJ,cutoff=hatC)

}


############################## Box-Cox method


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
  ROC=cbind(ss,1-pnorm(qnorm(1-ss,mx,sx),my,sy))
  AUC=mean(ROC[,2])
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


############################## Zhou XH, Lin H, 2008 (ZL)



datafreq1=function(x,y)
{
  Tvalue=c(x,y)
  group=c(rep(0,length(x)),rep(1,length(y)))

  tt=sort(unique(Tvalue))
  slr=tt
  skr=tt
  for(i in 1:length(tt))
  {
    slr[i]=sum((Tvalue==tt[i])&(group==1))
    skr[i]=sum((Tvalue==tt[i])&(group==0))
  }

  sd=as.numeric(skr>0)+as.numeric((skr>0)&(slr>0))
  dsd=diff(sd)
  ind=c( ((dsd==0)&(sd[-1]<=1)),FALSE)

  fy=c(min(tt)-1,tt[!ind])

  ni=length(fy)-1

  kr=rep(0,ni)
  lr=rep(0,ni)

  for(i in 1:ni)
  {
    kr[i]=sum((x>fy[i])&(x<=fy[i+1]))
    lr[i]=sum((y>fy[i])&(y<=fy[i+1]))
  }

  cbind(kr,lr)
}


loglik=function(x,y,alp,cvec)
{
  alp0=alp[1]
  alp1=alp[2]
  cvec1=c(-1000,cvec)
  cvec2=c(cvec,1000)

  out=datafreq1(x,y)
  kr=out[,1]
  lr=out[,2]

  part1=kr*log(pnorm(cvec2)-pnorm(cvec1)+1e-10)
  part2=lr*log(pnorm(-alp0+alp1*cvec2)-pnorm(-alp0+alp1*cvec1)+1e-10)

  -(sum(part1)+sum(part2))
}

loglik2=function(x,y,theta)
{

  out=datafreq1(x,y)
  kr=out[,1]
  lr=out[,2]

  alp=theta[1:2]
  alp0=alp[1]
  alp1=alp[2]

  cvec=cumsum(theta[-(1:2)])
  cvec1=c(-1000,cvec)
  cvec2=c(cvec,1000)

  part1=kr*log(pnorm(cvec2)-pnorm(cvec1)+1e-10)
  part2=lr*log(pnorm(-alp0+alp1*cvec2)-pnorm(-alp0+alp1*cvec1)+1e-10)

  -(sum(part1)+sum(part2))
}


maxite=function(x,y,theta)
{
  alp=theta[1:2]
  cvec=theta[-(1:2)]
  ni1=length(cvec)

  for(i in 1:ni1)
  {
    Extcvec=c(-10,cvec,10)
    cveci=cvec[-i]
    pln1=function(Ci)
    {

      newcvec=sort(c(cveci,Ci))
      out=loglik(x,y,alp,newcvec)
      out
    }
    output=nlminb(start=cvec[i],pln1,lower=c(Extcvec[i]),upper=c(Extcvec[i+2]))
    cvec[i]=output$par
  }

  pln2=function(alpha)
  {
    loglik(x,y,alpha,cvec)
  }
  output=nlminb(alp,pln2,lower=c(-10,0),upper=c(10,10))
  list(new=c(output$par,cvec),lik=output$objective)
}


newmax=function(x,y,theta)
{
  ni=nrow(datafreq1(x,y))
  ite=0
  out=maxite(x,y,theta)
  oldtheta=out$new
  oldlik=out$lik
  err=1
  newtheta=c(out$new[1:3],diff(out$new[-(1:2)]))


  while(err>1e-4)
  {
    output=optim(par=newtheta,fn=loglik2,
                 method="L-BFGS-B",x=x,y=y,lower=c(-10,0,-20,rep(0,ni-2)),
                 upper=c(10,10,rep(20,ni-1)),control=list(maxit=30))
    newtheta=output$par
    newlik=output$value
    err=abs(newlik-oldlik)
    oldlik=newlik
    ite=ite+1
  }
  list(new=newtheta,lik=newlik,ite=ite)
}


mle3=function(x,y)
{

  ni=nrow(datafreq1(x,y))
  out=newmax(x,y,c(1,1,sort(runif(ni-1,-5,5)) ))

  newtheta=c(out$new[1:2],cumsum(out$new[-(1:2)]))

  list(alp=newtheta,cvec=newtheta[-c(1:2)])
}


estzl=function(x,y,nss=10^4)
{
  oxy=mle3(x,y)$alp
  alp=oxy[c(1,2)]
  cvec=oxy[-c(1,2)]

  Tvalue=c(x,y)
  group=c(rep(0,length(x)),rep(1,length(y)))

  tt=sort(unique(Tvalue))
  slr=tt
  skr=tt
  for(i in 1:length(tt))
  {
    slr[i]=sum((Tvalue==tt[i])&(group==1))
    skr[i]=sum((Tvalue==tt[i])&(group==0))
  }

  sd=as.numeric(skr>0)+as.numeric((skr>0)&(slr>0))
  dsd=diff(sd)
  ind=c( ((dsd==0)&(sd[-1]<=1)),FALSE)

  index=c(1:length(tt))[!ind]
  fy=tt[index]
  status=sd[index]

  alp0=alp[1]
  alp1=alp[2]
  hatc=tt
  F0cvec=c(0,pnorm(cvec),1)
  F1cvec=c(0,pnorm(-alp0+alp1*cvec),1)

  for(i in 1:length(fy))
  {
    if(status[i]==2)
    {
      hatc[index[i]]=cvec[i]
    }

    if(status[i]==1)
    {
      low=ifelse(i==1,1,index[i-1]+1)
      upp=index[i]

      tot=F0cvec[i+1]-F0cvec[i]
      prob=skr[low:upp]*tot/sum(skr[low:upp])
      hatc[low:upp]=qnorm(F0cvec[i]+cumsum(prob))
    }

    if(status[i]==0)
    {
      low=ifelse(i==1,1,index[i-1]+1)
      upp=index[i]

      tot=F1cvec[i+1]-F1cvec[i]
      prob=slr[low:upp]*tot/sum(slr[low:upp])
      hatc[low:upp]=(qnorm(F1cvec[i]+cumsum(prob))+alp0)/alp1
    }


  }
  hatF0= pnorm(hatc)
  hatF1=pnorm(-alp0+alp1*hatc)

  roc=function(u)
  {
    pnorm(alp[1]+alp[2]*qnorm(u))-u
  }

  opt.res=optimize(roc,c(0.001,0.999),maximum=T)
  youden=opt.res$objective
  uHat=qnorm(1-opt.res$maximum)
  cutoff=tt[which.min(abs(uHat-hatc))]

  ss=((1:nss)-0.5)/nss
  ROC=cbind(ss,pnorm(alp[1]+alp[2]*qnorm(ss)))
  auc=mean(pnorm(alp[1]+alp[2]*qnorm(ss)))

  list(ROC=ROC,AUC=auc,Youden=youden,cutoff=cutoff)
}



############################## Empirical CDF (ECDF)


estemp= function(x,y,nss=10^4)
{
  hatF0 = as.stepfun(ecdf(x))
  hatF1 = as.stepfun(ecdf(y))

  tt = sort(unique(c(x,y)) )


  estF0=hatF0(tt)
  estF1=hatF1(tt)
  ss=((1:nss)-0.5)/nss

  ROC=cbind(ss,as.numeric(sapply(ss,disroc,hatF0=estF0,hatF1=estF1)))
  AUC=mean(ROC[,2])
  n=length(tt)

  dif = hatF0(tt)-hatF1(tt)
  hatJ = max(dif)

  index = c(1:length(tt))[abs(dif-hatJ)<1/n^2 ]

  hatC = mean(tt[index])

  list(ROC=ROC,AUC=AUC,Youden=hatJ,cutoff=hatC)
}



############################## MNLE method

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
  ROC=cbind(ss,as.numeric(sapply(ss,disroc,hatF0=estF0,hatF1=estF1)))
  AUC=mean(ROC[,2])
  hatF1=cumsum(hatp1)
  hatF0=cumsum(hatp0)

  dif = hatF0-hatF1
  hatJ = max(dif)

  index = c(1: length(tt) )[abs(dif-hatJ)<1/n^2 ]

  hatC = mean(tt[index])

  list(ROC=ROC,AUC=AUC,Youden=hatJ,cutoff=hatC)
}

############################## Liu et al., 2011 (LZL)



estsinica=function(x,y,nss=10^4){
  tt=sort(unique(c(x,y)))
  n0=length(x)
  n1=length(y)
  vec=rep(0,n1)
  ROC=function(u){
    for(i in 1:n1){
      vec[i]=(mean(x>=y[i])>=u)
    }
    1-mean(vec)
  }
  ss=((1:nss)-0.5)/nss
  roc=cbind(ss,sapply(ss,FUN = ROC))
  AUC=mean(roc[,2])
  hatJ=max(roc[,2]-ss)
  hatU=mean(ss[(roc[,2]-ss)>=hatJ])
  hatF0=as.stepfun(ecdf(x))
  hatC=tt[which.min(abs(hatU+hatF0(tt)-1))]
  list(ROC=roc,AUC=AUC,Youden=hatJ,cutoff=hatC)
}


############################## MSLE method

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
  ROC=cbind(ss,as.numeric(sapply(ss,disroc,hatF0=estF0,hatF1=estF1)))
  AUC=mean(ROC[,2])

  hatdif=abs(lcmfit$slope.knots-lam)

  len=length(hatdif)
  index=(1:len)[hatdif<=min(hatdif)]

  ind=index[1]
  hatW=lcmfit$x.knots[ind]
  hatF1=lcmfit$y.knots[ind]/lam
  hatJ=(hatW-hatF1)/(1-lam)
  indC=(1:length(tt))[estW==hatW]
  hatC=mean(tt[indC])

  list(ROC=ROC,AUC=AUC,Youden=hatJ,cutoff=hatC)

}


############################## Kernel method

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
  ROC=cbind(ss,as.numeric(sapply(ss,disroc,hatF0=estF0,hatF1=estF1)))
  AUC=mean(ROC[,2])
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
