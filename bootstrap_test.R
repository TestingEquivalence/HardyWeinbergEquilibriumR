source("distance.R")
source("asymptotic_test.R")

randomPoint<-function(tab){
  #erzeuge hier diagonal matrix 
  size=nrow(tab)*ncol(tab)
  x=runif(size,0,1)
  x=x/sum(x)
  m=matrix(data=x,nrow=nrow(tab), ncol=ncol(tab),byrow = TRUE)
  m=triangle(m)
  return(m)
}

randomExteriorPoint<-function(tab, eps, distance){
  repeat{
    m = randomPoint(tab)
    t= distance(m)
    if (t>eps) return(m)
  }
}

linComb<-function(x,y,a){
  return((1-a)*x+a*y) 
}

linearBoundaryPoint<-function(p,q,eps,distance){
  aim<-function(a){
    lc=linComb(p,q,a)
    dst=distance(lc)
    return(dst-eps)
  }
  
  aMin=uniroot(aim, c(0,1))
  return(linComb(p,q,aMin$root))
}



protoBstTest<-function(tab,n,distance,eps,exteriorPoints,nSimulation){
  #calculate test statistic
  t=distance(tab)
  
  #estimate closest boundary point
  rp=tab
  
  df<-function(x){
    r=x-tab
    rv=as.vector(r)
    srv=rv*rv
    sse=sum(srv)
    return(sqrt(sse))
  }
  
  if (t<eps){
    #function(p,q,eps,distance)
    bps=lapply(X=exteriorPoints,FUN=linearBoundaryPoint,
               q=tab,eps=eps, distance=distance)
    dst=lapply(bps, df)
    pos=which.min(dst)
    rp=bps[[pos]]
  }
  
  
  #simulate bootstrap sample
  i=c(1:nSimulation)
  f<-function(k){
    vrp=as.vector(rp)
    v=rmultinom(n=1,size=n,prob=vrp)
    v=v/sum(v)
    m=matrix(dat=v,nrow=nrow(tab), ncol=ncol(tab))
    dst=distance(m)
    return(dst)
  }
  sample=lapply(i,f)
  
  #bootstrap test
  pValue=sum(sample<t)/nSimulation
  return(pValue)
}

#' The bootstrap test is based on the re-sampling method called bootstrap.
#' The bootstrap test is more precise and reliable than the asymptotic test.
#' However, it should be used carefully because the test is approximate 
#' and may be anti-conservative. 
#' In order to obtain a conservative test reducing of alpha
#' (usually halving) or slight shrinkage of the tolerance parameter epsilon
#' may be appropriate. We prefer the slight shrinkage of the tolerance parameter 
#' because it is more effective and the significance level remains unchanged.
#'  
#' The test statistic is scaled Euclidian distance between the contingency table
#' and the product measure of the marginal distributions.
#' \code{bootstrap_test_absolute} bootstrap test for approximate row column independence
#' in two way contingency tables. 
#' The test statistic is scaled Euclidian distance between the counting frequencies
#' and the product measure of the marginal distributions.
#' @param tab contingency table containing the counts of events
#' @param alpha significance level
#' @param nSimulation number of bootstrap samples, default 10000 
#' @param nExteriorPoints number of random directions to search for a boundary point,
#' default is (nrow(tab)+ncol(tab))*50
#' @return test returns the minimum tolerance parameter epsilon,
#' for which the approximate independence can be shown

bootstrap_test_conditional<-function(tab, alpha, 
                         nSimulation=10000, 
                         nExteriorPoints=0,
                         eps=0){
  #find start value for min eps
  #use for this purpose the asymptotic test with 
  #small safety margin
  #or set it eps if provided
  beps=eps
  if (eps==0)
    beps=asymptotic_test_conditional(tab,alpha)*1.1
  
  
  
  n=sum(tab)
  tab=tab/n
  
  #number of search directions and seed
  if (nExteriorPoints==0) 
    nExteriorPoints=nrow(tab)*50*4
  
  set.seed(10071977)
  
  distance<-function(x){
    dst=cond_l22(x)
    return(sqrt(dst))
  }
  
  #calculate exterior points
  f<-function(x){
    randomExteriorPoint(tab,beps,distance)
  }
  
  i=c(1:nExteriorPoints)
  exteriorPoints=lapply(i, f)
  
  # if epsilon given, make short-cut and
  # calculate the logical value only
  if (eps>0){
    pval=protoBstTest(tab,n,distance,eps,exteriorPoints,nSimulation)
    l=pval<=alpha
    ls=list(p_value=pval,result=l)
    return(ls)
  }
  
  
  #calculate min epsilon
  ff<-function(x){
    set.seed(01012019)
    pval=protoBstTest(tab,n,distance,eps = x,exteriorPoints,nSimulation)
    pval-alpha
  }
  
  #check boundary values
  #check lower bound
  lb=ff(0)
  if (lb<0) return(0)
  
  #check upper bound
  ub=ff(beps)
  if (ub>0) return(NA)
  
  res=uniroot(ff,c(0,beps))
  return(res$root)
}

bootstrap_test_minimum<-function(tab, alpha, 
                                 nSimulation=10000, 
                                 nExteriorPoints=0,
                                 eps=0){
  #find start value for min eps
  #use for this purpose the asymptotic test with 
  #small safety margin
  #or set it eps if provided
  beps=eps
  if (eps==0)
    beps=asymptotic_test_minimum(tab,alpha)*1.1
  
  n=sum(tab)
  tab=tab/n
  
  #number of search directions and seed
  if (nExteriorPoints==0) 
    nExteriorPoints=nrow(tab)*50*4
  
  set.seed(10071977)
  
  distance<-function(x){
    res=min_l22(x)
    return(sqrt(res$val))
  }
  
  #calculate exterior points
  f<-function(x){
    randomExteriorPoint(tab,beps,distance)
  }
  
  i=c(1:nExteriorPoints)
  exteriorPoints=lapply(i, f)
  
  # if epsilon given, make short-cut and
  # calculate the logical value only
  if (eps>0){
    pval=protoBstTest(tab,n,distance,eps,exteriorPoints,nSimulation)
    l=pval<=alpha
    return(list(p_value=pval,result=l))
  }
  
  #calculate min epsilon
  ff<-function(x){
    set.seed(01012019)
    pval=protoBstTest(tab,n,distance,eps = x,exteriorPoints,nSimulation)
    pval-alpha
  }
  
  #check boundary values
  #check lower bound
  lb=ff(0)
  if (lb<0) return(0)
  
  #check upper bound
  ub=ff(beps)
  if (ub>0) return(NA)
  
  res=uniroot(ff,c(0,beps))
  return(res$root)
}


