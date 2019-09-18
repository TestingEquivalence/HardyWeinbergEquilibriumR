source("distance.R")
source("asymptotic_test.R")

protoBstTest<-function(tab,n,distance,eps,nSimulation){
  #calculate test statistic
  #t=distance(tab)
  
  #simulate bootstrap sample
  i=c(1:nSimulation)
  f<-function(k){
    vrp=as.vector(tab)
    v=rmultinom(n=1,size=n,prob=vrp)
    v=v/sum(v)
    m=matrix(dat=v,nrow=nrow(tab), ncol=ncol(tab))
    dst=distance(m)
    return(dst)
  }
  sample=lapply(i,f)
  
  #bootstrap test
  pValue=sum(sample>eps)/nSimulation
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
                         eps=0){
  #find start value for min eps
  #use for this purpose the asymptotic test with 
  #small safety margin
  #or set it eps if provided
  beps=eps
  if (eps==0)
    beps=asymptotic_test_conditional(tab,alpha)*1.2
  
  
  
  n=sum(tab)
  tab=tab/n
  
  set.seed(10071977)
  
  distance<-function(x){
    cond_l2(x)
  }
  
  # if epsilon given, make short-cut and
  # calculate the logical value only
  if (eps>0){
    pval=protoBstTest(tab,n,distance,eps,nSimulation)
    l=pval<=alpha
    ls=list(p_value=pval,result=l)
    return(ls)
  }
  
  
  #calculate min epsilon
  ff<-function(x){
    set.seed(01012019)
    pval=protoBstTest(tab,n,distance,eps = x,nSimulation)
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
  
  set.seed(10071977)
  
  distance<-function(x){
    res=min_l22(x)
    return(sqrt(res$val))
  }
  
  # if epsilon given, make short-cut and
  # calculate the logical value only
  if (eps>0){
    pval=protoBstTest(tab,n,distance,eps,nSimulation)
    l=pval<=alpha
    return(list(p_value=pval,result=l))
  }
  
  #calculate min epsilon
  ff<-function(x){
    set.seed(01012019)
    pval=protoBstTest(tab,n,distance,eps = x,nSimulation)
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


