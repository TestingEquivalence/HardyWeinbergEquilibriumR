source("distance.R")

tests=c(1,2,3,4)
names(tests)=c("asymptotic_test_conditional",
                 "resampling_test_conditional",
                 "asymptotic_test_minimum",
                 "resampling_test_minimum")

asympt_stdev<-function(p,derivative){
  vec = derivative
  vnsq_1  = sum(p*vec*vec)
  
  k=length(p)
  vnsq_2=0
  for (j1 in 1:k)
    for (j2 in 1:k)
      vnsq_2 = vnsq_2 + vec[j1] * vec[j2] * p[j1] * p[j2]
  
  
  vnsq  = vnsq_1 - vnsq_2
  return (sqrt(vnsq))
}


#' The asymptotic test is based on the asymptotic distribution of the test statistic. 
#' The test statistic is scaled Euclidian distance between the contingency table
#' and the product measure of the marginal distributions.
#' The asymptotic test needs some sufficiently large number of the observations
#' in any cell of the contingency table.
#' It should be used carefully because the test is approximate 
#' and may be anti-conservative at some points. 
#' In order to obtain a conservative test reducing of alpha  (usually halving) or
#' slight shrinkage of the tolerance parameter epsilon may be appropriate.
#' We prefer the slight shrinkage of the tolerance parameter 
#' because it is more effective and the significance level remains unchanged.
#' \code{asymptotic_test_absolute} asymptotic test for approximate row column independence
#' in two way contingency tables. 
#' The test statistic is scaled Euclidian distance between the contingency table
#' and the product measure of the marginal distributions.
#' @param tab contingency table containing the counts of events
#' @param alpha significance level
#' @return test returns the minimum tolerance parameter epsilon,
#' for which the approximate independence can be shown

asymptotic_test_conditional<-function(tab, alpha){
  n=sum(tab)
  tab=tab/n
  
  vtab=as.vector(t(tab))
  der=derivative_cond_l22(tab)
  vder=as.vector(t(der))
  
  vol = asympt_stdev(vtab,vder) / sqrt(n)
  qt=qnorm(1-alpha,0,1)
  t= cond_l22(tab)
  eps = t + qt*vol
  eps=sqrt(eps)
  return(eps)
}

asymptotic_test_minimum<-function(tab, alpha){
  #normalize tab
  n=sum(tab)
  tab=tab/n
  
  #calculate minimum distance
  res=min_l22(tab)
  prod=product(res$par)
  q=triangle(prod)
  
  
  vtab=as.vector(t(tab))
  der=l22_first_derivative(tab, q)
  vder=as.vector(t(der))
  
  vol = asympt_stdev(vtab,vder) / sqrt(n)
  qt=qnorm(1-alpha,0,1)
  t= res$value
  eps = t + qt*vol
  eps=sqrt(eps)
  return(eps)
}

resampling_stdev<-function(p,T,n, nSimulation){
  i=c(1:nSimulation)
  f<-function(k){
    vp=as.vector(p)
    v=rmultinom(n=1,size=n,prob=vp)
    v=v/sum(v)
    m=matrix(dat=v,nrow=nrow(p), ncol=ncol(p))
    return(T(m))
  }
  sample=sapply(i,f)
  return(sqrt(var(sample)))
}

resampling_test_conditional<-function(tab, alpha, nSimulation=2000){
  n=sum(tab)
  tab=tab/n
  set.seed(01012020)
  
  vol =resampling_stdev(p=tab,T=cond_l2,nSimulation = nSimulation,n=n)
  qt=qnorm(1-alpha,0,1)
  t= cond_l2(tab)
  eps = t + qt*vol
  return(eps)
}

resampling_test_minimum<-function(tab, alpha, nSimulation=2000){
  n=sum(tab)
  tab=tab/n
  set.seed(01012020)
  
  vol =resampling_stdev(p=tab,T=min_l2,nSimulation = nSimulation,n=n)
  qt=qnorm(1-alpha,0,1)
  t= min_l2(tab)
  eps = t + qt*vol
  return(eps)
}