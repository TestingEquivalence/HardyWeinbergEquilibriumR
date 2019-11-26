source("distance.R")

#vector of test numbers and names
#test numbers are used in research.R to choose,
#which tests should be simulated
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
#' The test statistic is Euclidian distance between the genotype distribution
#' and Hardy Weinberg Equilibrium, which is implied by allele distribution.
#' The test should be used carefully because it is approximate 
#' and may be anti-conservative at some points. 
#' In order to obtain a conservative test reducing of alpha  (usually halving) or
#' slight shrinkage of the tolerance parameter epsilon may be appropriate.
#' \code{asymptotic_test_conditional} asymptotic equivalence test for Hardy Weinberg Equilibrium
#' @param tab genotype distribution (lower triangle)
#' @param alpha significance level
#' @return test returns the minimum tolerance parameter,
#' for which equivalence could be shown
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


#' The asymptotic test is based on the asymptotic distribution of the test statistic. 
#' The test statistic is the minimum Euclidean distance between the genotype distribution 
#' and the family of Hardy Weinberg Equilibriums.
#' The test should be used carefully because it is approximate 
#' and may be anti-conservative at some points. 
#' In order to obtain a conservative test reducing of alpha  (usually halving) or
#' slight shrinkage of the tolerance parameter may be appropriate.
#' \code{asymptotic_test_minimum} asymptotic equivalence test for Hardy Weinberg Equilibrium
#' @param tab genotype distribution (lower triangle)
#' @param alpha significance level
#' @return test returns the minimum tolerance parameter,
#' for which equivalence could be shown
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

#' The resampling test uses bootstrap to estimate the variance of the test statistics.
#' The test statistic is Euclidian distance between the genotype distribution
#' and Hardy Weinberg Equilibrium, which is implied by allele distribution.
#' The test should be used carefully because it is approximate 
#' and may be anti-conservative at some points. 
#' In order to obtain a conservative test reducing of alpha  (usually halving) or
#' slight shrinkage of the tolerance parameter epsilon may be appropriate.
#' \code{resampling_test_conditional} equivalence test for Hardy Weinberg Equilibrium
#' @param tab genotype distribution (lower triangle)
#' @param alpha significance level
#' @param nSimulation number of resampling simulationen
#' @return test returns the minimum tolerance parameter,
#' for which equivalence could be shown
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

#' The resampling test uses bootstrap to estimate the variance of the test statistics.
#' The test statistic is the minimum Euclidean distance between the genotype distribution 
#' and the family of Hardy Weinberg Equilibriums.
#' The test should be used carefully because it is approximate 
#' and may be anti-conservative at some points. 
#' In order to obtain a conservative test reducing of alpha  (usually halving) or
#' slight shrinkage of the tolerance parameter may be appropriate.
#' \code{resampling_test_minimum} equivalence test for Hardy Weinberg Equilibrium
#' @param tab genotype distribution (lower triangle)
#' @param alpha significance level
#' @param nSimulation number of resampling simulationen
#' @return test returns the minimum tolerance parameter,
#' for which equivalence could be shown

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