source("distance.R")
source("asymptotic_test.R")
source("bootstrap_test.R")
source("data_sets.R")
source("simulation.R")

rp<-function(i,p,n){
  rp=rmultinom(n=1,size=n,prob=p)
  rp=rp/sum(rp)
  res=triangle(product(rp))
  return(res)
}

powerAtHWE<-function(p,n,eps){
  alpha=0.05
  nSamples=10
  hwe=p2hwe(p)
  
  #power of asymptotic test, conditional distance
  test<-function(tab){
    minEps=asymptotic_test_conditional(tab,alpha)
    return(minEps<=eps)
  }

  set.seed(11072019)
  p1=power(test,n,hwe,nSamples)
  
  #power of bootstrap test, conditional distance
  test<-function(tab){
    res=bootstrap_test_conditional(tab,alpha, eps=eps)
    return(res$result)
  }
  
  set.seed(11072019)
  p2=power(test,n,hwe,nSamples)
  
  #power of asymptotic test, minimum distance
  test<-function(tab){
    minEps=asymptotic_test_minimum(tab,alpha)
    return(minEps<=eps)
  }
  
  set.seed(11072019)
  p3=power(test,n,hwe,nSamples)
  
  #power of bootstrap test, conditional distance
  test<-function(tab){
    res=bootstrap_test_minimum(tab,alpha, eps=eps)
    return(res$result)
  }
  
  set.seed(11072019)
  p4=power(test,n,hwe,nSamples)
  
  
  
  return(c(eps,p1,p2,p3,p4))
}
