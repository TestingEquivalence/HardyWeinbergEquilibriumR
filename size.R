source("distance.R")
source("asymptotic_test.R")
source("bootstrap_test.R")
source("data_sets.R")

rp<-function(i,p,n){
  rp=rmultinom(n=1,size=n,prob=p)
  rp=rp/sum(rp)
  res=triangle(product(rp))
  return(res)
}

size<-function(table, test){
  n=sum(tab)
  tab=tab/n
  
  p=startValue(tab)
  
  hwe_points=c(triangle(product(p)))
  i=c(1:100)
  set.seed(10072019)
  vec=lapply(i, rp, n=n, p=p)
  hwe_points=c(hwe_points,vec)
  
  lapply(hwe_points, test)
}

alpha=0.05
test_asymp_cond_dst<-function(tab){
  asymptotic_test_conditional(tab,alpha)
}

size_asymp_cond_dst=size(  
