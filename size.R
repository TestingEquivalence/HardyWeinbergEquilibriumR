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

pList<-function(tab){
  
}

alpha=0.05
eps=0.13

f<-function(i){
  min_eps=asymptotic_test_conditional(example1,alpha)
  return(min_eps<=eps)
}

power(f,1000)
