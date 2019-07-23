source("distance.R")
source("asymptotic_test.R")
source("bootstrap_test.R")
source("data_sets.R")
source("simulation.R")

closeRandomPoint<-function(tab, eps, distance){
  n=sum(tab)
  tab=tab/n
  
  repeat{
    rtab=as.vector(tab)
    v=rmultinom(n=1,size=n,prob=rtab)
    v=v/n
    m=matrix(data=v,nrow=nrow(tab), ncol=ncol(tab))
    t= distance(m)
    if (t>eps) return(m)
  }
  
}

closeBoundaryPoint<-function(i,tab,eps,distance){
  p=closeRandomPoint(tab,eps,distance)
  n=sum(tab)
  tab=tab/n
  
  if (identical(distance,cond_l2)){
    q=p2triangle(startValue(tab))
  }
  
  if (identical(distance,min_l2)){
    q=min_l22(tab)$par
    q=p2triangle(q)
  }
  
  res=linearBoundaryPoint(p,q,eps,distance)
  return(res)
}

powerAtBoundary<-function(tab, nSamples,distance,bootstrap){
  i=c(1:100)
  set.seed(01082019)
  alpha=0.05
  n=sum(tab)
  
  eps=cond_l2(tab/n)
  
  #generate close boundary points
  boundaryPoints=lapply(i, closeBoundaryPoint,tab,eps,distance)
  
  if (identical(distance,cond_l2)){
    if (bootstrap){
      #power of the bootstrap test for conditional distance
      test<-function(point){
        bootstrap_test_conditional(tab=point,alpha,eps)
      }
    }
    else{
      #power of the asymptotic test for conditional distance
      test<-function(point){
        minEps=asymptotic_test_conditional(point,alpha)
        return(minEps<=eps)
      }
    }
  }
  
  if (identical(distance,min_l2)){
    if (bootstrap){
      
    }
    else{
      
    }
    
  }
  
  
  sapply(boundaryPoints,power,test,n, nSamples)
}




boundaryPower<-function(tab, nSamples, selector){
  asy_cond=rep(NA,100)
  bst_cond=rep(NA,100)
  asy_min=rep(NA,100)
  bst_min=rep(NA,100)
  
  
  if (selector[1]){
    asy_cond=powerConditionalDistanceAsymptotic(tab,nSamples)
    print("asympt. cond. done!")
  }
  
  if (selector[2]){
    bst_cond=powerConditionalDistanceBootstrap(tab,nSamples)
    print("bst. cond. done!")
  }
  
  if (selector[3]){
    asy_min=powerMinimumDistanceAsymptotic(tab, nSamples)
    print("asympt. min. done!")
  }
  
  if (selector[4])
  bst_min=powerMinimumDistanceBootstrap(tab,nSamples)
  print("asympt. min. done!")
  
  df=data.frame(asy_cond, bst_cnd, asy_min, bst_min)
  return(df)
}

powerExample1=boundaryPower(example1,1000,c(TRUE,FALSE,FALSE,FALSE))


