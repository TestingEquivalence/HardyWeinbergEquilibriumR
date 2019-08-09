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

powerAtBoundaryConditionalAsympt<-function(tab, nSamples, cl){
  i=c(1:100)
  set.seed(01082019)
  alpha=0.05
  n=sum(tab)
  
  #generate close boundary points
  eps=cond_l2(tab/n)
  boundaryPoints=lapply(i, closeBoundaryPoint,tab,eps,cond_l2)
  
  #power of the asymptotic test for conditional distance
  test<-function(point){
    minEps=asymptotic_test_conditional(point,alpha)
    return(minEps<=eps)
  }
  
  j=1
  res=rep(NA,100)
  
  for (point in boundaryPoints){
    res[j]=power(point,test,n,nSamples,cl)
    print(paste(j, " point done"))
    j=j+1
  }
  
  return(res)
}

powerAtBoundaryMinAsympt<-function(tab, nSamples, cl){
  i=c(1:100)
  set.seed(01082019)
  alpha=0.05
  n=sum(tab)
  
  #generate close boundary points
  eps=min_l2(tab/n)
  boundaryPoints=lapply(i, closeBoundaryPoint,tab,eps,min_l2)
  
  #power of the asymptotic test for conditional distance
  test<-function(point){
    minEps=asymptotic_test_minimum(point,alpha)
    return(minEps<=eps)
  }
  
  j=1
  res=rep(NA,100)
  
  for (point in boundaryPoints){
    res[j]=power(point,test,n,nSamples,cl)
    print(paste(j, " point done"))
    j=j+1
  }
  
  return(res)
}
 

powerAtBoundaryConditionalBst<-function(tab, nSamples, nSim,cl){
  i=c(1:100)
  set.seed(01082019)
  alpha=0.05
  n=sum(tab)
  nSim=nSim/1
  
  #generate close boundary points
  eps=cond_l2(tab/n)
  boundaryPoints=lapply(i, closeBoundaryPoint,tab,eps,cond_l2)
  
  j=1
  res=rep(NA,100)
  
  for (point in boundaryPoints){
    #power of the asymptotic test for conditional distance
    test<-function(point){
      res=bootstrap_test_conditional(tab = point,alpha=alpha,eps = eps,nSimulation =nSim)
      return(res)
    }
    
    res[j]=power(point,test,n,nSamples,cl)
    print(paste(j, " point done"))
    j=j+1
  }
  
  return(res)
}
 
powerAtBoundaryMinBst<-function(tab, nSamples, nSim,cl){
  i=c(1:100)
  set.seed(01082019)
  alpha=0.05
  n=sum(tab)
  nSim=nSim/1
  
  #generate close boundary points
  eps=min_l2(tab/n)
  boundaryPoints=lapply(i, closeBoundaryPoint,tab,eps,min_l2)
  
  #power of the asymptotic test for conditional distance
  test<-function(point){
    res=bootstrap_test_minimum(tab = point,alpha=alpha,eps = eps,nSimulation =nSim)
    return(res)
  }
  
  j=1
  res=rep(NA,100)
  
  for (point in boundaryPoints){
    res[j]=power(point,test,n,nSamples,cl)
    print(paste(j, " point done"))
    j=j+1
  }
  
  return(res)
}



boundaryPower<-function(tab, nSamples, selector,nSimulation){
  asy_cond=rep(NA,100)
  bst_cond=rep(NA,100)
  asy_min=rep(NA,100)
  bst_min=rep(NA,100)
  
  cl=getCluster()
  
  
  if (selector[1]){
    asy_cond=powerAtBoundaryConditionalAsympt(tab,nSamples,cl)
    print("asympt. cond. done!")
  }
  
  if (selector[2]){
    bst_cond=powerAtBoundaryConditionalBst(tab,nSamples,nSimulation,cl)
    print("bst. cond. done!")
  }
  
  if (selector[3]){
    asy_min=powerAtBoundaryMinAsympt(tab,nSamples,cl)
    print("asympt. min. done!")
  }
  
  if (selector[4]){
    bst_min=powerAtBoundaryMinBst(tab,nSamples,nSimulation,cl)
    print("bst. min. done!")
  }
  
  stopCluster(cl)
  
  df=data.frame(asy_cond, bst_cond, asy_min, bst_min)
  return(df)
}





