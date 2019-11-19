source("tests.R")
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



powerAtBoundaryConditionalAsympt<-function(tab, nSamples, cl,alpha, scaleFactor){
  i=c(1:100)
  set.seed(01082019)
  n=sum(tab)
  alpha=alpha*1
  
  #generate close boundary points
  eps=cond_l2(tab/n)
  boundaryPoints=lapply(i, closeBoundaryPoint,tab,eps,cond_l2)
  
  #power of the asymptotic test for conditional distance
  eps=eps* scaleFactor
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

powerAtBoundaryMinAsympt<-function(tab, nSamples, cl,alpha, scaleFactor){
  i=c(1:100)
  set.seed(01082019)
  n=sum(tab)
  
  #generate close boundary points
  eps=min_l2(tab/n)
  boundaryPoints=lapply(i, closeBoundaryPoint,tab,eps,min_l2)
  
  #power of the asymptotic test for the minimum distance
  eps=eps*scaleFactor
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
 
powerAtBoundaryConditionalBst<-function(tab, nSamples, nSim,cl,alpha, scaleFactor){
  i=c(1:100)
  set.seed(01082019)
  n=sum(tab)
  nSim=nSim/1
  
  #generate close boundary points
  eps=cond_l2(tab/n)
  boundaryPoints=lapply(i, closeBoundaryPoint,tab,eps,cond_l2)
  
  
  j=1
  res=rep(NA,100)
  
  eps=eps* scaleFactor
  for (point in boundaryPoints){
    #power of the bootstrap test for conditional distance
    test<-function(point){
      res=bootstrap_test_conditional(tab = point,alpha=alpha,eps = eps,nSimulation =nSim)
      return(res$result)
    }
    
    res[j]=power(point,test,n,nSamples,cl)
    print(paste(j, " point done"))
    j=j+1
  }
  
  return(res)
}
 
powerAtBoundaryMinBst<-function(tab, nSamples, nSim,cl,alpha, scaleFactor){
  i=c(1:100)
  set.seed(01082019)
  n=sum(tab)
  nSim=nSim/1
  
  #generate close boundary points
  eps=min_l2(tab/n)
  boundaryPoints=lapply(i, closeBoundaryPoint,tab,eps,min_l2)
  
  #power of the bootstrap test for minimum distance
  eps=eps* scaleFactor
  test<-function(point){
    res=bootstrap_test_minimum(tab = point,alpha=alpha,eps = eps,nSimulation =nSim)
    return(res$result)
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

powerAtBoundaryConditionalResampling<-function(tab, nSamples, cl,alpha, scaleFactor){
  i=c(1:100)
  set.seed(01082019)
  n=sum(tab)
  alpha=alpha*1
  
  #generate close boundary points
  eps=cond_l2(tab/n)
  boundaryPoints=lapply(i, closeBoundaryPoint,tab,eps,cond_l2)
  
  #power of the asymptotic test for conditional distance
  eps=eps* scaleFactor
  test<-function(point){
    minEps=resampling_test_conditional(point,alpha)
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

powerAtBoundaryMinResampling<-function(tab, nSamples, cl,alpha, scaleFactor){
  i=c(1:100)
  set.seed(01082019)
  n=sum(tab)
  alpha=alpha*1
  
  #generate close boundary points
  eps=min_l2(tab/n)
  boundaryPoints=lapply(i, closeBoundaryPoint,tab,eps,min_l2)
  
  #power of the asymptotic test for conditional distance
  eps=eps* scaleFactor
  test<-function(point){
    minEps=resampling_test_minimum(point,alpha)
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



boundaryPower<-function(tab, nSamples, selector,nSimulation,alpha, scaleFactor){
  asy_cond=rep(NA,100)
  bst_cond=rep(NA,100)
  res_cond=rep(NA,100)
  asy_min=rep(NA,100)
  bst_min=rep(NA,100)
  res_min=rep(NA,100)
  
  
  cl=getCluster()
  
  
  if (selector[1]){
    asy_cond=powerAtBoundaryConditionalAsympt(tab,nSamples,cl,alpha, scaleFactor)
    print("asympt. cond. done!")
  }
  
  if (selector[2]){
    bst_cond=powerAtBoundaryConditionalBst(tab,nSamples,nSimulation,cl,alpha, scaleFactor)
    print("bst. cond. done!")
  }
  
  if (selector[3]){
    asy_min=powerAtBoundaryMinAsympt(tab,nSamples,cl,alpha, scaleFactor)
    print("asympt. min. done!")
  }
  
  if (selector[4]){
    bst_min=powerAtBoundaryMinBst(tab,nSamples,nSimulation,cl,alpha, scaleFactor)
    print("bst. min. done!")
  }
  
  if (selector[5]){
    res_cond=powerAtBoundaryConditionalResampling(tab,nSamples,cl,alpha, scaleFactor)
    print("res. cond. done!")
  }
  
  if (selector[6]){
    res_min=powerAtBoundaryMinResampling(tab,nSamples,cl,alpha,scaleFactor)
    print("res. min. done!")
  }
  
  stopCluster(cl)
  
  df=data.frame(asy_cond, bst_cond, res_cond, asy_min, bst_min,res_min)
  return(df)
}





