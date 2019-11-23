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

boundaryPower<-function(tab, nSamples, testsToDo){
  alpha=0.05
  asy_cond=rep(NA,100)
  res_cond=rep(NA,100)
  asy_min=rep(NA,100)
  res_min=rep(NA,100)
  
  i=c(1:100)
  n=sum(tab)
  
  #generate close boundary points for min. distance
  set.seed(01082019)
  eps=min_l2(tab/n)
  boundaryPointsMin=lapply(i, closeBoundaryPoint,tab,eps,min_l2)
  
  #generate close boundary points for conditional distance
  set.seed(01082019)
  eps=cond_l2(tab/n)
  boundaryPointsCond=lapply(i, closeBoundaryPoint,tab,eps,cond_l2)
  
  
  
  cl=getCluster()
  
  if (1 %in% testsToDo){
    
    test<-function(tab){
      minEps=asymptotic_test_conditional(tab,alpha)
      return(minEps<=eps)
    }
    
    set.seed(11072019)
    asy_cond= sapply(boundaryPointsCond, power, test,n,nSamples,cl)
  }
  
  p2=NA
  if (2 %in% testsToDo){
    
    test<-function(tab){
      minEps=resampling_test_conditional(tab,alpha)
      return(minEps<=eps)
    }
    
    set.seed(11072019)
    res_cond= sapply(boundaryPointsCond, power, test,n,nSamples,cl)
  }
  
  #power of asymptotic test, minimum distance
  p3=NA
  if (3 %in% testsToDo){
    test<-function(tab){
      minEps=asymptotic_test_minimum(tab,alpha)
      return(minEps<=eps)
    }
    
    set.seed(11072019)
    asy_min= sapply(boundaryPointsMin, power, test,n,nSamples,cl)
  }
  
  p4=NA
  if (4 %in% testsToDo){
    test<-function(tab){
      minEps=resampling_test_minimum(tab,alpha)
      return(minEps<=eps)
    }
    
    set.seed(11072019)
    res_min= sapply(boundaryPointsMin, power, test,n,nSamples,cl)
  }
  
  stopCluster(cl)
  
  df=data.frame(asy_cond, res_cond, asy_min, res_min)
  return(df)
}





