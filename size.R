
powerAtHWE<-function(p,n,eps,nSamples, testsToDo, cl){
  alpha=0.05
  tab=p2hwe(p)
  
  #power of asymptotic test, conditional distance
  p1=NA
  
  if (1 %in% testsToDo){
    
    test<-function(tab){
      minEps=asymptotic_test_conditional(tab,alpha)
      return(minEps<=eps)
    }
    
    set.seed(11072019)
    p1=power(tab,test,n,nSamples,cl)
  }
  
  p2=NA
  if (2 %in% testsToDo){
    
    test<-function(tab){
      minEps=resampling_test_conditional(tab,alpha)
      return(minEps<=eps)
    }
    
    set.seed(11072019)
    p2=power(tab,test,n,nSamples,cl)
  }
  
  #power of asymptotic test, minimum distance
  p3=NA
  if (3 %in% testsToDo){
    test<-function(tab){
      minEps=asymptotic_test_minimum(tab,alpha)
      return(minEps<=eps)
    }
    
    set.seed(11072019)
    p3=power(tab,test,n,nSamples,cl)
  }
  
  p4=NA
  if (4 %in% testsToDo){
    test<-function(tab){
      minEps=resampling_test_minimum(tab,alpha)
      return(minEps<=eps)
    }
    
    set.seed(11072019)
    p4=power(tab,test,n,nSamples,cl)
  }
  
  return(c(eps,p1,p2,p3,p4))
}

# power of sample
powerAtPoint<-function(tab, eps, nSamples, testNumber){
  n=sum(tab)
  p=startValue(tab/n)
  cl=getCluster()
  res=powerAtHWE(p,n,eps,nSamples,testNumber,cl)
  stopCluster(cl)
  return(res)
}


rp<-function(i,p,n){
  rp=rmultinom(n=1,size=n,prob=p)
  rp=rp/sum(rp)
}


# power sensitivity
powerSensitivity<-function(tab, eps, nSamples,testsToDo){
  n=sum(tab)
  p=startValue(tab/n)
  
  nPoints=100
  i=c(1:nPoints)
  
  set.seed(01082019)
  points=lapply(i,rp,p,n)
  
  j=1
  res=matrix(data = NA, nrow=100, ncol=5)
  
  cl=getCluster()
  
  for (point in points){
    res[j,]=powerAtHWE(p=point,n,eps,nSamples,testsToDo,cl)
    print(paste(j, " point done"))
    j=j+1
  }

  stopCluster(cl)
  colnames(res)=c("eps","asy_cond","res_cond","asy_min","res_min")
  return(res)
}

