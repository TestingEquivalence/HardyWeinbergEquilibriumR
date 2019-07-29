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

powerAtBoundary<-function(tab, nSamples,distance,bootstrap, nSimulation,cl){
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
        bootstrap_test_conditional(tab=point,alpha,eps, nSimulation)
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
      
      #power of the bootstrap test for minimum distance
      test<-function(point){
        bootstrap_test_minimum(tab=point,alpha,eps, nSimulation)
      }
    }
    else{
      #power of the asymptotic test for minimum distance
      test<-function(point){
        minEps=asymptotic_test_minimum(point,alpha)
        return(minEps<=eps)
      }
    }
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
    asy_cond=powerAtBoundary(tab,nSamples,distance = cond_l2,bootstrap = FALSE, nSimulation, cl)
    print("asympt. cond. done!")
  }
  
  if (selector[2]){
    bst_cond=powerAtBoundary(tab,nSamples,distance = cond_l2,bootstrap = TRUE, nSimulation, cl)
    print("bst. cond. done!")
  }
  
  if (selector[3]){
    asy_min=powerAtBoundary(tab,nSamples,distance =min_l2,bootstrap = FALSE, nSimulation,cl)
    print("asympt. min. done!")
  }
  
  if (selector[4]){
    bst_min=powerAtBoundary(tab,nSamples,distance = min_l2,bootstrap = TRUE, nSimulation,cl)
    print("bst. min. done!")
  }
  
  stopCluster(cl)
  
  df=data.frame(asy_cond, bst_cond, asy_min, bst_min)
  return(df)
}

# powerExample1=boundaryPower(tab=example1, nSamples =  1000,
#                             selector = c(TRUE,FALSE,FALSE,FALSE), nSimulation =  1000)
# write.table(powerExample1,"powerExample1.txt")
# 
# powerExample2=boundaryPower(tab=example2, nSamples =  1000,
#                             selector=c(TRUE,FALSE,FALSE,FALSE),nSimulation =  1000)
# write.table(powerExample2,"powerExample2.txt")
# 
# powerExample3=boundaryPower(tab=example3, nSamples =  1000,
#                             selector =  c(TRUE,FALSE,FALSE,FALSE), nSimulation =  1000)
# write.table(powerExample3,"powerExample3.txt")
# 
# powerExample4=boundaryPower(tab=example4, nSamples =  1000,
#                             selector =  c(TRUE,FALSE,FALSE,FALSE), nSimulation = 1000)
# write.table(powerExample4,"powerExample4.txt")




