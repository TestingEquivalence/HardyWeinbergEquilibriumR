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

closeBoundaryPoint<-function(tab,eps,distance){
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

# closeRandomPoint(example1,0.1, cond_l2)
# closeBoundaryPoint(example1,0.1,cond_l2)
