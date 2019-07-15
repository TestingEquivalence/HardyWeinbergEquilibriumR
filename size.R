source("distance.R")
source("asymptotic_test.R")
source("bootstrap_test.R")
source("data_sets.R")
source("simulation.R")

powerAtHWE<-function(p,n,eps,nSamples,selector){
  alpha=0.05
  hwe=p2hwe(p)
  
  #power of asymptotic test, conditional distance
  p1=0
  
  if (selector[1]){
    test<-function(tab){
      minEps=asymptotic_test_conditional(tab,alpha)
      return(minEps<=eps)
    }
    
    set.seed(11072019)
    p1=power(test,n,hwe,nSamples)
  }
  
  #power of bootstrap test, conditional distance
  p2=0
  
  if (selector[2]){
    test<-function(tab){
      res=bootstrap_test_conditional(tab,alpha, eps=eps)
      return(res$result)
    }
    
    set.seed(11072019)
    p2=power(test,n,hwe,nSamples)
  }
  
  #power of asymptotic test, minimum distance
  p3=0
  
  if (selector[3]){
    test<-function(tab){
      minEps=asymptotic_test_minimum(tab,alpha)
      return(minEps<=eps)
    }
    
    set.seed(11072019)
    p3=power(test,n,hwe,nSamples)
  }
  
  #power of bootstrap test, conditional distance
  p4=0
  
  if (selector[4]){
    test<-function(tab){
      res=bootstrap_test_minimum(tab,alpha, eps=eps)
      return(res$result)
    }
    
    set.seed(11072019)
    p4=power(test,n,hwe,nSamples)
  }
  
  
  return(c(eps,p1,p2,p3,p4))
}

# power of sample
powerAtPoint<-function(tab, eps, nSamples, selector){
  n=sum(tab)
  p=startValue(tab/n)
  powerAtHWE(p,n,eps,nSamples,selector)
}


rp<-function(i,p,n){
  rp=rmultinom(n=1,size=n,prob=p)
  rp=rp/sum(rp)
}


# power sensitivity
powerSensitivity<-function(tab, eps, nSamples,selector){
  n=sum(tab)
  p=startValue(tab/n)
  
  nPoints=100
  i=c(1:nPoints)
  
  set.seed(01082019)
  points=lapply(i,rp,p,n)

  sapply(points, powerAtHWE,eps=eps,nSamples=nSamples,selector=selector, n=n)
}

# power at sample1
#sizeSample1=matrix(data=NA, nrow=10, ncol=5)
#colnames(sizeSample1)=c("eps","asy_cond","bst_cnd","asy_min","bst_min")

#sizeSample1[1,]=powerAtPoint(example1,0.12,1000,c(TRUE,TRUE,TRUE,TRUE))
#sizeSample1[2,]=powerAtPoint(example1,0.10,1000,c(TRUE,TRUE,TRUE,TRUE)) 
#sizeSample1[3,]=powerAtPoint(example1,0.09,1000,c(TRUE,TRUE,TRUE,TRUE)) 
#sizeSample1[4,]=powerAtPoint(example1,0.08,1000,c(TRUE,TRUE,TRUE,TRUE)) 
#sizeSample1[5,]=powerAtPoint(example1,0.07,1000,c(TRUE,TRUE,TRUE,TRUE)) 

#write.table(sizeSample1, "sizeSample1.txt")

# power at sample2
#sizeSample2=matrix(data=NA, nrow=10, ncol=5)
#colnames(sizeSample2)=c("eps","asy_cond","bst_cnd","asy_min","bst_min")

#sizeSample2[1,]=powerAtPoint(example2,0.12,1000,c(TRUE,TRUE,TRUE,TRUE))
#sizeSample2[2,]=powerAtPoint(example2,0.1,1000,c(TRUE,TRUE,TRUE,TRUE))
#sizeSample2[3,]=powerAtPoint(example2,0.09,1000,c(TRUE,TRUE,TRUE,TRUE))
#sizeSample2[4,]=powerAtPoint(example2,0.08,1000,c(TRUE,TRUE,TRUE,TRUE))
#sizeSample2[5,]=powerAtPoint(example2,0.07,1000,c(TRUE,TRUE,TRUE,TRUE))

#write.table(sizeSample2, "sizeSample2.txt")

# power at sample3
#sizeSample3=matrix(data=NA, nrow=10, ncol=5)
#colnames(sizeSample3)=c("eps","asy_cond","bst_cnd","asy_min","bst_min")

# sizeSample3[1,]=powerAtPoint(example3,0.018,1000,c(TRUE,TRUE,TRUE,TRUE))
# sizeSample3[2,]=powerAtPoint(example3,0.016,1000,c(TRUE,TRUE,TRUE,TRUE))
# sizeSample3[3,]=powerAtPoint(example3,0.014,1000,c(TRUE,TRUE,TRUE,TRUE))
# sizeSample3[4,]=powerAtPoint(example3,0.012,1000,c(TRUE,TRUE,TRUE,TRUE))
# sizeSample3[5,]=powerAtPoint(example3,0.010,1000,c(TRUE,TRUE,TRUE,TRUE))
# 
# 
# write.table(sizeSample3, "sizeSample3.txt")


# power at sample4
# sizeSample4=matrix(data=NA, nrow=10, ncol=5)
# colnames(sizeSample4)=c("eps","asy_cond","bst_cnd","asy_min","bst_min")
# 
# sizeSample4[1,]=powerAtPoint(example4,0.06,1000,c(TRUE,TRUE,TRUE,TRUE))
# sizeSample4[2,]=powerAtPoint(example4,0.05,1000,c(TRUE,TRUE,TRUE,TRUE))
# sizeSample4[3,]=powerAtPoint(example4,0.04,1000,c(TRUE,TRUE,TRUE,TRUE))
# sizeSample4[4,]=powerAtPoint(example4,0.03,1000,c(TRUE,TRUE,TRUE,TRUE))
# 
# 
# write.table(sizeSample4, "sizeSample4.txt")

sensi_example1=powerSensitivity(tab=example1,eps=0.1,nSamples = 1000, selector = c(TRUE,FALSE,FALSE,FALSE))
sensi_example2=powerSensitivity(tab=example1,eps=0.1,nSamples = 1000, selector = c(TRUE,FALSE,FALSE,FALSE))
