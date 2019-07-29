library(parallel)
source("distance.R")
source("asymptotic_test.R")
source("bootstrap_test.R")


# Calculate the number of cores
getCluster<-function(){
  no_cores <- detectCores() - 1

  # Initiate cluster
  cl <- makeCluster(no_cores,'SOCK')
  clusterExport(cl,c("asymptotic_test_conditional","cond_l22_da","asympt_stdev","cond_l22", "l22",
                   "derivative_cond_l22", "bootstrap_test_conditional", "randomExteriorPoint",
                   "startValue","triangle","product","cond_l22_db","randomPoint","protoBstTest",
                   "linearBoundaryPoint", "linComb","asymptotic_test_minimum","min_l22","fn",
                   "p2triangle","l22_first_derivative","bootstrap_test_minimum"))
                
  return(cl)
}

sample<-function(i,tab,n){
  vtab=as.vector(tab)
  v=rmultinom(n=1,size=n,prob=vtab)
  m=matrix(dat=v,nrow=nrow(tab), ncol=ncol(tab))
  return(m)
}

power<-function(tab, test, n,  nSamples, cl ){
  source("distance.R")
  source("asymptotic_test.R")
  source("bootstrap_test.R")
  i=c(1:nSamples)
  sampleList=lapply(i, sample, tab,n)
  v=parSapply(cl,sampleList, test)
  #v=sapply(sampleList, test)
  return(sum(v==TRUE)/nSamples)
}

