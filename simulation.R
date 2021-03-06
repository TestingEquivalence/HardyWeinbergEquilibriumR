library(parallel)

# Calculate the number of cores
getCluster<-function(){
  no_cores <- detectCores() - 1

  # Initiate cluster
  cl <- makeCluster(no_cores,'SOCK')
  clusterExport(cl,c("asymptotic_test_conditional","cond_l22_da","asympt_stdev",
                     "cond_l22", "l22","startValue",
                   "derivative_cond_l22","triangle","product","cond_l22_db",
                   "asymptotic_test_minimum","min_l22","fn","cond_l2","min_l2",
                   "p2triangle","l22_first_derivative",
                   "closeRandomPoint",
                   "resampling_test_conditional","resampling_stdev","resampling_test_minimum"))
                
  return(cl)
}

sample<-function(i,tab,n){
  vtab=as.vector(tab)
  v=rmultinom(n=1,size=n,prob=vtab)
  m=matrix(dat=v,nrow=nrow(tab), ncol=ncol(tab))
  return(m)
}

power<-function(tab, test, n,  nSamples, cl){
  i=c(1:nSamples)
  sampleList=lapply(i, sample, tab,n)
  v=parSapply(cl,sampleList, test)
  #v=sapply(sampleList, test)
  return(sum(v==TRUE)/nSamples)
}

