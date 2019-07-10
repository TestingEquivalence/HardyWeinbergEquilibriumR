sample<-function(i,tab,n){
  vtab=as.vector(tab)
  v=rmultinom(n=1,size=n,prob=vtab)
  m=matrix(dat=v,nrow=nrow(tab), ncol=ncol(tab))
  return(m)
}

power<-function(test, n, tab, nSamples){
  i=c(1:n)
  sampleList=lapply(list, function)
  return(sum(v==TRUE)/n)
}