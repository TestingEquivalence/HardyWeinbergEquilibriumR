product<-function(p){
  res=p %*% t(p)
  return(res)
}

triangle<-function(tab){
  #set upper triangle to zero
  us=upper.tri(tab,diag=FALSE)
  ttab=tab
  ttab[us]=0
  
  #add upper triangle
  ls=lower.tri(tab,diag=TRUE)
  rtab=tab
  rtab[ls]=0
  rtab=t(rtab)
  ttab=ttab+rtab
  
  return(ttab)
}

l22<-function(tab1,tab2){
  t=tab1-tab2
  t=t*t
  s=sum(t)
  return(s)
}

l22_first_derivative<-function(tab1,tab2){
  res=tab1-tab2
  res=2*res
  return(res)
}

startValue<-function(tab){
  r=rowSums(tab)
  c=colSums(tab)
  p=(r+c)/2
  return(p)
}

p2triangle<-function(p){
  #normalize
  p=p/sum(p)
  
  prod=product(p)
  t=triangle(prod)
  return(t)
}

fn<-function(p, table){
  t=p2triangle(p)
  dst=l22(table,t)
  return(dst)
}

min_l22<-function(table){
   start=startValue(table)
   l=length(start)
   lowerBound=rep(0,l)
   upperBound=rep(1,l)
   opRes=optim(par=start,fn=fn, gr=NULL, method="L-BFGS-B", lower=lowerBound,
               upper=upperBound, table=table)
   return(opRes)               
}

cond_l22<-function(tab){
  p=startValue(tab)
  prod=product(p)
  t=triangle(prod)
  
  dst=l22(tab,t)
  return(dst)
}

cond_l22_da<-function(p,i,j,k,m){
  res=0

  if (j>i) return(res)
    
  if ((i==k) & (j!=m)) res=-p[j]
  if ((i!=k) & (j==m)) res=-p[i]
  if ((i==k) & (j==m)) res=1-(p[i]+p[j])  
  
  return(res)
}

cond_l22_db<-function(tab, tprod, p, k, m){
  diff=tab-tprod
  s= 0
  for (i in  1:nrow(tab))
    for (j in 1:ncol(tab)) 
      s = s + 2 * diff[i, j] * cond_l22_da(p,i, j, k,m)
  return(s) 
}

derivative_cond_l22<-function(tab){
  p=startValue(tab)
  tprod=triangle(product(p))
  
  d=matrix(data=NA,nrow=nrow(tab),ncol=ncol(tab))
  for (k in  1:nrow(tab))
    for (m in 1:ncol(tab)) 
      d[k,m]=cond_l22_db(tab,tprod,p,k,m)
  return(d)
}
