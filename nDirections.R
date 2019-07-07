source("distance.R")
source("asymptotic_test.R")
source("bootstrap_test.R")

# significance level
alpha=0.05

#table of tests results for the conditional distance
nDir_cond_dst=matrix(data=NA, nrow=7,ncol=4)
colnames(nDir_cond_dst)=c("example 1","example 2","example 3", "example 4")
rownames(nDir_cond_dst)=c(5,10,20,50,100,200,500)

d=nrow(example1)
n