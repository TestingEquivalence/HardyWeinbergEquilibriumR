source("distance.R")
source("asymptotic_test.R")
source("bootstrap_test.R")
source("data_sets.R")

# significance level
alpha=0.05
nBstSamples=10000

#table of tests results for the conditional distance
nDir_cond_dst=matrix(data=NA, nrow=7,ncol=4)
colnames(nDir_cond_dst)=c("example 1","example 2","example 3", "example 4")
rownames(nDir_cond_dst)=c(5,10,20,50,100,200,500)

d=nrow(example1)
nDir_cond_dst[1,1]=bootstrap_test_conditional(example1,alpha,nBstSamples,d*5)
nDir_cond_dst[2,1]=bootstrap_test_conditional(example1,alpha,nBstSamples,d*10)
nDir_cond_dst[3,1]=bootstrap_test_conditional(example1,alpha,nBstSamples,d*20)
nDir_cond_dst[4,1]=bootstrap_test_conditional(example1,alpha,nBstSamples,d*50)
nDir_cond_dst[5,1]=bootstrap_test_conditional(example1,alpha,nBstSamples,d*100)
nDir_cond_dst[6,1]=bootstrap_test_conditional(example1,alpha,nBstSamples,d*200)
nDir_cond_dst[7,1]=bootstrap_test_conditional(example1,alpha,nBstSamples,d*500)

d=nrow(example2)
nDir_cond_dst[1,2]=bootstrap_test_conditional(example2,alpha,nBstSamples,d*5)
nDir_cond_dst[2,2]=bootstrap_test_conditional(example2,alpha,nBstSamples,d*10)
nDir_cond_dst[3,2]=bootstrap_test_conditional(example2,alpha,nBstSamples,d*20)
nDir_cond_dst[4,2]=bootstrap_test_conditional(example2,alpha,nBstSamples,d*50)
nDir_cond_dst[5,2]=bootstrap_test_conditional(example2,alpha,nBstSamples,d*100)
nDir_cond_dst[6,2]=bootstrap_test_conditional(example2,alpha,nBstSamples,d*200)
nDir_cond_dst[7,2]=bootstrap_test_conditional(example2,alpha,nBstSamples,d*500)

d=nrow(example3)
nDir_cond_dst[1,3]=bootstrap_test_conditional(example3,alpha,nBstSamples,d*5)
nDir_cond_dst[2,3]=bootstrap_test_conditional(example3,alpha,nBstSamples,d*10)
nDir_cond_dst[3,3]=bootstrap_test_conditional(example3,alpha,nBstSamples,d*20)
nDir_cond_dst[4,3]=bootstrap_test_conditional(example3,alpha,nBstSamples,d*50)
nDir_cond_dst[5,3]=bootstrap_test_conditional(example3,alpha,nBstSamples,d*100)
nDir_cond_dst[6,3]=bootstrap_test_conditional(example3,alpha,nBstSamples,d*200)
nDir_cond_dst[7,3]=bootstrap_test_conditional(example3,alpha,nBstSamples,d*500)


d=nrow(example4)
nDir_cond_dst[1,4]=bootstrap_test_conditional(example4,alpha,nBstSamples,d*5)
nDir_cond_dst[2,4]=bootstrap_test_conditional(example4,alpha,nBstSamples,d*10)
nDir_cond_dst[3,4]=bootstrap_test_conditional(example4,alpha,nBstSamples,d*20)
nDir_cond_dst[4,4]=bootstrap_test_conditional(example4,alpha,nBstSamples,d*50)
nDir_cond_dst[5,4]=bootstrap_test_conditional(example4,alpha,nBstSamples,d*100)
nDir_cond_dst[6,4]=bootstrap_test_conditional(example4,alpha,nBstSamples,d*200)
nDir_cond_dst[7,4]=bootstrap_test_conditional(example4,alpha,nBstSamples,d*500)

write.table(nDir_cond_dst,file="nDir_cond_dst.txt")

#table of tests results for the minimum distance
nDir_min_dst=matrix(data=NA, nrow=7,ncol=4)
colnames(nDir_min_dst)=c("example 1","example 2","example 3", "example 4")
rownames(nDir_min_dst)=c(5,10,20,50,100,200,500)

d=nrow(example1)
nDir_min_dst[1,1]=bootstrap_test_minimum(example1,alpha,nBstSamples,d*5)
nDir_min_dst[2,1]=bootstrap_test_minimum(example1,alpha,nBstSamples,d*10)
nDir_min_dst[3,1]=bootstrap_test_minimum(example1,alpha,nBstSamples,d*20)
nDir_min_dst[4,1]=bootstrap_test_minimum(example1,alpha,nBstSamples,d*50)
nDir_min_dst[5,1]=bootstrap_test_minimum(example1,alpha,nBstSamples,d*100)
nDir_min_dst[6,1]=bootstrap_test_minimum(example1,alpha,nBstSamples,d*200)
nDir_min_dst[7,1]=bootstrap_test_minimum(example1,alpha,nBstSamples,d*500)


d=nrow(example2)
nDir_min_dst[1,2]=bootstrap_test_minimum(example2,alpha,nBstSamples,d*5)
nDir_min_dst[2,2]=bootstrap_test_minimum(example2,alpha,nBstSamples,d*10)
nDir_min_dst[3,2]=bootstrap_test_minimum(example2,alpha,nBstSamples,d*20)
nDir_min_dst[4,2]=bootstrap_test_minimum(example2,alpha,nBstSamples,d*50)
nDir_min_dst[5,2]=bootstrap_test_minimum(example2,alpha,nBstSamples,d*100)
nDir_min_dst[6,2]=bootstrap_test_minimum(example2,alpha,nBstSamples,d*200)
nDir_min_dst[7,2]=bootstrap_test_minimum(example2,alpha,nBstSamples,d*500)

d=nrow(example3)
nDir_min_dst[1,3]=bootstrap_test_minimum(example3,alpha,nBstSamples,d*5)
nDir_min_dst[2,3]=bootstrap_test_minimum(example3,alpha,nBstSamples,d*10)
nDir_min_dst[3,3]=bootstrap_test_minimum(example3,alpha,nBstSamples,d*20)
nDir_min_dst[4,3]=bootstrap_test_minimum(example3,alpha,nBstSamples,d*50)
nDir_min_dst[5,3]=bootstrap_test_minimum(example3,alpha,nBstSamples,d*100)
nDir_min_dst[6,3]=bootstrap_test_minimum(example3,alpha,nBstSamples,d*200)
nDir_min_dst[7,3]=bootstrap_test_minimum(example3,alpha,nBstSamples,d*500)

d=nrow(example4)
nDir_min_dst[1,4]=bootstrap_test_minimum(example4,alpha,nBstSamples,d*5)
nDir_min_dst[2,4]=bootstrap_test_minimum(example4,alpha,nBstSamples,d*10)
nDir_min_dst[3,4]=bootstrap_test_minimum(example4,alpha,nBstSamples,d*20)
nDir_min_dst[4,4]=bootstrap_test_minimum(example4,alpha,nBstSamples,d*50)
nDir_min_dst[5,4]=bootstrap_test_minimum(example4,alpha,nBstSamples,d*100)
nDir_min_dst[6,4]=bootstrap_test_minimum(example4,alpha,nBstSamples,d*200)
nDir_min_dst[7,4]=bootstrap_test_minimum(example4,alpha,nBstSamples,d*500)

write.table(nDir_min_dst,file="nDir_min_dst.txt")
