source("tests.R")
source("data_sets.R")
source("simulation.R")
source("size.R")
source("power.R")

# We simulate the test power at the Hardy Weinberg equilibrium points, 
# which allele distribution of alleles equals these of the given data set.
# The power is simulated for different values of tolerance parameter epsilon to
# obtain the overview of the II. type error and effective sample size
# at the Hardy Weinberg equilibrium conditionally on the alleles distribution.
# Note, that the II. type error equals 1-test power.

# number of samples to estimate test power
nSamples=1000
# tests to consider
testsToDo=tests

#power at example 1
sizeSample1=matrix(data=NA, nrow=5, ncol=5)
colnames(sizeSample1)=c("eps","asy_cond","res_cond","asy_min","res_min")

sizeSample1[1,]=powerAtPoint(example1,0.12,nSamples,testsToDo)
sizeSample1[2,]=powerAtPoint(example1,0.10,nSamples,testsToDo)
sizeSample1[3,]=powerAtPoint(example1,0.09,nSamples,testsToDo)
sizeSample1[4,]=powerAtPoint(example1,0.08,nSamples,testsToDo)
sizeSample1[5,]=powerAtPoint(example1,0.07,nSamples,testsToDo)

write.table(sizeSample1, "sizeSample1.txt")

# power at example 2
sizeSample2=matrix(data=NA, nrow=5, ncol=5)
colnames(sizeSample2)=c("eps","asy_cond","res_cond","asy_min","res_min")

sizeSample2[1,]=powerAtPoint(example2,0.12,nSamples,testsToDo)
sizeSample2[2,]=powerAtPoint(example2,0.10,nSamples,testsToDo)
sizeSample2[3,]=powerAtPoint(example2,0.09,nSamples,testsToDo)
sizeSample2[4,]=powerAtPoint(example2,0.08,nSamples,testsToDo)
sizeSample2[5,]=powerAtPoint(example2,0.07,nSamples,testsToDo)

write.table(sizeSample2, "sizeSample2.txt")

# power at example 3
sizeSample3=matrix(data=NA, nrow=5, ncol=5)
colnames(sizeSample3)=c("eps","asy_cond","res_cond","asy_min","res_min")

sizeSample3[1,]=powerAtPoint(example3,0.018,nSamples,testsToDo)
sizeSample3[2,]=powerAtPoint(example3,0.016,nSamples,testsToDo)
sizeSample3[3,]=powerAtPoint(example3,0.014,nSamples,testsToDo)
sizeSample3[4,]=powerAtPoint(example3,0.012,nSamples,testsToDo)
sizeSample3[5,]=powerAtPoint(example3,0.010,nSamples,testsToDo)


write.table(sizeSample3, "sizeSample3.txt")


# power at example 4
sizeSample4=matrix(data=NA, nrow=5, ncol=5)
colnames(sizeSample4)=c("eps","asy_cond","res_cond","asy_min","res_min")

sizeSample4[1,]=powerAtPoint(example4,0.018,nSamples,testsToDo)
sizeSample4[2,]=powerAtPoint(example4,0.016,nSamples,testsToDo)
sizeSample4[3,]=powerAtPoint(example4,0.014,nSamples,testsToDo)
sizeSample4[4,]=powerAtPoint(example4,0.012,nSamples,testsToDo)
sizeSample4[5,]=powerAtPoint(example4,0.010,nSamples,testsToDo)


write.table(sizeSample4, "sizeSample4.txt")

# Next we investigate if the type II error is sensitive to the allele distribution,
# because the true allele distribution is unknown and the estimator is subject to
# sampling error.
# For this porpose, we estimate the allele distribution from the counting data.
# Then the random samples from the allele distribution are generated and the sample size
# equals these of the original data set. We generate 100 random samples which is sufficient
# for the considered data sets, because the test power vary very little from point to point.
# Then the test power is computed for each random sample 
# under assumption of Hardy Weinberg equilibrium.

sensi_example1=powerSensitivity(tab=example1,eps=0.1,nSamples,testsToDo)
write.table(sensi_example1,"sensi_example1.txt")

sensi_example2=powerSensitivity(tab=example2,eps=0.1,nSamples,testsToDo)
write.table(sensi_example2,"sensi_example2.txt")

sensi_example3=powerSensitivity(tab=example3,eps=0.016,nSamples,testsToDo)
write.table(sensi_example3,"sensi_example3.txt")

# The I. type  error is an important measure of the test quality and 
# it should not exceed the confidence level alpha in the best case.
# The I. type error is the maximum of the test power over the boundary of H0.
# In our case the boundary of H0 is a very complex set and 
# the I. type error can not be estimated exactly.
# Instead, we generate random boundary points, which are close to the sample distribution.
# The test power at these boundary points provides a good insight 
# into the test behavior at the bondary of H0 and
# provides a best practice estimate for the I.type error.
# The boundary points are generated as follows:
# - Random sample from original sata set is generated usind common bootstrap resampling.
# - If the distance between random sample and Hardy Weinberg equilibrium is smaller
#   than epsilon, we reject the random sample and repeat previous step.
# - We build a linear combination between the random sample and 
#   corresponding Hardy Weinberg equilibrium equilibrium point. 
#   The corresponding Hardy Weinberg equilibrium point depends on distance,
#   which is used for testing.
# - The weight of the linear combination is adjusted such that 
#   the linear combination  is at the boundary of H0.

powerExample1=boundaryPower(tab=example1,nSamples,testsToDo)
write.table(powerExample1,"powerExample1.txt")
 
powerExample2=boundaryPower(tab=example2,nSamples,testsToDo)
write.table(powerExample2,"powerExample2.txt")
 
powerExample3=boundaryPower(tab=example3,nSamples,testsToDo)
write.table(powerExample3,"powerExample3.txt")
 

