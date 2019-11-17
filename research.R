source("distance.R")
source("asymptotic_test.R")
source("bootstrap_test.R")
source("data_sets.R")
source("simulation.R")
source("size.R")
source("power.R")

# We simulate the test power at the Hardy Weinberg equilibrium points, 
# which rand distribution (distribution of alleles) equals the given data set.
# The power is simulated for different values of tolerance parameter epsilon to
# obtain the overview of the II. type error at the Hardy Weinberg equilibrium
# conditionally on the alleles distribution.
# Note that the II. type error equals 1-power.

# number of samples to estimate test power
nSamples=1000

# Selector contains logical values, which select the evaluation 
# for the different tests as follows: 
# selector[1] for the asymptotic test based on the conditional distance
# selector[2] for the bootstrap test based on the conditional distance
# selector[3] for the asymptotic test based on the minimum distance
# selector[4] for the bootstrap test based on the minimum distance
# Note that the power calculation needs a very long time (days) for the bootstrap tests.
# Power calculation uses parallel processing. Thus a strong parallel hardware is advantageous.
selector=c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)


# The number of bootstrap simulations is set to 2000 instead of 10000 in order to
# make power calcualation numerically feasible.
nSimulation=2000


#power at sample1
sizeSample1=matrix(data=NA, nrow=10, ncol=7)
colnames(sizeSample1)=c("eps","asy_cond","bst_cnd","asy_min","bst_min","res_cond","res_min")

sizeSample1[1,]=powerAtPoint(example1,0.12,nSamples,selector,nSimulation)
sizeSample1[2,]=powerAtPoint(example1,0.10,nSamples,selector,nSimulation)
sizeSample1[3,]=powerAtPoint(example1,0.09,nSamples,selector,nSimulation)
sizeSample1[4,]=powerAtPoint(example1,0.08,nSamples,selector,nSimulation)
sizeSample1[5,]=powerAtPoint(example1,0.07,nSamples,selector,nSimulation)

write.table(sizeSample1, "sizeSample1.txt")

# power at sample2
sizeSample2=matrix(data=NA, nrow=10, ncol=7)
colnames(sizeSample2)=c("eps","asy_cond","bst_cnd","asy_min","bst_min","res_cond","res_min")

sizeSample2[1,]=powerAtPoint(example2,0.12,nSamples,selector,nSimulation)
sizeSample2[2,]=powerAtPoint(example2,0.10,nSamples,selector,nSimulation)
sizeSample2[3,]=powerAtPoint(example2,0.09,nSamples,selector,nSimulation)
sizeSample2[4,]=powerAtPoint(example2,0.08,nSamples,selector,nSimulation)
sizeSample2[5,]=powerAtPoint(example2,0.07,nSamples,selector,nSimulation)

write.table(sizeSample2, "sizeSample2.txt")

# power at sample3
sizeSample3=matrix(data=NA, nrow=10, ncol=7)
colnames(sizeSample3)=c("eps","asy_cond","bst_cnd","asy_min","bst_min","res_cond","res_min")

sizeSample3[1,]=powerAtPoint(example3,0.018,nSamples,selector,nSimulation)
sizeSample3[2,]=powerAtPoint(example3,0.016,nSamples,selector,nSimulation)
sizeSample3[3,]=powerAtPoint(example3,0.014,nSamples,selector,nSimulation)
sizeSample3[4,]=powerAtPoint(example3,0.012,nSamples,selector,nSimulation)
sizeSample3[5,]=powerAtPoint(example3,0.010,nSamples,selector,nSimulation)


write.table(sizeSample3, "sizeSample3.txt")


# power at sample4
sizeSample4=matrix(data=NA, nrow=10, ncol=7)
colnames(sizeSample4)=c("eps","asy_cond","bst_cnd","asy_min","bst_min","res_cond","res_min")

sizeSample4[1,]=powerAtPoint(example4,0.06,nSamples,selector,nSimulation)
sizeSample4[2,]=powerAtPoint(example4,0.05,nSamples,selector,nSimulation)
sizeSample4[3,]=powerAtPoint(example4,0.04,nSamples,selector,nSimulation)
sizeSample4[4,]=powerAtPoint(example4,0.03,nSamples,selector,nSimulation)


write.table(sizeSample4, "sizeSample4.txt")

# power at sample5
sizeSample5=matrix(data=NA, nrow=10, ncol=7)
colnames(sizeSample5)=c("eps","asy_cond","bst_cnd","asy_min","bst_min","res_cond","res_min")

sizeSample5[1,]=powerAtPoint(example5,0.018,nSamples,selector,nSimulation)
sizeSample5[2,]=powerAtPoint(example5,0.016,nSamples,selector,nSimulation)
sizeSample5[3,]=powerAtPoint(example5,0.014,nSamples,selector,nSimulation)
sizeSample5[4,]=powerAtPoint(example5,0.012,nSamples,selector,nSimulation)
sizeSample5[5,]=powerAtPoint(example5,0.010,nSamples,selector,nSimulation)

write.table(sizeSample5, "sizeSample5.txt")

# Next we investigate if the type II error is sensitive to the allele distribution,
# because the true allele distribution is unknown and the estimator is subject to
# sampling error.
# For this porpose we estimate the allele distribution from the counting data.
# Then the random samples from the allele distribution are generated and the sample size
# equals these of the original data set. We generate 100 random samples which is sufficient
# for the considered data sets, because the test power vary very little from point to point.
# Then the test power is computed for each random sample under assumption of Hardy Weinberg equilibrium.

sensi_example1=powerSensitivity(tab=example1,eps=0.1,nSamples,selector,nSimulation)
write.table(sensi_example1,"sensi_example1.txt")

sensi_example2=powerSensitivity(tab=example2,eps=0.1,nSamples,selector,nSimulation)
write.table(sensi_example2,"sensi_example2.txt")

sensi_example3=powerSensitivity(tab=example3,eps=0.016,nSamples,selector,nSimulation)
write.table(sensi_example3,"sensi_example3.txt")

sensi_example4=powerSensitivity(tab=example4,eps=0.05,nSamples,selector,nSimulation)
write.table(sensi_example4,"sensi_example4.txt")

sensi_example5=powerSensitivity(tab=example5,eps=0.016,nSamples,selector,nSimulation)
write.table(sensi_example5,"sensi_example5.txt")


# The I. type  error is an important measure of the test quality and in the best case it should 
# not exceed the confidence level alpha.
# The I. type eroor is the maximum of the test power over the boundary of H0.
# In our case the boundary of H0 is a very complex set and the I. type error can not be estimated exactly.
# Instead we generate the random boundary points, which are close to the sample distribution.
# The test power at these boundary points provides a good insight into the test behavior at the bondary 
# of H0 and is a best practice estimate for the I.type error.
# The boundary points are generated as follows:
# - Random sample from original sata set is generated usind common bootstrap resampling.
# - If the distance between random sample and Hardy Weinberg equilibrium is smaller than epsilon,
#   we reject the random sample and repeat previous step.
# - We build a linear combination between the random sample and corresponding Hardy Weinberg equilibrium
#   point. 
#   The corresponding Hardy Weinberg equilibrium point depends on distance, which is used for testing.
# - The weight of the linear combination is adjusted such that 
#   the linear combination  is at the boundary of H0.

nSamples=1000
nSimulation=2000
#selector=c(TRUE,TRUE,TRUE,TRUE)
selector=c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE)
alpha=0.05
scaleFactor=1

powerExample1=boundaryPower(tab=example1,nSamples,selector,nSimulation,alpha, scaleFactor)
write.table(powerExample1,"powerExample1.txt")
 
powerExample2=boundaryPower(tab=example2,nSamples,selector,nSimulation,alpha, scaleFactor)
write.table(powerExample2,"powerExample2.txt")
 
powerExample3=boundaryPower(tab=example3,nSamples,selector,nSimulation,alpha, scaleFactor)
write.table(powerExample3,"powerExample3.txt")
 
powerExample4=boundaryPower(tab=example4,nSamples,selector,nSimulation,alpha, scaleFactor)
write.table(powerExample4,"powerExample4.txt")

powerExample5=boundaryPower(tab=example5,nSamples,selector,nSimulation,alpha,  scaleFactor)
write.table(powerExample5,"powerExample5.txt")

