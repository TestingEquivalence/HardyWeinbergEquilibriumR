source("distance.R")
source("asymptotic_test.R")
source("bootstrap_test2.R")
source("data_sets.R")

# significance level
alpha=0.05

#table of tests results for the conditional distance
results_conditional_distance=matrix(data=NA, nrow=4,ncol=3)
rownames(results_conditional_distance)=c("example 1","example 2","example 3", "example 4")
colnames(results_conditional_distance)=c("distance","asymptotic test","bootstrap test")

#table of tests results for the minimum distance
results_minimum_distance=matrix(data=NA, nrow=4,ncol=3)
rownames(results_minimum_distance)=c("example 1","example 2","example 3", "example 4")
colnames(results_minimum_distance)=c("distance","asymptotic test","bootstrap test")



# 1. example
# rheumatoid arthritis (RA) study
# Wordsworth P, Pile KD, Buckley JD, Lanchbury JSS, Ollier B, Lathrop M and Bell JI (1992)
# HLA heterozygosity contributes to susceptibility to rheumatoid arthritis. 
# Am J Hum Genet 51:585-591.
print(example1)

results_conditional_distance[1,1]=cond_l2(example1)
results_conditional_distance[1,2]=asymptotic_test_conditional(example1,alpha)
results_conditional_distance[1,3]=bootstrap_test_conditional(example1,alpha)

results_minimum_distance[1,1]=min_l2(example1)
results_minimum_distance[1,2]=asymptotic_test_minimum(example1,alpha)
results_minimum_distance[1,3]=bootstrap_test_minimum(example1,alpha)


# 2. example
# Rousset, F., 2008 genepop’007: 
# A complete re-implementation of the genepop software for Windows and Linux. 
# Mol. Ecol. Res. 8: 103–106.
print(example2)

results_conditional_distance[2,1]=cond_l2(example2)
results_conditional_distance[2,2]=asymptotic_test_conditional(example2,alpha)
results_conditional_distance[2,3]=bootstrap_test_conditional(example2,alpha)

results_minimum_distance[2,1]=min_l2(example2)
results_minimum_distance[2,2]=asymptotic_test_minimum(example2,alpha)
results_minimum_distance[2,3]=bootstrap_test_minimum(example2,alpha)


# 3. example
# Genotype frequency data at Rhesus locus.
# Cavalli-Sforza, L. L. and Bodmer, W. F. (1971). 
# The Genetics of Hutman Populations. San Francisco: W. H. Freeman

print(example3)

results_conditional_distance[3,1]=cond_l2(example3)
results_conditional_distance[3,2]=asymptotic_test_conditional(example3,alpha)
results_conditional_distance[3,3]=bootstrap_test_conditional(example3,alpha)

results_minimum_distance[3,1]=min_l2(example3)
results_minimum_distance[3,2]=asymptotic_test_minimum(example3,alpha)
results_minimum_distance[3,3]=bootstrap_test_minimum(example3,alpha)


# 4. example
# HFE genotype of the whole population referred for a personal or family
# history of iron overload (patients) 
# Györffy B., Kocsis I. and Vasarhely B.
# Biallelic genotype distributions in papers published in Gut between 1998 and 2003: 
# altered conclusions after recalculating the Hardy-Weinberg equilibrium
# Gut 2004;53:614–616

print(example4)

results_conditional_distance[4,1]=cond_l2(example4)
results_conditional_distance[4,2]=asymptotic_test_conditional(example4,alpha)
results_conditional_distance[4,3]=bootstrap_test_conditional(example4,alpha)

results_minimum_distance[4,1]=min_l2(example4)
results_minimum_distance[4,2]=asymptotic_test_minimum(example4,alpha)
results_minimum_distance[4,3]=bootstrap_test_minimum(example4,alpha)
 
write.table(results_conditional_distance,"results_conditional_distance.txt")
write.table(results_minimum_distance,"results_minimum_distance.txt")

