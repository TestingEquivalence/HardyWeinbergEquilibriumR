source("tests.R")
source("data_sets.R")

# significance level
alpha=0.05

#table of tests results for the conditional distance
results_conditional_distance=matrix(data=NA, nrow=4,ncol=3)
rownames(results_conditional_distance)=c("example 1","example 2","example 3", "example 4")
colnames(results_conditional_distance)=c("distance","asymptotic test","resampling test")

#table of tests results for the minimum distance
results_minimum_distance=matrix(data=NA, nrow=4,ncol=3)
rownames(results_minimum_distance)=c("example 1","example 2","example 3", "example 4")
colnames(results_minimum_distance)=c("distance","asymptotic test","resampling test")



# 1. example
# rheumatoid arthritis study
# Wordsworth P, Pile KD, Buckley JD, Lanchbury JSS, Ollier B, Lathrop M and Bell JI (1992)
# HLA heterozygosity contributes to susceptibility to rheumatoid arthritis. 
# Am J Hum Genet 51:585-591.
print(example1)

results_conditional_distance[1,1]=cond_l2(example1)
results_conditional_distance[1,2]=asymptotic_test_conditional(example1,alpha)
results_conditional_distance[1,3]=resampling_test_conditional(example1,alpha)

results_minimum_distance[1,1]=min_l2(example1)
results_minimum_distance[1,2]=asymptotic_test_minimum(example1,alpha)
results_minimum_distance[1,3]=resampling_test_minimum(example1,alpha)


# 2. example
# Rousset, F., 2008 genepop’007: 
# A complete re-implementation of the genepop software for Windows and Linux. 
# Mol. Ecol. Res. 8: 103–106.
print(example2)

results_conditional_distance[2,1]=cond_l2(example2)
results_conditional_distance[2,2]=asymptotic_test_conditional(example2,alpha)
results_conditional_distance[2,3]=resampling_test_conditional(example2,alpha)

results_minimum_distance[2,1]=min_l2(example2)
results_minimum_distance[2,2]=asymptotic_test_minimum(example2,alpha)
results_minimum_distance[2,3]=resampling_test_minimum(example2,alpha)


# 3. example
# Genotype frequency data at Rhesus locus.
# Cavalli-Sforza, L. L. and Bodmer, W. F. (1971). 
# The Genetics of Hutman Populations. San Francisco: W. H. Freeman

print(example3)

results_conditional_distance[3,1]=cond_l2(example3)
results_conditional_distance[3,2]=asymptotic_test_conditional(example3,alpha)
results_conditional_distance[3,3]=resampling_test_conditional(example3,alpha)

results_minimum_distance[3,1]=min_l2(example3)
results_minimum_distance[3,2]=asymptotic_test_minimum(example3,alpha)
results_minimum_distance[3,3]=resampling_test_minimum(example3,alpha)
