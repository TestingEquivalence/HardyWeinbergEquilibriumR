# 1. example
# rheumatoid arthritis (RA) study
# Wordsworth P, Pile KD, Buckley JD, Lanchbury JSS, Ollier B, Lathrop M and Bell JI (1992)
# HLA heterozygosity contributes to susceptibility to rheumatoid arthritis. 
# Am J Hum Genet 51:585-591.
example1 = matrix(data=c(5, 0, 0, 0,
                         40, 12, 0, 0,
                         6, 32, 2, 0,
                         30, 55, 15, 33),
                  nrow=4, ncol=4, byrow=TRUE)

# 2. example
# Rousset, F., 2008 genepop’007: 
# A complete re-implementation of the genepop software for Windows and Linux. 
# Mol. Ecol. Res. 8: 103–106.
example2 = matrix(data=c(2, 0, 0, 0,
                         12, 24, 0, 0,
                         30, 34, 54, 0,
                         22, 21, 20, 10),
                  nrow=4, ncol=4, byrow=TRUE)

# 3. example
# Genotype frequency data at Rhesus locus. 
# Cavalli-Sforza, L. L. and Bodmer, W. F. (1971). 
# The Genetics of Human Populations. San Francisco: W. H. Freeman

example3= matrix(data=c(1236,   0,   0,   0,   0,   0,   0,  0,  0,
                        120,    3,   0,   0,   0,   0,   0,  0,  0,
                        18,     0,   0,   0,   0,   0,   0,  0,  0,
                        982,   55,   7, 249,   0,   0,   0,  0,  0,
                        32,     1,   0,  12,   0,   0,   0,  0,  0,
                        2582, 132,  20,1162,  29,1312,   0,  0,  0,
                        6,      0,   0,   4,   0,   4,   0,  0,  0, 
                        2,      0,   0,   0,   0,   0,   0,  0,  0,
                        115,    5,   2,   53,  1, 149,   0,  0,  4),
                        nrow=9, ncol=9, byrow=TRUE)

# 4. example
# Genotype frequency data from Louis and Dempster (1987)
# Louis, E. J. and Dempster, E. R. (1987). 
# An exact test for Hardy-Weinberg and multiple alleles.
# Biometrics 43, 805-811. 
example4 = matrix(data=c( 0,  0,  0,  0,
                          3,  1,  0,  0,
                          5, 18,  1,  0,
                          3,  7,  5,  2),
                  nrow=4, ncol=4, byrow=TRUE)
