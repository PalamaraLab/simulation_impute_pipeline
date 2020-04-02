library(stringr)

pairs <- dir("/homes/shi/workspace/dev/mypsmc/analysis/simdecode/exp03_chip0.6_sim20k_ref1k_ref1000/test_outpairs/pairs")
id1 <- str_split(pairs, "-", simplify=TRUE)[, 1]
id2 <- str_replace(str_split(pairs, "-", simplify=TRUE)[, 2], ".geno", "")

hap1 <- str_split(id1, "_", simplify=TRUE)[, 3]
hap2 <- str_split(id2, "_", simplify=TRUE)[, 3]
id1 <- str_split(id1, "_", simplify=TRUE)[, 2]
id2 <- str_split(id2, "_", simplify=TRUE)[, 2]


