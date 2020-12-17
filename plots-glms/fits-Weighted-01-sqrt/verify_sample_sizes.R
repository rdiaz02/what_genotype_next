## Verify the sample sizes (e.g., no filtering on sampledProp)
library(data.table)
setDTthreads(threads = 0)

load("../../wide_g.RData")

tmp_data <- wide_g

## This is common for all

max_g_7 <- 6
max_g_10 <- 9

tmp_data <- tmp_data[
        !is.na(min_js) 
        & (
            ((numGenes == 7) & (sourceGenotype_nMut <= max_g_7))
            |
            ((numGenes == 10) & (sourceGenotype_nMut <= max_g_10))
        )
       , 
    ]

nrow(tmp_data) ## 2389559

sum(tmp_data$numGenes == 10) ## 1715676
sum(tmp_data$numGenes == 7) ## 673883

sum((tmp_data$numGenes == 10) & (tmp_data$detect == "uniform") & (tmp_data$size_split == 4000)) ## 195825
sum((tmp_data$numGenes == 7) & (tmp_data$detect == "uniform") & (tmp_data$size_split == 4000)) ## 76085


## Which are the same as

system('grep "nrow data set" *.Rout')


## Yes, there is a stron argument in favor of not removing based on
## sampleProp: we could have lots of genotypes, each rarely observed, that
## together show a common a pattern pattern. So any one of the genotypes
## is uncommon, but the class of genotypes is common. That is why weighting
## makes sense, filtering does not.
