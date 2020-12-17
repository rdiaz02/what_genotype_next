library(data.table)
library(parallel)
library(spatstat)
library(memisc)

## library(RhpcBLASctl)
## cores_blas_omp <- detectCores()
## RhpcBLASctl::blas_set_num_threads(cores_blas_omp)
## RhpcBLASctl::omp_set_num_threads(cores_blas_omp)

## From explore-reduce-method-comp.R
load("../../method_comp_fit.RData")

localdata <- method_comp_fit

max_g_7 <- 6
max_g_10 <- 9

localdata <- localdata[
(
            ((numGenes == 7) & (sourceGenotype_nMut <= max_g_7))
            |
            ((numGenes == 10) & (sourceGenotype_nMut <= max_g_10))
        ),
    ]

cat("\n nrow(localdata) post trimming last genotype")
nrow(localdata)


## Remember this is how we use weights in the fits
## Here we just want to use this overall. 
mfw <- nrow(localdata)/sum(localdata$sampledProp)
cat("\n Multiplication factor for weights is ", mfw, "\n")
localdata[, `:=`( sampledProp2 = mfw * sampledProp)]

## They are not 1, as we have removed genotypes: the last one
tmpweights <- localdata[, sum(sampledProp),
                        by = .(id, replicate, size_split, detect)]

summary(tmpweights$V1)
quantile(tmpweights$V1, probs = seq(0, 1, 0.05))
## > quantile(tmpweights$V1, probs = seq(0, 1, 0.05))
##    0%    5%   10%   15%   20%   25%   30%   35%   40%   45%   50%   55%   60% 
## 0.610 0.851 0.901 0.940 0.974 0.991 0.997 0.999 1.000 1.000 1.000 1.000 1.000 
##   65%   70%   75%   80%   85%   90%   95%  100% 
## 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000

tmpweights2 <- localdata[, sum(sampledProp2),
                        by = .(id, replicate, size_split, detect)]

summary(tmpweights2$V1)
quantile(tmpweights2$V1, probs = seq(0, 1, 0.05))


## This we did also in plots-glms/fits-Weighted-01-sqrt/min_js_by_detect_sample_size.R
## The file was created in plots-glms/data_for_weighted_glmertree_plots.R
load("../data_for_weighted_glmertree_plots.RData")

dt <- data_for_weighted_glmertree_plots
dim(dt)
dim(localdata) ## same number observations
summary(dt$MCCBN_js) ## recall 15 that are missing from MCCBN
summary(localdata$min_js)
summary(dt$min_js)

## Same as in plots-glms/fits-Weighted-01-sqrt/plotting-lmertree-fits-weighted-01-sqrt.R
spatstat::weighted.quantile(dt$min_js, dt$sampledProp,
                            probs = c(0.1, 0.2, 0.25, 0.5))


d1 <- dt[, weighted.mean(min_js, sampledProp),
    by = .(numGenes, sample_size, detect)]
setkey(d1, numGenes, sample_size, detect)
d1

d2 <- localdata[, weighted.mean(min_js, sampledProp),
          by = .(numGenes, size_split, detect)]

d3 <- localdata[, weighted.mean(min_js, sampledProp2),
          by = .(numGenes, size_split, detect)]

## So, yes, we are dealing with the same data and the weigthing resales
## changes nothing in this case
summary(d1$V1 - d2$V1)
summary(d2$V1 - d3$V1)

########

## define a few variables after a left join

## dt[, `:=` (good_TD = (pmin(MHN_td_js, CBN_td_js) <= 0.0677),
##            good_CE = (pmin(MHN_js, CBN_js, MCCBN_js) <= 0.0677)
##                )]
##

## Use pmax!! All are good
dt[, `:=` (good_TDm = (pmax(MHN_td_js, CBN_td_js) <= 0.0677),
           good_CEm = (pmax(MHN_js, CBN_js) <= 0.0677)
           )]

dt_tmp <- dt[, .(id, sourceGenotype, replicate, sample_size, detect,
                 min_js, good_TDm, good_CEm)]

dt_tmp <- dt_tmp[, `:=` (size_split = sample_size, old_min_js = min_js)]

dt_tmp <- dt_tmp[, .(id, sourceGenotype,
                     replicate, size_split,
                     detect,
                     old_min_js,
                     good_TDm, good_CEm)]
## old_min_js is just for checking

local2 <- merge(localdata, dt_tmp, all.x = TRUE)

## checks
dim(local2)
dim(localdata)
summary(local2$old_min_js - local2$min_js)


local2 <- local2[, `:=` (similar_TD = (TD <= 0.0677),
                         similar_CE = (CE <= 0.0677)
                       , different_CE_TD = (TD_comp_CE >= 0.7)
                         )]


options(digits = 3)

## Ignoring different_CE_TD


toLatex(ftable(prop.table(xtabs(wsum ~ good_TDm + similar_TD, 
                                local2[,
                                       .(wsum = sum(sampledProp2)),
                                       by = .(good_TDm, similar_TD)]))))


toLatex(ftable(prop.table(xtabs(wsum ~ good_CEm + similar_CE, 
                                local2[,
                                       .(wsum = sum(sampledProp2)),
                                       by = .(good_CEm, similar_CE)]))))




toLatex(ftable(prop.table(xtabs(wsum ~ different_CE_TD + similar_TD + good_TDm , 
                                local2[,
                                       .(wsum = sum(sampledProp2)),
                                       by = .(good_TDm, similar_TD, different_CE_TD)]))))


toLatex(ftable(prop.table(xtabs(wsum ~ different_CE_TD + similar_CE + good_CEm , 
                                local2[,
                                       .(wsum = sum(sampledProp2)),
                                       by = .(good_CEm, similar_CE, different_CE_TD)]))))



toLatex(ftable(prop.table(xtabs(wsum ~ different_CE_TD + similar_TD + similar_CE + good_CEm , 
                                local2[,
                                       .(wsum = sum(sampledProp2)),
                                       by = .(good_CEm, similar_TD, similar_CE, different_CE_TD)]))))


toLatex(ftable(prop.table(xtabs(wsum ~ different_CE_TD + similar_TD + similar_CE + good_TDm , 
                                local2[,
                                       .(wsum = sum(sampledProp2)),
                                       by = .(good_TDm, similar_TD, similar_CE, different_CE_TD)]))))




## xtabs(wsum ~ different_CE_TD + similar_CE + good_CEm , 
##                                 local2[,
##                                        .(wsum = sum(sampledProp2)),
##                                        by = .(good_CEm, similar_CE, different_CE_TD)]))))

