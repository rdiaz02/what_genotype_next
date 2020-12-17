library(spatstat)
library(data.table)
library(memisc)
## In when_good_conditions_0677.R we also do tables like this

## Note that the last genotype is no longer present
load("../data_for_weighted_glmertree_plots.RData")


dt <- data_for_weighted_glmertree_plots


dt[, weighted.mean(min_js, sampledProp),
   by = .(numGenes, sample_size, detect)]
##     numGenes sample_size  detect     V1
##  1:        7          50   large 0.3749
##  2:        7         200   large 0.2847
##  3:        7        4000   large 0.2175
##  4:        7          50   small 0.3256
##  5:        7         200   small 0.2883
##  6:        7        4000   small 0.2736
##  7:        7          50 uniform 0.2603
##  8:        7         200 uniform 0.2269
##  9:        7        4000 uniform 0.2154
## 10:       10          50   large 0.5432
## 11:       10         200   large 0.4293
## 12:       10        4000   large 0.2972
## 13:       10          50   small 0.4075
## 14:       10         200   small 0.3409
## 15:       10        4000   small 0.3229
## 16:       10          50 uniform 0.3237
## 17:       10         200 uniform 0.2713
## 18:       10        4000 uniform 0.2600


dt[, spatstat::weighted.median(min_js, sampledProp),
   by = .(numGenes, sample_size, detect)]
##     numGenes sample_size  detect     V1
##  1:        7          50   large 0.2668
##  2:        7         200   large 0.1878
##  3:        7        4000   large 0.1425
##  4:        7          50   small 0.2473
##  5:        7         200   small 0.2199
##  6:        7        4000   small 0.2108
##  7:        7          50 uniform 0.1912
##  8:        7         200 uniform 0.1655
##  9:        7        4000 uniform 0.1581
## 10:       10          50   large 0.5795
## 11:       10         200   large 0.3353
## 12:       10        4000   large 0.2215
## 13:       10          50   small 0.3351
## 14:       10         200   small 0.2874
## 15:       10        4000   small 0.2775
## 16:       10          50 uniform 0.2585
## 17:       10         200 uniform 0.2213
## 18:       10        4000 uniform 0.2138


dtx1 <- dt[, .(wmean = weighted.mean(min_js, sampledProp),
               wmedian = weighted.median(min_js, sampledProp)),
           by = .(typeLandscape, numGenes, sample_size, detect)]

dtx1$typeLandscape <- factor(dtx1$typeLandscape, levels = c("Represent.", "Local maxima", "RMF"))

setkey(dtx1, typeLandscape, numGenes, sample_size, detect)

dtx1

options(digits = 3)
toLatex(ftable(xtabs(wmean ~ numGenes + sample_size + detect +   typeLandscape , data = dtx1)))


toLatex(ftable(xtabs(wmedian ~ numGenes + sample_size + detect +   typeLandscape , data = dtx1)))


### Recheck the output for overall quantiles, as shown in
##  plotting-lmertree-fits-weighted-01-sqrt.R

spatstat::weighted.quantile(x = dt$min_js,
                            w = dt$sampledProp,
                            probs = c(0, .05, .075, 0.08, 0.085, 0.09, 0.095,
                                      .1, .15, .2, .25, .3,
                                      .5, .6, .75, .90, 1.0)
                            )

##  0.0%    5.0%    7.5%    8.0%    8.5%    9.0%    9.5%   10.0%   15.0%   20.0% 
## 0.00000 0.00000 0.00462 0.00708 0.00907 0.01163 0.01396 0.01668 0.04216 0.06773 
##   25.0%   30.0%   50.0%   60.0%   75.0%   90.0%  100.0% 
## 0.09426 0.11935 0.23652 0.30603 0.46663 0.78483 1.00000 
## > 


## Sure there are more efficient ways to do this

fp <- function(data, pp = c(0, .05, .1, .15, .2, .25, .4, .5, .75, 1.0)) {
    mm <- c("min_js", "MHN_js", "CBN_js", "MHN_td_js", "CBN_td_js",
            "MCCBN_js", "MCCBN_td_js",
            "CAPRI_AIC_js", "CAPRI_BIC_js",
            "CBN_uw_js" , "MCCBN_uw_js",
            "CAPRESE_js", "OT_js", "OT_uw_js")

    mout <- matrix(NA, nrow = length(pp), ncol = length(mm))
    colnames(mout) <- mm
    rownames(mout) <- as.character(pp)
    
    for(i in 1:length(mm)) {
        mout[, i] <- spatstat::weighted.quantile(x = data[[mm[i]]],
                            w = data$sampledProp,
                            probs = pp
                            )
    }
    return(mout)
}




qall <- fp(dt)


toLatex(qall)

toLatex(fp(dt[(sample_size == 4000) & (detect == "uniform") & (typeLandscape == "Represent."),  ]))

toLatex(fp(dt[(sample_size == 4000) & (detect == "uniform") & (typeLandscape == "Local maxima"),  ]))

toLatex(fp(dt[(sample_size == 4000) & (detect == "uniform") & (typeLandscape == "RMF"),  ]))
