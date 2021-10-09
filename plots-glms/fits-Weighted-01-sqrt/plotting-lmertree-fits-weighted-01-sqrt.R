
source("../common_lmertree_fits.R")
source("../lmertree-grapcon-functions.R")
source("../common-lmertree-plotting-checking-pruning-functions.R")

library(ggplot2)
library(cowplot)
old <- theme_set(theme_cowplot())

#####################################

## Load data for plots
## created in /plots-glms/data_for_weighted_glmertree_plots.R
load("../data_for_weighted_glmertree_plots.RData.gz")
gc()
dim(data_for_weighted_glmertree_plots) ## 2389559     108

data_for_weighted_glmertree_plots[, `:=`( typeLandscape = typeLandscape_f)]
levels(data_for_weighted_glmertree_plots$typeLandscape_f) <- c("Rep", "LM", "RMF")

with(data_for_weighted_glmertree_plots, table(typeLandscape_f, typeLandscape))


## Recall we transformed min_js
data_for_weighted_glmertree_plots[, `:=`( min_js = sqrt(min_js))]


data_for_specific_weighted_10 <-
    data_for_weighted_glmertree_plots[(numGenes == 10) &
                             (sample_size == 4000) &
                             (detect == "uniform")
                             , ]

data_for_specific_weighted_7 <-
    data_for_weighted_glmertree_plots[(numGenes == 7) &
                             (sample_size == 4000) &
                             (detect == "uniform")
                             , ]

data_for_all_weighted_10 <-
    data_for_weighted_glmertree_plots[(numGenes == 10) , ]

data_for_all_weighted_7 <-
    data_for_weighted_glmertree_plots[(numGenes == 7) , ]

## The "nrow data set" for the fits are (for the checks)
##  7, specific:   76085
## 10, specific:  195825
##  7, all     :  673883
## 10, all     : 1715676



## Load fitted glmertrees

load("./keep-fit-RData-W-01-sqrt/lmertree_sqrt_W_Obs_minsize_01_fl_specific_7genes_gamma_epist_withFL_Inf_NM_fits.RData.gz")
W_specific_7 <- lmertree_sqrt_W_Obs_minsize_01_fl_specific_7genes_gamma_epist_withFL_Inf_NM_allFL_Inf


load("./keep-fit-RData-W-01-sqrt/lmertree_sqrt_W_Obs_minsize_01_fl_specific_10genes_gamma_epist_withFL_Inf_NM_fits.RData.gz")
W_specific_10 <- lmertree_sqrt_W_Obs_minsize_01_fl_specific_10genes_gamma_epist_withFL_Inf_NM_allFL_Inf

load("./keep-fit-RData-W-01-sqrt/lmertree_sqrt_W_Obs_minsize_01_fLandscape-rsign-peaks-recoded_7genes_Inf_NM_fits.RData.gz")
W_all_7 <- get("lmertree_sqrt_W_Obs_minsize_01_fLandscape-rsign-peaks-recoded_7genes_Inf_NM_allFL_Inf")

load("./keep-fit-RData-W-01-sqrt/lmertree_sqrt_W_Obs_minsize_01_fLandscape-rsign-peaks-recoded_10genes_Inf_NM_fits.RData.gz")
W_all_10 <- get("lmertree_sqrt_W_Obs_minsize_01_fLandscape-rsign-peaks-recoded_10genes_Inf_NM_allFL_Inf")





## Check (transformation done correctly and we are fitting what we think)


stopifnot(isTRUE(all.equal(
    weighted_summary_data(data_for_specific_weighted_7),
    weighted_summary(data_party(W_specific_7$tree, 1)$min_js,
                     data_party(W_specific_7$tree, 1)[["(weights)"]])
)))


stopifnot(isTRUE(all.equal(
    weighted_summary_data(data_for_specific_weighted_10),
    weighted_summary(data_party(W_specific_10$tree, 1)$min_js,
                     data_party(W_specific_10$tree, 1)[["(weights)"]])
)))


stopifnot(isTRUE(all.equal(
    weighted_summary_data(data_for_all_weighted_7),
    weighted_summary(data_party(W_all_7$tree, 1)$min_js,
                     data_party(W_all_7$tree, 1)[["(weights)"]])
)))


stopifnot(isTRUE(all.equal(
    weighted_summary_data(data_for_all_weighted_10),
    weighted_summary(data_party(W_all_10$tree, 1)$min_js,
                     data_party(W_all_10$tree, 1)[["(weights)"]])
)))



weighted_summary_data(data_for_specific_weighted_7)
weighted_summary_data2(data_for_specific_weighted_7) ## a lot faster

how_many_less_than_js(data_for_specific_weighted_7) # 0.48

Hmisc::wtd.quantile(data_for_specific_weighted_7$min_js,
             weights = data_for_specific_weighted_7$sampledProp,
             normwt = TRUE,
             probs = 0.4823749)

weighted_summary_data2(data_for_specific_weighted_10)
how_many_less_than_js(data_for_specific_weighted_10) # 0.38

weighted_summary_data2(data_for_all_weighted_7)
how_many_less_than_js(data_for_all_weighted_7) # 0.42

weighted_summary_data2(data_for_all_weighted_10)
how_many_less_than_js(data_for_all_weighted_10) # 0.29


how_many_less_than_js(data_for_specific_weighted_7,
                      threshold = sqrt(0.1)) ## .37
how_many_less_than_js(data_for_specific_weighted_10, threshold = sqrt(0.1)) ## .27
how_many_less_than_js(data_for_all_weighted_7, threshold = sqrt(0.1)) ## .31
how_many_less_than_js(data_for_all_weighted_10, threshold = sqrt(0.1)) ## .21


## > weighted_summary_data2(data_for_specific_weighted_7) ## a lot faster
##        0%       10%       15%       20%       25%       40%       50%       75%      100% 
## 0.0000000 0.0000000 0.1212250 0.1841334 0.2320607 0.3344770 0.3975865 0.5577398 1.0000000 
## > weighted_summary_data2(data_for_specific_weighted_10)
##        0%       10%       15%       20%       25%       40%       50%       75%      100% 
## 0.0000000 0.1140453 0.1932585 0.2475516 0.2949771 0.4011804 0.4623327 0.6147035 1.0000000 
## > weighted_summary_data2(data_for_all_weighted_7)
##         0%        10%        15%        20%        25%        40%        50%        75%       100% 
## 0.00000000 0.07873796 0.16899052 0.22232117 0.26891747 0.37609455 0.44457776 0.63379873 1.00000000 
## > weighted_summary_data2(data_for_all_weighted_10)
##        0%       10%       15%       20%       25%       40%       50%       75%      100% 
## 0.0000000 0.1789731 0.2499359 0.3051627 0.3497459 0.4636725 0.5290890 0.7276738 1.0000000 
## > 

weighted_summary_data2(data_for_specific_weighted_7, c(0, .1, .15, .2, .25, .3), TRUE)
weighted_summary_data2(data_for_specific_weighted_10, c(0, .1, .15, .2, .25, .3), TRUE)
weighted_summary_data2(data_for_all_weighted_7, c(0, .1, .15, .2, .25, .3), TRUE)
weighted_summary_data2(data_for_all_weighted_10, c(0, .1, .15, .2, .25, .3), TRUE)


## > weighted_summary_data2(data_for_specific_weighted_7, c(0, .1, .15, .2, .25, .3), TRUE)
##         0%        10%        15%        20%        25%        30% 
## 0.00000000 0.00000000 0.01469551 0.03390511 0.05385216 0.07367471 
## > weighted_summary_data2(data_for_specific_weighted_10, c(0, .1, .15, .2, .25, .3), TRUE)
##         0%        10%        15%        20%        25%        30% 
## 0.00000000 0.01300633 0.03734885 0.06128178 0.08701148 0.11048488 
## > weighted_summary_data2(data_for_all_weighted_7, c(0, .1, .15, .2, .25, .3), TRUE)
##          0%         10%         15%         20%         25%         30% 
## 0.000000000 0.006199666 0.028557796 0.049426703 0.072316608 0.094292732 
## > weighted_summary_data2(data_for_all_weighted_10, c(0, .1, .15, .2, .25, .3), TRUE)
##         0%        10%        15%        20%        25%        30% 
## 0.00000000 0.03203138 0.06246794 0.09312424 0.12232222 0.15474827 
## > 

threshold_fits <- c(
    weighted_summary_data2(data_for_specific_weighted_7, c(.25), TRUE),
    weighted_summary_data2(data_for_specific_weighted_10, c(.25), TRUE),
    weighted_summary_data2(data_for_all_weighted_7, c(.25), TRUE),
    weighted_summary_data2(data_for_all_weighted_10, c(.25), TRUE)
)


weighted_summary_data2(data_for_weighted_glmertree_plots, sqrt_inv = TRUE)
##        0%        10%        15%        20%        25%        40%        50%        75%       100% 
## 0.00000000 0.01668478 0.04215511 0.06773060 0.09426294 0.17554359 0.23651955 0.46662663 1.00000000 
## >
 

## > weighted_summary_data2(data_for_weighted_glmertree_plots, sqrt_inv = TRUE)
##         0%        10%        15%        20%        25%        30%        35%        40%        50%        75% 
## 0.00000000 0.01668478 0.04215511 0.06773060 0.09426294 0.11934654 0.14776386 0.17554359 0.23651955 0.46662663 
##       100% 
## 1.00000000

how_many_less_than_js(data_for_weighted_glmertree_plots, sqrt(0.13)) ## 0.32
how_many_less_than_js(data_for_weighted_glmertree_plots, sqrt(0.18)) ## 0.4

## Much better in general
##  set last argument to TRUE if you want to see the plot (not by default,
##  as it can break some non-interactive work)
scatterplot_median_intercept(W_specific_7$tree, TRUE)
## seems much better: cors: 0.94, 0.94, slope .86
scatterplot_median_intercept(W_specific_10$tree, TRUE) ## .92, .91, .97
scatterplot_median_intercept(W_all_7$tree, TRUE, FALSE) ## .85, .90, .63: ## yes, not great
scatterplot_median_intercept(W_all_10$tree, TRUE, FALSE)
## 0.88, 0.87, slope 0.8, bu 1.05 if through origin




do_prune_again <- FALSE ## pruned trees are precomputed for
## repeated runs where we change grpahical stuff. No need to redo this again
## If you've never pruned, well prune here

if(do_prune_again) {
W_specific_7_prunedb_sq <- fast_prune_sqrt(W_specific_7$tree, maxFit = 0.06773060,
                                minN = 100,  mindiff = 0.06773060/4,
                                statistic = "intercept", weighted = TRUE)


W_specific_10_prunedb_sq <- fast_prune_sqrt(W_specific_10$tree, maxFit = 0.06773060,
                                minN = 100,  mindiff = 0.06773060/4,
                                statistic = "intercept", weighted = TRUE)


W_all_7_prunedb_sq <- fast_prune_sqrt(W_all_7$tree, maxFit = 0.06773060,
                              minN = 100,  mindiff = 0.06773060/4,
                              statistic = "intercept", weighted = TRUE)

W_all_10_prunedb_sq <- fast_prune_sqrt(W_all_10$tree, maxFit = 0.06773060,
                              minN = 100,  mindiff = 0.06773060/4,
                              statistic = "intercept", weighted = TRUE)


## For faster loading 
save(file = "pruned_trees_sqrt.RData", list = ls(pattern = glob2rx("W_*_prunedb_sq")),
     compress = FALSE)
## compress it for permanent storage
}

load("pruned_trees_sqrt.RData.gz")


scatterplot_median_intercept(W_specific_7_prunedb_sq, TRUE, plot = FALSE)
## .94, .94
scatterplot_median_intercept(W_specific_10_prunedb_sq, TRUE, plot = FALSE)
## .95, .93
scatterplot_median_intercept(W_all_7_prunedb_sq, TRUE, plot = FALSE)
## .84, .81 but slope is low: 0.56, but 1.16 if forced to 0, so OK
scatterplot_median_intercept(W_all_10_prunedb_sq, TRUE, plot = FALSE)
## .83, .87, slope .63 and 1.03 if through origin


print(warnings())

all_vars <- c("MHN"
            , "MHN_td"
            , "CBN"
            , "CAPRESE"
            , "OT"
            , "CBN_td"
            , "MCCBN"
            , "CAPRI_BIC"
            , "CAPRI_AIC"
            , "OT_uw"
            , "MCCBN_td"
            , "MCCBN_uw"
            , "CBN_uw"
              )


all_vars_except_uw <- c("MHN"
                      , "MHN_td"
                      , "CBN"
                      , "CAPRESE"
                      , "OT"
                      , "CBN_td"
                      , "MCCBN"
                      , "CAPRI_BIC"
                      , "CAPRI_AIC"
                      ## , "OT_uw"
                      , "MCCBN_td"
                      ## , "MCCBN_uw"
                      ## , "CBN_uw"
              )


## Yes, weighted boxplots and means.
##    weighted argument passed to created_terminal_plot_data_2
##    which returns a vector of weights
##    inside my_node_bp5 we use colWeightedMeans (from matrixStats)
##    That for the weighted means of performance
##   And the box plot uses the weights: x <- boxplot(rep.int(yn, wn), plot = FALSE)
## And yes, back transformed, so not the sqrt transformation
##     we pass the sqrt_transf = TRUE
##     used by my_node_bp5
##     if(sqrt_transf) yn <- yn^2





## Following reviewer's request, we produce figures excluding the unweighted,
## to unclutter. Figures with all methods are also provided in the suppl mat.

data_and_plot(W_specific_7_prunedb_sq,
              data_for_specific_weighted_7, 76085L,
              main = "weighted-7-specific-no-split-01-sqrt-all_methods_except_uw-js",
              weighted = TRUE,
              size_y_text = 12,
              height = 22,
              ## width = 30,
              best_vars = all_vars_except_uw,
              sqrt_transf  = TRUE, show_bottom = "js")

data_and_plot(W_specific_10_prunedb_sq,
              data_for_specific_weighted_10, 195825L,
              main = "weighted-10-specific-no-split-01-sqrt-all_methods_except_uw-js",
              weighted = TRUE,
              size_y_text = 12,
              height = 22,
              best_vars = all_vars_except_uw,
              sqrt_transf  = TRUE, show_bottom = "js")


data_and_plot(W_all_7_prunedb_sq,
              data_for_all_weighted_7, 673883L,
              main = "weighted-7-ALL-no-split-01-sqrt-all_methods_except_uw-js",
              weighted = TRUE,
              size_y_text = 12,
              height = 22,
              best_vars = all_vars_except_uw,
              sqrt_transf  = TRUE, show_bottom = "js")

data_and_plot(W_all_10_prunedb_sq,
              data_for_all_weighted_10, 1715676L,
              main = "weighted-10-ALL-no-split-01-sqrt-all_methods_except_uw-js",
              weighted = TRUE,
              size_y_text = 12,
              height = 22,
              best_vars = all_vars_except_uw,
              sqrt_transf  = TRUE, show_bottom = "js")



## These are for supplementary material only and include all methods.  As seen
## above, we also create figures, one for main text, the rest for suppl mat,
## which contain only a subset of methods: all except the unweighted.

## Since these have not changed, we do not rerun them for revised
## version of ms.
if(FALSE) {
    data_and_plot(W_specific_7_prunedb_sq,
                  data_for_specific_weighted_7, 76085L,
                  main = "weighted-7-specific-no-split-01-sqrt-all_methods-js",
                  weighted = TRUE,
                  size_y_text = 8,
                  height = 22,
                  ## width = 30,
                  best_vars = all_vars,
                  sqrt_transf  = TRUE, show_bottom = "js")
    data_and_plot(W_specific_10_prunedb_sq,
                  data_for_specific_weighted_10, 195825L,
                  main = "weighted-10-specific-no-split-01-sqrt-all_methods-js",
                  weighted = TRUE,
                  size_y_text = 8,
                  height = 22,
                  best_vars = all_vars,
                  sqrt_transf  = TRUE, show_bottom = "js")


    data_and_plot(W_all_7_prunedb_sq,
                  data_for_all_weighted_7, 673883L,
                  main = "weighted-7-ALL-no-split-01-sqrt-all_methods-js",
                  weighted = TRUE,
                  size_y_text = 8,
                  height = 22,
                  best_vars = all_vars,
                  sqrt_transf  = TRUE, show_bottom = "js")

    data_and_plot(W_all_10_prunedb_sq,
                  data_for_all_weighted_10, 1715676L,
                  main = "weighted-10-ALL-no-split-01-sqrt-all_methods-js",
                  weighted = TRUE,
                  size_y_text = 8,
                  height = 22,
                  best_vars = all_vars,
                  sqrt_transf  = TRUE, show_bottom = "js")
}

list_fits <- ls(pattern = glob2rx("W_*"))

minint <- lapply(list_fits, function(x) {
    if (!inherits(get(x), "modelparty")) {
         obj <- get(x)$tree
    } else {
        obj <- get(x)
    }
    min(unlist(nodeapply(obj, ids = nodeids(obj, terminal = TRUE),
                     FUN = function(n)
                         n$info$coefficients[["(Intercept)"]]))
        )})

lapply(minint, function(x) x^2)

q()


## [[1]]
## [1] 0.009550728

## [[2]]
## [1] 0.009550728

## [[3]]
## [1] 0.0163032

## [[4]]
## [1] 0.0163032

## [[5]]
## [1] 0.01227104

## [[6]]
## [1] 0.01227104

## [[7]]
## [1] 7.286324e-05

## [[8]]
## [1] 0.0005834646

## list_fits
## [1] "W_all_10"                 "W_all_10_prunedb_sq"      "W_all_7"                  "W_all_7_prunedb_sq"      
## [5] "W_specific_10"            "W_specific_10_prunedb_sq" "W_specific_7"             "W_specific_7_prunedb_sq"



### Other presentations; not used. 


## ## Paper main figure
## data_and_plot(W_specific_7_prunedb,
##               data_for_specific_weighted_7, 76085L,
##               main = "weighted-7-specific-no-split-01-sqrt-best",
##               weighted = TRUE,
##               sqrt_transf = TRUE,
##               size_y_text = 10,
##               height = 22)

## ## Suppl mat: all methods
## data_and_plot(W_specific_7_prunedb,
##               data_for_specific_weighted_7, 76085L,
##               main = "weighted-7-specific-no-split-01-sqrt-all_methods-best",
##               weighted = TRUE,
##               size_y_text = 8,
##               height = 22,
##               best_vars = all_vars,
##               sqrt_transf = TRUE)


## data_and_plot(W_specific_10_prunedb,
##               data_for_specific_weighted_10, 195825L,
##               main = "weighted-10-specific-no-split-01-sqrt-all_methods-best",
##               weighted = TRUE,
##               size_y_text = 8,
##               height = 22,
##               best_vars = all_vars,
##               sqrt_transf = TRUE)



## data_and_plot(W_all_7_prunedb,
##               data_for_all_weighted_7, 673883L,
##               main = "weighted-7-ALL-no-split-01-sqrt-all_methods-best",
##               weighted = TRUE,
##               size_y_text = 8,
##               height = 22,
##               best_vars = all_vars,
##               sqrt_transf = TRUE)


## data_and_plot(W_all_10_prunedb,
##               data_for_all_weighted_10, 1715676L,
##               main = "weighted-10-ALL-no-split-01-sqrt-all_methods-best",
##               weighted = TRUE,
##               size_y_text = 8,
##               height = 22,
##               best_vars = all_vars,
##               sqrt_transf = TRUE)









## data_and_plot(W_specific_7_prunedb,
##               data_for_specific_weighted_7, 76085L,
##               main = "weighted-7-specific-no-split-01-sqrt-js",
##               weighted = TRUE,
##               sqrt_transf = TRUE,
##               size_y_text = 10,
##               height = 22,
##               show_bottom = "js")






## data_and_plot(W_specific_7_prunedb,
##               data_for_specific_weighted_7, 76085L,
##               main = "weighted-7-specific-no-split-01-sqrt-all_methods-js",
##               weighted = TRUE,
##               size_y_text = 8,
##               height = 22,
##               best_vars = all_vars,
##               sqrt_transf  = TRUE, show_bottom = "js")







## data_and_plot(W_specific_7_prunedb,
##               data_for_specific_weighted_7, 76085L,
##               main = "weighted-7-specific-no-split-01-sqrt-all_methods-average-rank",
##               weighted = TRUE,
##               size_y_text = 8,
##               height = 22,
##               best_vars = all_vars,
##               sqrt_transf  = TRUE, show_bottom = "rank")


## data_and_plot(W_specific_10_prunedb,
##               data_for_specific_weighted_10, 195825L,
##               main = "weighted-10-specific-no-split-01-sqrt-all_methods-average-rank",
##               weighted = TRUE,
##               size_y_text = 8,
##               height = 22,
##               best_vars = all_vars,
##               sqrt_transf  = TRUE, show_bottom = "rank")


## data_and_plot(W_all_7_prunedb,
##               data_for_all_weighted_7, 673883L,
##               main = "weighted-7-ALL-no-split-01-sqrt-all_methods-average-rank",
##               weighted = TRUE,
##               size_y_text = 8,
##               height = 22,
##               best_vars = all_vars,
##               sqrt_transf  = TRUE, show_bottom = "rank")

## data_and_plot(W_all_10_prunedb,
##               data_for_all_weighted_10, 1715676L,
##               main = "weighted-10-ALL-no-split-01-sqrt-all_methods-average-rank",
##               weighted = TRUE,
##               size_y_text = 8,
##               height = 22,
##               best_vars = all_vars,
##               sqrt_transf  = TRUE, show_bottom = "rank")






## ## Show intercept (means of fit). No back transform the sqrt, therefore
## data_and_plot(W_specific_7_prunedb,
##               data_for_specific_weighted_7, 76085L,
##               main = "weighted-7-specific-no-split-01-sqrt-show-mean",
##               weighted = TRUE,
##               size_y_text = 10,
##               height = 22,
##               ## best_vars = all_vars,
##               sqrt_transf = FALSE,
##               show_mean = TRUE)


## data_and_plot(W_specific_10_prunedb,
##               data_for_specific_weighted_10, 195825L,
##               main = "weighted-10-specific-no-split-01-sqrt-show-mean",
##               weighted = TRUE,
##               size_y_text = 10,
##               height = 22,
##               ## best_vars = all_vars,
##               sqrt_transf = FALSE,
##               show_mean = TRUE)


## data_and_plot(W_all_7_prunedb,
##               data_for_all_weighted_7, 673883L,
##               main = "weighted-7-ALL-no-split-01-sqrt-show-mean",
##               weighted = TRUE,
##               size_y_text = 10,
##               height = 22,
##               sqrt_transf = FALSE,
##               show_mean = TRUE)


## data_and_plot(W_all_10_prunedb,
##               data_for_all_weighted_10, 1715676L,
##               main = "weighted-10-ALL-no-split-01-sqrt-show-mean",
##               weighted = TRUE,
##               size_y_text = 10,
##               height = 22,
##               sqrt_transf = FALSE,
##               show_mean = TRUE)







## From the ones in the original scale
## > scatterplot_median_intercept(W_specific_10_02$tree, TRUE)

## 	Pearson's product-moment correlation

## data:  medians and intercepts
## t = 11.07, df = 32, p-value = 1.782e-12
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.7902730 0.9442848
## sample estimates:
##       cor 
## 0.8904673 


## 	Spearman's rank correlation rho

## data:  medians and intercepts
## S = 692, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.8942704 


## Call:
## lm(formula = intercepts ~ medians)

## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.15597 -0.05499  0.00933  0.05382  0.12250 

## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.09874    0.02533   3.898 0.000466 ***
## medians      0.97253    0.08786  11.070 1.78e-12 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 0.07472 on 32 degrees of freedom
## Multiple R-squared:  0.7929,	Adjusted R-squared:  0.7865 
## F-statistic: 122.5 on 1 and 32 DF,  p-value: 1.782e-12

## > scatterplot_median_intercept(W_specific_10$tree, TRUE)

## 	Pearson's product-moment correlation

## data:  medians and intercepts
## t = 16.127, df = 67, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.8302592 0.9317471
## sample estimates:
##       cor 
## 0.8917127 


## 	Spearman's rank correlation rho

## data:  medians and intercepts
## S = 5492.6, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.8996611 


## Call:
## lm(formula = intercepts ~ medians)

## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.151664 -0.041985 -0.007124  0.040800  0.156907 

## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.06541    0.01804   3.627 0.000555 ***
## medians      1.05220    0.06525  16.127  < 2e-16 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 0.07379 on 67 degrees of freedom
## Multiple R-squared:  0.7952,	Adjusted R-squared:  0.7921 
## F-statistic: 260.1 on 1 and 67 DF,  p-value: < 2.2e-16

## Warning message:
## In cor.test.default(medians, intercepts, method = "spearman") :
##   Cannot compute exact p-value with ties
## > scatterplot_median_intercept(W_specific_7_02$tree, TRUE)

## 	Pearson's product-moment correlation

## data:  medians and intercepts
## t = 11.022, df = 32, p-value = 1.99e-12
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.7888256 0.9438669
## sample estimates:
##       cor 
## 0.8896688 


## 	Spearman's rank correlation rho

## data:  medians and intercepts
## S = 786.56, p-value = 7.26e-12
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.8798228 


## Call:
## lm(formula = intercepts ~ medians)

## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.09601 -0.05305  0.00749  0.04678  0.14378 

## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.07165    0.02031   3.528  0.00129 ** 
## medians      1.06140    0.09630  11.022 1.99e-12 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 0.06362 on 32 degrees of freedom
## Multiple R-squared:  0.7915,	Adjusted R-squared:  0.785 
## F-statistic: 121.5 on 1 and 32 DF,  p-value: 1.99e-12

## Warning message:
## In cor.test.default(medians, intercepts, method = "spearman") :
##   Cannot compute exact p-value with ties
## > scatterplot_median_intercept(W_specific_7$tree, TRUE)

## 	Pearson's product-moment correlation

## data:  medians and intercepts
## t = 16.82, df = 55, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.8593166 0.9492516
## sample estimates:
##       cor 
## 0.9150053 


## 	Spearman's rank correlation rho

## data:  medians and intercepts
## S = 2781.5, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##      rho 
## 0.909854 


## Call:
## lm(formula = intercepts ~ medians)

## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.17145 -0.03454  0.00358  0.03562  0.15518 

## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.06559    0.01461    4.49 3.69e-05 ***
## medians      1.07062    0.06365   16.82  < 2e-16 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 0.06281 on 55 degrees of freedom
## Multiple R-squared:  0.8372,	Adjusted R-squared:  0.8343 
## F-statistic: 282.9 on 1 and 55 DF,  p-value: < 2.2e-16

## Warning message:
## In cor.test.default(medians, intercepts, method = "spearman") :
##   Cannot compute exact p-value with ties
## > 
## > scatterplot_median_intercept(W_all_10_02$tree, TRUE)

## 	Pearson's product-moment correlation

## data:  medians and intercepts
## t = 7.2533, df = 35, p-value = 1.805e-08
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.6020806 0.8783880
## sample estimates:
##       cor 
## 0.7749233 


## 	Spearman's rank correlation rho

## data:  medians and intercepts
## S = 1976, p-value = 3.69e-07
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.7657658 


## Call:
## lm(formula = intercepts ~ medians)

## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.25644 -0.05665  0.01678  0.06038  0.19314 

## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.28352    0.03018   9.394 4.24e-11 ***
## medians      0.58348    0.08044   7.253 1.80e-08 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 0.09694 on 35 degrees of freedom
## Multiple R-squared:  0.6005,	Adjusted R-squared:  0.5891 
## F-statistic: 52.61 on 1 and 35 DF,  p-value: 1.805e-08

## > scatterplot_median_intercept(W_all_10$tree, TRUE)

## 	Pearson's product-moment correlation

## data:  medians and intercepts
## t = 18.446, df = 66, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.8656290 0.9469695
## sample estimates:
##       cor 
## 0.9151735 


## 	Spearman's rank correlation rho

## data:  medians and intercepts
## S = 5202.5, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.9007033 


## Call:
## lm(formula = intercepts ~ medians)

## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.17255 -0.02865  0.01019  0.03493  0.12189 

## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.26309    0.01274   20.65   <2e-16 ***
## medians      0.63014    0.03416   18.45   <2e-16 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 0.05629 on 66 degrees of freedom
## Multiple R-squared:  0.8375,	Adjusted R-squared:  0.8351 
## F-statistic: 340.3 on 1 and 66 DF,  p-value: < 2.2e-16

## Warning message:
## In cor.test.default(medians, intercepts, method = "spearman") :
##   Cannot compute exact p-value with ties
## > 
## > scatterplot_median_intercept(W_all_7_02$tree, TRUE)

## 	Pearson's product-moment correlation

## data:  medians and intercepts
## t = 9.6173, df = 31, p-value = 8.071e-11
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.7428233 0.9318718
## sample estimates:
##       cor 
## 0.8654317 


## 	Spearman's rank correlation rho

## data:  medians and intercepts
## S = 412, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.9311497 


## Call:
## lm(formula = intercepts ~ medians)

## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.18506 -0.05520 -0.01038  0.04002  0.17463 

## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.12199    0.02560   4.765 4.21e-05 ***
## medians      0.89603    0.09317   9.617 8.07e-11 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 0.08083 on 31 degrees of freedom
## Multiple R-squared:  0.749,	Adjusted R-squared:  0.7409 
## F-statistic: 92.49 on 1 and 31 DF,  p-value: 8.071e-11

## > scatterplot_median_intercept(W_all_7$tree, TRUE)

## 	Pearson's product-moment correlation

## data:  medians and intercepts
## t = 10.226, df = 62, p-value = 6.257e-15
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.678690 0.868901
## sample estimates:
##       cor 
## 0.7923182 


## 	Spearman's rank correlation rho

## data:  medians and intercepts
## S = 8255.6, p-value = 4.544e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.8109983 


## Call:
## lm(formula = intercepts ~ medians)

## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.19926 -0.07593  0.01140  0.07006  0.17493 

## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.12837    0.02042   6.288 3.61e-08 ***
## medians      0.77999    0.07628  10.226 6.26e-15 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 0.08992 on 62 degrees of freedom
## Multiple R-squared:  0.6278,	Adjusted R-squared:  0.6218 
## F-statistic: 104.6 on 1 and 62 DF,  p-value: 6.257e-15

## Warning message:
## In cor.test.default(medians, intercepts, method = "spearman") :
##   Cannot compute exact p-value with ties
## > 
## > 
## > 
## rm(list = ls(pattern = "_02$"))



