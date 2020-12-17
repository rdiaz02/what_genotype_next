## Runs the glmertree over genotypes for the min_js variable, single for
## all fitness landscapes

## Here, I split the 7 and 10 gene fits, since the meaning of the number
## of genes variable is different.


date()
print(Sys.getpid())

source("common_lmertree_fits.R")

RhpcBLASctl::blas_set_num_threads(4)
RhpcBLASctl::omp_set_num_threads(4)




## Models for min_js: plotted as
## partykit:::plot.lmtree(fit$tree, tp_args = list(id = FALSE), gp )

## use plot_best in common_lmertree_fits.R



## I rename same variables for easier to read trees
## Create also new variables for the crossed random effects
## for easier specification
create_local_data <- function(data) {
    print(nrow(data))
    localdata <- data
    cat("\n post sampledProp trimming\n ")
    print(nrow(localdata))
    setnames(localdata,
             old = c("sourceGenotype_nMut",
                     "size_split",
                     "freq_most_freq_mean_no450"),
             new = c("nMut",
                     "sample_size",
                     "SSWM"))
    ## localdata[, `:=` (
    ##     idgenot = paste0(id, "_", sourceGenotype),
    ##     idrep   = paste0(id, "_", replicate))]
    return(localdata)
}


rm(localformula)
rm(localdata)

## created in create_wide_g_long_g.R
load("../../wide_g.RData")


print(nrow(wide_g))
localdata <- create_local_data(wide_g)
rm(wide_g)
gc()

localdata[, `:=`( typeLandscape = typeLandscape_f)]
levels(localdata$typeLandscape_f) <- c("Rep", "LM", "RMF")

with(localdata, table(typeLandscape_f, typeLandscape))

localdata[, `:=`( numGenes = as.factor(numGenes))]

localdata <- localdata[ (numGenes == 10) 
                     , ]

cat("\n nrow(localdata) ")
nrow(localdata)

## So as to get the multiplication factor right,
## do here the additional trimming of data

localdata <- localdata[
    !is.na(min_js) 
    & (
        ((numGenes == 10) & (nMut <= 9))
    )
  ,
]

cat("\n nrow(localdata) post-trimming before mfw ")
nrow(localdata)


mfw <- nrow(localdata)/sum(localdata$sampledProp)
cat("\n Multiplication factor for weights is ", mfw, "\n")
localdata[, `:=`( sampledProp = mfw * sampledProp)]

localdata <- localdata[, `:=`(min_js = sqrt(min_js))]

## Remember we convert nMut and sample_size into factors
## inside tree_best

localformula <- formula(min_js ~ 1 |
                            (1 | id) +
                            (1 | id:sourceGenotype) +
                            (1 | id:replicate)
                        |
                        nMut ##  sourceGenotype_nMut
                        + detect
                        + sample_size ## + size_split
                        ## + numGenes
                        + SSWM ##  + freq_most_freq_mean_no450
                        + gamma
                        + epistRSign
                        + numObservedPeaks
                         + diff_obs_prop
                         + observedProp
                        + fitnessRank
                        + propLocalMax)


tree_best_time_print_save("allFL", Inf, "lmertree_sqrt_W_Obs_minsize_01_fLandscape-rsign-peaks-recoded_10genes_Inf_BBY",
                          lmer.control = lmerControl(optimizer = "bobyqa",
                                                     optCtrl=list(maxfun=5e5))
                          , weights = TRUE,minsizef = 0.01)

print(warnings())


tree_best_time_print_save("allFL", Inf, "lmertree_sqrt_W_Obs_minsize_01_fLandscape-rsign-peaks-recoded_10genes_Inf_NM",
                          lmer.control = lmerControl(optimizer = "Nelder_Mead")                          
                          , weights = TRUE,minsizef = 0.01)

print(warnings())


tree_best_time_print_save("allFL", Inf, "lmertree_sqrt_W_Obs_minsize_01_fLandscape-rsign-peaks-recoded_10genes_Inf_default"
                          , weights = TRUE,minsizef = 0.01)

print(warnings())


tree_best_time_print_save("allFL", Inf, "lmertree_sqrt_W_Obs_minsize_01_fLandscape-rsign-peaks-recoded_10genes_Inf_NLM",
                          lmer.control = lmerControl(optimizer = "optimx",
                                                     optCtrl=list(method= "nlminb"))
                          , weights = TRUE,minsizef = 0.01)


print(warnings())


tree_best_time_print_save("allFL", Inf, "lmertree_sqrt_W_Obs_minsize_01_fLandscape-rsign-peaks-recoded_10genes_Inf_LBF",
                          lmer.control = lmerControl(optimizer = "optimx",
                                                     optCtrl=list(method = "L-BFGS-B"))
                          , weights = TRUE,minsizef = 0.01)

print(warnings())

date()

