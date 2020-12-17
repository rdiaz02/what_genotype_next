## Prepare data for plots of the trees
source("common_lmertree_fits.R")

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


load("../wide_g.RData")


print(nrow(wide_g))
localdata <- create_local_data(wide_g)
## rm(wide_g)
gc()

## Filter the max_g_7 and max_g_10
max_g_7 <- 6
max_g_10 <- 9

localdata <- localdata[
    ((numGenes == 7) & (nMut <= max_g_7))
    |
    ((numGenes == 10) & (nMut <= max_g_10)),
]


## This is wide_g without the genotypes with 7 and 10 mutations when
## number of loci are 7 and 10, respectively. But we still have NAs in
## min_js.
print(nrow(localdata))

## From create_wide_g_long_g.R
## This contains the rank of each method, the actual value of each method, etc
load("../wide_genotype_no_average_rank.RData")



wide_genotype_no_average_rank[,
                              `:=`(
                                  which_single_best_js_f = factor(which_single_best_js),
                                  typeLandscape_f =
                                      factor(typeLandscape,
                                             levels = c("Represent.", "Local maxima", "RMF")),
                                  detect = factor(detect),
                                  id = factor(id),
                                  sourceGenotype = factor(sourceGenotype)
                                  ## accessible is a method property!!
                                  ## sourceGenotype_accessible = as.integer(sourceGenotype_accessible)
                              )]

min(wide_genotype_no_average_rank$propLocalMax, na.rm = TRUE)
min(wide_genotype_no_average_rank$observedProp, na.rm = TRUE)
min(wide_genotype_no_average_rank$sampledProp, na.rm = TRUE)

summary(wide_genotype_no_average_rank[, list(propLocalMax, observedProp, sampledProp)])

## These NAs should be 0
wide_genotype_no_average_rank[is.na(propLocalMax), propLocalMax := 0]
wide_genotype_no_average_rank[is.na(observedProp), observedProp := 0]
wide_genotype_no_average_rank[is.na(sampledProp), sampledProp := 0]

wide_genotype_no_average_rank[, diff_obs_prop := sampledProp - observedProp]

## Yes, now needed as thee are still cases with sampledProp == 0
nrow(wide_genotype_no_average_rank)
wide_genotype_no_average_rank <- wide_genotype_no_average_rank[sampledProp > 0, ]
nrow(wide_genotype_no_average_rank)

## We do this to rename some variables. Note nrow does not change
localdata2 <- create_local_data(wide_genotype_no_average_rank)

## Again, filter the max_g_7 and max_g_10
localdata2 <- localdata2[
    ((numGenes == 7) & (nMut <= max_g_7))
    |
    ((numGenes == 10) & (nMut <= max_g_10)),
]

print(nrow(localdata2))

## same as
print(nrow(localdata))
stopifnot(nrow(localdata)== nrow(localdata2))

## The two data sets contain the same values in the exact same order.

## We proceed this way because we will extract a matrix where columns are
## all methods, with a 1 if that method was the best (matrix.best,
## below). This is exctracted from "all_best_js" column in localdata2. In
## fact, localdata (wide_g) is never used for anything, except for these
## paranoid checks that we have selected the exact same cases.
stopifnot(identical(localdata$min_js, localdata2$min_js))
stopifnot(identical(localdata$id, localdata2$id))
stopifnot(identical(localdata$sourceGenotype, localdata2$sourceGenotype))
stopifnot(identical(localdata$detect, localdata2$detect))
stopifnot(identical(localdata$sample_size, localdata2$sample_size))
stopifnot(identical(localdata$replicate, localdata2$replicate))

## A few values, for the hell of showing that we get some with multiple best
localdata2$all_best_js[c(1, 1000, 2000, 10000, nrow(localdata2))]

## used to check a check below
uu <- localdata2$all_best_js[c(1, 1000, 2000, 10000, nrow(localdata2))]

## Brutal, but a check too: just look at levels: all methods should be there
levels <- sort(unique(unlist(strsplit(localdata2$all_best_js, ", "))))
names(levels) <- levels

among_best <- function(x, thelevels = levels) {
    tmp <- as.integer(thelevels %in% unlist(strsplit(x, split = ", ")))
    names(tmp) <- names(thelevels)
    return(tmp)
}

t(vapply(uu, FUN = function(z) among_best(z, thelevels = levels),
       FUN.VALUE = integer(length(levels)),
       USE.NAMES = TRUE))


matrix.best <- t(vapply(localdata2$all_best_js,
                      FUN = function(z) among_best(z, thelevels = levels),
                      FUN.VALUE = integer(length(levels)),
                      USE.NAMES = TRUE))

dim(matrix.best)
row.names(matrix.best) <- NULL

## Prepare for merging, which is just like a cbind
matrix.best <- data.table(matrix.best)
matrix.best[, rownum := .I]
localdata2[, rownum := .I]

colnames(matrix.best) <-
    gsub("_js", "_is_best", colnames(matrix.best))


localdata3 <- merge(localdata2, matrix.best, by = "rownum")

localdata3 <- localdata3[!is.na(min_js), ]

## Check: rsmb comes from matrix.best, number_best_js was in localdata2
rsmb <- rowSums(matrix.best[, !"rownum"])
rsmb <- rsmb[rsmb > 0]

stopifnot( (rsmb - localdata3$number_best_js) == 0)

data_for_weighted_glmertree_plots <- localdata3

print(nrow(localdata3))

save(file = "data_for_weighted_glmertree_plots.RData",
     data_for_weighted_glmertree_plots,
     compress = FALSE)
