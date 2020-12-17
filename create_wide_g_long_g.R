## As the name says. Minimal, concise and with factors,
## for the models

library(data.table)
setDTthreads(threads = 0)



long_genot <- data.table(readRDS("data_no_any_no_null.rds")$data)

long_genot[,
           `:=`(
               typeLandscape_f =
                   factor(typeLandscape,
                          levels = c("Represent.", "Local maxima", "RMF")),
               detect = factor(detect),
               id = factor(id),
               method = factor(method),
               sourceGenotype = factor(sourceGenotype),
               sourceGenotype_accessible = as.integer(sourceGenotype_accessible)
           )]

## changes:
## propLocalMax: all NAs are 0
## observedProp: all NAs are 0
## sampledProp:  all NAs are 0

min(long_genot$propLocalMax, na.rm = TRUE)
min(long_genot$observedProp, na.rm = TRUE)
min(long_genot$sampledProp, na.rm = TRUE)

## sourceGenotype_accessible: this is irrelevant:
with(long_genot, any(is.na(sourceGenotype_accessible) & !(is.na(js))))
## but why NA? because it is NA, as these are failures


summary(long_genot[, list(propLocalMax, observedProp, sampledProp)])

long_genot[is.na(propLocalMax), propLocalMax := 0]
long_genot[is.na(observedProp), observedProp := 0]
long_genot[is.na(sampledProp), sampledProp := 0]

long_genot[, diff_obs_prop := sampledProp - observedProp]

long_g <- long_genot[
    ## !is.na(min_js)
   ,
    list(
        js
      , method
      , id
      , sourceGenotype        
      , detect
      , sourceGenotype_nMut
      , size_split, freq_most_freq_mean_no450
      , numObservedPeaks, gamma, numGenes, typeLandscape_f
      , w2, w3, epistRSign
      , lod_h
      , diff_obs_prop
      , observedProp
      , fitnessRank
      , fitnessRankNoZero
      , propLocalMax
      , sampledGenotypesNumber
      , numMuts_mean
      , sourceGenotype_accessible
      , sampledProp
      , replicate
    )
]

summary(long_g)

## created in create_rank_data.R
load(file = "wide_genotype_no_average_rank.RData")

## simpler for models below
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

wide_genotype_no_average_rank[is.na(propLocalMax), propLocalMax := 0]
wide_genotype_no_average_rank[is.na(observedProp), observedProp := 0]
wide_genotype_no_average_rank[is.na(sampledProp), sampledProp := 0]

wide_genotype_no_average_rank[, diff_obs_prop := sampledProp - observedProp]

wide_g <- wide_genotype_no_average_rank[
   ,
    list(
        min_js
      , id
      , sourceGenotype
      , detect
      , sourceGenotype_nMut
      , size_split, freq_most_freq_mean_no450
      , numObservedPeaks, gamma, numGenes, typeLandscape_f
      , w2, w3, epistRSign
      , lod_h
      , diff_obs_prop
      , observedProp
      , fitnessRank
      , fitnessRankNoZero
      , propLocalMax
      , sampledGenotypesNumber
      , numMuts_mean
        ##      , sourceGenotype_accessible ## NOPE: this is a method property
      , sampledProp
      , replicate
    )
]

## Only those present
long_g <- long_g[sampledProp > 0, ]
wide_g <- wide_g[sampledProp > 0, ]

save(file = "long_g.RData", long_g, compress = FALSE)
save(file = "wide_g.RData", wide_g, compress = FALSE)



