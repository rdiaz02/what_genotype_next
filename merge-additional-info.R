date()
cat("Replicates version: Combining data into single table\n")

#### Change paths as needed

rm(list = ls())
gc()

load("pre-table-replicates.RData") ## takes a couple of minutes
flDirectory <- "/home/jdc/next-genotype/fitland"
load(file.path(flDirectory, "fl_sampl.RData"))
load("outObsFreqs.RData")
load("outRanksFitness.RData")
load("weightsAll.RData")
load("outLocalMax.RData")


source("oncoFunctions.R")


df_r <- dplyr::bind_rows(df)

cat("Replicates version: formatting and columns\n")

## change "" genotype into "WT" or "root" and set accessibility to TRUE
df_r$sourceGenotype[df_r$sourceGenotype==""] <- "root"

## format "size_split" column
df_r$size_split[df_r$size_split=="size_split_50"] <- "50"
df_r$size_split[df_r$size_split=="size_split_200"] <- "200"
df_r$size_split[df_r$size_split=="size_split_4000"] <- "4000"
df_r$size_split <- as.numeric(df_r$size_split)

## format "detect" column
df_r$detect[df_r$detect=="detect_large"] <- "large"
df_r$detect[df_r$detect=="detect_small"] <- "small"
df_r$detect[df_r$detect=="detect_unif"] <- "uniform"

######## Merging sampling characteristics
fl_sampl$Detection <- as.character(fl_sampl$Detection)
fl_sampl$Detection[fl_sampl$Detection=="unif"] <- "uniform"


#### To compare merge and left_join
compare_merge_left_join <- FALSE

if(compare_merge_left_join) {
    system.time(
        df_r_merge <- merge(df_r, fl_sampl,
                            by.x=c("id","detect"),by.y=c("ID","Detection"),
                            all=T,suffixes=NULL)
    ) #### 164 s
}

#### This is 10x to 16x faster
system.time(
    df_r_join <- dplyr::left_join(df_r, fl_sampl,
                                  by = c("id" = "ID", "detect" = "Detection"))
) 

#### Checks
##   Null model not under any sampling. We will change this below.
stopifnot(all(is.na(unique(df_r$detect[df_r$cpm == "null"]))))
##   The rest are what they should be
stopifnot(unique(df_r$detect[df_r$cpm != "null"]) == c("large", "small", "uniform"))

stopifnot(all(is.na(unique(df_r_join$detect[df_r_join$cpm == "null"]))))
stopifnot(unique(df_r_join$detect[df_r_join$cpm != "null"]) == c("large", "small", "uniform"))

stopifnot(nrow(df_r) == nrow(df_r_join))


if(compare_merge_left_join) {
    #### checks
    stopifnot(identical(nrow(df_r_join), nrow(df_r_merge)))
    stopifnot(identical(ncol(df_r_join), ncol(df_r_merge)))
    #### order of colnames changes, but columns are the same
    stopifnot(all(colnames(df_r_join) %in% colnames(df_r_merge)))
    stopifnot(all(colnames(df_r_merge) %in% colnames(df_r_join)))

    df_r_join <- df_r_join[, colnames(df_r_merge)]

    #### still, not identical
    which_no_i <-
        sapply(1:ncol(df_r_merge),
               function(i) !identical(df_r_merge[, i], df_r_join[, i]))
    #### what gives here?
    #### ordering differs, starting with row 191
    df_r_merge[190:194, c(1:4)]
    df_r_join[190:194, c(1:4)]

    #### Lets fix that to compare
    #### (yes, this resorting is slow)
    df_r_m_2 <- df_r_merge[order(df_r_merge$id,
                                 df_r_merge$detect,
                                 df_r_merge$size_split,
                                 df_r_merge$cpm,
                                 df_r_merge$replicate,
                                 df_r_merge$sourceGenotype),
                           ]

    df_r_j_2 <- df_r_join[order(df_r_join$id,
                                df_r_join$detect,
                                df_r_join$size_split,
                                df_r_join$cpm,
                                df_r_join$replicate,
                                df_r_join$sourceGenotype),
                          ]

    #### This apparently fails
    identical(df_r_j_2, df_r_m_2)

    #### But this does not
    which_no_i_2 <-
        sapply(1:ncol(df_r_m_2),
               function(i) !identical(df_r_m_2[, i], df_r_j_2[, i]))

    #### As usual, the culprit are attributes
    identical(attributes(df_r_m_2)$row.names,
              attributes(df_r_j_2)$row.names)
    #### And that happened in one of join/merge
    #### (this should not be identical, unless one or both set reset
    #### rownames after join/merge)
    head(row.names(df_r_join)) 
    head(row.names(df_r_merge))
    #### Reset row names
    row.names(df_r_m_2) <- NULL
    row.names(df_r_j_2) <- NULL


    #### And now, both are identical
    stopifnot(identical(df_r_j_2, df_r_m_2))

    #### game over for checks
    rm(df_r_merge, df_r_m_2, df_r_j_2, which_no_i, which_no_i_2)
    gc()
}

#### df_r_join left around helps to backtrack errors
df_r <- df_r_join

#### Adding weights (computed externally with makeTableWeights.R)
#### recall to change unif to uniform

nrow(weightsAll)
weightsAll$detect[weightsAll$detect == "unif"] <- "uniform"
weightsAll$Genotype[weightsAll$Genotype == "WT"] <- "root"

## The null can break the left_joins in the sense that it
## has not detect (is NA)

## That is irrelevant for most data (and was irrelevant in the merge
## before for fl characters). But we must have entries for the weights of
## genotypes according to "true" and "sampled". True is the same
## regardless of sampling, but sample is not. We could add True weights to
## null easily, but not sampled, as we need different weights for each
## case. So replicate null.


nullp <- dplyr::filter(df_r, cpm == "null")
df_r_nn <- dplyr::filter(df_r, cpm != "null")

nullpL <- nullpS <- nullpU <- nullp
nullpL$detect <- "large"
nullpS$detect <- "small"
nullpU$detect <- "uniform"
nullpA <- rbind(nullpL, nullpS, nullpU)
stopifnot(length(table(nullpA$detec, useNA = "ifany")) == 3)
stopifnot(table(nullpA$detec, useNA = "ifany") == nrow(nullp))
df_r2 <- rbind(df_r_nn, nullpA)
stopifnot(nrow(df_r2) == ( sum(df_r$cpm != "null") + 3 * sum(df_r$cpm == "null") ))



df_r <- df_r2
rm(df_r2)
gc()



(nr_df_r_before_join <- nrow(df_r))
df_r <- dplyr::left_join(df_r, weightsAll,
                         by = c("id" = "ID",
                                "detect",
                                "sourceGenotype" = "Genotype"))

nrow(df_r) #### same as before

stopifnot(nrow(df_r) == nr_df_r_before_join)


############## Anything in POM not in True? Yes, but ignore ##############
####  (how? in C++ code POM is updated even if no creation of
####   output table because it happens in between keepEverys)
df_r_no_any <- dplyr::filter(df_r, sourceGenotype != "any")
df_r_no_any <- dplyr::filter(df_r_no_any, cpm != "null")

not_in_true <- which(is.na(df_r_no_any$TrueProp))
#### There are some
length(not_in_true) ## 56745


#### but very tiny frequencies : 3 simulations at most
summary(df_r_no_any$sourceGenotype_freqInPOM[not_in_true])
max(df_r_no_any$sourceGenotype_freqInPOM[not_in_true])
20000 * max(df_r_no_any$sourceGenotype_freqInPOM[not_in_true])

#### And in 95% of the cases in just one sample
20000 * quantile(
            df_r_no_any$sourceGenotype_freqInPOM[not_in_true],
            probs = c(0.95, 0.96, 0.97, 0.99))

#### Affects 82 IDs
length(unique(df_r_no_any$id[not_in_true]))

## Do not break anything on reruns
stopifnot(length(unique(df_r_no_any$id[not_in_true])) == 82)
stopifnot(length(not_in_true) == 56745)
stopifnot(all.equal(20000 * max(df_r_no_any$sourceGenotype_freqInPOM[not_in_true]) , 3))


#### Check if anything in True not in POM
####    This should not happen, as True is a subset of POM

nrow(weightsAll)
nrow(df_r_join)

df_tnt <- dplyr::left_join(weightsAll, df_r_join,
                         by = c("ID" = "id",
                                "detect",
                                "Genotype" = "sourceGenotype"))
nrow(df_tnt)

not_in_pom <- which(is.na(df_tnt$sourceGenotype_freqInPOM))
stopifnot(length(not_in_pom) == 0) #### this is 0

#### This true, but this does not verify anything about not in POMs
not_in_pomb <- which(is.na(df_r_no_any$sourceGenotype_freqInPOM))
stopifnot(length(not_in_pomb) == 0) #### this is 0

rm(df_tnt)
rm(df_r_no_any)
gc()

##### Observed frequencies in sampling, fitness ranks, local max
outObsFreqs$Genotype[outObsFreqs$Genotype == "WT"] <- "root"
outObsFreqs$detect[outObsFreqs$detect == "unif"] <- "uniform"
outRanksFitness$Genotype[outRanksFitness$Genotype == "WT"] <- "root"
outLocalMax$Genotype[outLocalMax$Genotype == "WT"] <- "root"


nra <- nrow(df_r)
df_r <- dplyr::left_join(df_r, outRanksFitness,
                         by = c("id" = "ID",
                                "sourceGenotype" = "Genotype"))
stopifnot(nrow(df_r) == nra)


####  Minimal verifications. Should all be none or 0
stopifnot(length(setdiff(unique(outObsFreqs$ID), unique(df_r$id))) == 0)
stopifnot(length(setdiff(unique(outObsFreqs$detect), unique(df_r$detect))) == 0)
stopifnot(length(setdiff(unique(outObsFreqs$replicate), unique(df_r$replicate))) == 0)
stopifnot(length(setdiff(unique(outObsFreqs$Genotype), unique(df_r$sourceGenotype))) == 0)
stopifnot(length(setdiff(unique(outObsFreqs$size_split), unique(df_r$size_split))) == 0)


df_r <- dplyr::left_join(df_r,
                          outObsFreqs,
                          by = c("id" = "ID",
                                 "detect",
                                 "replicate",
                                 "size_split",
                                 "sourceGenotype" = "Genotype"))

stopifnot(nrow(df_r) == nra)

df_r <- dplyr::left_join(df_r,
                          outLocalMax,
                          by = c("id" = "ID",
                                 "sourceGenotype" = "Genotype"))

stopifnot(nrow(df_r) == nra)

rm(df_r_join)
gc()




#### Can be handy when dealing with replicates
#### and, yes, arrange much faster than [order(...
#### About 25 seconds
system.time(df_r <- dplyr::arrange(df_r,
                                   id,
                                   detect,
                                   size_split,
                                   cpm,
                                   sourceGenotype,
                                   replicate
                                   ))


## Avoid this: dangerous (could exchange column contents) and forces a
## copy of complete data.frame, and unclear (which four columns of df_r
## are we dropping?)

## setdiff(colnames(df_r),
##                    c(
##         "id",
##         "cpm", "size_split",
##         "sqDiff", "sqDiff_eq", "js", "js_eq", "hellinger",
##         "hellinger_eq", "spearman", "spearman_pval",
##         "flags",
##         "sourceGenotype", "sourceGenotype_nMut",
##         "sourceGenotype_freqInPOM","sourceGenotype_accessible",
##         "detect", "sampledGenotypesDiversity", "sampledGenotypesNumber",
##         "sampledGenotypesGenesAbove_0.01", "sampledGenotypesGenesAbove_0.1",
##         "Mean_muts", "Median_muts", "Var_muts",
##         "Stdev_muts", "Kurtosis_muts","Skewness_muts",
##         "pom_h", "lod_h",
##         "freq_most_freq_mean_no450", "how_many_gt_5p_mean_no450",
##         "num_accessible_genots", "Init_Size", "Mutation",
##         "rnst", "nGenes", "typeLandscape",
##         "num_local_peaks", "num_observed_peaks",
##         "epist_magn", "epist_sign", "epist_rsign",
##         "w.1.", "w.2.", "w.3..","gamma",
##         "replicate",
##         "SampledFreq", "SampledProp", "TrueFreq", "TrueProp",
##         "fitnessRank", "fitnessRankNoZero", "ObservedFreq", "ObservedProp",
##         "FreqLocalMax", "PropLocalMax"
##     ))

## But type_Landscape is handy for plots, as more explicit
## table(df_r$typeLandscape, df_r$type_Landscape)
## unique(df_r$typeLandscape)
## We will rm AFTER checking all OK (below) as we rename levels

## ## and diversity_observed_peaks is a measure of diversity of evol
## with(df_r, table(detect, Sampling, useNA = "always"))
## with(df_r, table(nGenes, Ngenes, useNA = "always"))
## unique(df_r$Ngenes)
## ## We will drop Sampling, Ngenes

## df_r <-
##     df_r[, c(
##         "id",
##         "cpm", "size_split",
##         "sqDiff", "sqDiff_eq", "js", "js_eq", "hellinger",
##         "hellinger_eq", "spearman", "spearman_pval",
##         "flags",
##         "sourceGenotype", "sourceGenotype_nMut",
##         "sourceGenotype_freqInPOM","sourceGenotype_accessible",
##         "detect", "sampledGenotypesDiversity", "sampledGenotypesNumber",
##         "sampledGenotypesGenesAbove_0.01", "sampledGenotypesGenesAbove_0.1",
##         "Mean_muts", "Median_muts", "Var_muts",
##         "Stdev_muts", "Kurtosis_muts","Skewness_muts",
##         "pom_h", "lod_h",
##         "freq_most_freq_mean_no450", "how_many_gt_5p_mean_no450",
##         "num_accessible_genots", "Init_Size", "Mutation",
##         "rnst", "nGenes", "typeLandscape",
##         "num_local_peaks", "num_observed_peaks",
##         "epist_magn", "epist_sign", "epist_rsign",
##         "w.1.", "w.2.", "w.3..","gamma",
##         "replicate",
##         "SampledFreq", "SampledProp", "TrueFreq", "TrueProp",
##         "fitnessRank", "fitnessRankNoZero", "ObservedFreq", "ObservedProp",
##         "FreqLocalMax", "PropLocalMax"
##     )]

## colnames(df_r) <-
##     c(
##         "id",
##         "method","size_split",
##         "sqDiff", "sqDiff_eq", "js", "js_eq", "hellinger",
##         "hellinger_eq", "spearman", "spearman_pval",
##         "flags",
##         "sourceGenotype","sourceGenotype_numMuts",
##         "sourceGenotype_freqInPOM","sourceGenotype_accessible",
##         "detect", "sampledGenotypesDiversity", "sampledGenotypesNumber",
##         "sampledGenotypesGenesAbove_0.01", "sampledGenotypesGenesAbove_0.1",
##         "numMuts_mean", "numMuts_median", "numMuts_var",
##         "numMuts_stDev", "numMuts_kurtosis","nMuts_skewness",
##         "pom_h", "lod_h",
##         "freq_most_freq_mean_no450", "how_many_gt_5p_mean_no450",
##         "numAccessibleGenotypes", "initSize", "mutationRate",
##         "rnst", "numGenes", "typeLandscape",
##         "numLocalPeaks", "numObservedPeaks",
##         "epistMagn", "epistSign", "epistRSign",
##         "w1", "w2", "w3","gamma",
##         "replicate",
##         "sampledFreq", "sampledProp", "trueFreq", "trueProp",
##         "fitnessRank", "fitnessRankNoZero", "observedFreq", "observedProp",
##         "freqLocalMax", "propLocalMax"
##     )


## checks typeLandscape type_Landscape equivalence
df_r_tl_na <- dplyr::filter(df_r, is.na(type_Landscape))
stopifnot(length(df_r_tl_na$cpm) > 10000) ## make sure something
stopifnot(unique(df_r_tl_na$cpm) == "null")
df_r_tl_nonull <- dplyr::filter(df_r, cpm != "null")
stopifnot(any(is.na(df_r$type_Landscape)))
stopifnot(!any(is.na(df_r_tl_nonull$type_Landscape)))

tlnn <- with(df_r_tl_nonull,
             table(typeLandscape, as.character(type_Landscape), useNA = "ifany"))
stopifnot(dim(tlnn) == c(3, 3))
stopifnot(sum(diag(tlnn)) == nrow(df_r_tl_nonull))
## in contrast
tldfo <- with(df_r,
              table(typeLandscape, as.character(type_Landscape), useNA = "ifany"))
stopifnot(dim(tldfo) == c(3, 4))
stopifnot(sum(diag(tldfo)) < nrow(df_r))
stopifnot( (sum(diag(tldfo)) + sum(tldfo[, 4])) == nrow(df_r))
rm(tldfo, tlnn)


## A check about replicates
stopifnot(!any(is.na(df_r_tl_nonull$replicate)))
stopifnot(all(is.na(df_r[df_r$cpm == "null", ]$replicate)))
## yes, some: the nulls
stopifnot(any(is.na(df_r[df_r$sourceGenotype == "any", ]$replicate)))
stopifnot(!any(is.na(df_r[(df_r$sourceGenotype == "any") & (df_r$cpm != "null"), ]$replicate)))

rm(df_r_tl_na)
rm(df_r_tl_nonull)

df_r <- df_r[, !(names(df_r) %in% c("Sampling", "Ngenes", "type_Landscape"))]


## dplyr::rename seems simpler than data.table::set.names, as we have pairs old = new
## I am using quoted arguments
df_r <- dplyr::rename(df_r,
                      "method" =           "cpm",
                      "numMuts_mean" =     "Mean_muts",
                      "numMuts_median" =   "Median_muts",
                      "numMuts_var" =      "Var_muts",
                      "numMuts_stDev" =    "Stdev_muts",
                      "numMuts_kurtosis" = "Kurtosis_muts",
                      "numMuts_skewness" =   "Skewness_muts",
                      "initSize" =         "Init_Size",
                      "mutationRate" =     "Mutation",
                      "numGenes" =         "nGenes",
                      "numLocalPeaks" =    "num_local_peaks",
                      "numObservedPeaks" = "num_observed_peaks",
                      "numAccessibleGenotypes" = "num_accessible_genots",
                      "diversityObservedPeaks" = "diversity_observed_peaks",
                      "epistMagn" =        "epist_magn",
                      "epistSign" =        "epist_sign",
                      "epistRSign" =       "epist_rsign",
                      "w1" =               "w.1.",
                      "w2" =               "w.2.",
                      "w3" =               "w.3..",
                      "sampledFreq" =      "SampledFreq",
                      "sampledProp" =      "SampledProp",
                      "trueFreq" =         "TrueFreq",
                      "trueProp" =         "TrueProp",
                      "observedFreq" =     "ObservedFreq",
                      "observedProp" =     "ObservedProp",
                      "freqLocalMax" =     "FreqLocalMax",
                      "propLocalMax" =     "PropLocalMax"
)


#### explanation of each column
columnsExplained <- list(
    identification=list(
        id="Data identifier",
        replicate = "Replicate number"
    ),
    cpm_and_input=list(
        method="Method used for next genotype prediction",
        size_split="Sample size: size of the input matrix given to the method"
    ),
    statistics=list(
        sqDiff="Square root of average of squared differences",
        sqDiff_eq="Square root of average of squared differences (equiprobabilized)",
        js="Jensen-Shannon distance (square root of the Jensen-Shannon divergence), in log base 2 units",
        js_eq="Jensen-Shannon distance (equiprobabilized)",
        hellinger="Hellinger distance",
        hellinger_eq="Hellinger distance (equiprobabilized)",
        spearman="Spearman's rank correlation value",
        spearman_pval="P-value of the rank correlation"
    ),
    flags=list(
        flags="Warnings shown when unfusing or rearranging genotype names"
    ),
    properties_of_source_genotype=list(
        sourceGenotype="Name of the source genotype (if set to 'any', the statistical parameters correspond to the averages across all source genotypes, weighted by their frequency in the POM)",
        sourceGenotype_nMut="Number of mutations of the source genotype",
        sourceGenotype_freqInPOM="Frequency of appearance of the source genotype in the POM",
        sourceGenotype_accessible="Fraction of the replicates in which the source genotype was accessible according to the method",
        sampledFreq = "Frequency of genotype in the 20000 samples corresponding to the actual detection regime used (so over a simulation*detection, all genotypes add to 20000).",
        sampledProp = "Proportion of genotype in the 20000 samples corresponding to the actual detection regime used (so over a simulation*detection, all genotypes add to 1).",
        trueFreq = "Frecuency of genotype as the most common genotype during the simulations, computed over all the regularly spaced full population samples of each simulation (i.e., from the pops.by.time object); scaled to give equal weight to all 20000 simulations, these add up to 20000 for every simulation. These numbers are the same for all detection regimes.",
        trueProp = "Like trueFreq, but proportion; these add up to 1 for every simulation.",
        fitnessRank = "The rank of genotypes' fitness (where 1 is largest fitness).",
        fitnessRankNoZero = "The rank of genotypes' fitness (where 1 is largest fitness), but NA for non-viable genotypes (fitness <= 1e-9).",
        observedFreq = "Observed frequency of the genotype in the specific sample from which the CPM was built.",
        observedProp = "Observed proportion of the genotype in the specific sample from which the CPM was built (so sums to 1 in the replicate).",
        freqLocalMax = "Frequency with which the genotype is a local maximum (end of LOD)",
        propLocalMax = "Proportion of times the genotype is a local maximum (end of LOD); sums to 1 over the landscape"
    ),
    sampling=list(
        detect="Detection regime: when tumors are sampled (large, small, uniform)",
        sampledGenotypesDiversity="Diversity of sampled genotypes",
        sampledGenotypesNumber="Number of unique genotypes in the sample",
        sampledGenotypesGenesAbove_0.01="Number of unique genotypes with a frequency above 0.01 present in the sample",
        sampledGenotypesGenesAbove_0.1="Number of unique genotypes with a frequency above 0.1 present in the sample",
        numMuts_mean="Average number of mutations in the sample",
        numMuts_median="Median number of mutations in the sample",
        numMuts_var="Variance of the number of mutations in the sample",
        numMuts_stDev="Standard deviation of the number of mutations in the sample",
        numMuts_kurtosis="Kurtosis of the distribution of the number of mutations in the sample",
        numMuts_skewness="Skewness of the distribution of the number of mutations in the sample",
        pom_h="POM diversity",
        lod_h="LOD diversity",
        freq_most_freq_mean_no450="Frequency of the most frequent genotype in the sample",
        how_many_gt_5p_mean_no450="How many genotypes have a frequency > 5% in the sample"
    ),
    simulation_and_evolutionary_process=list(
        numAccessibleGenotypes="Number of accesible genotypes in the fitness landscape",
        initSize="Initial number of wild-type cells in the simulation",
        mutationRate="Mutation rate regime"
    ),
    fitness_landscape=list(
        rnst="Fitness landscape identifier (redundant check: one-to-one between id and rnst)",
        numGenes="Number of driver genes (7 or 10)",
        typeLandscape="Type of fitness landscape",
        numLocalPeaks="Number of local peaks (maxima) in the fitness landscape under the no-back mutation assumption",
        numObservedPeaks="Number of local peaks in the landscape that are actually visited in the evolutionary simulations",
        diversityObservedPeaks = "Diversity of observed local peaks"
    ),
    epistasis=list(
        epistMagn="Fraction of pairs of loci with magnitude epistasis in the landscape",
        epistSign="Fraction of pairs of loci with sign epistasis",
        epistRSign="Fraction of pairs of loci with reciprocal sign epistasis",
        w1="Fourier expansion of the landscape: fraction of coefficients of order 1",
        w2="Fourier expansion of the landscape: fraction of coefficients of order 2",
        w3="Fourier expansion of the landscape: fraction of coefficients of order 3 or higher",
        gamma="Correlation in fitness effects between genotypes that differ by one locus (Ferretti et al., 2016)"
    )
)


## fix formatting
cat("\n")
cat("Applying format")
cat("\n")

df_r$method[df_r$method == "caprese"] <- "CAPRESE"
df_r$method[df_r$method == "capri_aic"] <- "CAPRI_AIC"
df_r$method[df_r$method == "capri_bic"] <- "CAPRI_BIC"
df_r$method[df_r$method == "cbn_ot"] <- "CBN"
df_r$method[df_r$method == "mccbn"] <- "MCCBN"
df_r$method[df_r$method == "mhn"] <- "MHN"
df_r$method[df_r$method == "ot"] <- "OT"
df_r$method[df_r$method == "td-cbn_ot"] <- "CBN_td"
df_r$method[df_r$method == "td-mccbn"] <- "MCCBN_td"
df_r$method[df_r$method == "td-mhn"] <- "MHN_td"
df_r$method[df_r$method == "uw-cbn_ot"] <- "CBN_uw"
df_r$method[df_r$method == "uw-mccbn"] <- "MCCBN_uw"
df_r$method[df_r$method == "uw-ot"] <- "OT_uw"

df_r$typeLandscape[df_r$typeLandscape == "Local"] <- "Local maxima"


## A few additional checks. Takes a few minutes
stopifnot(min(df_r$sourceGenotype_freqInPOM, na.rm = TRUE) > 0)

u <- table(df_r[, c("id", "detect")])
uu <- as.data.frame(u)

v <- table(df_r[, c("id", "size_split")])
vv <- as.data.frame(v)

stopifnot(summary(as.vector(table(uu$id))) == 3) #### all a 3
stopifnot(table(uu$detect) == 1260) #### all 1260

stopifnot(summary(as.vector(table(vv$id))) == 3) #### all a 3
stopifnot(table(vv$size_split) == 1260) #### all 1260

## null now has detect
u2 <- table(df_r[, c("id", "detect")], useNA = "ifany")
stopifnot(identical(u2, u))

v2 <- table(df_r[, c("id", "size_split")], useNA = "ifany")
vv2 <- as.data.frame(v2)
stopifnot(length(unique(vv2$size_split)) == 4)
stopifnot(summary(as.vector(table(vv2$id))) == 4) ## the NA for null
stopifnot(table(vv2$size_split) == 1260) #### all 1260
stopifnot(sum(is.na(vv2$size_split)) == 1260)

x <- table(df_r[, c("detect", "size_split")])
xx <- as.data.frame(x)
stopifnot(length(unique(x)) == 1) ## all the same
stopifnot( (sum(x) + sum(df_r$method == "null")) == nrow(df_r) )

x2 <- table(df_r[, c("detect", "size_split")], useNA = "ifany")
xx2 <- as.data.frame(x2)
stopifnot(length(unique(x2)) == 2) ## all the same + the NA
stopifnot(length(unique(x2[x2 > 0])) == 2) ## all the same + the NA + 0 counts
stopifnot(sum(x2) == nrow(df_r))

w <- table(df_r[, c("id", "size_split", "detect")])
ww <- as.data.frame(w)

stopifnot(summary(as.vector(table(ww$id))) == 9) #### all a 9
stopifnot(table(ww$size_split) == (1260 * 3)) #### all 1260 * 3
stopifnot(table(ww$detect) == (1260 * 3)) #### all 1260 * 3

w2 <- table(df_r[, c("id", "size_split", "detect")], useNA = "ifany")
ww2 <- as.data.frame(w2)
stopifnot(summary(as.vector(table(ww2$id))) == 12) #### 9 + 3
## 3: NA for size * three detect for the null
ww3 <- ww2[ww2$Freq > 0, ]
stopifnot(identical(ww3, ww2))
stopifnot(table(ww2$size_split, useNA = "ifany") == (1260 * 3)) #### all 1260 * 3
stopifnot(table(ww2$detect, useNA = "ifany") == (1260 * 4)) #### all 1260 * (3 + 1)
## where the 4 are the three size split + the NA size split

stopifnot(table(ww2$size_split, useNA = "ifany") == c(3780, 3780, 3780, 3780)) 
stopifnot(table(ww2$detect, useNA = "ifany") == c(5040, 5040, 5040))

which_no_repl <- which(is.na(df_r$replicate))

t1 <- table(df_r[which_no_repl, "method"]) #### all should be null
stopifnot(names(t1) == "null")
stopifnot(length(unique(table(df_r$replicate))) == 1) #### identical
summary(df_r$replicate) #### 1 to 5, and NAs as given above
stopifnot(length(unique(table(df_r$replicate, useNA = "ifany"))) == 2)
stopifnot(summary(na.omit(df_r$replicate)) == c(1, 2, 3, 3, 4, 5))

tt <- table(df_r[, c("typeLandscape")])
stopifnot(length(tt) == 3) #### the three landscapes

ss <- table(df_r[, c("id", "typeLandscape")])
sss <- as.data.frame(ss)

stopifnot(summary(as.vector(table(sss$id))) == 3)
stopifnot(summary(as.vector(table(sss[sss$Freq > 0, ]$id))) == 1) #### all a 1
stopifnot(table(sss$typeLandscape) == 1260) #### all 1260
stopifnot(table(sss[sss$Freq > 0, ]$typeLandscape) == 420) #### all 420
## remember: 1260 landscapes, 1/3 in each type


                                        ## save output
cat("\n")
cat("Saving output")
cat("\n")


#### outFile_with_replicates <- file.path(saveDirectory,"table-replicates.rds")
outFile_with_replicates <- "table-replicates.rds"

system.time(
    saveRDS(list(data = df_r,
                 columnsExplained = columnsExplained),
            file = outFile_with_replicates,
            compress = FALSE)
)

## for the sake of curiosity
gc()






#### df_r_no_any <- dplyr::filter(df_r, sourceGenotype != "any")
#### df_r_no_any <- dplyr::filter(df_r_no_any, cpm != "null")

#### nrow(df_r_no_any)
date()
