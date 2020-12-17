## Compute the statistics averaged over genotypes
## weighting by POM frequency

## Beware: need to deal with NAs. A lot of boring stuff related to it
## below.

date()

rm(list = ls())

library(data.table)
library(dplyr)

## slow!
system.time(df <- readRDS("table-replicates.rds"))

## Simpler for now
data <- df$data

## We take averages, etc, and it is much simpler to average over the relevant
## columns of the relevant variables, but the constant ID-level or sample-level
## variables take that info from fl_sampl.
## So we will be merging info. Load here the file and
## rename variables. Leave ID and Detection as those
## are not output later
flDirectory <- "/home/jdc/next-genotype/fitland"
load(file.path(flDirectory, "fl_sampl.RData"))
fl_sampl$Detection[fl_sampl$Detection=="unif"] <- "uniform"
## Redundant
stopifnot(with(fl_sampl, all(Sampling == Detection)))
fl_sampl <- fl_sampl[, !(names(fl_sampl) %in% c("Sampling"))]
fl_sampl <- dplyr::rename(fl_sampl,
                      "numMuts_mean" =     "Mean_muts",
                      "numMuts_median" =   "Median_muts",
                      "numMuts_var" =      "Var_muts",
                      "numMuts_stDev" =    "Stdev_muts",
                      "numMuts_kurtosis" = "Kurtosis_muts",
                      "numMuts_skewness" =   "Skewness_muts",
                      "initSize" =         "Init_Size",
                      "mutationRate" =     "Mutation",
                      "numGenes" =         "Ngenes",
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
                      "typeLandscape" =    "type_Landscape"
                      )





#### Explanation of each column
##   Common to a lot of code below
columnsExplained <- list(
    note = list(note = "Not all columns need to be present, depending on what was averaged over"),
    identification=list(
        id="Data identifier",
        replicate = "Replicate number"
    ),
    cpm_and_input=list(
        method="Method used for next genotype prediction",
        size_split="Sample size: size of the input matrix given to the method"
    ),
    statistics=list(
        js_w_real = "Weighted average JS distance, genotype weights proportional to true frequency during simulations",
        js_w_sampl = "Weighted average JS distance, genotype weights proportional to frequency during simulations on sampling regime used",
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
        ## Ngenes = "Number of driver genes (7 or 10)", ## redundant
        typeLandscape="Type of fitness landscape",
        ## type_Landscape="Type of fitness landscape",
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

## This helps to verify the column names

allv <- sort(unname(unlist(lapply(columnsExplained, function(x) names(x)))))


###############################################
#####
#####   Prechecks of NAs and nrows
#####
###############################################

data_any <- dplyr::filter(data, sourceGenotype == "any")
data_any_null <- dplyr::filter(data_any, method == "null")
data_any_no_null <- dplyr::filter(data_any, method != "null")
data_no_any <- dplyr::filter(data, sourceGenotype != "any")
nrow(data_no_any)


##### Check number of replicates with valid data
nrow(data_any)
stopifnot(nrow(data_any) == (1260 * 9 * 5 * (length(unique(data$method)) - 1)+ 1260 * 3))

## There should be no NAs in the null model, ever
summary(data_any_null$js)
stopifnot(!any(is.na(data_any_null$js)))

## There should be NAs in the non-null. Those with CPM crashes
summary(data_any_no_null$js) ## some NAs
stopifnot(any(is.na(data_any_no_null$js)))

## Verify those are all NAs 
which_crash_any <- grep("ERROR_CPM_ANALYSIS", data_any$flags)
stopifnot(all(is.na(data_any[which_crash_any, "js"])))

## But there are NAs in some of those without CPM crashes
stopifnot(any(is.na(data_any[-which_crash_any, "js"])))
which(is.na(data_any[-which_crash_any, "js"]))
## and 6 cases total
sum(is.na(data_any[-which_crash_any, "js"]))
## Check one example
dplyr::filter(data, (id == "Dv8B7RHWDguRzH1W") &
                    (size_split == 50) &
                    (detect == "large") &
                    (replicate == 1) &
                    (method == "OT")
              )[, c(1, 4, 5, 6, 7, 8, 10, 14, 51, 52, 55, 56, 20)]
## js returned a NaN

## Do we have many NaNs?
sum(is.nan(data_no_any$js)) ## 6, as above

stopifnot(sum(is.nan(data_no_any$js)) == 6)
## over all, 12
sum(is.nan(data$js)) ## 12, as we add the 6 for genotypes to those for the corresponding array
stopifnot(sum(is.nan(data$js)) == 12)

## See them all. Of course, any is irrelevant
dplyr::filter(data, (is.nan(js)) &
                    (sourceGenotype != "any")
              )[, c(1, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 51, 52, 55, 56, 20)]
## the 12 is shown to compare
## and check
## Those NaNs are to be turned to 0s. See all gory details in why-the-NaN.R
## basically, a numerical issue in borderline cases (hey, 6 out of > 39 million)
## with KL divergence.


## Do!
data$js[is.nan(data$js)] <- 0

## Recreate and redo checks
data_any <- dplyr::filter(data, sourceGenotype == "any")
data_any_null <- dplyr::filter(data_any, method == "null")
data_any_no_null <- dplyr::filter(data_any, method != "null")
data_no_any <- dplyr::filter(data, sourceGenotype != "any")
nrow(data_no_any)

##### Check number of replicates with valid data
nrow(data_any)
stopifnot(nrow(data_any) == (1260 * 9 * 5 * (length(unique(data$method)) - 1)+ 1260 * 3))

## There should be no NAs in the null model, ever
summary(data_any_null$js)
stopifnot(!any(is.na(data_any_null$js)))

## There should be NAs in the non-null. Those with CPM crashes
summary(data_any_no_null$js) ## some NAs
stopifnot(any(is.na(data_any_no_null$js)))

## Verify those are all NAs 
which_crash_any <- grep("ERROR_CPM_ANALYSIS", data_any$flags)
stopifnot(all(is.na(data_any[which_crash_any, "js"])))

## Now, no NAs when no CPM crashes
stopifnot(!any(is.na(data_any[-which_crash_any, "js"])))
which(is.na(data_any[-which_crash_any, "js"]))




##################################################
######
######     Array-level averages of statistics
######
##################################################

unique_combs <- c("id", "replicate", "method", "size_split", "detect")


#### Verify no need to scale, as the props (true and sampled) sum to 1.

## tmp files to check scaling
tmp_true <-
    dplyr::group_by_at(data_no_any[, c("trueProp", unique_combs)],
                       vars(one_of(unique_combs)))
tmp_sampled <-
    dplyr::group_by_at(data_no_any[, c("sampledProp", unique_combs)],
                       vars(one_of(unique_combs)))

ds_true_nonarm <- dplyr::summarize(tmp_true,
                            sum = sum(.data[["trueProp"]],
                                      na.rm = FALSE))
ds_true <- dplyr::summarize(tmp_true,
                            sum = sum(.data[["trueProp"]],
                                      na.rm = TRUE))

ds_sampled_nonarm <- dplyr::summarize(tmp_sampled,
                            sum = sum(.data[["sampledProp"]],
                                      na.rm = FALSE))
ds_sampled <- dplyr::summarize(tmp_sampled,
                            sum = sum(.data[["sampledProp"]],
                                      na.rm = TRUE))

## There are NAs here
summary(ds_true_nonarm[, "sum"])
## No NAs
summary(ds_true[, "sum"])
## Both add up to one
stopifnot(isTRUE(all.equal(ds_sampled$sum, rep(1, nrow(ds_sampled)))))
stopifnot(isTRUE(all.equal(ds_true$sum, rep(1, nrow(ds_true)))))

### </Verify no need to scale


####### Matrix-level
###      Logic
###           Subset and keep only 
###           Multiply by weights the statistic, then sum
###           Then, merge the data with fl_statistics

## Yes, keep sourceGenotype initially, before merging, so we can
## check if needed
matsum <- data_no_any[, c(unique_combs, "sourceGenotype",
                          "trueProp", "sampledProp", "js")]

matsum$js_real <- with(matsum, js * trueProp)
matsum$js_samp <- with(matsum, js * sampledProp)

matsum2 <- dplyr::group_by_at(matsum,
                              vars(one_of(unique_combs)))

## None of these are good ideas. See below when we explain mysum
## matsum3 <- dplyr::summarize(matsum2,
##                             js_w_real = sum(.data[["js_real"]],
##                                             na.rm = TRUE),
##                             js_w_sampl = sum(.data[["js_samp"]],
##                                              na.rm = TRUE)
##                            )

## matsum4 <- dplyr::summarize(matsum2,
##                             js_w_real = sum(.data[["js_real"]],
##                                             na.rm = FALSE),
##                             js_w_sampl = sum(.data[["js_samp"]],
##                                              na.rm = FALSE)
##                             )

## See below for the reason of mysum
mysum <- function(x) if (all(is.na(x))) NA else sum(x, na.rm=TRUE)
matsum5 <- dplyr::summarize(matsum2,
                            js_w_real = mysum(.data[["js_real"]]),
                            js_w_sampl = mysum(.data[["js_samp"]])
                            )

matsum5b <- dplyr::summarize(matsum2,
                            js_w_real = mysum(js_real),
                            js_w_sampl = mysum(js_samp)
                            )
stopifnot(identical(matsum5, matsum5b))
rm(matsum5b)
####  What about NAs?
##    Recall: NAs can arise from NAs in weights or NAs in js

### NAs in js
## Remember many js as NA...
which_na_js <- which(is.na(data_no_any$js))
## ... but all when ERROR_CPM_ANALYSIS
error_cpm <- grep("ERROR_CPM_ANALYSIS", data_no_any$flag, fixed = TRUE)
stopifnot(!any(is.na(data_no_any[-error_cpm, "js"])))

## An example
data_no_any[error_cpm[1:3], c(1, 4, 5, 6, 7, 8, 10, 14, 51, 52, 55, 56, 20)]


### NAs in weights:
nas_weight_true <- which(is.na(data_no_any$trueProp))
## Yes, a bunch
head(data_no_any[nas_weight_true, c(1, 4, 5, 6, 7, 8, 10, 14, 51, 52, 55, 56)])
nrow(data_no_any[nas_weight_true, c(1, 4, 5, 6, 7, 8, 10, 14, 51, 52, 55, 56)])

## For example, genotype D,I: it appears once in POM but not in the true genotypes
dplyr::filter(data_no_any, (id == "10aJ10EAqKH7VA3ccp") &
                     (method == "CAPRESE") &
                     (size_split == 50) &
                     (detect == "large") &
                     (replicate == 1)
              )[, c(1, 4, 5, 6, 7, 8, 10, 14, 51, 52, 55, 56)]

## This was documented in merge-additional-info.R

tmp_true_Prop_na <-
    data_no_any[nas_weight_true, c(1, 4, 5, 6, 7, 8, 10, 14, 51, 52, 55, 56)]
tmp_id_genot_prop_na <-
    paste0(tmp_true_Prop_na$id,"_", tmp_true_Prop_na$sourceGenotype)

length(unique(tmp_id_genot_prop_na)) ## 97 different instances
## out of the
length(
    unique(paste0(data_no_any[data_no_any$method != "null", "id"],
                  "_",
                  data_no_any[data_no_any$method != "null", "sourceGenotype"])))
## 66864
## With
max(tmp_true_Prop_na$sourceGenotype_freqInPOM) ## 0.00015
## so 3 cases out of the 20000 simulations as max
## Of course, this happens much more with sampled. Sure.


## We want to keep all level combinations, even if NA for all (e.g., those
## crashed) and for those with missing values in some genotypes (e.g., one
## of the props) the sum with those removed.
## 

## Example of crash
## Definitely want the NAs here, and we want to keep the level combination
dplyr::filter(matsum5, (id == "10H8qiQAOoVrnY9T3") &
                     (method == "CAPRESE") &
                     (size_split == 50) &
                     (detect == "small") &
                     (replicate == 2)
              )

## An example where a genotype did not have a value for trueProp
## We want an average over the non-NA rows
dplyr::filter(matsum5, (id == "10aJ10EAqKH7VA3ccp") &
                     (method == "CAPRESE") &
                     (size_split == 50) &
                     (detect == "large") &
                     (replicate == 1)
              )

## Check. The two statistics 0 <= x <= 1
stopifnot(na.omit(matsum5$js_w_real) <= 1)
stopifnot(na.omit(matsum5$js_w_real) >= 0)
stopifnot(na.omit(matsum5$js_w_sampl) <= 1)
stopifnot(na.omit(matsum5$js_w_sampl) >= 0)

## Dimensions
stopifnot(
    nrow(matsum5) == (1260 * 9 * 5 * (length(unique(data_no_any$method)) - 1)+ 1260 * 3)
)

stopifnot(nrow(matsum5[matsum5$method == "null",]) == 1260 * 3)


## Average over replicates
unique_combs_over_repl <- c("id", "method", "size_split", "detect")
## same as setdiff(unique_combs, "replicate")

matsum6 <- dplyr::group_by_at(matsum5,
                              vars(one_of(unique_combs_over_repl)))

array_statistics_over_replicate <-
    dplyr::summarize(matsum6,
                     js_w_real = mean(.data[["js_w_real"]], na.rm = TRUE),
                     js_w_sampl = mean(.data[["js_w_sampl"]], na.rm = TRUE)
                     )

## Check specific cases and dimensions
stopifnot(nrow(array_statistics_over_replicate) ==
          ((1260 * 9 * (length(unique(matsum6$method)) - 1) ) + (1260 * 3)))
stopifnot(
    nrow(
        array_statistics_over_replicate[
            array_statistics_over_replicate$method == "null", ]) ==
          (1260 * 3))



## this crashed in one replicate
dplyr::filter(array_statistics_over_replicate, (id == "10H8qiQAOoVrnY9T3") &
                                              (method == "CAPRESE") &
                                              (size_split == 50) &
                                              (detect == "small"))

## trueProp missing in one case
dplyr::filter(array_statistics_over_replicate, (id == "10aJ10EAqKH7VA3ccp") &
                                              (method == "CAPRESE") &
                                              (size_split == 50) &
                                              (detect == "large"))



array_statistics <- dplyr::left_join(matsum5, fl_sampl,
                                     by = c("id" = "ID", "detect" = "Detection"))

array_statistics_over_replicate <-
    dplyr::left_join(array_statistics_over_replicate, fl_sampl,
                     by = c("id" = "ID", "detect" = "Detection"))


stopifnot(
    nrow(array_statistics) ==
    (1260 * 9 * 5 * (length(unique(data_no_any$method)) - 1)+ 1260 * 3)
)

stopifnot(nrow(array_statistics[array_statistics$method == "null",]) == 1260 * 3)

stopifnot(nrow(array_statistics_over_replicate) ==
          ((1260 * 9 * (length(unique(matsum6$method)) - 1) ) + (1260 * 3)))

stopifnot(
    nrow(
        array_statistics_over_replicate[
            array_statistics_over_replicate$method == "null", ]) ==
          (1260 * 3))




system.time(saveRDS(list(data = array_statistics,
             columnsExplained = columnsExplained),
        file = "array_statistics.rds",
        compress = FALSE
        ))

system.time(saveRDS(list(data = array_statistics_over_replicate,
             columnsExplained = columnsExplained),
        file = "array_statistics_over_replicate.rds",
        compress = FALSE
        ))


### Checking COLUMN NAMES array column names

setdiff(allv, sort(colnames(array_statistics)))
setdiff(allv, sort(colnames(array_statistics_over_replicate)))

setdiff(sort(colnames(array_statistics)), allv)
setdiff(sort(colnames(array_statistics_over_replicate)), allv)



##################################################
######
######     Full data for models, without null
######
##################################################

## Just for the the hell of it. Well, to really drive that
## data.table is much better for me.
system.time(
    data_no_any_tibble <- as_tibble(data_no_any)
) ## 0.004
## note: not using setDT

system.time(
    data_no_any_dt <- as.data.table(data_no_any)
) ## 8.97


system.time(
    data_no_any_no_null_tibble <- dplyr::filter(data_no_any_tibble,
                                                method != "null")
) ## 12.4

system.time(
    data_no_any_no_null_dt <- data_no_any_dt[method != "null", ]
) ## 5.4


system.time(saveRDS(list(data = data_no_any_no_null_dt,
             columnsExplained = columnsExplained),
        file = "data_no_any_no_null.rds",
        compress = FALSE
        ))

## Merge over replicates
data_no_any_no_null_merged <-
    data_no_any_no_null_dt[,
                           .(js = mean(js, na.rm = TRUE),
                             sampledFreq = mean(sampledFreq, na.rm = TRUE),
                             sampledFreq_unique = length(unique(sampledFreq)),
                             sampledProp = mean(sampledProp, na.rm = TRUE),
                             sampledProp_unique = length(unique(sampledProp)),
                             trueFreq = mean(trueFreq, na.rm = TRUE),
                             trueFreq_unique = length(unique(trueFreq)),
                             trueProp = mean(trueProp, na.rm = TRUE),
                             trueProp_unique = length(unique(trueProp)),
                             observedFreq = mean(observedFreq, na.rm = TRUE),
                             observedFreq_unique = length(unique(observedFreq)),
                             observedProp = mean(observedProp, na.rm = TRUE),
                             observedProp_unique = length(unique(observedProp)),
                             fitnessRank = mean(fitnessRank, na.rm = TRUE),
                             fitnessRank_unique = length(unique(fitnessRank)),
                             fitnessRankNoZero = mean(fitnessRankNoZero, na.rm = TRUE),
                             fitnessRankNoZero_unique = length(unique(fitnessRankNoZero)),
                             freqLocalMax = mean(freqLocalMax, na.rm = TRUE),
                             freqLocalMax_unique = length(unique(freqLocalMax)),
                             propLocalMax = mean(propLocalMax, na.rm = TRUE),
                             propLocalMax_unique = length(unique(propLocalMax))
                           ## , flags = paste(unique(flags))
                           ## , flags_unique = length(unique(flags))
                           , sourceGenotype_accessible =
                                 mean(sourceGenotype_accessible, na.rm = TRUE)
                           , sourceGenotype_accessible_unique =
                                 length(unique(sourceGenotype_accessible))
                             ),
                           by = list(id, sourceGenotype, method, size_split, detect,
                                     ## this is not a unique comb. but we want it in
                                     ## output
                                     sourceGenotype_nMut)]

nrow(data_no_any_no_null_merged)

## checks
stopifnot(
    all(data_no_any_no_null_merged$sampledFreq_unique ==1) 
)
stopifnot(
    all(data_no_any_no_null_merged$sampledProp_unique ==1)
)
stopifnot(
    all(data_no_any_no_null_merged$trueFreq_unique ==1)
)
stopifnot(
    all(data_no_any_no_null_merged$trueProp_unique ==1)
)
stopifnot(
    !all(data_no_any_no_null_merged$observedFreq_unique ==1)
)
stopifnot(
    !all(data_no_any_no_null_merged$observedProp_unique ==1)
)
stopifnot(
    all(data_no_any_no_null_merged$fitnessRank_unique ==1)
)
stopifnot(
    all(data_no_any_no_null_merged$fitnessRankNoZero_unique ==1)
)
stopifnot(
    all(data_no_any_no_null_merged$freqLocalMax_unique ==1)
)
stopifnot(
    all(data_no_any_no_null_merged$propLocalMax_unique ==1)
)




## Those columns are not fully useless. Any column with a non-1 value
## helps remember those are not constant over replicates




fl_sampl_dt <- as.data.table(fl_sampl)

data_no_any_no_null_over_replicate <- merge(
    data_no_any_no_null_merged,
    fl_sampl,
    all.x = TRUE,
    by.x = c("id", "detect"),
    by.y = c("ID", "Detection")
)


nrow(data_no_any_no_null_merged)

nrow(data_no_any_no_null_dt)/nrow(data_no_any_no_null_merged)


system.time(saveRDS(list(data = data_no_any_no_null_over_replicate,
             columnsExplained = columnsExplained),
        file = "data_no_any_no_null_over_replicate.rds",
        compress = FALSE
        ))



### Checking COLUMN NAMES full data column names

setdiff(allv, sort(colnames(data_no_any_no_null_dt)))
setdiff(allv, sort(colnames(data_no_any_no_null_over_replicate)))

setdiff(sort(colnames(data_no_any_no_null_dt)), allv)
setdiff(sort(colnames(data_no_any_no_null_over_replicate)), allv)


rm(data_no_any_dt)
rm(data_no_any_no_null_tibble)
rm(data_no_any_tibble)
rm(data_no_any_no_null_dt)
rm(data_no_any_no_null_merged)
rm(data_no_any_no_null_over_replicate)

gc()



##################################################
######
######     Full data for models
######
##################################################

## Just for the the hell of it. Well, to really drive that
## data.table is much better for me.

## Note: I want to preserve the null
## Sure, there are no replicates of the null
##   so the merged is as good as the non-merged if we just want
##   to use the null
## Yes, many steps from above repeated

system.time(
    data_no_any_tibble <- as_tibble(data_no_any)
) ## 0.004

## Note: not using setDT
system.time(
    data_no_any_dt <- as.data.table(data_no_any)
) ## 8.97



system.time(saveRDS(list(data = data_no_any_dt,
             columnsExplained = columnsExplained),
        file = "data_no_any.rds",
        compress = FALSE
        ))

## Merge over replicates
data_no_any_merged <-
    data_no_any_dt[,
                   .(js = mean(js, na.rm = TRUE),
                     sampledFreq = mean(sampledFreq, na.rm = TRUE),
                     sampledFreq_unique = length(unique(sampledFreq)),
                     sampledProp = mean(sampledProp, na.rm = TRUE),
                     sampledProp_unique = length(unique(sampledProp)),
                     trueFreq = mean(trueFreq, na.rm = TRUE),
                     trueFreq_unique = length(unique(trueFreq)),
                     trueProp = mean(trueProp, na.rm = TRUE),
                     trueProp_unique = length(unique(trueProp)),
                     observedFreq = mean(observedFreq, na.rm = TRUE),
                     observedFreq_unique = length(unique(observedFreq)),
                     observedProp = mean(observedProp, na.rm = TRUE),
                     observedProp_unique = length(unique(observedProp)),
                     fitnessRank = mean(fitnessRank, na.rm = TRUE),
                     fitnessRank_unique = length(unique(fitnessRank)),
                     fitnessRankNoZero = mean(fitnessRankNoZero, na.rm = TRUE),
                     fitnessRankNoZero_unique = length(unique(fitnessRankNoZero)),
                     freqLocalMax = mean(freqLocalMax, na.rm = TRUE),
                     freqLocalMax_unique = length(unique(freqLocalMax)),
                     propLocalMax = mean(propLocalMax, na.rm = TRUE),
                     propLocalMax_unique = length(unique(propLocalMax))
                   ## , flags = paste(unique(flags))
                   ## , flags_unique = length(unique(flags))
                   , sourceGenotype_accessible =
                         mean(sourceGenotype_accessible, na.rm = TRUE)
                   , sourceGenotype_accessible_unique =
                         length(unique(sourceGenotype_accessible))                     
                     ),
                   by = list(id, sourceGenotype, method, size_split, detect,
                             sourceGenotype_nMut)]

nrow(data_no_any_merged)

## checks
stopifnot(
    all(data_no_any_merged$sampledFreq_unique ==1)
)
stopifnot(
    all(data_no_any_merged$sampledProp_unique ==1)
)
stopifnot(
    all(data_no_any_merged$trueFreq_unique ==1)
)
stopifnot(
    all(data_no_any_merged$trueProp_unique ==1)
)
stopifnot(
    !all(data_no_any_merged$observedFreq_unique ==1)
)
stopifnot(
    !all(data_no_any_merged$observedProp_unique ==1)
)
stopifnot(
    all(data_no_any_merged$fitnessRank_unique ==1)
)
stopifnot(
    all(data_no_any_merged$fitnessRankNoZero_unique ==1)
)
stopifnot(
    all(data_no_any_merged$freqLocalMax_unique ==1)
)
stopifnot(
    all(data_no_any_merged$propLocalMax_unique ==1)
)

## Those columns are not fully useless. Any column with a non-1 value
## helps remember those are not constant over replicates


fl_sampl_dt <- as.data.table(fl_sampl)

data_no_any_over_replicate <- merge(
    data_no_any_merged,
    fl_sampl,
    all.x = TRUE,
    by.x = c("id", "detect"),
    by.y = c("ID", "Detection")
)


nrow(data_no_any_merged)
nrow(data_no_any_over_replicate)

nrow(data_no_any_dt)/nrow(data_no_any_over_replicate)


system.time(saveRDS(list(data = data_no_any_over_replicate,
             columnsExplained = columnsExplained),
        file = "data_no_any_over_replicate.rds",
        compress = FALSE
        ))

### Checking COLUMN NAMES full data column names

setdiff(allv, sort(colnames(data_no_any_dt)))
setdiff(allv, sort(colnames(data_no_any_over_replicate)))

setdiff(sort(colnames(data_no_any_dt)), allv)
setdiff(sort(colnames(data_no_any_over_replicate)), allv)



rm(data_no_any_dt)
rm(data_no_any_tibble)
rm(data_no_any_merged)

gc()






##################################################
######
######     Number of mutations averages of statistics
######
##################################################


unique_combs_muts <- c("id", "replicate", "sourceGenotype_nMut",
                       "method", "size_split", "detect")

unique_combs_muts_over_repl <- c("id", "sourceGenotype_nMut",
                                 "method", "size_split", "detect")

mutsum <- data_no_any[, c(unique_combs_muts,
                          "sourceGenotype",
                          "trueProp", "sampledProp", "js")]


## Do now use weighted.mean. See help: Missing values in ‘w’ are
## nothandled specially and so give a missing value as the result.
## See, e.g., https://stackoverflow.com/questions/40269022/weighted-average-in-r-using-na-weights
## See examples below

## my_weighted_mean <- function(x, w){
##     df_omit <- na.omit(data.frame(x, w))
##     weighted.mean(df_omit$x, df_omit$w)
## }

my_weighted_mean2 <- function(x, w){
    df_omit <- na.omit(cbind(x, w))
    weighted.mean(df_omit[, 1], df_omit[, 2])
}


## Use data table

dtmutsum <- as.data.table( data_no_any[, c(unique_combs_muts,
                          "sourceGenotype",
                          "trueProp", "sampledProp", "js")])
system.time(
    setkeyv(dtmutsum, cols = unique_combs_muts)
) ## 1.7

system.time(
    dtmutsum4 <- dtmutsum[ ,
                          .(js_w_real = my_weighted_mean2(js, trueProp),
                            js_w_sampl = my_weighted_mean2(js, sampledProp)),
                          by = unique_combs_muts]
) ## 174


stopifnot(nrow(dtmutsum4) <=
          (1260 * 9 * 5 * (length(unique(data$method)) - 1) * 0.5 * (8 + 11)
           + 1260 * 3 * 0.5 * (8 + 11)))

## Some rows are expected to be missing, as not all simulations have the
## genotype with all mutations. For example:

maxminmut <- dtmutsum[ ,
                   .(max_muts = max(sourceGenotype_nMut),
                    min_muts = min(sourceGenotype_nMut)),
                   by = "id"]

stopifnot(all(maxminmut$min_muts == 0))
stopifnot(all(maxminmut$max_muts <= 10))
## This is not just 7 or 10.
table(maxminmut$max_muts)


maxmut_full <- dtmutsum[ ,
                        .(max_muts = max(sourceGenotype_nMut)),
                        by = list(id, replicate, method, size_split, detect)]

## Why? Because not all intermediate levels present
nrow(dtmutsum4) < sum(1 + maxmut_full$max_muts)

## Check null for ease
Null_dtmutsum4 <- dtmutsum4[ method == "null",  ]
Null_maxmut_full <- maxmut_full[ method == "null", ]
Null_dtmutsum <- dtmutsum[ method == "null", ]

stopifnot(nrow(Null_dtmutsum4) < sum(1 + Null_maxmut_full$max_muts))

table(Null_maxmut_full$max_muts)

## No monotonic decrease. Missing intermediate genotypes
table(Null_dtmutsum4$sourceGenotype_nMut)

missing_muts <- function(x) {
    maxm <- max(x)
    missing <- setdiff(1:maxm, unique(x))
    if(length(missing)) {
        return(paste(missing, collapse = ", "))
    } else {
        return("None")
    }
}

Null_which_missing <- Null_dtmutsum[ ,
    .(which_missing_mut = missing_muts(sourceGenotype_nMut)),
    by = list(id, replicate)
]

## Yes, missing intermediate genotypes
table(Null_which_missing$which_missing_mut)


## For instance
Null_which_missing[which(Null_which_missing$which_missing_mut == "2, 6"), ]

dtmutsum[ (id == "utFV4f4ZwaTtF6F7") & (method == "null") ]
sum(dtmutsum[ (id == "utFV4f4ZwaTtF6F7") & (method == "null") & (detect == "large")]$trueProp)


## We are missing these many: total, minus None and add two for the two
## cases that miss two
sum(table(Null_which_missing$which_missing_mut)) - 1206 + 2 ## 56

sum(1 + Null_maxmut_full$max_muts) - nrow(Null_dtmutsum4) ## = 3 * 56
## the 3 * is because we had the three detect in the null but no sizes
## yes, yes, in the Null_which_missing, we do not look over detect, either.

stopifnot(sum(1 + Null_maxmut_full$max_muts) - nrow(Null_dtmutsum4) == (3 * 56))


## With the non null
No_null_dtmutsum4 <- dtmutsum4[ method != "null",  ]
No_null_maxmut_full <- maxmut_full[ method != "null", ]
No_null_dtmutsum <- dtmutsum[ method != "null", ]

stopifnot(nrow(No_null_dtmutsum4) < sum(1 + No_null_maxmut_full$max_muts))
sum(1 + No_null_maxmut_full$max_muts) - nrow(No_null_dtmutsum4) ## 32760

table(No_null_maxmut_full$max_muts)
table(No_null_dtmutsum4$sourceGenotype_nMut)

table(No_null_dtmutsum4$sourceGenotype_nMut)/
    table(Null_dtmutsum4$sourceGenotype_nMut)

## a constant of 195 = 3 * 5 * (length(unique(data$method)) - 1)
## where the 3 instead of 9 is because we have 3 detect and no sizes

stopifnot(all(
    (table(No_null_dtmutsum4$sourceGenotype_nMut)/
     table(Null_dtmutsum4$sourceGenotype_nMut)) == 195))

## Repeat missing combinations. Just a multipe of the null.
No_null_which_missing <- No_null_dtmutsum[ ,
    .(which_missing_mut = missing_muts(sourceGenotype_nMut)),
    by = list(id, replicate, method, size_split, detect)
]

table(No_null_which_missing$which_missing_mut)

table(No_null_which_missing$which_missing_mut)/table(Null_which_missing$which_missing_mut)
## 585

## these are the combinations for each ID
9 * 5 * (length(unique(data$method)) - 1) ## 585

stopifnot(
    all(
    (table(No_null_which_missing$which_missing_mut)/
     table(Null_which_missing$which_missing_mut)) == 585)
)

## Missing these many total
sum(table(No_null_which_missing$which_missing_mut)) - 705510 + 2 * 585
## 32760

sum(1 + No_null_maxmut_full$max_muts) - nrow(No_null_dtmutsum4)
## 32760

stopifnot(
    (sum(table(No_null_which_missing$which_missing_mut)) - 705510 + 2 * 585) ==
    (sum(1 + No_null_maxmut_full$max_muts) - nrow(No_null_dtmutsum4))
)

sum(1 + maxmut_full$max_muts) - nrow(dtmutsum4) ## 32928
32760 + 3 * 56 ## 32928

stopifnot(
(sum(1 + maxmut_full$max_muts) - nrow(dtmutsum4) ) ==
32760 + 3 * 56
)


## Average over replicates
system.time(
    dtmutsum5 <- dtmutsum4[ ,
                          .(js_w_real = mean(js_w_real, na.rm = TRUE),
                            js_w_sampl = mean(js_w_sampl, na.rm = TRUE)),
                          by = unique_combs_muts_over_repl]
) 


## Merge the data for sample.
## There must surely be a data table equivalent but


num_mut_statistics <-
    dplyr::left_join(as.data.frame(dtmutsum4),
                     fl_sampl,
                     by = c("id" = "ID", "detect" = "Detection"))

stopifnot(nrow(dtmutsum4) == nrow(num_mut_statistics))

num_mut_statistics_over_replicate <-
    dplyr::left_join(as.data.frame(dtmutsum5),
                     fl_sampl,
                     by = c("id" = "ID", "detect" = "Detection"))

stopifnot(nrow(dtmutsum5) == nrow(num_mut_statistics_over_replicate))

print(object.size(num_mut_statistics_over_replicate),
      units = "MB")



system.time(saveRDS(list(data = num_mut_statistics,
             columnsExplained = columnsExplained),
        file = "num_mut_statistics.rds",
        compress = FALSE
        ))

system.time(saveRDS(list(data = num_mut_statistics_over_replicate,
             columnsExplained = columnsExplained),
        file = "num_mut_statistics_over_replicate.rds",
        compress = FALSE
        ))




### Checking COLUMN NAMES num mut data column names

setdiff(allv, sort(colnames(num_mut_statistics)))
setdiff(allv, sort(colnames(num_mut_statistics_over_replicate)))

setdiff(sort(colnames(num_mut_statistics)), allv)
setdiff(sort(colnames(num_mut_statistics_over_replicate)), allv)


##################################################
######
######     What are the NAs?
######
##################################################

##########  Summary

#### Array stats
##      All NAs are from cases where the CPM crashed


#### Array stats, merged over replicate
####    Two cases, all failed in the CPM. They are NaNs


#### Num mut stats (for the number of mutations statistics)
##      The NaNs in JS with real weight are due to having all the JS be NA
##      in those combinations
##      
##      The NaNs in JS with sample weights are all due to either all the
##      sample weights being NAs or all the JS being NAs


#### Num mut stats (for the number of mutations statistics) over replicate
##      If there are NaNs in the original data (Num mut stats) and they
##      affect all replicates, necessarily we have NaNs.
##
##      All NAs are NaNs, but they do not affect at both types of sampling
##      identically (143 for JS with real weights, 61517 for JS with
##      sample weights)


##########  </Summary


##########   Details

#### Array stats, over replicate
## All NAs are NaNs, and affect at both sampling identically
na_js_r <- which(is.na(array_statistics_over_replicate$js_w_real))
na_js_s <- which(is.na(array_statistics_over_replicate$js_w_sampl))

nan_js_r <- which(is.nan(array_statistics_over_replicate$js_w_real))
nan_js_s <- which(is.nan(array_statistics_over_replicate$js_w_sampl))

stopifnot(length(setdiff(na_js_r, nan_js_r)) == 0)
stopifnot(length(setdiff(na_js_s, nan_js_s)) == 0)
stopifnot(length(setdiff(nan_js_s, nan_js_r)) == 0)

length(nan_js_r)
length(nan_js_s)

## Two simulations
as.data.frame(array_statistics_over_replicate[nan_js_r, ])[, c(1:6)]

tmp1 <- data[data$id %in% c("CsQCk5oitWE3VWiq", "gHFY6JsOn10PnBYf3"), ]

tmp2 <- dplyr::filter(tmp1, (method != "any") &
                            (detect == "small") &
                            (size_split == 50)
                      )
## They were all NAs
stopifnot(all(is.na(tmp2$js)))
## ... because all of these failed in CPM
stopifnot(all(tmp2$flags == "not a matrix: ERROR_CPM_ANALYSIS"))
## and the average of NAs is a NaN if we na.rm = TRUE
mean(c(NA, NA), na.rm = TRUE)

#### Array stats
##   All NAs are from cases where the CPM crashed
na_js_r <- which(is.na(array_statistics$js_w_real))
na_js_s <- which(is.na(array_statistics$js_w_sampl))

## No NaNs
nan_js_r <- which(is.nan(array_statistics$js_w_real))
nan_js_s <- which(is.nan(array_statistics$js_w_sampl))
stopifnot(length(nan_js_r) == 0)
stopifnot(length(nan_js_s) == 0)

## Missing in same places for both statistics
stopifnot(length(setdiff(na_js_r, na_js_s)) == 0)

array_statistics[na_js_r, ]

## Looks like those that crash, as seen from this
which_crash_any <- grep("ERROR_CPM_ANALYSIS", data_any$flags)
stopifnot(length(which_crash_any) == length(na_js_r))

## More carefully
## (see also analysis in lines 243 and ff. But here all of the entries
## for an array result in NA)

data_no_any_crash <- dplyr::filter(data_no_any,
                                   flags == "not a matrix: ERROR_CPM_ANALYSIS")

## Verify all crashed and obtain conditions of crash
data_no_any_crash <- as.data.table(data_no_any_crash)

dnac <-
    data_no_any_crash[ ,
                      .(count_js = length(js),
                        count_js_na = sum(is.na(js)),
                        count_error_cpm = sum(flags == "not a matrix: ERROR_CPM_ANALYSIS")),
                      by = list(id, replicate, method, size_split, detect)]

stopifnot(dnac$count_js == dnac$count_js_na)
stopifnot(dnac$count_js == dnac$count_error_cpm)
## Makes sense
stopifnot(nrow(dnac) == length(na_js_r))
## Strong check
crashed_o <- sort(with(dnac, paste(id, replicate, method, size_split, detect)))
array_crashed <- array_statistics[na_js_r, ]
crashed_a <- sort(with(array_crashed, paste(id, replicate, method, size_split, detect)))
stopifnot(identical(crashed_o, crashed_a))




#### Num mut stats
## All NAs are NaNs, but they do not affect at both types of sampling identically

na_js_r <- which(is.na(num_mut_statistics$js_w_real))
na_js_s <- which(is.na(num_mut_statistics$js_w_sampl))

nan_js_r <- which(is.nan(num_mut_statistics$js_w_real))
nan_js_s <- which(is.nan(num_mut_statistics$js_w_sampl))

stopifnot(length(setdiff(na_js_r, nan_js_r)) == 0)
stopifnot(length(setdiff(na_js_s, nan_js_s)) == 0)
stopifnot(identical(na_js_r, nan_js_r))
stopifnot(identical(na_js_s, nan_js_s))

stopifnot(length(setdiff(nan_js_s, nan_js_r)) != 0)

## Either all JS were NA or NaN or all weights were NA or NaN
my_weighted_mean2(c(1, 2), c(NA, 1))
my_weighted_mean2(c(NA, 2), c(1, 1))
my_weighted_mean2(c(NA, NA), c(2, 1))
my_weighted_mean2(c(1, 2), c(NA, NA))

my_weighted_mean2(c(1, 2), c(NaN, 1))
my_weighted_mean2(c(NaN, 2), c(1, 1))
my_weighted_mean2(c(NaN, NaN), c(2, 1))
my_weighted_mean2(c(1, 2), c(NaN, NaN))


dtmut_nas <- dtmutsum[ ,
                      .(
                          all_js_na = all(is.na(js)),
                          all_js_nan = all(is.nan(js)),
                          all_s_na = all(is.na(sampledProp)),
                          all_s_nan = all(is.nan(sampledProp)),
                          all_r_na = all(is.na(trueProp)),
                          all_r_nan = all(is.nan(trueProp)),
                          sum_js_na = sum(is.na(js)),
                          sum_js_nan = sum(is.nan(js)),
                          sum_s_na = sum(is.na(sampledProp)),
                          sum_s_nan = sum(is.nan(sampledProp)),
                          sum_r_na = sum(is.na(trueProp)),
                          sum_r_nan = sum(is.nan(trueProp))
                          ),
                      by = unique_combs_muts]

## Yes, there are NaNs in JS, 6 cases. See file "why-the-NaN.R"

## We are only interested in cases where any of the components of the
## calculation are always NA or NaN

dtmut_nas_allnan <- dtmut_nas[all_js_na | all_js_nan |
                              all_s_na | all_s_nan |
                              all_r_na | all_r_nan,
                              ]

## js is never NaN in this set
stopifnot(all(dtmut_nas_allnan$sum_js_nan == 0))

## "real" are never all NA or NaNs in this set
stopifnot(all(dtmut_nas_allnan$all_r_na == FALSE))
stopifnot(all(dtmut_nas_allnan$all_r_nan == FALSE))
## sampling are never NaN
stopifnot(all(dtmut_nas_allnan$all_s_nan == FALSE))


summary(dtmut_nas_allnan)

## So cause of NAs in the array are all the JS were NA (NA, not NaN) or
## all the "sample" weights were all NA (again, NA, not NaN)

## JS, with real weights
js_na <- dtmut_nas_allnan[all_js_na == TRUE, ]

nan_js_o <- num_mut_statistics[nan_js_r, ]

id_nan_js_o <- sort(with(nan_js_o, paste(id, replicate, sourceGenotype_nMut,
                                       method, size_split, detect)))
id_na_js_x <- sort(with(js_na, paste(id, replicate, sourceGenotype_nMut,
                                       method, size_split, detect)))

## The NaNs in JS with real weight are due to having all the JS be NA
## in those combinations
stopifnot(identical(id_nan_js_o, id_na_js_x))

## JS with "sample" weights?
ws_na <- dtmut_nas_allnan[all_s_na == TRUE, ]
na_js_s_o <- num_mut_statistics[nan_js_s, ]

id_na_js_s_o <- sort(with(na_js_s_o, paste(id, replicate, sourceGenotype_nMut,
                                       method, size_split, detect)))
id_ws_na <- sort(with(ws_na, paste(id, replicate, sourceGenotype_nMut,
                                       method, size_split, detect)))
## Nope!
stopifnot(!identical(id_na_js_s_o, id_ws_na))
## Why? Because we are missing elements! If all JS is NA you also have
## JS with sampled weight be NaN. This shows we are missing things.
stopifnot(nrow(ws_na) < nrow(na_js_s_o))


## Now, add that condition
ws_na_or_js_na <- dtmut_nas_allnan[(all_s_na == TRUE) | (all_js_na == TRUE), ]
stopifnot(nrow(ws_na) < nrow(ws_na_or_js_na))
## There 8415 cases where we had weights but JS was NA
nrow(ws_na_or_js_na) - nrow(ws_na) 
##
with(dtmut_nas_allnan, sum( (all_js_na == TRUE) & (all_s_na == FALSE)))

id_ws_na_or_js_na <- sort(with(ws_na_or_js_na, paste(id, replicate, sourceGenotype_nMut,
                                       method, size_split, detect)))



stopifnot(identical(id_na_js_s_o, id_ws_na_or_js_na))
##  So NaNs in JS with sample weights are all due to either all the sample
##  weights being NAs or all the JS being NAs


#### Num mut stats, over replicate
##    All NAs are NaNs, but they do not affect at both types of sampling
##    identically

na_js_r <- which(is.na(num_mut_statistics_over_replicate$js_w_real))
na_js_s <- which(is.na(num_mut_statistics_over_replicate$js_w_sampl))

nan_js_r <- which(is.nan(num_mut_statistics_over_replicate$js_w_real))
nan_js_s <- which(is.nan(num_mut_statistics_over_replicate$js_w_sampl))

stopifnot(length(setdiff(na_js_r, nan_js_r)) == 0)
stopifnot(length(setdiff(na_js_s, nan_js_s)) == 0)
stopifnot(identical(na_js_r, nan_js_r))
stopifnot(identical(na_js_s, nan_js_s))

stopifnot(length(setdiff(nan_js_s, nan_js_r)) != 0)

length(nan_js_r)
length(nan_js_s)



##################################################
######
######     Appendix: Why we want to use mysum
######
##################################################

## Suppose this
matsum3 <- dplyr::summarize(matsum2,
                            js_w_real = sum(.data[["js_real"]],
                                            na.rm = TRUE),
                            js_w_sampl = sum(.data[["js_samp"]],
                                             na.rm = TRUE)
                           )

matsum3b <- dplyr::summarize(matsum2,
                            js_w_real = sum(js_real,
                                            na.rm = TRUE),
                            js_w_sampl = sum(js_samp,
                                             na.rm = TRUE)
                           )
stopifnot(identical(matsum3, matsum3b))

matsum4 <- dplyr::summarize(matsum2,
                            js_w_real = sum(.data[["js_real"]],
                                            na.rm = FALSE),
                            js_w_sampl = sum(.data[["js_samp"]],
                                             na.rm = FALSE)
                            )

dplyr::filter(matsum3, (id == "10H8qiQAOoVrnY9T3") &
                     (method == "CAPRESE") &
                     (size_split == 50) &
                     (detect == "small") &
                     (replicate == 2)
              )

dplyr::filter(matsum4, (id == "10H8qiQAOoVrnY9T3") &
                     (method == "CAPRESE") &
                     (size_split == 50) &
                     (detect == "small") &
                     (replicate == 2)
              )
## Definitely want the NAs here, not a 0. So na.rm = FALSE or equivalent
## when creating those (matsum4)
## But why does na.rm = TRUE gives a 0?
dplyr::filter(matsum2, (id == "10H8qiQAOoVrnY9T3") &
                     (method == "CAPRESE") &
                     (size_split == 50) &
                     (detect == "small") &
                     (replicate == 2)
              )

## What gives here?
## See this thread. A little exploration
df1 <- data.frame(v1 = c("a", "a", "a", "b", "b"),
                  v2 = c("X", "X", "Y", "Y", "Y"),
                  v3 = c(NA, NA, 3, 4, NA))

gbdf1 <- dplyr::group_by(df1, v1, v2)

## The first case makes no sense?
dplyr::summarize(gbdf1, v3_sum = sum(v3, na.rm = TRUE))
dplyr::summarize(gbdf1, v3_sum = sum(v3, na.rm = FALSE))

## Contrast with
aggregate(v3 ~ v1 + v2, function(x) sum(x), data = df1)
aggregate(v3 ~ v1 + v2, function(x) sum(x), data = df1, na.action = na.pass)
aggregate(v3 ~ v1 + v2, function(x) sum(x, na.rm = TRUE),
          data = df1, na.action = na.pass)
## the mean seems not to be affected
aggregate(v3 ~ v1 + v2, function(x) mean(x), data = df1)
aggregate(v3 ~ v1 + v2, function(x) mean(x), data = df1, na.action = na.pass)
aggregate(v3 ~ v1 + v2, function(x) mean(x, na.rm = TRUE),
          data = df1, na.action = na.pass)

## but look at tapply
with(df1, tapply(v3, list(v2, v1), sum, na.rm = TRUE))

aggregate(v3 ~ v1 + v2, function(x) mean(x, na.rm = TRUE), data = df1)
aggregate(v3 ~ v1 + v2, function(x) mean(x, na.rm = FALSE), data = df1)

## Use Dalgaards suggestion
mysum <- function(x) if (all(is.na(x))) NA else sum(x,na.rm=T)
## sanity is back
dplyr::summarize(gbdf1, v3_sum = mysum(v3))
## I am missing a combination
aggregate(v3 ~ v1 + v2, mysum, data = df1)
## sanity back
aggregate(v3 ~ v1 + v2, mysum, data = df1, na.action = na.pass)
## extra level
with(df1, tapply(v3, list(v2, v1), mysum))
## NULL for extra level
with(df1, tapply(v3, list(v2, v1), mysum, simplify = FALSE))

## Using anything else would lead to differences here

## An example where CPM crashed
## Definitely want the NAs here, not a 0. So na.rm = FALSE or equivalent
dplyr::filter(matsum3, (id == "10H8qiQAOoVrnY9T3") &
                     (method == "CAPRESE") &
                     (size_split == 50) &
                     (detect == "small") &
                     (replicate == 2)
              )

dplyr::filter(matsum4, (id == "10H8qiQAOoVrnY9T3") &
                     (method == "CAPRESE") &
                     (size_split == 50) &
                     (detect == "small") &
                     (replicate == 2)
              )
## An example where a genotype did not have a value for trueProp
## We want an average over the non-NA rows
dplyr::filter(matsum3, (id == "10aJ10EAqKH7VA3ccp") &
                     (method == "CAPRESE") &
                     (size_split == 50) &
                     (detect == "large") &
                     (replicate == 1)
              )

## An example where a genotype did not have a value for trueProp
## We want an average over the non-NA rows
dplyr::filter(matsum4, (id == "10aJ10EAqKH7VA3ccp") &
                     (method == "CAPRESE") &
                     (size_split == 50) &
                     (detect == "large") &
                     (replicate == 1)
              )


##################################################
######
######     Appendix: Why we want to use my_weighted_mean
######
##################################################


df2 <- data.frame(v1 = c("a", "a", "a", "b", "b", "c", "c", "d", "d", "e", "e"),
                  v2 = c("X", "X", "Y", "Y", "Y", "U", "U", "V", "V", "W", "W"),
                  v3 = c(NA, NA, 3, 4, NA, 1, 2, 3,  4,  5, 6),
                  v4 = c(1,  2,  3, 4, 5, NA, 7, NA, NA, 4, 5))

gbdf2 <- dplyr::group_by(df2, v1, v2)

## Do now use weighted.mean. See help: Missing values in ‘w’ are
## nothandled specially and so give a missing value as the result.
## See, e.g., https://stackoverflow.com/questions/40269022/weighted-average-in-r-using-na-weights


df3 <-
     as.data.frame( dplyr::summarize(gbdf2, v3_wg_s = weighted.mean(v3, w
     = v4, na.rm = TRUE)))
## see value for combination c, U


df4 <-
    as.data.frame(
        dplyr::summarize(gbdf2,
                         v3_wg_s = my_weighted_mean2(v3, w = v4)))

stopifnot(identical(df4[, 3],
                    c(
                        NaN,
                        3.0,
                        4.0,
                        2.0,
                        NaN,
                        (20 + 30)/(9)
                    )))

## data.table

dtdf2 <- as.data.table(df2)

dtdf4 <- dtdf2[, .(v3_wg_s = my_weighted_mean2(v3, v4),
                   v3_m2 = mean(v3)), ## for the hell of it
               by = list(v1, v2)]

gb <- c("v1", "v2")

dtdf5 <- dtdf2[, .(v3_wg_s = my_weighted_mean2(v3, v4),
                   v3_m2 = mean(v3)), ## for the hell of it
               by = gb]

stopifnot(identical(dtdf4, dtdf5))

stopifnot(identical(as.data.frame(dtdf4[, 3])[, 1],
                    c(
                        NaN,
                        3.0,
                        4.0,
                        2.0,
                        NaN,
                        (20 + 30)/(9)
                    )))



## ## Using data.table now
## mutsum2 <- dplyr::group_by_at(mutsum,
##                               vars(one_of(unique_combs_muts)))


## ## Not sure this does what I intend, as no "by"
## system.time({
##     require(dtplyr)
##     dtmutsum2 <- dtplyr::lazy_dt(mutsum2)
##     dtmutsum3 <- as_tibble(dplyr::summarize(dtmutsum2,
##                             js_w_real = my_weighted_mean2(js, trueProp),
##                             js_w_sampl = my_weighted_mean2(js, sampledProp)))
## })



## ## Because dplyr sucks so much but so smart
## options(dplyr.width = Inf)

## ## Wow, slow
## system.time(mutsum3 <- dplyr::summarize(mutsum2,
##                             js_w_real = my_weighted_mean2(js, trueProp),
##                             js_w_sampl = my_weighted_mean2(js, sampledProp))
##             )

date()
