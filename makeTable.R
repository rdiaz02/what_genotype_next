# ----------

# This script takes transition matrices generated with the scripts:
#   next-genotype_transitionMatrix-sim.R
#   next-genotype_transitionMatrix-cpm.R
# and structured into formatted objects by the script:
#   structData.R
# And performs matrix similarity tests. It then generates a data frame
# containing all the data, and saves it to the specified location.

# ----------

source("oncoFunctions.R")

# read arguments from command line
args <- commandArgs(trailingOnly=T)
possibleArgs <- c("-scratch","-mode=b","-mode=pr")
if (!all(args %in% possibleArgs)) {
  
  cat("WARNING: invalid arguments")
  cat("\n")
  cat("Possible arguments are:")
  cat("\n")
  cat("     ")
  cat(paste(possibleArgs,collapse=", "))
  cat("\n")
  stop()
  
} else {
  
  ### list here arguments and variables they toggle:
  
  # format-only mode
  onlyFormat <- T
  if ("-scratch" %in% args) {
    onlyFormat <- F
  } else {
    onlyFormat <- T
  }
  if (onlyFormat) {
    cat("[makeTable] Running formatting-only routine")
    cat("\n")
  } else {
    cat("[makeTable] Making table from scratch")
    cat("\n")
  }
  
  # accessibility mismatches
  modePR <- F
  if ("-mode=pr" %in% args & !("-mode=b" %in% args)) {
    modePR <- T
  } else if (!("-mode=pr" %in% args) & "-mode=b" %in% args) {
    modePR <- F
  } else if ("-mode=pr" %in% args & "-mode=b" %in% args) {
    stop("Only one option for -mode can be passed")
  }
  if (modePR) {
    cat("[makeTable] Using 'punishment/reward' mode to handle accessibility mismatches")
    cat("\n")
  } else {
        cat("[makeTable] Using 'benevolent' mode to handle accessibility mismatches")
    cat("\n")
  }
  
}
cat("\n")

# directories to load/save from
if (!onlyFormat) loadDirectory <- askDir(defaultDir="./data",
                                         message="Enter directory where data are stored.")
flDirectory <- askDir(defaultDir="./fitland",
                      message="Enter directory where fitness landscape & sampling stats are stored.")
saveDirectory <- askDir(defaultDir="./",
                        message="Enter directory to save output to.")

# create save directory (if it doesn't exist)
dir.create(saveDirectory,showWarnings=F)
  
if (!onlyFormat) {
  
  # list all files
  cat("Checking provided directory")
  cat("\n")
  cat(paste("  > ",loadDirectory,sep=""))
  cat("\n")
  files <- list.files(loadDirectory,full.names=T,rec=F)
  cat(paste("  > Found",length(files),"files"))
  cat("\n")
  cat("\n")
  
  # loop through files
  cat("Processing files")
  cat("\n")
  pboptions(type="txt")
  df <- pblapply(files,
                function(file) {
  
                  # load file
                  load(file)
                  
                  # initialize output data frame
                  df <- data.frame(id=character(),
                                   nGenes=numeric(),
                                   typeLandscape=character(),
                                   cpm=character(),
                                   size_split=character(),
                                   detect=character(),
                                   replicate=numeric(),
                                   sourceGenotype=character(),
                                   sourceGenotype_nMut=numeric(),
                                   sourceGenotype_freqInPOM=numeric(),
                                   sourceGenotype_accessible=logical(),
                                   sqDiff=numeric(),
                                   sqDiff_fix=numeric(),
                                   sqDiff_eq=numeric(),
                                   sqDiff_eq_fix=numeric(),
                                   js=numeric(),
                                   js_fix=numeric(),
                                   js_eq=numeric(),
                                   js_eq_fix=numeric(),
                                   hellinger=numeric(),
                                   hellinger_fix=numeric(),
                                   hellinger_eq=numeric(),
                                   hellinger_eq_fix=numeric(),
                                   spearman=numeric(),
                                   spearman_pval=numeric(),
                                   flags=character())
                  
                  # compare matrices (for all CPMs, size_splits, detection regimes and replicates)
                  for (j in 1:length(data[["cpm"]])) { # CPM index
                    for (k in 1:length(data[["cpm"]][[1]])) { # size_split index  
                      for (l in 1:length(data[["cpm"]][[1]][[1]])) { # detect index
                        for (m in 1:length(data[["cpm"]][[1]][[1]][[1]])) { # replicate index
                          print(c(j,k,l,m))
                          t <-
                            compareMatrices(data[["cpm"]][[j]][[k]][[l]][[m]][["transitionMatrix"]],
                                            data[["sim"]][["transitionMatrix"]],
                                            rowWeights=data[["sim"]][["timesInPOM"]],
                                            threshold=c(0,0))
                          t <- data.frame(id=data[["ID"]],
                                          nGenes=data[["nGenes"]],
                                          typeLandscape=data[["typeLandscape"]],
                                          cpm=names(data[["cpm"]])[j],
                                          size_split=names(data[["cpm"]][[j]])[k],
                                          detect=names(data[["cpm"]][[j]][[k]])[l],
                                          replicate=m,
                                          sourceGenotype=rownames(t),
                                          sourceGenotype_nMut=nMut(rownames(t)),
                                          sourceGenotype_freqInPOM=c(
                                            data[["sim"]][["timesInPOM"]]/max(
                                              data[["sim"]][["timesInPOM"]]
                                            ),NA),
                                          sourceGenotype_accessible=c(
                                            isAccessible(
                                              data[["cpm"]][[j]][[k]][[l]][[m]][["transitionMatrix"]]
                                            ),NA),
                                          sqDiff=t$sqDiff,
                                          sqDiff_fix=t$sqDiff_fix,
                                          sqDiff_eq=t$sqDiff_eq,
                                          sqDiff_eq_fix=t$sqDiff_eq_fix,
                                          js=t$js,
                                          js_fix=t$js_fix,
                                          js_eq=t$js_eq,
                                          js_eq_fix=t$js_eq_fix,
                                          hellinger=t$hellinger,
                                          hellinger_fix=t$hellinger_fix,
                                          hellinger_eq=t$hellinger_eq,
                                          hellinger_eq_fix=t$hellinger_eq_fix,
                                          spearman=t$spearman,
                                          spearman_pval=t$spearman_pval,
                                          flags=data[["cpm"]][[j]][[k]][[l]][[m]][["flags"]])
                          df <- rbind(df,t)
                        }
                      }
                    }
                  }
                  
                  # compare matrices with null model
                  t <- compareMatrices(data[["null"]][["transitionMatrix"]],
                                       data[["sim"]][["transitionMatrix"]],
                                       rowWeights=data[["sim"]][["timesInPOM"]],
                                       threshold=c(0,0))
                  t <- data.frame(id=data[["ID"]],
                                  nGenes=data[["nGenes"]],
                                  typeLandscape=data[["typeLandscape"]],
                                  cpm="null",
                                  size_split=NA,
                                  detect=NA,
                                  replicate=NA,
                                  sourceGenotype=rownames(t),
                                  sourceGenotype_nMut=nMut(rownames(t)),
                                  sourceGenotype_freqInPOM=c(
                                    data[["sim"]][["timesInPOM"]]/max(
                                      data[["sim"]][["timesInPOM"]]
                                    ),NA),
                                  sourceGenotype_accessible=c(
                                    isAccessible(
                                      data[["null"]][["transitionMatrix"]]
                                    ),NA),
                                  sqDiff=t$sqDiff,
                                  sqDiff_fix=t$sqDiff_fix,
                                  sqDiff_eq=t$sqDiff_eq,
                                  sqDiff_eq_fix=t$sqDiff_eq_fix,
                                  js=t$js,
                                  js_fix=t$js_fix,
                                  js_eq=t$js_eq,
                                  js_eq_fix=t$js_eq_fix,
                                  hellinger=t$hellinger,
                                  hellinger_fix=t$hellinger_fix,
                                  hellinger_eq=t$hellinger_eq,
                                  hellinger_eq_fix=t$hellinger_eq_fix,
                                  spearman=t$spearman,
                                  spearman_pval=t$spearman_pval,
                                  flags="")
                  
                  df <- rbind(df,t)
                  
                  # fiter out genots. unaccessible by the simulations & the methods
                  ### FIXME: do NOT filter this beforehand so we can later choose
                  ### between the "penalize/reward" or the "benevolent" strategies
                  if(F) {
                    n <- which(df$sourceGenotype_freqInPOM==0 &
                                 df$sourceGenotype_accessible==F)
                    if (length(n)) df <- df[-n,]
                  }
                  
                  # sanity check: if something is wrong with equiprobabilizing 
                  # and/or fixing, stop and warn (see sanityStats() in
                  # oncoFunctions source file)
                  sanityStats(df)
                  
                  return(df)
  
                },
                cl=detectCores())
  
  # save temporary file
  cat("\n")
  cat("Saving temporary file (unformatted table)")
  cat("\n")
  outFile <- file.path(saveDirectory,"makeTable_tmp.rds")
  saveRDS(df,file=outFile)

} else { # if we only want to format a previously generated table, load it
  
  cat("Attempting to load unformatted table from current directory")
  cat("\n")
  df <- readRDS("makeTable_tmp.rds")
  
}


## To rerun from here:
##    - load makeTable_tmp.rds
##    - modePR to whatever (FALSE most time, unless punishment/reward)
##  - then run this block and the local block

  


# choose an approach and adjust statistics accordingly
cat("\n")
cat("Adjusting statistics")
cat("\n")
pboptions(type="txt")
df <- pblapply(df,
               function(df) {
                 
                 stats <- c("sqDiff","js","hellinger")
                 
                 if (modePR) {
                   # if 'punishment/reward': use stats as is (remove _fix ones)
                   df[,paste(stats,"_fix",sep="")] <- NULL
                   df[,paste(stats,"_eq_fix",sep="")] <- NULL
                 } else {
                   # if 'benevolent': first use the _fix versions of stats
                   df[,stats] <-
                     df[,paste(stats,"_fix",sep="")]
                   df[,paste(stats,"_eq",sep="")] <-
                     df[,paste(stats,"_eq_fix",sep="")]
                   df[,paste(stats,"_fix",sep="")] <- NULL
                   df[,paste(stats,"_eq_fix",sep="")] <- NULL
                   
                   # then remove table elements of unaccesible genotypes
                   # (in the simuls)
                   df <- df[df$sourceGenotype_freqInPOM>0 | df$sourceGenotype=="any",]
                 }
                 
                 return(df)
                 
               },
               cl=detectCores())

save(file = "pre-table-replicates.RData", df, compress = FALSE)

### No longer used. We stop here and prepare large table and combine
### additional info in merge-additional-info.R


## ## Since loading and running previous step takes a while


## ## RDU: code for creating, testing, saving replicate
## ## file. All inside local. No stuff should exist afterwords
## ## to affect the rest

## local({
##     ## Watch out: up to here, using about 25% of RAM in Draco (~ 90 GB)
##     ## I don't quite get it, since the object size is 6 GB
##     cat("Replicates version: Combining data into single table\n")

##     df_r <- dplyr::bind_rows(df)

##     cat("Replicates version: formatting and columns")

##                                         # change "" genotype into "WT" or "root" and set accessibility to TRUE
##     df_r$sourceGenotype[df_r$sourceGenotype==""] <- "root"

##                                         # format "size_split" column
##     df_r$size_split[df_r$size_split=="size_split_50"] <- "50"
##     df_r$size_split[df_r$size_split=="size_split_200"] <- "200"
##     df_r$size_split[df_r$size_split=="size_split_4000"] <- "4000"
##     df_r$size_split <- as.numeric(df_r$size_split)

##                                         # format "detect" column
##     df_r$detect[df_r$detect=="detect_large"] <- "large"
##     df_r$detect[df_r$detect=="detect_small"] <- "small"
##     df_r$detect[df_r$detect=="detect_unif"] <- "unif"
    
    
##                                         # load fitness landscape and sampling characteristics, set all "unif" to "uniform"
##     load("./fitland/fl_sampl.RData")

##     df_r$detect <- as.character(df_r$detect)
##     df_r$detect[df_r$detect=="unif"] <- "uniform"
##     fl_sampl$Detection <- as.character(fl_sampl$Detection)
##     fl_sampl$Detection[fl_sampl$Detection=="unif"] <- "uniform"

##                                         # merge tables
##     cat("\n")
##     cat("Attaching fitness landscape & sampling statistics")
##     cat("\n")


##     ## To compare merge and left_join
##     compare_merge_left_join <- FALSE

##     if(compare_merge_left_join) {
##         system.time(
##             df_r_merge <- merge(df_r, fl_sampl,
##                                 by.x=c("id","detect"),by.y=c("ID","Detection"),
##                                 all=T,suffixes=NULL)
##         ) ## 164 s
##     }

##     ## This is 10x to 16x faster
##     system.time(
##         df_r_join <- dplyr::left_join(df_r, fl_sampl,
##                                       by = c("id" = "ID", "detect" = "Detection"))
##     ) 

##     if(compare_merge_left_join) {
##         ## checks
##         stopifnot(identical(nrow(df_r_join), nrow(df_r_merge)))
##         stopifnot(identical(ncol(df_r_join), ncol(df_r_merge)))
##         ## order of colnames changes, but columns are the same
##         stopifnot(all(colnames(df_r_join) %in% colnames(df_r_merge)))
##         stopifnot(all(colnames(df_r_merge) %in% colnames(df_r_join)))

##         df_r_join <- df_r_join[, colnames(df_r_merge)]

##         ## still, not identical
##         which_no_i <-
##             sapply(1:ncol(df_r_merge),
##                    function(i) !identical(df_r_merge[, i], df_r_join[, i]))
##         ## what gives here?
##         ## ordering differs, starting with row 191
##         df_r_merge[190:194, c(1:4)]
##         df_r_join[190:194, c(1:4)]

##         ## Lets fix that to compare
##         ## (yes, this resorting is slow)
##         df_r_m_2 <- df_r_merge[order(df_r_merge$id,
##                                      df_r_merge$detect,
##                                      df_r_merge$size_split,
##                                      df_r_merge$cpm,
##                                      df_r_merge$replicate,
##                                      df_r_merge$sourceGenotype),
##                                ]

##         df_r_j_2 <- df_r_join[order(df_r_join$id,
##                                     df_r_join$detect,
##                                     df_r_join$size_split,
##                                     df_r_join$cpm,
##                                     df_r_join$replicate,
##                                     df_r_join$sourceGenotype),
##                               ]

##         ## This apparently fails
##         identical(df_r_j_2, df_r_m_2)

##         ## But this does not
##         which_no_i_2 <-
##             sapply(1:ncol(df_r_m_2),
##                    function(i) !identical(df_r_m_2[, i], df_r_j_2[, i]))

##         ## As usual, the culprit are attributes
##         identical(attributes(df_r_m_2)$row.names,
##                   attributes(df_r_j_2)$row.names)
##         ## And that happened in one of join/merge
##         ## (this should not be identical, unless one or both set reset
##         ## rownames after join/merge)
##         head(row.names(df_r_join)) 
##         head(row.names(df_r_merge))
##         ## Reset row names
##         row.names(df_r_m_2) <- NULL
##         row.names(df_r_j_2) <- NULL


##         ## And now, both are identical
##         stopifnot(identical(df_r_j_2, df_r_m_2))

##         ## game over for checks
##         rm(df_r_merge, df_r_m_2, df_r_j_2, which_no_i, which_no_i_2)
##         gc()
##     }

##     df_r <- df_r_join

##     ## Adding weights (computed externally with makeTableWeights.R)
##     ## recall to change unif to uniform
##     load("weightsAll.RData")
##     nrow(weightsAll)
##     weightsAll$detect[weightsAll$detect == "unif"] <- "uniform"

##     ## FIXME: wrong: I called it WT, should be called root
##     weightsAll$Genotype[weightsAll$Genotype == "WT"] <- "root"

    
##     (nr_df_r_before_join <- nrow(df_r))
##     df_r <- dplyr::left_join(df_r, weightsAll,
##                              by = c("id" = "ID",
##                                     "detect",
##                                     "sourceGenotype" = "Genotype"))
    
##     nrow(df_r) ## same as before

##     stopifnot(nrow(df_r) == nr_df_r_before_join)

    
##     ####### Anything in POM not in True? Yes, but ignore #######
##     ##  (how? in C++ code POM is updated even if no creation of
##     ##   output table because it happens in between keepEverys)
##     df_r_no_any <- dplyr::filter(df_r, sourceGenotype != "any")
##     df_r_no_any <- dplyr::filter(df_r_no_any, method != "null")

##     not_in_true <- which(is.na(df_r_no_any$trueProp))
##     ## There are some
##     length(not_in_true)
##     ## but very tiny frequencies : 3 simulations at most
##     summary(df_r_no_any$sourceGenotype_freqInPOM[not_in_true])
##     max(df_r_no_any$sourceGenotype_freqInPOM[not_in_true])
##     20000 * max(df_r_no_any$sourceGenotype_freqInPOM[not_in_true])

##     ## And in 95% of the cases in just one sample
##     20000 * quantile(
##                 df_r_no_any$sourceGenotype_freqInPOM[not_in_true],
##                 probs = c(0.95, 0.96, 0.97, 0.99))
    
##     ## Affects 82 IDs
##     length(unique(df_r_no_any$id[not_in_true]))


##     ## FIXME: check if anything in True not in POM
##     ##   a left join, with weightsAll on left, or a full join
##     ##   using minimal set of columns?



##     rm(df_r_join)
##     gc()

    
##     ## Can be handy when dealing with replicates
##     ## and, yes, arrange much faster than [order(...
##     ## About 25 seconds
##     system.time(df_r <- dplyr::arrange(df_r,
##                                        id,
##                                        detect,
##                                        size_split,
##                                        cpm,
##                                        sourceGenotype,
##                                        replicate
##                                        ))

##     df_r <-
##         df_r[, c(
##             "id",
##             "cpm", "size_split",
##             "sqDiff", "sqDiff_eq", "js", "js_eq", "hellinger", "hellinger_eq", "spearman", "spearman_pval",
##             "flags",
##             "sourceGenotype", "sourceGenotype_nMut", "sourceGenotype_freqInPOM","sourceGenotype_accessible",
##             "detect", "sampledGenotypesDiversity", "sampledGenotypesNumber",
##             "sampledGenotypesGenesAbove_0.01", "sampledGenotypesGenesAbove_0.1",
##             "Mean_muts", "Median_muts", "Var_muts",
##             "Stdev_muts", "Kurtosis_muts","Skewness_muts",
##             "pom_h", "lod_h",
##             "freq_most_freq_mean_no450", "how_many_gt_5p_mean_no450",
##             "num_accessible_genots", "Init_Size", "Mutation",
##             "rnst", "nGenes", "typeLandscape",
##             "num_local_peaks", "num_observed_peaks",
##             "epist_magn", "epist_sign", "epist_rsign",
##             "w.1.", "w.2.", "w.3..","gamma",
##             "replicate",
##             "SampledFreq", "SampledProp", "TrueFreq", "TrueProp"
##         )]

##     colnames(df_r) <-
##         c(
##             "id",
##             "method","size_split",
##             "sqDiff", "sqDiff_eq", "js", "js_eq", "hellinger", "hellinger_eq", "spearman", "spearman_pval",
##             "flags",
##             "sourceGenotype","sourceGenotype_numMuts","sourceGenotype_freqInPOM","sourceGenotype_accessible",
##             "detect", "sampledGenotypesDiversity", "sampledGenotypesNumber",
##             "sampledGenotypesGenesAbove_0.01", "sampledGenotypesGenesAbove_0.1",
##             "numMuts_mean", "numMuts_median", "numMuts_var",
##             "numMuts_stDev", "numMuts_kurtosis","nMuts_skewness",
##             "pom_h", "lod_h",
##             "freq_most_freq_mean_no450", "how_many_gt_5p_mean_no450",
##             "numAccessibleGenotypes", "initSize", "mutationRate",
##             "rnst", "numGenes", "typeLandscape",
##             "numLocalPeaks", "numObservedPeaks",
##             "epistMagn", "epistSign", "epistRSign",
##             "w1", "w2", "w3","gamma",
##             "replicate",
##             "sampledFreq", "sampledProp", "trueFreq", "trueProp"
##         )


##     ## explanation of each column
##     columnsExplained <- list(
##         identification=list(
##             id="Data identifier",
##             replicate = "Replicate number"
##         ),
##         cpm_and_input=list(
##             method="Method used for next genotype prediction",
##             size_split="Sample size: size of the input matrix given to the method"
##         ),
##         statistics=list(
##             sqDiff="Square root of average of squared differences",
##             sqDiff_eq="Square root of average of squared differences (equiprobabilized)",
##             js="Jensen-Shannon distance (square root of the Jensen-Shannon divergence), in log base 2 units",
##             js_eq="Jensen-Shannon distance (equiprobabilized)",
##             hellinger="Hellinger distance",
##             hellinger_eq="Hellinger distance (equiprobabilized)",
##             spearman="Spearman's rank correlation value",
##             spearman_pval="P-value of the rank correlation"
##         ),
##         flags=list(
##             flags="Warnings shown when unfusing or rearranging genotype names"
##         ),
##         properties_of_source_genotype=list(
##             sourceGenotype="Name of the source genotype (if set to 'any', the statistical parameters correspond to the averages across all source genotypes, weighted by their frequency in the POM)",
##             sourceGenotype_numMuts="Number of mutations of the source genotype",
##             sourceGenotype_freqInPOM="Frequency of appearance of the source genotype in the POM",
##             sourceGenotype_accessible="Fraction of the replicates in which the source genotype was accessible according to the method",
##             sampledFreq = "Frequency of genotype in the 20000 samples corresponding to the actual detection regime used (so over a simulation*detection, all genotypes add to 20000).",
##             sampledProp = "Proportion of genotype in the 20000 samples corresponding to the actual detection regime used (so over a simulation*detection, all genotypes add to 1).",
##             trueFreq = "Frecuency of genotype as the most common genotype during the simulations, computed over all the regularly spaced full population samples of each simulation (i.e., from the pops.by.time object); scaled to give equal weight to all 20000 simulations, these add up to 20000 for every simulation. These numbers are the same for all detection regimes.",
##             trueProp = "Like trueFreq, but proportion; these add up to 1 for every simulation."
##         ),
##         sampling=list(
##             detect="Detection regime: when tumors are sampled (large, small, uniform)",
##             sampledGenotypesDiversity="Diversity of sampled genotypes",
##             sampledGenotypesNumber="Number of unique genotypes in the sample",
##             sampledGenotypesAbove_0.01="Number of unique genotypes with a frequency above 0.01 present in the sample",
##             sampledGenotypesAbove_0.01="Number of unique genotypes with a frequency above 0.1 present in the sample",
##             numMuts_mean="Average number of mutations in the sample",
##             numMuts_median="Median number of mutations in the sample",
##             numMuts_var="Variance of the number of mutations in the sample",
##             numMuts_stDev="Standard deviation of the number of mutations in the sample",
##             numMuts_kurtosis="Kurtosis of the distribution of the number of mutations in the sample",
##             numMuts_skewness="Skewness of the distribution of the number of mutations in the sample",
##             pom_h="POM diversity",
##             lod_h="LOD diversity",
##             freq_most_freq_mean_no450="Frequency of the most frequent genotype in the sample",
##             how_many_gt_5p_mean_no450="How many genotypes have a frequency > 5% in the sample"
##         ),
##         simulation_and_evolutionary_process=list(
##             numAccessibleGenotypes="Number of accesible genotypes in the fitness landscape",
##             initSize="Initial number of wild-type cells in the simulation",
##             mutationRate="Mutation rate regime"
##         ),
##         fitness_landscape=list(
##             rnst="Fitness landscape identifier (redundant check: one-to-one between id and rnst)",
##             numGenes="Number of driver genes (7 or 10)",
##             typeLandscape="Type of fitness landscape",
##             numLocalPeaks="Number of local peaks (maxima) in the fitness landscape under the no-back mutation assumption",
##             numObservedPeaks="Number of local peaks in the landscape that are actually visited in the evolutionary simulations"
##         ),
##         epistasis=list(
##             epistMagn="Fraction of pairs of loci with magnitude epistasis in the landscape",
##             epistSign="Fraction of pairs of loci with sign epistasis",
##             epistRSign="Fraction of pairs of loci with reciprocal sign epistasis",
##             w1="Fourier expansion of the landscape: fraction of coefficients of order 1",
##             w2="Fourier expansion of the landscape: fraction of coefficients of order 2",
##             w3="Fourier expansion of the landscape: fraction of coefficients of order 3 or higher",
##             gamma="Correlation in fitness effects between genotypes that differ by one locus (Ferretti et al., 2016)"
##         )
##     )


##                                         # fix formatting
##     cat("\n")
##     cat("Applying format")
##     cat("\n")
##     df_r$method <- gsub("^caprese$", "CAPRESE", df_r$method)
##     df_r$method <- gsub("^capri_aic$", "CAPRI_AIC", df_r$method)
##     df_r$method <- gsub("^capri_bic$", "CAPRI_BIC", df_r$method)
##     df_r$method <- gsub("^cbn_ot$", "CBN", df_r$method)
##     df_r$method <- gsub("^mccbn$", "MCCBN", df_r$method)
##     df_r$method <- gsub("^mhn$", "MHN", df_r$method)
##     df_r$method <- gsub("^ot$", "OT", df_r$method)
##     df_r$method <- gsub("^td-cbn_ot$", "CBN_td", df_r$method)
##     df_r$method <- gsub("^td-mccbn$", "MCCBN_td", df_r$method)
##     df_r$method <- gsub("^td-mhn$", "MHN_td", df_r$method)
##     df_r$method <- gsub("^uw-cbn_ot$", "CBN_uw", df_r$method)
##     df_r$method <- gsub("^uw-mccbn$", "MCCBN_uw", df_r$method)
##     df_r$method <- gsub("^uw-ot$", "OT_uw", df_r$method)

##     df_r$typeLandscape <- gsub("^Local$", "Local maxima", df_r$typeLandscape)


    
##                                         # a few checks
##     if(FALSE){
##         u <- table(df_r[, c("id", "detect")])
##         uu <- as.data.frame(u)

##         v <- table(df_r[, c("id", "size_split")])
##         vv <- as.data.frame(v)

##         summary(as.vector(table(uu$id))) ## all a 3
##         table(uu$detect) ## all 1260

##         summary(as.vector(table(vv$id))) ## all a 3
##         table(vv$size_split) ## all 1260

##         w <- table(df_r[, c("id", "size_split", "detect")])
##         ww <- as.data.frame(w)

##         summary(as.vector(table(ww$id))) ## all a 9
##         table(ww$size_split) ## all 1260 * 3
##         table(ww$detect) ## all 1260 * 3

##         which_no_repl <- which(is.na(df_r$replicate))

##         table(df_r[which_no_repl, "method"]) ## all should be null
##         table(df_r$replicate) ## identical
##         summary(df_r$replicate) ## 1 to 5, and NAs as given above

##         tt <- table(df_r[, c("typeLandscape")])
##         tt ## the three landscapes

##         ss <- table(df_r[, c("id", "typeLandscape")])
##         sss <- as.data.frame(ss)

##         summary(as.vector(table(sss$id))) ## all a 3
##         summary(as.vector(table(sss[sss$Freq > 0, ]$id))) ## all a 1
##         table(sss$typeLandscape) ## all 1260
##         table(sss[sss$Freq > 0, ]$typeLandscape) ## all 420

##     }

##                                         # save output
##     cat("\n")
##     cat("Saving output")
##     cat("\n")

##     outFile_with_replicates <- file.path(saveDirectory,"table-replicates.rds")
##     system.time(
##         saveRDS(list(data = df_r,
##                      columnsExplained = columnsExplained),
##                 file = outFile_with_replicates,
##                 compress = FALSE)
##     )

## })




###  No longer used. Averages and array statistics computed separately in
###   average-and-array-statistics.R






## # aggregate replicates
## cat("\n")
## cat("Aggregating replicate statistics")
## cat("\n")
## pboptions(type="txt")
## df <- pblapply(df,
##                function(df) {
                 
##                  # aggregate stats for methods with replicates (all except "null")
##                  x <- df[df$cpm!="null",]
##                  x <- aggregate(x,by=list(x$id,
##                                           x$cpm,
##                                           x$size_split,
##                                           x$detect,
##                                           x$sourceGenotype),
##                                 FUN=function(x) {
##                                   if(all(is.na(x))) {
##                                     return(NA)
##                                   } else {
##                                     if(is.numeric(x)) return(mean(x,na.rm=T))
##                                     if(is.character(x)) return(paste(unique(x),collapse=" | "))
##                                     if(is.logical(x)) return(mean(x,na.rm=T))
##                                   }
##                                 })
##                  x <- x[,colnames(x) %in% colnames(df)]
                 
##                  # attach null model (size_split and detect columns are NA)
##                  x <- rbind(x,df[df$cpm=="null",])

##                  x$replicate <- NULL
##                  return(x)
                 
##                },
##                cl=detectCores())

## # combine data frames into single table
## cat("\n")
## cat("Combining data into single table")
## cat("\n")
## df <- do.call(rbind,df)

## # change "" genotype into "WT" or "root" and set accessibility to TRUE
## df$sourceGenotype[df$sourceGenotype==""] <- "root"

## # format "size_split" column
## df$size_split[df$size_split=="size_split_50"] <- "50"
## df$size_split[df$size_split=="size_split_200"] <- "200"
## df$size_split[df$size_split=="size_split_4000"] <- "4000"
## df$size_split <- as.numeric(df$size_split)

## # format "detect" column
## df$detect[df$detect=="detect_large"] <- "large"
## df$detect[df$detect=="detect_small"] <- "small"
## df$detect[df$detect=="detect_unif"] <- "unif"

## # load fitness landscape and sampling characteristics, set all "unif" to "uniform"
## load("./fitland/fl_sampl.RData")

## df$detect <- as.character(df$detect)
## df$detect[df$detect=="unif"] <- "uniform"
## fl_sampl$Detection <- as.character(fl_sampl$Detection)
## fl_sampl$Detection[fl_sampl$Detection=="unif"] <- "uniform"

## # merge tables
## cat("\n")
## cat("Attaching fitness landscape & sampling statistics")
## cat("\n")
## df <- merge(df,fl_sampl,
##             by.x=c("id","detect"),by.y=c("ID","Detection"),
##             all=T,suffixes=NULL)

## df <-
##     df[, c(
##         "id",
##         "cpm", "size_split",
##         "sqDiff", "sqDiff_eq", "js", "js_eq", "hellinger", "hellinger_eq", "spearman", "spearman_pval",
##         "flags",
##         "sourceGenotype", "sourceGenotype_nMut", "sourceGenotype_freqInPOM","sourceGenotype_accessible",
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
##         "w.1.", "w.2.", "w.3..","gamma"
##     )]

## colnames(df) <-
##     c(
##         "id",
##         "method","size_split",
##         "sqDiff", "sqDiff_eq", "js", "js_eq", "hellinger", "hellinger_eq", "spearman", "spearman_pval",
##         "flags",
##         "sourceGenotype","sourceGenotype_numMuts","sourceGenotype_freqInPOM","sourceGenotype_accessible",
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
##         "w1", "w2", "w3","gamma"
##     )

## ## explanation of each column
## columnsExplained <- list(
##   identification=list(
##     id="Data identifier"
##     ),
##   cpm_and_input=list(
##     method="Method used for next genotype prediction",
##     size_split="Sample size: size of the input matrix given to the method"
##     ),
##   statistics=list(
##     sqDiff="Square root of average of squared differences",
##     sqDiff_eq="Square root of average of squared differences (equiprobabilized)",
##     js="Jensen-Shannon distance (square root of the Jensen-Shannon divergence), in log base 2 units",
##     js_eq="Jensen-Shannon distance (equiprobabilized)",
##     hellinger="Hellinger distance",
##     hellinger_eq="Hellinger distance (equiprobabilized)",
##     spearman="Spearman's rank correlation value",
##     spearman_pval="P-value of the rank correlation"
##   ),
##   flags=list(
##     flags="Warnings shown when unfusing or rearranging genotype names"
##   ),
##   properties_of_source_genotype=list(
##     sourceGenotype="Name of the source genotype (if set to 'any', the statistical parameters correspond to the averages across all source genotypes, weighted by their frequency in the POM)",
##     sourceGenotype_numMuts="Number of mutations of the source genotype",
##     sourceGenotype_freqInPOM="Frequency of appearance of the source genotype in the POM",
##     sourceGenotype_accessible="Fraction of the replicates in which the source genotype was accessible according to the method"
##   ),
##   sampling=list(
##     detect="Detection regime: when tumors are sampled (large, small, uniform)",
##     sampledGenotypesDiversity="Diversity of sampled genotypes",
##     sampledGenotypesNumber="Number of unique genotypes in the sample",
##     sampledGenotypesAbove_0.01="Number of unique genotypes with a frequency above 0.01 present in the sample",
##     sampledGenotypesAbove_0.01="Number of unique genotypes with a frequency above 0.1 present in the sample",
##     numMuts_mean="Average number of mutations in the sample",
##     numMuts_median="Median number of mutations in the sample",
##     numMuts_var="Variance of the number of mutations in the sample",
##     numMuts_stDev="Standard deviation of the number of mutations in the sample",
##     numMuts_kurtosis="Kurtosis of the distribution of the number of mutations in the sample",
##     numMuts_skewness="Skewness of the distribution of the number of mutations in the sample",
##     pom_h="POM diversity",
##     lod_h="LOD diversity",
##     freq_most_freq_mean_no450="Frequency of the most frequent genotype in the sample",
##     how_many_gt_5p_mean_no450="How many genotypes have a frequency > 5% in the sample"
##   ),
##   simulation_and_evolutionary_process=list(
##     numAccessibleGenotypes="Number of accesible genotypes in the fitness landscape",
##     initSize="Initial number of wild-type cells in the simulation",
##     mutationRate="Mutation rate regime"
##   ),
##   fitness_landscape=list(
##     rnst="Fitness landscape identifier (redundant check: one-to-one between id and rnst)",
##     numGenes="Number of driver genes (7 or 10)",
##     typeLandscape="Type of fitness landscape",
##     numLocalPeaks="Number of local peaks (maxima) in the fitness landscape under the no-back mutation assumption",
##     numObservedPeaks="Number of local peaks in the landscape that are actually visited in the evolutionary simulations"
##   ),
##   epistasis=list(
##     epistMagn="Fraction of pairs of loci with magnitude epistasis in the landscape",
##     epistSign="Fraction of pairs of loci with sign epistasis",
##     epistRSign="Fraction of pairs of loci with reciprocal sign epistasis",
##     w1="Fourier expansion of the landscape: fraction of coefficients of order 1",
##     w2="Fourier expansion of the landscape: fraction of coefficients of order 2",
##     w3="Fourier expansion of the landscape: fraction of coefficients of order 3 or higher",
##     gamma="Correlation in fitness effects between genotypes that differ by one locus (Ferretti et al., 2016)"
##   )
## )

## # fix formatting
## cat("\n")
## cat("Applying format")
## cat("\n")
## df$method <- gsub("^caprese$", "CAPRESE", df$method)
## df$method <- gsub("^capri_aic$", "CAPRI_AIC", df$method)
## df$method <- gsub("^capri_bic$", "CAPRI_BIC", df$method)
## df$method <- gsub("^cbn_ot$", "CBN", df$method)
## df$method <- gsub("^mccbn$", "MCCBN", df$method)
## df$method <- gsub("^mhn$", "MHN", df$method)
## df$method <- gsub("^ot$", "OT", df$method)
## df$method <- gsub("^td-cbn_ot$", "CBN_td", df$method)
## df$method <- gsub("^td-mccbn$", "MCCBN_td", df$method)
## df$method <- gsub("^td-mhn$", "MHN_td", df$method)
## df$method <- gsub("^uw-cbn_ot$", "CBN_uw", df$method)
## df$method <- gsub("^uw-mccbn$", "MCCBN_uw", df$method)
## df$method <- gsub("^uw-ot$", "OT_uw", df$method)

## df$typeLandscape <- gsub("^Local$", "Local maxima", df$typeLandscape)

## # a few checks
## if(F){
##   u <- table(df[, c("id", "detect")])
##   uu <- as.data.frame(u)
##   summary(as.vector(table(uu$id))) ## all a 3
##   table(uu$detect) ## all 1260
## }

## # save output
## cat("\n")
## cat("Saving output")
## cat("\n")
## data <- df
## outFile <- file.path(saveDirectory,"table.rds")
## saveRDS(list(data=data,
##             columnsExplained=columnsExplained),
##         file=outFile)



