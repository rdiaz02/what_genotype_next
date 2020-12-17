## Create data set in "wide" format, for easier:

##  a) finding out best method
##  b) finding out best possible performance
##  c) paired t-tests et al.
##  d) pairs plots (see pairs-plots.R)

rm(list = ls())
## library(dplyr)
library(data.table)

## Yes, we could load each, transform, rm, gc(), but except for
## one, the rest are not that big and makes my life a lot simpler.

## Created in makeTable-to-diff-wrt-null-2.R
load("genotype_diff_null.RData")
load("array_diff_null.RData")
load("genotype_diff_null_no_average.RData")
load("num_mut_diff_null.RData")


## Convert to data.table by reference
setDT(genotype_diff_null)
setDT(genotype_diff_null_no_average)
setDT(array_diff_null)
setDT(num_mut_diff_null)



## To have a common set of columns that are used below
## Verify first
stopifnot(is.null(genotype_diff_null$replicate))
genotype_diff_null$replicate <- NA


stopifnot(is.null(array_diff_null$sourceGenotype)) 
stopifnot(is.null(array_diff_null$sourceGenotype_nMut)) 
stopifnot(is.null(array_diff_null$sourceGenotype_accessible)) 
stopifnot(is.null(array_diff_null$replicate)) 

array_diff_null$sourceGenotype <- "any"
array_diff_null$sourceGenotype_nMut <- NA
array_diff_null$sourceGenotype_accessible <- NA
array_diff_null$replicate <- NA

## partial matching, though we could use more idiomatic DT expressions
stopifnot(is.null(num_mut_diff_null[["sourceGenotype"]])) 
stopifnot(is.null(num_mut_diff_null$sourceGenotype_accessible)) 
stopifnot(is.null(num_mut_diff_null$replicate)) 

num_mut_diff_null$sourceGenotype <- NA
num_mut_diff_null$sourceGenotype_accessible <- NA
num_mut_diff_null$replicate <- NA



## df where columns are response vars in each method
## Helpful for paired t-tests and paired plots of method results
f_wide_out <- function(x) {
    ## This is VERT hackish, but things are tailored for the different
    ## input data

    nn <- deparse(substitute(x))

    if(nn == "array_diff_null") input <- "array"
    if(nn == "genotype_diff_null") input <- "genotype"
    if(nn == "genotype_diff_null_no_average") input <- "genotype_no_average"
    if(nn == "num_mut_diff_null") input <- "num_mut"
 
    
    if(input == "array") {
        stopifnot(is.na(x$sourceGenotype_nMut))
        stopifnot(is.na(x$sourceGenotype_accessible))
    }

    if(input == "genotype") {
        stopifnot(!is.na(x$sourceGenotype_nMut))
        ## accessible is sometimes NA
        ## double check this
        naaccess <- which(is.na(x$sourceGenotype_accessible))
        ## We are dealing with the data we think
        stopifnot(length(naaccess) < nrow(x))
        ## When accessible is na the js is na too. Those are
        ## failures in analysis
        nanjs <- which(is.nan(x$js))
        stopifnot(identical(nanjs, naaccess))
    }
    
    
    
    ## ## Figure out the data we are dealing with
    ## ## but check our guess is correct
    ## ## More checks later when we look at column names
    ## if(all(x$sourceGenotype == "any")) {
    ##     input <- "array"
    ##     stopifnot(is.na(x$sourceGenotype_nMut))
    ##     stopifnot(is.na(x$sourceGenotype_accessible))
    ## } else if (all(is.na(x$sourceGenotype))) {
    ##     input <- "num_mut"
    ## } else if (all(is.na(x$sourceGenotype)) &&
    ##            ) {
    ##     input <- "genotype"
    ##     stopifnot(!is.na(x$sourceGenotype_nMut))
    ##     ## accessible is sometimes NA
    ##     ## double check this
    ##     naaccess <- which(is.na(x$sourceGenotype_accessible))
    ##     ## We are dealing with the data we think
    ##     stopifnot(length(naaccess) < nrow(x))
    ##     ## When accessible is na the js is na too. Those are
    ##     ## failures in analysis
    ##     nanjs <- which(is.nan(x$js))
    ##     stopifnot(identical(nanjs, naaccess))
    ## }
    
    common_columns <- c(
        ## Columns to match
        "id", "sourceGenotype",
        "detect", "size_split",
        "sourceGenotype_nMut",
        "replicate",
        ## Redundant, but keep to check
        "rnst",
        ## Descriptors
        "numGenes", "initSize",
        "typeLandscape", "mutationRate",
        ## Genotype descriptors
        ## "sourceGenotype_accessible", ## FIXME: wrong!!! accessible is a method property ## see bottom
        ## Landscape and sample stats
        "sampledGenotypesDiversity",
        "sampledGenotypesNumber",
        "sampledGenotypesGenesAbove_0.01",
        "sampledGenotypesGenesAbove_0.1",
        "numMuts_mean",
        "numMuts_median", "numMuts_var",
        "numMuts_stDev", "numMuts_kurtosis",
        "numMuts_skewness", "pom_h",
        "lod_h", "freq_most_freq_mean_no450",
        "how_many_gt_5p_mean_no450",
        "numAccessibleGenotypes",
        "numLocalPeaks",
        "numObservedPeaks",
        "diversityObservedPeaks",
        "epistMagn",
        "epistSign", "epistRSign",
        "w1", "w2",
        "w3", "gamma")
    
    if(input == "array") {
        response_columns <- c("js_w_real", "js_w_sampl",
                              "js_w_real_diff", "js_w_sampl_diff")
        common_columns <- c(common_columns,
                            "js_w_real_null",
                            "js_w_sampl_null"
                            )
        stopifnot(length(setdiff(response_columns, colnames(x))) == 0)
        stopifnot(length(setdiff(common_columns, colnames(x))) == 0)
    }
    if(input == "num_mut") {
        response_columns <- c("js_w_real", "js_w_sampl",
                              "js_w_real_diff", "js_w_sampl_diff")
        common_columns <- c(common_columns,
                            "js_w_real_null",
                            "js_w_sampl_null"
                            )
        stopifnot(length(setdiff(response_columns, colnames(x))) == 0)
        stopifnot(length(setdiff(common_columns, colnames(x))) == 0)
    }
    if(input == "genotype") {
        response_columns <- c("js", "js_diff")
        common_columns <- c(common_columns,
                            "sampledFreq",                     
                            "sampledFreq_unique",              
                            "sampledProp",                     
                            "sampledProp_unique",              
                            "trueFreq",                        
                            "trueFreq_unique",                 
                            "trueProp",                        
                            "trueProp_unique",                 
                            "observedFreq",                    
                            "observedFreq_unique",             
                            "observedProp",                    
                            "observedProp_unique",             
                            "fitnessRank",                     
                            "fitnessRank_unique",              
                            "fitnessRankNoZero",               
                            "fitnessRankNoZero_unique",        
                            "freqLocalMax",                    
                            "freqLocalMax_unique",             
                            "propLocalMax",                    
                            "propLocalMax_unique",
                            "sourceGenotype_accessible_unique",
                            "js_null"
                            )
        stopifnot(length(setdiff(response_columns, colnames(x))) == 0)
        stopifnot(length(setdiff(common_columns, colnames(x))) == 0)
    }
    if(input == "genotype_no_average") {
        response_columns <- c("js", "js_diff")
        common_columns <- c(common_columns,
                            "sampledFreq",                     
                            "sampledProp",                     
                            "trueFreq",                        
                            "trueProp",                        
                            "observedFreq",                    
                            "observedProp",                    
                            "fitnessRank",                     
                            "fitnessRankNoZero",               
                            "freqLocalMax",                    
                            "propLocalMax",                    
                            "js_null"
                            )
        stopifnot(length(setdiff(response_columns, colnames(x))) == 0)
        stopifnot(length(setdiff(common_columns, colnames(x))) == 0)
    }

    
    x <- x[method != "null", ]
    x$method <- factor(x$method,
                       levels=c(
                           "MHN_td",
                           "CBN_td",
                           "MCCBN_td",
                           "MHN",
                           "CBN",
                           "MCCBN",
                           "CBN_uw",
                           "MCCBN_uw",
                           "CAPRI_AIC",
                           "CAPRI_BIC",
                           "OT",
                           "OT_uw",
                           "CAPRESE"
                       )
                       )
    data_tmp_s <- split(x, x$method)

    ## Instead of split, we could process using by? Too much of a mess?
    
    ## Obtain dfs with only the response and identifiers
    ##   but rename columns

    ## js, js_eq, js_diff, js_eq_diff
    ## id, sourceGenotype,
    ## detect, size_split
    ##  method
    ##  To check: rnst

    ## Unique output
    unique_out_cols <- c(
        ## Method
        "method",
        ## Columns to match
        "id", "sourceGenotype", "detect", "size_split",
        "sourceGenotype_nMut",
        "replicate",
        ## Redundant, for check
        "rnst",
        ## Responses
        response_columns
    )
    
    data_tmp_s_unique <-
        lapply(data_tmp_s,
               function(x) {
                   tmp <- x[, ..unique_out_cols]
                   ## Rename response columns
                   colnames(tmp)[c(9:(9 + length(response_columns) - 1))] <-
                       paste0(unique(tmp$method), "_",
                              colnames(tmp)[c(9:(9 + length(response_columns) - 1))])
                   ## No need for method anymore
                   tmp[, method:=NULL]
                   return(tmp)
               }
               )

    ## Extract common columns. Any method will do for that
    data_tmp_common <- data_tmp_s[[1]][, ..common_columns]

    ## Join all method tables
   
    all_m <- Reduce(
        function(...) merge.data.table(..., by = c("id", "sourceGenotype",
                                                   "detect", "size_split",
                                                   "sourceGenotype_nMut",
                                                   "replicate",
                                                   "rnst" ## ok, using it
                                                   ## to avoid warnings from DT
                                                   )),
        data_tmp_s_unique)


    ## Checks
    dim(all_m)
    nrows <- lapply(data_tmp_s_unique, nrow)
    stopifnot(length(unique(nrows)) == 1)
    stopifnot(unique(nrows) == nrow(all_m))

    ## ## More checks: use rnst
    ## dummy2 <- all_m[, grep("rnst", colnames(all_m), fixed = TRUE)]
    ## for(i in 2:ncol(dummy2)) stopifnot(dummy2[, 1] == dummy2[, i])
    ## rm(dummy2)
    ## gc()

    ## rm the rnst columns
    ## all_m <- all_m[, -(grep("rnst", colnames(all_m), fixed = TRUE)[-1])]
    stopifnot(nrow(all_m) == nrow(data_tmp_common))

    ## Join the method-specific output with the common columns
    all_m_c <- merge(all_m, data_tmp_common,
                     by = c("id", "sourceGenotype",
                            "detect", "size_split",
                            "sourceGenotype_nMut",
                            "replicate",
                            "rnst"
                            ))
    
    stopifnot(nrow(all_m) == nrow(all_m_c))
    ## dummy2 <- all_m_c[, grep("rnst", colnames(all_m_c), fixed = TRUE)]
    ## stopifnot(dummy2[, 1] == dummy2[, 2])
    ## rm(dummy2)
    gc()
    return(all_m_c)
}

wide_array <- f_wide_out(array_diff_null)
wide_genotype   <- f_wide_out(genotype_diff_null)
wide_num_mut <- f_wide_out(num_mut_diff_null)
wide_genotype_no_average <- f_wide_out(genotype_diff_null_no_average)



save(file = "wide_array.RData", wide_array, compress = FALSE)
save(file = "wide_genotype.RData", wide_genotype, compress = FALSE)
save(file = "wide_genotype_no_average.RData",
     wide_genotype_no_average, compress = FALSE)
save(file = "wide_num_mut.RData", wide_num_mut, compress = FALSE)



## This gives and error. No problem I want it to stop
stop()






## ## Former dplyr code

## ## To have a common set of columns that are used below
## ## Verify first
## stopifnot(is.null(genotype_diff_null$replicate))
## genotype_diff_null$replicate <- NA
## stopifnot(is.null(array_diff_null$sourceGenotype)) 
## stopifnot(is.null(array_diff_null$sourceGenotype_nMut)) 
## stopifnot(is.null(array_diff_null$sourceGenotype_accessible)) 
## stopifnot(is.null(array_diff_null$replicate)) 

## array_diff_null$sourceGenotype <- "any"
## array_diff_null$sourceGenotype_nMut <- NA
## array_diff_null$sourceGenotype_accessible <- NA
## array_diff_null$replicate <- NA

## ## partial matching
## stopifnot(is.null(num_mut_diff_null[["sourceGenotype"]])) 
## stopifnot(is.null(num_mut_diff_null$sourceGenotype_accessible)) 
## stopifnot(is.null(num_mut_diff_null$replicate)) 

## num_mut_diff_null$sourceGenotype <- NA
## num_mut_diff_null$sourceGenotype_accessible <- NA
## num_mut_diff_null$replicate <- NA





## f_wide_out <- function(x) {

##     ## This is VERT hackish, but things are tailored for the different
##     ## input data

##     nn <- deparse(substitute(x))

##     if(nn == "array_diff_null") input <- "array"
##     if(nn == "genotype_diff_null") input <- "genotype"
##     if(nn == "genotype_diff_null_no_average") input <- "genotype_no_average"
##     if(nn == "num_mut_diff_null") input <- "num_mut"
 
    
##     if(input == "array") {
##         stopifnot(is.na(x$sourceGenotype_nMut))
##         stopifnot(is.na(x$sourceGenotype_accessible))
##     }

##     if(input == "genotype") {
##         stopifnot(!is.na(x$sourceGenotype_nMut))
##         ## accessible is sometimes NA
##         ## double check this
##         naaccess <- which(is.na(x$sourceGenotype_accessible))
##         ## We are dealing with the data we think
##         stopifnot(length(naaccess) < nrow(x))
##         ## When accessible is na the js is na too. Those are
##         ## failures in analysis
##         nanjs <- which(is.nan(x$js))
##         stopifnot(identical(nanjs, naaccess))
##     }
    
    
    
##     ## ## Figure out the data we are dealing with
##     ## ## but check our guess is correct
##     ## ## More checks later when we look at column names
##     ## if(all(x$sourceGenotype == "any")) {
##     ##     input <- "array"
##     ##     stopifnot(is.na(x$sourceGenotype_nMut))
##     ##     stopifnot(is.na(x$sourceGenotype_accessible))
##     ## } else if (all(is.na(x$sourceGenotype))) {
##     ##     input <- "num_mut"
##     ## } else if (all(is.na(x$sourceGenotype)) &&
##     ##            ) {
##     ##     input <- "genotype"
##     ##     stopifnot(!is.na(x$sourceGenotype_nMut))
##     ##     ## accessible is sometimes NA
##     ##     ## double check this
##     ##     naaccess <- which(is.na(x$sourceGenotype_accessible))
##     ##     ## We are dealing with the data we think
##     ##     stopifnot(length(naaccess) < nrow(x))
##     ##     ## When accessible is na the js is na too. Those are
##     ##     ## failures in analysis
##     ##     nanjs <- which(is.nan(x$js))
##     ##     stopifnot(identical(nanjs, naaccess))
##     ## }
    
##     common_columns <- c(
##         ## Columns to match
##         "id", "sourceGenotype",
##         "detect", "size_split",
##         "sourceGenotype_nMut",
##         "replicate",
##         ## Redundant, but keep to check
##         "rnst",
##         ## Descriptors
##         "numGenes", "initSize",
##         "typeLandscape", "mutationRate",
##         ## Genotype descriptors
##         "sourceGenotype_accessible",
##         ## Landscape and sample stats
##         "sampledGenotypesDiversity",
##         "sampledGenotypesNumber",
##         "sampledGenotypesGenesAbove_0.01",
##         "sampledGenotypesGenesAbove_0.1",
##         "numMuts_mean",
##         "numMuts_median", "numMuts_var",
##         "numMuts_stDev", "numMuts_kurtosis",
##         "numMuts_skewness", "pom_h",
##         "lod_h", "freq_most_freq_mean_no450",
##         "how_many_gt_5p_mean_no450",
##         "numAccessibleGenotypes",
##         "numLocalPeaks",
##         "numObservedPeaks",
##         "diversityObservedPeaks",
##         "epistMagn",
##         "epistSign", "epistRSign",
##         "w1", "w2",
##         "w3", "gamma")
    
##     if(input == "array") {
##         response_columns <- c("js_w_real", "js_w_sampl",
##                               "js_w_real_diff", "js_w_sampl_diff")
##         common_columns <- c(common_columns,
##                             "js_w_real_null",
##                             "js_w_sampl_null"
##                             )
##         stopifnot(length(setdiff(response_columns, colnames(x))) == 0)
##         stopifnot(length(setdiff(common_columns, colnames(x))) == 0)
##     }
##     if(input == "num_mut") {
##         response_columns <- c("js_w_real", "js_w_sampl",
##                               "js_w_real_diff", "js_w_sampl_diff")
##         common_columns <- c(common_columns,
##                             "js_w_real_null",
##                             "js_w_sampl_null"
##                             )
##         stopifnot(length(setdiff(response_columns, colnames(x))) == 0)
##         stopifnot(length(setdiff(common_columns, colnames(x))) == 0)
##     }
##     if(input == "genotype") {
##         response_columns <- c("js", "js_diff")
##         common_columns <- c(common_columns,
##                             "sampledFreq",                     
##                             "sampledFreq_unique",              
##                             "sampledProp",                     
##                             "sampledProp_unique",              
##                             "trueFreq",                        
##                             "trueFreq_unique",                 
##                             "trueProp",                        
##                             "trueProp_unique",                 
##                             "observedFreq",                    
##                             "observedFreq_unique",             
##                             "observedProp",                    
##                             "observedProp_unique",             
##                             "fitnessRank",                     
##                             "fitnessRank_unique",              
##                             "fitnessRankNoZero",               
##                             "fitnessRankNoZero_unique",        
##                             "freqLocalMax",                    
##                             "freqLocalMax_unique",             
##                             "propLocalMax",                    
##                             "propLocalMax_unique",
##                             "sourceGenotype_accessible",                            
##                             "sourceGenotype_accessible_unique",
##                             "diversityObservedPeaks",     
##                             "js_null"
##                             )
##         stopifnot(length(setdiff(response_columns, colnames(x))) == 0)
##         stopifnot(length(setdiff(common_columns, colnames(x))) == 0)
##     }
##     if(input == "genotype_no_average") {
##         response_columns <- c("js", "js_diff")
##         common_columns <- c(common_columns,
##                             "sampledFreq",                     
##                             "sampledProp",                     
##                             "trueFreq",                        
##                             "trueProp",                        
##                             "observedFreq",                    
##                             "observedProp",                    
##                             "fitnessRank",                     
##                             "fitnessRankNoZero",               
##                             "freqLocalMax",                    
##                             "propLocalMax",                    
##                             "sourceGenotype_accessible",
##                             "diversityObservedPeaks",     
##                             "js_null"
##                             )
##         stopifnot(length(setdiff(response_columns, colnames(x))) == 0)
##         stopifnot(length(setdiff(common_columns, colnames(x))) == 0)
##     }

    
##     x <- dplyr::filter(x, method != "null")
##     x$method <- factor(x$method,
##                        levels=c(
##                            "MHN_td",
##                            "CBN_td",
##                            "MCCBN_td",
##                            "MHN",
##                            "CBN",
##                            "MCCBN",
##                            "CBN_uw",
##                            "MCCBN_uw",
##                            "CAPRI_AIC",
##                            "CAPRI_BIC",
##                            "OT",
##                            "OT_uw",
##                            "CAPRESE"
##                        )
##                        )
##     data_tmp_s <- split(x, x$method)

##     ## Obtain dfs with only the response and identifiers
##     ##   but rename columns

##     ## js, js_eq, js_diff, js_eq_diff
##     ## id, sourceGenotype,
##     ## detect, size_split
##     ##  method
##     ##  To check: rnst

##     ## Unique output
##     data_tmp_s_unique <-
##         lapply(data_tmp_s,
##                function(x) {
##                    tmp <- x[, c(
##                        ## Method
##                        "method",
##                        ## Columns to match
##                        "id", "sourceGenotype", "detect", "size_split",
##                        "sourceGenotype_nMut",
##                        "replicate",
##                        ## Redundant, for check
##                        "rnst",
##                        ## Responses
##                        response_columns
##                    )]
##                    ## Rename response columns
##                    colnames(tmp)[c(9:(9 + length(response_columns) - 1))] <-
##                        paste0(unique(tmp$method), "_",
##                               colnames(tmp)[c(9:(9 + length(response_columns) - 1))])
##                    ## No need for method anymore
##                    tmp <- tmp[, -1]
##                    return(tmp)
##                }
##                )

##     ## Extract common columns. Any method will do for that
##     data_tmp_common <- data_tmp_s[[1]][, common_columns]

##     ## Join all method tables
##     all_m <- Reduce(function(...) left_join(..., by = c("id", "sourceGenotype",
##                                                         "detect", "size_split",
##                                                         "sourceGenotype_nMut",
##                                                         "replicate"
##                                                         )),
##                     data_tmp_s_unique)
##     ## Checks
##     dim(all_m)
##     nrows <- lapply(data_tmp_s_unique, nrow)
##     stopifnot(length(unique(nrows)) == 1)
##     stopifnot(unique(nrows) == nrow(all_m))
##     ## More checks: use rnst
##     dummy2 <- all_m[, grep("rnst", colnames(all_m), fixed = TRUE)]
##     for(i in 2:ncol(dummy2)) stopifnot(dummy2[, 1] == dummy2[, i])
##     rm(dummy2)
##     gc()

##     ## rm the rnst columns
##     all_m <- all_m[, -(grep("rnst", colnames(all_m), fixed = TRUE)[-1])]
##     stopifnot(nrow(all_m) == nrow(data_tmp_common))

##     ## Join the method-specific output with the common columns
##     all_m_c <- dplyr::full_join(all_m, data_tmp_common,
##                                 by = c("id", "sourceGenotype",
##                                        "detect", "size_split",
##                                        "sourceGenotype_nMut",
##                                        "replicate"
##                                        ))
    
##     stopifnot(nrow(all_m) == nrow(all_m_c))
##     dummy2 <- all_m_c[, grep("rnst", colnames(all_m_c), fixed = TRUE)]
##     stopifnot(dummy2[, 1] == dummy2[, 2])
##     rm(dummy2)
##     gc()

##     return(all_m_c)
## }



## FIXME: the "sourceGenotype_accessible" issue

## It turns out that for the wide, we have been using the value from
## MHN. These can be verified this way

## load("genotype_diff_null_no_average.RData")
## load("wide_genotype_no_average_rank.RData")
## aggregate(sourceGenotype_accessible ~ method,
##           FUN = function(x) sum(x, na.rm = TRUE),
##           data = genotype_diff_null_no_average)

## ##       method sourceGenotype_accessible
## ## 1    CAPRESE                   1449137
## ## 2  CAPRI_AIC                    757353
## ## 3  CAPRI_BIC                    939277
## ## 4        CBN                   1348033
## ## 5     CBN_td                   1348033
## ## 6     CBN_uw                   1348033
## ## 7      MCCBN                   1043920
## ## 8   MCCBN_td                   1043920
## ## 9   MCCBN_uw                   1043920
## ## 10       MHN                   2618524
## ## 11    MHN_td                   2618524
## ## 12      null                    200592
## ## 13        OT                   1207031
## ## 14     OT_uw                   1207031

## sum(wide_genotype_no_average_rank$sourceGenotype_accessible, na.rm = TRUE)
## ## [1] 2618524
