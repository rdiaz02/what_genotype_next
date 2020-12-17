## Create data with the differences w.r.t. null incorporated as
## variables
date()
rm(list = ls())

library(dplyr)

###########################################
##
## Genotype-by-genotype data, over replicate
##
###########################################


df <- readRDS("data_no_any_over_replicate.rds")
data <- df$data

columnsExplained <- df$columnsExplained

## Getting the difference w.r.t null

nullm <- dplyr::filter(df$data, method == "null")

## Note
table(df$data$method)/nrow(nullm)
## 3
## where 3 = 3 size_split
## Recall that in "merge-additional-info.R" we replicated the null for the three
## detections. See around line 156
nullmljs <- nullm[nullm$detect == "large", "js"]
nullmsjs <- nullm[nullm$detect == "small", "js"]
nullmujs <- nullm[nullm$detect == "uniform", "js"]
stopifnot(identical(nullmljs, nullmsjs))
stopifnot(identical(nullmljs, nullmujs))


## detect and size_split make no difference for null model
## recall no value for size_split
stopifnot(is.na(nullm$size_split))

nullm2 <- nullm[nullm$detect == "large",
                c("id", "sourceGenotype",
                    ## Next are for checking
                    "typeLandscape",
                    "numGenes", 
                    "sourceGenotype_nMut",
                  ## "sourceGenotype_freqInPOM",
                    "method", ## rename this
                    ## The response to subtract
                    "js")]

colnames(nullm2)[3:ncol(nullm2)] <- paste0(colnames(nullm2)[3:ncol(nullm2)], "_null")
## rm(nullm)
rm(df)
gc()

data2 <- dplyr::full_join(data, nullm2, by = c("id", "sourceGenotype"))
tmp <- dplyr::left_join(data, nullm2, by = c("id", "sourceGenotype"))


## checks
stopifnot(nrow(data) == nrow(data2))

stopifnot(identical(tmp, data2))
rm(tmp)
gc()
stopifnot(table(data2$method_null) == nrow(data))
stopifnot(nrow(data2) == nrow(data))
stopifnot(identical(data2$typeLandscape, data2$typeLandscape_null))
stopifnot(identical(data2$numGenes, data2$numGenes_null))
stopifnot(identical(data2$sourceGenotype_numMuts, data2$sourceGenotype_numMuts_null))
## silly: these columns not present
## stopifnot(identical(data2$sourceGenotype_freqInPOM, data2$sourceGenotype_freqInPOM_null))
stopifnot(data2$method_null == "null")

data2$js_diff <- data2$js_null - data2$js
## rm duplicated columns just for checking
data2 <- data2[, -c(60:63)]

stopifnot(dplyr::filter(data2, method == "null")$js_diff == 0)

genotype_diff_null <- data2

save(file = "genotype_diff_null.RData", genotype_diff_null,
     compress = FALSE)


rm(nullm, df, data, nullm2)
rm(list = ls())
gc()


###############################################
##
##       Array data
##
###############################################
rm(list = ls())
dfa <- readRDS("array_statistics_over_replicate.rds")
da <- dfa$data


nullm <- dplyr::filter(da, method == "null")

## Note
table(da$method)/nrow(nullm)
## 3
## where 3 = 3 size_split


## detect makes no difference for js_w_real, but makes a difference
## for js_w_sampl

## recall no value for size_split
stopifnot(is.na(nullm$size_split))

## checks
nullmljs <- nullm[nullm$detect == "large", "js_w_real"]
nullmsjs <- nullm[nullm$detect == "small", "js_w_real"]
nullmujs <- nullm[nullm$detect == "uniform", "js_w_real"]
stopifnot(identical(nullmljs, nullmsjs))
stopifnot(identical(nullmljs, nullmujs))

nullmljs <- nullm[nullm$detect == "large", "js_w_sampl"]
nullmsjs <- nullm[nullm$detect == "small", "js_w_sampl"]
nullmujs <- nullm[nullm$detect == "uniform", "js_w_sampl"]
stopifnot(!identical(nullmljs, nullmsjs))
stopifnot(!identical(nullmljs, nullmujs))




nullm2 <- nullm[,
                c("id", "detect",
                    ## Next are for checking
                    "typeLandscape",
                    "numGenes", 
                    "method", ## rename this
                    ## The response to subtract
                    "js_w_real", "js_w_sampl")]

colnames(nullm2)[3:ncol(nullm2)] <-
    paste0(colnames(nullm2)[3:ncol(nullm2)], "_null")

dataA <- dplyr::full_join(da, nullm2, by = c("id", "detect"))
tmp <- dplyr::left_join(da, nullm2, by = c("id", "detect"))
## checks
stopifnot(nrow(da) == nrow(dataA))

stopifnot(identical(tmp, dataA))
rm(tmp)
gc()
stopifnot(table(dataA$method_null) == nrow(da))
stopifnot(nrow(dataA) == nrow(da))
stopifnot(identical(dataA$typeLandscape, dataA$typeLandscape_null))
stopifnot(identical(dataA$numGenes, dataA$numGenes_null))
stopifnot(dataA$method_null == "null")


dataA$js_w_real_diff <- dataA$js_w_real_null - dataA$js_w_real
dataA$js_w_sampl_diff <- dataA$js_w_sampl_null - dataA$js_w_sampl

stopifnot(dplyr::filter(dataA, method == "null")$js_w_real_diff == 0)
stopifnot(dplyr::filter(dataA, method == "null")$js_w_sampl_diff == 0)

## rm duplicated columns just for checking
dataA <- dataA[-c(37:39)]

array_diff_null <- dataA

save(file = "array_diff_null.RData", array_diff_null,
     compress = FALSE)


## Yes, the "js_null" or "js_sampl_null" column is redundant for those which
## are null. And could be obtained for the rest by diff + js but this it is
## simpler to leave it.





###########################################
##
## Genotype-by-genotype data, not averaged over replicate
##
###########################################

rm(data2, tmp, nullm, nullm2)
rm(list = ls())

df <- readRDS("data_no_any.rds")
data <- df$data

columnsExplained <- df$columnsExplained

## Getting the difference w.r.t null

nullm <- dplyr::filter(df$data, method == "null")

## Note
table(df$data$method)/nrow(nullm)
## 15
## where 3 = 3 size_split* 5 replicates
## Recall that in "merge-additional-info.R" we replicated the null for the three
## detections. See around line 156

## detect and size_split and replicate make no difference for null model

## detect and size_split make no difference for null model
## recall no value for size_split or replicate
stopifnot(is.na(nullm$size_split))
stopifnot(is.na(nullm$replicate))

nullmljs <- nullm[nullm$detect == "large", "js"]
nullmsjs <- nullm[nullm$detect == "small", "js"]
nullmujs <- nullm[nullm$detect == "uniform", "js"]
stopifnot(identical(nullmljs, nullmsjs))
stopifnot(identical(nullmljs, nullmujs))





nullm2 <- nullm[nullm$detect == "large",
                c("id", "sourceGenotype",
                    ## Next are for checking
                    "typeLandscape",
                    "numGenes", 
                    "sourceGenotype_nMut",
                  ## "sourceGenotype_freqInPOM",
                    "method", ## rename this
                    ## The response to subtract
                    "js")]

colnames(nullm2)[3:ncol(nullm2)] <- paste0(colnames(nullm2)[3:ncol(nullm2)], "_null")
## rm(nullm)
rm(df)
gc()

data2 <- dplyr::full_join(data, nullm2, by = c("id", "sourceGenotype"))
tmp <- dplyr::left_join(data, nullm2, by = c("id", "sourceGenotype"))


## checks
stopifnot(nrow(data) == nrow(data2))
stopifnot(identical(tmp, data2))
rm(tmp)
gc()
stopifnot(table(data2$method_null) == nrow(data))
stopifnot(nrow(data2) == nrow(data))
stopifnot(identical(data2$typeLandscape, data2$typeLandscape_null))
stopifnot(identical(data2$numGenes, data2$numGenes_null))
stopifnot(identical(data2$sourceGenotype_numMuts, data2$sourceGenotype_numMuts_null))
## stopifnot(identical(data2$sourceGenotype_freqInPOM, data2$sourceGenotype_freqInPOM_null))
stopifnot(data2$method_null == "null")

data2$js_diff <- data2$js_null - data2$js
## rm duplicated columns just for checking
data2 <- data2[, -c(59:62)]

stopifnot(dplyr::filter(data2, method == "null")$js_diff == 0)

genotype_diff_null_no_average <- data2

save(file = "genotype_diff_null_no_average.RData",
     genotype_diff_null_no_average,
     compress = FALSE)


rm(nullm, df, data, nullm2)
rm(list = ls())
gc()



###########################################
##
## Number mutations, over replicate
##
###########################################

rm(list = ls())
dfa <- readRDS("num_mut_statistics_over_replicate.rds")
da <- dfa$data


nullm <- dplyr::filter(da, method == "null")

## Note
table(da$method)/nrow(nullm)
## 3
## where 3 = 3 size_split


## detect makes no difference for js_w_real, but makes a difference
## for js_w_sampl
## recall no value for size_split

## checks
stopifnot(is.na(nullm$size_split))

nullmljs <- nullm[nullm$detect == "large", "js_w_real"]
nullmsjs <- nullm[nullm$detect == "small", "js_w_real"]
nullmujs <- nullm[nullm$detect == "uniform", "js_w_real"]
stopifnot(identical(nullmljs, nullmsjs))
stopifnot(identical(nullmljs, nullmujs))

nullmljs <- nullm[nullm$detect == "large", "js_w_sampl"]
nullmsjs <- nullm[nullm$detect == "small", "js_w_sampl"]
nullmujs <- nullm[nullm$detect == "uniform", "js_w_sampl"]
stopifnot(!identical(nullmljs, nullmsjs))
stopifnot(!identical(nullmljs, nullmujs))


nullm2 <- nullm[,
                  c("id", "sourceGenotype_nMut",
                    "detect",
                    ## Next are for checking
                    "typeLandscape",
                    "numGenes", 
                    "method", ## rename this
                    ## The response to subtract
                    "js_w_real", "js_w_sampl")]

colnames(nullm2)[4:ncol(nullm2)] <-
    paste0(colnames(nullm2)[4:ncol(nullm2)], "_null")

dataA <- dplyr::full_join(da, nullm2, by = c("id", "sourceGenotype_nMut", "detect"))
tmp <- dplyr::left_join(da, nullm2, by = c("id", "sourceGenotype_nMut", "detect"))
## checks
stopifnot(identical(tmp, dataA))
rm(tmp)
gc()
stopifnot(table(dataA$method_null) == nrow(da))
stopifnot(nrow(dataA) == nrow(da))
stopifnot(identical(dataA$typeLandscape, dataA$typeLandscape_null))
stopifnot(identical(dataA$numGenes, dataA$numGenes_null))
stopifnot(dataA$method_null == "null")


dataA$js_w_real_diff <- dataA$js_w_real_null - dataA$js_w_real
dataA$js_w_sampl_diff <- dataA$js_w_sampl_null - dataA$js_w_sampl

stopifnot(dplyr::filter(dataA, method == "null")$js_w_real_diff == 0)
stopifnot(na.omit(dplyr::filter(dataA, method == "null")$js_w_sampl_diff) == 0)

## rm duplicated columns just for checking
dataA <- dataA[-c(38:40)]

num_mut_diff_null <- dataA

save(file = "num_mut_diff_null.RData", num_mut_diff_null,
     compress = FALSE)


date()



