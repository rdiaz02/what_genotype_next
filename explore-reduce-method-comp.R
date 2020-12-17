library(data.table)
setDTthreads(threads = 0)

## Created in merge-compsfixNulls-wide_g.R
load("./plots-glms/wide_g_method_comp.RData")

df <- wide_g_method_comp

colnames(df)

## This is dangerous
colsnum <- colnames(df)[25:52]
df[, (colsnum) := lapply(.SD, as.numeric), .SDcols = colsnum]
## All of them generate
## 1: In lapply(.SD, as.numeric) : NAs introduced by coercion

## 1460 NAs

sum(is.na(df$"CBN -vs- MHN"))

100 * sum(is.na(df$"CBN -vs- MHN"))/nrow(df) ## 0.06

lapply(df[, 25:52], function(x) sum(is.na(x)))

somenas <- which(is.na(df$"CBN -vs- MHN"))[c(1, 15, 1000, 1460)]

df[somenas, 1:27]

## OK, only those from the 16 cases where MCCBN had an NA
dfnona <- df[!is.na(df$min_js), ]
lapply(dfnona[, 25:52], function(x) sum(is.na(x)))

## 0.77
cor.test(dfnona$"CBN -vs- MHN", dfnona$"MCCBN -vs- MHN", method = "spearman")

## -0.0024
cor.test(dfnona$"CBN -vs- MCCBN", dfnona$"MCCBN -vs- MHN", method = "spearman")

## -0.1853
cor.test(dfnona$"CBN -vs- MHN", dfnona$"CBN -vs- MCCBN", method = "spearman")


## 0.863
cor.test(dfnona$"CAPRESE -vs- MHN", dfnona$"MHN -vs- OT", method = "spearman")

## 0.8237
cor.test(dfnona$"CAPRESE -vs- CBN", dfnona$"CBN -vs- OT", method = "spearman")

## - 0.25
cor.test(dfnona$"CAPRESE -vs- MHN", dfnona$"CAPRESE -vs- OT", method = "spearman")

## 0.45
cor.test(dfnona$"MHN -vs- MHN_td", dfnona$"CBN -vs- CBN_td", method = "spearman")

## 0.6211
cor.test(dfnona$"MCCBN -vs- MCCBN_td", dfnona$"CBN -vs- CBN_td", method = "spearman")

## 0.35
cor.test(dfnona$"MHN -vs- MHN_td", dfnona$"MCCBN -vs- MCCBN_td", method = "spearman")

## 0.3353
cor.test(dfnona$"CBN_td -vs- MHN", dfnona$"CBN -vs- MHN_td", method = "spearman")



## 0.19
cor.test(dfnona$"MHN -vs- MHN_td", dfnona$"CBN -vs- MHN", method = "spearman")


## Need to rename columns
colnames(dfnona) <- gsub(" -vs- ", "..vs..", colnames(dfnona))

## TD_comp_CE: time-discretized compared to CompetingExponentials

## Based on the figures, MCCBN_td is not as good a representative of the
## "_td" cases. Its average performance, when MHN_th and CBN_td do a great
## job, is generally poorer.

dfnona[, `:=`(CBN_MHN = (1/3) * (CBN..vs..MHN + CBN..vs..MCCBN + MCCBN..vs..MHN),
              OT_CAPRESE   = CAPRESE..vs..OT,
              ## TD =   (1/3) * (CBN_td..vs..MHN_td + CBN_td..vs..MCCBN_td + MCCBN_td..vs..MHN_td),
              TD =   CBN_td..vs..MHN_td,              
              ## TD_comp_CE = (1/3) * (CBN..vs..CBN_td + MCCBN..vs..MCCBN_td + MHN..vs..MHN_td),
              TD_comp_CE = (1/2) * (CBN..vs..CBN_td + MHN..vs..MHN_td),              
              CBN_comp_OT = (1/6) * (CBN..vs..OT + MCCBN..vs..OT + MHN..vs..OT +
                                     CAPRESE..vs..CBN + CAPRESE..vs..MCCBN + CAPRESE..vs..MHN),
              CE = CBN..vs..MHN,
              MHN_CE_TD = MHN..vs..MHN_td,
              CBN_CE_TD = CBN..vs..CBN_td
              ) ]

nrow(dfnona)
summary(dfnona)

## ## remove the 16 missing. Nope, don't
## dfnona_clean <- dfnona[!is.na(CBN_MHN), ]

dfnona_clean <- dfnona

nrow(dfnona_clean)
summary(dfnona_clean)

## Keep only variables we will need

dfnona_clean[, `:=`( typeLandscape = typeLandscape_f)]
levels(dfnona_clean$typeLandscape_f) <- c("Rep", "LM", "RMF")

with(dfnona_clean, table(typeLandscape_f, typeLandscape))



dfnona_clean <- dfnona_clean[, .(id, sourceGenotype, replicate,
                                 size_split, detect, 
                                 min_js, propLocalMax,
                                 sourceGenotype_nMut, numGenes,
                                 sampledProp,
                                 CE,
                                 MHN_CE_TD,
                                 CBN_CE_TD,
                                 CBN_MHN, OT_CAPRESE, TD,
                                 TD_comp_CE, CBN_comp_OT,
                                 typeLandscape_f
                                 )]


method_comp_fit <- dfnona_clean

save(file = "method_comp_fit.RData", method_comp_fit,
     compress = FALSE)

## All of them being similar is captured by, well, all of them being small

## 0.06
cor.test(method_comp_fit$CBN_MHN, method_comp_fit$OT_CAPRESE, method = "spearman")

## 0.7
cor.test(method_comp_fit$CBN_MHN, method_comp_fit$TD, method = "spearman")

## -0.05
cor.test(method_comp_fit$OT_CAPRESE, method_comp_fit$TD, method = "spearman")

## -0.15
cor.test(method_comp_fit$TD_comp_CE, method_comp_fit$TD, method = "spearman")

## 0.08
cor.test(method_comp_fit$TD_comp_CE, method_comp_fit$CBN_MHN, method = "spearman")

## 0.44
cor.test(method_comp_fit$TD_comp_CE, method_comp_fit$OT_CAPRESE, method = "spearman")

## 0.17
cor.test(method_comp_fit$TD_comp_CE, method_comp_fit$CBN_comp_OT, method = "spearman")


## 0.73
cor.test(method_comp_fit$CBN_MHN, method_comp_fit$CBN_comp_OT, method = "spearman")

## 0.26
cor.test(method_comp_fit$OT_CAPRESE, method_comp_fit$CBN_comp_OT, method = "spearman")




## In case we wanted all vars, which we dont
## dfnona_all_vars <- dfnona
## dfnona_all_vars <- dfnona_all_vars[!is.na(TD), ]

## colsToDelete <- colnames(dfnona_all_vars)[c(8:10, 12:20, 22:23)]

## dfnona_all_vars[, (colsToDelete) := NULL]

## method_comp_rf_fit <- dfnona_all_vars

## save(file = "method_comp_rf_fit.RData", method_comp_rf_fit,
##      compress = FALSE)










## ###############
## #####   Older stuff, pasting other info

## load("data_for_weighted_glmertree_plots.RData")

## nrow(data_for_weighted_glmertree_plots)
## nrow(df)

## df2 <- df
## max_g_7 <- 6
## max_g_10 <- 9

## df2 <- df2[
##     ((numGenes == 7) & (sourceGenotype_nMut <= max_g_7))
##     |
##     ((numGenes == 10) & (sourceGenotype_nMut <= max_g_10)),
##     ]

## nrow(df2)
## any(is.na(df2$min_js))

## df2 <- df2[!is.na(df2$min_js), ]
## nrow(df2)

## ## Make sure same names
## setnames(data_for_weighted_glmertree_plots,
##          old = c("sample_size"),
##          new = c("size_split"))

## samecols <- intersect(colnames(df2), colnames(data_for_weighted_glmertree_plots))

## df3 <- merge.data.table(df2, data_for_weighted_glmertree_plots,
##                         by = samecols,
##                         all.x = TRUE, all.y = TRUE)

## dim(df3)
## dim(df2)


## sum(is.na(df3$"CBN -vs- MHN"))
## sum(is.na(df3$"CBN -vs- MHN"))/nrow(df3) ## 34%

## somenas <- which(is.na(df3$"CBN -vs- MHN"))[c(1, 15, 1000, 200000, 815000, 815398)]

## df3[somenas, ]
## ## with(df, cor.test("CBN -vs- MCCBN", "CBN -vs- MHN"))

## df3somenas <- df3[somenas, ]

## save(file = "df3_somenas.RData", df3somenas)

## summary(df3[, 25:52])
