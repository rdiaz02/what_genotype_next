## Merge output from method differences and wide_g for glmertree fits.

## This takes a while (~ 5 minutes)
dd <- readRDS("comps-fixNulls.rds")

## steps
## - clean up the extra stuff in the columns
## - do a merge

library(data.table)
setDTthreads(threads = 0)

ddt <- as.data.table(dd)

load("wide_g.RData")

## source genotype???

any(unique(wide_g$sourceGenotype) == "WT") ## FALSE
any(unique(wide_g$sourceGenotype) == "root") ## TRUE

uddtsg <- unique(ddt$sourceGenotype)
any(uddtsg == "root") ## FALSE
any(uddtsg == "WT") ## FALSE
any(uddtsg == " ") ## FALSE
any(uddtsg == "") ## TRUE
## We assume the last are the root. We will see.




## ddt[, `:=`(size_split_clean = gsub("size_split_", "", size_split),
##            detect_clean = gsub("detect_", "", detect),
##            sourceGenotype_clean = sourceGenotype)]
## ddt$sourceGenotype_clean[ddt$sourceGenotype_clean == ""] <- "root"

ddt[, `:=`(size_split = gsub("size_split_", "", size_split),
           detect = gsub("detect_", "", detect))
    ]

ddt$detect[ddt$detect == "unif"] <- "uniform"
unique(ddt$detect)

ddt$sourceGenotype[ddt$sourceGenotype == ""] <- "root"


## Same class. Reset later to factors
wide_g[, `:=`(
    id = as.character(id),
    detect = as.character(detect),
    sourceGenotype = as.character(sourceGenotype)
) ]

## Leave as numeric.
ddt[, `:=`(size_split = as.numeric(size_split))]


wide_g_method_comp <- merge(wide_g, ddt,
                              by = c("id", "size_split", "detect", "sourceGenotype", "replicate"),
                              all.x = TRUE, all.y = FALSE)

dim(wide_g)
dim(wide_g_method_comp)

head(wide_g_method_comp)

stopifnot(nrow(wide_g) == nrow(wide_g_method_comp))


## Reset id, detect, sourceGenoytpe to factors
## size_split will later be converted in the fitting routine
wide_g_method_comp[, `:=`(id = as.factor(id),
                          detect = as.factor(detect),
                          sourceGenotype = as.factor(sourceGenotype))]


save(file = "wide_g_method_comp.RData", wide_g_method_comp, compress = FALSE)


## quick check
