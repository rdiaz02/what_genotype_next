##### As it says, this obtains a table for genotypes' weights-
## Launch as
## nohup R --vanilla --slave -f input_name.R &> input_name.Rout &
rm(list = ls())
source("oncoFunctions.R")

## library(dplyr)
## library(parallel)
## library(OncoSimulR)

## Change paths as needed

dirsSampled <- "/home/ramon/No-Backuppc/cpm-analyzed-all-selected"

## ## For playing
## dirsSampled <- "~/tmp"
## dirsSims <- "~/tmp"


filesSampled <- list.files(dirsSampled,
                           pattern = glob2rx("ANALYZED__ID*.rds"),
                           full.names = TRUE)

length(filesSampled)
## When we parallelized, do not have patterns of all easy early, etc.
filesSampled <- sample(filesSampled)

date()
## Warm up compiler and check
null <- dplyr::bind_rows(
                   lapply(filesSampled[1:5],
                          obs_freqs_in_sample))
null

rm(null)
gc()


pboptions(type = "txt")
cat("\n Doing sampled \n")
system.time(
outObsFreqs <- pbmclapply(filesSampled,
                          obs_freqs_in_sample,
                          ## No: just reorder input
                          ## mc.preschedule = FALSE,
                          mc.cores = detectCores())
)

save(file = "outObsFreqs.RData", outObsFreqs,
     compress = FALSE)
date()

outObsFreqs <- dplyr::bind_rows(outObsFreqs)
save(file = "outObsFreqs.RData", outObsFreqs,
     compress = FALSE)
date()
