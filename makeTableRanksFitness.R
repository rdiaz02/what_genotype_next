##### As it says, this obtains a table for genotypes' weights-
## Launch as
## nohup R --vanilla --slave -f input.R &> input.Rout &

## This is so fast (< 5 seconds in Draco) that can run interactively

rm(list = ls())
source("oncoFunctions.R")

## Change paths as needed

dirsSims <- c("/Disk2/ramon/No-Backuppc/selected-simulations-A",
             "/Disk2/ramon/No-Backuppc/selected-simulations-B")

## ## For playing
## dirsSampled <- "~/tmp"
## dirsSims <- "~/tmp"


## filesSim give the "true probabilities" of genotypes
## filesSample4000 are used to obtain the probabilities of genotypes under
##  the sampling regime (4000 * 5 = 20000, and that is where we get
##  the probabilities from)

filesSim <- list.files(dirsSims, pattern = glob2rx("landscape_ID_*.rds"),
                       full.names = TRUE)
length(filesSim)




date()
## Warm up compiler and check

null <- dplyr::bind_rows(
                   lapply(filesSim[1:5],
                          fitness_rank_genotypes)
               )

null
rm(null)
gc()

## Do for real
pboptions(type = "txt")
cat("\n Doing sampled \n")

outRanksFitness <- pbmclapply(filesSim,
                              fitness_rank_genotypes,
                              mc.cores = detectCores())

save(file = "outRanksFitness.RData", outRanksFitness,
     compress = FALSE)
date()

outRanksFitness <- dplyr::bind_rows(outRanksFitness)

save(file = "outRanksFitness.RData", outRanksFitness,
     compress = FALSE)
date()

