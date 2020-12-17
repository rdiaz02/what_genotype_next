##### As it says, this obtains a table for genotypes' weights-
## Launch as
## nohup R --vanilla --slave -f makeTableWeigths.R &> makeTableWeights.Rout &
rm(list = ls())
source("oncoFunctions.R")

library(codetools)
checkUsageEnv(env = .GlobalEnv)

## library(dplyr)
## library(parallel)
## library(OncoSimulR)

## Change paths as needed

dirsSims <- c("/Disk2/ramon/No-Backuppc/selected-simulations-A",
             "/Disk2/ramon/No-Backuppc/selected-simulations-B")
dirsSampled <- "/home/ramon/No-Backuppc/cpm-analyzed-all-selected"

## ## For playing
## dirsSampled <- "~/tmp"
## dirsSims <- "~/tmp"


## filesSim give the "true probabilities" of genotypes
## filesSample4000 are used to obtain the probabilities of genotypes under
##  the sampling regime (4000 * 5 = 20000, and that is where we get
##  the probabilities from)

filesSim <- list.files(dirsSims, pattern="simobj", full.names = TRUE)

filesSampled <- list.files(dirsSampled,
                           pattern = glob2rx("ANALYZED__ID*.rds"),
                           full.names = TRUE)
filesSampled4000 <- grep("_size_split_4000\\.rds", filesSampled,
                         value = TRUE)

stopifnot(length(filesSim) == 1260)
stopifnot(length(filesSampled4000) == (1260 * 3))


date()
## Warm up compiler and check
null <- dplyr::bind_rows(
                   lapply(filesSampled4000[1:5],
                          freqs_from_sampling_full))
null

null <- dplyr::bind_rows(
                   lapply(filesSim[1:5],
                          true_Freqs_full_2))

null
rm(null)
gc()

#### Sampled is much faster
## Do for real
pboptions(type = "txt")
cat("\n Doing sampled \n")
outFreqsSampled <- pbmclapply(filesSampled4000,
                            freqs_from_sampling_full,
                            mc.cores = detectCores())
save(file = "outFreqsSampled.RData", outFreqsSampled,
     compress = FALSE)
date()

outFreqsSampled <- dplyr::bind_rows(outFreqsSampled)
save(file = "outFreqsSampled.RData", outFreqsSampled,
     compress = FALSE)
date()



#### Simulations

date()
## Do for real
cat("\n Doing simulations \n")

## This is very slow. I think it is I/O bound (load average relative to
## %CPU, R processes waiting, etc, and files leave in an HDD, not the
## SSD). If rerunning, it could make sense to run this with lapply
## but parallelize inside true_Freqs. Place also timing sentinels
## in true_Freq_nulls.
## I tried it, with little improvement with the _2 versions
## in limited testing.
outFreqsTrue <- pbmclapply(filesSim,
                         true_Freqs_full_2,
                         mc.cores = 10) ## detectCores())
save(file = "outFreqsTrue.RData", outFreqsTrue,
     compress = FALSE)

## I was getting out of RAM errors
any_error <- unlist(lapply(outFreqsTrue, function(x) inherits(x, "try-error")))
if(sum(any_error)) {
    print(outFreqsTrue[any_error])
    sum(any_error)
    }


date()
outFreqsTrue <- dplyr::bind_rows(outFreqsTrue)
save(file = "outFreqsTrue.RData", outFreqsTrue,
     compress = FALSE)
gc()

date()



## BEWARE here!!
## Had we done this

## weightsAll <- dplyr::full_join(outFreqsSampled,
##                                outFreqsTrue)

## in weightsAll, some genotypes, that appear in the "true"
## never appeared in one or more of the three detect. But they
## could been in a POM. Fix that. 


## How can we fix that? Easy. Create a dummies and replicate.
## The True is the same in all sampling regimes
## and we ensure those genotypes always present as they should
outFreqsTrueL <- outFreqsTrueS <- outFreqsTrueU <- outFreqsTrue

outFreqsTrueL$detect <- "large"
outFreqsTrueS$detect <- "small"
outFreqsTrueU$detect <- "unif"

outFreqsTrue3 <- rbind(outFreqsTrueL,
                       outFreqsTrueS,
                       outFreqsTrueU)

nrow(outFreqsSampled)
nrow(outFreqsTrue)
nrow(outFreqsTrue3)


weightsAll <- dplyr::full_join(outFreqsSampled,
                               outFreqsTrue3)

nrow(weightsAll)

## Remember: NAs mean some genotypes not present
## We expect no missing values in TrueFreq and TrueProp
## except for the genotypes with the splits. But none here
## So no missing there
summary(weightsAll)

save(file = "weightsAll.RData", weightsAll,
     compress = FALSE)


date()







