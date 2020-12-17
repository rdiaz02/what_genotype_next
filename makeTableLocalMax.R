##### As it says, this obtains a table for genotypes' frequency as local max
## Launch as
## nohup R --vanilla --slave -f input.R &> input.Rout &

rm(list = ls())
source("oncoFunctions.R")

library(codetools)
checkUsageEnv(env = .GlobalEnv)

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

filesSim <- list.files(dirsSims, pattern = "simobj",
                       full.names = TRUE)
stopifnot(length(filesSim) == 1260)




date()
## Warm up compiler and check

## system.time(
## null <- dplyr::bind_rows(
##                    lapply(filesSim[c(1, 15, 2, 3)],
##                           local_max)
##                )
## )

## system.time(
## null <- dplyr::bind_rows(
##                    mclapply(filesSim[c(1, 15, 2, 3)],
##                             local_max,
##                             mc.cores = 4)
##                )
## )


system.time(
null <- dplyr::bind_rows(
                   lapply(filesSim[c(1, 15)],
                          local_max_np)
               )
)

## check compiled
local_max_np

## system.time(
## null <- dplyr::bind_rows(
##                    mclapply(filesSim[c(1, 15, 2, 3)],
##                             local_max_np,
##                             mc.cores = 4)
##                )
## )



null
rm(null)
gc()

pboptions(type = "txt")
outLocalMax <- pbmclapply(filesSim,
                          local_max_np,
                          mc.cores = 10) ## detectCores())

cat("\n Save before bind_rows \n")
save(file = "outLocalMax.RData", outLocalMax,
     compress = FALSE)
date()

stopifnot(length(outLocalMax) == 1260)

## I was getting out of RAM errors
any_error <- unlist(lapply(outLocalMax, function(x) inherits(x, "try-error")))
if(sum(any_error)) {
    print(outLocalMax[any_error])
    sum(any_error)
    }

outLocalMax <- dplyr::bind_rows(outLocalMax)

cat("\n Save after bind_rows \n")
save(file = "outLocalMax.RData", outLocalMax,
     compress = FALSE)
date()






