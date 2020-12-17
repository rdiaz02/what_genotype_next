## Copyright 2020 Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


## This is used to obtain Schill et al.,'s MHN output from the simulated
## data.

date()
machine <- Sys.info()["nodename"]
print(machine)


library(parallel)
library(stringr)

######################################################################

##   Load common code for analysis

######################################################################

source("schill-trans-mat.R")
source("pre-process.R")


## From run-all-methods.R
add_pseudosamples <- function(x, n00 = "auto3") {
    if(n00 == "auto") {
        if(nrow(x) <= 500) {
            n00 <- round(nrow(x) * 0.10)
        } else {
            n00 <- round(nrow(x) * 0.05)
        }
    } else if(n00 == "auto2") {
        ## add only if max. freq. of any gene is > 95%
        fmax <- max(colSums(x))/nrow(x)
        if(fmax > 0.95)
            n00 <- round(nrow(x) * 0.05)
        else
            n00 <- 0
    } else if(n00 == "auto3") {
        ## add only if any gene is 100%
        ## add just 1
        fmax <- max(colSums(x))/nrow(x)
        if(fmax == 1)
            n00 <- 1
        else
            n00 <- 0
    }    
    return(rbind(x,
                 matrix(0L, nrow = n00, ncol = ncol(x))
                 ))
}




## From function all_methods in file run-all-methods.R in paper with Claudia

## by default, only write file
mhn_on_analyzed <- function(thefile, return_object = FALSE) {
    x <- readRDS(thefile)
    outn <- x[1:12]
    ## muvar <- str_split_fixed(
    ##     str_split_fixed(thefile, "_muvar_", n = 2)[1, 2], "_", n = 2)[1, 1]
    ##string <- unlist(x[c(7, 2, 3, 9, 12)])
    ## I added another level of nesting somehow
    string <- as.data.frame(x[c(7, 8, 2, 3, 9, 12, 1)], stringsAsFactors = FALSE)
    string <- paste(paste(names(string), string, sep = "_"), collapse = "__")
    num_splits <- length(x$allM)

    cat("\n    Doing string = ", string, "\n")


    out_w <- lapply(1:num_splits,
                    function(n) {
                        list(MHN = mhn_on_split(x$allM[[n]]))
                    })

    theout <- c(outn,
                out_paths_genots = list(out_w)
                )

    saveRDS(theout, file = paste0("MHN_out_", string, ".rds"))
    cat("\n       done string = ", string, "\n")
    if(return_object) return(theout)
}


## An allM component from an ANALYZED__ID -> Schill's MHN output
##       with additional error handling
mhn_on_split <- function(y, string = NULL) {
    x <- y$input_data_pseudosamples
    x_tmp <- pre_process(x, remove.constant = FALSE, min.freq = 0)


    if(ncol(x_tmp) >= 2) {
        time_mhn <- system.time(MHN <- try(do_MHN(x_tmp, lambda = 1/nrow(x_tmp))))["elapsed"]
    } else {
        time_mhn <- NA
        cat("\n         Cannot run:  ncol < 2   \n")
    }

    ## Names of components and how they match those from CBN, etc
    ## theta: there is no equivalent
    ## transitionRateMatrix: weighted_fgraph (but only for CBN and MCCBN, not OT)
    ## transitionMatrixCompExp: trans_mat_genots
    ## transitionMatrixTimeDiscretized: transitionMatrixTimeDiscretized
    ##                          from new code for CBN and MCCBN
    ##                          see time-discretized-CBN.R
    
    ## Add appropriate error indicators
    if(ncol(x_tmp) < 2) {
        MHN <- list(theta = "ERROR_CPM_ANALYSIS",
                    transitionRateMatrix = "ERROR_CPM_ANALYSIS",
                    transitionMatrixTimeDiscretized = "ERROR_CPM_ANALYSIS",
                    transitionMatrixCompExp = "ERROR_CPM_ANALYSIS",
                    time_mhn = "NA",
                    likely_error = "ncol_x")
    } else if(inherits(MHN, "try-error")) {
        MHN <- list(theta = "ERROR_CPM_ANALYSIS",
                    transitionRateMatrix = "ERROR_CPM_ANALYSIS",
                    transitionMatrixTimeDiscretized = "ERROR_CPM_ANALYSIS",
                    transitionMatrixCompExp = "ERROR_CPM_ANALYSIS",
                    time_mhn = time_mhn,
                    likely_error = "Error_in_run")
        cat("\n         Error_in_run:   \n")
    } else {
        MHN <- c(MHN,
                 time_mhn = time_mhn,
                 likely_error = "No_error")
    }
    cat("\n MHN time: ", time_mhn, "\n")

    return(MHN)
    ## out <- list(MHN = MHN,
    ##             time_mhn = time_mhn
    ##             )
}





library(codetools)
checkUsageEnv(env = .GlobalEnv)


######################################################################

##   Select the files to use

######################################################################


## ff <- dir(
##     path = "../selected-simulations-splitted-data",
##     pattern = glob2rx("^n_splitted_data__ID*"), full.names = TRUE)

## ## do 7, then 10
## ffdo7 <- grep("_ng_7_fl_ng_7_", ff, value = TRUE)
## ffdo10 <- grep("_ng_10_fl_ng_10_", ff, value = TRUE)

## rm(ff)

## ff <- c(ffdo7, ffdo10)

## ## for checks
## lff <- length(ff)
## lffdo7 <- length(ffdo7)
## lffdo10 <- length(ffdo10)
## stopifnot(lff == (2 * lffdo7))
## stopifnot(lffdo7 == lffdo10)
## stopifnot(length(ff) == 3780)

## Numbers are
## 3780 total, 1890 each ngenes, with those begin 630 * 3 (types of sampling)
## and 630 = 35 * 3 * 2 * 3 (replicates, type landscape, muvar, initsize)



## The ANALYZED files do have the pseudodata, etc.
aa <- dir(
    path = "../cpm-analyzed-all-selected",
    pattern = glob2rx("ANALYZED__ID*.rds"), full.names = TRUE)

## 11340 = 1260 * 3 * 3
## 11340 = 3780 * 3
stopifnot(length(aa) == 11340)


## a single example run
cucu <- mhn_on_analyzed(aa[[1]])
names(cucu)
length(cucu$out_paths_genots)
names(cucu$out_paths_genots[[1]])
names(cucu$out_paths_genots[[1]][["MHN"]])

mclapply(aa,
         mhn_on_analyzed,
         mc.cores = detectCores(),
         mc.preschedule = FALSE,
         mc.allow.recursive = FALSE)

date()
gc()











## ## From function all_methods in file run-all-methods.R in paper with Claudia
## mhn_method <- function(x,
##                        n00 = "auto3", 
##                        min.freq = 0) {


##     x000 <- x
##     x <- add_pseudosamples(x, n00 = n00)
##     ## remove.constant makes no difference IFF we add pseudosamples, as
##     ## there can be no constant column when we add pseudosamples
##     x_tmp <- pre_process(x, remove.constant = FALSE, min.freq = min.freq)


##     if(ncol(x_tmp) >= 2) {
##         time_mhn <- system.time(MHN <- try(do_MHN(x_tmp, lambda = 0.01)))["elapsed"]
##     } else {
##         time_mhn <- NA
##         MHN <- NA
##         cat("\n         Cannot run:  ncol < 2   \n")
##     }
##     cat("\n MHN time: ", time_mhn, "\n")

##     out <- list(MHN = MHN,
##                 time_mhn = time_mhn
##                 )

##     ## FIXME prepare the appropriate structure here
   
    
## }




## ## From run_on_list_splitted, in new-run-CPMs.R (paper with Claudia)
## ## Main change: insteadll of all_methods, calling mhn_method

## ## name of file (rds) with splitted data, size_split -> output of simulations
## ##   based on "split_run_methods" but here data are already splitted
## run_MHN_on_list_splitted <- function(x, save = TRUE) {
##     inittime <- Sys.time()
##     cat("\n       doing file ", x)
##     u <- readRDS(x)

##     ## Extract things for naming
##     pos <- c(6, 27, 8, 17, 15, 23, 25, 21, 37, 39)
##     pieces <- str_split_fixed(x, "_", n = 40)
##     v <- pieces[pos]
##     names(v) <- pieces[pos - 1]
    
##     if(v["RMF"] == "TRUE") {
##         typeLandscape <- "RMF" 
##     } else if( (v["RMF"] == "FALSE" ) && (v["nrnds"] == "0") ) {
##         typeLandscape <- "Represent."
##     } else if( (v["RMF"] == "FALSE" ) && (v["nrnds"] == "50") ) {
##         typeLandscape <- "Local_peaks"
##     } else {
##         stop("unknown option")
##     }
##     fname <- strsplit(x, "_betaA_")[[1]][1]
##     name <- paste0("ID_", v["ID"],
##                    "_rnst_", v["rnst"],
##                    "_typeLandscape_", typeLandscape,
##                    "_ng_", v["ng"],
##                    "_muvar_", ifelse(v["var"] == "NA", FALSE, TRUE),
##                    "_initS_", v["initS"],
##                    "_beta_A_", v["betaA"],
##                    "_beta_B_", v["betaB"])

##     ## This looks much more cumbersome than the much simpler looping of,
##     ## say, new-new-diversity-samples.R. Here I inherited the way of doing
##     ## the splitting and analysis from older scripts. I no longer sample,
##     ## but I keep the idea of cycling over the replicates given a split
##     ## size.
    
##     size_splits <- c(50, 200, 4000)

##     outs <- lapply(size_splits, function(size_split) {
##         inittime_i <- Sys.time()
##         sspl <- which(unlist(lapply(u, function(v) v$size_split == size_split)))
##         uspl <- u[sspl]
##         ## allM <- lapply(uspl, function(d) am(d$data))
##         ## allM <- lapply(uspl, function(d) all_methods(d$data))
##         allM <- lapply(uspl, function(d) mhn_method(d$data))

##         ## FIXME: these need not be inside this function! It should have
##         ## been out of the function. This is common to all.
##         if( (v["betaA"] == "1") && (v["betaB"] == "1") ) {
##             detect <- "unif"
##         } else if ( (v["betaA"] == "3") && (v["betaB"] == "5") ) {
##             detect <- "small"
##         } else if ( (v["betaA"] == "5") && (v["betaB"] == "3") ) {
##             detect <- "large"
##         } else {
##             stop("Invalid beta options")
##         }
##         out <- list(
##             detect = detect,
##             ## N = 20000,
##             ngenes = v["ng"],
##             initSize = v["initS"],
##             mu = 1e-5,
##             ## next will not work. so this is a useless variable
##             ## FIXME: so why did I not do something useful?
##             mumaxprodvar = ifelse(v["var"], 25, 1),
##             nrounds = v["nrnds"],
##             ## muvar = v["muvar"],
##             ID = v["ID"],
##             rnst = v["rnst"],
##             typeLandscape = typeLandscape,
##             is_dag_a_tree = v["tree"],
##             fname = fname,
##             size_split = size_split,
##             ## beta_a = v["betaA"],
##             ## beta_b = v["betaB"],
##             allM = allM            
##         )
##         ## We save the data in the rds files, but not for the output
##         ## Saving in the rds is also redundant but anyway
##         if(save) {
##             nameout <- paste0(name, "_size_split_", size_split)
##             saveRDS(out, file = paste0("ANALYZED__", nameout, ".rds"))
##         }
##         cat("\n              doing CPM size_split ", size_split, " took(ssp) ",
##             difftime(Sys.time(), inittime_i, units = "min"),
##             " mins")
##         for(ll in 1:length(out$allM)) {
##             out$allM[[ll]]$input_data <- NULL
##             out$allM[[ll]]$input_data_pseudosamples <- NULL
##             ## these might not exist, if no model was fitted
##             try(out$allM[[ll]]$CAPRESE$TRONCO_model$genotypes <- NULL, silent = TRUE)
##             try(out$allM[[ll]]$CAPRI_BIC$TRONCO_model$genotypes <- NULL, silent = TRUE)
##             try(out$allM[[ll]]$CAPRI_AIC$TRONCO_model$genotypes <- NULL, silent = TRUE)
##         }
##         return(out)
##     })

##     cat("\n      done file ", x, ".  Took ",
##         difftime(Sys.time(), inittime, units = "min"), " mins \n")
##     ## single, three element list, to be analyzed with other scripts
##     return(do.call(c, list(outs)))
## }
