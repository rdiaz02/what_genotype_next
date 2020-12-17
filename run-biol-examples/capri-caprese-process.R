## Copyright 2016, 2017, 2018 Ramon Diaz-Uriarte

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


## rm(list = ls())
library(TRONCO)
library(help = TRONCO)
library(parallel)

## Recall caprese and capri crash if only one column of data.

## f_caprese <- function(x, cores.ratio = 1, silent = FALSE,
##                       nboot = 100) {
##     ## Like above, but says what file it is processing and is verbose
##     if(!silent)
##         cat("\n .... Starting with data ", x$name, "\n")
##     datax <- x$out$popSample
##     datax <- pre_process(datax, remove.constant = TRUE)
##     if(ncol(datax) < 2) {
##         return(list(name = x$name, params = x$params,
##                     nocols = TRUE,
##                     time = NA,
##                     res = NA))
##     } else {
##         datax <- tronco_common(datax)
##         time <- system.time({
##             otout <- try(caprese_proc(genotData = datax, nboot = nboot,
##                                       silent = silent,
##                                       cores.ratio = cores.ratio))
##         })
##         return(list(name = x$name, params = x$params,
##                     nocols = FALSE,
##                     time = time,
##                     res = otout))
##     }
## }

## f_capri <- function(x, silent = FALSE, nboot = 100,
##                     nbootWilcox = 100,
##                     regularization = "bic",
##                     cores.ratio = 1) {
##     ## nboot is for the bootstrap per se
##     ## nbootWilcox is the "nboot" argument in tronco.capri.
##     if(!silent)
##         cat("\n .... Starting with data ", x$name, "\n")
##     datax <- x$out$popSample
##     datax <- pre_process(datax, remove.constant = TRUE)
##     if(ncol(datax) < 2) {
##         return(list(name = x$name, params = x$params,
##                     nocols = TRUE,
##                     time = NA,
##                     res = NA))
##     } else {
##         datax <- tronco_common(datax)
##         time <- system.time({
##             otout <- try(capri_proc(genotData = datax, nboot = nboot,
##                                     nbootWilcox = nbootWilcox,
##                                     regularization = regularization,
##                                     command = "hc", pvalue = 0.05,
##                                     cores.ratio = cores.ratio))
##         })
##         return(list(name = x$name, params = x$params,
##                     nocols = FALSE,
##                     time = time,
##                     res = otout))
##     }
## }






tronco_common <- function(x) {
    if(is.null(rownames(x)))
        rownames(x) <- 1:(nrow(x))
    return(TRONCO::import.genotypes(data.frame(x)))
}




caprese_capri_common <- function(original.adj.mat,
                                 boot.adj.mat, bt, pr, tp, hg,
                                 genotData) {
    message("          Starting caprese_capri_common ", date())
    ## Why do we do all of this, instead of just using the bt and other
    ## confidences? Because that would mean we loose all information about
    ## how often there is a Root -> some node in the the boostrap
    
    ## We spend most of the complex code just dealing with their naming
    ## stuff.
    ## Paranoid checks
    ## Note that this need not be true: stopifnot(identical(boot.adj.mat[-1, -1], bt))
    ## as first are in terms of all edges, and second only in terms of edges in
    ## reconstructed tree; so boot.adj.mat can have non-zero entries for relationships
    ## that are not in the reconstructed tree.
    ## But this should be true
    stopifnot(all(boot.adj.mat[-1, -1] >=  bt))
    ## And the next too, of course
    stopifnot(identical(rownames(boot.adj.mat[-1, -1]),
                        rownames(bt)))
    stopifnot(identical(colnames(boot.adj.mat[-1, -1]),
                        colnames(bt)))
    stopifnot(identical(rownames(boot.adj.mat),
                        colnames(boot.adj.mat)))
    

    ## Complicated because of what they do with names and indices. Since
    ## such a mess, I check below. I do it two ways.

    ## First (and the one that leads to out) is:
    ## - 1. start with large adjacency matrix (includes Root, or None):
    ##    create a data frame with pairs of Rows/Columns of non-zero entries.
    ## - 2. remove from that entries not in original tree, but preserve root:
    ##     - create a matrix of pairs Row/Column based on original matrix
    ##       (this has no Root)
    ##     - add a "Root -> Something" node to the pairs of rows and columns
    ##         to keep
    ##     - use this to keep from 1. only entries in original tree.
    ## - 3. Add tp. hg. pr.
    ##     - Use the Rows/Cols from 2. to index into the tp, hg, pr matrices.
    ##     - As "None" is not a valid row name, make it "NA" in the index.

    ## At the end, we check the results in out are the same as those obtained
    ## by accessing in a very different way.

    
    ## Get the indices from the bootstrap results, including their None
    ## Some of this is redundant, but helps see what happens because of
    ## the mess of their naming scheme, jumping to "G?", etc.
    indices <- which(boot.adj.mat > 0, arr.ind = TRUE)
    bootNamesResults <- data.frame(Row = indices[, 1],
                                   Col = indices[, 2],
                                   RowName = rownames(boot.adj.mat)[indices[, 1]],
                                   ColName = colnames(boot.adj.mat)[indices[, 2]],
                                   bootfreq = boot.adj.mat[indices],
                                   stringsAsFactors = FALSE)
    
    ## But some of these might not correspond to edges in the original
    ## tree. Find those in the original, and remove the remaining.  Yes, a
    ## lot of this would be easier just adding a 1 but we cannot bet on
    ## it?
    ## We need to do extra work because the original adjacency matrix
    ## has removed the root node.
    original.root.edge <- which(colSums(original.adj.mat) == 0)
    original.adj.mat.ind <- which(original.adj.mat == 1, arr.ind = TRUE)
    original.adj.mat.row.name <-
        rownames(original.adj.mat)[original.adj.mat.ind[, 1]]
    original.adj.mat.col.name <-
        colnames(original.adj.mat)[original.adj.mat.ind[, 2]]
    fromRoot <- colnames(original.adj.mat)[original.root.edge]
    originalAMRowsCols <- data.frame(RowName = c(rep("None", length(fromRoot)),
                                         original.adj.mat.row.name),
                             ColName = c(fromRoot, original.adj.mat.col.name),
                             stringsAsFactors = FALSE)

    ## Keep only those in original model
    bootNamesResults <- dplyr::inner_join(bootNamesResults, originalAMRowsCols,
                                       by = c("RowName", "ColName")) 

    ## We also want access to tp, pr, etc.
    bootNamesResults$RowNameNoRoot <- bootNamesResults$RowName
    bootNamesResults$RowNameNoRoot[bootNamesResults$RowName == "None"] <- NA
    ## There should be none of these?
    bootNamesResults$ColNameNoRoot <- bootNamesResults$ColName
    stopifnot(!any(bootNamesResults$ColName == "None"))

    ## Many indices are NA, so they will return an NA value.
    pr.tp.pos.get <- cbind(bootNamesResults$RowNameNoRoot,
                           bootNamesResults$ColNameNoRoot)

    ## Create lookupNames object to return in terms of original names
    lookupNames <- genotData$annotations[, "event"]
    ## Check it is doing what we think: they just add "G" to the events,
    ## numbered 1 to ... All next three must be true
    nln <- sub("^G", "", names(lookupNames))
    if( ! ( all(as.numeric(nln) == seq_along(lookupNames))   &&
            all(nln == as.character(seq_along(lookupNames))) &&
            all(names(lookupNames) == paste0("G", seq_along(lookupNames)))) )
        stop("Unexpected names of genotype data annotations")
    lookupNames <- c(lookupNames, "None" = "Root")
    bootNamesResults$RowNameAsOriginal <- lookupNames[bootNamesResults$RowName]
    bootNamesResults$ColNameAsOriginal <- lookupNames[bootNamesResults$ColName]
    
    out <- data.frame(From = bootNamesResults$RowNameAsOriginal,
                      To = bootNamesResults$ColNameAsOriginal,
                      edge = paste(bootNamesResults$RowNameAsOriginal,
                                   bootNamesResults$ColNameAsOriginal,
                                   sep = " -> "),
                      edgeBootFreq = bootNamesResults$bootfreq,
                      pr = pr[pr.tp.pos.get],
                      tp = tp[pr.tp.pos.get],
                      hg = hg[pr.tp.pos.get],
                      stringsAsFactors = FALSE
                      )

    ## Check EVERYTHING but accessing elements in a different way.
    ## The row and column names in results
    mat.ind.res <- cbind(bootNamesResults[, "RowName"],
                         bootNamesResults[, "ColName"])
    ## Check the bootstrap. Set to zero any not in the original tree.
    ## Above we removed from the the results object, the bootNameResults.
    ## Here we set them to zero. We check that possible reordering in output
    ## does not screw match.
    not.in.original <- which(original.adj.mat == 0, arr.ind = TRUE)
    zeroed.boot.adj.mat <- boot.adj.mat
    zeroed.boot.adj.mat[not.in.original + 1] <- 0
    stopifnot(all(zeroed.boot.adj.mat[mat.ind.res] == out$edgeBootFreq))
    ## same
    stopifnot(all(boot.adj.mat[mat.ind.res] == bootNamesResults$bootfreq))
    ## and make sure no na in bootstrap results
    stopifnot(!any(is.na(out$edgeBootFreq)))
    
    ## Check pr, tp, hg.
    ## First, check consistent col and rownames
    stopifnot(identical(colnames(pr), rownames(pr)))
    stopifnot(identical(colnames(tp), rownames(tp)))
    stopifnot(identical(colnames(hg), rownames(hg)))
    stopifnot(identical(colnames(bt), rownames(bt)))    
    stopifnot(identical(colnames(tp), colnames(pr)))
    stopifnot(identical(colnames(tp), colnames(hg)))
    stopifnot(identical(colnames(tp), colnames(bt)))
    ## Now, create matrices of pr, tp, hg, bt but with a "None"
    pr.augmented <- cbind(NA, rbind(NA, pr))
    colnames(pr.augmented) <- rownames(pr.augmented) <- c("None", colnames(pr))
    tp.augmented <- cbind(NA, rbind(NA, tp))
    colnames(tp.augmented) <- rownames(tp.augmented) <- c("None", colnames(tp))
    hg.augmented <- cbind(NA, rbind(NA, hg))
    colnames(hg.augmented) <- rownames(hg.augmented) <- c("None", colnames(hg))
    ## Bootstrap now extending the bt object; entries from Root will be missing
    bt.augmented <- cbind(NA, rbind(NA, bt))
    colnames(bt.augmented) <- rownames(bt.augmented) <- c("None", colnames(bt))
    ## Above we accessed with a matrix with indices that had NA.
    ## Here we access with a matrix of indices that have no NAs, but the object,
    ## tp.augmented, etc, have NAs in any None -> Something node.
    ## We na.omit each separately just in case; we do not want to get it
    ## right because of NAs in different positions in each column
    stopifnot(all(na.omit(pr.augmented[mat.ind.res]) == na.omit(out$pr)))
    stopifnot(all(na.omit(tp.augmented[mat.ind.res]) == na.omit(out$tp)))
    stopifnot(all(na.omit(hg.augmented[mat.ind.res]) == na.omit(out$hg)))
    
    ## bt is messier since bt has NA for those from None;
    ## the na.omit only uses bt.augmented, since results have no NA for boot;
    ## as checked above
    btc <- na.omit(bt.augmented[mat.ind.res] == out$edgeBootFreq)
    stopifnot(all(btc))
    message("         Ending caprese_capri_common ", date())
    return(out)
}



caprese_capri_common_no_boot <- function(original.adj.mat,
                                         pr, tp, hg,
                                         genotData) {
    message("         Starting caprese_capri_common_no_boot ", date())
    ## Same as function above, but for no bootstrap.
    ## We use same logic and names
    original.root.edge <- which(colSums(original.adj.mat) == 0)
    original.adj.mat.ind <- which(original.adj.mat == 1, arr.ind = TRUE)
    original.adj.mat.row.name <-
        rownames(original.adj.mat)[original.adj.mat.ind[, 1]]
    original.adj.mat.col.name <-
        colnames(original.adj.mat)[original.adj.mat.ind[, 2]]
    fromRoot <- colnames(original.adj.mat)[original.root.edge]
    originalAMRowsCols <- data.frame(RowName = c(rep("None", length(fromRoot)),
                                                 original.adj.mat.row.name),
                                     ColName = c(fromRoot, original.adj.mat.col.name),
                                     stringsAsFactors = FALSE)
    ## This is just silly, but to keep same names
    bootNamesResults <- originalAMRowsCols

    ## We also want access to tp, pr, etc.
    bootNamesResults$RowNameNoRoot <- bootNamesResults$RowName
    bootNamesResults$RowNameNoRoot[bootNamesResults$RowName == "None"] <- NA
    ## There should be none of these?
    bootNamesResults$ColNameNoRoot <- bootNamesResults$ColName
    stopifnot(!any(bootNamesResults$ColName == "None"))

    ## Many indices are NA, so they will return an NA value.
    pr.tp.pos.get <- cbind(bootNamesResults$RowNameNoRoot,
                           bootNamesResults$ColNameNoRoot)
    
     ## Create lookupNames object to return in terms of original names
    lookupNames <- genotData$annotations[, "event"]
    ## Check it is doing what we think: they just add "G" to the events,
    ## numbered 1 to ... All next three must be true
    nln <- sub("^G", "", names(lookupNames))
    if( ! ( all(as.numeric(nln) == seq_along(lookupNames))   &&
            all(nln == as.character(seq_along(lookupNames))) &&
            all(names(lookupNames) == paste0("G", seq_along(lookupNames)))) )
        stop("Unexpected names of genotype data annotations")
    lookupNames <- c(lookupNames, "None" = "Root")
    bootNamesResults$RowNameAsOriginal <- lookupNames[bootNamesResults$RowName]
    bootNamesResults$ColNameAsOriginal <- lookupNames[bootNamesResults$ColName]
    
    out <- data.frame(From = bootNamesResults$RowNameAsOriginal,
                      To = bootNamesResults$ColNameAsOriginal,
                      edge = paste(bootNamesResults$RowNameAsOriginal,
                                   bootNamesResults$ColNameAsOriginal,
                                   sep = " -> "),
                      edgeBootFreq = NA,
                      pr = pr[pr.tp.pos.get],
                      tp = tp[pr.tp.pos.get],
                      hg = hg[pr.tp.pos.get],
                      stringsAsFactors = FALSE
                      )

    ## Check EVERYTHING but accessing elements in a different way.
    ## The row and column names in results
    mat.ind.res <- cbind(bootNamesResults[, "RowName"],
                         bootNamesResults[, "ColName"])
    ## Check the bootstrap. Set to zero any not in the original tree.
    ## Above we removed from the the results object, the bootNameResults.
    ## Here we set them to zero. We check that possible reordering in output
    ## does not screw match.
    ## not.in.original <- which(original.adj.mat == 0, arr.ind = TRUE)
    
    ## Check pr, tp, hg.
    ## First, check consistent col and rownames
    stopifnot(identical(colnames(pr), rownames(pr)))
    stopifnot(identical(colnames(tp), rownames(tp)))
    stopifnot(identical(colnames(hg), rownames(hg)))
    stopifnot(identical(colnames(tp), colnames(pr)))
    stopifnot(identical(colnames(tp), colnames(hg)))
    ## Now, create matrices of pr, tp, hg, bt but with a "None"
    pr.augmented <- cbind(NA, rbind(NA, pr))
    colnames(pr.augmented) <- rownames(pr.augmented) <- c("None", colnames(pr))
    tp.augmented <- cbind(NA, rbind(NA, tp))
    colnames(tp.augmented) <- rownames(tp.augmented) <- c("None", colnames(tp))
    hg.augmented <- cbind(NA, rbind(NA, hg))
    colnames(hg.augmented) <- rownames(hg.augmented) <- c("None", colnames(hg))
    ## Above we accessed with a matrix with indices that had NA.
    ## Here we access with a matrix of indices that have no NAs, but the object,
    ## tp.augmented, etc, have NAs in any None -> Something node.
    ## We na.omit each separately just in case; we do not want to get it
    ## right because of NAs in different positions in each column
    stopifnot(all(na.omit(pr.augmented[mat.ind.res]) == na.omit(out$pr)))
    stopifnot(all(na.omit(tp.augmented[mat.ind.res]) == na.omit(out$tp)))
    stopifnot(all(na.omit(hg.augmented[mat.ind.res]) == na.omit(out$hg)))
    
    message("         Ending caprese_capri_common_no_boot ", date())    
    return(out)
}

useless_extra <- function(model, boot) {
    ## as.selective.advantage.relations(model) ## already in out
    ## as.bootstrap.scores(boot) ## already in out
    ## Just in case, since all this is so terribly confusing
    if(is.null(boot)) {
        return(list(selective_advantage = as.selective.advantage.relations(model),
                    bootstrap_scores = NA))
    } else  {
        return(list(selective_advantage = as.selective.advantage.relations(model),
             bootstrap_scores = as.bootstrap.scores(boot)))
    }
}

estimated_probs <- function(model) {
    ## We provide a redundant thing, as I never know with CAPRI and
    ## CAPRESE...  Actually, this is all rather useless, as these are the
    ## observed, not the predicted. And the model, inside in
    ## model$capri_bic$probabilities$probabilities.fit$estimated.conditional.probs
    ## and similar, only contains NAs.
    mp <- as.marginal.probs(model)
    cp <- as.conditional.probs(model)
    jp <- as.joint.probs(model)
    mp1 <- data.frame(mp[[1]])
    cp1 <- data.frame(cp[[1]])
    jp1 <- jp[[1]]
    mp1$Gene <- model$annotations[rownames(mp1), "event"]
    cp1$Gene <- model$annotations[rownames(cp1), "event"]
    colnames(jp1) <- model$annotations[colnames(jp1), "event"]
    rownames(jp1) <- model$annotations[rownames(jp1), "event"]
    
    list(marginal = mp,
         conditional = cp,
         joint = jp,
         marginal_clean = mp1,
         conditional_clean = cp1,
         joint_clean = jp1)
}


caprese_capri_cv <- function(model, evaleloss, evalprederr, evalposterr,
                             cores.ratio = 1e-6) {
    message("  Starting caprese_capri_cv ", date())
    ## overall eloss
    if(evaleloss) {
        model <- tronco.kfold.eloss(model, runs = 10, k = 10)
        eloss <- as.kfold.eloss(model)
    } else {
        eloss <- NA
    }
    ## These are per child node, not per relationship
    if(evalprederr) {
        model <- tronco.kfold.prederr(model, runs = 10, k = 10,
                                      cores.ratio = cores.ratio) ## all child nodes
        prederr <- as.kfold.prederr(model)
    } else {
        prederr <- NA
    }
    if(evalposterr) {
        model <- tronco.kfold.posterr(model, runs = 10, k = 10,
                                      cores.ratio = cores.ratio)
        posterr <- as.kfold.posterr(model)
    } else {
        posterr <- NA
    }
    message("  Ending caprese_capri_cv ", date())
    return(list(eloss = eloss,
                prederr = prederr,
                posterr = posterr))
}

caprese_proc <- function(genotData = NULL, nboot = 200, datax = NULL,
                         silent = TRUE,
                         evaleloss   = TRUE,
                         evalprederr = TRUE,
                         evalposterr = TRUE,
                         cores.ratio = 1e-6) {
    if(is.null(datax) && is.null(genotData))
        stop("one of data or genotData needed")
    if( !is.null(datax) && !is.null(genotData) )
        warning("data and genotData passed; we will process data")
    if(!is.null(datax))
        genotData <- import.genotypes(datax)
    message("  Starting tronco.caprese ", date())
    model <- tronco.caprese(genotData, silent = silent)
    ## message("  Starting tronco.bootstrap, in caprese ", date())
    ## boot <- tronco.bootstrap(model, nboot = nboot, silent = silent,
    ##                          cores.ratio = cores.ratio)
    ## message("  Ending tronco.bootstrap, in capri ", date())
    ## original.adj.mat <- as.adj.matrix(model)$caprese
    ## boot.adj.mat <- boot$bootstrap$caprese$npb$bootstrap.adj.matrix$frequency
    ## confidences <- as.confidence(boot, conf = c("pr", "tp", "npb", "hg"))
    ## bt <- confidences$npb$caprese ## only used for checking

    if(nboot > 0) {
        message("  Starting tronco.bootstrap, in caprese ", date())
        boot <- tronco.bootstrap(model, nboot = nboot, silent = silent,
                                 cores.ratio = cores.ratio)
        message("  Ending tronco.bootstrap, in capri ", date())
        original.adj.mat <- as.adj.matrix(model)$caprese
        boot.adj.mat <- boot$bootstrap$caprese$npb$bootstrap.adj.matrix$frequency
        confidences <- as.confidence(boot, conf = c("pr", "tp", "npb", "hg"))
        bt <- confidences$npb$caprese ## only used for checking
    } else {
        boot <- NA
        original.adj.mat <- as.adj.matrix(model)$caprese
        confidences <- as.confidence(model, conf = c("pr", "tp", "hg"))
    }

    pr <- confidences$pr
    tp <- confidences$tp
    hg <- confidences$hg
    ## eloss <- confidences$eloss
    edgeBootFreq <- "somedummytoshutupnovisiblebinding"
    if(nboot > 0)
        out <- caprese_capri_common(original.adj.mat, boot.adj.mat, bt, pr, tp, hg,
                                    genotData)
    else
        out <- caprese_capri_common_no_boot(original.adj.mat, pr, tp, hg,
                                            genotData)
    cvp <- caprese_capri_cv(model, evaleloss, evalprederr, evalposterr,
                            cores.ratio = cores.ratio)
    out <- dplyr::rename(out,
                  CAPRESE_edgeBootFreq = edgeBootFreq,
                  CAPRESE_pr = pr,
                  CAPRESE_tp = tp,
                  CAPRESE_hg = hg
                  )
    
    ## Could use with(cvp, get0(whatever, ifnotfound = NA))
    ## but more cumbersome
        ## Since CAPRI and CAPRESE are so tricky we better return the whole object
    ## with data removed

    ## pos.genot.model <- which(names(model) == "genotypes")
    ## pos.genot.boot <- which(names(boot) == "genotypes")

    if(nboot == 0) boot <- NULL

    return(list(edges = out,
                CAPRESE_eloss = if(evaleloss){cvp$eloss} else {NA},
                CAPRESE_prederr = if(evalprederr){cvp$prederr$caprese} else {NA},
                CAPRESE_posterr = if(evalposterr){cvp$posterr$caprese} else {NA},
                ## CAPRESE_probs = estimated_probs(model),
                CAPRESE_useless_extra = useless_extra(model, boot),
                ## CAPRESE_whole_model = model[-pos.genot.model],
                ## CAPRESE_whole_boot = boot[-pos.genot.boot],
                TRONCO_model_boot = boot,
                TRONCO_model = model
                ))
}




capri_proc <- function(genotData = NULL, nboot = 100,
                       nbootWilcox = 100, datax = NULL,
                       regularization = "bic", command = "hc",
                       pvalue = 0.05, silent = TRUE,
                       evaleloss   = TRUE,
                       evalprederr = TRUE,
                       evalposterr = TRUE,
                       cores.ratio = 1e-6) {
    if(is.null(datax) && is.null(genotData))
        stop("one of data or genotData needed")
    if( !is.null(datax) && !is.null(genotData) )
        warning("datax and genotData passed; we will process datax")
    if(!is.null(datax))
        genotData <- import.genotypes(datax)
    message("  Starting tronco.capri ", date())
    model <- tronco.capri(genotData, command = command,
                          regularization = regularization,
                          nboot = nbootWilcox, pvalue = pvalue,
                          silent = silent)


    
    ## message("  Starting tronco.bootstrap, in capri ", date())
    ## boot <- tronco.bootstrap(model, nboot = nboot, silent = silent)
    ## message("  Ending tronco.bootstrap, in capri ", date())
    ## confidences <- as.confidence(boot, conf = c("pr", "tp", "npb", "hg"))
    ## if (regularization == "bic") {
    ##     original.adj.mat <- as.adj.matrix(model, type = "fit")$capri_bic
    ##     boot.adj.mat <- boot$bootstrap$capri_bic$npb$bootstrap.adj.matrix$frequency
    ##     bt <- confidences$npb$capri_bic
    ## } else if (regularization == "aic") {
    ##     original.adj.mat <- as.adj.matrix(model, type = "fit")$capri_aic
    ##     boot.adj.mat <- boot$bootstrap$capri_aic$npb$bootstrap.adj.matrix$frequency
    ##     bt <- confidences$npb$capri_aic
    ## }
    if(nboot > 0) {
        message("  Starting tronco.bootstrap, in capri ", date())
        boot <- tronco.bootstrap(model, nboot = nboot, silent = silent)
        message("  Ending tronco.bootstrap, in capri ", date())
        confidences <- as.confidence(boot, conf = c("pr", "tp", "npb", "hg"))
        if (regularization == "bic") {
            original.adj.mat <- as.adj.matrix(model, type = "fit")$capri_bic
            boot.adj.mat <- boot$bootstrap$capri_bic$npb$bootstrap.adj.matrix$frequency
            bt <- confidences$npb$capri_bic
        } else if (regularization == "aic") {
            original.adj.mat <- as.adj.matrix(model, type = "fit")$capri_aic
            boot.adj.mat <- boot$bootstrap$capri_aic$npb$bootstrap.adj.matrix$frequency
            bt <- confidences$npb$capri_aic
        }
    } else {
        boot <- NA
        confidences <- as.confidence(model, conf = c("pr", "tp", "hg"))
        if (regularization == "bic") {
            original.adj.mat <- as.adj.matrix(model, type = "fit")$capri_bic
        } else if (regularization == "aic") {
            original.adj.mat <- as.adj.matrix(model, type = "fit")$capri_aic
        }
    }

    pr <- confidences$pr
    tp <- confidences$tp
    hg <- confidences$hg
    edgeBootFreq <- "somedummytoshutupnovisiblebinding"
    if(nboot > 0)
        out <- caprese_capri_common(original.adj.mat, boot.adj.mat, bt, pr, tp, hg,
                                    genotData)
    else
        out <- caprese_capri_common_no_boot(original.adj.mat, pr, tp, hg,
                                    genotData)
    cvp <- caprese_capri_cv(model, evaleloss, evalprederr, evalposterr,
                            cores.ratio = cores.ratio)

    out <- dplyr::rename(out,
                  CAPRI_edgeBootFreq = edgeBootFreq,
                  CAPRI_pr = pr,
                  CAPRI_tp = tp,
                  CAPRI_hg = hg
                  )
    ## Since CAPRI and CAPRESE are so tricky we better return the whole object
    ## with data removed
    ## pos.genot.model <- which(names(model) == "genotypes")
    ## pos.genot.boot <- which(names(boot) == "genotypes")

    if(nboot == 0) boot <- NULL
    
    return(list(edges = out,
                CAPRI_eloss = if(evaleloss){cvp$eloss} else {NA},
                CAPRI_prederr = if(evalprederr){cvp$prederr$capri} else {NA},
                CAPRI_posterr = if(evalposterr){cvp$posterr$capri} else {NA},
                ## CAPRI_probs = estimated_probs(model),
                CAPRI_useless_extra = useless_extra(model, boot),
                ## CAPRI_whole_model = model[-pos.genot.model],
                ## CAPRI_whole_boot = boot[-pos.genot.boot],
                TRONCO_model_boot = boot,
                TRONCO_model = model
                ))
}


library(codetools)
checkUsageEnv(env = .GlobalEnv)



