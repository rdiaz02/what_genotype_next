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


options(mc.cores = 1L)
options(boot.parallel = "no")
options(boot.ncpus = 1L)

## rm(list = ls())
library(Oncotree)
library(help = Oncotree)


## Also return the frequency of the original tree
## see
## > print.boottree
## function (x, ...) 
## {
##     orig.string <- paste(x$original$parent.num, collapse = ".")
##     boot.idx <- match(orig.string, as.character(x$tree.list$Tree))
##     cat("Out of the", sum(x$tree.list$Freq), "replicates", "there are", 
##         nrow(x$tree.list), "unique trees with frequencies from", 
##         max(x$tree.list$Freq), "down to", min(x$tree.list$Freq), 
##         "\n")
##     cat("The bootstrap process found the original tree", x$tree.list$Freq[boot.idx], 
##         "times\n")
##     invisible(x)
## }

f_ot <- function(x, nboot = 100) {
    datax <- x$out$popSample
    datax <- pre_process(datax, remove.constant = FALSE)
    if(ncol(datax) < 2) {
        return(list(name = x$name, params = x$params,
                    nocols = TRUE,
                    time = NA,
                    res = NA))
    } else {
        time <- system.time({
            otout <- try(ot_proc(datax, nboot = nboot))
        })
        return(list(name = x$name, params = x$params,
                    nocols = FALSE,
                    time = time,
                    res = otout))
    }
}

ot_consensus_sb <- function(x) {
    ## Taking the code from plot.boottree
    child <- x$original$child
    parent.num <- as.numeric(x$consensus)
    parent <- character()
    nmut <- length(child)
    for (i in 1:nmut) {
        number <- parent.num[i]
        if (i == 1) {
            parent[i] <- ""
        }
        else {
            parent[i] <- child[number]
        }
    }

    ## why was this here??
    ## mostfreq <- list(child = child, parent = parent,
    ##                  parent.num = parent.num)

    
    ## Now, extract the From and To in my usual way
    edges.matrix <- cbind(parent = parent,
                          child = child)
    ## Just removing the "" -> Root
    if(identical(edges.matrix[1, ], c(parent = "", child = "Root"))){
        edges.matrix <- edges.matrix[-1, , drop = FALSE]
    } else {
        stop("this is not expected; where is 1, 2 as Root?")
    }
    
    return(data.frame(From = edges.matrix[, "parent"],
                      To = edges.matrix[, "child"],
                      edge = paste(edges.matrix[, "parent"],
                                   edges.matrix[, "child"],
                                   sep = " -> "),
                      stringsAsFactors = FALSE))
}


## You can verify this function by running examples
## and doing plot(boot, draw.consensus = TRUE, minfreq = 100)
ot_consensus_mine <- function(fit, boot) {
    cn <- colnames(fit$data)
    child <- boot$original$child
    stopifnot(child[1] == "Root")
    child <- child[-1]
    consensus <- boot$consensus
    stopifnot(consensus[1] == 0)
    stopifnot(sum(consensus == 0) == 1)
    parent <-cn[consensus[-1]]
    return(data.frame(
        From = parent,
        To = child,
        edge = paste(parent, child, sep = " -> "),
        stringsAsFactors = FALSE
    ))
}
## Some checks
check_ot_consensus <- function(data, nb) {
    otf <- oncotree.fit(data)
    otb <- bootstrap.oncotree(otf, R = nb, type = "nonparametric")
    par(mfrow = c(2, 1))
    plot(otb, draw.consensus = TRUE, minfreq = 100)
    cat("\n mine\n")
    print(mine <- ot_consensus_mine(otf, otb))
    cat("\n sb\n")
    print(otsb <- ot_consensus_sb(otb))
    identical(otsb, mine)
}
## and do, for instance
## check_ot_consensus(ov.cgh, 5)


ot_consensus <- ot_consensus_sb

ot_proc <- function(datax, nboot = 1000, distribution.oncotree = TRUE) {

    ## Like ot_proc, but gives frequency of original tree
    error.fun <- "std"
    message(" Starting ot.fit ", date())
    ot.fit <- try(oncotree.fit(datax))
    if(inherits(ot.fit, "try-error")) {
        error.fun <- "tryNULL"
        ot.fit <- try(oncotree.fit(datax, error.fun = NULL))
    }
    message(" Done ot.fit ", date())
    edges.matrix <- cbind(parent = ot.fit$parent$parent,
                          child = ot.fit$parent$child)
    edge.weights <- ot.fit$parent$est.weight
    
    if(identical(edges.matrix[1, ], c(parent = "", child = "Root"))){
        edges.matrix <- edges.matrix[-1, , drop = FALSE]
        edge.weights <- edge.weights[-1]
        if(is.null(edge.weights)) {
            ## no est.weight, probably because error.fun is NULL
            ## use observed
            ## Was not in version for Baseline, but was never needed as
            ## we never faced this issue.
            edge.weights <- ot.fit$parent$obs.weight[-1]
            error.fun <- paste(error.fun, "observed.weights", sep = ",")
        }
    } else {
        stop("this is not expected; where is 1, 2 as Root?")
    }

    ## Marginal freqs, observed and fitted.
    ## We will add them in the output associated with the "To" node,
    ## as these are the predicted frequencies of the receiving node.
    obs_marginal <- colSums(datax)/nrow(datax)
    pred_marginal <- try(marginal.distr(ot.fit))
    if(inherits(pred_marginal, "try-error")) {
        pred_marginal <- marginal.distr(ot.fit, with.errors = FALSE)
        error.fun <- paste(error.fun, "pred.marginal.no.error", sep = ",")
    }
    if(nboot > 0) {
        message(" Starting bootstrap.oncotree ", date())
        ot.boot <- bootstrap.oncotree(ot.fit,
                                      R = nboot, 
                                      type = 'nonparametric')
        message(" Done bootstrap.oncotree ", date())
        ot.boot.freqs <- ot.boot$parent.freq
        ## From print.boottree, to return freq of original tree
        orig.string <- paste(ot.boot$original$parent.num, collapse = ".")
        boot.idx <- match(orig.string, as.character(ot.boot$tree.list$Tree))
        ot.boot.original <- ot.boot$tree.list$Freq[boot.idx]/nboot

        ## Two paranoid checks
        if(! all(edges.matrix[, 1, drop = FALSE] %in% colnames(ot.boot.freqs)))
            stop("colnames ot.boot.freqs weird")
        if(! all(edges.matrix[, 2, drop = FALSE] %in% rownames(ot.boot.freqs)))
            stop("rownames ot.boot.freqs weird")
        boot.freq <- ot.boot.freqs[edges.matrix]/nboot
        consensus <- try(ot_consensus(ot.boot))
    } else {
        boot.freq <- consensus <- ot.boot.original <- NA
    }
        
    
    if(distribution.oncotree) {
        ## with many genotypes and/or large trees can lead to unexpectedly
        ## very long computing times. Removed for now.
        ## Well, I think that only happens if with.errors = TRUE
        ## which we probably don't want anyway.
        message(" Starting distribution.oncotree ", date())
        ## Observed and expected frequencies of genotypes
        est_genots <- distribution.oncotree(ot.fit, with.probs = TRUE,
                                            with.errors = FALSE) ## TRUE)
        tt <- as.data.frame(datax)
        obs_genots <- aggregate(tt, by = tt, length)[1:(ncol(tt) + 1)]
        colnames(obs_genots)[ncol(obs_genots)] <- "Counts"
        message(" Ending distribution.oncotree ", date())

        message(" Starting observed vs expected, oncotree ", date())
        ## Observed and expected, 2 events. From vignette
        est2way <- t(data.matrix(est_genots[2:(ncol(est_genots) - 1)])) %*% diag(est_genots$Prob) %*%
            data.matrix(est_genots[2:(ncol(est_genots) - 1)])
        obs2way <-t(ot.fit$data[,-1]) %*% ot.fit$data[,-1]/nrow(ot.fit$data)
        message(" Ending observed vs expected, oncotree ", date())
    } else {
        est_genots <- NA
        est2way <- NA
        obs_genots <- NA
        obs2way <- NA
    }
    
    return(list(edges = data.frame(From = edges.matrix[, "parent"],
                      To = edges.matrix[, "child"],
                      edge = paste(edges.matrix[, "parent"],
                                   edges.matrix[, "child"],
                                   sep = " -> "),
                      OT_edgeBootFreq = boot.freq,
                      OT_edgeWeight = edge.weights,
                      OT_obsMarginal = obs_marginal[edges.matrix[, "child"]],
                      OT_predMarginal = pred_marginal[edges.matrix[, "child"]],
                      stringsAsFactors = FALSE),
                consensus = consensus, ## ot_consensus(ot.boot),
                OT_error.fun  = error.fun,
                ot.boot.original = ot.boot.original, ## Frequency of original tree among boot
                genots_predicted = est_genots,
                genots_observed = obs_genots,
                two_way_predicted = est2way,
                two_way_observed = obs2way
                      ))
}

## For genotypes, obs and pred
## then match dd and the tt table, but I might as well
## paste and tabulate and compare, create a long table with indicator,
## do a chi-square, and test and find largest diffs


## require(gtools)
## require(dplyr)
## freq_combs <- function(x, joint = 3) {
##     ## First columns are genes, last is Prob or Freq
##     ## Beware to remove the 1, if coming from OT expected

##     ## When using this, make sure the observed and predicted have columns
##     ## ordered the same way.

    
##     ## joint: 3 for three-way, etc.
##     ng <- ncol(x) - 1
##     cc <- combinations(ng, joint)

##     yes_comb <- function(z, s = joint) {
##         ## these have the given joint combination
##         which(rowSums(z[, , drop = FALSE]) == s)
##     }


##     these_cols_comb <- function(cols, data = x) {
##         joint <- length(cols)
##         ## data has the genotypes and the last column is freq
##         rr <- yes_comb(data[, cols, drop = FALSE], s = joint)
##         s <- 0
##         if(length(rr) > 0) {
##             s <- sum(data[rr, ncol(data)])
##         }
##         nn <- paste(colnames(data)[cols], collapse = ",")
##         return(data.frame(comb = nn, freq = s, stringsAsFactors = FALSE))
##     }
##     browser()
##     return(dplyr::bind_rows(apply(cc, 1,
##                                   function(u) these_cols_comb(u, x))))
## }




## The "drop = FALSE" here but not in old versions: if we never have a
## single column of data, then edges.matrix always has at least two rows.

## In former versions, pre_process would not allow passing a constant
## column. So any column would have at least two states. And all the
## functions only attempt anything if the number of columns is at least
## 2. So with two columns, one constant, we would never call any analysis
## function.


## pre_process now allows passing a constant column. So we have two
## columns. Analysis is allowed. But then OT removes the constant column.


## Why do we allow now data with only a single column? Because we now are
## trying to also recover the Root -> something. If you know you are going
## to exclude those, then data sets with only 1 column make no sense.




## FIXME: to add?

## rank order of worst predicted genes?
## small bootstrap
## discrepancy consensus and original


library(codetools)
checkUsageEnv(env = .GlobalEnv)
