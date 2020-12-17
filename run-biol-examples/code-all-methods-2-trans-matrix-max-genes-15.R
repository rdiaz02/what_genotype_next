## Copyright 2016, 2017, 2018, 2020 Ramon Diaz-Uriarte

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



## Given a data set (patients as rows, genes as columns) return the
## transition matrices between genotypes according to all methods.

## Done by function: all_methods_2_trans_mat
## almost at the bottom.


## You need to have
## - Schills code
## - MCCBN installed: from git, then R CMD INSTALL: https://github.com/cbg-ethz/MC-CBN
## - CBN: use my version with fixes, in the repo:
##     cd to ct-cbn-0.1.04b-with-rdu-bug-fix-write-lambda-likelihood
##     then:
##       ./configure
##       make
##       and cp, mv, ln the two binaries to a place in the path
##      BEWARE! sometimes weird things can happen (make not working, etc)
##        if you have stale files. 
##        If things break, the simplest is to checkout a new copy from the repo
##        and do the ./configure , make dance there.
##        Probably you need to have autoconf-archive installed.
## - OT, and CAPRIand CAPRESE (the corresponding packages Oncotree and TRONCO)

date()

## Since we are not parallelizing here, you might want to set
thiscores <- 1
## export OPENBLAS_NUM_THREADS=thiscores
## export OMP_NUM_THREADS=thiscores
## with the 36 substituted by the number of your CPUs
## in the script that launches analyses (unless those runs are parallelized)
## I also do it here for R itself
library(RhpcBLASctl)
RhpcBLASctl::blas_set_num_threads(thiscores)
RhpcBLASctl::omp_set_num_threads(thiscores)


library(compiler)
enableJIT(3)
setCompilerOptions(optimize = 3)
## ## library(parallel)

library(dplyr)
library(Oncotree)
library(readr)
library(data.table)
library(igraph)
library(TRONCO)
library(parallel)
library(foreach)
library(doParallel)
library(Rgraphviz)
library(stringr)
library(pryr)
library(OncoSimulR)
library(testthat)
library(Matrix)

source("schill-trans-mat.R") ## yes, a few seconds because of the testing

## source("LOD-POM.R", echo = FALSE, max.deparse.length = 0)
source("pre-process.R", echo = FALSE, max.deparse.length = 0)
source("capri-caprese-process.R", echo = FALSE, max.deparse.length = 0)
source("ot-process.R", echo = FALSE, max.deparse.length = 0)
## source("dip-process.R", echo = FALSE, max.deparse.length = 0)  ## No longer using it
source("cbn-process.R", echo = FALSE, max.deparse.length = 0)
source("mccbn-process.R", echo = FALSE, max.deparse.length = 0)

registerDoSEQ() ## for TRONCO

boot_data_index <- function(x, boot) {
    ## Used by CBN and DiP
    ## boot is an integer. 0 means no boot
    ## that is because I reuse boot for two purposes
    boot <- as.logical(boot)
    if(boot) {
        ind <- sample(nrow(x), nrow(x), replace = TRUE)
        return(x[ind, , drop = FALSE])
    } else {
        return(x)
    }
}

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
        if(fmax == 1) {
            cat("\n  Added one pseudosample \n ")
            n00 <- 1
        } else { 
            n00 <- 0
        }
    }
    return(rbind(x,
                 matrix(0L, nrow = n00, ncol = ncol(x))
                 ))
    ## cn <- colnames(x)
    
    ## tmp <- rbind(x,
    ##              matrix(0L, nrow = n00, ncol = ncol(x))
    ##              )
    ## colnames(tmp) <- cn
    ## return(tmp)
}

## NOTE: MCCBN allowed to run with arbitrary number of columns
all_methods <- function(x, nboot = 0, nboot_caprese_capri = 0,
                        nboot_cbn = 0, nboot_dip = 0,
                        n00 = "auto3", caprese_capri_minimal = TRUE,
                        caprese_capri_cores.ratio = 0,
                        distribution_oncotree = FALSE,
                        min.freq = 0,
                        cores_cbn = 1) { ## I think we want to keep min.freq to 0?
    
    if(caprese_capri_minimal) {
        ## Most of the time we only want the graphs, that's it
        evaleloss <- FALSE
        evalprederr <- FALSE
        evalposterr <- FALSE
    } else {
        evaleloss <- TRUE
        evalprederr <- TRUE
        evalposterr <- TRUE
    }
    x000 <- x
    x <- add_pseudosamples(x, n00 = n00)
    ## remove.constant makes no difference IFF we add pseudosamples, as
    ## there can be no constant column when we add pseudosamples
    x_tmp <- pre_process(x, remove.constant = FALSE, min.freq = min.freq)

    cat("\n     Doing OT")
    if(ncol(x_tmp) >= 2) {
        time_ot <- system.time(OT <- try(
                                   suppressMessages(
                                       ot_proc(x_tmp,
                                               nboot = nboot,
                                               distribution.oncotree = distribution_oncotree))))["elapsed"]
    } else {
        OT <- NA
        time_ot <- NA
    }

    ## remove.constant I think is here because if we do not add
    ## pseudosamples, caprese and capri will fail (fail, not just give a
    ## tree with a silly edge) if we pass n00 = 0 to add
    ## pseudosamples. And in some cases we have run this without adding
    ## pseudosamples. Not with the simulations, though.
    x_tmp <- pre_process(x, remove.constant = TRUE, min.freq = min.freq)
    if(ncol(x_tmp) >= 2) {
        x_tmp <- tronco_common(x_tmp)
        cat("\n     Doing CAPRESE")
        time_caprese <- system.time(CAPRESE <- try(suppressMessages(
                                        caprese_proc(genotData = x_tmp,
                                    nboot = nboot_caprese_capri,
                                    silent = TRUE,
                                    cores.ratio = caprese_capri_cores.ratio,
                                    evaleloss = evaleloss,
                                    evalprederr = evalprederr,
                                    evalposterr = evalposterr))))["elapsed"]
        cat("\n     Doing CAPRI_BIC")
        time_capri_bic <- system.time(CAPRI_BIC <- try(suppressMessages(
                                      capri_proc(genotData = x_tmp,
                                    nboot = nboot_caprese_capri,
                                    nbootWilcox = 100,
                                    regularization = "bic",
                                    command = "hc",
                                    pvalue = 0.05,
                                    cores.ratio = caprese_capri_cores.ratio,
                                    evaleloss = evaleloss,
                                    evalprederr = evalprederr,
                                    evalposterr = evalposterr))))["elapsed"]
        cat("\n     Doing CAPRI_AIC")
        time_capri_aic <- system.time(CAPRI_AIC <- try(suppressMessages(
                                      capri_proc(genotData = x_tmp,
                                    nboot = nboot_caprese_capri,
                                    nbootWilcox = 100,
                                    regularization = "aic",
                                    command = "hc",
                                    pvalue = 0.05,
                                    cores.ratio = caprese_capri_cores.ratio,
                                    evaleloss = evaleloss,
                                    evalprederr = evalprederr,
                                    evalposterr = evalposterr))))["elapsed"]
    } else {
        CAPRESE <- CAPRI_AIC <- CAPRI_BIC <- NA
        time_caprese <- time_capri_aic <- time_capri_bic <- NA
    }

    ## The code (ct-cbn.h) has
    ## 	if ((n < 1) || (n > 25))
    ## {
    ## 	fprintf(stderr, "Error:  Number of events is %d.  Supported range is {1, ..., 14}.\n", n);
    ## 	exit(1);
    ## }
    ## So despite the message, we can do up to 25?
    x_tmp <- pre_process(x, remove.constant = FALSE, min.freq = min.freq,
                         max.cols = 25)
    if(ncol(x_tmp) >= 2) {
        ## CBN_linear <- try(cbn_proc(x_tmp, addname = "tmpl",
        ##                            init.poset = "linear",
        ##                            nboot = nboot_cbn, parall = TRUE))
        ## CBN_linear <- NA ## disabled
        ## And note the parall argument has no effect as I am not using mclapply
        ## This was changed in cbn-process.R on commit b9027d5, on 2016-09-26.
        ## Yes, ugly!
        
        cat("\n     Doing CBN")
        time_cbn_ot <- system.time(CBN_ot <- try(cbn_proc(x_tmp, addname = "tmpo",
                                   init.poset = "OT",
                                   nboot = nboot_cbn, parall = TRUE,
                                   cores = cores_cbn)))["elapsed"]
    } else {
       time_cbn_ot <- NA
       CBN_ot <- NA 
    }

    x_tmp <- pre_process(x, remove.constant = FALSE, min.freq = min.freq,
                         max.cols = NULL)
    cat("\n     Doing MCCBN")
    if(ncol(x_tmp) >= 2) {
        ## MCCBN should be able to run with many more columns. See max.cols above
        time_mccbn <- system.time(MCCBN <- try(mccbn_proc(x_tmp)))["elapsed"]
    } else {
        ## MCCBN <- CBN_linear <- CBN_ot <- NA
        time_mccbn <- NA
        MCCBN <- NA
    }
    ## x_tmp <- pre_process(x, remove.constant = FALSE, min.freq = min.freq,
    ##                      max.cols = 20)
    ## if(ncol(x_tmp) >= 2) {
    ##     ## if we are really going to use boottstrap,
    ##     ## it is more reasonable to parallelize over boot
    ##     ## than over p-range
    ##     if(nboot_dip == 0) {
    ##         ## When many runs from shell, we do not want any parall
    ##         ## cpr <- detectCores()
    ##         cpr <- 1
    ##         cpb <- 1
    ##         pb <- FALSE
    ##     } else {
    ##         cpr <- 1
    ##         cpb <- detectCores()
    ##         pb <- TRUE
    ##     }
    ##     time_dip_mpn <- system.time(DIP_MPN <- try(dip_proc(x_tmp, method = "MPN",
    ##                             addname = "dm",
    ##                             nboot = nboot_dip,
    ##                             cores_for_p_range = cpr,
    ##                             cores_for_boot = cpb,
    ##                             parall_boot = pb)))["elapsed"]
    ##     time_dip_smpn <- system.time(DIP_SMPN <- try(dip_proc(x_tmp, method = "SMPN", addname = "ds",
    ##                              nboot = nboot_dip,
    ##                              cores_for_p_range = cpr,
    ##                              cores_for_boot = cpb,
    ##                              parall_boot = pb)))["elapsed"]
    ## } else {
    ##     DIP_MPN <- DIP_SMPN <- NA
    ##     time_dip_mpn <- time_dip_smpn <- NA
    ## }

    cat("\n                  _times_cpm: ot", time_ot,
        " cbn ", time_cbn_ot,
        " mccbn ", time_mccbn,
        " caprese ", time_caprese,
        " capri_aic ", time_capri_aic,
        " capri_bic ", time_capri_bic,
        "\n")
    out <- list(
        OT = OT,
        CAPRESE = CAPRESE,
        CAPRI_BIC = CAPRI_BIC,
        CAPRI_AIC = CAPRI_AIC,        
        ## DIP_MPN = DIP_MPN,
        ## DIP_SMPN = DIP_SMPN,
        CBN_ot = CBN_ot,
        MCCBN = MCCBN,
        time_ot = time_ot,
        time_caprese = time_caprese,
        time_capri_bic = time_capri_bic,
        time_capri_aic = time_capri_aic,        
        ## time_dip_mpn = time_dip_mpn,
        ## time_dip_smpn = time_dip_smpn,
        time_cbn_ot = time_cbn_ot,
        time_mccbn = time_mccbn,
        input_data = x000, ## yes, return this!!!
        input_data_pseudosamples = x
    )
    return(out)
}




## This is coming from
## Cancer_Data_sets/CPM-weighted-paths-biol.R
## in the supplementary material for Diaz-Uriarte and Vasallo.
## We leave in there things we don't really need. Simpler.


## DAG of restrictions (as data frame) -> vector of accessible genotypes and graph of DAG of restrictions
## return all the accessible genotypes
##     from a DAG of genes plus the DAG as igraph object
df_2_access_genots_and_graph <- function(x) {
    
    ## minor detail: if x is a sparse adjacency mat.
    ##    from igraph this still works. But don't do that.
    g <- igraph::graph_from_data_frame(x, directed = TRUE)
    
    ## g: an igraph object, with a "Root"
    ## returns a list
    children <- sort(setdiff(V(g)$name, "Root"))
    node_depth <- unlist(lapply(children,
                                function(node)
                                    max(unlist(lapply(all_simple_paths(g,
                                                                       from = "Root",
                                                                       to = node),
                                                      length)))
                                ))
    
    names(node_depth) <- children
    node_depth <- sort(node_depth)
    ## pre-allocate a list.
    ## FIXME: Could be smarter as a function of dim(x)?
    all_gty <- vector(mode = "list", length = 100) 
    i <- 1
    for(j in seq_along(node_depth)) {
        tmp_gty_1 <- sort(setdiff(names(subcomponent(g, v = names(node_depth)[j],
                                                     mode = "in")),
                                  "Root"))
        all_gty[[i]] <- tmp_gty_1
        
        ## only do union of those not contained in the genotype
        ## to_unite <-  which(unlist(lapply(all_gty[1:(i-1)],
        ##                                  function(z) length(setdiff(z, tmp_gty_1)) > 0)))
        
        if(i > 1) {
            to_unite <-  which(unlist(lapply(all_gty[1:(i-1)],
                                             function(z) !all(z %in% tmp_gty_1))))
        } else {
            to_unite <- vector(length = 0)
        }
        
        if(length(to_unite)) {
            ## we need unique as some sets are the same
            tmp_gty_2 <- unique(lapply(all_gty[to_unite],
                                       function(u) sort(union(u, tmp_gty_1))))
            ## check remaining space in preallocated list.
            ## expand if needed.
            if(length(all_gty) < (i + 1 + length(tmp_gty_2))) {
                all_gty <- c(all_gty,
                             vector(mode = "list", length = max(length(all_gty),
                                                                2 + length(tmp_gty_2))))
            }
            all_gty[(i + 1):(i + length(tmp_gty_2))] <- tmp_gty_2
            i <- i + length(tmp_gty_2)
        }
        i <- i + 1
    }
    all_gty <- all_gty[1:(i - 1)]
    ng <- unlist(lapply(all_gty, length))
    all_gty <- all_gty[order(ng)]
    return(list(accessible_genots = all_gty,
                graph = g))
}


## From function of same name in ruggify-functions.R

## vector of accessible genotypes -> adjacency matrix of genotypes (fitness graph)
## return maximally connected fitness graph for a given set of accessible genotypes
unrestricted_fitness_graph <- function(gacc, plot = FALSE) {
    
    gs <- unlist(lapply(gacc, function(g) paste0(g, collapse = ", ")))
    gs <- c("WT", gs)
    nmut <- c(0, vapply(gacc, length, 1))
    
    adjmat <- matrix(0L, nrow = length(gs), ncol = length(gs))
    rownames(adjmat) <- colnames(adjmat) <- gs
    
    adjmat["WT", gs[which(nmut == 1)]] <- 1L
    
    for(m in 2:max(nmut)){
        g <- gs[which(nmut == m)]
        for (gn in g) {
            parents <- gacc[which(nmut == m-1)-1]
            gns <- unlist(strsplit(gn, ", "))
            parents <- parents[which(unlist(lapply(parents,
                                                   function(p)
                                                       length(setdiff(gns, p)))) == 1)]
            for (p in parents){
                adjmat[paste0(p, collapse = ", "), gn] <- 1L
            }
        }
    }
    if (plot)
        mccbn::plot_poset(adjmat) ## , title = "G0 (unrestricted)")

    stopifnot(all(adjmat %in% c(0L, 1L) ))
    storage.mode(adjmat) <- "integer"

    return(adjmat)
}



## Trying sparse matrices
## https://stackoverflow.com/questions/23107837/r-sparse-matrix-from-list-of-dimension-names
## https://stackoverflow.com/questions/26207850/create-sparse-matrix-from-a-data-frame
## https://cmdlinetips.com/2019/05/introduction-to-sparse-matrices-in-r/
## https://www.gormanalysis.com/blog/sparse-matrix-construction-and-use-in-r/

## vector of accessible genotypes -> adjacency matrix of genotypes (fitness graph)
## return maximally connected fitness graph for a given set of accessible genotypes
unrestricted_fitness_graph_sparseM <- function(gacc, plot = FALSE) {
    gs <- unlist(lapply(gacc, function(g) paste0(g, collapse = ", ")))
    gs <- c("WT", gs)
    nmut <- c(0, vapply(gacc, length, 1))
    
    ## adjmat <- matrix(0L, nrow = length(gs), ncol = length(gs))
    ## rownames(adjmat) <- colnames(adjmat) <- gs
    
    ## ## adjmat <- Matrix(0L, nrow = length(gs), ncol = length(gs),
    ## ##                  sparse = TRUE, dimnames = list(gs, gs))
    ## ## This works but I don't feel comfortable
    ## ## adjmat["WT", gs[which(nmut == 1)]] <- 1L
    ## adjmat[i = rep(ii, length(jj)), j = jj] <- 1L

    jj <- match(gs[which(nmut == 1)], gs)
    ii <- rep.int(match("WT", gs), length(jj))
    adjmat <- sparseMatrix(i = ii, j = jj, x = 1L,
                           dims = c(length(gs), length(gs)), 
                           dimnames = list(gs, gs))

    for(m in 2:max(nmut)){
        g <- gs[which(nmut == m)]
        for (gn in g) {
            parents <- gacc[which(nmut == m-1)-1]
            gns <- unlist(strsplit(gn, ", "))
            parents <- parents[which(unlist(lapply(parents,
                                                   function(p)
                                                       length(setdiff(gns, p)))) == 1)]
            for (p in parents){
                ## Works but better via indices, I think
                ## adjmat[paste0(p, collapse = ", "), gn] <- 1L
                jjj <- match(gn, gs)
                iii <- rep.int(match(paste0(p, collapse = ", "), gs),
                               length(jjj))
                adjmat[iii, jjj] <- 1L
            }
        }
    }
    ## Probably will not work
    if (plot)
        mccbn::plot_poset(adjmat) ## , title = "G0 (unrestricted)")
    ## stopifnot(all(adjmat %in% c(0L, 1L) ))
    ## storage.mode(adjmat) <- "integer"
    return(adjmat)
}



## list of accessible genotypes -> global maximum
##    Beware: just a single maximum
get_global_max <- function(ag) {
    muts <- unlist(lapply(ag, length))
    max_muts <- max(muts)
    ind_max_muts <- which(muts == max_muts)
    if(length(ind_max_muts) != 1) stop("eh??!! ind_max_muts")
    return(paste0(ag[[ind_max_muts]], collapse = ", "))
}


## < /from CPMs-paths-genotypes-and-comb.R >

## output of CPM analysis, string -> accessible genotypes and paths
##                  string: just to identify errors
##    the _w: weights, so we add probs.
##   Modification of function of same name, without _w in
##   CPMs-paths-genotypes-and-comb.R

cpm_access_genots_paths_w <- function(x, string = NULL,
                                    names_weights_paths =
                                        c("rerun_lambda",
                                          "lambda",
                                          "OT_edgeWeight")) {
    if(inherits(x, "try-error") || is.na(x) || is.null(x)) {
        ## The CPM analysis produced no edges component, so
        ## nothing can be done
        if(inherits(x, "try-error")) likely_error <- "Error_in_run"
        if(is.na(x)) likely_error <- "ncol_x"
        if(is.null(x)) likely_error <- "other_error"
        return(list(accessible_genots = "ERROR_CPM_ANALYSIS",
                    num_accessible_genots = "ERROR_CPM_ANALYSIS",
                    CPM_DAG_as_igraph = "ERROR_CPM_ANALYSIS",
                    fgraph = "ERROR_CPM_ANALYSIS",
                    num_paths_to_max = "ERROR_CPM_ANALYSIS",
                    paths_error = "ERROR_CPM_ANALYSIS",
                    paths_to_max = "ERROR_CPM_ANALYSIS",
                    weighted_fgraph = "ERROR_CPM_ANALYSIS",
                    trans_mat_genots = "ERROR_CPM_ANALYSIS",
                    unweighted_paths_to_max = "ERROR_CPM_ANALYSIS",
                    weighted_paths_to_max = "ERROR_CPM_ANALYSIS",
                    diversity_weighted_paths_to_max = "ERROR_CPM_ANALYSIS",
                    likely_error = likely_error
                    
                ))
    }
                         
    x <- x$edges
    tmp <- try(df_2_access_genots_and_graph(x[, c("From", "To")]))
    if(inherits(tmp, "try-error")) {
        stop("how is this happening? there was edges component!")
    } else {
        accessible_genots <- tmp$accessible_genots
        fgraph <- unrestricted_fitness_graph(accessible_genots)
        fgraphi <- igraph::graph_from_adjacency_matrix(fgraph)
        gmax <- get_global_max(accessible_genots)
    }

    ## based on n_pahts_single_max in compute-numpaths-clonal-int-stats.R
    ## I still need this to get the paths from the CPM
    lpaths <- NA
    paths_to_max <- NA
    paths_error <- FALSE
    paths <- try(igraph::all_simple_paths(fgraphi,
                                          from = "WT",
                                          to = gmax,
                                          mode = "out"))
    if(inherits(paths, "try-error")) {
        cat("\n     ERROR_in_paths_in_calling_string = ", string, "\n")
        lpaths <- -99
        paths_error <- TRUE
        paths_to_max <- "PATHS_ERROR" ## this will never appear anywhere FIXME
        ## will have to look for the paths_error variable
    } else {
        lpaths <- length(paths)
        paths_to_max <-  unlist(lapply(paths,
                                       function(x) paste(igraph::as_ids(x),
                                                         collapse = " -> ")))
    }

    
    ## Logic of obtaining weighted fitness graph
    
    ##  - obtain all accessible genotypes and fitness graph (unrestricted
    ##      fitness graph) from CPM

    ##  - from the lambdas/probs. of genes given genes obtain
    ##       probabilities/lambdas of genotypes given genotypes. The call
    ##       to "transition_fg" (where weights is the conditional
    ##       prob./lambda) of descendant gene given parent gene

    ##  - when using CBN, might not be probabilities. Make sure they are
    ##    transition probabilities between genotypes: sweep

    ##  - find the probability of each path: call to
    ##    do_weighted_paths_to_max, that uses the list of all paths_to_max. 

    
    which_col_weights <- which(colnames(x) %in% names_weights_paths)
    if(length(which_col_weights) > 1) {
        stop("more than one column with weights")
    } else if (length(which_col_weights) == 1) {
        stopifnot(colnames(x)[2] == "To")
        weights <- unique(x[, c(2, which_col_weights)])
        if(any(duplicated(weights[, "To"]))) {
            stop("Different lambda/weight for same destination gene.",
                 " This is not allowed with conjunctive DAGs")
        }
        rownames(weights) <- weights[, "To"]
        weighted_fgraph <- transition_fg(fgraph, weights)
    } else {
        ## why would we return something? It is NA
        ## weighted_fgraph <- fgraph
        weighted_fgraph <- NA
    }

    if (length(which_col_weights) == 1) {
        ## weighted_fgraph need not have each row sum to 1, if they are lambdas
        ## from CBN for instance. So make sure they are transition matrices
        ## between genotypes.
        trans_mat_genots <- sweep(weighted_fgraph, 1,
                                  rowSums(weighted_fgraph), FUN = "/")
        trans_mat_genots[is.nan(trans_mat_genots)] <- 0
    } else {
        trans_mat_genots <- NA
    }


    if(length(which_col_weights) == 1) {
        weighted_paths_to_max <- do_weighted_paths_to_max(paths_to_max,
                                                          trans_mat_genots)
        diversity_weighted_paths_to_max <-
            OncoSimulR:::shannonI(weighted_paths_to_max[, "probability"])

    } else {
        weighted_paths_to_max <- NA
        diversity_weighted_paths_to_max <- NA
    }
    
    return(list(accessible_genots = accessible_genots,
                num_accessible_genots = length(accessible_genots),
                CPM_DAG_as_igraph = tmp$graph,
                ## fgraph_AM = fgraph,
                fgraph = fgraph,
                num_paths_to_max = lpaths,
                paths_error = paths_error,
                weighted_fgraph = weighted_fgraph,
                trans_mat_genots = trans_mat_genots,
                unweighted_paths_to_max = paths_to_max,
                weighted_paths_to_max = weighted_paths_to_max,
                diversity_weighted_paths_to_max =
                    diversity_weighted_paths_to_max,
                likely_error = "No_error"
                ))
}





## < /from CPMs-paths-genotypes-and-comb.R >

## output of CPM analysis, string -> accessible genotypes and paths
##                  string: just to identify errors
##    the _w: weights, so we add probs.
##   Modification of function of same name, without _w in
##   CPMs-paths-genotypes-and-comb.R

## Like the one above, but only with necessary output
##  for both speed and size and using sparse matrices.
cpm_access_genots_paths_w_simplified <- function(x, string = NULL,
                                    names_weights_paths =
                                        c("rerun_lambda",
                                          "lambda",
                                          "OT_edgeWeight")) {
    if(inherits(x, "try-error") || is.na(x) || is.null(x)) {
        ## The CPM analysis produced no edges component, so
        ## nothing can be done
        ## if(inherits(x, "try-error")) likely_error <- "Error_in_run"
        ## if(is.na(x)) likely_error <- "ncol_x"
        ## if(is.null(x)) likely_error <- "other_error"
        return(list(## accessible_genots = "ERROR_CPM_ANALYSIS",
                    ## num_accessible_genots = "ERROR_CPM_ANALYSIS",
                    ## CPM_DAG_as_igraph = "ERROR_CPM_ANALYSIS",
                    fgraph = "ERROR_CPM_ANALYSIS",
                    ## num_paths_to_max = "ERROR_CPM_ANALYSIS",
                    ## paths_error = "ERROR_CPM_ANALYSIS",
                    ## paths_to_max = "ERROR_CPM_ANALYSIS",
                    weighted_fgraph = "ERROR_CPM_ANALYSIS",
                    trans_mat_genots = "ERROR_CPM_ANALYSIS"
                    ## unweighted_paths_to_max = "ERROR_CPM_ANALYSIS",
                    ## weighted_paths_to_max = "ERROR_CPM_ANALYSIS",
                    ## diversity_weighted_paths_to_max = "ERROR_CPM_ANALYSIS",
                    ## likely_error = likely_error
                    
                ))
    }
                         
    x <- x$edges
    tmp <- try(df_2_access_genots_and_graph(x[, c("From", "To")]))
    if(inherits(tmp, "try-error")) {
        stop("how is this happening? there was edges component!")
    } else {
        accessible_genots <- tmp$accessible_genots
        fgraph <- unrestricted_fitness_graph_sparseM(accessible_genots)
        ## fgraphi <- igraph::graph_from_adjacency_matrix(fgraph)
        ## gmax <- get_global_max(accessible_genots)
    }

    
    ## ## based on n_pahts_single_max in compute-numpaths-clonal-int-stats.R
    ## ## I still need this to get the paths from the CPM
    ## lpaths <- NA
    ## paths_to_max <- NA
    ## paths_error <- FALSE
    ## paths <- try(igraph::all_simple_paths(fgraphi,
    ##                                       from = "WT",
    ##                                       to = gmax,
    ##                                       mode = "out"))
    ## ## if(inherits(paths, "try-error")) {
    ## ##     cat("\n     ERROR_in_paths_in_calling_string = ", string, "\n")
    ## ##     lpaths <- -99
    ## ##     paths_error <- TRUE
    ## ##     paths_to_max <- "PATHS_ERROR" ## this will never appear anywhere FIXME
    ## ##     ## will have to look for the paths_error variable
    ## ## } else {
    ## ##     lpaths <- length(paths)
    ## ##     paths_to_max <-  unlist(lapply(paths,
    ## ##                                    function(x) paste(igraph::as_ids(x),
    ## ##                                                      collapse = " -> ")))
    ## ## }

    
    ## Logic of obtaining weighted fitness graph
    
    ##  - obtain all accessible genotypes and fitness graph (unrestricted
    ##      fitness graph) from CPM

    ##  - from the lambdas/probs. of genes given genes obtain
    ##       probabilities/lambdas of genotypes given genotypes. The call
    ##       to "transition_fg" (where weights is the conditional
    ##       prob./lambda) of descendant gene given parent gene

    ##  - when using CBN, might not be probabilities. Make sure they are
    ##    transition probabilities between genotypes: sweep

    ##  - find the probability of each path: call to
    ##    do_weighted_paths_to_max, that uses the list of all paths_to_max. 

    which_col_weights <- which(colnames(x) %in% names_weights_paths)
    if(length(which_col_weights) > 1) {
        stop("more than one column with weights")
    } else if (length(which_col_weights) == 1) {
        stopifnot(colnames(x)[2] == "To")
        weights <- unique(x[, c(2, which_col_weights)])
        if(any(duplicated(weights[, "To"]))) {
            stop("Different lambda/weight for same destination gene.",
                 " This is not allowed with conjunctive DAGs")
        }
        rownames(weights) <- weights[, "To"]
        weighted_fgraph <- transition_fg_sparseM(fgraph, weights)
    } else {
        ## why would we return something? It is NA
        ## weighted_fgraph <- fgraph
        weighted_fgraph <- NA
    }

    if (length(which_col_weights) == 1) {
        ## weighted_fgraph need not have each row sum to 1, if they are lambdas
        ## from CBN for instance. So make sure they are transition matrices
        ## between genotypes.
        
        ## trans_mat_genots <- sweep(weighted_fgraph, 1,
        ##                           rowSums(weighted_fgraph), FUN = "/")

        trans_mat_genots <- rowScaleMatrix(weighted_fgraph)
    } else {
        trans_mat_genots <- NA
    }


    ## if(length(which_col_weights) == 1) {
    ##     weighted_paths_to_max <- do_weighted_paths_to_max(paths_to_max,
    ##                                                       trans_mat_genots)
    ##     diversity_weighted_paths_to_max <-
    ##         OncoSimulR:::shannonI(weighted_paths_to_max[, "probability"])

    ## } else {
    ##     weighted_paths_to_max <- NA
    ##     diversity_weighted_paths_to_max <- NA
    ## }
    
    return(list(
        ## accessible_genots = accessible_genots,
        ## num_accessible_genots = length(accessible_genots),
        ## CPM_DAG_as_igraph = tmp$graph,
        ## fgraph_AM = fgraph,
        fgraph = fgraph,
        ## num_paths_to_max = lpaths,
        ## paths_error = paths_error,
        weighted_fgraph = weighted_fgraph,
        trans_mat_genots = trans_mat_genots
        ## unweighted_paths_to_max = paths_to_max,
        ## weighted_paths_to_max = weighted_paths_to_max,
        ## diversity_weighted_paths_to_max =  diversity_weighted_paths_to_max,
        ## likely_error = "No_error"
    ))
}







## adjacency matrix genotypes, row, column of that matrix,
##      lambdas/weights of descendant gene given parent gene ->
##     the lambda of descendant genotype
## Return the lambda/weight of a single genotype coming from a single genotype
get_single_lambda <- function(x, row, col, weights) {
    genot_to <- unlist(strsplit(colnames(x)[col], ", ", fixed = TRUE))
    genot_from <- unlist(strsplit(rownames(x)[row], ", ", fixed = TRUE))
    ## WT stuff
    if((length(genot_from) == 1) && (genot_from == "WT")) {
        genot_from == ""
    }
    added_gene <- setdiff(genot_to, genot_from)
    return(weights[added_gene, 2])
}



## fitness graph, weights of probs/lambdas descendant gene given parent gene -> weighted fitness graph
##             the fitness graph with weights (relative to each node)
##             of jumping to each of the descendant genotypes

## ## No longer used. Turning into transition matrix done outside
## ## 
## ##          transition = TRUE: return transition matrix
## ##          if they sumto zero by row this is really the transition matrix
## ##          of genotypes
transition_fg <- function(x, weights) { ## , transition = TRUE) {
    pos_do <- which(x == 1, arr.ind = TRUE)
    wfg <- x
    wfg[] <- 0
    tmp <- unlist(
        Map(function(r, c)
             get_single_lambda(x, r, c, weights),
            pos_do[, 1], pos_do[, 2]))
    wfg[pos_do] <- tmp
        
    ## ## Could use Map but let's use a loop instead
    ## for(p in 1:nrow(pos_do)) {
    ##     wfg[pos_do[p, ] ] <-  get_single_lambda(x,
    ##                                             pos_do[p, 1], pos_do[p, 2],
    ##                                             weights)
    ## }
    return(wfg)
}


## fitness graph, weights of probs/lambdas descendant gene given parent gene -> weighted fitness graph
##             the fitness graph with weights (relative to each node)
##             of jumping to each of the descendant genotypes

## ## No longer used. Turning into transition matrix done outside
## ## 
## ##          transition = TRUE: return transition matrix
## ##          if they sumto zero by row this is really the transition matrix
## ##          of genotypes
transition_fg_sparseM <- function(x, weights) { ## , transition = TRUE) {
    ## pos_do <- which(x == 1, arr.ind = TRUE)
    pos_do <- as.matrix(summary(x)[, c("i", "j")])
    wfg <- x
    wfg[] <- 0
    tmp <- unlist(
        Map(function(r, c)
             get_single_lambda(x, r, c, weights),
            pos_do[, 1], pos_do[, 2]))
    wfg[pos_do] <- tmp
    ## ## Could use Map but let's use a loop instead
    ## for(p in 1:nrow(pos_do)) {
    ##     wfg[pos_do[p, ] ] <-  get_single_lambda(x,
    ##                                             pos_do[p, 1], pos_do[p, 2],
    ##                                             weights)
    ## }
    return(wfg)
}

## vector of paths to max, transition matrix genotypes ->
##     weighted paths to max
do_weighted_paths_to_max <- function(paths, trans_mat) {
    probs <- vapply(paths, function(u) prob_single_path(u, trans_mat),
                    FUN.VALUE = 0.0)
    stopifnot(isTRUE(all.equal(sum(probs), 1)))
    df <- data.frame(path = paths,
                     probability = probs, stringsAsFactors = FALSE)

    rownames(df) <- NULL
    return(df)
}

## path to max, transition matrix genotypes ->
##    probability of path
prob_single_path <- function(path, trans_mat) {
    ## the indices of the sequence of genotypes
    ii <- which(row.names(trans_mat) %in%
                strsplit(path, " -> ", fixed = TRUE)[[1]])
    prod(trans_mat[cbind(ii[-length(ii)], ii[-1])])
}



## ## file with a single entry of cpm output -> accessible genotypes and paths
## ##    file is an ANALYZED__ID* file, as produced by run-CPMs.R
## cpm_out_to_paths_genots_w <- function(thefile,
##                                       methods = c("CBN", "MCCBN")
##                                       ## methods = c("OT", "CAPRESE",
##                                       ##               "CAPRI_BIC", "CAPRI_AIC",
##                                       ##           "CBN_ot", "MCCBN")
##                                       ) {
##     out <- NULL
##     load(thefile)
##     ## this creates the out, that codetools complaints about
##     outn <- out[1:5]
##     if(is.null(out$bootstrap)) {
##         outn$bootstap <- FALSE
##     } else {
##         outn$bootstrap <- out$bootstrap
##     }

##     string0 <- as.data.frame(outn, stringsAsFactors = FALSE)

##     string <- paste(paste(names(string0), string0, sep = "_"), collapse = "__")

##     out_w <- cpm_access_genots_paths_w(out$cbn_out, string = string)

##     theout <- c(outn,
##                 out_paths_genots = list(out_w)
##                 )
    
##     ## saveRDS(theout, file = paste0("paths_genots_index_", index, "_", string, ".rds"))
##     ## naming more consistent with that of true fitness graphs
##     saveRDS(theout, file = paste0("paths_", string, ".rds"))
##     cat("\n       done string = ", string, "\n")

##     likely_error <- out_w$likely_error
##     paths_error <- out_w$paths_error
    
##     m1 <- data.frame(num_paths_to_max = out_w$num_paths_to_max,
##                      num_accessible_genots = out_w$num_accessible_genots,
##                      diversity_weighted_paths_to_max = out_w$diversity_weighted_paths_to_max)
                     
##     m1$replicate <- outn$iter
##     likely_error <- as.data.frame(likely_error)
##     paths_error <- as.data.frame(paths_error)
    
##     m1 <- cbind(m1, likely_error)
##     m1 <- cbind(m1, paths_error)
    
##     rm(theout)
##     dfoutn <- as.data.frame(outn)
##     row.names(dfoutn) <- NULL
    
##     retout <- cbind(dfoutn, m1)
##     return(retout)
## }




## a simple check
any_constant_col <- function(x) {
    nr <- nrow(x)
    mcs <- max(colSums(x))
    any(mcs == nr)
}

## convert data frame to a matrix with 0L and 1L
df_2_mat_integer <- function(x) {
    x1 <- as.matrix(x)
    if(max(abs(x - x1)) != 0) stop("failed conversion to matrix")
    x2 <- x1
    storage.mode(x2) <- "integer"
    if(max(abs(x2 - x1)) != 0) stop("Not in 0L, 1L")
    if(max(abs(x - x2)) != 0) stop("df not in 0L, 1L") ## paranoia
    return(x2)
}

## To make it explicit
## but do not set the last row to NaNs or similar.
## Following same logic as in trans_rate_to_trans_mat (in MHN dir)
rowScaleMatrix <- function(x) {
    tm <- x
    sx <- rowSums(x)
    ii <- which(sx > 0)
    for(i in ii) {
        tm[i, ] <- tm[i, ]/sx[i]
    }
    tm
}



######################################################################
######################################################################

####################### Tests #################

## Not using a test_that block here. But doing it locally
local({


ex0 <- list(edges = data.frame(
                From = c("Root", "B", "C", "D"),
                To =   c("B", "C", "D", "A"),
                rerun_lambda = c(1, 2, 3, 4)
            ))

    
ex1 <- list(edges = data.frame(
                From = c("Root", "Root", "A", "B", "C"),
                To =   c("A", "B", "C", "C", "D"),
                rerun_lambda = c(1, 2, 3, 3, 4)
            ))

ex2 <- list(edges = data.frame(
                From = c("Root", "Root", "A", "B"),
                To =   c("A", "B", "C", "D"),
                rerun_lambda = c(10, 11, 12, 14)
            ))

ex3 <- list(edges = data.frame(
                From = c("Root", "Root", "A", "C"),
                To =   c("A", "B", "C", "D"),
                rerun_lambda = c(2, 3, 4, 5)
            ))

ex4 <- list(edges = data.frame(
                From = c("Root", "Root", "A", "A"),
                To =   c("A", "B", "C", "D"),
                rerun_lambda = c(1, 2, 3, 4)
            ))

ex5 <- list(edges = data.frame(
                From = c("Root", "Root", "Root", "A", "B", "C"),
                To =   c("A", "B", "C", "D", "D", "D"),
                rerun_lambda = c(1, 2, 3, 4, 4, 4)
            ))
## yes. this is wrong. Tested below.
ex5e <- list(edges = data.frame(
                From = c("Root", "Root", "Root", "A", "B", "C"),
                To =   c("A", "B", "C", "D", "D", "D"),
                rerun_lambda = c(1, 2, 3, 4, 4, 5)
            ))

ex6 <- list(edges = data.frame(
                From = c("Root", "Root", "Root", "A", "B"),
                To =   c("A", "B", "D", "C", "C"),
                rerun_lambda = c(1, 2, 3, 5, 5)
            ))

ex7 <- list(edges = data.frame(
                From = c("Root", "E", "E", "E", "A", "B"),
                To =   c("E",    "A", "B", "D", "C", "C"),
                rerun_lambda = c(1, 2, 3, 4, 5, 5)
            ))

ex8 <- list(edges = data.frame(
                From = c("Root", "Root", "A", "B", "C", "C"),
                To =   c("A",     "B",   "C", "C", "D",  "E"),
                rerun_lambda = c(1, 2, 3, 3, 4, 5)
            ))

ex9 <- list(edges = data.frame(
                From = c("Root", "A", "B", "B", "B"),
                To =   c("A",    "B", "C", "D", "E"),
                rerun_lambda = c(6, 7, 8, 9, 10)
            ))
## wrong, triggers error
ex10e <- list(edges = data.frame(
                From = c("Root", "A", "A", "B", "C", "D"),
                To =   c("A",    "B", "C", "D", "D", "E"),
                rerun_lambda = c(1, 2, 3, 4, 5, 5)
              ))

ex10 <- list(edges = data.frame(
                From = c("Root", "A", "A", "B", "C", "D"),
                To =   c("A",    "B", "C", "D", "D", "E"),
                rerun_lambda = c(1, 2, 3, 4, 4, 5)
            ))


ex11 <- list(edges = data.frame(
                From = c(rep("Root", 4), "A", "B", "C", "D"),
                To =   c("A",  "B", "C", "D", rep("E", 4)),
                rerun_lambda = c(3, 4, 5, 6, rep(7, 4))
            ))



oex0 <- cpm_access_genots_paths_w(ex0)
oex1 <- cpm_access_genots_paths_w(ex1)
oex2 <- cpm_access_genots_paths_w(ex2)
oex3 <- cpm_access_genots_paths_w(ex3)
oex4 <- cpm_access_genots_paths_w(ex4)
oex5 <- cpm_access_genots_paths_w(ex5)
expect_error(cpm_access_genots_paths_w(ex5e), "Different lambda", fixed = TRUE)

oex6 <- cpm_access_genots_paths_w(ex6)
oex7 <- cpm_access_genots_paths_w(ex7)
oex8 <- cpm_access_genots_paths_w(ex8)
oex9 <- cpm_access_genots_paths_w(ex9)
expect_error(cpm_access_genots_paths_w(ex10e), "Different lambda", fixed = TRUE)
oex10 <- cpm_access_genots_paths_w(ex10)
oex11 <- cpm_access_genots_paths_w(ex11)


expect_equivalent(oex0$weighted_paths[, 2], 1)

expect_equivalent(oex1$weighted_paths[, 2],
                  c(1/3 * 1 * 1 * 1, 2/3 * 1 * 1 * 1))

expect_equivalent(oex2$weighted_paths[, 2],
                  c(10/21 * 11/23 * 12/26 * 1,
                    10/21 * 11/23 * 14/26 * 1,
                    10/21 * 12/23 * 1     * 1,
                    11/21 * 10/24 * 12/26 * 1,
                    11/21 * 10/24 * 14/26 * 1,
                    11/21 * 14/24 * 1     * 1))

expect_equivalent(oex3$weighted_paths[, 2],
                  c(2/5 * 3/7 * 1 * 1,
                    2/5 * 4/7 * 3/8 * 1,
                    2/5 * 4/7 * 5/8 * 1,
                    3/5 * 1 *   1   * 1
                  ))

expect_equivalent(oex4$weighted_paths[, 2],
                  c(1/3 * 2/9 * 3/7 * 1,
                    1/3 * 2/9 * 4/7 * 1,
                    1/3 * 3/9 * 2/6 * 1,
                    1/3 * 3/9 * 4/6 * 1,
                    1/3 * 4/9 * 2/5 * 1,
                    1/3 * 4/9 * 3/5 * 1,
                    2/3 * 1 *   3/7 * 1,
                    2/3 * 1 *   4/7 * 1
                  ))

expect_equivalent(oex5$weighted_paths[, 2],
                  c(1/6 * 2/5 * 1 * 1,
                    1/6 * 3/5 * 1 * 1,
                    2/6 * 1/4 * 1 * 1,
                    2/6 * 3/4 * 1 * 1,
                    3/6 * 1/3 * 1 * 1,
                    3/6 * 2/3 * 1 * 1
                    ))

expect_equivalent(oex6$weighted_paths[, 2],
                    c(1/6 * 2/5 * 3/8 * 1,
                      1/6 * 2/5 * 5/8 * 1,
                      1/6 * 3/5 * 1 * 1,
                      2/6 * 1/4 * 3/8 * 1,
                      2/6 * 1/4 * 5/8 * 1,
                      2/6 * 3/4 * 1 * 1,
                      3/6 * 1/3 * 1 * 1,
                      3/6 * 2/3 * 1 * 1
                    ))


expect_equivalent(oex7$weighted_paths[, 2],
                  c(1 * 2/9 * 3/7 * 4/9,
                    1 * 2/9 * 3/7 * 5/9,
                    1 * 2/9 * 4/7 * 1,
                    1 * 3/9 * 2/6 * 4/9,
                    1 * 3/9 * 2/6 * 5/9,
                    1 * 3/9 * 4/6 * 1,
                    1 * 4/9 * 2/5 * 1,
                    1 * 4/9 * 3/5 * 1
                  ))
                      
expect_equivalent(oex8$weighted_paths[, 2],
                  c(1/3 * 1 * 1 * 4/9,
                    1/3 * 1 * 1 * 5/9,
                    2/3 * 1 * 1 * 4/9,
                    2/3 * 1 * 1 * 5/9
                  ))
                      
 
expect_equivalent(oex9$weighted_paths[, 2],
                  c(6/6 * 7/7 * 8/27 * 9/19 * 10/10,
                    6/6 * 7/7 * 8/27 * 10/19 * 9/9,                    
                    6/6 * 7/7 * 9/27 * 8/18 * 10/10,
                    6/6 * 7/7 * 9/27 * 10/18 * 8/8,                   
                    6/6 * 7/7 * 10/27 * 8/17 * 9/9,
                     6/6 * 7/7 * 10/27 * 9/17 * 8/8
                    ))

expect_equivalent(oex10$weighted_paths[, 2],
                  c(1/1 * 2/5 * 3/3 * 4/4 * 5/5,
                    1/1 * 3/5 * 2/2 * 4/4 * 5/5                                        
                    ))

expect_equivalent(oex11$weighted_paths[, 2],
                  c(3/18 * 4/15 * 5/11 * 7/7,
                    3/18 * 4/15 * 6/11 * 7/7,
                    3/18 * 5/15 * 4/10 * 7/7,
                    3/18 * 5/15 * 6/10 * 7/7,
                    3/18 * 6/15 * 4/9 * 7/7,
                    3/18 * 6/15 * 5/9 * 7/7,
                    4/18 * 3/14 * 5/11 * 7/7,                    
                    4/18 * 3/14 * 6/11 * 7/7,
                    4/18 * 5/14 * 3/9 * 7/7,
                    4/18 * 5/14 * 6/9 * 7/7,
                    4/18 * 6/14 * 3/8 * 7/7,
                    4/18 * 6/14 * 5/8 * 7/7,                    
                    5/18 * 3/13 * 4/10 * 7/7,
                    5/18 * 3/13 * 6/10 * 7/7,
                    5/18 * 4/13 * 3/9 * 7/7,
                    5/18 * 4/13 * 6/9 * 7/7,
                    5/18 * 6/13 * 3/7 * 7/7,
                    5/18 * 6/13 * 4/7 * 7/7,                    
                    6/18 * 3/12 * 4/9 * 7/7,
                    6/18 * 3/12 * 5/9 * 7/7,
                    6/18 * 4/12 * 3/8 * 7/7,
                    6/18 * 4/12 * 5/8 * 7/7,
                    6/18 * 5/12 * 3/7 * 7/7,
                    6/18 * 5/12 * 4/7 * 7/7                    
                    ))




## Test the simplified code
oex0_simplified <- cpm_access_genots_paths_w_simplified(ex0)
oex1_simplified <- cpm_access_genots_paths_w_simplified(ex1)
oex2_simplified <- cpm_access_genots_paths_w_simplified(ex2)
oex3_simplified <- cpm_access_genots_paths_w_simplified(ex3)
oex4_simplified <- cpm_access_genots_paths_w_simplified(ex4)
oex5_simplified <- cpm_access_genots_paths_w_simplified(ex5)
expect_error(cpm_access_genots_paths_w_simplified(ex5e), "Different lambda", fixed = TRUE)

## We are now using sparseMatrices: fuller testing component by component below

oex6_simplified <- cpm_access_genots_paths_w_simplified(ex6)
oex7_simplified <- cpm_access_genots_paths_w_simplified(ex7)
oex8_simplified <- cpm_access_genots_paths_w_simplified(ex8)
oex9_simplified <- cpm_access_genots_paths_w_simplified(ex9)
expect_error(cpm_access_genots_paths_w_simplified(ex10e), "Different lambda", fixed = TRUE)
oex10_simplified <- cpm_access_genots_paths_w_simplified(ex10)
oex11_simplified <- cpm_access_genots_paths_w_simplified(ex11)

expect_equal(oex0[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
             lapply(oex0_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))

expect_equal(oex1[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
             lapply(oex1_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))

expect_equal(oex2[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
             lapply(oex2_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))

expect_equal(oex3[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
             lapply(oex3_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))

expect_equal(oex4[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
             lapply(oex4_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))

expect_equal(oex5[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
             lapply(oex5_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))

expect_equal(oex6[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
             lapply(oex6_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))

expect_equal(oex7[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
             lapply(oex7_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))

expect_equal(oex8[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
             lapply(oex8_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))

expect_equal(oex9[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
             lapply(oex9_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))

expect_equal(oex10[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
             lapply(oex10_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))

expect_equal(oex11[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
             lapply(oex11_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))

gacc1 <- list("A", "B", "D", c("A", "B"), c("A", "D"), c("B", "D"), c("A", 
"C"), c("A", "B", "D"), c("A", "B", "C"), c("A", "C", "D"), c("A", 
"B", "C", "D"))

am1 <- unrestricted_fitness_graph(gacc1)
am1M <- unrestricted_fitness_graph_sparseM(gacc1)
am1Mm <- as.matrix(am1M)
storage.mode(am1Mm) <- "integer"
expect_identical(am1, am1Mm)

wg <- structure(list(To = c("A", "B", "C", "D"), OT_edgeWeight = c(0.525915054637741, 
0.101508072999909, 0.108065026764796, 0.302542959038882)), row.names = c("A", 
"B", "C", "D"), class = "data.frame")

tm1 <- transition_fg(am1, wg)
tm1M <- transition_fg_sparseM(am1M, wg)
tm1M <- as.matrix(tm1M)
expect_identical(tm1, tm1M)

## To test with a full run. Checking sparse matrix implementation

cpm_out_others1 <- list(OT = list(edges = structure(list(From = c("Root", "Root", 
"A", "Root"), To = c("A", "B", "C", "D"), edge = c("Root -> A", 
"Root -> B", "A -> C", "Root -> D"), OT_edgeBootFreq = c(NA, 
NA, NA, NA), OT_edgeWeight = c(0.525915054637741, 0.101508072999909, 
0.108065026764796, 0.302542959038882), OT_obsMarginal = c(0.56, 
0.18, 0.14, 0.36), OT_predMarginal = c(0.56, 0.18, 0.139999436433903, 
0.36)), class = "data.frame", row.names = c("A", "B", "C", "D"
)), consensus = NA, OT_error.fun = "std", ot.boot.original = NA, 
    genots_predicted = NA, genots_observed = NA, two_way_predicted = NA, 
    two_way_observed = NA), CAPRESE = list(edges = structure(list(
    From = c("Root", "Root", "Root", "A"), To = c("A", "B", "D", 
    "C"), edge = c("Root -> A", "Root -> B", "Root -> D", "A -> C"
    ), CAPRESE_edgeBootFreq = c(NA, NA, NA, NA), CAPRESE_pr = c(NA, 
    NA, NA, 0.223188719001567), CAPRESE_tp = c(NA, NA, NA, 1), 
    CAPRESE_hg = c(NA, NA, NA, 0.0948328267477203)), class = "data.frame", row.names = c(NA, 
-4L)), CAPRESE_eloss = NA, CAPRESE_prederr = NA, CAPRESE_posterr = NA, 
    CAPRESE_useless_extra = list(selective_advantage = list(caprese = structure(list(
        SELECTS = "variant A", SELECTED = "variant C", OBS.SELECTS = 28, 
        OBS.SELECTED = 7, TEMPORAL.PRIORITY = 1, PROBABILITY.RAISING = 0.223188719001567, 
        HYPERGEOMETRIC = 0.0948328267477203), row.names = "1", class = "data.frame")), 
        bootstrap_scores = NA), TRONCO_model_boot = NULL, TRONCO_model = list(
        genotypes = structure(c(0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 
        1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 
        0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 
        0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 
        0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 
        0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 
        0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 
        0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 
        0L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 
        1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L), .Dim = c(50L, 
        4L), .Dimnames = list(c("1", "2", "3", "4", "5", "6", 
        "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", 
        "17", "18", "19", "20", "21", "22", "23", "24", "25", 
        "26", "27", "28", "29", "30", "31", "32", "33", "34", 
        "35", "36", "37", "38", "39", "40", "41", "42", "43", 
        "44", "45", "46", "47", "48", "49", "50"), c("G1", "G2", 
        "G3", "G4"))), annotations = structure(c("variant", "variant", 
        "variant", "variant", "A", "B", "C", "D"), .Dim = c(4L, 
        2L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("type", 
        "event"))), types = structure("Darkgreen", .Dim = c(1L, 
        1L), .Dimnames = list("variant", "color")), hypotheses = NA, 
        confidence = structure(list(structure(c(0, 0.321428571428571, 
        0.25, 0.642857142857143, 1, 0, 0.777777777777778, 1, 
        1, 1, 0, 1, 1, 0.5, 0.388888888888889, 0), .Dim = c(4L, 
        4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("G1", 
        "G2", "G3", "G4"))), structure(c(-1, -0.125943176525384, 
        0.132346627230226, -0.00509762193462174, -0.171557211613104, 
        -1, -0.123178689550371, -0.281478505423865, 0.223188719001567, 
        -0.125943176525384, -1, -0.14193770830939, -0.00649653637701445, 
        -0.255488752279943, -0.123178689550371, -1), .Dim = c(4L, 
        4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("G1", 
        "G2", "G3", "G4"))), structure(c(0, 0.657728939304999, 
        0.0948328267477203, 0.403066884258277, 0.657728939304999, 
        0, 0.369777142376588, 0.707693654795176, 0.0948328267477203, 
        0.369777142376588, 0, 0.494537285101578, 0.403066884258277, 
        0.707693654795176, 0.494537285101578, 0), .Dim = c(4L, 
        4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("G1", 
        "G2", "G3", "G4")))), .Dim = c(3L, 1L), .Dimnames = list(
            c("temporal priority", "probability raising", "hypergeometric test"
            ), "confidence")), model = list(caprese = list(probabilities = list(
            probabilities.observed = list(marginal.probs = structure(c(0.56, 
            0.18, 0.14, 0.36), .Dim = c(4L, 1L), .Dimnames = list(
                c("G1", "G2", "G3", "G4"), "marginal probability")), 
                joint.probs = structure(c(0.56, 0.08, 0.1, 0.2, 
                0.08, 0.18, 0.02, 0.04, 0.1, 0.02, 0.14, 0.04, 
                0.2, 0.04, 0.04, 0.36), .Dim = c(4L, 4L), .Dimnames = list(
                  c("G1", "G2", "G3", "G4"), c("G1", "G2", "G3", 
                  "G4"))), conditional.probs = structure(c(1, 
                1, 0.178571428571429, 1), .Dim = c(4L, 1L), .Dimnames = list(
                  c("G1", "G2", "G3", "G4"), "conditional probability"))), 
            probabilities.fit = list(estimated.marginal.probs = NA, 
                estimated.joint.probs = NA, estimated.conditional.probs = NA)), 
            parents.pos = structure(c(-1, -1, 1, -1), .Dim = c(4L, 
            1L), .Dimnames = list(c("G1", "G2", "G3", "G4"), 
                "parents")), error.rates = list(error.fp = NA, 
                error.fn = NA), adj.matrix = list(adj.matrix.fit = structure(c(0, 
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), .Dim = c(4L, 
            4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), 
                c("G1", "G2", "G3", "G4")))), logLik = -110.377237990856)), 
        parameters = list(algorithm = "CAPRESE", lambda = 0.5, 
            silent = TRUE, error.rates = list(epos = 0, eneg = 0)), 
        execution.time = structure(c(user.self = 0.00199999999998113, 
        sys.self = 0, elapsed = 0.00100000000020373, user.child = 0, 
        sys.child = 0), class = "proc_time"))), CAPRI_BIC = list(
    edges = structure(list(From = c("Root", "Root", "Root", "Root"
    ), To = c("A", "B", "C", "D"), edge = c("Root -> A", "Root -> B", 
    "Root -> C", "Root -> D"), CAPRI_edgeBootFreq = c(NA, NA, 
    NA, NA), CAPRI_pr = c(NA_real_, NA_real_, NA_real_, NA_real_
    ), CAPRI_tp = c(NA_real_, NA_real_, NA_real_, NA_real_), 
        CAPRI_hg = c(NA_real_, NA_real_, NA_real_, NA_real_)), class = "data.frame", row.names = c(NA, 
    -4L)), CAPRI_eloss = NA, CAPRI_prederr = NA, CAPRI_posterr = NA, 
    CAPRI_useless_extra = list(selective_advantage = list(capri_bic = NULL), 
        bootstrap_scores = NA), TRONCO_model_boot = NULL, TRONCO_model = list(
        genotypes = structure(c(0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 
        1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 
        0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 
        0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 
        0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 
        0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 
        0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 
        0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 
        0L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 
        1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L), .Dim = c(50L, 
        4L), .Dimnames = list(c("1", "2", "3", "4", "5", "6", 
        "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", 
        "17", "18", "19", "20", "21", "22", "23", "24", "25", 
        "26", "27", "28", "29", "30", "31", "32", "33", "34", 
        "35", "36", "37", "38", "39", "40", "41", "42", "43", 
        "44", "45", "46", "47", "48", "49", "50"), c("G1", "G2", 
        "G3", "G4"))), annotations = structure(c("variant", "variant", 
        "variant", "variant", "A", "B", "C", "D"), .Dim = c(4L, 
        2L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("type", 
        "event"))), types = structure("Darkgreen", .Dim = c(1L, 
        1L), .Dimnames = list("variant", "color")), hypotheses = NA, 
        adj.matrix.prima.facie = structure(c(0, 0, 0, 0, 0, 0, 
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0), .Dim = c(4L, 4L), .Dimnames = list(
            c("G1", "G2", "G3", "G4"), c("G1", "G2", "G3", "G4"
            ))), confidence = structure(list(structure(c(1, 1, 
        1, 1, 9.39901768425988e-35, 1, 0.999999998490081, 1.54941877188073e-32, 
        9.94334396623681e-35, 1.53273973984097e-09, 1, 2.65146581110054e-34, 
        4.00518485843201e-31, 1, 1, 1), .Dim = c(4L, 4L), .Dimnames = list(
            c("G1", "G2", "G3", "G4"), c("G1", "G2", "G3", "G4"
            ))), structure(c(1, 0.999999999923137, 3.71813560969053e-15, 
        0.846484514790407, 0.999999999999868, 1, 0.999886735833476, 
        1, 8.73588365187455e-17, 0.999958011279254, 1, 0.999997831094051, 
        0.781442879272848, 1, 0.999999534542517, 1), .Dim = c(4L, 
        4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("G1", 
        "G2", "G3", "G4"))), structure(c(1, 0.872936649650717, 
        0.0948328267477203, 0.635507694679078, 0.872936649650717, 
        1, 0.369777142376588, 0.707693654795176, 0.0948328267477203, 
        0.369777142376588, 1, 0.494537285101578, 0.635507694679078, 
        0.707693654795176, 0.494537285101578, 1), .Dim = c(4L, 
        4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("G1", 
        "G2", "G3", "G4")))), .Dim = c(3L, 1L), .Dimnames = list(
            c("temporal priority", "probability raising", "hypergeometric test"
            ), "confidence")), model = list(capri_bic = list(
            probabilities = list(probabilities.observed = list(
                marginal.probs = structure(c(0.5534, 0.183, 0.1418, 
                0.3628), .Dim = c(4L, 1L), .Dimnames = list(c("G1", 
                "G2", "G3", "G4"), "marginal probability")), 
                joint.probs = structure(c(0.5534, 0.0776, 0.1008, 
                0.1978, 0.0776, 0.183, 0.0204, 0.0422, 0.1008, 
                0.0204, 0.1418, 0.041, 0.1978, 0.0422, 0.041, 
                0.3628), .Dim = c(4L, 4L), .Dimnames = list(c("G1", 
                "G2", "G3", "G4"), c("G1", "G2", "G3", "G4"))), 
                conditional.probs = structure(list(1, 1, 1, 1), .Dim = c(4L, 
                1L), .Dimnames = list(c("G1", "G2", "G3", "G4"
                ), "conditional probability"))), probabilities.fit = list(
                estimated.marginal.probs = NA, estimated.joint.probs = NA, 
                estimated.conditional.probs = NA)), parents.pos = structure(list(
                -1, -1, -1, -1), .Dim = c(4L, 1L), .Dimnames = list(
                c("G1", "G2", "G3", "G4"), "parents")), error.rates = list(
                error.fp = NA, error.fn = NA), adj.matrix = list(
                adj.matrix.pf = structure(c(0, 0, 0, 0, 0, 0, 
                0, 0, 1, 0, 0, 0, 0, 0, 0, 0), .Dim = c(4L, 4L
                ), .Dimnames = list(c("G1", "G2", "G3", "G4"), 
                  c("G1", "G2", "G3", "G4"))), adj.matrix.fit = structure(c(0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), .Dim = c(4L, 
                4L), .Dimnames = list(c("G1", "G2", "G3", "G4"
                ), c("G1", "G2", "G3", "G4")))), score = -118.609294356862, 
            logLik = -110.785248346005)), parameters = list(algorithm = "CAPRI", 
            command = "hc", regularization = "bic", do.boot = TRUE, 
            nboot = 100, pvalue = 0.05, min.boot = 3, min.stat = TRUE, 
            boot.seed = NULL, silent = TRUE, error.rates = list(
                epos = 0, eneg = 0), restart = 100), execution.time = structure(c(user.self = 0.0769999999999982, 
        sys.self = 0, elapsed = 0.077000000001135, user.child = 0, 
        sys.child = 0), class = "proc_time"))), CAPRI_AIC = list(
    edges = structure(list(From = c("Root", "Root", "Root", "Root"
    ), To = c("A", "B", "C", "D"), edge = c("Root -> A", "Root -> B", 
    "Root -> C", "Root -> D"), CAPRI_edgeBootFreq = c(NA, NA, 
    NA, NA), CAPRI_pr = c(NA_real_, NA_real_, NA_real_, NA_real_
    ), CAPRI_tp = c(NA_real_, NA_real_, NA_real_, NA_real_), 
        CAPRI_hg = c(NA_real_, NA_real_, NA_real_, NA_real_)), class = "data.frame", row.names = c(NA, 
    -4L)), CAPRI_eloss = NA, CAPRI_prederr = NA, CAPRI_posterr = NA, 
    CAPRI_useless_extra = list(selective_advantage = list(capri_aic = NULL), 
        bootstrap_scores = NA), TRONCO_model_boot = NULL, TRONCO_model = list(
        genotypes = structure(c(0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 
        1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 
        0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 
        0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 
        0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 
        0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 
        0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 
        0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 
        0L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 
        1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L), .Dim = c(50L, 
        4L), .Dimnames = list(c("1", "2", "3", "4", "5", "6", 
        "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", 
        "17", "18", "19", "20", "21", "22", "23", "24", "25", 
        "26", "27", "28", "29", "30", "31", "32", "33", "34", 
        "35", "36", "37", "38", "39", "40", "41", "42", "43", 
        "44", "45", "46", "47", "48", "49", "50"), c("G1", "G2", 
        "G3", "G4"))), annotations = structure(c("variant", "variant", 
        "variant", "variant", "A", "B", "C", "D"), .Dim = c(4L, 
        2L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("type", 
        "event"))), types = structure("Darkgreen", .Dim = c(1L, 
        1L), .Dimnames = list("variant", "color")), hypotheses = NA, 
        adj.matrix.prima.facie = structure(c(0, 0, 0, 0, 0, 0, 
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0), .Dim = c(4L, 4L), .Dimnames = list(
            c("G1", "G2", "G3", "G4"), c("G1", "G2", "G3", "G4"
            ))), confidence = structure(list(structure(c(1, 1, 
        1, 1, 9.81725428393167e-35, 1, 0.999990867323265, 8.96819397729642e-32, 
        9.74566638078606e-35, 9.23440871142156e-06, 1, 3.31407732958239e-34, 
        1.34669152568532e-32, 1, 1, 1), .Dim = c(4L, 4L), .Dimnames = list(
            c("G1", "G2", "G3", "G4"), c("G1", "G2", "G3", "G4"
            ))), structure(c(1, 0.999999998167757, 2.73780178422517e-19, 
        0.463492715389974, 0.999999986477575, 1, 0.998529402441165, 
        1, 5.49487631358828e-17, 0.99916699252184, 1, 0.999998976978537, 
        0.525815417121827, 1, 0.999999778739998, 1), .Dim = c(4L, 
        4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("G1", 
        "G2", "G3", "G4"))), structure(c(1, 0.657728939304999, 
        0.0948328267477203, 0.506373911111967, 0.657728939304999, 
        1, 0.369777142376588, 0.890399114532546, 0.0948328267477203, 
        0.369777142376588, 1, 0.445575084798027, 0.506373911111967, 
        0.890399114532546, 0.445575084798027, 1), .Dim = c(4L, 
        4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("G1", 
        "G2", "G3", "G4")))), .Dim = c(3L, 1L), .Dimnames = list(
            c("temporal priority", "probability raising", "hypergeometric test"
            ), "confidence")), model = list(capri_aic = list(
            probabilities = list(probabilities.observed = list(
                marginal.probs = structure(c(0.5678, 0.1802, 
                0.1472, 0.3494), .Dim = c(4L, 1L), .Dimnames = list(
                  c("G1", "G2", "G3", "G4"), "marginal probability")), 
                joint.probs = structure(c(0.5678, 0.0864, 0.1066, 
                0.1978, 0.0864, 0.1802, 0.0222, 0.0374, 0.1066, 
                0.0222, 0.1472, 0.0406, 0.1978, 0.0374, 0.0406, 
                0.3494), .Dim = c(4L, 4L), .Dimnames = list(c("G1", 
                "G2", "G3", "G4"), c("G1", "G2", "G3", "G4"))), 
                conditional.probs = structure(list(1, 1, 1, 1), .Dim = c(4L, 
                1L), .Dimnames = list(c("G1", "G2", "G3", "G4"
                ), "conditional probability"))), probabilities.fit = list(
                estimated.marginal.probs = NA, estimated.joint.probs = NA, 
                estimated.conditional.probs = NA)), parents.pos = structure(list(
                -1, -1, -1, -1), .Dim = c(4L, 1L), .Dimnames = list(
                c("G1", "G2", "G3", "G4"), "parents")), error.rates = list(
                error.fp = NA, error.fn = NA), adj.matrix = list(
                adj.matrix.pf = structure(c(0, 0, 0, 0, 0, 0, 
                0, 0, 1, 0, 0, 0, 0, 0, 0, 0), .Dim = c(4L, 4L
                ), .Dimnames = list(c("G1", "G2", "G3", "G4"), 
                  c("G1", "G2", "G3", "G4"))), adj.matrix.fit = structure(c(0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), .Dim = c(4L, 
                4L), .Dimnames = list(c("G1", "G2", "G3", "G4"
                ), c("G1", "G2", "G3", "G4")))), score = -114.785248346005, 
            logLik = -110.785248346005)), parameters = list(algorithm = "CAPRI", 
            command = "hc", regularization = "aic", do.boot = TRUE, 
            nboot = 100, pvalue = 0.05, min.boot = 3, min.stat = TRUE, 
            boot.seed = NULL, silent = TRUE, error.rates = list(
                epos = 0, eneg = 0), restart = 100), execution.time = structure(c(user.self = 0.0769999999999982, 
        sys.self = 0, elapsed = 0.0769999999993161, user.child = 0, 
        sys.child = 0), class = "proc_time"))), CBN_ot = list(
    edges = structure(list(From = c("Root", "C", "D", "A", "Root"
    ), To = c("A", "B", "B", "C", "D"), edge = c("Root -> A", 
    "C -> B", "D -> B", "A -> C", "Root -> D"), init_lambda = c(1.414067, 
    0.009031, 0.009031, 0.012027, 0.381915), final_lambda = c(1.416377, 
    0.109282, 0.109282, 0.014903, 0.380714), rerun_lambda = c(1.416373, 
    0.110483, 0.110483, 0.01489, 0.380711), CBN_edgeBootFreq = c(NA, 
    NA, NA, NA, NA)), class = "data.frame", row.names = c(NA, 
    -5L)), nboot = 0, init.poset = "OT"), MCCBN = list(edges = structure(list(
    From = c("Root", "Root", "Root", "A"), To = c("A", "B", "D", 
    "C"), edge = c("Root -> A", "Root -> B", "Root -> D", "A -> C"
    ), lambda = c(1.17159022808865, 0.217522208241225, 0.460587464855744, 
    0.25995329754024)), row.names = c(NA, -4L), class = "data.frame")), 
    time_ot = c(elapsed = 0.00399999999899592), time_caprese = c(elapsed = 0.011000000000422), 
    time_capri_bic = c(elapsed = 0.0860000000011496), time_capri_aic = c(elapsed = 0.0869999999995343), 
    time_cbn_ot = c(elapsed = 0.199999999998909), time_mccbn = c(elapsed = 0.514999999999418), 
    input_data = structure(c(0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 0L, 
    0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 
    0L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 
    1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 
    0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 
    1L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 
    0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L), .Dim = c(50L, 
    4L), .Dimnames = list(c("1", "2", "3", "4", "5", "6", "7", 
    "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", 
    "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", 
    "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", 
    "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", 
    "48", "49", "50"), c("A", "B", "C", "D"))), input_data_pseudosamples = structure(c(0L, 
    0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 
    0L, 0L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 
    1L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 
    0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 
    0L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 
    0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 
    1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 
    0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 
    0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 
    0L, 1L, 0L, 0L), .Dim = c(50L, 4L), .Dimnames = list(c("1", 
    "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", 
    "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", 
    "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", 
    "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", 
    "43", "44", "45", "46", "47", "48", "49", "50"), c("A", "B", 
    "C", "D"))))

mm <- c("OT", "CAPRESE", "CAPRI_BIC", "CAPRI_AIC",
                 "CBN_ot", "MCCBN")

mm0  <- lapply(cpm_out_others1[mm], cpm_access_genots_paths_w)
mmSM <- lapply(cpm_out_others1[mm], cpm_access_genots_paths_w_simplified)

## Identical unweighted transition matrices (fitness graphs)
uw0 <- lapply(mm0, function(x) rowScaleMatrix(x$fgraph))
uwSM <- lapply(mmSM, function(x) rowScaleMatrix(x$fgraph))
uwSMm <- lapply(uwSM, as.matrix)
expect_identical(uw0, uwSMm)

## Identical Weighted transition matrices
wg0 <- lapply(mm0[c("OT", "MCCBN", "CBN_ot")],
              function(x) x$trans_mat_genots)
wgSM <- lapply(mmSM[c("OT", "MCCBN", "CBN_ot")],
              function(x) x$trans_mat_genots)
wgSMm <- lapply(wgSM, as.matrix)
expect_identical(wg0, wgSMm)

## Diagonal
td0 <- lapply(mm0[c("MCCBN", "CBN_ot")],
              function(x)
                  trans_rate_to_trans_mat(x$weighted_fgraph,
                                          method = "uniformization")) 
tdSM <- lapply(mmSM[c("MCCBN", "CBN_ot")],
              function(x)
                  trans_rate_to_trans_mat(x$weighted_fgraph,
                                          method = "uniformization")) 
tdSMm <- lapply(tdSM, as.matrix)
expect_identical(td0, tdSMm)

## recheck competing exponentials done differently
wg2 <- lapply(mm0[c("OT", "MCCBN", "CBN_ot")],
                  function(x) trans_rate_to_trans_mat(x$weighted_fgraph,
                                                      method = "competingExponentials"))
wg2SM <- lapply(mmSM[c("OT", "MCCBN", "CBN_ot")],
                  function(x) as.matrix(trans_rate_to_trans_mat(x$weighted_fgraph,
                                                      method = "competingExponentials")))
expect_identical(wg2, wg2SM)
## tripwire to check tests run
## expect_identical(1, 3)

rm(mm, mm0, mmSM)
rm(cpm_out_others1)

test_others <- function(data) {
    data <- as.matrix(data)
    data <- df_2_mat_integer(data)
    cpm_out_others2 <- all_methods(data)
    mm <- c("OT", "CAPRESE", "CAPRI_BIC", "CAPRI_AIC",
            "CBN_ot", "MCCBN")

    mm0  <- lapply(cpm_out_others2[mm], cpm_access_genots_paths_w)
    mmSM <- lapply(cpm_out_others2[mm], cpm_access_genots_paths_w_simplified)

    ## Identical unweighted transition matrices (fitness graphs)
    uw0 <- lapply(mm0, function(x) rowScaleMatrix(x$fgraph))
    uwSM <- lapply(mmSM, function(x) rowScaleMatrix(x$fgraph))
    uwSMm <- lapply(uwSM, as.matrix)
    expect_identical(uw0, uwSMm)

    ## Identical Weighted transition matrices
    wg0 <- lapply(mm0[c("OT", "MCCBN", "CBN_ot")],
                  function(x) x$trans_mat_genots)
    wgSM <- lapply(mmSM[c("OT", "MCCBN", "CBN_ot")],
                   function(x) x$trans_mat_genots)
    wgSMm <- lapply(wgSM, as.matrix)
    expect_equal(wg0, wgSMm)

    ## Diagonal
    td0 <- lapply(mm0[c("MCCBN", "CBN_ot")],
                  function(x)
                      trans_rate_to_trans_mat(x$weighted_fgraph,
                                              method = "uniformization")) 
    tdSM <- lapply(mmSM[c("MCCBN", "CBN_ot")],
                   function(x)
                       trans_rate_to_trans_mat(x$weighted_fgraph,
                                               method = "uniformization")) 
    tdSMm <- lapply(tdSM, as.matrix)
    expect_equal(td0, tdSMm)

    ## recheck competing exponentials done differently
    wg2 <- lapply(mm0[c("OT", "MCCBN", "CBN_ot")],
                  function(x) trans_rate_to_trans_mat(x$weighted_fgraph,
                                                      method = "competingExponentials"))
    wg2SM <- lapply(mmSM[c("OT", "MCCBN", "CBN_ot")],
                    function(x) as.matrix(trans_rate_to_trans_mat(x$weighted_fgraph,
                                                                  method = "competingExponentials")))
    expect_equal(wg2, wg2SM)
    
}

## Additional testing


Dat1 <- readRDS(file="./MHN/data/BreastCancer.rds") [1:50, 1:4]
test_others(Dat1)

Dat1 <- readRDS(file="./MHN/data/ColorectalCancer.rds")[1:40, 1:6]
test_others(Dat1)

Dat1 <- readRDS(file="./MHN/data/RenalCellCarcinoma.rds")[1:30,2:6]
test_others(Dat1)

Dat1 <- readRDS(file="./MHN/data/Glioblastoma.rds")[1:20, 1:5]
test_others(Dat1)

rm(Dat1)

}) ## end of testing



## Pass a data set as a matrix with subjects as rows and genes as columns

all_methods_2_trans_mat <- function(x, cores_cbn = 1) {
    x <- df_2_mat_integer(x)

    cat("\n  Number of genes before limiting = ", ncol(x))
    ## Always limit to 15
    x <- pre_process(x, remove.constant = FALSE,
                     min.freq = 0, max.cols = 15)
    
    cat("\n  Number of genes after limiting = ", ncol(x), "\n")

    methods <- c("OT", "CAPRESE", "CAPRI_BIC", "CAPRI_AIC",
                 "CBN_ot", "MCCBN")

    cat("\n     Doing MHN")
    time_schill <- system.time(
        out_schill <- do_MHN2(x, lambda = 1/nrow(x)))["elapsed"]

    cat("\n  time MHN = ", time_schill)
    cpm_out_others <- all_methods(x, cores_cbn = cores_cbn)

    pre_trans_mat_others <- lapply(cpm_out_others[methods], cpm_access_genots_paths_w_simplified)

    cat("\n    getting transition matrices for all non-mhn methods ")
    ## Unweighted
    uw <- lapply(pre_trans_mat_others, function(x) rowScaleMatrix(x$fgraph))
    ## Weighted
    wg <- lapply(pre_trans_mat_others[c("OT", "MCCBN", "CBN_ot")],
                 function(x) x$trans_mat_genots)
    ## Diagonal
    td <- lapply(pre_trans_mat_others[c("MCCBN", "CBN_ot")],
                 function(x) trans_rate_to_trans_mat(x$weighted_fgraph,
                                                     method = "uniformization"))
    ## ## Paranoid check
    ## wg2 <- lapply(pre_trans_mat_others[c("OT", "MCCBN", "CBN_ot")],
    ##               function(x) trans_rate_to_trans_mat(x$weighted_fgraph,
    ##                                                   method = "competingExponentials"))
    
    ## stopifnot(identical(wg, wg2))
    
    
    return(list(OT = wg$OT, OT_u = uw$OT,
         CBN = wg$CBN_ot, CBN_td = td$CBN_ot, CBN_uw = uw$CBN_ot,
         MCCBN = wg$MCCBN, MCCBN_td = td$MCCBN, MCCBN_uw = uw$MCCBN,
         CAPRESE = uw$CAPRESE,
         CAPRI_BIC = uw$CAPRI_BIC,
         CAPRI_AIC = uw$CAPRI_AIC,
         MHN = out_schill$transitionMatrixCompExp,
         MHN_td = out_schill$transitionMatrixTimeDiscretized
         ))
}


library(codetools)
checkUsageEnv(env = .GlobalEnv)



#####################################################################

## ## mini example
## Dat1 <- readRDS(file="./MHN/data/BreastCancer.rds") [1:50, 1:4]
## colnames(Dat1) <- LETTERS[1:ncol(Dat1)]

## Dat1 <- as.matrix(Dat1)

## simplerun <- all_methods_2_trans_mat(Dat1)

## rm(simplerun)





## Getting an idea of sizes of data if we did not use sparse matrices
if(FALSE) {
## 8 bytes per float as can be checked doing
u <- runif((2^10) * (2^10))
print(object.size(u), units = "b")
8 * length(u)

uu <-  runif((2^15) * (2^15))
length(uu)
print(object.size(uu), units = "b")
print(object.size(uu), units = "GB") ## 8 GB

length(uu) * 8 /(1024 * 1024 * 1024)

## so an object for 2^20 * 2^20 would take

((2^20)^2) * 8 / (1024 * 1024 * 1024)

## or 8192 GB which is what R complaints about if I try this
uuu <-  runif((2^20) * (2^20))

## Yes, there are no memory limits in the machines
## install.packages("devtools", dependencies = TRUE)
##devtools::install_github("krlmlr/ulimit")
library(ulimit)
ulimit::memory_limit()

}



## MHN does not run if we use mclapply after setting the threads for OMP >
## 1

## all_data_out <- mclapply(all_data, all_methods_2_trans_mat, mc.cores
## = 36)

## Yes, you can use mclapply, just make sure to set the threads of OMP to 1
## The following shows it. It just did not make any sense above.
if(FALSE) {

    lapply(all_data, dim)

    RhpcBLASctl::omp_set_num_threads(36)

    system.time(cucu3 <- lapply(all_data[c(3, 4, 20, 21:24)], do_MHN))

    RhpcBLASctl::omp_set_num_threads(1)
    system.time(cucu4 <- mclapply(all_data[c(3, 4, 20, 21:24)], do_MHN,
                                  mc.cores = 36))
    stopifnot(identical(cucu3, cucu4))
}


## If you wanted to launch a single process, you could do this
## ## Since I have a lot fewer than 36 data sets, and easier to catch issues
## ## as they show up, do not parallelize. But use multiple threads
## all_data_out <- lapply(1:length(all_data),
##                        function(i) {
##                            cat("\n #################################")
##                            cat("\n #################################")
##                            cat("\n\n   Doing data set ", names(all_data)[i])
##                            cat(" with dim = ", dim(all_data[[i]]))
##                            return(all_methods_2_trans_mat(all_data[[i]]))
##                        }
##                        )
