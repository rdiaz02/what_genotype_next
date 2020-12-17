## Recall weighted_fgraph in all_paths_CPMs*.rds
## keeps probs/rates. Only for CBN_ot and MCCBN

## See file: new-CPMs-weighted-paths-genots.R
## see function: cpm_access_genots_paths_w

## where comments say

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


## OT. Nope, cannot do it. What OT returns are kind of conditional probabilities.
## Recall, in paper with Claudia: "It is also possible to perform a
## similar operation with the output of OT, and use the edge weights from the
## fits of OT to obtain the probabilities of transition to each descendant
## genotype and, from them, the probabilities of the different paths to the
## single fitness maximum.  It must be noted that these probabilities are not
## really returned by the model, since the OTs used are untimed oncogenetic
## trees (Desper et al., 1999; Szabo and Boucher, 2008)."

date()  ## in Draco, runs in about 2 minutes



## transition rate matrix -> transition matrix
##    if method == uniformization, we assume no diagonal entry
##    because we compute it
trans_rate_to_trans_mat <- function(x,
                                    method = c("competingExponentials",
                                               "uniformization"),
                                    paranoidCheck = TRUE) {
    method <- match.arg(method)
    if(method == "competingExponentials") {
        ## Conditional on there being a transition, thus diagonals are
        ## zero
        sx <- rowSums(x)
        tm <- sweep(x, 1, sx, "/")
        tm[nrow(tm), ] <- 0 ## last row is 0
        if(paranoidCheck) {
            stopifnot(isTRUE(all.equal(c(rep(1, nrow(x) - 1), 0),
                                       rowSums(tm),
                                       check.attributes = FALSE)))  
        }
    } else if (method == "uniformization") {
        ## The time-discretized version. Diagonals can be non-zero
        ## Using uniformization method, as in Schill et al., 2020,
        ## "Modelling cancer progression using Mutual Hazard Networks",
        ## Bioinformatics
        ## Fig.6 legend, who cite Grassman, 1977
        ## (Wikipedia also has an entry:
        ## https://en.wikipedia.org/wiki/Uniformization_(probability_theory )
        dd <- diag(x)
        if(!isTRUE(all(dd == 0))) stop("Diagonal of x is not 0")
        diag(x) <- -1 * rowSums(x)
        gamma <- max(abs(diag(x)))
        tm <- diag(nrow(x)) + x/gamma
        if(paranoidCheck) {
            stopifnot(isTRUE(all.equal(rep(1, nrow(x)),
                                       rowSums(tm),
                                       check.attributes = FALSE)))
        }
    }
    return(tm)
}


## list with weighted fitness graph (transition rate matrix from CBN and MCCBN) ->
##    transition matrix, with non-zero diagonals
##   a wrapper to handle error conditions
##   And be verbose
unif_trans_mat <- function(x) {
    wft <- x[["weighted_fgraph"]]
    if( (length(wft) == 1) &&
        (wft == "ERROR_CPM_ANALYSIS")) {
        cat("\n     ************* No weighted_fgraph: ERROR_CPM_ANALYSIS  ****** \n")
        return("ERROR_CPM_ANALYSIS")
    } else {       
        trans_rate_to_trans_mat(wft, method = "uniformization")
    }
}


## single CPM output file name ->
##      transition matrices with non-zero diagonals
do_unif_trans_mat <- function(f) {
    x <- readRDS(f)
    landscape <- x[["ID"]]
    detect <- x[["detect"]]
    size_split <- x[["size_split"]]
    typeLandscape <- x[["typeLandscape"]]
    cat("\n Doing landscape ", landscape,
        " detect ", detect,
        " size_split", size_split,
        " typeLandscape ", typeLandscape,
        "\n")
    ngenes <- as.integer(as.character(x[["ngenes"]]))
    if(! (ngenes %in% c(7, 10) )) stop("ngenes not 7 or 10")

    xo <- x[["out_paths_genots"]]


    methods <- c("CBN_ot", "MCCBN")
    ## reps <- length(xo)
    ## xg <- expand.grid(methods, 1:reps, stringsAsFactors = FALSE)
    ## colnames(xg) <- c("Method", "Replicate")
    ## xg$fname <- paste0("tpaths",
    ##                "_ID_", landscape,
    ##                "_detect_", detect,
    ##                "_size_split_", size_split,
    ##                "_typeLandscape_", typeLandscape,
    ##                "_rep_", xg[, "Replicate"],
    ##                "__m__", xg[, "Method"], ## because of CBN_ot naming
    ##                ".rds"
    ##                )

    ## So as to preserve original structure
    m1 <- lapply(xo,
           function(x) { ## loop over replicates
               lapply(x[methods],
                      function(z) {
                          ## same name as from Schill's MHN'
                          return(list(transitionMatrixTimeDiscretized =
                                          unif_trans_mat(z)))
                      }
                      )
    })

    
    ## m1 <- Map(function(replicate, method)
    ##     trans_rate_to_trans_mat(xo[[replicate]][[method]][["weighted_fgraph"]],
    ##                             method = "uniformization"),
    ##     xg[, "Replicate"],
    ##     xg[, "Method"]
    ##     )

    
    theout <- c(
        x[1:12],
        ## trans_matrix_diagonal = list(m1) ## NOPE! use same names
        ## as the component underneath
        ## has the right name: uniformiz_trans_matrix
        out_paths_genots = list(m1)        
    )
    ## Preserve all original name, but replace
    ## all_paths_CPMs_genots by transition_matrices_diag
    fname <- paste0(
        "trans_mat_time_discr_",
        strsplit(f, "all_paths_CPMs_genots_")[[1]][2])
    saveRDS(theout, fname)
    cat("\n  done  ", fname , "\n")

}

## How do names of CBN components match those from MHN, from Schill
## in the lists? From "run-Schill-MHN-trans-mat.R"

    ## Names of components and how they match those from CBN, etc
    ## theta: there is no equivalent
    ## transitionRateMatrix: weighted_fgraph (but only for CBN and MCCBN, not OT)
    ## transitionMatrixCompExp: trans_mat_genots
    ## transitionMatrixTimeDiscretized: transitionMatrixTimeDiscretized




library(codetools)
checkUsageEnv(env = .GlobalEnv)


run_example <- FALSE

if(run_example){
    ## weighted_fgraph keeps probs/rates.
    u <- readRDS("all_paths_CPMs_genots_ID_10fs6bc54hwzLal8K__rnst_w33ijb81haFWU1cN__ngenes_7__initSize_1e+06__typeLandscape_RMF__size_split_4000__detect_unif.rds")

    ## Note output probs, not rates. Can't use them
    u$out_paths_genots[[1]][["OT"]][["weighted_fgraph"]]

    ## rates
    u$out_paths_genots[[1]][["CBN_ot"]][["weighted_fgraph"]]
    u$out_paths_genots[[1]][["MCCBN"]][["weighted_fgraph"]]

    unif_trans_mat(u$out_paths_genots[[1]][["MCCBN"]][["weighted_fgraph"]])
    unif_trans_mat(u$out_paths_genots[[1]][["CBN_ot"]][["weighted_fgraph"]])
}



cpmf <- dir(path = "../cpm-analyzed-all-selected",
            pattern = glob2rx("all_paths_CPMs*.rds"),
            full.names = TRUE)
## 35  replicates * 3 initSize * 2 mutations * 3 fitnes landscape
##         types * 2 ngenes
## And we do 3 sampling schemes and 3 sizes of sampling
stopifnot(length(cpmf) == (3 * 3 * 1260))
## and we will produce this many files
## 2 methods, 5 replicates. Oh man.
length(cpmf) * 2 * 5


## get compiler to run
for(i in 1:4) {
    do_unif_trans_mat(cpmf[i])
}

library(parallel)

mclapply(cpmf,
         function(y)
             do_unif_trans_mat(y),
         mc.cores = detectCores(),
         mc.preschedule = TRUE,
         mc.allow.recursive = FALSE)

date()
gc()





## launch as
## nohup R-devel --vanilla --slave -f time-discretized-CBN.R &> time-discretized-CBN.Rout &

## check number of missing fgraphs: 211
## grep "No weighted_fgraph" uniformiz.Rout | wc -l


