## Copyright 2020 Ramon Diaz-Uriarte, Juan Diaz Colunga

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


library(matrixStats)
library(parallel)
library(spatstat) ## weighted stuff
library(spatstat.geom) ## Where weighted.quantile now is
library(Hmisc) ## wtd.quantile

## What a hack!!
## Could have used Hmisc::wtd.quantile
## or spatstat::weighted.quantile
weighted_summary <- function(x, w) {
    wn <- round(w/min(w))
    y <- rep.int(x, wn)
    return(quantile(y, c(0, .1, .15, .2, .25, .4, .5, .75, 1.0)))
}

weighted_summary_data <- function(data) {
    x <- data$min_js
    w <- data$sampledProp
    weighted_summary(x, w)
}


weighted_summary_data2 <- function(data,
                                   probs = c(0, .1, .15, .2, .25, .3, .35, .4, .5, .75, 1.0),
                                   sqrt_inv = FALSE) {
    x <- data$min_js
    w <- data$sampledProp
    if(sqrt_inv) x <- x^2
    spatstat.geom::weighted.quantile(x, w,
                                probs = probs)
}


how_many_less_than_js <- function(data, threshold = sqrt(0.15)) {
    x <- data$min_js
    w <- data$sampledProp
    sum_w_less <- sum(w[x <= threshold])
    sum_w <- sum(w)
    return(c(prop_less = sum_w_less/sum_w, sum_w = sum_w))
}

## a fitted glmertree and a set of variables to show in histogram ->
##    list with data in a form usable by my_node_bp4 plotting function
##   variables all finish as "_is_best"
##   omit that part, which will be rm
create_terminal_plot_data <- function(tree,
                                      data = data_for_glmertree_plots,
                                      best_vars = c("MHN", "MHN_td", "CBN",
                                                    "CBN_td", "OT"),
                                      weighted = FALSE) {
    ## Sometimes, there can be rounding error issues
    rules_nodes_2 <- my.list.rules.party(tree,
                                         i = nodeids(tree, terminal = TRUE))
    ## Needed if we pass a tree like tree[2]
    ## NO!!! This would mask problems like not using the first rule of the split
    ## names(rules_nodes_2) <- nodeids(tree, terminal = TRUE)
    terminal_plot_data <- vector(mode = "list", length = length(nodeids(tree)))
    names(terminal_plot_data) <- 1:length(terminal_plot_data)
    best_vars <- paste0(best_vars, "_is_best")
    if(weighted) {
        mfw <- nrow(data)/sum(data$sampledProp)
        data$weights <- mfw * data$sampledProp
    }
    for(i in names(rules_nodes_2)) {
        terminal_plot_data[[i]]$min_js <-
            data[eval(parse(text = rules_nodes_2[i] )), min_js]
        tmp <- data[eval(parse(text = rules_nodes_2[i] )), ..best_vars]
        colnames(tmp) <- gsub("_is_best", "", colnames(tmp))
        terminal_plot_data[[i]]$is_best <- tmp
        rm(tmp)
        if(weighted) {
            terminal_plot_data[[i]]$weights <-
                data[eval(parse(text = rules_nodes_2[i] )), weights]
        } 
    }
    gc()
    return(terminal_plot_data)
}

## if average_rank = TRUE, show average rank
## if average_js = TRUE, show average js
## yes, had I thought this, could be cleaner
## show: one of "is_best", "rank", "js"
create_terminal_plot_data_2 <- function(tree,
                                      data = data_for_glmertree_plots,
                                      best_vars = c("MHN", "MHN_td", "CBN",
                                                    "CBN_td", "OT"),
                                      weighted = FALSE,
                                      show_bottom = "is_best") {
    ## Sometimes, there can be rounding error issues
    rules_nodes_2 <- my.list.rules.party(tree,
                                         i = nodeids(tree, terminal = TRUE))
    ## Needed if we pass a tree like tree[2]
    ## NO!!! This would mask problems like not using the first rule of the split
    ## names(rules_nodes_2) <- nodeids(tree, terminal = TRUE)
    terminal_plot_data <- vector(mode = "list", length = length(nodeids(tree)))
    names(terminal_plot_data) <- 1:length(terminal_plot_data)
    if(show_bottom == "is_best")
        best_vars <- paste0(best_vars, "_is_best")
    else if(show_bottom == "rank")
        best_vars <- paste0("rank_", best_vars, "_js")
    else if(show_bottom == "js")
        best_vars <- paste0(best_vars, "_js")
    
    if(weighted) {
        mfw <- nrow(data)/sum(data$sampledProp)
        data$weights <- mfw * data$sampledProp
    }
    for(i in names(rules_nodes_2)) {
        terminal_plot_data[[i]]$min_js <-
            data[eval(parse(text = rules_nodes_2[i] )), min_js]
        tmp <- data[eval(parse(text = rules_nodes_2[i] )), ..best_vars]
        if(show_bottom == "is_best")
            colnames(tmp) <- gsub("_is_best$", "", colnames(tmp))
        else if(show_bottom == "rank")
            colnames(tmp) <- gsub("^rank_", "", gsub("_js$", "", colnames(tmp)))
        else if(show_bottom == "js")
            colnames(tmp) <- gsub("_js$", "", colnames(tmp))
        terminal_plot_data[[i]]$is_best <- tmp
        rm(tmp)
        if(weighted) {
            terminal_plot_data[[i]]$weights <-
                data[eval(parse(text = rules_nodes_2[i] )), weights]
        } 
    }
    gc()
    return(terminal_plot_data)
}



## spit out the total number of observations in the terminal
## nodes, and verify they sum 1379997, the total number of observations used
## The number to check against, default, should also be identical to the number
## in the Rout with the fits. "nrow data set"
node_counts_check <- function(object, default = 1379997L) {
   ## checks:
    termin_nodes_n0 <- lapply(nodeids(object, terminal = TRUE),
                             function(x) nrow(data_party(object, x)))
    termin_nodes_n <- nodeapply(object,
                                ids = nodeids(object, terminal = TRUE),
                                FUN = function(n) n$info$nobs)
    stopifnot(identical(sum(unlist(termin_nodes_n0)),default))
    if(is.null(object$fitted[["(weights)"]])) {
        stopifnot(identical(sum(unlist(termin_nodes_n)),
                            sum(unlist(termin_nodes_n0))))
        print(termin_nodes_n)
    } else {
        stopifnot(isTRUE(all.equal(
            sum(unlist(termin_nodes_n)),
            sum(object$fitted[["(weights)"]]) )
            ))
        print(c(termin_nodes_n, termin_nodes_n0))
    }
}


## Like node_counts_check, but with more detail.
## Remember that matching the split rules sometimes leads to minor
## disadjustements because of rounding issues. This looks into it. 
node_counts_check2 <- function(object, data,
                               default = 1379997L) {
   ## checks:
    termin_nodes_n0 <- lapply(nodeids(object, terminal = TRUE),
                             function(x) nrow(data_party(object, x)))
    termin_nodes_n <- nodeapply(object,
                                ids = nodeids(object, terminal = TRUE),
                                FUN = function(n) n$info$nobs)
    stopifnot(identical(sum(unlist(termin_nodes_n0)),default))

    if(is.null(object$fitted[["(weights)"]])) {
        stopifnot(identical(sum(unlist(termin_nodes_n)),
                            sum(unlist(termin_nodes_n0))))
    } else {
        stopifnot(isTRUE(all.equal(
            sum(unlist(termin_nodes_n)),
            sum(object$fitted[["(weights)"]]) )
            ))
    }
  

    nrb <- unlist(lapply(data, function(x) nrow(x$is_best)))
    nb <- unlist(lapply(data, function(x) length(x$min_js)))
    nrb <- nrb[nrb > 0]
    nb <- nb[nb > 0]
    stopifnot(identical(nrb, nb))
    cat("\n You can expect minor mismatches\n")
    if(is.null(object$fitted[["(weights)"]])) {
        print(cbind(nrb, unlist(termin_nodes_n)))
        print(nrb - unlist(termin_nodes_n))
        ## print(termin_nodes_n)
    } else {
        print(cbind(nrb,
                    nodes_n = unlist(termin_nodes_n),
                    nodes_n0 = unlist(termin_nodes_n0)))
        cat("\nnrb - n0\n")
        print(nrb - unlist(termin_nodes_n0))
        cat("nrb - n\n")
        ## of coyurse, here possibly major mismatches
        ## as these are different things
        print(nrb - unlist(termin_nodes_n))
        ## print(termin_nodes_n)
        ## print(termin_nodes_n0)
    }
}


#### A lot of the code below is for pruning

## If weight is not null, use a weighted median
node_medians <- function(obj, weighted, cores = detectCores()) {
    if(!weighted) {
        medians <- mclapply(nodeids(obj, terminal = TRUE),
                            function(x) median(data_party(obj, x)$min_js),
                            mc.cores = cores)
    } else {
        medians <- mclapply(nodeids(obj, terminal = TRUE),
                          function(x) weightedMedian(data_party(obj, x)$min_js,
                                                     data_party(obj, x)[["(weights)"]]),
                          mc.cores = cores)
        
    }
    names(medians) <- nodeids(obj, terminal = TRUE)
    return(medians)
}

## as it says: show median and intercept
scatterplot_median_intercept <- function(obj, weighted, plot = FALSE) {
    medians <- unlist(node_medians(obj, weighted))
    intercepts <- unlist(nodeapply(obj, ids = nodeids(obj, terminal = TRUE),
                   FUN = function(n)
                       n$info$coefficients[["(Intercept)"]]))
    if(plot) {
        plot(intercepts ~ medians)
        abline(0, 1, lty = 1)
        abline(lm(intercepts ~ medians), lty = 2)
    }
    print(cor.test(medians, intercepts))
    print(cor.test(medians, intercepts, method = "spearman"))
    print(summary(lm(intercepts ~ medians)))
    print(summary(lm(intercepts ~ -1 + medians)))
}

## Note: transformations are handled externally. For example, with sqrt
## you might want to pass maxFit = sqrt(0.25), mindiff = sqrt(0.015), etc
## But that handling of transformations can be mistaken!!!
## sqrt(a) - sqrt(b) compared to sqrt(c)!!!
## Eh!!! sqrt(0.165) - sqrt(0.15) = 0.019, but sqrt(0.015) = 0.122!!!

should_we_prune <- function(obj, maxFit = 0.25, minN = 100,
                            mindiff = 0.005,
                            statistic = "median",
                            weighted = FALSE) {
    message("\n BEWARE of transformations. This might not do what you want!!")

    if(statistic == "median") {
        statistic <- node_medians(obj, weighted)
    } else if (statistic == "intercept") {
         statistic <-
             nodeapply(obj, ids = nodeids(obj, terminal = TRUE),
                   FUN = function(n)
                       n$info$coefficients[["(Intercept)"]])
    }

    termin_nodes_n <- nodeapply(obj,
                                ids = nodeids(obj, terminal = TRUE),
                                FUN = function(n) n$info$nobs)

    ## We do not care about these nodes
    node_go <- names(which((
        (statistic > maxFit) | (termin_nodes_n < minN))
        ))

    ## We prune a node whose two kids must be removed.
    ## Find parents of nodes in node_go. A node to be pruned
    ## is duplicated as a parent.
    
    parents <- parent_kids(obj)
    candidate_parents <- parents[node_go]
    which_prune <- unname(candidate_parents[duplicated(candidate_parents)])

    ## ## But we also want to prune a node if one its children has n < minN,
    ## ## regardless of its intercept
    ## ## NOPE
    ## which_prune <- unique(c(which_prune,
    ##                         parents[names(which(termin_nodes_n < minN))]))

    ## Finally collapse two children node if they differ less than mindiff
    ## in their statistic
    ## Yeah, could have used same logic above. Oh well.
    parents_of_terminal <- parents[names(statistic)]
    parents_with_two <- parents_of_terminal[duplicated(parents_of_terminal)]
    for(p in parents_with_two) {
        sister_stat <- statistic[names(parents_of_terminal[parents_of_terminal == p])]
        if(abs(diff(unlist(sister_stat))) < mindiff) {
            which_prune <- unique(c(which_prune, p))
        }
    }
    return(list(which_prune = which_prune))
}


## We assume data are sqrt transformed, so we backtransform them
should_we_prune_sqrt <- function(obj,
                                 maxFit = 0.15, minN = 100,
                            mindiff = 0.015,
                            statistic = "median",
                            weighted = FALSE) {
    

    if(statistic == "median") {
        statistic <- node_medians(obj, weighted)
    } else if (statistic == "intercept") {
         statistic <-
             nodeapply(obj, ids = nodeids(obj, terminal = TRUE),
                   FUN = function(n)
                       n$info$coefficients[["(Intercept)"]])
    }

    statistic <- lapply(statistic, function(z) z^2)
    
    termin_nodes_n <- nodeapply(obj,
                                ids = nodeids(obj, terminal = TRUE),
                                FUN = function(n) n$info$nobs)

    ## We do not care about these nodes
    node_go <- names(which((
        (statistic > maxFit) | (termin_nodes_n < minN))
        ))


    ## We prune a node whose two kids must be removed.
    ## Find parents of nodes in node_go. A node to be pruned
    ## is duplicated as a parent.

    parents <- parent_kids(obj)
    candidate_parents <- parents[node_go]
    which_prune <- unname(candidate_parents[duplicated(candidate_parents)])

    if(length(which_prune)) {
        message(" ... prunning by maxFit, nodes ", paste(which_prune, collapse = " "), "\n")
    } 
    
    ## ## But we also want to prune a node if one its children has n < minN,
    ## ## regardless of its intercept
    ## ## NOPE!
    ## which_prune <- unique(c(which_prune,
    ##                         parents[names(which(termin_nodes_n < minN))]))

    ## Finally collapse two children node if they differ less than mindiff
    ## in their statistic
    ## Yeah, could have used same logic above. Oh well.
    parents_of_terminal <- parents[names(statistic)]
    parents_with_two <- parents_of_terminal[duplicated(parents_of_terminal)]
    for(p in parents_with_two) {
        sister_stat <- statistic[names(parents_of_terminal[parents_of_terminal == p])]
        if(abs(diff(unlist(sister_stat))) < mindiff) {
            if( !(p %in% which_prune) ) {
                message("  .... prunning by mindiff; p = ", p, ". Only by mindiff")
            } else {
                message("  .... prunning by mindiff; p = ", p)
            }
            which_prune <- unique(c(which_prune, p))
        }
    }
    return(list(which_prune = which_prune))
}



my_prune <- function(obj, maxFit = 0.25, minN = 100, mindiff = 0.005,
                     statistic = "median", weighted = FALSE) {
    prune <- TRUE
    i <- 0
    while(prune) {
        i <- i + 1
        to_prune <- should_we_prune(obj, maxFit, minN, mindiff, statistic, weighted)
        if(length(to_prune$which_prune) == 0) break
        print(to_prune$which_prune)
        cat("\n      ... Iteration ", i, ". Prunning id = ")
        ## Prune one by one! Because of the relabelling
        id <- sort(to_prune$which_prune)[1]
        cat(id, "\n")
        ## obj <- prune_relabel(obj, id)
        obj <- nodeprune(obj, id)


        ## Intermediate check: catch possible issues ASAP
        ## and do not allow possible errors to cancel each other
        nr_2i <- unlist(nodeapply(obj,
                                 ids = nodeids(obj, terminal = TRUE),
                                 FUN = function(n) n$info$nobs))
        if(!weighted) {
            nr_dpi <- unlist(lapply(nodeids(obj, terminal = TRUE),
                                   function(x) nrow(data_party(obj, x))))
        } else {
            nr_dpi <- unlist(lapply(nodeids(obj, terminal = TRUE),
                                   function(x) sum(data_party(obj, x)[["(weights)"]])))
        }
        ## if(!isTRUE(all.equal(nr_dpi, unname(nr_2i)))) browser()
        stopifnot(isTRUE(all.equal(nr_dpi, unname(nr_2i))))
        rm(nr_dpi, nr_2i)
    }
    ## Final check
    nr_2 <- unlist(nodeapply(obj,
                                 ids = nodeids(obj, terminal = TRUE),
                                 FUN = function(n) n$info$nobs))
    if(!weighted) {
        nr_dp <- unlist(lapply(nodeids(obj, terminal = TRUE),
                               function(x) nrow(data_party(obj, x))))
    } else {
        nr_dp <- unlist(lapply(nodeids(obj, terminal = TRUE),
                               function(x) sum(data_party(obj, x)[["(weights)"]])))
    }
    ## if(!isTRUE(all.equal(nr_dp, unname(nr_2)))) browser()
    stopifnot(isTRUE(all.equal(nr_dp, unname(nr_2))))
    return(obj)
}


## Note: transformations are handled externally. For example, with sqrt
## you might want to pass maxFit = sqrt(0.25), mindiff = sqrt(0.015), etc
## But this would not work!!! sqrt(a) - sqrt(b) compared to sqrt(c)!!!
## Eh!!! sqrt(0.165) - sqrt(0.15) = 0.019, but sqrt(0.015) = 0.122!!!

fast_prune <- function(obj, maxFit = 0.25, minN = 100, mindiff = 0.005,
                       statistic = "median", weighted = FALSE) {
    message("\n BEWARE of transformations!!!! This might not do what you want!!!")
    prune <- TRUE
    i <- 0
    while(prune) {
        i <- i + 1
        to_prune <- should_we_prune(obj, maxFit, minN, mindiff, statistic, weighted)
        if(length(to_prune$which_prune) == 0) break
        print(to_prune$which_prune)
        cat("\n      ... Iteration ", i, ". Prunning ids = ")
        ## No need to prune one by one anymore
        id <- sort(to_prune$which_prune)
        cat(id, "\n")
        ## obj <- prune_relabel(obj, id)
        obj <- nodeprune(obj, id)


        ## Intermediate check: catch possible issues ASAP
        ## and do not allow possible errors to cancel each other
        nr_2i <- unlist(nodeapply(obj,
                                 ids = nodeids(obj, terminal = TRUE),
                                 FUN = function(n) n$info$nobs))
        if(!weighted) {
            nr_dpi <- unlist(lapply(nodeids(obj, terminal = TRUE),
                                   function(x) nrow(data_party(obj, x))))
        } else {
            nr_dpi <- unlist(lapply(nodeids(obj, terminal = TRUE),
                                   function(x) sum(data_party(obj, x)[["(weights)"]])))
        }
        ## if(!isTRUE(all.equal(nr_dpi, unname(nr_2i)))) browser()
        stopifnot(isTRUE(all.equal(nr_dpi, unname(nr_2i))))
        rm(nr_dpi, nr_2i)
    }
    ## Final check
    nr_2 <- unlist(nodeapply(obj,
                                 ids = nodeids(obj, terminal = TRUE),
                                 FUN = function(n) n$info$nobs))
    if(!weighted) {
        nr_dp <- unlist(lapply(nodeids(obj, terminal = TRUE),
                               function(x) nrow(data_party(obj, x))))
    } else {
        nr_dp <- unlist(lapply(nodeids(obj, terminal = TRUE),
                               function(x) sum(data_party(obj, x)[["(weights)"]])))
    }
    ## if(!isTRUE(all.equal(nr_dp, unname(nr_2)))) browser()
    stopifnot(isTRUE(all.equal(nr_dp, unname(nr_2))))
    return(obj)
}


## maxFit and mindiff are given in the original scale!!!
fast_prune_sqrt <- function(obj, maxFit = 0.15, minN = 100, mindiff = 0.015,
                     statistic = "median", weighted = FALSE) {
    prune <- TRUE
    i <- 0
    while(prune) {
        i <- i + 1
        to_prune <- should_we_prune_sqrt(obj, maxFit, minN, mindiff, statistic, weighted)
        if(length(to_prune$which_prune) == 0) break
        print(to_prune$which_prune)
        cat("\n Iteration ", i, ". Prunning ids = ")
        ## No need to prune one by one anymore
        id <- sort(to_prune$which_prune)
        cat(id, "\n\n")
        ## obj <- prune_relabel(obj, id)
        obj <- nodeprune(obj, id)


        ## Intermediate check: catch possible issues ASAP
        ## and do not allow possible errors to cancel each other
        nr_2i <- unlist(nodeapply(obj,
                                 ids = nodeids(obj, terminal = TRUE),
                                 FUN = function(n) n$info$nobs))
        if(!weighted) {
            nr_dpi <- unlist(lapply(nodeids(obj, terminal = TRUE),
                                   function(x) nrow(data_party(obj, x))))
        } else {
            nr_dpi <- unlist(lapply(nodeids(obj, terminal = TRUE),
                                   function(x) sum(data_party(obj, x)[["(weights)"]])))
        }
        ## if(!isTRUE(all.equal(nr_dpi, unname(nr_2i)))) browser()
        stopifnot(isTRUE(all.equal(nr_dpi, unname(nr_2i))))
        rm(nr_dpi, nr_2i)
    }
    ## Final check
    nr_2 <- unlist(nodeapply(obj,
                                 ids = nodeids(obj, terminal = TRUE),
                                 FUN = function(n) n$info$nobs))
    if(!weighted) {
        nr_dp <- unlist(lapply(nodeids(obj, terminal = TRUE),
                               function(x) nrow(data_party(obj, x))))
    } else {
        nr_dp <- unlist(lapply(nodeids(obj, terminal = TRUE),
                               function(x) sum(data_party(obj, x)[["(weights)"]])))
    }
    ## if(!isTRUE(all.equal(nr_dp, unname(nr_2)))) browser()
    stopifnot(isTRUE(all.equal(nr_dp, unname(nr_2))))
    return(obj)
}



## Seems to be fixed as of 2020-07-22.
## Remember issue in https://stackoverflow.com/questions/62614893/partykit-data-party-on-nodepruned-tree-how-when-should-nodeprune-rename-nod
## This prunes and fixes the ids in fitted
prune_relabel <- function(obj, id) {
    stopifnot(length(id) == 1)
    objb <- nodeprune(obj, id)
    children <- rep(id, 2) + c(1, 2)
    rename_ids <- nodeids(obj, terminal = TRUE)
    rename_ids <- rename_ids[rename_ids > max(children)]
    ## names(rename_ids) <- rename_ids
    ## rename_ids <- rename_ids - 2
    rename_ids <- sort(rename_ids)
    ## simpler than lookup table since not all ids in rename_ids
    for (ri in rename_ids) 
        objb$fitted[["(fitted)"]][objb$fitted[["(fitted)"]] == ri] <- (ri - 2)
    return(objb)
}

## Based on https://stackoverflow.com/a/53085559/3670562
## Returns vector with names kids and entry parents
## So a lookup vector with children as names
parent_kids <- function(obj) {
    tmp <- t(sapply(as.list(obj$node), function(n) {
        if(is.null(n$kids)) c(n$id, NA, NA) else c(n$id, n$kids)
    }))
    ## colnames(tmp) <- c("parent", "kid1", "kid2")
    ## tmp
    kid_to_parent <- cbind(parent = c(tmp[, 1], tmp[, 1]),
                       kid = c(tmp[, 2], tmp[, 3]))
    kid_to_parent <- na.omit(kid_to_parent)
    kp <- kid_to_parent[, 1]
    names(kp) <- kid_to_parent[, 2]
    return(kp)
}


## Wrappers for plotting. 

## Plot the tree. Since I do not splits here, this is a way to check
## and the plotting includes checks too.
suppressWarnings(try(rm(terminal_plot_data), silent = TRUE))
data_and_plot <- function(obj, data, counts = 106000L, main = NULL,
                          best_vars =  c("MHN"
                                        , "MHN_td"
                                        , "CBN"
                                        , "CAPRESE"
                                       , "CBN_td"),
                          size_x_text = 11,
                          size_y_text = 12,
                          angle_x_text = 45,
                          width = 36,
                          height = 18,
                          weighted = FALSE,
                          sqrt_transf = FALSE,
                          show_mean = FALSE,
                          show_bottom = "is_best") {
    if(is.null(main)) main <- deparse(substitute(obj))

    suppressWarnings(try(rm(terminal_plot_data), silent = TRUE))
    suppressWarnings(try(rm(terminal_plot_data, envir = .GlobalEnv), silent = TRUE))
    datatmp <- create_terminal_plot_data_2(obj
                                       , data = data
                                       , best_vars = best_vars
                                       , weighted = weighted  
                                       , show_bottom = show_bottom)
    ## node_counts_check(obj, counts)
    node_counts_check2(obj, datatmp, counts)
    terminal_plot_data <<- datatmp
    pdf(file = paste0(main, ".pdf"),
        width = width, height = height)
    my_tree_plot(obj,
             inner_panel = my_node_inner2,
             terminal_panel = my_node_bp5,
             tp_args = list(size_x_text = size_x_text,
                            angle_x_text = angle_x_text,
                            size_y_text = size_y_text,
                            sqrt_transf = sqrt_transf,
                            show_mean = show_mean,
                            show_bottom = show_bottom),
             main = main
             )
    dev.off()
    suppressWarnings(try(rm(terminal_plot_data), silent = TRUE))
    suppressWarnings(try(rm(terminal_plot_data, envir = .GlobalEnv), silent = TRUE))
}

## Fancier plots, split. Always check with above.
split_tree_plot <- function(obj, data,
                            width_2 = 8.5, 
                            height_2 = 8.5, 
                            yunit_2 = 2.8, 
                            width_rest = 32, 
                            height_rest = 16, 
                            size_x_text = 11,
                            size_y_text = 12,
                            angle_x_text = 45, 
                            width_inner = 12, 
                            height_inner = 9, 
                            hspace = "95mm", 
                            vspace = "-15mm", 
                            main = NULL,
                            best_vars = c("MHN"
                                        , "MHN_td"
                                        , "CBN"
                                        , "CAPRESE"
                                        , "CBN_td"),
                            weighted = FALSE,
                            sqrt_transf = FALSE
                            ) {
    rm(terminal_plot_data)
    rm(terminal_plot_data, envir = .GlobalEnv)
    if(is.null(main)) main <- deparse(substitute(obj))
    
    pk <- parent_kids(obj)
    upper_children <- names(pk[pk == 1])
    obj_branches <- nodeprune(obj, as.integer(upper_children))
    ## ## Not needed anymore the next line!
    ## ## The problem in https://stackoverflow.com/q/62614893/3670562
    ## obj_branches$fitted[["(fitted)"]] <-
    ##     c(2, 3)[match(obj_branches$fitted[["(fitted)"]], as.integer(upper_children))]
    
    ## Rename
    names(obj_branches) <- c("1", upper_children)
    
    ## Checks
    dim(data_party(obj, 1))
    dim(data_party(obj_branches, 1))
    dim(data_party(obj_branches, 2))
    dim(data_party(obj_branches, 3))
    stopifnot(identical(dim(data_party(obj, 1)),
                        dim(data_party(obj_branches, 1))))
    stopifnot(nrow(data_party(obj, 1)) ==
              (nrow(data_party(obj_branches, 2)) +
               nrow(data_party(obj_branches, 3))
              ))
    
    rm(terminal_plot_data)
    rm(terminal_plot_data, envir = .GlobalEnv)
    
    datatmp <- create_terminal_plot_data(obj,
                                         data = data,
                                         best_vars = best_vars,
                                         weighted = weighted)
    
    system("rm t_first_split.pdf")
    system("rm t_right_split.pdf")
    system("rm t_left_split.pdf")
    system("rm t_first_split-crop.pdf")
    system("rm t_right_split-crop.pdf")
    system("rm t_left_split-crop.pdf")
    system("rm t_tmp.pdf")
    system("rm t_tmp.tex")
    system("rm t_tmp.aux")
    system("rm t_tmp.log")
        
    pdf(file = "t_first_split.pdf", width = width_2, height = height_2)
    my_tree_plot(obj_branches,
                 terminal_panel = my_node_number2,
                 inner_panel = my_node_inner2,
                 tp_args = list(yunit = yunit_2))
    dev.off()

    terminal_plot_data <<- datatmp

    pdf(file = "t_left_split.pdf", width = width_rest, height = height_rest)
    my_tree_plot(obj[as.integer(upper_children)[1]],
                 inner_panel = my_node_inner2,
                 terminal_panel = my_node_bp5,
                 tp_args = list(size_x_text = size_x_text,
                                size_y_text = size_y_text,
                                angle_x_text = angle_x_text,
                                Ntot = nrow(obj$fitted),
                                sqrt_transf = sqrt_transf),
                 ip_args = list(width_inner = width_inner,
                                height_inner = height_inner
                                ))
    dev.off()
    
    pdf(file = "t_right_split.pdf", width = width_rest, height = height_rest)
    my_tree_plot(obj[as.integer(upper_children)[2]],
                 inner_panel = my_node_inner2,
                 terminal_panel = my_node_bp5,
                 tp_args = list(size_x_text = size_x_text,
                                size_y_text = size_y_text,
                                angle_x_text = angle_x_text,
                                Ntot = nrow(obj$fitted),
                                sqrt_transf = sqrt_transf),
                 ip_args = list(width_inner = width_inner,
                                height_inner = height_inner
                                ))
    dev.off()
    
    
    system("pdfcrop t_first_split.pdf")
    system("pdfcrop t_right_split.pdf")
    system("pdfcrop t_left_split.pdf")
    
    ## From idea in https://tex.stackexchange.com/a/140837
    thelatexfile <- paste0("
\\documentclass{article}
\\usepackage{graphicx}
\\usepackage{subfig}
\\begin{document}
\\thispagestyle{empty}
\\begin{figure}
\\begin{tabular}{c}
% (a) \\\\
\\hspace*{", hspace, "}
\\vspace*{", vspace, "}
\\includegraphics[width=20mm]{t_first_split-crop.pdf} \\\\[2pt]
% (b) \\\\[2pt]
\\includegraphics[width=150mm]{t_left_split-crop.pdf} \\\\[2pt]
% (c) \\\\[2pt]
\\includegraphics[width=150mm]{t_right_split-crop.pdf} \\\\
\\end{tabular}
\\caption*{",
main, "}
\\end{figure}
\\end{document}
")

    write(file = "t_tmp.tex", thelatexfile)
    system("texi2pdf t_tmp.tex")
    system("pdfcrop t_tmp.pdf")

    system("rm t_first_split.pdf")
    system("rm t_right_split.pdf")
    system("rm t_left_split.pdf")
    system("rm t_first_split-crop.pdf")
    system("rm t_right_split-crop.pdf")
    system("rm t_left_split-crop.pdf")
    system("rm t_tmp.pdf")
    system("rm t_tmp.tex")
    system("rm t_tmp.aux")
    system("rm t_tmp.log")
    system(paste("mv t_tmp-crop.pdf", paste0(" ", main, ".pdf")))
    rm(terminal_plot_data)
    rm(terminal_plot_data, envir = .GlobalEnv)
}




