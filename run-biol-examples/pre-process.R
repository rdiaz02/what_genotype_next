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



## Why do we allow now data with only a single column after using
## pre_process, in all the analysis scripts? Because we now are trying to
## also recover the Root -> something. If you know you are going to
## exclude those, as in the Baseline analyses, then data sets with only 1
## column make no sense: you only get "Root -> something"


merge_all_dups <- function(x) {
    ## Return a "column deduplicated" matrix
    while(TRUE) {
        dups <- which(duplicated(x, MARGIN = 2))
        if(length(dups))
            x <- merge_the_ident(x, dups[1])
        else(break)
    }
    return(x)
}

## In the future:
## use another character, not a "_". Something like "@@!" or something
## that will be unique, for whose pressence we test before starting, and
## that will also be used in the unfusing code.
merge_the_ident <- function(x, column) {
    ## Merges a column, column, identified as duplicated
    thedup <- x[, column]
    cnpre <- colnames(x)
    idt <- which( colSums(x == thedup) == nrow(x))
    idt <- setdiff(idt, column)
    cn <- colnames(x)[-c(idt, column)]
    x <- cbind(x[, -c(idt, column)], thedup)
    colnames(x) <- c(cn, paste(cnpre[c(column, idt)], collapse = "_"))
    return(x)
}


pre_process <- function(x, remove.constant, min.freq = 0.05,
                        max.cols = NULL) {
    ## As it says: preprocess with usual configuration:
    ## - no columns with tiny frequency (default is min of 0.05),
    ## - duplicated columns are merged
    ## - max.cols, in case we are using CBN et al.
    ## - if remove.constant, then no constant columns (those always present)
    ##          - those always absent are always removed
    ##          - No default for remove.constant: you must be explicit.
    if(!identical(sort(unique(as.vector(x))), c(0L, 1L)))
        stop("Values in x not in 0, 1")
    ## if(!isTRUE(all.equal(sort(unique(as.vector(x))), c(0L, 1L))))
    ##     stop("Values in x not in 0, 1")
    
    if(min.freq < 0) stop("min.freq has to be positive or 0")
    if(is.null(colnames(x))) stop("colnames must exist")
    nsubj <- nrow(x)
    cs1 <- colSums(x)
    ## Anything constant or with freq less than min.freq.
    if(remove.constant)
        rm1 <- which( (cs1 == 0) | (cs1 == nsubj) | (cs1 < (min.freq * nsubj)))
    else
        rm1 <- which( (cs1 == 0) | (cs1 < (min.freq * nsubj)))
    if(length(rm1))
        x <- x[, -rm1, drop = FALSE]
    ## Any duplicated cols
    x <- merge_all_dups(x)
    ## For CBN
    if(!is.null(max.cols)) {
        nc <- ncol(x)
        if(nc > max.cols) {
            keep <- order(colSums(x), decreasing = TRUE)[seq.int(max.cols)]
            x <- x[, keep, drop = FALSE]
        }
    }
    return(x)
}






library(codetools)
checkUsageEnv(env = .GlobalEnv)




