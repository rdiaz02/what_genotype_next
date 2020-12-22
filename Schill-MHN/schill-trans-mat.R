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


## Code to run Schill et al.'s MHN and obtain transition matrices between
## genotypes. Several implementations are provided, some using sparse
## matrices. In general, the best option is to use do_MHN2 for moderate to
## large numbers of features.


library(gtools)



## Use these wrappers
## observations (rows as patients, columns genes) -> transition matrix genotypes

##          Remember my transition matrices between genotypes have origin
##          in rows, destination in columns. Transposed w.r.t. to their
##          Figure 2, left.

do_MHN <- function(x,  lambda = 1/nrow(x)) {
    ## lambda 0.01 is what they use by default in their bioRxiv, p. 7 paper.
    ## In the paper it is 1/nrow(x). See paper and emails.
    ## What we did: for simulations we used lambda = 0.01.
    ## (see, e.g., function mhn_on_split, in file run-Schill-MHN-trans-mat.R)
    ## but for biological data we used lambda = 1/nrow(x).
    ## Differences are often minor. lmabda = 0.01 seems an easier setting
    ## for reproducibility with simulations since number of features
    ## can actually change in sampled data sets.
    ## See Suppl Mat.
    mhnd <- Data.to.pD(x)
    cat("\n      MHN: done Data.to.pD ", date(), "\n")
    theta <- Learn.MHN(mhnd, lambda = lambda)
    cat("\n      MHN: done Learn.MHN ", date(), "\n")
    colnames(theta) <- rownames(theta) <- colnames(x)
    trm <- theta_to_trans_rate_3(theta,
                                 inner_transition = inner_transitionRate_3_1)
    cat("\n      MHN: done theta_to_trans_rate_3 ", date(), "\n")

    return(list(
        theta = theta,
        transitionRateMatrix = trm,
        transitionMatrixTimeDiscretized =
            trans_rate_to_trans_mat(trm,
                                    method = "uniformization",
                                    paranoidCheck = TRUE),
        transitionMatrixCompExp =
            trans_rate_to_trans_mat(trm,
                                    method = "competingExponentials",
                                    paranoidCheck = TRUE)
                ))
}


## Identical to do_MHN, but with an argument for sparse, and corresponding
## additional code
do_MHN2 <- function(x,  lambda = 1/nrow(x), sparse = TRUE) {
    ## lambda 0.01 is what they use by default (see bioRxiv, p. 7)
    ## Nope, it is 1/nrow(x). See paper and emails.
    mhnd <- Data.to.pD(x)
    cat("\n      MHN: done Data.to.pD ", date(), "\n")
    theta <- Learn.MHN(mhnd, lambda = lambda)
    cat("\n      MHN: done Learn.MHN ", date(), "\n")
    colnames(theta) <- rownames(theta) <- colnames(x)
    if(!sparse) {
        trm <- theta_to_trans_rate_3(theta,
                                     inner_transition = inner_transitionRate_3_1)
        cat("\n      MHN: done theta_to_trans_rate_3 ", date(), "\n")
    } else {
        trm <- theta_to_trans_rate_3_SM(theta,
                                        inner_transition = inner_transitionRate_3_1)
        cat("\n      MHN: done theta_to_trans_rate_3_SM ", date(), "\n")
    }
    
    return(list(
        theta = theta,
        transitionRateMatrix = trm,
        transitionMatrixTimeDiscretized =
            trans_rate_to_trans_mat(trm,
                                    method = "uniformization",
                                    paranoidCheck = TRUE),
        transitionMatrixCompExp =
            trans_rate_to_trans_mat(trm,
                                    method = "competingExponentials",
                                    paranoidCheck = TRUE)
                ))
}




## transition rate matrix -> transition matrix
##    if method == uniformization, we assume no diagonal entry
##    because we compute it
##   similar code in time-discretized-CBN.R, but w/o ability to handle
##   sparse matrices

trans_rate_to_trans_mat <- function(x,
                                    method = c("competingExponentials",
                                               "uniformization"),
                                    paranoidCheck = TRUE) {
    method <- match.arg(method)
    if(method == "competingExponentials") {
        ## Conditional on there being a transition, thus diagonals are
        ## zero
        sx <- rowSums(x)
        if( inherits(x, "dgCMatrix") ) { ## sparse matrices
            tm <- x
            ii <- which(sx > 0)
            for(i in ii) {
                tm[i, ] <- tm[i, ]/sx[i]
            }
        } else {
            tm <- sweep(x, 1, sx, "/")
            tm[nrow(tm), ] <- 0 ## last row is 0
        }
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
        ## In p. 243 they say
        ## "diagonal elements are defined as Qxx = -\sum y;x so that
        ## columns sum to zero"
        ## They talk of columns, because they transpose it relative
        ## to what I do
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


## theta from Learn.MHN -> transition rate matrix
##      Note that the diagonal is not added
##         This computes the products of the Theta, as in Fig.2 right
theta_to_trans_rate_1 <- function(theta) {
    Theta <- exp(theta)
    geneNames <- colnames(theta)

    k <- ncol(theta)
    genots <- c(list(rep(0, k)), allGenotypes_former(k))
    numGenots <- length(genots)
    genotNames <- unlist(
        lapply(genots,
               function(x)
                   paste(geneNames[which(x == 1L)], sep = "", collapse = ", "))
        )
    genotNames[genotNames == ""] <- "WT"
    ## this single call is the one that takes most time
    ## t1 <- Sys.time()
    ## These calls take tiny time
    dim1 <- rep(genots, numGenots)
    dim2 <- rep(genots, rep(numGenots, numGenots))
    ## t2 <- Sys.time()
    tr2 <- function(g1, g2) transitionRateB(g1, g2, Theta)

    TRM <- mapply(tr2,
                  dim1,
                  dim2)
    ## TRM <- mapply(function(g1, g2) transitionRate(g1, g2, Theta),
    ##               dim1,
    ##               dim2)

    ## TRM <- mapply(function(g1, g2) transitionRate(g1, g2, Theta),
    ##            rep(genots, numGenots),
    ##            rep(genots, rep(numGenots, numGenots)))

    TRM <- matrix(TRM, ncol = numGenots, byrow = FALSE)
    colnames(TRM) <- rownames(TRM) <- genotNames
    ## t3 <- Sys.time()

    ## cat("\nt2 - t1 :", t2 - t1)
    ## cat("\nt3 - t2 :", t3 - t2)
    ## cat("\n")
    return(TRM)
}



## integer (number of genes) -> all genotypes as 0,1 vectors
allGenotypes_2 <- function(k) {
    ## From OncoSimulR
    f1 <- function(n) {
        lapply(seq.int(n), function(x) combinations(n = n, r = x))}

    list.of.vectors <- function(y) {
        ## there's got to be a simpler way
        lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}),
                      recursive = FALSE),
               function(m) m[[1]])
    }
   
    mutated <- list.of.vectors(f1(k))
    num_mutated <- lapply(mutated, length)
    
    ## number of genes, mutated positions -> binary genotype as vector of
    ## 0, 1
    binary_genotype <- function(x, k) {
        y <- rep(0L, k)
        y[x] <- 1L
        return(y)
    }
    bin_genot <- lapply(mutated, function(x) binary_genotype(x, k = k))
    ## Map(function(nm, m) list(num_mutated = nm, genot = m), num_mutated, bin_genot)
    df <- data.frame(num_mut = unlist(num_mutated),
                     mutated = I(mutated),  ## not needed?
                     bin_genotype = I(bin_genot))
    ## Add WT
    df <- rbind(data.frame(num_mut = 0,
                           mutated = NA,
                           bin_genotype = I(list(rep(0, k)))),
                df)
    return(df)
}




## integer (number of genes) -> all genotypes as 0,1 vectors
allGenotypes_3 <- function(k) {
    ## From OncoSimulR
    f1 <- function(n) {
        lapply(seq.int(n), function(x) combinations(n = n, r = x))}

    list.of.vectors <- function(y) {
        ## there's got to be a simpler way
        lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}),
                      recursive = FALSE),
               function(m) m[[1]])
    }
   
    mutated <- list.of.vectors(f1(k))
    num_mutated <- lapply(mutated, length)
    
    ## number of genes, mutated positions -> binary genotype as vector of
    ## 0, 1
    binary_genotype <- function(x, k) {
        y <- rep(0L, k)
        y[x] <- 1L
        return(y)
    }
    bin_genot <- lapply(mutated, function(x) binary_genotype(x, k = k))
    
    return(list(num_mut = c(0, unlist(num_mutated)),
                mutated = c(list(NA), mutated),
                bin_genotype = c(list(rep(0, k)), bin_genot)))

}





## ## integer (number of genes) -> all genotypes as 0,1 vectors
## allGenotypes_1 <- function(k) {
##     ## From OncoSimulR
##     f1 <- function(n) {
##         lapply(seq.int(n), function(x) combinations(n = n, r = x))}

##     list.of.vectors <- function(y) {
##         ## there's got to be a simpler way
##         lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}),
##                       recursive = FALSE),
##                function(m) m[[1]])
##     }
   
##     mutated <- list.of.vectors(f1(k))
##     num_mutated <- lapply(mutated, length)
    
##     ## number of genes, mutated positions -> binary genotype as vector of
##     ## 0, 1
##     binary_genotype <- function(x, k) {
##         y <- rep(0L, k)
##         y[x] <- 1L
##         return(y)
##     }
##     bin_genot <- lapply(mutated, function(x) binary_genotype(x, k = k))
##     Map(function(nm, m) list(num_mutated = nm, genot = m), num_mutated, bin_genot)
##     ## FIXME: add WT
## }


## integer (number of genes) -> all genotypes as 0,1 vectors
allGenotypes_former <- function(k) {
    ## From OncoSimulR
    f1 <- function(n) {
        lapply(seq.int(n), function(x) combinations(n = n, r = x))}


    list.of.vectors <- function(y) {
        ## there's got to be a simpler way
        lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}),
                      recursive = FALSE),
               function(m) m[[1]])
    }
   
    mutated <- list.of.vectors(f1(k))
    
    ## number of genes, mutated positions -> binary genotype as vector of
    ## 0, 1
    binary_genotype <- function(x, k) {
        y <- rep(0L, k)
        y[x] <- 1L
        return(y)
    }
    lapply(mutated, function(x) binary_genotype(x, k = k))
}


## ## two genotypes as vector of 0,1 : can we transition from x -> y?
## ## only if exactly a difference of 1
## canTransition <- function(x, y) {
##     sum(y - x) == 1
## }


## genotype, genotype, Theta (as exp(theta)) -> transition rate x -> y
transitionRate <- function(x, y, Theta) {
    if(sum(y) != (sum(x) + 1) ) {
        return(0)
    } else {
        posy <- which(y != x)
        if(length(posy) != 1) {
            return(0)
        } else {
            posx <- which(x == 1L)
            if(length(posx) == 0) {
                ret <- return(Theta[posy, posy])
            } else {
                ret <- (Theta[posy, posy] * prod(Theta[posy, posx]))
            }
            if(length(ret) > 1) {
                cat("\n here")
                stop()
                browser()
            }
            else (return(ret))
            ## if(length(posx) == 0) return(Theta[posy, posy])
            ## else return(Theta[posy, posy] * cumprod(Theta[posy, posx]))
        }
    }
}

## this is faster
## genotype, genotype, Theta (as exp(theta)) -> transition rate x -> y
transitionRateB <- function(x, y, Theta) {
    if( (sum(y) != (sum(x) + 1)) || (sum(y != x) != 1) ) {        
        return(0)
    } else {
        posy <- which(y != x)
        posx <- which(x == 1L)
        if(length(posx) == 0) {
            return(Theta[posy, posy])
        } else {
            return(Theta[posy, posy] * prod(Theta[posy, posx]))
        }
    }
}


## this is faster
## row number of genotype, genotypes data frame
##            Theta (as exp(theta))
##            maximum number of mutations -> transition rate x -> y
transitionRateC <- function(i, genotypes,  Theta, maxmut,
                            inner_transition) {

    num_genots <- nrow(genotypes)
    ## Last genotype
    if(i == num_genots) return(rep(0, num_genots)) 

    ## t11 <- Sys.time()
    ## All genotypes we are sure we cannot transition to
    nmuts <- genotypes[i, "num_mut"]
    ## Necessarily same of fewer mutations
    tmp1 <- rep(0, length.out = (i - 1))
    ## Necessarily not reachable
    ## t12 <- Sys.time()
    if( (nmuts + 2) <= maxmut ) {
        mi2 <- match(nmuts + 2, genotypes$num_mut)
        tmp3 <- rep(0, length.out = num_genots - mi2 + 1)
        upper_gg <- (mi2 - 1)
    } else {
        tmp3 <- double(0)
        upper_gg <- num_genots
    }
    ## t13 <- Sys.time()
    ## All those with same number of mutations are also
    ## necessarily not reachable.
    mi1 <- match(nmuts + 1, genotypes$num_mut)
    tmp2 <- rep(0, length.out = mi1 - i)
    ##  t14 <- Sys.time()
    ## The rest
    ## gg <- genotypes[c(i, mi1:upper_gg), ]

    ## t15 <- Sys.time()

    ## ff <- function(x) inner_transition(1, x, gg, Theta)
    ## qs <- vapply(2:nrow(gg),
    ##              ff,
    ##              double(1)
    ##              )

    ff <- function(x) inner_transition(i, x, genotypes, Theta)
    qs <- vapply(mi1:upper_gg,
                 ff,
                 double(1)
                 )
    ## t16 <- Sys.time()
    ## alltimes <- c(t12 - t11, t13 - t12, t14 - t13, t15 - t14, t16 - t15)
    ## cat("\n    all times: ", alltimes, "\n")
    return(c(tmp1, tmp2, qs, tmp3))
}



## this is faster
## row number of genotype, genotypes data frame
##            Theta (as exp(theta))
##            maximum number of mutations -> transition rate x -> y
transitionRateC3 <- function(i, genotypes,  Theta, maxmut,
                            inner_transition) {

    num_genots <- length(genotypes$num_mut)
    ## Last genotype
    if(i == num_genots) return(rep(0, num_genots)) 

    ## t11 <- Sys.time()
    ## All genotypes we are sure we cannot transition to
    nmuts <- genotypes$num_mut[i]
    ## Necessarily same of fewer mutations
    tmp1 <- rep(0, length.out = (i - 1))
    ## Necessarily not reachable
    ## t12 <- Sys.time()
    if( (nmuts + 2) <= maxmut ) {
        mi2 <- match(nmuts + 2, genotypes$num_mut)
        tmp3 <- rep(0, length.out = num_genots - mi2 + 1)
        upper_gg <- (mi2 - 1)
    } else {
        tmp3 <- double(0)
        upper_gg <- num_genots
    }
    ## t13 <- Sys.time()
    ## All those with same number of mutations are also
    ## necessarily not reachable.
    mi1 <- match(nmuts + 1, genotypes$num_mut)
    tmp2 <- rep(0, length.out = mi1 - i)
    ##  t14 <- Sys.time()
    ## The rest
    ## gg <- genotypes[c(i, mi1:upper_gg), ]

    ## t15 <- Sys.time()

    ## ff <- function(x) inner_transition(1, x, gg, Theta)
    ## qs <- vapply(2:nrow(gg),
    ##              ff,
    ##              double(1)
    ##              )

    ff <- function(x) inner_transition(i, x, genotypes, Theta)
    qs <- vapply(mi1:upper_gg,
                 ff,
                 double(1)
                 )
    ## t16 <- Sys.time()
    ## alltimes <- c(t12 - t11, t13 - t12, t14 - t13, t15 - t14, t16 - t15)
    ## cat("\n    all times: ", alltimes, "\n")
    return(c(tmp1, tmp2, qs, tmp3))
}



## row number of genotype, genotypes data frame
##            Theta (as exp(theta))
##            maximum number of mutations ->
##      transition rate x -> y [only those that are not 0 by construction]
##   output to be used for populating a sparse matrix
transitionRateC3_SM <- function(i, genotypes,  Theta, maxmut,
                            inner_transition) {

    num_genots <- length(genotypes$num_mut)
    ## Last genotype
    ## But this is silly
    if(i == num_genots) return(cbind(j = num_genots, x = 0))

    ## t11 <- Sys.time()
    ## All genotypes we are sure we cannot transition to
    nmuts <- genotypes$num_mut[i]
    ## Necessarily same of fewer mutations
    ## Necessarily not reachable
    if( (nmuts + 2) <= maxmut ) {
        mi2 <- match(nmuts + 2, genotypes$num_mut)
        upper_gg <- (mi2 - 1)
    } else {
        upper_gg <- num_genots
    }

    ## All those with same number of mutations are also
    ## necessarily not reachable.
    mi1 <- match(nmuts + 1, genotypes$num_mut)

    ff <- function(x) inner_transition(i, x, genotypes, Theta)
    qs <- vapply(mi1:upper_gg,
                 ff,
                 double(1)
                 )
    return(cbind(j = mi1:upper_gg, 
                 x = qs))
}


## genotype position (row), genotype position,
##    genotype data frame
##    Theta (as exp(theta)) -> transition rate x -> y
inner_transitionRate_3_1 <- function(i, j, genotypes, Theta) {
    
    if( genotypes$num_mut[j] != (genotypes$num_mut[i] + 1) ) return(0)

    x <- genotypes$bin_genotype[[i]]
    y <- genotypes$bin_genotype[[j]]

    if( length(posy <- which(y != x)) != 1 ) return(0)

    posx <- which(x == 1L)
    if(length(posx) == 0) return(Theta[posy, posy])

    return(Theta[posy, posy] * prod(Theta[posy, posx]))

}

## genotype position (row), genotype position,
##    genotype data frame
##    Theta (as exp(theta)) -> transition rate x -> y
inner_transitionRate_1 <- function(i, j, genotypes, Theta) {

    x <- genotypes[[i, "bin_genotype"]]
    y <- genotypes[[j, "bin_genotype"]]
    
    if( (genotypes[j, "num_mut"] != (genotypes[i, "num_mut"] + 1) ) ||
        (length(posy <- which(y != x)) != 1) ){
        ## (sum(y != x) != 1) ) {
        return(0)
    } else {
        ## posy <- which(y != x)
        posx <- which(x == 1L)
        if(length(posx) == 0) {
            return(Theta[posy, posy])
        } else {
            return(Theta[posy, posy] * prod(Theta[posy, posx]))
        }
    }
}


## genotype position (row), genotype position,
##    genotype data frame
##    Theta (as exp(theta)) -> transition rate x -> y
inner_transitionRate_3_2 <- function(i, j, genotypes, Theta) {

    if( genotypes$num_mut[j] != (genotypes$num_mut[i] + 1) ) return(0)

    x <- genotypes$mutated[[i]]
    y <- genotypes$mutated[[j]]

    if( length(posy <- setdiff(y, x)) != 1) return(0)

    if(all(is.na(x))) return(Theta[posy, posy])

    return(Theta[posy, posy] * prod(Theta[posy, x]))

}


## genotype position (row), genotype position,
##    genotype data frame
##    Theta (as exp(theta)) -> transition rate x -> y
inner_transitionRate_2 <- function(i, j, genotypes, Theta) {

    if(genotypes[j, "num_mut"] != (genotypes[i, "num_mut"] + 1) ) return(0)

    x <- genotypes[[i, "mutated"]]
    y <- genotypes[[j, "mutated"]]

    if( length(posy <- setdiff(y, x)) != 1) return(0)

    if(all(is.na(x))) return(Theta[posy, posy])

    return(Theta[posy, posy] * prod(Theta[posy, x]))

}



## ## genotype position (row), genotype position,
## ##    genotype data frame
## ##    Theta (as exp(theta)) -> transition rate x -> y
## inner_transitionRate_2 <- function(i, j, genotypes, Theta) {

##     x <- genotypes[[i, "mutated"]]
##     y <- genotypes[[j, "mutated"]]
    
##     if( (genotypes[j, "num_mut"] != (genotypes[i, "num_mut"] + 1) ) ||
##         (length(posy <- setdiff(y, x)) != 1) ) {
##         ## (length(setdiff(y, x)) != 1) ) {
##         return(0)
##     } else {
##         ## posy <- setdiff(y, x)
##         if(all(is.na(x))) { ## x is WT
##             return(Theta[posy, posy])
##         } else {
##             return(Theta[posy, posy] * prod(Theta[posy, x]))
##         }
##     }
## }




## theta from Learn.MHN
##      function used to compute theta -> transition rate matrix
theta_to_trans_rate <- function(theta,
                                inner_transition = inner_transitionRate_1) {

    ## t1 <- Sys.time()
    Theta <- exp(theta)
    geneNames <- colnames(theta)
    
    k <- ncol(theta)
    genots <- allGenotypes_2(k)
    numGenots <- nrow(genots)

    genotNames <- unlist(
        lapply(genots$bin_genotype,
               function(x)
                   paste(geneNames[which(x == 1L)], sep = "", collapse = ", "))
        )
    genotNames[genotNames == ""] <- "WT"

    ## This single call is the one that takes most time
    ## t2 <- Sys.time()

    TRM <- vapply(seq.int(numGenots),
           function(x)
               transitionRateC(x, genotypes = genots,
                               Theta = Theta, maxmut = k,
                               inner_transition = inner_transition),
           double(numGenots))
    TRM <- matrix(TRM, ncol = numGenots, byrow = TRUE)
    colnames(TRM) <- rownames(TRM) <- genotNames
    ## t3 <- Sys.time()

    ## cat("\nt2 - t1 :", t2 - t1)
    ## cat("\nt3 - t2 :", t3 - t2)
    ## cat("\n")
    return(TRM)
}



## theta from Learn.MHN
##      function used to compute theta -> transition rate matrix
##      Note that the diagonal is not added
##         This computes the products of the Theta, as in Fig.2 right
theta_to_trans_rate_3 <- function(theta,
                                inner_transition = inner_transitionRate_3_1) {

    ## t1 <- Sys.time()
    Theta <- exp(theta)
    geneNames <- colnames(theta)
    
    k <- ncol(theta)
    genots <- allGenotypes_3(k)
    numGenots <- length(genots$num_mut)

    genotNames <- unlist(
        lapply(genots$bin_genotype,
               function(x)
                   paste(geneNames[which(x == 1L)], sep = "", collapse = ", "))
        )
    genotNames[genotNames == ""] <- "WT"

    ## This single call is the one that takes most time
    ## t2 <- Sys.time()

    TRM <- vapply(seq.int(numGenots),
           function(x)
               transitionRateC3(x, genotypes = genots,
                               Theta = Theta, maxmut = k,
                               inner_transition = inner_transition),
           double(numGenots))
    TRM <- matrix(TRM, ncol = numGenots, byrow = TRUE)
    colnames(TRM) <- rownames(TRM) <- genotNames
    ## t3 <- Sys.time()

    ## cat("\nt2 - t1 :", t2 - t1)
    ## cat("\nt3 - t2 :", t3 - t2)
    ## cat("\n")
    return(TRM)
}



## theta from Learn.MHN
##      function used to compute trans rate matrix -> transition rate matrix
theta_to_trans_rate_3_SM <- function(theta,
                                     inner_transition = inner_transitionRate_3_1) {

    ## t1 <- Sys.time()
    Theta <- exp(theta)
    geneNames <- colnames(theta)
    
    k <- ncol(theta)
    genots <- allGenotypes_3(k)
    numGenots <- length(genots$num_mut)

    genotNames <- unlist(
        lapply(genots$bin_genotype,
               function(x)
                   paste(geneNames[which(x == 1L)], sep = "", collapse = ", "))
        )
    genotNames[genotNames == ""] <- "WT"

    ## Initialize sparseMatrix in first call
    trmv <- transitionRateC3_SM(1, genotypes = genots,
                                Theta = Theta, maxmut = k,
                                inner_transition = inner_transition)
    TRM <- sparseMatrix(i = rep(1, length(trmv[, "j"])),
                        j = trmv[, "j"],
                        x = trmv[, "x"],
                        dims = c(numGenots, numGenots),
                        dimnames = list(genotNames, genotNames))
    ## Can skip last one
    for(i in seq.int(2, numGenots - 1)) {
        trmv <- transitionRateC3_SM(i, genotypes = genots,
                                    Theta = Theta, maxmut = k,
                                    inner_transition = inner_transition)
        TRM[i = i, j = trmv[, "j"]] <- trmv[, "x"]
    }
    return(TRM)
}

######################################################################
######################################################################
local({

    ## Minimal tests
    
    ## from rounding the first 4 x 4 entries of Schill's
    ## cancer example
    me1 <- structure(c(0.5, 0, 0, 0,
                       0, -1.2, 0, 0,
                       0, 0, -1.1, 0,
                       -0.1, 0, 0.3, 0.1),
                     .Dim = c(4L, 4L),
                     .Dimnames =
                         list(c("A", "B", "C", "D"),
                              c("A", "B", "C", "D")))

    eme1 <- exp(me1)

    ## off-diagonal terms first
    trm0 <- matrix(0, ncol = 16, nrow = 16)

    geneNames <- colnames(eme1)
    k <- 4
    genots <- c(list(rep(0, k)), allGenotypes_former(k))
    genotNames <- unlist(
        lapply(genots,
               function(x)
                   paste(geneNames[which(x == 1L)], sep = "", collapse = ", "))
    )
    genotNames[genotNames == ""] <- "WT"

    colnames(trm0) <- rownames(trm0) <- genotNames


    trm0["WT", "A"] <- eme1["A", "A"]
    trm0["WT", "B"] <- eme1["B", "B"]
    trm0["WT", "C"] <- eme1["C", "C"]
    trm0["WT", "D"] <- eme1["D", "D"]    

    trm0["A", "A, B"] <- eme1["B", "B"] * eme1["B", "A"]
    trm0["A", "A, C"] <- eme1["C", "C"] * eme1["C", "A"]
    trm0["A", "A, D"] <- eme1["D", "D"] * eme1["D", "A"]

    trm0["B", "A, B"] <- eme1["A", "A"] * eme1["A", "B"]
    trm0["B", "B, C"] <- eme1["C", "C"] * eme1["C", "B"]
    trm0["B", "B, D"] <- eme1["D", "D"] * eme1["D", "B"]

    trm0["C", "A, C"] <- eme1["A", "A"] * eme1["A", "C"]
    trm0["C", "B, C"] <- eme1["B", "B"] * eme1["B", "C"]
    trm0["C", "C, D"] <- eme1["D", "D"] * eme1["D", "C"]
    
    trm0["D", "A, D"] <- eme1["A", "A"] * eme1["A", "D"]
    trm0["D", "B, D"] <- eme1["B", "B"] * eme1["B", "D"]
    trm0["D", "C, D"] <- eme1["C", "C"] * eme1["C", "D"]

    trm0["A, B", "A, B, C"] <- eme1["C", "C"] * eme1["C", "A"] * eme1["C", "B"]
    trm0["A, B", "A, B, D"] <- eme1["D", "D"] * eme1["D", "A"] * eme1["D", "B"]

    trm0["A, C", "A, B, C"] <- eme1["B", "B"] * eme1["B", "A"] * eme1["B", "C"]
    trm0["A, C", "A, C, D"] <- eme1["D", "D"] * eme1["D", "A"] * eme1["D", "C"]

    trm0["A, D", "A, B, D"] <- eme1["B", "B"] * eme1["B", "A"] * eme1["B", "D"]
    trm0["A, D", "A, C, D"] <- eme1["C", "C"] * eme1["C", "A"] * eme1["C", "D"]

    trm0["B, C", "A, B, C"] <- eme1["A", "A"] * eme1["A", "B"] * eme1["A", "C"]
    trm0["B, C", "B, C, D"] <- eme1["D", "D"] * eme1["D", "B"] * eme1["D", "C"]    
    
    trm0["B, D", "A, B, D"] <- eme1["A", "A"] * eme1["A", "B"] * eme1["A", "D"]
    trm0["B, D", "B, C, D"] <- eme1["C", "C"] * eme1["C", "B"] * eme1["C", "D"]
    
    trm0["C, D", "A, C, D"] <- eme1["A", "A"] * eme1["A", "C"] * eme1["A", "D"]
    trm0["C, D", "B, C, D"] <- eme1["B", "B"] * eme1["B", "C"] * eme1["B", "D"]

    trm0["A, B, C", "A, B, C, D"] <-
        eme1["D", "D"] * eme1["D", "A"] * eme1["D", "B"] * eme1["D", "C"]

    trm0["A, B, D", "A, B, C, D"] <-
        eme1["C", "C"] * eme1["C", "A"] * eme1["C", "B"] * eme1["C", "D"]

    trm0["A, C, D", "A, B, C, D"] <-
        eme1["B", "B"] * eme1["B", "A"] * eme1["B", "C"] * eme1["B", "D"]

    trm0["B, C, D", "A, B, C, D"] <-
        eme1["A", "A"] * eme1["A", "B"] * eme1["A", "C"] * eme1["A", "D"]


    ## Adding the diagonal
    ## Q_ii = - \Sum Q_ji
    trm1 <- trm0
    diag(trm1) <- -1 * rowSums(trm0)

    ## paranoid check
    stopifnot(isTRUE(all.equal(rep(0, 16),
                     rowSums(trm1),
                     check.attributes = FALSE)))


    ## Verify that transition rates between different genotypes are the
    ## same from time discretized and using the competing exponentials
    ## approach

    
    ## Uniformization
    gamma <- max(abs(diag(trm1)))
    time_discr_trm <- diag(16) + trm1/gamma

    trans_mat_time_discr <- time_discr_trm
    diag(trans_mat_time_discr) <- 0
    s_tmtd <- rowSums(trans_mat_time_discr)
    trans_mat_time_discr <- sweep(trans_mat_time_discr, 1, s_tmtd, "/")
    ## note: last row is NaN
    trans_mat_time_discr[16, ] <- 0
    ## paranoid check
    stopifnot(isTRUE(all.equal(c(rep(1, 15), 0),
                               rowSums(trans_mat_time_discr),
                               check.attributes = FALSE)))


    ## Competing exponentials
    trans_mat_comp_exp <- trm0
    s_tmce <- rowSums(trans_mat_comp_exp)
    trans_mat_comp_exp <- sweep(trans_mat_comp_exp, 1, s_tmce, "/")
    trans_mat_comp_exp[16, ] <- 0 ## last row is 0
    ## paranoid check
    stopifnot(isTRUE(all.equal(c(rep(1, 15), 0),
                               rowSums(trans_mat_comp_exp),
                               check.attributes = FALSE)))

    ## Verify identity, as it should be
    stopifnot(isTRUE(all.equal(trans_mat_comp_exp,
                               trans_mat_time_discr)))

    ## For us, now, it is faster and simpler if we just use the competing
    ## exponentials one.

    tm_u <- trans_rate_to_trans_mat(trm0, method = "uniformization")
    tm_ce <- trans_rate_to_trans_mat(trm0, method = "competingExponentials")

    stopifnot(tm_u[nrow(tm_u), ncol(tm_u)] == 1)
    stopifnot(isTRUE(all.equal(tm_u, time_discr_trm)))
    stopifnot(isTRUE(all.equal(tm_ce, trans_mat_comp_exp)))
    
    ## Verify code to obtain transition rate matrix
    ttr_1 <- theta_to_trans_rate_1(me1)
    ttr_2 <- theta_to_trans_rate(me1, inner_transition = inner_transitionRate_1)
    ttr_22 <- theta_to_trans_rate(me1, inner_transition = inner_transitionRate_2)
    ttr_31 <- theta_to_trans_rate_3(me1, inner_transition = inner_transitionRate_3_1)
    ttr_32 <- theta_to_trans_rate_3(me1, inner_transition = inner_transitionRate_3_2)            

    stopifnot(identical(ttr_1, ttr_2))
    stopifnot(identical(ttr_2, ttr_22))
    stopifnot(identical(ttr_2, ttr_31))
    stopifnot(identical(ttr_2, ttr_32))
    
    stopifnot(identical(ttr_1, trm0))
    stopifnot(identical(ttr_2, trm0))
    stopifnot(identical(ttr_31, trm0))
    stopifnot(identical(ttr_32, trm0))

    ## In particular, note the asymmetrical changes
    ## Figure 2, right
    ## From D mutated to both C and D mutated
    ## theta_33 * theta_34
    stopifnot(all.equal(eme1[3, 3] * eme1[3, 4], trm0["D", "C, D"]))
    ## This is different, of course
    stopifnot(eme1[3, 3] * eme1[4, 3] != trm0["D", "C, D"])
    
})

######################################################################
######################################################################






######################################################################
## From Schill's ExampleApplications.R

pwd <- getwd()
setwd("./MHN")
source("UtilityFunctions.R")
source("ModelConstruction.R")
source("Likelihood.R")
source("RegularizedOptimization.R")
setwd(pwd)
rm(pwd)


## additional tests and example calls
local({

    pwd2 <- getwd()
    setwd("./MHN")

    ## Example of calls from original Schill code
    ##Cancer Progression Data----------------
    
    Dat <- readRDS(file="data/BreastCancer.rds") 

    pD <- Data.to.pD(Dat)
    Theta.BC <- Learn.MHN(pD, lambda=0.01)

    colnames(Theta.BC) <- colnames(Dat)
    rownames(Theta.BC) <- colnames(Theta.BC)

######################################################################
    setwd(pwd2)
    rm(pwd2)

    ## smaller:
    e1 <- Theta.BC[1:4, 1:4]

    system.time(te1_1 <- theta_to_trans_rate_1(e1))
    system.time(te1_2 <- theta_to_trans_rate(e1))
    system.time(te1_3 <- theta_to_trans_rate_3(e1))

    identical(te1_1, te1_2)
    identical(te1_1, te1_3)


    tt31 <- theta_to_trans_rate_3(Theta.BC, inner_transition = inner_transitionRate_3_1)

    ## slow, and the faster one is slower!!
    system.time(tt1 <- theta_to_trans_rate_1(Theta.BC))
    ## these two are very slow
    ## system.time(tt2 <- theta_to_trans_rate(Theta.BC))
    ## system.time(tt22 <- theta_to_trans_rate(Theta.BC, inner_transition = inner_transitionRate_2))
    ## but see this! Culprit was accessing a data frame?
    system.time(tt31 <- theta_to_trans_rate_3(Theta.BC, inner_transition = inner_transitionRate_3_1))
    ## but this is slower by x2. Why? setdiff, I guess
    system.time(tt32 <- theta_to_trans_rate_3(Theta.BC, inner_transition = inner_transitionRate_3_2))
    ## stopifnot(identical(tt1, tt2))
    ## stopifnot(identical(tt1, tt22))
    stopifnot(identical(tt1, tt31))
    stopifnot(identical(tt31, tt32))

    ## If profiling
    ## library(proftools)
    ## https://cran.r-project.org/web/packages/proftools/vignettes/proftools.pdf

    
    ## pd <- profileExpr(theta_to_trans_rate(Theta.BC, inner_transition = inner_transitionRate_2))
    ## pd1 <- profileExpr(theta_to_trans_rate(Theta.BC, inner_transition = inner_transitionRate_1))
    ## funSummary(pd)
    ## callSummary(pd)
    
    ## funSummary(pd1)
    ## callSummary(pd1)



theta0 <- structure(c(0.26, -0.57, 0, -0.97, -0.43, -1.18, 0, -0.29, 0, 
0, -1.84, 0, -0.23, -0.5, -0.04, -0.24), .Dim = c(4L, 4L), .Dimnames = list(
    c("A", "B", "C", "D"), c("A", "B", "C", "D")))


trm0 <- theta_to_trans_rate_3(theta0,
                              inner_transition = inner_transitionRate_3_1)
trmSM <- theta_to_trans_rate_3_SM(theta0,
                                  inner_transition = inner_transitionRate_3_1)
trmSMm <- as.matrix(trmSM)

stopifnot(identical(trm0, trmSMm))


## Use all four data sets
Dat1 <- readRDS(file="./MHN/data/BreastCancer.rds") [1:50, 1:4]
colnames(Dat1) <- LETTERS[1:ncol(Dat1)]
Dat1 <- as.matrix(Dat1)

mhn0 <- do_MHN(Dat1)
mhnSM <- do_MHN2(Dat1)

stopifnot(identical(mhn0$theta, mhnSM$theta))
stopifnot(identical(mhn0$transitionRateMatrix, as.matrix(mhnSM$transitionRateMatrix)))
stopifnot(isTRUE(all.equal(mhn0$transitionMatrixTimeDiscretized,
                    as.matrix(mhnSM$transitionMatrixTimeDiscretized))))
stopifnot(isTRUE(all.equal(mhn0$transitionMatrixCompExp,
                    as.matrix(mhnSM$transitionMatrixCompExp))))




Dat1 <- readRDS(file="./MHN/data/ColorectalCancer.rds")[1:300, 1:9] 
# colnames(Dat1) <- LETTERS[1:ncol(Dat1)]
Dat1 <- as.matrix(Dat1)

mhn0 <- do_MHN(Dat1)
mhnSM <- do_MHN2(Dat1)

stopifnot(identical(mhn0$theta, mhnSM$theta))
stopifnot(identical(mhn0$transitionRateMatrix, as.matrix(mhnSM$transitionRateMatrix)))
stopifnot(isTRUE(all.equal(mhn0$transitionMatrixTimeDiscretized,
                    as.matrix(mhnSM$transitionMatrixTimeDiscretized))))
stopifnot(isTRUE(all.equal(mhn0$transitionMatrixCompExp,
                    as.matrix(mhnSM$transitionMatrixCompExp))))


Dat1 <- readRDS(file="./MHN/data/RenalCellCarcinoma.rds")[, 1:9]
# colnames(Dat1) <- LETTERS[1:ncol(Dat1)]
Dat1 <- as.matrix(Dat1)

mhn0 <- do_MHN(Dat1)
mhnSM <- do_MHN2(Dat1)

stopifnot(identical(mhn0$theta, mhnSM$theta))
stopifnot(identical(mhn0$transitionRateMatrix, as.matrix(mhnSM$transitionRateMatrix)))
stopifnot(isTRUE(all.equal(mhn0$transitionMatrixTimeDiscretized,
                    as.matrix(mhnSM$transitionMatrixTimeDiscretized))))
stopifnot(isTRUE(all.equal(mhn0$transitionMatrixCompExp,
                    as.matrix(mhnSM$transitionMatrixCompExp))))


Dat11 <- readRDS(file="./MHN/data/Glioblastoma.rds")


for(i in 1:10) {
    Dat1 <- Dat11[, sample(1:ncol(Dat11), 6)]
    Dat1 <- as.matrix(Dat1)
    mhn0 <- do_MHN(Dat1)
    mhnSM <- do_MHN2(Dat1)
    stopifnot(identical(mhn0$theta, mhnSM$theta))
    stopifnot(identical(mhn0$transitionRateMatrix, as.matrix(mhnSM$transitionRateMatrix)))
    stopifnot(isTRUE(all.equal(mhn0$transitionMatrixTimeDiscretized,
                               as.matrix(mhnSM$transitionMatrixTimeDiscretized))))
    stopifnot(isTRUE(all.equal(mhn0$transitionMatrixCompExp,
                               as.matrix(mhnSM$transitionMatrixCompExp))))
}


## tripwire to check tests run
## stopifnot(3 == 2)

})



## A couple of paranoid checks that we are using the indices of theta
## correctly

local({
    library(OncoSimulR)
    ## 
    N <- 200
    na <- N
    nc <- N + 21 ## round( 10 * runif(1))
    nab <- N + 5 ## round( 10 * runif(1)) + 
    nbc <- N + 7 + round( 10 * runif(1))
    n00 <- N/10 + round( 10 * runif(1))
    dB <- matrix(
        c(
            rep(c(1, 0, 0, 0), na) 
          , rep(c(0, 0, 1, 0), nc)
      , rep(c(1, 1, 0, 0), nab)
      , rep(c(0, 1, 1, 0), nbc)
      , rep(c(0, 0, 0, 0), n00)
        ), ncol = 4, byrow = TRUE
    )
    colnames(dB) <- LETTERS[1:4]
    sampledGenotypes(dB)
    try(do_MHN(dB), silent = TRUE) ## this crashes. Of course: don't do silly things like passing
    ## a data frame with one or more colSums == 0
    dB2 <- dB[, -4]
    mm1 <- do_MHN(dB2)   
    ## From A to A, B:
    mm1$transitionMatrixCompExp["A", "A, B"]
    mm1$transitionRateMatrix["A", "A, B"]
    exp(mm1$theta["B", "B"]) * exp(mm1$theta["B", "A"])
    
    ## same, of course
    stopifnot(all.equal( mm1$transitionRateMatrix["A", "A, B"],
                        exp(mm1$theta["B", "B"]) * exp(mm1$theta["B", "A"])))

    ## not the same, though the first occasionally is identical if
    ## symmetry and thetas are identical A, B and B, A, if na = nab. That
    ## is why we make sure different N above
    ## stopifnot(mm1$transitionRateMatrix["A", "A, B"] != 
    ##           exp(mm1$theta["B", "B"]) * exp(mm1$theta["A", "B"]))

    ## stopifnot(mm1$transitionRateMatrix["A", "A, B"] != 
    ##           exp(mm1$theta["A", "A"]) * exp(mm1$theta["A", "B"]))

    ## stopifnot(mm1$transitionRateMatrix["A", "A, B"] != 
    ##           exp(mm1$theta["A", "A"]) * exp(mm1$theta["B", "A"]))



    if(isTRUE(all.equal(
        mm1$transitionRateMatrix["A", "A, B"], 
        exp(mm1$theta["B", "B"]) * exp(mm1$theta["A", "B"])))) {
        cat("\n \n Here they are identical too 0_1")
        print(exp(mm1$theta["A", "B"]))
        print(exp(mm1$theta["B", "A"]))
        print(exp(mm1$theta["B", "B"]))
        print(sampledGenotypes(dB))
    }

    if(isTRUE(all.equal(
        mm1$transitionRateMatrix["A", "A, B"], 
        exp(mm1$theta["A", "A"]) * exp(mm1$theta["A", "B"])))) {
        cat("\n \n Here they are identical too 0_2")
        print(exp(mm1$theta["A", "B"]))
        print(exp(mm1$theta["B", "A"]))
        print(exp(mm1$theta["A", "A"]))
        print(sampledGenotypes(dB))
    }


    if(isTRUE(all.equal(
        mm1$transitionRateMatrix["A", "A, B"], 
        exp(mm1$theta["A", "A"]) * exp(mm1$theta["B", "A"])))) {
        cat("\n \n Here they are identical too 0_3")
        print(exp(mm1$theta["A", "B"]))
        print(exp(mm1$theta["B", "A"]))
        print(exp(mm1$theta["A", "A"]))
        print(sampledGenotypes(dB))
    }

    
    ## And output matches the observed
    ## From A we move mostly to A, B
    ## From B, we can move to either A, B or B, C
})



local({
    ## 
    N <- 200
    na <- N
    nc <- N + 3 + round( 10 * runif(1))
    nab <- N + 5 + round( 10 * runif(1))
    ncd <- N + 7 + round( 10 * runif(1))
    n00 <- N/10 + round( 10 * runif(1))
    dB <- matrix(
        c(
            rep(c(1, 0, 0, 0), na) 
          , rep(c(0, 0, 1, 0), nc)
          , rep(c(1, 1, 0, 0), nab)
          , rep(c(0, 0, 1, 1), ncd)
          , rep(c(0, 0, 0, 0), n00)
        ), ncol = 4, byrow = TRUE
    )
    colnames(dB) <- LETTERS[1:4]
    sampledGenotypes(dB)
    mm4 <- do_MHN(dB)

    round(mm4$transitionMatrixCompExp, 3)
    
    ## From A to A, B:
    mm4$transitionMatrixCompExp["A", "A, B"]
    mm4$transitionRateMatrix["A", "A, B"]
    exp(mm4$theta["B", "B"]) * exp(mm4$theta["B", "A"])

    ## same, of course
    stopifnot(all.equal( mm4$transitionRateMatrix["A", "A, B"],
                        exp(mm4$theta["B", "B"]) * exp(mm4$theta["B", "A"])))

    
    ## not the same most of the time. Occasionally identical, but not a problem
    ## See code below and aobve
    ## stopifnot(mm4$transitionRateMatrix["A", "A, B"] != 
    ##           exp(mm4$theta["B", "B"]) * exp(mm4$theta["A", "B"]))

    ## stopifnot(mm4$transitionRateMatrix["A", "A, B"] != 
    ##           exp(mm4$theta["A", "A"]) * exp(mm4$theta["A", "B"]))

    ## stopifnot(mm4$transitionRateMatrix["A", "A, B"] != 
    ##           exp(mm4$theta["A", "A"]) * exp(mm4$theta["B", "A"]))



    
    ## same, of course
    stopifnot(all.equal( mm4$transitionRateMatrix["B", "A, B"],
                        exp(mm4$theta["A", "A"]) * exp(mm4$theta["A", "B"])))

    
    stopifnot(all.equal( mm4$transitionRateMatrix["C", "C, D"],
                        exp(mm4$theta["D", "D"]) * exp(mm4$theta["D", "C"])))


    stopifnot(all.equal( mm4$transitionRateMatrix["D", "C, D"],
                        exp(mm4$theta["C", "C"]) * exp(mm4$theta["C", "D"])))


    stopifnot(all.equal( mm4$transitionRateMatrix["D", "B, D"],
                        exp(mm4$theta["B", "B"]) * exp(mm4$theta["B", "D"])))


    stopifnot(all.equal( mm4$transitionRateMatrix["C, D", "A, C, D"],
                        exp(mm4$theta["A", "A"]) *
                        exp(mm4$theta["A", "D"]) *
                        exp(mm4$theta["A", "C"])))


    stopifnot(all.equal( mm4$transitionRateMatrix["B, C", "A, B, C"],
                        exp(mm4$theta["A", "A"]) *
                        exp(mm4$theta["A", "B"]) *
                        exp(mm4$theta["A", "C"])))

    stopifnot(all.equal( mm4$transitionRateMatrix["B, C", "B, C, D"],
                        exp(mm4$theta["D", "D"]) *
                        exp(mm4$theta["D", "B"]) *
                        exp(mm4$theta["D", "C"])))
    

    stopifnot(all.equal( mm4$transitionRateMatrix["B, C, D", "A, B, C, D"],
                        exp(mm4$theta["A", "A"]) *
                        exp(mm4$theta["A", "B"]) *
                        exp(mm4$theta["A", "D"]) *
                        exp(mm4$theta["A", "C"])))


    round(mm4$transitionRateMatrix, 3)
    
    ## And output matches the observed
    ## From A we move mostly to A, B
    ## From B, we can move to either A, B or B, C
    ## From C, we move to C, D
   
})




local({
    ## 

    N <- 200
    na <- N
    nc <- 5 * N + 13 + round( 10 * runif(1))
    nd <- round(1.5 * N) + round( 10 * runif(1))
    nab <- 2 * N + round( 10 * runif(1))
    nabc <- 3 * N + round(10 * runif(1))
    nabd <- 5 + round(10 * runif(1))    
    n00 <- N/10 + round( 10 * runif(1))
    dB <- matrix(
        c(
            rep(c(1, 0, 0, 0), na) 
          , rep(c(0, 0, 1, 0), nc)
          , rep(c(0, 0, 0, 1), nd)            
          , rep(c(1, 1, 0, 0), nab)
          , rep(c(1, 1, 1, 0), nabc)
          , rep(c(1, 1, 0, 1), nabd)            
          , rep(c(0, 0, 0, 0), n00)
        ), ncol = 4, byrow = TRUE
    )
    colnames(dB) <- LETTERS[1:4]
    sampledGenotypes(dB)
    mm77 <- do_MHN(dB)
    options(width = 300)
    
    round(mm77$transitionMatrixCompExp, 3)
    round(mm77$transitionMatrixTimeDiscretized, 3)
          
    ## From A to A, B:
    mm77$transitionMatrixCompExp["A", "A, B"]
    mm77$transitionRateMatrix["A", "A, B"]
    exp(mm77$theta["B", "B"]) * exp(mm77$theta["B", "A"])

    ## same, of course
    stopifnot(all.equal( mm77$transitionRateMatrix["A", "A, B"],
                        exp(mm77$theta["B", "B"]) * exp(mm77$theta["B", "A"])))

    ## not the same most of the time, though occasionally identical. Not a problem
    if(isTRUE(all.equal(
        mm77$transitionRateMatrix["A", "A, B"], 
        exp(mm77$theta["B", "B"]) * exp(mm77$theta["A", "B"])))) {
        cat("\n \n Here they are identical too")
        print(exp(mm77$theta["A", "B"]))
        print(exp(mm77$theta["B", "A"]))
        print(exp(mm77$theta["B", "B"]))
        print(sampledGenotypes(dB))
    }


    ## Same as above: this could, be chance, be equal
    ## stopifnot(!isTRUE(all.equal(
    ##     mm77$transitionRateMatrix["A", "A, B"], 
    ##     exp(mm77$theta["A", "A"]) * exp(mm77$theta["A", "B"]))))
    
    ## stopifnot(!isTRUE(all.equal(
    ##     mm77$transitionRateMatrix["A", "A, B"], 
    ##     exp(mm77$theta["A", "A"]) * exp(mm77$theta["B", "A"]))))


    ## same, of course
    stopifnot(all.equal( mm77$transitionRateMatrix["A", "A, C"],
                        exp(mm77$theta["C", "C"]) * exp(mm77$theta["C", "A"])))
    ## not necessarily always different if theta[a,c] != theta[c, a]
    ## stopifnot(!isTRUE(all.equal(
    ##     mm77$transitionRateMatrix["A", "A, C"], 
    ##     exp(mm77$theta["C", "C"]) * exp(mm77$theta["A", "C"]))))
    
    if(isTRUE(all.equal(
        mm77$transitionRateMatrix["A", "A, C"], 
        exp(mm77$theta["C", "C"]) * exp(mm77$theta["A", "C"])))) {
        cat("\n \n Here they are identical too 2")
        print(exp(mm77$theta["A", "C"]))
        print(exp(mm77$theta["C", "A"]))
        print(exp(mm77$theta["C", "C"]))
        print(sampledGenotypes(dB))
    }

    
    ## stopifnot(!isTRUE(all.equal(
    ##     mm77$transitionRateMatrix["A", "A, C"], 
    ##     exp(mm77$theta["A", "A"]) * exp(mm77$theta["A", "C"]))))
    
    ## stopifnot(!isTRUE(all.equal(
    ##     mm77$transitionRateMatrix["A", "A, C"], 
    ##     exp(mm77$theta["A", "A"]) * exp(mm77$theta["C", "A"]))))


    ## same, of course
    stopifnot(all.equal( mm77$transitionRateMatrix["A, B", "A, B, C"],
                        exp(mm77$theta["C", "C"]) *
                        exp(mm77$theta["C", "A"]) *
                        exp(mm77$theta["C", "B"]) 
                        ))

    stopifnot(all.equal( mm77$transitionRateMatrix["A, B", "A, B, D"],
                        exp(mm77$theta["D", "D"]) *
                        exp(mm77$theta["D", "A"]) *
                        exp(mm77$theta["D", "B"]) 
                        ))

    ## A very extreme one, from a non-existing genotype
    stopifnot(all.equal( mm77$transitionRateMatrix["A, C", "A, B, C"],
                        exp(mm77$theta["B", "B"]) *
                        exp(mm77$theta["B", "A"]) *
                        exp(mm77$theta["B", "C"]) 
                        ))

    
    stopifnot(all.equal( mm77$transitionRateMatrix["B, C, D", "A, B, C, D"],
                        exp(mm77$theta["A", "A"]) *
                        exp(mm77$theta["A", "B"]) *
                        exp(mm77$theta["A", "D"]) *
                        exp(mm77$theta["A", "C"])))

    stopifnot(all.equal( mm77$transitionRateMatrix["A, B, C", "A, B, C, D"],
                        exp(mm77$theta["D", "D"]) *
                        exp(mm77$theta["D", "B"]) *
                        exp(mm77$theta["D", "A"]) *
                        exp(mm77$theta["D", "C"])))
    
    stopifnot(all.equal( mm77$transitionRateMatrix["A, B, D", "A, B, C, D"],
                        exp(mm77$theta["C", "C"]) *
                        exp(mm77$theta["C", "B"]) *
                        exp(mm77$theta["C", "D"]) *
                        exp(mm77$theta["C", "A"])))
    
    stopifnot(all.equal( mm77$transitionRateMatrix["A, C, D", "A, B, C, D"],
                        exp(mm77$theta["B", "B"]) *
                        exp(mm77$theta["B", "A"]) *
                        exp(mm77$theta["B", "D"]) *
                        exp(mm77$theta["B", "C"])))


    

    round(mm77$transitionRateMatrix, 3)
    
    ## And output matches the observed
    ## In particular: A, B to A,B,C vs A, B, D
})


## trypwire
## stopifnot(2 == 3)



local({

    ## And a six gene example. Also testing the sparse matrix implementation
    N <- 200
    na <- N
    nc <- 3 * N + 2 + round( 10 * runif(1))
    nd <- round(1.5 * N) + round( 10 * runif(1))
    nab <- 2 * N + round( 10 * runif(1))
    nabc <- 3 * N + round(10 * runif(1))
    nabd <- 5 * N + round(10 * runif(1))
    nabce <- 5 * N + round(10 * runif(1))
    nabcf <- 3 * N + round(10 * runif(1))
    nabcde <- 3 * N + round(10 * runif(1))
    nabcef <- 1 * N + round(10 * runif(1))
    nabcdef <- 2 * N + round(10 * runif(1))
    n00 <- N/10 + round( 10 * runif(1))
    dB <- matrix(
        c(
            rep(c(1, 0, 0, 0, 0, 0), na) 
          , rep(c(0, 0, 1, 0, 0, 0), nc)
          , rep(c(0, 0, 0, 1, 0, 0), nd)            
          , rep(c(1, 1, 0, 0, 0, 0), nab)
          , rep(c(1, 1, 1, 0, 0, 0), nabc)
          , rep(c(1, 1, 0, 1, 0, 0), nabd)            
          , rep(c(0, 0, 0, 0, 0, 0), n00)
          , rep(c(1, 1, 1, 0, 1, 0), nabce)
          , rep(c(1, 1, 1, 0, 0, 1), nabcf)
          , rep(c(1, 1, 1, 1, 1, 0), nabcde)
          , rep(c(1, 1, 1, 0, 1, 1), nabcef)            
          , rep(c(1, 1, 1, 1, 1, 1), nabcdef)            
        ), ncol = 6, byrow = TRUE
    )
    colnames(dB) <- LETTERS[1:6]
    sampledGenotypes(dB)
    mm2 <- do_MHN(dB)
    mm22 <- do_MHN2(dB)
    options(width = 350)
    
    round(mm2$transitionMatrixCompExp, 2)
    round(mm2$transitionMatrixTimeDiscretized, 2)


    ## Yes, sparse ones
    stopifnot(inherits(mm22$transitionRateMatrix, "dgCMatrix"))
    stopifnot(inherits(mm22$transitionMatrixTimeDiscretized, "dgCMatrix"))
    stopifnot(inherits(mm22$transitionMatrixCompExp, "dgCMatrix"))

    ## To compare same values, need to convert into matrices
    stopifnot(identical(mm2$transitionRateMatrix,
                        as.matrix(mm22$transitionRateMatrix)))
    ## Identical would fail
    stopifnot(all.equal(mm2$transitionMatrixCompExp,
                        as.matrix(mm22$transitionMatrixCompExp)))
    
    stopifnot(all.equal(mm2$transitionMatrixTimeDiscretized,
                        as.matrix(mm22$transitionMatrixTimeDiscretized)))


    ## This would fail, of course unless we turned them into matrices
    stopifnot(!identical(mm2$transitionRateMatrix,
                        mm22$transitionRateMatrix))
    stopifnot(!identical(mm2$transitionMatrixTimeDiscretized,
                        mm22$transitionMatrixTimeDiscretized))
    stopifnot(!identical(mm2$transitionMatrixCompExp,
                        mm22$transitionMatrixCompExp))

    
    ## Does not seem great for how we go
    ## ABCE to ABCDE and ABCEF. But probably a hard model anyway,
    
    ## From A to A, B:
    mm2$transitionMatrixCompExp["A", "A, B"]
    mm2$transitionRateMatrix["A", "A, B"]
    exp(mm2$theta["B", "B"]) * exp(mm2$theta["B", "A"])

    ## same, of course
    stopifnot(all.equal( mm2$transitionRateMatrix["A", "A, B"],
                        exp(mm2$theta["B", "B"]) * exp(mm2$theta["B", "A"])))

    ## not the same, but not always
    ## stopifnot(!isTRUE(all.equal(
    ##     mm2$transitionRateMatrix["A", "A, B"], 
    ##     exp(mm2$theta["B", "B"]) * exp(mm2$theta["A", "B"]))))

    ## stopifnot(!isTRUE(all.equal(
    ##     mm2$transitionRateMatrix["A", "A, B"], 
    ##     exp(mm2$theta["A", "A"]) * exp(mm2$theta["A", "B"]))))
    
    ## stopifnot(!isTRUE(all.equal(
    ##     mm2$transitionRateMatrix["A", "A, B"], 
    ##     exp(mm2$theta["A", "A"]) * exp(mm2$theta["B", "A"]))))

    if(isTRUE(all.equal(
        mm2$transitionRateMatrix["A", "A, B"], 
        exp(mm2$theta["B", "B"]) * exp(mm2$theta["A", "B"])))) {
        cat("\n \n Here they are identical too 3_1")
        print(exp(mm2$theta["A", "B"]))
        print(exp(mm2$theta["B", "A"]))
        print(exp(mm2$theta["B", "B"]))
        print(sampledGenotypes(dB))
    }
    if(isTRUE(all.equal(
        mm2$transitionRateMatrix["A", "A, B"], 
        exp(mm2$theta["A", "A"]) * exp(mm2$theta["A", "B"])))) {
        cat("\n \n Here they are identical too 3_2")
        print(exp(mm2$theta["A", "B"]))
        print(exp(mm2$theta["B", "A"]))
        print(exp(mm2$theta["A", "A"]))
        print(sampledGenotypes(dB))
    }
    if(isTRUE(all.equal(
        mm2$transitionRateMatrix["A", "A, B"], 
        exp(mm2$theta["A", "A"]) * exp(mm2$theta["B", "A"])))) {
        cat("\n \n Here they are identical too 3_3")
        print(exp(mm2$theta["A", "B"]))
        print(exp(mm2$theta["B", "A"]))
        print(exp(mm2$theta["A", "A"]))
        print(sampledGenotypes(dB))
    }



    

    ## ## same, of course
    ## stopifnot(all.equal( mm2$transitionRateMatrix["A", "A, C"],
    ##                     exp(mm2$theta["C", "C"]) * exp(mm2$theta["C", "A"])))
    ## not the same, but not always
    if(isTRUE(all.equal(
        mm2$transitionRateMatrix["A", "A, C"], 
        exp(mm2$theta["C", "C"]) * exp(mm2$theta["A", "C"])))) {
        cat("\n \n Here they are identical too 3_4")
        print(exp(mm2$theta["A", "C"]))
        print(exp(mm2$theta["C", "A"]))
        print(exp(mm2$theta["C", "C"]))
        print(sampledGenotypes(dB))
    }

    if(isTRUE(all.equal(
        mm2$transitionRateMatrix["A", "A, C"], 
        exp(mm2$theta["A", "A"]) * exp(mm2$theta["A", "C"])))) {
        cat("\n \n Here they are identical too 3_5")
        print(exp(mm2$theta["A", "C"]))
        print(exp(mm2$theta["C", "A"]))
        print(exp(mm2$theta["A", "A"]))
        print(sampledGenotypes(dB))
    }

    if(isTRUE(all.equal(
        mm2$transitionRateMatrix["A", "A, C"], 
        exp(mm2$theta["A", "A"]) * exp(mm2$theta["C", "A"])))) {
        cat("\n \n Here they are identical too 3_6")
        print(exp(mm2$theta["A", "C"]))
        print(exp(mm2$theta["C", "A"]))
        print(exp(mm2$theta["A", "A"]))
        print(sampledGenotypes(dB))
    }

        
    ## stopifnot(!isTRUE(all.equal(
    ##     mm2$transitionRateMatrix["A", "A, C"], 
    ##     exp(mm2$theta["C", "C"]) * exp(mm2$theta["A", "C"]))))

    ## Same story as above.
    ## stopifnot(!isTRUE(all.equal(
    ##     mm2$transitionRateMatrix["A", "A, C"], 
    ##     exp(mm2$theta["A", "A"]) * exp(mm2$theta["A", "C"]))))
    
    ## stopifnot(!isTRUE(all.equal(
    ##     mm2$transitionRateMatrix["A", "A, C"], 
    ##     exp(mm2$theta["A", "A"]) * exp(mm2$theta["C", "A"]))))


    ## same, of course
    stopifnot(all.equal( mm2$transitionRateMatrix["A, B", "A, B, C"],
                        exp(mm2$theta["C", "C"]) *
                        exp(mm2$theta["C", "A"]) *
                        exp(mm2$theta["C", "B"]) 
                        ))

    stopifnot(all.equal( mm2$transitionRateMatrix["A, B", "A, B, D"],
                        exp(mm2$theta["D", "D"]) *
                        exp(mm2$theta["D", "A"]) *
                        exp(mm2$theta["D", "B"]) 
                        ))

    ## A very extreme one, from a non-existing genotype
    stopifnot(all.equal( mm2$transitionRateMatrix["A, C", "A, B, C"],
                        exp(mm2$theta["B", "B"]) *
                        exp(mm2$theta["B", "A"]) *
                        exp(mm2$theta["B", "C"]) 
                        ))

    
    stopifnot(all.equal( mm2$transitionRateMatrix["B, C, D", "A, B, C, D"],
                        exp(mm2$theta["A", "A"]) *
                        exp(mm2$theta["A", "B"]) *
                        exp(mm2$theta["A", "D"]) *
                        exp(mm2$theta["A", "C"])))

    stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C", "A, B, C, D"],
                        exp(mm2$theta["D", "D"]) *
                        exp(mm2$theta["D", "B"]) *
                        exp(mm2$theta["D", "A"]) *
                        exp(mm2$theta["D", "C"])))
    
    stopifnot(all.equal( mm2$transitionRateMatrix["A, B, D", "A, B, C, D"],
                        exp(mm2$theta["C", "C"]) *
                        exp(mm2$theta["C", "B"]) *
                        exp(mm2$theta["C", "D"]) *
                        exp(mm2$theta["C", "A"])))
    
    stopifnot(all.equal( mm2$transitionRateMatrix["A, C, D", "A, B, C, D"],
                        exp(mm2$theta["B", "B"]) *
                        exp(mm2$theta["B", "A"]) *
                        exp(mm2$theta["B", "D"]) *
                        exp(mm2$theta["B", "C"])))


    ## 5 genes
    stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C, E", "A, B, C, D, E"],
                        exp(mm2$theta["D", "D"]) *
                        exp(mm2$theta["D", "A"]) *
                        exp(mm2$theta["D", "B"]) *
                        exp(mm2$theta["D", "C"]) *                        
                        exp(mm2$theta["D", "E"])))

    stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C, E", "A, B, C, E, F"],
                        exp(mm2$theta["F", "F"]) *
                        exp(mm2$theta["F", "A"]) *
                        exp(mm2$theta["F", "B"]) *
                        exp(mm2$theta["F", "C"]) *                        
                        exp(mm2$theta["F", "E"])))

    stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C, F", "A, B, C, E, F"],
                        exp(mm2$theta["E", "E"]) *
                        exp(mm2$theta["E", "A"]) *
                        exp(mm2$theta["E", "B"]) *
                        exp(mm2$theta["E", "C"]) *                        
                        exp(mm2$theta["E", "F"])))

    stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C, F", "A, B, C, D, F"],
                        exp(mm2$theta["D", "D"]) *
                        exp(mm2$theta["D", "A"]) *
                        exp(mm2$theta["D", "B"]) *
                        exp(mm2$theta["D", "C"]) *                        
                        exp(mm2$theta["D", "F"])))

    ## six genes
    stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C, D, E", "A, B, C, D, E, F"],
                        exp(mm2$theta["F", "F"]) *
                        exp(mm2$theta["F", "A"]) *
                        exp(mm2$theta["F", "B"]) *
                        exp(mm2$theta["F", "C"]) *
                        exp(mm2$theta["F", "D"]) *                        
                        exp(mm2$theta["F", "E"])))

    stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C, D, F", "A, B, C, D, E, F"],
                        exp(mm2$theta["E", "E"]) *
                        exp(mm2$theta["E", "A"]) *
                        exp(mm2$theta["E", "B"]) *
                        exp(mm2$theta["E", "C"]) *
                        exp(mm2$theta["E", "D"]) *                        
                        exp(mm2$theta["E", "F"])))


    stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C, E, F", "A, B, C, D, E, F"],
                        exp(mm2$theta["D", "D"]) *
                        exp(mm2$theta["D", "A"]) *
                        exp(mm2$theta["D", "B"]) *
                        exp(mm2$theta["D", "C"]) *
                        exp(mm2$theta["D", "E"]) *                        
                        exp(mm2$theta["D", "F"])))

    ## six genes, with sparse. Totally unnecessary since they were identical to non-sparse
    stopifnot(all.equal( mm22$transitionRateMatrix["A, B, C, D, E", "A, B, C, D, E, F"],
                        exp(mm22$theta["F", "F"]) *
                        exp(mm22$theta["F", "A"]) *
                        exp(mm22$theta["F", "B"]) *
                        exp(mm22$theta["F", "C"]) *
                        exp(mm22$theta["F", "D"]) *                        
                        exp(mm22$theta["F", "E"])))

    stopifnot(all.equal( mm22$transitionRateMatrix["A, B, C, D, F", "A, B, C, D, E, F"],
                        exp(mm22$theta["E", "E"]) *
                        exp(mm22$theta["E", "A"]) *
                        exp(mm22$theta["E", "B"]) *
                        exp(mm22$theta["E", "C"]) *
                        exp(mm22$theta["E", "D"]) *                        
                        exp(mm22$theta["E", "F"])))


    stopifnot(all.equal( mm22$transitionRateMatrix["A, B, C, E, F", "A, B, C, D, E, F"],
                        exp(mm22$theta["D", "D"]) *
                        exp(mm22$theta["D", "A"]) *
                        exp(mm22$theta["D", "B"]) *
                        exp(mm22$theta["D", "C"]) *
                        exp(mm22$theta["D", "E"]) *                        
                        exp(mm22$theta["D", "F"])))

    

    round(mm2$transitionRateMatrix, 3)
    
    ## And output matches the observed
    ## In particular: A, B to A,B,C vs A, B, D
})





## ## Seems hard to make it any faster
## pd3 <- profileExpr(theta_to_trans_rate_3(Theta.BC, inner_transition = inner_transitionRate_3_1))
## ## pd32 <- profileExpr(theta_to_trans_rate_3(Theta.BC, inner_transition = inner_transitionRate_3_2))
## options(width = 150)
## funSummary(pd3)
## callSummary(pd3)

## ## funSummary(pd32)
## ## callSummary(pd32)


library(codetools)
checkUsageEnv(env = .GlobalEnv)



