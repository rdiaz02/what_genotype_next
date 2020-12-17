## Examples of local maxima

source("schill-trans-mat.R")
library(OncoSimulR)



options(width = 200)

options(digits = 2)

## 
N <- 200
na <- N
nc <- N + round( 10 * runif(1))
nab <- N + round( 10 * runif(1))
nbc <- N + round( 10 * runif(1))
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
do_MHN(dB) ## this crashes. Of course: don't do silly things like passing
           ## a data frame with one or more colSums == 0
dB2 <- dB[, -4]
do_MHN(dB2) 


## 
N <- 200
na <- N
nc <- N + round( 10 * runif(1))
nab <- N + round( 10 * runif(1))
ncd <- N + round( 10 * runif(1))
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
do_MHN(dB) 




## Playing with different thresholds to call something a local max
play <- function(mm) {
    ## mm <- do.call(rbind, lapply(1:(reps/2),
    ##                             function(x) do.call(rbind, list(m1, m1))))

    ## Can we detect local peaks here?
    out1 <- do_MHN(mm)
    ## ## The next could never catch local peaks
    ## out2 <- do_MHN(mm, "competingExponentials")

    cat("\n theta \n")
    print(round(out1$theta, 2))
    
    cat("\n transition matrix \n")
    print(round(out1$transitionMatrixTimeDiscretized, 2))
    ## round(out2$transitionMatrixTimeDiscretized, 2)
    cat("\n diagonal \n")
    print(round(diag(out1$transitionMatrixTimeDiscretized), 2))
    cat("\n*******************\n\n")
    cat("\n diag >= 0.9 \n")
    ## Way too many?
    print(which(diag(out1$transitionMatrixTimeDiscretized) >= 0.9))
    cat("\n*******************\n\n")
    cat("\n diag >= 0.95 \n")
    ## Way too many?
    print(which(diag(out1$transitionMatrixTimeDiscretized) >= 0.95))
    cat("\n*******************\n\n")
    cat("\n diag >= 0.975 \n")
    ## Way too many?
    print(which(diag(out1$transitionMatrixTimeDiscretized) >= 0.975))
    cat("\n*******************\n\n")    
    cat("\n diag >= 0.99 \n")
    ## Hummm... we miss A
    print(which(diag(out1$transitionMatrixTimeDiscretized) >= 0.99))
    cat("\n*******************\n\n")
    cat("\n diag >= 0.999 \n")
    ## Too far
    print(which(diag(out1$transitionMatrixTimeDiscretized) >= 0.999))
    
}



## From Diaz-Uriarte & Vasallo, Fig. 1, B
## hard example, because 1,1,0,0 is local max
## but there is also a 1,1,1,0. 
N <- 200
na <- N
nb <- N
nc <- N
nab <- N
nac <- N
nbc <- N
nabc <- N
nbcd <- N
## nabd <- 50
nabcd <- N
n00 <- N/10
dB <- matrix(
    c(
        rep(c(1, 0, 0, 0), na) 
      , rep(c(0, 1, 0, 0), nb)
      , rep(c(0, 0, 1, 0), nc)
      , rep(c(1, 1, 0, 0), nab)
      , rep(c(1, 0, 1, 0), nac)
      , rep(c(0, 1, 1, 0), nbc)
      , rep(c(1, 1, 1, 0), nabc)
        ##      , rep(c(1, 1, 0, 1), nabd)
      , rep(c(0, 1, 1, 1), nbcd)
      , rep(c(1, 1, 1, 1), nabcd)
      , rep(c(0, 0, 0, 0), n00)
    ), ncol = 4, byrow = TRUE
)
colnames(dB) <- LETTERS[1:4]
play(dB)


## From Diaz-Uriarte & Vasallo, Fig. 1, B
## hard example, because 1,1,0,0 is local max
## but there is also a 1,1,1,0. 
N <- 200
na <- N
nb <- N
nc <- N
nab <- N * 2
nac <- N * 2
nbc <- N
nabc <- N
nbcd <- N
## nabd <- 50
nabcd <- N
n00 <- N/10
dB <- matrix(
    c(
        rep(c(1, 0, 0, 0), na) 
      , rep(c(0, 1, 0, 0), nb)
      , rep(c(0, 0, 1, 0), nc)
      , rep(c(1, 1, 0, 0), nab)
      , rep(c(1, 0, 1, 0), nac)
      , rep(c(0, 1, 1, 0), nbc)
      , rep(c(1, 1, 1, 0), nabc)
        ##      , rep(c(1, 1, 0, 1), nabd)
      , rep(c(0, 1, 1, 1), nbcd)
      , rep(c(1, 1, 1, 1), nabcd)
      , rep(c(0, 0, 0, 0), n00)
    ), ncol = 4, byrow = TRUE
)
colnames(dB) <- LETTERS[1:4]
play(dB)





## From Diaz-Uriarte & Vasallo, Fig. 1, C
## misses A as local max, unless is proportionately large
N <- 50
na <- N * 10
nb <- N
nc <- N
## nab <- 0
## nac <- 0
nbc <- N
nabc <- N * 2
nbcd <- N * 2
## nabd <- 50
## nabcd <- 0
n00 <- N/10
dB <- matrix(
    c(
        rep(c(1, 0, 0, 0), na) 
      , rep(c(0, 1, 0, 0), nb)
      , rep(c(0, 0, 1, 0), nc)
        ## , rep(c(1, 1, 0, 0), nab)
        ## , rep(c(1, 0, 1, 0), nac)
      , rep(c(0, 1, 1, 0), nbc)
      , rep(c(1, 1, 1, 0), nabc)
        ##      , rep(c(1, 1, 0, 1), nabd)
      , rep(c(0, 1, 1, 1), nbcd)
        ##      , rep(c(1, 1, 1, 1), nabcd)
      , rep(c(0, 0, 0, 0), n00)
    ), ncol = 4, byrow = TRUE
)
colnames(dB) <- LETTERS[1:4]
play(dB)




## Much simpler, because a two-way interaction model (first order epistasis)
## can model this.
## Wt -> A
## Wt -> B -> BC
## WT -> D -> DE -> DEF
N <- 200
na <- N
nb <- N
nbc <- N
nd <- N
nde <- N
ndef <- N
n00 <- N
m1 <- matrix(
    c(
        rep(c(1, 0, 0, 0, 0, 0), na)
      , rep(c(0, 1, 0, 0, 0, 0), nb)
      , rep(c(0, 1, 1, 0, 0, 0), nbc)
      , rep(c(0, 0, 0, 1, 0, 0), nd)
      , rep(c(0, 0, 0, 1, 1, 0), nde)
      , rep(c(0, 0, 0, 1, 1, 1), ndef)
      , rep(c(0, 0, 0, 0, 0, 0), n00)  
    ), ncol = 6, byrow = TRUE
)
colnames(m1) <- LETTERS[1:6]
sampledGenotypes(m1)
play(m1)



##   Misses A and BC unless large frequency
## Wt -> A
## Wt -> B -> BC           [cannot tell if can reach BC from C]
## Wt -> B -> BD -> BDE
## Wt -> C -> CE -> CEF
## Wt -> D -> DE -> DEF
## "    "     "  -> BDE   [so convergence with 3rd path]
N <- 200
na <- N 
nb <- N
nbc <- N 
nbd <- N
nbde <- N
nc <- N
nce <- N
ncef <- N 
nd <- N
nde <- N
ndef <- N
n00 <- N
m1 <- matrix(
    c(
        rep(c(1, 0, 0, 0, 0, 0), na)
      , rep(c(0, 1, 0, 0, 0, 0), nb)
      , rep(c(0, 1, 1, 0, 0, 0), nbc)
      , rep(c(0, 1, 0, 1, 0, 0), nbd)
      , rep(c(0, 1, 0, 1, 1, 0), nbde)       
      , rep(c(0, 0, 1, 0, 0, 0), nc)
      , rep(c(0, 0, 1, 0, 1, 0), nce)
      , rep(c(0, 0, 1, 0, 1, 1), ncef)              
      , rep(c(0, 0, 0, 1, 0, 0), nd)
      , rep(c(0, 0, 0, 1, 1, 0), nde)
      , rep(c(0, 0, 0, 1, 1, 1), ndef)
      , rep(c(0, 0, 0, 0, 0, 0), n00)  
    ), ncol = 6, byrow = TRUE
)
colnames(m1) <- LETTERS[1:6]
play(m1)



## Ehh??? either makes all single local max or too few local max.
## Why?
## and the transition matrix is weird, very weird
## Wt -> A
## Wt -> B -> BC           [cannot tell if can reach BC from C]
## Wt -> C -> BC
## Wt -> D -> DE -> ADE -> ADEF
## Wt -> F -> EF -> AEF -> ADEF

N <- 200
na <- N * 3
nb <- N
nc <- N
nbc <- N * 2
nd <- N
nde <- N * 2
nade <- N * 2.5
nadef <- N * 1.5
nf <- N
nef <- N
naef <- N * 1.5
n00 <- 0
m1 <- matrix(
    c(
        rep(c(1, 0, 0, 0, 0, 0), na)
      , rep(c(0, 1, 0, 0, 0, 0), nb)
      , rep(c(0, 0, 1, 0, 0, 0), nc)
      , rep(c(0, 1, 1, 0, 0, 0), nbc)
      , rep(c(0, 0, 0, 1, 0, 0), nd)
      , rep(c(0, 0, 0, 1, 1, 0), nde)
      , rep(c(1, 0, 0, 1, 1, 0), nade)
      , rep(c(1, 0, 0, 1, 1, 1), nadef)               
      , rep(c(0, 0, 0, 0, 0, 1), nf)
      , rep(c(0, 0, 0, 0, 1, 1), nef)
      , rep(c(1, 0, 0, 0, 1, 1), naef)       
      , rep(c(0, 0, 0, 0, 0, 0), n00)  
    ), ncol = 6, byrow = TRUE
)
colnames(m1) <- LETTERS[1:6]
sampledGenotypes(m1) ## check
play(m1)



## Ehh??!!
## Like former, but without convergence in ADEF
## Wt -> A
## Wt -> B -> BC           [cannot tell if can reach BC from C]
## Wt -> C -> BC
## Wt -> D -> DE -> ADE 
## Wt -> F -> EF -> AEF 
N <- 200
na <- N
nb <- N
nc <- N
nbc <- N
nd <- N
nde <- N
nade <- N
nf <- N
nef <- N
naef <- N
n00 <- N
m1 <- matrix(
    c(
        rep(c(1, 0, 0, 0, 0, 0), na)
      , rep(c(0, 1, 0, 0, 0, 0), nb)
      , rep(c(0, 0, 1, 0, 0, 0), nc)
      , rep(c(0, 1, 1, 0, 0, 0), nbc)
      , rep(c(0, 0, 0, 1, 0, 0), nd)
      , rep(c(0, 0, 0, 1, 1, 0), nde)
      , rep(c(1, 0, 0, 1, 1, 0), nade)
      , rep(c(0, 0, 0, 0, 0, 1), nf)
      , rep(c(0, 0, 0, 0, 1, 1), nef)
      , rep(c(1, 0, 0, 0, 1, 1), naef)       
      , rep(c(0, 0, 0, 0, 0, 0), n00)  
    ), ncol = 6, byrow = TRUE
)
colnames(m1) <- LETTERS[1:6]
play(m1)



## Here, the time discretized gives 0 for WT-WT
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
mm1


## Here, the time discretized gives 0 for WT-WT and for A-A
N <- 200
na <- 5
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
mm1

## Here, WT-WT is not 0 in time discretized
N <- 200
na <- 5
nc <- N + 21 ## round( 10 * runif(1))
nab <- N + 5 ## round( 10 * runif(1)) + 
nbc <- N + 7 + round( 10 * runif(1))
n00 <- 50
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
mm1

