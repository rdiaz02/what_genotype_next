### LIBRARIES & SETTINGS
###

library(compiler)
enableJIT(3)
setCompilerOptions(optimize = 3)
options(stringsAsFactors = FALSE)

library(OncoSimulR)
library(parallel)
library(pryr)
library(dplyr)
library(testthat)
library(stringr)
library(combinat)
library(pbapply)
library(pbmcapply)









### FUNCTIONS
###

# prompt: where to get files from/save files to?
askDir <- function(defaultDir=".",message="Choose a directory.") {
  
  # ask for directory
  cat(message)
  cat("\n")
  cat(paste("Leave blank for default: ",defaultDir,sep=""))
  cat("\n")
  cat("  > Enter path: ")
  input <- readLines("stdin",1)
  cat("\n")
  
  # set default directory if input was blank
  if (input == "") input <- defaultDir
  
  # change dot (.) by current directory
  if (grepl(".",input)) {
    input <- sub(".",getwd(),input,fixed=T)
  }
  
  # show a warning if directory doesn't exist
  if (!dir.exists(input)) {
    cat("WARNING: Provided directory doesn't exist")
    cat("\n")
    cat("Continue? [y/n] ")
    check <- readLines("stdin",1)
    cat("\n")
    if(check!="y") stop()
  }
  
  return(input)
  
}

# get number of genes from a file name
numGenesFromFile <- function(filename) {
  case1 <- as.numeric(str_match(filename, "ng_(.*?)_fl")[2])
  case2 <- as.numeric(str_match(filename, "ngenes_(.*?)__")[2])
  out <- c(case1,case2)
  out <- out[!is.na(out)]
  return(out)
}

# get fitness landscape ID from file names
flIDFromFile <- function(filename) {
  str_match(filename, "ID_(.*?)_")[2]
}

# get size_split from file names
sizeSplitFromFile <- function(filename) {
  paste("size_split_",
        str_match(filename, "size_split_(.*?)_")[2],
        sep="")
}

# get detection regime from file names
detectFromFile <- function(filename) {
  paste("detect_",
        str_match(filename, "_detect_(.*?).RData")[2],
        sep="")
}

# get landscape type from file names
flTypeFromFile <- function(filename) {
  str_match(filename, "_typeLandscape_(.*?)_")[2]
}

# find number of mutations of a genotype given its name
nMut <- function(genotype) {
  
  f <- function(genotype) {
    out <- length(strsplit(genotype,",")[[1]])
    if (genotype=="none" | genotype=="end" | genotype=="any") out <- NA
    if (tolower(genotype)=="wt" | tolower(genotype)=="root") out <- 0
    return(out)
  }
  
  out <- sapply(genotype,f)
  names(out) <- NULL
  
  return(out)
  
}

# generate genotype names, IDs, etc
makeGenotypes <- function(numGenes) {
  
  # undo repeated number of genes
  numGenes <- unique(numGenes)
  
  # function: convert genotype IDs to names
  genotypeIDtoName <- function(ID) {
    paste(LETTERS[which(strsplit(ID,"")[[1]]==1)],collapse=", ")
  }
  
  # function: get gentoype names and IDs from given number of genes
  allGenotypes <- function(numGenes) {
    
    # full list of genotypes
    numGenotypes <- 2^numGenes
    
    genotypeIDs <- do.call(expand.grid,rep(list(0:1),numGenes))
    genotypeIDs <- apply(genotypeIDs,1,paste,collapse="")
    genotypeNames <- sapply(genotypeIDs,genotypeIDtoName)
    genotypeMuts <- nMut(genotypeNames)
    names(genotypeMuts) <- genotypeNames
    names(genotypeIDs) <- genotypeNames
    
    out <- list(numGenes=numGenes,numGenotypes=numGenotypes,
                genotypeIDs=genotypeIDs,genotypeNames=genotypeNames,
                genotypeMuts=genotypeMuts)
    
    return(out)
    
  }
  
  aux <- vector(mode="list",length=length(numGenes))
  names(aux) <- numGenes
  
  # get all required genotype names and IDs for all possible numGenes in the files
  for (i in 1:length(numGenes)) {
    aux[[as.character(numGenes[i])]] <- allGenotypes(numGenes[i])
  }
  
  return(aux)
  
}

# find chain of (n+1)-genotypes given a POM (i.e. chain of n-genotypes)
nextInLOD <- function(lod,pom) {

  # find (n+1)-genotype in a LOD given a n-genotype
  nextInLOD_n <- function(genotype,lod) {
    nextGenotype <- lod[nMut(lod)==nMut(genotype)+1]
    
    # return "end" or "none" if there is no next genotype
    if (length(nextGenotype)==0) {
      if(genotype==lod[length(lod)]) {
        nextGenotype <- "end"
      } else {
        nextGenotype <- "none"
      }
    }
    
    return(nextGenotype)
  }
  
  # apply to full LOD
  out <- sapply(pom,nextInLOD_n,lod=lod)
  names(out) <- NULL
  return(out)
  
}

# output of nextInLOD() in transition matrix form
nextInLOD_transitionMatrix <- function(lod,pom,allGenotypes) {
  numGenotypes <- length(allGenotypes)
  t <- matrix(0,nrow=numGenotypes,ncol=numGenotypes+2)
  rownames(t) <- allGenotypes
  colnames(t) <- c(allGenotypes,"none","end")
  rows <- pom
  cols <- nextInLOD(lod,pom)
  for (j in 1:length(rows)) {
    t[rows[j]==rownames(t),cols[j]==colnames(t)] <-
      t[rows[j]==rownames(t),cols[j]==colnames(t)] + 1
  }
  
  # row normalization
  tNorm <- matrix(rep(rowSums(t),numGenotypes+2),
                  nrow=numGenotypes,ncol=numGenotypes+2)
  tNorm[tNorm==0] <- 1
  t <- t/tNorm
  
  # return t or update transitionMatrix
  return(t)
  
}

# check which genotypes out of a given array are in a POM
whichInPOM <- function(genotypes,pom) {
  
  # check if a certain genotype is in a given POM
  isInPOM <- function(genotype,pom) {
    as.numeric(genotype %in% pom)
  }
  
  out <- sapply(genotypes,isInPOM,pom=pom)
  names(out) <- NULL
  return(out)
  
}

# unfuse columns/rows of the kind X_Y_Z_... (X, Y, Z, ... = gene names)
unfuseGenotype <- function(genotype) {
  
  # remove everything after first "_"
  f <- function(genotype) {
    paste(gsub("_.*","",
               strsplit(genotype,", ")[[1]]),
          collapse=", ")
  }
  unfusedGenotype <- sapply(genotype,f)
  names(unfusedGenotype) <- NULL
  
  # flag if something was done
  n <- unfusedGenotype == genotype
  if (sum(n) != length(n)) {
    flag <- paste("unfused genotype(s):",paste(genotype[!n],collapse=" ; "))
  } else {
    flag <- NULL
  }
  
  out <- list(unfusedGenotype=unfusedGenotype,
              flags=flag)
  
  return(out)
  
}

# adjust the order of the letters that appear in a genotype name
orderGenotype <- function(genotype) {
  
  f <- function(genotype) {
    paste(sort(strsplit(genotype,", ")[[1]]),collapse=", ")
  }
  orderedGenotype <- sapply(genotype,f)
  names(orderedGenotype) <- NULL
  
  # return flag if something was done
  n <- genotype != orderedGenotype
  if (sum(n)>0) {
    flag <- paste("ordered genotype(s):",paste(genotype[n],collapse=" ; "))
  } else {
    flag <- NULL
  }
  
  out <- list(orderedGenotype=orderedGenotype,
              flags=flag)
  
  return(out)
  
}

# take a matrix from CPM output, format it into a transition matrix
structTransitionMatrix <- function(transitionMatrix_in,allGenotypes) {
  
  numGenotypes <- length(allGenotypes)
  
  # proceed only if there is a matrix
  if (is.matrix(transitionMatrix_in)) {
    
    # initialize output (formatted matrix)
    out <- transitionMatrix_in
    
    # rename "WT" genotype to "" in matrix
    colnames(out)[colnames(out)=="WT"] <- ""
    rownames(out)[rownames(out)=="WT"] <- ""
    
    # unfuse genotype names in rows/columns
    unfG <- unfuseGenotype(colnames(out))
    flags <- unfG$flags
    colnames(out) <- unfG$unfusedGenotype
    rownames(out) <- unfG$unfusedGenotype
    
    # adjust order of letters in genotype names
    ordG <- orderGenotype(colnames(out))
    flags <- c(flags,ordG$flags)
    colnames(out) <- ordG$orderedGenotype
    rownames(out) <- ordG$orderedGenotype
    
    # flags to return
    flags <- paste(flags,collapse=" | ")
    
    # add "end" column
    genotypeEnd <- which(rowSums(out)==0)
    out <- cbind(out,end=0)
    out[genotypeEnd,"end"] <- 1
    
    # diagonal elements interpreted as "end" entries (remaining in the same state)
    out[,"end"] <- out[,"end"] + diag(out)
    diag(out) <- 0
    
    # how many genotypes are missing? (not accessible according to CPM)
    accGenotypes <- dim(out)[1]
    missGenotypes <- numGenotypes - accGenotypes
    
    # if >0 gentoypes are missing, add those to the matrix (as "zero" entries)
    if (missGenotypes>0) {
      out <- cbind(out,matrix(0,ncol=missGenotypes,nrow=nrow(out)))
      out <- rbind(out,matrix(0,nrow=missGenotypes,ncol=ncol(out)))
      
      missGenotypeNames <- allGenotypes[!(allGenotypes %in% rownames(out))]
      rownames(out)[(accGenotypes+1):nrow(out)] <- missGenotypeNames
      colnames(out)[(accGenotypes+2):ncol(out)] <- missGenotypeNames
    }
    
    # rearrange rows and columns
    newColOrder <- c(match(allGenotypes,colnames(out)),accGenotypes+1)
    newRowOrder <- match(allGenotypes,rownames(out))
    out <- out[newRowOrder,newColOrder]
    
    # add "none" column
    out <- cbind(out[,-ncol(out)],none=0,end=out[,ncol(out)])
    
    # row-normalize (probably not needed at this point?)
    # this IS needed since some of the un-structured CPM output matrices are
    # not row-normalized to begin with
    if(T) {
      tNorm <- matrix(rep(rowSums(out),ncol(out)),
                      nrow=nrow(out),ncol=ncol(out))
      tNorm[tNorm==0] <- 1
      out <- out/tNorm
    }
    
  } else {
    
    # flag if input is not a matrix
    flags <- paste("not a matrix:",transitionMatrix_in)
    
    # return NULL if input is not a matrix
    out <- NULL
    
  }
  
  # return restructured matrix and flags
  out <- list(transitionMatrix=out,flags=flags)
  return(out)
  
}

# make a null matrix from an input matrix with row and column names
nullMatrix <- function(A) {
  
  # not-normalized null matrix
  muts <- nMut(rownames(A))
  N <- matrix(0,ncol=ncol(A),nrow=nrow(A))
  colnames(N) <- colnames(A)
  rownames(N) <- rownames(A)
  for (i in 0:(max(muts)-1)){
    N[muts==i,c(muts,Inf,Inf)==(i+1)] <- 1
  }
  N[,"end"] <- 1
  N[!rownames(N)=="","none"] <- 1
  
  # normalized null matrix
  Nn <- do.call(cbind,replicate(ncol(N),
                                rowSums(N),simplify=F))
  Nn <- N/Nn
  
  return(list(raw=N,norm=Nn))
  
}

# check what genotypes are accessible given a transition matrix
isAccessible <- function(A) {
  
  # if A is not a matrix, return NULL
  if(!is.matrix(A)) {
    return(NULL)
  } else {
    out <- colSums(A[,1:nrow(A)])>0
    out[names(out) %in% c("","WT","root")] <- TRUE
    return(out)
  }
  
}

# average square of differences (square-rooted)
sqDiff <- function(p, q, fix=F) {
  
  pq <- c(p, q)
  if(any(is.na(pq))) return(NA)
  if(any(pq > 1)) stop("some prob. > 1")
  if(any(pq < 0)) stop("some prob. < 0")
  stopifnot(isTRUE(all.equal(sum(p), 1.0)) | isTRUE(all.equal(sum(p), 0.0)))
  stopifnot(isTRUE(all.equal(sum(q), 1.0)) | isTRUE(all.equal(sum(q), 0.0)))
  stopifnot(identical(length(p), length(q)))
  
  if (all(pq==0)) return(0) # minimum possible value
  
  if ((any(p>0) & all(q==0)) | (all(p==0) & any(q>0))) {
    if(fix) {
      if(any(p>0)) q <- rep(1/length(q),length(q))
      if(all(p==0)) p <- rep(1/length(p),length(p))
      return(sqDiff(p,q))
    } else {
      return(1) # max possible value
    }
  }
  
  if(any(p>0) & any(q>0)) return(sqrt(mean((p-q)^2)))
  
}

## See below for improved code: my_kldiv_2
# Kullback-Leibler distance:  Dkl = D(x || y)
my_kldiv <- function(x, y) {
  
  ## Based on https://stackoverflow.com/a/27953616
  ## with modification for removal of 0 in numerator of Kullback-Leibler as per
  ## https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence#Definition
    ## Using log2 so that max is 1. See wikipedia or
    ##  Lin, J. (1991). "Divergence measures based on the shannon entropy"
    ## IEEE Transactions on Information Theory. 37 (1): 145–151. 
    
  ## remove zeros in x
  rmx <- which(x == 0)
  if(length(rmx))
    return(sum(x[-rmx] * log2(x[-rmx]/y[-rmx])))
  else
    return(sum(x * log2(x/y)))
  
}
## See below for improved code. jensen_shannon_2

# Jensen-Shannon entropy
jensen_shannon <- function(p, q, fix=F) {
  
  ## FIXME: could have checked none in m is 0?
  ## well, one could pass a silly c(a, b, 0), c(e, f, 0)
  ## and should we not allow for that?
  ## KL will not crash (we remove those, as we should, and they are really
  ## irrelevant anyway)
  
  ## FIXME: generalized function so it can handle all-zeros vectors as input
  
  pq <- c(p, q)
  if(any(is.na(pq))) return(NA)
  if(any(pq > 1)) stop("some prob. > 1")
  if(any(pq < 0)) stop("some prob. < 0")
  stopifnot(isTRUE(all.equal(sum(p), 1.0)) | isTRUE(all.equal(sum(p), 0.0)))
  stopifnot(isTRUE(all.equal(sum(q), 1.0)) | isTRUE(all.equal(sum(q), 0.0)))
  stopifnot(identical(length(p), length(q)))
  m <- 0.5 * (p + q)
  
  if (all(pq==0)) return(0) # minimum possible value
  
  if ((any(p>0) & all(q==0)) | (all(p==0) & any(q>0))) {
    # return(1) #max. value
    
    ### FIXME: the output of this (p is all 0 and q is finite or vice-versa)
    ### now depends on the "eq" flag: if F it will return 1, if T it will
    ### equiprobabilize the all 0 input and calculate jensen_shannon as usual
    
    if(fix) {
      if(any(p>0)) q <- rep(1/length(q),length(q))
      if(all(p==0)) p <- rep(1/length(p),length(p))
      return(jensen_shannon(p,q))
    } else {
      return(1) # max. possible value
    }
  }
  
  if(any(p>0) & any(q>0)) return(0.5 *  (my_kldiv(p, m) + my_kldiv(q, m)))

    ## FIXME: yes, FIXME for the future.
    ## how do we know we never fail to fall into those cases?
    ## stop("JS undefined")
    
  ## rmp <- which(p == 0)
  ## rmq <- which(q == 0)
  ## kl1 <- sum(p[-rmp] * log(p[-rmp]/m[-rmp]))
  ## kl2 <- sum(q[-rmq] * log(q[-rmq]/m[-rmq]))
  ## return(0.5 * (kl1 + kl2))
  ## <- 0.5 * (sum(p * log(p / m)) + sum(q * log(q / m)))
  
}

# simplified Jensen-Shannon
easyJS <- function(p,q,fixNulls=F) {
  
  # sanity checks
  pq <- c(p, q)
  if(any(is.na(pq))) return(NA)
  if(any(pq > 1)) stop("some prob. > 1")
  if(any(pq < 0)) stop("some prob. < 0")
  stopifnot(isTRUE(all.equal(sum(p), 1.0)) | isTRUE(all.equal(sum(p), 0.0)))
  stopifnot(isTRUE(all.equal(sum(q), 1.0)) | isTRUE(all.equal(sum(q), 0.0)))
  stopifnot(identical(length(p), length(q)))
  
  # if fixNulls, use null model instead of all-zero vectors
  if(fixNulls) {
    if(all(p==0)) p <- rep(1/length(p),length(p))
    if(all(q==0)) q <- rep(1/length(q),length(q))
  }
  
  # "standard" case: both p and q are non-zero
  m <- 0.5 * (p + q)
  if(any(p>0) & any(q>0))
    return(0.5 *  (my_kldiv(p, m) + my_kldiv(q, m)))
  else
    return(NA)
  
}

# Hellinger distance
hellinger_distance <- function(p, q, fix=F) {
  
  pq <- c(p, q)
  if(any(is.na(pq))) return(NA)
  if(any(pq > 1)) stop("some prob. > 1")
  if(any(pq < 0)) stop("some prob. < 0")
  stopifnot(isTRUE(all.equal(sum(p), 1.0)) | isTRUE(all.equal(sum(p), 0.0)))
  stopifnot(isTRUE(all.equal(sum(q), 1.0)) | isTRUE(all.equal(sum(q), 0.0)))
  stopifnot(identical(length(p), length(q)))
  
  if (all(pq==0)) return(0) # minimum possible value
  
  if ((any(p>0) & all(q==0)) | (all(p==0) & any(q>0))) {
    if(fix) {
      if(any(p>0)) q <- rep(1/length(q),length(q))
      if(all(p==0)) p <- rep(1/length(p),length(p))
      return(hellinger_distance(p,q))
    } else {
      return(1) # max possible value
    }
  }
  
  if(any(p>0) & any(q>0)) return( (1/sqrt(2)) * sqrt(sum((sqrt(p) - sqrt(q))^2)) )
  
}

# Spearman correlation
spearman_corr <- function(p, q) {
  
  pq <- c(p, q)
  if(any(is.na(pq))) return(c(rho=NA,
                              pval=NA))
  if(any(pq > 1)) stop("some prob. > 1")
  if(any(pq < 0)) stop("some prob. < 0")
  stopifnot(isTRUE(all.equal(sum(p), 1.0)) | isTRUE(all.equal(sum(p), 0.0)))
  stopifnot(isTRUE(all.equal(sum(q), 1.0)) | isTRUE(all.equal(sum(q), 0.0)))
  stopifnot(identical(length(p), length(q)))
  
  if (all(pq==0)) return(c(rho=1, # maximum possible value
                           pval=1))
  
  if ((any(p>0) & all(q==0)) | (all(p==0) & any(q>0))) {
    return(c(rho=0, # min possible value
             pval=1))
  }
  
  if(any(p>0) & any(q>0)) {
  
    ct <- suppressWarnings(cor.test(p, q, method = "spearman"))
    
    # if output is NA (p and/or q equiprobable), return 0 instead
    if (is.na(ct$estimate[["rho"]])) {
      ct$estimate[["rho"]] <- 0
      ct$p.value <- 1
    }
    
    return(c(rho=ct$estimate[["rho"]],
             pval=ct$p.value))
    
  }
  
}

eq <- function(p,threshold=0) {
  
  if(any(is.na(p))) return(NA)
  if(any(p > 1)) stop("some prob. > 1")
  if(any(p < 0)) stop("some prob. < 0")
  stopifnot(isTRUE(all.equal(sum(p), 1.0)) | isTRUE(all.equal(sum(p), 0.0)))
  
  p[p<=threshold] <- 0
  p[p>threshold] <- 1
  if (any(p>0)) p <- p/sum(p)
  
  return(p)
}

# wrapper: all statistics
allStats <- function(x, y, threshold=c(0,0)) {
    c(sqDiff = sqDiff(x, y),
      sqDiff_fix = sqDiff(x, y, fix=T),
      sqDiff_eq = sqDiff(eq(x,threshold=threshold[1]),eq(y,threshold=threshold[2])),
      sqDiff_eq_fix = sqDiff(eq(x,threshold=threshold[1]),eq(y,threshold=threshold[2]),fix=T),
      JS = sqrt(jensen_shannon(x, y)),
      JS_fix = sqrt(jensen_shannon(x, y, fix=T)),
      JS_eq = sqrt(jensen_shannon(eq(x,threshold=threshold[1]),eq(y,threshold=threshold[2]))),
      JS_eq_fix = sqrt(jensen_shannon(eq(x,threshold=threshold[1]),eq(y,threshold=threshold[2]),fix=T)),
      He = hellinger_distance(x, y),
      He_fix = hellinger_distance(x, y, fix=T),
      He_eq = hellinger_distance(eq(x,threshold=threshold[1]),eq(y,threshold=threshold[2])),
      He_eq_fix = hellinger_distance(eq(x,threshold=threshold[1]),eq(y,threshold=threshold[2]),fix=T),
      rank_corr = spearman_corr(x, y)[["rho"]],
      rank_corr_p = spearman_corr(x, y)[["pval"]])
}

# compare matrices
compareMatrices <- function(A, B, rowWeights=NULL, threshold=c(0,0)) {
  
  # output NAs in certain situations (see conditions below)
  statsNA <- allStats(NA,NA)
  
  # return NAs if at least one input is not a matrix
  if(!is.matrix(A) | !is.matrix(B)) {
    
    if (is.matrix(A)) X <- A
    if (is.matrix(B)) X <- B
    
    out <- as.data.frame(matrix(NA,nrow=nrow(X)+1,ncol=length(statsNA)))
    colnames(out) <- names(statsNA)
    rownames(out) <- c(rownames(X),"any")
    
  } else {
    
    # check matrix dimensions and row/column names match
    if (ncol(A) != ncol(B) | nrow(A) != nrow(B)) {
      stop("matrix dimensions must agree")
    }
    if (any(colnames(A) != colnames(B)) | any(rownames(A) != rownames(B))) {
      stop("matrix row and column names must agree")
    }
    
    # if there was no input for rowWeights, make them uniform
    if (!length(rowWeights)) rowWeights <- rep(1,nrow(A))
    if (length(rowWeights) != nrow(A)) {
      stop("length of weights vector and number of matrix rows differ")
    }
    
    # null matrix (for non-trivial element access)
    N <- nullMatrix(A)$raw
    
    # row-wise stats
    rowStats <- lapply(1:nrow(A),
                       function(row) {
                         p <- A[row,N[row,]==1]
                         q <- B[row,N[row,]==1]
                         return(allStats(p,q,
                                         threshold=threshold))
                      })
    rowStats <- matrix(unlist(rowStats),nrow=length(rowStats),byrow=T)
    colnames(rowStats) <- names(statsNA)
    rownames(rowStats) <- rownames(A)
    
    # identify trivial rows (all zeros in A and B, or only one non-trivial element)
    trivialRows <- (!rowSums(A) & !rowSums(B)) | rowSums(N)==1
    
    # adjust weights
    rowWeights[trivialRows] <- 0
    rowWeights <- rowWeights/sum(rowWeights)
    rowWeights <- do.call(cbind,replicate(ncol(rowStats),rowWeights,simplify=F))
    
    # get full matrix stats (average)
    matrixStats <- colSums(rowStats[!trivialRows,]*rowWeights[!trivialRows,])
    
    # final output
    out <- as.data.frame(rbind(rowStats,any=matrixStats))
    
  }
  
  names(out) <- c("sqDiff","sqDiff_fix","sqDiff_eq","sqDiff_eq_fix",
                  "js","js_fix","js_eq","js_eq_fix",
                  "hellinger","hellinger_fix","hellinger_eq","hellinger_eq_fix",
                  "spearman","spearman_pval")
  return(out)
  
}

# simplified compareMatrices
easyCompare <- function(A,B,numGenes,fixNulls=F) {
  
  # return NAs if at least one input is not a matrix
  if(!is.matrix(A) | !is.matrix(B)) {
    
    g <- makeGenotypes(numGenes)[[1]]$genotypeNames
    out <- rep("err_no_matrix",length(g))
    names(out) <- g
    
  } else {
    
    # null matrix (for non-trivial element access)
    N <- nullMatrix(A)$raw
    
    # row-wise JS
    out <- sapply(1:nrow(A),
                  function(row) {
                    p <- A[row,N[row,]==1]
                    q <- B[row,N[row,]==1]
                    return(easyJS(p,q,fixNulls))
                  })
    names(out) <- rownames(A)
    
    }
  
  return(out)
  
}

## x: oncoSimulIndiv output -> table of most abundant genotype at each
##       sampling time.
##     correct_length: if TRUE, all individuals contribute the same to
##     the estimates  (i.e., individuals with longer runs do not contribute
##     more)
## called from trueFreqs, below.
maxg_freqs2 <- function(x, correct_length) {
    ## Yes, loop much faster here
    tmp <- x$pops.by.time[, -1]
    r <- nrow(tmp)
    out <- rep(-9, r)
    for(i in seq(r)) out[i] <- which.max(tmp[i, ])
    tt <- table(x$GenotypesLabels[out])
    if(correct_length) tt <- tt/sum(tt)
    data.frame(v1 = names(tt),
               v2 = c(tt),
               row.names = NULL,
               fix.empty.names = FALSE,
               check.names = FALSE,
               stringsAsFactors = FALSE)
}

## x: a simulation output object, simobj -> frequencies of most abundant
##      genotypes
##     if correct_length = TRUE, all individuals contribute the same
##     to the estimates (i.e., individuals with longer runs do not contribute
##     more)
true_Freqs <- function(x, correct_length = TRUE) {
    cat("\n      starting true_Freqs")
    lex <- length(x)
    stopifnot(lex <= 20000)
    if(lex > 20000) warning("\n More than 20000 \n")
    ## FIXME: if we get failures or warnings:
    ##   pass filename (i.e., file)
    ##   and return it under each condition
    ##   simpler and faster than tryCatch
    t22 <- Sys.time()
    ## all_maxg_freqs <- lapply(x,  maxg_freqs2, correct_length)
    all_maxg_freqs <- mclapply(x,  maxg_freqs2,
                                correct_length, mc.cores = detectCores())
    cat("\n          done maxg_freqs2; took ", Sys.time() - t22, "\n")

    ## all_maxg_freqs <- mlp(x)
    all_maxg_freqs_C <- dplyr::bind_rows(all_maxg_freqs)
    ## dplgb <- dplyr::group_by(all_maxg_freqs_C, v1)
    ## true_freqs <- as.data.frame(dplyr::summarize(dplgb, Freq = sum(v2)))
    ## colnames(true_freqs)[1] <- "Genotype"
    true_freqs <- aggregate( v2 ~ v1, data = all_maxg_freqs_C, FUN = sum)
    colnames(true_freqs) <- c("Genotype", "TrueFreq")
    true_freqs$TrueProp <- true_freqs$TrueFreq/sum(true_freqs$TrueFreq)
    true_freqs
}

## true frequencies, but you pass a file and also get
##   all the info: ID.
true_Freqs_full <- function(x, correct_length = TRUE) {
    if(exists("simo")) stop("A simo here!")
    ## try(rm(simo), silent = TRUE)
    cat("\n Doing simul file ", x, "...")
    cat(" ... loading file;")
    t1 <- Sys.time()
    simo <- loadRData(x)
    t2 <- Sys.time()
    cat("loaded; took", t2 - t1, "\n")
    ID <- strsplit(strsplit(x, "simobj_ID_")[[1]][2],
                   "_ng_")[[1]][1]
    ## stopifnot(nchar(ID) == 16)
    gg <- true_Freqs(simo)
    rm(simo)
    whichwt <- (gg$Genotype == "")
    stopifnot(sum(whichwt) == 1)
    gg$Genotype[whichwt] <- "WT"
    return(
        cbind(ID = ID,
              gg,
              stringsAsFactors = FALSE, row.names = NULL)
    )
}

## Possibly faster versions of the above?
## load inside the inner function, so do not pass big things
## and parallelize inside
## x: a simulation output object, simobj -> frequencies of most abundant
##      genotypes
##     if correct_length = TRUE, all individuals contribute the same
##     to the estimates (i.e., individuals with longer runs do not contribute
##     more)
true_Freqs_2 <- function(file, correct_length = TRUE) {
    if(exists("simo")) stop("A simo here!")
    ## try(rm(simo), silent = TRUE)
    cat(" ... loading file;")
    t1 <- Sys.time()
    simo <- loadRData(file)
    t2 <- Sys.time()
    cat("loaded; took", t2 - t1, "\n")
    lex <- length(simo)
    stopifnot(lex <= 20000)
    if(lex > 20000) warning("\n More than 20000 \n")
    ## FIXME: if we get failures or warnings:
    ##   pass filename (i.e., file)
    ##   and return it under each condition
    ##   simpler and faster than tryCatch
    ## mclapply, called from another mclapply, blows up ram 
    all_maxg_freqs <- lapply(simo,  maxg_freqs2, correct_length)
    ## all_maxg_freqs <- mclapply(simo,  maxg_freqs2,
    ##                            correct_length, mc.cores = detectCores())
    cat("\n          done maxg_freqs2; took ", Sys.time() - t2, "\n")
    rm(simo)
    ## all_maxg_freqs <- mlp(x)
    all_maxg_freqs_C <- dplyr::bind_rows(all_maxg_freqs)
    ## dplgb <- dplyr::group_by(all_maxg_freqs_C, v1)
    ## true_freqs <- as.data.frame(dplyr::summarize(dplgb, Freq = sum(v2)))
    ## colnames(true_freqs)[1] <- "Genotype"
    true_freqs <- aggregate( v2 ~ v1, data = all_maxg_freqs_C, FUN = sum)
    colnames(true_freqs) <- c("Genotype", "TrueFreq")
    true_freqs$TrueProp <- true_freqs$TrueFreq/sum(true_freqs$TrueFreq)
    true_freqs
}

## true frequencies, but you pass a file, and also get
##   all the info: ID.
true_Freqs_full_2 <- function(x, correct_length = TRUE) {
    cat("\n Doing simul file ", x, "...")
    ID <- strsplit(strsplit(x, "simobj_ID_")[[1]][2],
                   "_ng_")[[1]][1]
    ## stopifnot(nchar(ID) == 16)
    gg <- true_Freqs_2(file = x)
    cat(" ; done true_Freqs_2\n")
    whichwt <- (gg$Genotype == "")
    stopifnot(sum(whichwt) == 1)
    gg$Genotype[whichwt] <- "WT"
    return(
        cbind(ID = ID,
              gg,
              stringsAsFactors = FALSE, row.names = NULL)
    )
}

## file of all_paths* and vector of ANALYZED_ files -> file from ANALYZED_
##      with matching sampling scheme, but always sample size of 4000
get_4000_file <- function(file, filesANALYZED) {
    ID <- strsplit(strsplit(file, "_ID_")[[1]][2], "__rnst_")[[1]][1]
    idf <- grep(paste0("ANALYZED__ID_", ID,
                       "_rnst_"),
                filesANALYZED, value = TRUE)
    idf4000 <- grep("_size_split_4000\\.rds", idf, value = TRUE)

    if(grepl("__detect_large\\.rds", file)) {
        idf4000_s <- grep("_beta_A_5_beta_B_3_size_split_", idf4000, value = TRUE)
    } else if(grepl("__detect_small\\.rds", file)) {
        idf4000_s <- grep("_beta_A_3_beta_B_5_size_split_", idf4000, value = TRUE)
    } else if(grepl("__detect_unif\\.rds", file)) {
        idf4000_s <- grep("_beta_A_1_beta_B_1_size_split_", idf4000, value = TRUE)
    }
    stopifnot(length(idf4000_s) == 1)
    return(idf4000_s)
}

## ANALYZED file -> genotype frequencies from 20000 samples
freqs_from_sampling <- function(file) {
    x <- readRDS(file)
    aid <- do.call(rbind, lapply(x$allM, function(x) x$input_data))
    ## Get genotypes
    sampledGenot <- OncoSimulR::sampledGenotypes(aid)
    colnames(sampledGenot) <- c("Genotype", "SampledFreq")
    sampledGenot$SampledProp <-
        sampledGenot$SampledFreq/sum(sampledGenot$SampledFreq)
    sampledGenot
}

## ANALYZED file -> genotype frequencies from 20000 samples and full info ID,
##    sampling, etc
freqs_from_sampling_full <- function(file) {
    ## cat("\n Doing sampling file ", file, "\n")
    x <- readRDS(file)
    stopifnot(x$size_split == 4000)
    aid <- do.call(rbind, lapply(x$allM, function(x) x$input_data))
    stopifnot(nrow(aid) == 20000)
    ## Get genotypes
    sampledGenot <- OncoSimulR::sampledGenotypes(aid)
    ## rm the sampledGenotypes class
    class(sampledGenot) <- "data.frame"
    colnames(sampledGenot) <- c("Genotype", "SampledFreq")
    sampledGenot$SampledProp <-
        sampledGenot$SampledFreq/sum(sampledGenot$SampledFreq)
    ## sampledGenot
    ## Only ID and detect are needed for matching. The rest are for
    ## paranoid checks. Since files are huge, removed.
    return(
        cbind(ID = x$ID,
              detect = x$detect,
              ## rnst = x$rnst,
              ## ngenes = x$ngenes,
              ## initSize = x$initSize,
              sampledGenot,
              stringsAsFactors = FALSE, row.names = NULL))
}



## Wrapper to OncoSimulR::sampledGenotypes
## with a check
obs_freq_dataset <- function(datat, size_split, replicate_num) {
    stopifnot(nrow(datat) == size_split)
    sG <- OncoSimulR::sampledGenotypes(datat)
    ## rm the sampledGenotypes class
    class(sG) <- "data.frame"
    colnames(sG) <- c("Genotype", "ObservedFreq")
    sumF <- sum(sG$ObservedFreq)
    stopifnot(sumF == size_split)
    sG$ObservedProp <-
        sG$ObservedFreq/sumF
    return(cbind(replicate = replicate_num,
                 sG,
                 stringsAsFactors = FALSE, row.names = NULL))
}


## ANALYZED file -> genotype frequencies from samples and full info ID,
##    sampling, etc
obs_freqs_in_sample <- function(file) {
    ## cat("\n Doing sampling file ", file, "\n")
    x <- readRDS(file)
    size_split <- x$size_split

    ## Don't : waiting, and too many processes'
    ## of <- mclapply(1:length(x$allM),
    ##              function(i) fsg(x$allM[[i]]$input_data,
    ##                              size_split,
    ##                              i),
    ##              mc.cores = detectCores()
    ##              )
    stopifnot(length(x$allM) == 5)
    of <- lapply(1:length(x$allM),
                 function(i) obs_freq_dataset(x$allM[[i]]$input_data,
                                 size_split,
                                 i)
                 )

    ofu <- dplyr::bind_rows(of)
    ## Only ID and detect are needed for matching. The rest are for
    ## paranoid checks. Since files are huge, removed.
    return(
        cbind(ID = x$ID,
              detect = x$detect,
              size_split = size_split,
              ## rnst = x$rnst,
              ## ngenes = x$ngenes,
              ## initSize = x$initSize,
              ofu,
              stringsAsFactors = FALSE, row.names = NULL))
}


## fitness landscape file -> rank of genotypes (two columns, second
##             with non-viable genots set to NA)
fitness_rank_genotypes <- function(file) {
    x <- readRDS(file)
    ID <- x$landscape
    ngenes <- ncol(x$fitness_landscape) - 1
    stopifnot(ngenes %in% c(7, 10) )
    stopifnot(x$fitness_landscape[, "Fitness"] > 0)
    geneNames <- LETTERS[1:ngenes]
    genots <- apply(x$fitness_landscape[, -(ngenes + 1)],
                    1,
                    function(u)
                       paste(sort(geneNames[as.logical(u)]), collapse = ", ") 
                    )
    genots[genots == ""] <- "WT"
    whichwt <- which(genots == "WT")
    stopifnot(x$fitness_landscape[whichwt, "Fitness"] == 1)
    ## largest fitness gets 1st rank
    fitnessRank <- rank( -(x$fitness_landscape[, "Fitness"]) )
    fitnessRank2 <- fitnessRank
    fitnessRank2[x$fitness_landscape[, "Fitness"] <= 1e-9] <- NA
    ## marank <- max(fitnessRank)
    mirank <- min(fitnessRank)
    wtr <- fitnessRank[whichwt]
    stopifnot(wtr > mirank) 
    ## stopifnot(wtr < marank) ## need not be true in some representable
    stopifnot(!is.na(fitnessRank2[wtr]))
    data.frame(ID = ID,
               Genotype = genots,
               fitnessRank = fitnessRank,
               fitnessRankNoZero = fitnessRank2,
               stringsAsFactors = FALSE,
               row.names = NULL
               )
}

## Modified from https://stackoverflow.com/a/25455968/3670562
loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    lsout <- ls()
    if(length(lsout) > 2) stop("More than one object in filename")
    get(lsout[lsout != "fileName"])
}



## file of a simulation output object, simobj -> frequencies each genotype is local max.
local_max <- function(file) {
    if(exists("simo")) stop("A simo here!")
    ## try(rm(simo), silent = TRUE)
    cat(" ... loading file;")
    t1 <- Sys.time()
    simo <- loadRData(file)
    t2 <- Sys.time()
    cat("loaded; took", t2 - t1, "\n")
    ID <- strsplit(strsplit(file, "simobj_ID_")[[1]][2],
                   "_ng_")[[1]][1]
    
    lex <- length(simo)
    stopifnot(lex <= 20000)
    if(lex > 20000) warning("\n More than 20000 \n")

    ## lod <- mclapply(simo, OncoSimulR:::LOD.internal,
    ##                        mc.cores = detectCores())
    ## lode <- mclapply(lod, function(x) x[length(x)],
    ##                  mc.cores = detectCores())

    lode <- mclapply(simo, function(x) {
        tmp <- OncoSimulR:::LOD.internal(x)
        return(tmp[length(tmp)])
    },
    mc.cores = detectCores())
    rm(simo)
    cat(" ; done lod last\n")
    stopifnot(length(lode) == 20000)
    lodt <- as.data.frame(table(unlist(lode)),
                          stringsAsFactors = FALSE)
    colnames(lodt) <- c("Genotype", "FreqLocalMax")
    slm <- sum(lodt$FreqLocalMax)
    stopifnot(slm == 20000)
    lodt$PropLocalMax <- lodt$FreqLocalMax/slm

    return(
        cbind(ID = ID,
              lodt,
              stringsAsFactors = FALSE, row.names = NULL)
    )
}

## file of a simulation output object, simobj -> frequencies each genotype
## is local max.
## this is the non-parallelized version
local_max_np <- function(file) {
    if(exists("simo")) stop("A simo here!")
    ## try(rm(simo), silent = TRUE)
    cat(" ... loading file;")
    t1 <- Sys.time()
    simo <- loadRData(file)
    t2 <- Sys.time()
    cat("loaded; took", t2 - t1, "\n")
    ID <- strsplit(strsplit(file, "simobj_ID_")[[1]][2],
                   "_ng_")[[1]][1]
    
    lex <- length(simo)
    stopifnot(lex <= 20000)
    if(lex > 20000) warning("\n More than 20000 \n")

    ## lod <- mclapply(simo, OncoSimulR:::LOD.internal,
    ##                        mc.cores = detectCores())
    ## lode <- mclapply(lod, function(x) x[length(x)],
    ##                  mc.cores = detectCores())

    lode <- lapply(simo, function(x) {
        tmp <- OncoSimulR:::LOD.internal(x)
        return(tmp[length(tmp)])
    })
    rm(simo)
    cat(" ; done lod last\n")
    stopifnot(length(lode) == 20000)
    lodt <- as.data.frame(table(unlist(lode)),
                          stringsAsFactors = FALSE)
    colnames(lodt) <- c("Genotype", "FreqLocalMax")
    slm <- sum(lodt$FreqLocalMax)
    stopifnot(slm == 20000)
    lodt$PropLocalMax <- lodt$FreqLocalMax/slm

    return(
        cbind(ID = ID,
              lodt,
              stringsAsFactors = FALSE, row.names = NULL)
    )
}

# sanity check stats
sanityStats <- function(df) {
  
  # remove "any" row
  df <- df[df$sourceGenotype!="any",]
  
  # if both matrix rows are all-zero:
  n <- df$sourceGenotype_freqInPOM==0 & df$sourceGenotype_accessible==F
  #   fix = T or F should be having no effect
  expect_equal(df$sqDiff[n],df$sqDiff_fix[n])
  expect_equal(df$sqDiff_eq[n],df$sqDiff_eq_fix[n])
  expect_equal(df$js[n],df$js_fix[n])
  expect_equal(df$js_eq[n],df$js_eq_fix[n])
  expect_equal(df$hellinger[n],df$hellinger_fix[n])
  expect_equal(df$hellinger_eq[n],df$hellinger_eq_fix[n])
  #   equiprobabilizing should be having no effect
  expect_equal(df$sqDiff[n],df$sqDiff_eq[n])
  expect_equal(df$sqDiff_fix[n],df$sqDiff_eq_fix[n])
  expect_equal(df$js[n],df$js_eq[n])
  expect_equal(df$js_fix[n],df$js_eq_fix[n])
  expect_equal(df$hellinger[n],df$hellinger_eq[n])
  expect_equal(df$hellinger_fix[n],df$hellinger_eq_fix[n])
  
  # if both matrix rows are NOT all-zero:
  n <- df$sourceGenotype_freqInPOM>0 & df$sourceGenotype_accessible==T
  #   fix = T or F should be having no effect
  expect_equal(df$sqDiff[n],df$sqDiff_fix[n])
  expect_equal(df$sqDiff_eq[n],df$sqDiff_eq_fix[n])
  expect_equal(df$js[n],df$js_fix[n])
  expect_equal(df$js_eq[n],df$js_eq_fix[n])
  expect_equal(df$hellinger[n],df$hellinger_fix[n])
  expect_equal(df$hellinger_eq[n],df$hellinger_eq_fix[n])
  #   equiprobabilizing should be having an effect (or not, in specific cases
  #   where the original vectors are already equiprobable). Also, sometimes
  #   equiprobabilizing will induce a decrease in performance
  #   (e.g. p<-c(0.9,0.05,0,0,0.05)  q<-c(0.9,0.05,0,0.05,0))
  #   Any alternative is possible in this case.
  
  # if matrix row is all-zero for the simulations but not for the method
  # or vice-versa (functions should be symmetrical):
  n <- (df$sourceGenotype_freqInPOM==0 & df$sourceGenotype_accessible==T) |
    (df$sourceGenotype_freqInPOM>0 & df$sourceGenotype_accessible==F)
  #   fix = T or F should be having an effect (consistently improving performance)
  stopifnot(all(df$sqDiff[n]>df$sqDiff_fix[n]))
  stopifnot(all(df$sqDiff_eq[n]>df$sqDiff_eq_fix[n]))
  stopifnot(all(df$js[n]>df$js_fix[n]))
  stopifnot(all(df$js_eq[n]>df$js_eq_fix[n]))
  stopifnot(all(df$hellinger[n]>df$hellinger_fix[n]))
  stopifnot(all(df$hellinger_eq[n]>df$hellinger_eq_fix[n]))
  #   equiprobabilizing should be having no effect with no fixing (all-zero
  #   vectors remain the same, thus output should still be maximal difference)
  expect_equal(df$sqDiff[n],df$sqDiff_eq[n])
  stopifnot(all(df$sqDiff[n]==1))
  expect_equal(df$js[n],df$js_eq[n])
  stopifnot(all(df$js[n]==1))
  expect_equal(df$hellinger[n],df$hellinger_eq[n])
  stopifnot(all(df$hellinger[n]==1))
  #   equiprobabilizing with fixing should yield consistently better performance
  #   than fixing without equiprobabilizing (or the same, if input non-zero vector
  #   was already equiprobable)
  #   This is not straightforward to see, but it is always true
  stopifnot(all(df$sqDiff_fix[n]>=df$sqDiff_eq_fix[n]))
  stopifnot(all(df$js_fix[n]>=df$js_eq_fix[n]))
  stopifnot(all(df$hellinger_fix[n]>=df$hellinger_eq_fix[n]))
  
}

# given a subset of genotypes in the LOD and a final (fixated) genotype, check if they are compatible
isCompatible <- function(genotypesLOD,genotypeFinal) {
  
  # function: convert genotype name (A, B, C...) to ID (11100...)
  genotypeNametoID <- function(genotype) {
    mutList <- match(strsplit(genotype,", ")[[1]],LETTERS)
    id <- rep(0,max(mutList))
    id[mutList] <- 1
    id <- paste(id,collapse="")
    return(id)
  }
  
  # same but for just one input (will sapply later)
  f <- function(genotypeLOD) {
    
    # check if inputs are names or IDs, convert them to IDs if necessary
    if (grepl(", ",genotypeLOD)) genotypeLOD <- genotypeNametoID(genotypeLOD)
    if (grepl(", ",genotypeFinal)) genotypeFinal <- genotypeNametoID(genotypeFinal)
    
    # make array of mutated (1) and non-mutated (0) loci
    g_in <- strsplit(genotypeLOD,"")[[1]]
    g_out <- strsplit(genotypeFinal,"")[[1]]
    
    # add non-mutated loci if necessary
    if (length(g_in)<length(g_out)) g_in <- c(g_in,rep(0,length(g_out)-length(g_in)))
    if (length(g_in)>length(g_out)) g_out <- c(g_out,rep(0,length(g_in)-length(g_out)))
    
    # check if mutated loci in g_in are also mutated in g_out
    return(all(g_out[g_in==1]==1))
    
  }
  
  return(sapply(genotypesLOD,f))
  
}

# prepare matrices for fast row-wise comparison (intended for use with biological data)
prepareBiolMatrices <- function(A,B,genots) {
  
  # A: matrix given by 1st CPM
  # B: matrix given by 2nd CPM
  # genots: observable genotypes (typically those with frequency > 0 in the sample)
  
  # sanity check: are matrices square and consistently named?
  stopifnot(ncol(A)==nrow(A))
  stopifnot(ncol(B)==nrow(B))
  stopifnot(all(colnames(A)==rownames(A)))
  stopifnot(all(colnames(B)==rownames(B)))
  
  ### FIXME: need to be careful with gene ordering in genotype names. E.g. genotype 
  ### "ATP2B2, TP53" can sometimes appear as "TP53, ATP2B2". The function orderGenotypes()
  ### fixes this (originally developed to deal with genotype naming such as "B, C_A" together
  ### with the unfuseGenotype() function)
  rownames(A) <- orderGenotype(rownames(A))$orderedGenotype
  colnames(A) <- rownames(A)
  rownames(B) <- orderGenotype(rownames(B))$orderedGenotype
  colnames(B) <- rownames(B)
  
  # add "end" column
  diagA <- diag(A)
  diag(A) <- 0
  A <- cbind(A,end=diagA)
  
  diagB <- diag(B)
  diag(B) <- 0
  B <- cbind(B,end=diagB)
  
  # trim matrices (keep observable rows only, then remove all-zero columns)
  A <- A[rownames(A) %in% genots,]
  B <- B[rownames(B) %in% genots,]
  
  A <- A[,colSums(A)>0]
  B <- B[,colSums(B)>0]
  
  # if observable genotypes are missing, add all-zero rows
  add_these <- matrix(0,ncol=ncol(A),nrow=sum(!genots %in% rownames(A)))
  rownames(add_these) <- genots[!genots %in% rownames(A)]
  A <- rbind(A,add_these)
  A <- A[genots,]
  
  add_these <- matrix(0,ncol=ncol(B),nrow=sum(!genots %in% rownames(B)))
  rownames(add_these) <- genots[!genots %in% rownames(B)]
  B <- rbind(B,add_these)
  B <- B[genots,]
  
  # if column names mismatch, add all-zero columns
  genots_union <- union(colnames(A),colnames(B))
  
  add_these <- matrix(0,nrow=nrow(A),ncol=sum(!genots_union %in% colnames(A)))
  colnames(add_these) <- genots_union[!genots_union %in% colnames(A)]
  A <- cbind(A,add_these)
  A <- A[,genots_union]
  
  add_these <- matrix(0,nrow=nrow(B),ncol=sum(!genots_union %in% colnames(B)))
  colnames(add_these) <- genots_union[!genots_union %in% colnames(B)]
  B <- cbind(B,add_these)
  B <- B[,genots_union]
  
  # sanity check
  rm(add_these)
  stopifnot(all(colnames(A)==colnames(B)) & all(rownames(A)==rownames(B)))
  
  # return formatted matrices
  return(list(A,B))
  
}

# jensen-shannon, optimized to work on sparse matrices (intended for use with biological data)
jsBiol <- function(p,q,nmut,ngenes) {
  
  # p: vector of probabilities
  # q: vector of probabilities
  # nmut: number of mutations of the source genotype
  # ngenes: total number of genes
  
  # sanity checks
  ### FIXED: these are performed in the jensen_shannon() function itself, no need to repeat
  ### here.
  ### In addition: with biological data, if instead of using
  ### stopifnot(isTRUE(all.equal(something,1.0)))
  ### we simply use
  ### stopifnot(something==1)
  ### for some reason the stop condition triggers when something is different to 1 by a small
  ### (~1e-16) amount. This is weird (both statements above should behave equally given
  ### numeric inputs). Ramon was the one to write the sanity checks in the isTRUE() form,
  ### so maybe he knows why the discrepancy?
  # stopifnot(sum(p)==1 | sum(p)==0)
  # stopifnot(sum(q)==1 | sum(q)==0)
  
  # base case: both methods give a prediction
  if (isTRUE(all.equal(sum(p),1)) & isTRUE(all.equal(sum(q),1))) {
    out <- jensen_shannon(p,q)
  }
  
  ### FIXME: if none of the methods give a prediction, return 0 (equivalent to using the same
  ### null model for both methods, therefore both give the same prediction and JS=0)
  if (isTRUE(all.equal(sum(p),0)) & isTRUE(all.equal(sum(q),0))) {
    out <- 0
  }
  
  # if one method gives a prediction and the other doesn't...
  if ((isTRUE(all.equal(sum(p),1)) & isTRUE(all.equal(sum(q),0)))
     |
     (isTRUE(all.equal(sum(p),0)) & isTRUE(all.equal(sum(q),1)))) {
    
    # we define p as the vector that has a prediction (is non-zero)
    p <- p+q
    p <- p[p>0]
    
    # number of potential genotypes to jump to from the sourceGenotype
    ngenots <- choose(ngenes,nmut+1) + 1 # + 1 because of the "end"
    
    # q (null prediction) is equiprobable
    q <- rep(1/ngenots,ngenots)
    
    # add zeros to p to match the length of q
    p <- c(p,rep(0,length(q)-length(p)))
    
    # calculate jensen-shannon distance normally
    out <- jensen_shannon(p,q)
    
  }
  
  return(out)
  
}





### TESTS
###

# numGenesFromFile
test_that("numGenesFromFile", {
  
  file_1 <- "some-random-filename_ng_4_fl_some-more-filename.some-extension"
  file_2 <- "some-random-filename_ngenes_9__$%&/.some-extension"
  
  expect_equal(numGenesFromFile(file_1),4)
  expect_equal(numGenesFromFile(file_2),9)
  
})

# nMut
test_that("nMut", {
  
  g <- c("","A","B","C","A, B","A, C","B, C","A, B, C")
  n <- c(0,1,1,1,2,2,2,3)
  
  expect_equal(nMut(g),n)
  
})

# makeGenotypes
test_that("makeGenotypes", {
  
  g <- c("","A","B","C","A, B","A, C","B, C","A, B, C")
  g_made <- makeGenotypes(3)$`3`$genotypeNames
  
  expect_setequal(g,g_made)
  
})

# structTransitionMatrix
test_that("structTransitionMatrix", {
  
  g <- c("","A","B","C","A, B","A, C","B, C","A, B, C")
  t <- matrix(c(c(0,0.25,0.75,0,0),
                c(0,0.05,0,0.95,0),
                c(0,0,0,1,0),
                c(0,0,0,0,1),
                c(0,0,0,0,1)),
              nrow=5,byrow=T)
  colnames(t) <- c("","A","C","A, C","A, B, C")
  rownames(t) <- colnames(t)
  
  t_struct <- matrix(c(c(0,.25,0,.75,0,0,0,0,0,0),
                       c(0,0,0,0,0,.95,0,0,0,.05),
                       c(0,0,0,0,0,0,0,0,0,0),
                       c(0,0,0,0,0,1,0,0,0,0),
                       c(0,0,0,0,0,0,0,0,0,0),
                       c(0,0,0,0,0,0,0,1,0,0),
                       c(0,0,0,0,0,0,0,0,0,0),
                       c(0,0,0,0,0,0,0,0,0,1)),
                     nrow=length(g),byrow=T)
  rownames(t_struct) <- g
  colnames(t_struct) <- c(g,"none","end")

  expect_equal(structTransitionMatrix(t,g)$transitionMatrix,t_struct)
  
})

# orderGenotype
test_that("orderGenotype", {
  
  p <- c("WT","A","B, A","B, D, A")
  p_ord <- c("WT","A","A, B","A, B, D")
  
  expect_equal(orderGenotype(p)$orderedGenotype,p_ord)
  
})

# unfuseGenotype
test_that("unfuseGenotype", {

  p <- c("WT","A_C","A_C, B","A_C, B, D_E_F")
  p_unf <- c("WT","A","A, B","A, B, D")
  
  expect_equal(unfuseGenotype(p)$unfusedGenotype,p_unf)

})

# eq
test_that("eq", {
  
  p1 <- c(0.1,0,0,0.7,0,0.2)
  p2 <- c(0,0,0,0,0,0)
  
  expect_equal(eq(p1),c(1,0,0,1,0,1)/3)
  expect_equal(eq(p2),p2)
  
})

# jensen_shannon
local({
  ## ## understand KLdiv output
  ## p <- abs(rnorm(15, 0, sd = 3)^2)
  ## p <- p/sum(p)
  ## q <- runif(15)
  ## q <- q/sum(q)
  ## KLdiv(cbind(p, q))
  ## my_kldiv(p, q) ## 1st row, 2nd col
  ## my_kldiv(p, q)
  
  p <- c(0.1, 0.3, 0.6, 0, 0)
  q <- c(0, 0, 0, 0.5, 0.5)
  ## stopifnot(isTRUE(all.equal(jensen_shannon(p, q), log(2))))
  ## Recall we are now using log2
  stopifnot(isTRUE(all.equal(jensen_shannon(p, q), 1)))
  
  stopifnot(isTRUE(all.equal(jensen_shannon(p, p), 0)))
  stopifnot(isTRUE(all.equal(jensen_shannon(q, q), 0)))
  
  ## example from stackoverflow
  ## https://stackoverflow.com/questions/11226627/jensen-shannon-divergence-in-r
  p <- c(0.00029421, 0.42837957, 0.1371827, 0.00029419, 0.00029419,
         0.40526004, 0.02741252, 0.00029422, 0.00029417, 0.00029418)
  
  q <- c(0.00476199, 0.004762, 0.004762, 0.00476202, 0.95714168,
         0.00476213, 0.00476212, 0.00476202, 0.00476202, 0.00476202)
  
  ## stopifnot(isTRUE(all.equal(jensen_shannon(p, q),
  ##                            0.6457538)))
  ## Again, log base 2
  stopifnot(isTRUE(all.equal(jensen_shannon(p, q),
                               0.6457538/log(2))))
    
  allStats(p, q)
  rm(p, q)
})

# jensen_shannon (with or without equiprobabilization, with or without fixing)
test_that("jensen_shannon (many scenarios)", {
  
  p <- c(0.1,0,0,0.7,0,0.2)
  q <- c(0.2,0,0,0.1,0,0.7)
  r <- c(0,0,0,0.9,0.1,0)
  z <- c(0,0,0,0,0,0)
  
  # general case: if both inputs are different and not all-zero, output should be
  # a number >0
  expect_gt(jensen_shannon(p,q),0)
  
  # if both inputs are all-zero, output should be 0
  expect_equal(jensen_shannon(z,z),0)
  
  # if one of the inputs is all-zero and the other isn't, with fix=F should return 1
  expect_equal(jensen_shannon(p,z),1)
  
  # previous case: performance should improve if fix=T
  expect_gt(jensen_shannon(p,z),jensen_shannon(p,z,fix=T))
  
  # if two vectors differ in the values but not in which elements are(n't) zero,
  # output should be 0 when equiprobabilizing
  expect_equal(jensen_shannon(eq(p),eq(q)),0)
  
})

# compareMatrices
test_that("compareMatrices", {
  
  g <- c("","A","B","C","A, B","A, C","B, C","A, B, C")
  
  A <- matrix(c(c(0,.25,0,.75,0,0,0,0,0,0),
                c(0,0,0,0,0,.95,0,0,0,.05),
                c(0,0,0,0,0,0,0,0,0,0),
                c(0,0,0,0,0,1,0,0,0,0),
                c(0,0,0,0,0,0,0,0,0,0),
                c(0,0,0,0,0,0,0,1,0,0),
                c(0,0,0,0,0,0,0,0,0,0),
                c(0,0,0,0,0,0,0,0,0,1)),
              nrow=length(g),byrow=T)
  rownames(A) <- g
  colnames(A) <- c(g,"none","end")
  
  B <- matrix(c(c(0,.35,0,.65,0,0,0,0,0,0),
                c(0,0,0,0,0,.85,0,0,0,.15),
                c(0,0,0,0,0,0,0,0,0,0),
                c(0,0,0,0,0,.5,.5,0,0,0),
                c(0,0,0,0,0,0,0,0,0,0),
                c(0,0,0,0,0,0,0,1,0,0),
                c(0,0,0,0,0,0,0,1,0,0),
                c(0,0,0,0,0,0,0,0,0,1)),
              nrow=length(g),byrow=T)
  rownames(B) <- g
  colnames(B) <- c(g,"none","end")
  
  s <- compareMatrices(A,B)
  
  
  A <- matrix(c(c(0,.25,0,.75,0,0,0,0,0,0),
                c(0,0,0,0,0,.95,0,0,0,.05),
                c(0,0,0,0,0,0,0,0,0,0),
                c(0,0,0,0,0,1,0,0,0,0),
                c(0,0,0,0,0,0,0,0,0,0),
                c(0,0,0,0,0,0,0,1,0,0),
                c(0,0,0,0,0,0,0,0,0,0),
                c(0,0,0,0,0,0,0,0,0,1)),
              nrow=length(g),byrow=T)
  rownames(A) <- g
  colnames(A) <- c(g,"none","end")
  
  B <- matrix(c(c(0,0,.35,.65,0,0,0,0,0,0),
                c(0,0,0,0,0,0,0,0,0,0),
                c(0,0,0,0,.85,0,0,0,0,.15),
                c(0,0,0,0,0,.5,.5,0,0,0),
                c(0,0,0,0,0,0,0,1,0,0),
                c(0,0,0,0,0,0,0,1,0,0),
                c(0,0,0,0,0,0,0,1,0,0),
                c(0,0,0,0,0,0,0,0,0,1)),
              nrow=length(g),byrow=T)
  rownames(B) <- g
  colnames(B) <- c(g,"none","end")
  
  s <- compareMatrices(A,B)
  
})


## FIXED code for extremely rare borderline

# Kullback-Leibler distance:  Dkl = D(x || y)
my_kldiv_2 <- function(x, y) {
  
  ## Based on https://stackoverflow.com/a/27953616
  ## with modification for removal of 0 in numerator of Kullback-Leibler as per
  ## https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence#Definition
    ## Using log2 so that max is 1. See wikipedia or
    ##  Lin, J. (1991). "Divergence measures based on the shannon entropy"
    ## IEEE Transactions on Information Theory. 37 (1): 145–151. 
    
  ## remove zeros in x
  rmx <- which(x == 0)
    if(length(rmx)) {
        tmp <- sum(x[-rmx] * log2(x[-rmx]/y[-rmx]))
    } else {
        tmp <- sum(x * log2(x/y))
    }
    if(tmp < 0) {
        warning(paste0("Kullback-Leibler < 0; value = ", tmp,
                       ". Setting to 0."))
        tmp <- 0
    }
    return(tmp)
}

# Jensen-Shannon entropy
jensen_shannon_2 <- function(p, q, fix=F) {
  
  ## FIXME: could have checked none in m is 0?
  ## well, one could pass a silly c(a, b, 0), c(e, f, 0)
  ## and should we not allow for that?
  ## KL will not crash (we remove those, as we should, and they are really
  ## irrelevant anyway)
  
  ## FIXME: generalized function so it can handle all-zeros vectors as input
  
  pq <- c(p, q)
  if(any(is.na(pq))) return(NA)
  if(any(pq > 1)) stop("some prob. > 1")
  if(any(pq < 0)) stop("some prob. < 0")
  stopifnot(isTRUE(all.equal(sum(p), 1.0)) | isTRUE(all.equal(sum(p), 0.0)))
  stopifnot(isTRUE(all.equal(sum(q), 1.0)) | isTRUE(all.equal(sum(q), 0.0)))
  stopifnot(identical(length(p), length(q)))
  m <- 0.5 * (p + q)
  
  if (all(pq==0)) return(0) # minimum possible value
  
  if ((any(p>0) & all(q==0)) | (all(p==0) & any(q>0))) {
    # return(1) #max. value
    
    ### FIXME: the output of this (p is all 0 and q is finite or vice-versa)
    ### now depends on the "eq" flag: if F it will return 1, if T it will
    ### equiprobabilize the all 0 input and calculate jensen_shannon as usual
    
    if(fix) {
      if(any(p>0)) q <- rep(1/length(q),length(q))
      if(all(p==0)) p <- rep(1/length(p),length(p))
      return(jensen_shannon_2(p,q))
    } else {
      return(1) # max. possible value
    }
  }
  
  if(any(p>0) & any(q>0)) return(0.5 *  (my_kldiv_2(p, m) + my_kldiv_2(q, m)))

    
    ## how do we know we never fail to fall into those cases?
    stop("JS undefined")
    
  ## rmp <- which(p == 0)
  ## rmq <- which(q == 0)
  ## kl1 <- sum(p[-rmp] * log(p[-rmp]/m[-rmp]))
  ## kl2 <- sum(q[-rmq] * log(q[-rmq]/m[-rmq]))
  ## return(0.5 * (kl1 + kl2))
  ## <- 0.5 * (sum(p * log(p / m)) + sum(q * log(q / m)))
  
} 


####
# jensen_shannon_2 tests
local({
  ## ## understand KLdiv output
  ## p <- abs(rnorm(15, 0, sd = 3)^2)
  ## p <- p/sum(p)
  ## q <- runif(15)
  ## q <- q/sum(q)
  ## KLdiv(cbind(p, q))
  ## my_kldiv(p, q) ## 1st row, 2nd col
  ## my_kldiv(p, q)
  
  p <- c(0.1, 0.3, 0.6, 0, 0)
  q <- c(0, 0, 0, 0.5, 0.5)
  ## stopifnot(isTRUE(all.equal(jensen_shannon_2(p, q), log(2))))
  ## Recall we are now using log2
  stopifnot(isTRUE(all.equal(jensen_shannon_2(p, q), 1)))
  
  stopifnot(isTRUE(all.equal(jensen_shannon_2(p, p), 0)))
  stopifnot(isTRUE(all.equal(jensen_shannon_2(q, q), 0)))
  
  ## example from stackoverflow
  ## https://stackoverflow.com/questions/11226627/jensen-shannon-divergence-in-r
  p <- c(0.00029421, 0.42837957, 0.1371827, 0.00029419, 0.00029419,
         0.40526004, 0.02741252, 0.00029422, 0.00029417, 0.00029418)
  
  q <- c(0.00476199, 0.004762, 0.004762, 0.00476202, 0.95714168,
         0.00476213, 0.00476212, 0.00476202, 0.00476202, 0.00476202)
  
  ## stopifnot(isTRUE(all.equal(jensen_shannon_2(p, q),
  ##                            0.6457538)))
  ## Again, log base 2
  stopifnot(isTRUE(all.equal(jensen_shannon_2(p, q),
                               0.6457538/log(2))))
    
  allStats(p, q)
  rm(p, q)
})

# jensen_shannon (with or without equiprobabilization, with or without fixing)
test_that("jensen_shannon (many scenarios)", {
  
  p <- c(0.1,0,0,0.7,0,0.2)
  q <- c(0.2,0,0,0.1,0,0.7)
  r <- c(0,0,0,0.9,0.1,0)
  z <- c(0,0,0,0,0,0)
  
  # general case: if both inputs are different and not all-zero, output should be
  # a number >0
  expect_gt(jensen_shannon_2(p,q),0)
  
  # if both inputs are all-zero, output should be 0
  expect_equal(jensen_shannon_2(z,z),0)
  
  # if one of the inputs is all-zero and the other isn't, with fix=F should return 1
  expect_equal(jensen_shannon_2(p,z),1)
  
  # previous case: performance should improve if fix=T
  expect_gt(jensen_shannon_2(p,z),jensen_shannon_2(p,z,fix=T))
  
  # if two vectors differ in the values but not in which elements are(n't) zero,
  # output should be 0 when equiprobabilizing
  expect_equal(jensen_shannon_2(eq(p),eq(q)),0)
  
})


test_that("jensen_shannon with the negative KL", {
    p <- c(0.51612903225806450180, 0.48387096774193544269)
    q <- c(0.5161290322580645018, 0.4838709677419354982)

    expect_equal(jensen_shannon_2(p, q), 0)
    expect_true(jensen_shannon(p, q) < 0 )
})
