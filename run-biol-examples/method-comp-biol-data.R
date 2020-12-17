# ----------

# Analyze biological data

## BEWARE:!!! Output is JS divergence (i.e., not its square root, not the
## distance)

## BEWARE: DO NOT use the output from the last genotype, as the TD vs CE
## comparisons are not performed correctly. This is irrelevant, since
## those are never used.

# ----------

source("../oncoFunctions.R")
library(Matrix)

# directories to load/save from
loadDirectory <- askDir(defaultDir="./run-biol-examples/RDatas-and-Routs-runs-genes-renamed",
                        message="Enter directory where data are stored.")
saveDirectory <- askDir(defaultDir="../biol-data-plots",
                        message="Enter directory to save output to.")

# create save directory (if it doesn't exist)
dir.create(saveDirectory,showWarnings=F)

# methods to compare
methods <- c("OT",
             "OT_uw",
             "CBN",
             "CBN_td",
             "CBN_uw",
             "MCCBN",
             "MCCBN_td",
             "MCCBN_uw",
             "CAPRESE",
             "CAPRI_BIC",
             "CAPRI_AIC",
             "MHN",
             "MHN_td")

# list files
files <- list.files(loadDirectory,pattern="transition_matrices_num",
                    full.names=T)

# get dataset numbers and names
dataNum <- as.numeric(str_match(files,"num\\s*(.*?)\\s*_dataname")[,2])
dataName <- str_match(files,"dataname_\\s*(.*?)\\s*_cores")[,2]

names(dataNum) <- files
names(dataName) <- files

# get info on datasets
load(file.path(loadDirectory,"..","genots.RData"))
genots_data <- genots
rm(genots)

# this function takes the name of a dataset as input and performs all method-to-method
# comparisons, then returns them as a data frame
# it is intended to be applied over all datasets
do_this_for_every_dataset <- function(file) {
  
  # load dataset  
  load(file)
  transmat <- get(ls(pattern="transition_matrices"))
  rm(list=ls(pattern="transition_matrices"))
  
  ### FIXME: typo: in list names "OT_u" should be "OT_uw"
  names(transmat)[names(transmat)=="OT_u"] <- "OT_uw"
  
  # list genes and genotypes in dataset
  genots <- genots_data[[dataName[file]]][,1]
  genes <- unique(unlist(strsplit(genots,", ")))
  genes <- genes[!genes=="WT"]
  
  ### FIXME: need to be careful with gene ordering in genotype names. E.g. genotype 
  ### "ATP2B2, TP53" can sometimes appear as "TP53, ATP2B2". The function orderGenotypes()
  ### fixes this (originally developed to deal with genotype naming such as "B, A_C" together
  ### with the unfuseGenotype() function)
  genots <- orderGenotype(genots)$orderedGenotype
  
  # check which methods provide a prediction and which don't
  givesPrediction <- function(method) {
    
    # load matrix for input method
    A <- as.matrix(transmat[[method]])
    
    # again, watch out for gene ordering in genotype names
    rownames(A) <- orderGenotype(rownames(A))$orderedGenotype
    
    # check which observable genotypes (genots, genotypes with freq>0 in the sample) are deemed
    # accessible by the method (appear in the transition matrix)
    return(genots %in% rownames(A)[rowSums(A)>0])
    
  }
  givesPrediction <- as.data.frame(sapply(methods,givesPrediction))
  
  # initialize data frame to store method comparisons
  comps <- givesPrediction
  colnames(comps) <- paste("givesPrediction_",colnames(comps),sep="")
  comps <- cbind(dataset_num=dataNum[file],
                 dataset_id=dataName[file],
                 sourceGenotype=genots,
                 sourceGenotype_nMut=nMut(genots),
                 comps,
                 givesPrediction_fraction=rowSums(comps)/length(methods),
                 row.names=NULL)
  
  # function to compare transition matrices of biological data for different methods
  compareMethodsBiol <- function(method1,method2) {
  
    # get matrices to compare (in regular, non-sparse matrix form)
    A <- as.matrix(transmat[[method1]])
    B <- as.matrix(transmat[[method2]])
    
    # preprocess matrices (keep observable genotype rows only, match column names, etc)
    x <- prepareBiolMatrices(A,B,genots)
    A <- x[[1]]
    B <- x[[2]]
    rm(x)
    
    # get jensen-shannon row by row (details of jsBiol() and special cases in oncoFunctions)
    js <- rep(NA,length(genots))
    names(js) <- genots
    for (r in 1:length(genots)) {
      js[r] <- jsBiol(A[r,],B[r,],nmut=nMut(genots[r]),ngenes=length(genes))
    }
    return(js)
  
  }
  
  # loop to compare methods
  for (i in 1:(length(methods)-1)) {
    for (j in (i+1):length(methods)) {
      js <- compareMethodsBiol(methods[i],methods[j])
      comps <- cbind(comps,
                     js)
      colnames(comps)[colnames(comps)=="js"] <- paste(methods[i],"-vs-",methods[j],sep="")
    }
  }
  
  # return results
  return(comps)
  
}

# get comparisons for every dataset
pboptions(type="txt")
comps <- pblapply(files,
                  do_this_for_every_dataset,
                  cl=detectCores())
comps <- do.call(rbind,comps)

# save
saveRDS(comps,file.path(saveDirectory,"biol-data.rds"))

