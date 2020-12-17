# ----------

# This script takes output files from CPMs and generates a transition matrix
# structured in such a way that it can be directly compared to the matrices
# obtained from OncoSimulR simulation files. It saves a .RData object containing
# five transition matrices (one per replicate), in a subdirectory within the
# one specified by saveDirectory. Each subdirectory corresponds to a different
# CPM.

# ----------

source("oncoFunctions.R")

# directories to load from/save to
loadDirectory <- askDir(defaultDir="/home/ramon/No-Backuppc/cpm-analyzed-all-selected",
                            message="Enter directory where .rds files containing evolutionary paths predicted by CPMs are.")
saveDirectory <- askDir(defaultDir="./cpm",
                        message="Enter directory to save output to.")

# list CPM output files
cat("Checking provided directory")
cat("\n")
cat(paste("  > ",loadDirectory,sep=""))
cat("\n")
files <- list.files(loadDirectory,
                    pattern="all_paths_CPMs(.*?).rds",full.names=T)
cat(paste("  > Found",length(files),"files"))
cat("\n")
cat("\n")


# generate all genotype names, IDs, etc
cat("Generating genotype names and IDs")
cat("\n")
cat("\n")
numGenes <- as.numeric(unlist(sapply(files,numGenesFromFile)))
aux <- makeGenotypes(numGenes)

# available CPMs (uw- indicates the unweighed transition matrix, fgraph, to be used)
cpm <- c("OT","CAPRESE","CAPRI_BIC","CAPRI_AIC","CBN_ot","MCCBN","OT","CBN_ot","MCCBN")
cpmDirNames <- c(tolower(cpm[1:6]),paste("uw-",c("ot","cbn_ot","mccbn"),sep=""))

# which matrix should be used for each cpm?
cpmMatrix <- c("fgraph","trans_mat_genots")[c(2,1,1,1,2,2,1,1,1)]

# create directories if not already there
dir.create(saveDirectory,showWarnings=F)
for (i in 1:length(cpm)) {
  dir <- paste(saveDirectory,cpmDirNames[i],sep="/")
  if (!dir.exists(dir)) dir.create(dir)
}

# sweep through files
cat("Processing files")
cat("\n")
pboptions(type="txt")
x <- pblapply(files,
              function(file, f4000 = filesANALYZED) {

                # load i-th file
                cpmObj <- readRDS(file)
                numReplicates <- length(cpmObj$out_paths_genots)

                  
                # full list of genotypes
                numGenes <- as.numeric(cpmObj$ngenes)
                numGenotypes <- aux[[as.character(numGenes)]]$numGenotypes
                genotypeNames <- aux[[as.character(numGenes)]]$genotypeNames
                genotypeIDs <- aux[[as.character(numGenes)]]$genotypeIDs
                genotypeMuts <- aux[[as.character(numGenes)]]$genotypeMuts
                
                # loop through all CPMs
                for (j in 1:length(cpmDirNames)) {
                  
                  # initialize transition matrix
                  transitionMatrix <- vector(mode="list",
                                             length=numReplicates)
                  
                  # loop through replicates
                  for (k in 1:numReplicates) {
                    t_jk <- cpmObj$out_paths_genots[[k]][[cpm[j]]][[cpmMatrix[j]]]
                    transitionMatrix[[k]] <-
                      structTransitionMatrix(t_jk,
                                             allGenotypes=genotypeNames)
                  }
                  
                  # save transition matrix in corresponding directory
                  outFile <- file.path(saveDirectory,cpmDirNames[j],
                                       paste("transMatrix_",
                                             str_match(basename(file),
                                                       "(.*?).rds")[2],
                                             ".RData",sep=""))
                  save(transitionMatrix,file=outFile)
                  
                }
                
                return(0)
              
              },
              cl=detectCores())

