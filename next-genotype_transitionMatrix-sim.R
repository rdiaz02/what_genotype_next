# ----------

# This script takes OncoSimulR output objects and generates the transition matrix
# corresponding to them. Each entry (i,j) of the matrix produced is an answer to the
# question: "given that a genotype i has been observed in the POM, what is the
# probability that the genotype in the LOD that has one mutation more than i is j?"

# Note that the answer to this question could be "there is no genotype in the LOD
# with one mutation more than the observed genotype in the POM". This could be
# because of one of two reasons: one, the observed genotype is the final genotype
# that has fixated, which corresponds to the "end" column in the transition matrix),
# or two, the genotype that ultimately fixates has less mutations than the observed
# one (i.e. the genotype that was the most abundant at the time of observation does
# not belong to the LOD and was ultimately outcompeted by a genotype with less
# mutations), which corresponds to the "none" column in the matrix.

# ----------

source("oncoFunctions.R")

# directories to load from/save to
loadDirectory <- askDir(defaultDir="/Disk2/ramon/No-Backuppc",
                            message="Enter parent directory of folders /selected-simulations-A and /selected-simulations-B containing OncoSimulR simulation files.")
saveDirectory <- askDir(defaultDir="./sim",
                        message="Enter directory to save output to.")

# create directory if needed
dir.create(saveDirectory,showWarnings=F)

# list simulation files
cat("Checking provided directory")
cat("\n")
cat(paste("  > ",loadDirectory,sep=""))
cat("\n")
dirs <- list.dirs(path=loadDirectory,full.names=T)

## FIXME: why grep after listing, instead of, by decree,
## pasting those two subdirectories?
dirs <- dirs[c(grep("selected-simulations-A",dirs),
               grep("selected-simulations-B",dirs))]
files <- list.files(dirs,pattern="simobj",full.names=T)
cat(paste("  > Found",length(files),"files"))
cat("\n")
cat("\n")

# generate all genotype names, IDs, etc
cat("Generating genotype names and IDs")
cat("\n")
cat("\n")
numGenes <- as.numeric(unlist(sapply(files,numGenesFromFile)))
aux <- makeGenotypes(numGenes)

# sweep through files
cat("Processing files")
cat("\n")
pboptions(type="txt")
x <- pblapply(files,
              function(file) {
  
                # load file, get POMs and LODs
                load(file)
                poms <- POM(simo)
                lods <- LOD(simo)

                 
                # number of genes and genotype names
                numGenes <- numGenesFromFile(file)
                numGenotypes <- aux[[as.character(numGenes)]]$numGenotypes
                genotypeNames <- aux[[as.character(numGenes)]]$genotypeNames
                genotypeIDs <- aux[[as.character(numGenes)]]$genotypeIDs
                genotypeMuts <- aux[[as.character(numGenes)]]$genotypeMuts
                
                # initialize transition matrix
                transitionMatrix <- matrix(0,nrow=numGenotypes,
                                           ncol=numGenotypes+2)
                rownames(transitionMatrix) <- genotypeNames
                colnames(transitionMatrix) <- c(genotypeNames,"none","end")
                
                # initialize frequency of genotypes
                timesInPOM <- rep(0,numGenotypes)
                names(timesInPOM) <- genotypeNames
                
                # extract transition matrix and genotype frequency from LODs/POMs
                for (j in 1:length(lods)) {
                  transitionMatrix <- transitionMatrix +
                    nextInLOD_transitionMatrix(lods[[j]],poms[[j]],
                                               allGenotypes=
                                                 genotypeNames)
                  timesInPOM <- timesInPOM +
                    whichInPOM(genotypeNames,poms[[j]])
                }
                
                # normalize by genotype appearances in POMs
                tNorm <- do.call(cbind,replicate(ncol(transitionMatrix),
                                                 timesInPOM,simplify=F))
                tNorm[tNorm==0] <- 1
                transitionMatrix <- transitionMatrix/tNorm
                
                # save data
                outFile <- file.path(saveDirectory,
                                     paste("transMatrix_",
                                           str_match(basename(file),
                                                     "(.*?).RData.gz")[2],
                                           ".RData",sep=""))
                save(transitionMatrix,timesInPOM,file=outFile)
                
                return(0)
                
              },
              cl=detectCores())

warnings()
