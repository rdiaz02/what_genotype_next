# ----------

# This script takes transition matrices generated with the scripts:
#   next-genotype_transitionMatrix-sim.R
#   next-genotype_transitionMatrix-cpm.R
#   next-genotype_transitionMatrix-timeDiscretizedCPMs.R
# And combines all the data in a single file (generates one file per fitness
# landscape ID).

# ----------

source("oncoFunctions.R")

# directories to load from/save to
simDirectory <- askDir(defaultDir="./sim",
                           message="Enter directory where .RData files containing transition matrices extracted from OncoSimulR simulations are.")
cpmDirectory <- askDir(defaultDir="./cpm",
                           message="Enter directory where .RData files containing formatted transition matrices extracted from CPM outputs are.")
saveDirectory <- askDir(defaultDir="./data",
                        message="Enter directory to save output to.")

# create save directory (if it doesn't exist)
dir.create(saveDirectory,showWarnings=F)

# list all CPM sub-directories under the cpmDirectory
cpm <- list.dirs(cpmDirectory,full.names=F,rec=F)

# list all files
cat("Checking provided directories")
cat("\n")
simFiles <- list.files(simDirectory,full.names=T,rec=F)
simFiles <- simFiles[!file.info(simFiles)$isdir]
cpmFiles <- vector(mode="list")
for (i in 1:length(cpm)) {
  files <- list.files(file.path(cpmDirectory,cpm[i]),
                      full.names=T,rec=F)
  files <- files[!file.info(files)$isdir]
  cpmFiles[[cpm[i]]] <- files
}
cpmFiles_listed <- unlist(cpmFiles)
names(cpmFiles_listed) <- NULL

# list all IDs
flIDs <- sapply(simFiles,flIDFromFile)
names(flIDs) <- NULL

# list all size_splits
size_splits <- sapply(cpmFiles_listed,sizeSplitFromFile)
names(size_splits) <- NULL
size_splits <- unique(size_splits)

# list all detect regimes
detects <- sapply(cpmFiles_listed,detectFromFile)
names(detects) <- NULL
detects <- unique(detects)

# match landscape IDs and types
flTypes <- data.frame(ID=sapply(cpmFiles_listed,flIDFromFile),
                typeLandscape=sapply(cpmFiles_listed,flTypeFromFile))
flTypes <- unique(flTypes)
rownames(flTypes) <- NULL

# check if IDs match or there are missing files
simIDs <- sapply(simFiles,flIDFromFile)
cpmIDs <- sapply(cpmFiles_listed,flIDFromFile)
cpmIDs <- as.data.frame(table(cpmIDs))

isMatched <- rep(NA,length(simIDs))
for (i in 1:length(simIDs)) {
  isMatched[i] <- simIDs[i] %in% cpmIDs$cpmIDs
}
if(sum(isMatched)==length(isMatched)) {
  cat("  > Landscape IDs matched")
} else {
  cat("  > WARNING: encountered landscape ID mismatch")
}
cat("\n")

if(sum(cpmIDs$Freq==length(size_splits)*length(detects)*length(cpm))==dim(cpmIDs)[1]) {
  cat("  > All CPM files located")
} else {
  stop("  > WARNING: missing files")
}
cat("\n")
cat("\n")

# group all matrices for a given ID under a same object
cat("Structuring data")
cat("\n")
pboptions(type="txt")
x <- pblapply(simIDs,
              function(ID) {

                # initialize output object
                out <- vector(mode="list")
                out[["ID"]] <- as.character(ID)
                out[["typeLandscape"]] <- flTypes$typeLandscape[flTypes$ID==ID]
                out[["nGenes"]] <- numGenesFromFile(simFiles[grepl(ID,simFiles)])
                
                # "sim" matrix
                load(simFiles[grepl(ID,simFiles)])
                out[["sim"]] <- list(transitionMatrix=transitionMatrix,
                                     timesInPOM=timesInPOM)
                
                # null model
                out[["null"]] <- list(transitionMatrix=
                                        nullMatrix(transitionMatrix)$norm)
                
                # loop through CPMs, size_splits and detection regimes (nest)
                cpm <- names(cpmFiles)
                for (i in 1:length(cpm)) { # cpm index
                  files <- cpmFiles[[cpm[i]]]
                  files <- files[grepl(ID,files)]
                  for (j in 1:length(size_splits)) { # "size_split" index
                    for (k in 1:length(detects)) { # "detect" index
                      n <- grepl(size_splits[j],files) & grepl(detects[k],files)
                      load(files[n])
                      out[["cpm"]][[cpm[i]]][[size_splits[j]]][[detects[k]]] <-
                        transitionMatrix
                    }
                  }
                }
                
                # save output
                data <- out
                outFile <- file.path(saveDirectory,
                                     paste(ID,".RData",sep=""))
                save(data,file=outFile)
                
                return(0)
                
              },
              cl=detectCores())







