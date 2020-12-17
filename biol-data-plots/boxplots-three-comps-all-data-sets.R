## Boxplots similar CE, similar TD.
source("../oncoFunctions.R")
library(ggplot2)
library(RColorBrewer)
library(magrittr)
library(Matrix)
library(dplyr)
library(data.table)
library(cowplot)

## If set to false, use the defaultDirs
##  (permits non-interactive work)
ask_user_for_dirs <- FALSE
# directories to load/save from
if(ask_user_for_dirs) {
loadDirectory <- askDir(defaultDir=".",
                        message="Enter directory where biol-data.rds table is stored.")

loadDirectory_g <- askDir(defaultDir="../run-biol-examples",
                        message="Enter directory where genots.RData file is stored.")

saveDirectory <- askDir(defaultDir="./plots",
                        message="Enter directory to save output to.")
} else {
    loadDirectory <- "."
    loadDirectory_g <- "../run-biol-examples"
    saveDirectory <- "./plots"
}


# load data
data <- readRDS(file.path(loadDirectory,"biol-data.rds"))

# create saveDirectory if it doesn't exist
dir.create(saveDirectory,showWarnings=F)


data00 <- data
## Recall JS divergence, not JS distance was used. We will transform it below.
## But rm columns we certainly will not use
data <- data[, c(1:17, 96, 42, 50, 60)]
data$"MHN-vs-MHN_td" <- sqrt(data$"MHN-vs-MHN_td")
data$"CBN-vs-CBN_td" <- sqrt(data$"CBN-vs-CBN_td")
data$"CBN-vs-MHN" <- sqrt(data$"CBN-vs-MHN")
data$"CBN_td-vs-MHN_td" <- sqrt(data$"CBN_td-vs-MHN_td")


## Load datasets info
load(file.path(loadDirectory_g,"genots.RData"))

## Get genotype frequencies and weights
tmp1 <- lapply(1:length(genots),
               function(i) {
                   tmpz <- cbind(genots[[i]], dataset_id = names(genots)[i])
                   tmpz$Nsamples <- sum(tmpz$Freq)
                   return(tmpz)
               })

genot_freqs <- do.call(rbind, tmp1)
## simpler merging below
colnames(genot_freqs)[colnames(genot_freqs) == "Genotype"] <- "sourceGenotype"

stopifnot(nrow(genot_freqs) == nrow(data))

## Used below 
ngenes <- unlist(lapply(genots,
       function(z) {
           tmp <- unique(unlist(strsplit(z$Genotype, ", ")))
           tmp <- tmp[!tmp == "WT"]
           return(length(tmp))
       }
       ))
ngenes <- data.frame(dataset_id = names(ngenes), nGenes = ngenes)



data3 <- dplyr::left_join(data, genot_freqs, by = c("dataset_id", "sourceGenotype"))

stopifnot(nrow(data3) == nrow(data))

data <- data3
## Weights so that all datasets have the same sum of weights
data$W_equi_dataset <-  data$Freq/data$Nsamples
tapply(data$W_equi_dataset, data$dataset_id, sum) ## Yes, 1

# fix gene naming
data$sourceGenotype <- gsub("^G_","",data$sourceGenotype)
data$sourceGenotype <- gsub(", G_",", ",data$sourceGenotype)
data$sourceGenotype <- gsub("_n_","-",data$sourceGenotype)
data$sourceGenotype <- gsub("_p_","+",data$sourceGenotype)
data$sourceGenotype <- gsub("_Pl_","(",data$sourceGenotype)
data$sourceGenotype <- gsub("_Pr_",")",data$sourceGenotype)



data0x <- data
data <- merge(data, ngenes, by="dataset_id", all=T)

stopifnot(nrow(data0x) == nrow(data))
## rm(data0x)
rm(ngenes)

data$sourceGenotype_nMutRel <- data$sourceGenotype_nMut/data$nGenes


# remove "schill" tag from dataset IDs
data$dataset_id <- gsub("gliob_schill_15","GBM_261",data$dataset_id)
data$dataset_id <- gsub("renal_schill","Renal_CGH",data$dataset_id)
data$dataset_id <- gsub("breast_schill","Breast_CGH",data$dataset_id)
data$dataset_id <- gsub("colon_schill","Colon_CGH",data$dataset_id)




## remove unused datasets

data <- data[ !data$dataset_id  %in% c("gliob_schill",
                                       "GBM_261"),
             ] 


## Double check genes
aggregate(nGenes ~ dataset_id, FUN = mean, data = data)
## Compare with ./run-biol-examples/existing-genotypes-and-sample-sizes.Rout

## Keep only columns we care about
data4 <- data[, c("dataset_id",
                  "dataset_num",
                  "sourceGenotype",
                  "sourceGenotype_nMut",
                  "sourceGenotype_nMutRel",                  
                  "Freq",
                  "Nsamples",
                  "nGenes",
                  "W_equi_dataset",
                  "CBN-vs-MHN",
                  "CBN_td-vs-MHN_td",
                  "MHN-vs-MHN_td"
                  )]

data4 <- data4[data4$sourceGenotype_nMutRel < 1, ]

## yes, hack, but faster
data5 <- rbind(cbind(data4, jsdiff = data4$"CBN-vs-MHN", comp = "CBN-vs-MHN"),
               cbind(data4, jsdiff = data4$"CBN_td-vs-MHN_td", comp = "CBN_td-vs-MHN_td"),
               cbind(data4, jsdiff = data4$"MHN-vs-MHN_td", comp = "MHN-vs-MHN_td"))

data5$comp <- factor(data5$comp, levels = c("CBN-vs-MHN", "CBN_td-vs-MHN_td", "MHN-vs-MHN_td"))

data5$dataset_id_r <- reorder(data5$dataset_id, desc(data5$dataset_id))

data5$dataset_id_r2 <- paste0(data5$dataset_id_r,
                              " (n = ", data5$Nsamples,
                              "; g = ", data5$nGenes, 
                              ")") 
data5$dataset_id_r2 <- reorder(data5$dataset_id_r2, desc(data5$dataset_id_r2))


pl1 <- brewer.pal(8, "Dark2")
pl2 <- colorRampPalette(pl1)(length(unique(data5$dataset_id)))
set.seed(2)
pl2 <- sample(pl2)
## pie(rep(1, length(pl2)), col = pl2 , main="")

gg1 <- ggplot(data = data5) +
    geom_boxplot(aes(y = dataset_id_r2, x = jsdiff, color = dataset_id)) +
    geom_jitter(aes(y = dataset_id_r2, x = jsdiff, color = dataset_id), width = 0) +
    facet_wrap(~comp) +
    xlab("JS between method predictions") + ylab("Data set") +
    scale_color_manual(values = pl2,
                        guide = "none") ## + theme_light()

ggsave(filename = file.path(saveDirectory,
                    "boxplots-three-comps_all_datasets_combined.pdf"),
       plot = gg1, dpi = 600,
       width = 11, height = 14)

