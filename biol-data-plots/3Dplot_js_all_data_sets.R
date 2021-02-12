## Boxplots similar CE, similar TD.
source("../oncoFunctions.R")
library(ggplot2)
library(RColorBrewer)
library(magrittr)
library(Matrix)
library(dplyr)
library(data.table)
library(cowplot)
library(plotly)
library(processx)

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
load(file.path(loadDirectory_g,"genots.RData.gz"))

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


## pl1 <- brewer.pal(8, "Dark2")
pl1 <- brewer.pal(8, "Set1")
pl1[6] <- "#f0bd00" # manually change yellow to a darker shade for clarity
pl2 <- colorRampPalette(pl1)(length(unique(data5$dataset_id)))
set.seed(2)
pl2 <- sample(pl2)
## pie(rep(1, length(pl2)), col = pl2 , main="")

gg1 <- ggplot(data = data5) +
    geom_boxplot(aes(y = dataset_id_r2, x = jsdiff, color = dataset_id),outlier.shape = NA) + ## outlier.shape set to NA to avoid plotting duplicates with geom_jitter
    geom_jitter(aes(y = dataset_id_r2, x = jsdiff, color = dataset_id), width = 0, size = 0.75) +
    facet_wrap(~comp) +
    ## ylab("Data set") + ## xlab("JS between method predictions") + 
    scale_color_manual(values = pl2,
                        guide = "none") +
    scale_x_continuous(name="JS between method predictions",
                     breaks=c(0,0.5,1), 
                     labels=c("0","0.5","1")) +
    scale_y_discrete(name="",
                     position="right") ## + theme_light()

gg1

data_to_plot <- data5
pp <- plot_ly(x=data_to_plot$`CBN-vs-MHN`,
              y=data_to_plot$`CBN_td-vs-MHN_td`,
              z=data_to_plot$`MHN-vs-MHN_td`,
              type="scatter3d",
              mode="markers",
              color=data_to_plot$dataset_id_r2,
              colors=pl2,
              marker=list(size=3))
pp <- layout(pp,
             scene=list(xaxis=list(title="JS<sub>CBN,MHN</sub>",
                                   range=c(-0.05,1.05)),
                        yaxis=list(title="JS<sub>CBN_td,MHN_td</sub>",
                                   range=c(-0.05,1.05)),
                        zaxis=list(title="JS<sub>MHN,MHN_td</sub>",
                                   range=c(-0.05,1.05)),
                        camera=list(eye=list(x=0.8,
                                             y=1.5,
                                             z=1.1))))
# pp

data_for_m <- data_to_plot[,c("dataset_id_r2","CBN-vs-MHN","CBN_td-vs-MHN_td","MHN-vs-MHN_td")]
data_for_m$dataset_id_r2 <- sapply(1:dim(data_for_m)[1],FUN=function(i) which(data_for_m$dataset_id_r2[i]==levels(data_for_m$dataset_id_r2)))
write.table(data_for_m,file="data_for_m.txt",sep="\t",row.names=F,col.names=F)

colors_for_m <- t(col2rgb(pl2))
write.table(colors_for_m,file="colors_for_m.txt",sep="\t",row.names=F,col.names=F)

# orca(pp,
#      file=file.path(saveDirectory,
#                     "3Dplot_js_all_data_sets.pdf"))

# ggsave(filename = file.path(saveDirectory,
#                     "boxplots-three-comps_all_datasets_combined.pdf"),
#        plot = gg1, dpi = 600,
#        width = 11, height = 14)

