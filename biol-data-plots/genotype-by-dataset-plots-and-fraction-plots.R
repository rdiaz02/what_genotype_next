## Individual data set plots for each genotype and fraction CE, fraction TD plots

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

## Same comparisons as in "when-good-conditions.R"

comps_ACROSS <- c("MHN-vs-MHN_td", "CBN-vs-CBN_td")
comps_CE <- c("CBN-vs-MHN") ## , "MCCBN-vs-MHN", "CBN-vs-MCCBN")
comps_TD <- c("CBN_td-vs-MHN_td")

## Recall JS divergence, not JS distance was used. We will transform it below.
## But rm columns we certainly will not use
data00 <- data
data <- data[, c(1:17, 96, 42, 50, 60)]
data$"MHN-vs-MHN_td" <- sqrt(data$"MHN-vs-MHN_td")
data$"CBN-vs-CBN_td" <- sqrt(data$"CBN-vs-CBN_td")
data$"CBN-vs-MHN" <- sqrt(data$"CBN-vs-MHN")
data$"CBN_td-vs-MHN_td" <- sqrt(data$"CBN_td-vs-MHN_td")



## restructure data for easy plotting
## This is like stacking by comparison
##  zz: a row (not row number) of dataframe
structRow2 <- function(zz,
                       comps_ce = comps_CE,
                       comps_td = comps_TD,
                       comps_across = comps_ACROSS) {
    zz <- zz[,!grepl("givesPrediction_",colnames(zz))]
  
  out <- data.frame(dataset_num = zz$dataset_num,
                    dataset_id = zz$dataset_id,
                    sourceGenotype = zz$sourceGenotype,
                    sourceGenotype_nMut = zz$sourceGenotype_nMut,
                    
                    comp = colnames(zz)[grep("-vs-",colnames(zz))],
                    js = as.numeric(zz[,grep("-vs-",colnames(zz))]))
  
  comp_type <- sapply(out$comp,
                      function(comp) {
                        if (comp %in% comps_ce) return("CE")
                        if (comp %in% comps_td) return("TD")
                        if (comp %in% comps_across) return("CE_TD")
                        else return(NA)
                      })
  stopifnot(is.vector(comp_type))
  out$comp_type <- comp_type
  ## keep relevant comparisons only (i.e. remove _uw methods)
  out <- out[!is.na(out$comp_type),]
  return(out)
}

dataStruct <- lapply(1:nrow(data), function(i) structRow2(data[i, ]))
dataStruct <- do.call(rbind,dataStruct)
data02 <- data
data <- dataStruct
## rm(dataStruct)


# remove NAs
data <- data[!is.na(data$js),]


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
stopifnot(nrow(genot_freqs) == nrow(data02))

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
tapply(data$W_equi_dataset, data$dataset_id, sum) ## Yes, 7: the seven comparisons

# fix gene naming
data$sourceGenotype <- gsub("^G_","",data$sourceGenotype)
data$sourceGenotype <- gsub(", G_",", ",data$sourceGenotype)
data$sourceGenotype <- gsub("_n_","-",data$sourceGenotype)
data$sourceGenotype <- gsub("_p_","+",data$sourceGenotype)
data$sourceGenotype <- gsub("_Pl_","(",data$sourceGenotype)
data$sourceGenotype <- gsub("_Pr_",")",data$sourceGenotype)


## break genotype names into multiple lines
## For plotting
trimGenotypeName <- function(genotype_name) {
  nchar_max <- 20

  genes <- strsplit(genotype_name,", ")[[1]]
  genes <- paste(genes,", ",sep="")
  genes[length(genes)] <- substr(genes[length(genes)],1,
                                 nchar(genes[length(genes)])-2)
  recursiPaste <- function(genes) {
    cumulative_nchar <- cumsum(nchar(genes))
    if (any(cumulative_nchar>nchar_max)) {
      if (any(cumulative_nchar<=nchar_max)) {
        line <- paste(genes[cumulative_nchar <= nchar_max],collapse="")
        line <- substr(line,1,nchar(line)-1)
        line <- paste(line,"\n",sep="")
        line <- c(line,recursiPaste(genes[cumulative_nchar > nchar_max]))
        return(line)
      } else {
        if (length(genes)>1) {
          line <- genes[1]
          line <- substr(line,1,nchar(line)-1)
          line <- paste(line,"\n",sep="")
          line <- c(line,recursiPaste(genes[2:length(genes)]))
          return(line)
        } else {
          line <- genes[1]
          return(line)
        }
      }
    } else {
      line <- paste(genes,collapse="")
      return(line)
    }
  }
  
  return(paste(recursiPaste(genes),collapse=""))
  
}

data$sourceGenotypeTrimmed <- sapply(data$sourceGenotype,trimGenotypeName)

# genotype names as factors (for ordering)
data <- data[order(data$sourceGenotype),]
data$sourceGenotype <- factor(data$sourceGenotype,
                              levels=unique(data$sourceGenotype[order(data$sourceGenotype_nMut)]))

## Is the difference indication of bug? Not yet, as we have not
## averaged over genes
table(data$comp_type)

# comparisons as factors (for ordering)
data$comp_type[data$comp_type=="CE"] <- "CE vs. CE"
data$comp_type[data$comp_type=="TD"] <- "TD vs. TD"
data$comp_type[data$comp_type=="CE_TD"] <- "CE vs. TD"

## data$comp_type <- factor(data$comp_type,
##                          levels=c("CE vs. CE","TD vs. TD","CE vs. TD"))

data$comp_type <- factor(data$comp_type,
                         levels=c("TD vs. TD","CE vs. CE","CE vs. TD"))



data0x <- data
data <- merge(data, ngenes, by="dataset_id", all=T)

stopifnot(nrow(data0x) == nrow(data))
rm(data0x)
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


## Ordering in the call seems strange, but allows comparing for double
## check below
ag1 <- aggregate(js ~ 
                     dataset_num + sourceGenotype_nMut +
                     sourceGenotypeTrimmed + 
                     Freq + Nsamples + W_equi_dataset +
                     nGenes + sourceGenotype_nMutRel +
                     comp_type + sourceGenotype + dataset_id,
                 FUN = mean, data = data)

## Double check 
ag2 <- aggregate(js ~ comp_type + sourceGenotype + dataset_id,
                 FUN = mean, data = data)

stopifnot(identical(nrow(ag1), nrow(ag2)))
stopifnot(identical(nrow(ag1$js), nrow(ag2$js)))

olddata <- data

data <- ag1

## Now, same number
table(data$comp_type)

nrow(data)
## remove the last genotype
data <- data[data$sourceGenotype_nMutRel < 1, ]


## This sucks badly, as we are undoing things we did above in structRow2
## but it is simpler. Put, for each genotype, the comparison in one column

dds <- split(data, f = data$comp_type)
stopifnot(identical(colnames(dds[[1]]), colnames(dds[[2]])))
stopifnot(identical(colnames(dds[[1]]), colnames(dds[[3]])))
stopifnot(colnames(dds[[1]])[9] == "comp_type")
stopifnot(colnames(dds[[1]])[12] == "js")
stopifnot(dds[[1]][, "comp_type"] == "TD vs. TD")
stopifnot(dds[[2]][, "comp_type"] == "CE vs. CE")
stopifnot(dds[[3]][, "comp_type"] == "CE vs. TD")


colnames(dds[[1]])[12] <- "JS_TD_vs_TD"
colnames(dds[[2]])[12] <- "JS_CE_vs_CE"
colnames(dds[[3]])[12] <- "JS_CE_vs_TD"

dds[[1]] <- dds[[1]][, -9]
dds[[2]] <- dds[[2]][, -9]
dds[[3]] <- dds[[3]][, -9]

dds12 <- dplyr::left_join(dds[[1]], dds[[2]])
dds123 <- dplyr::left_join(dds12, dds[[3]])

colnames(dds123)
stopifnot(nrow(dds123) == nrow(dds[[1]]))
stopifnot(ncol(dds123) == 2 + ncol(dds[[1]]))

dt <- data.table(dds123)
dt <- dt[, `:=` (similar_TD1 = (JS_TD_vs_TD <= 0.0677),
                 similar_CE1 = (JS_CE_vs_CE <= 0.0677),
                 different_CE_TD1 = (JS_CE_vs_TD >= 0.7))
                 ]
dt <- dt[, `:=` (all_three_conds = (similar_TD1) & (similar_CE1) & (different_CE_TD1))]


dt[, count := .N, by = .(dataset_id)]

dt3 <- dt[,
          .(
               sum_similar_CE1 = sum(similar_CE1 & different_CE_TD1),
               sum_similar_TD1 = sum(similar_TD1 & different_CE_TD1),               
               ## sure, same as counting
               ## sum_different_CE_TD = sum(different_CE_TD),
               sum_different_CE_TD1 = sum(different_CE_TD1),
               sum_all_three = sum(all_three_conds),
               Ngenotypes = mean(count),
               sampleSize = mean(Nsamples)),
           by = .(dataset_id)]

## Ngenotypes does NOT include the last genotype. Counting done after excluding it
dt4 <- dt3[, `:=` (fraction_similar_CE = sum_similar_CE1/Ngenotypes,
                   fraction_similar_TD = sum_similar_TD1/Ngenotypes,
                   fraction_different_CE_TD1 = sum_different_CE_TD1/Ngenotypes,
                   fraction_all_three = sum_all_three/Ngenotypes)]


## Compare Ngenotypes with this: number of observed genotypes,
lapply(genots, nrow)

## E.g., Pan_pa has 9 Ngenotypes but 10 genotypes. Not always the case, of
## course, since sometimes the last genotype not observed.

dt5 <- dt4
dt5$fraction_similar_CE <- round(100 * dt5$fraction_similar_CE)
dt5$fraction_similar_TD <- round(100 * dt5$fraction_similar_TD)
dt5$fraction_all_three <- round(100 * dt5$fraction_all_three)

dt5$dataset_id_r <- paste0(dt5$dataset_id, " (GT = ", dt5$Ngenotypes, ")")
dt5$dataset_id_r <- reorder(dt5$dataset_id_r, desc(dt5$dataset_id_r))


## I am not sure if the fact that I can use stupid things for "color"
## shows that ggplot is great or that it is a huge hack :)
## and I always end up having to double check the plot output
## It is sooo dangerous!!! 

gg1 <- ggplot(data = dt5) +
    geom_jitter(aes(x = fraction_similar_CE, y = dataset_id_r, color = "A"), width = 0) +
    geom_jitter(aes(x = fraction_similar_TD, y = dataset_id_r, color = "B"), width = 0) +
    geom_point(aes(x = fraction_all_three, y = dataset_id_r, color = "C")) +
    scale_color_manual(name = " ",
                       labels = c("Similar CE", "Similar TD", "Similar CE\nAND\nSimilar TD"),
                       values = c("red", "#56B4E9", "grey43")) +
    xlab("Fraction of genotypes relative to total number of genotypes") +
    ylab("Data set") + theme_light() +
    theme(legend.key.size = unit(1, "cm"))
##gg1
filename <- file.path(saveDirectory,"fraction_CE_TD_similar_all_datasets.pdf")

ggsave(filename = filename,
       plot = gg1, device="pdf",dpi=600,
       width=7.5, height=6)



## pl1 <- brewer.pal(8, "Dark2")
## pl2 <- colorRampPalette(pl1)(length(unique(dt5$dataset_id)))
## set.seed(2)
## pl2 <- sample(pl2)
## ## pie(rep(1, length(pl2)), col = pl2 , main="")

## gg6 <- ggplot(data = dt5) +
##     geom_point(aes(x = sampleSize, y = fraction_similar_CE,
##                    color = dataset_id)) +
##     xlab("Sample size") + ylab("Percentage similar CE") +
##     theme(legend.position = "none") +
##     scale_color_manual(values = pl2)
## ## gg6
## gg7 <- ggplot(data = dt5) +
##     geom_point(aes(x = sampleSize, y = fraction_similar_TD,
##                    color = dataset_id)) +
##     xlab("Sample size") + ylab("Percentage similar TD ") +
##     theme(legend.title = element_blank()) +
##     guides(color = guide_legend(ncol = 1)) +
##     scale_color_manual(values = pl2)
## ## gg7
## gg7 <- gg7 + theme_light() + theme(legend.title = element_blank())
## gg6 <- gg6 + theme_light() + theme(legend.position = "none")

## gg8 <- plot_grid(gg6, gg7, nrow = 1, ncol = 2,
##                  rel_widths = c(0.85, 1))



## ggsave(filename = file.path(saveDirectory,
##                             "fraction_CE_TD_similar_all_datasets_sample_size.pdf"),
##        plot = gg8, dpi = 600,
##        width = 12, height = 6.5)





## gg8b <- plot_grid(gg6, gg7, nrow = 1, ncol = 2,
##                   rel_widths = c(0.85, 1),
##                   labels = c("B", "C"))

## gg9 <- plot_grid(gg4,
##                  NULL,
##                  gg8b,
##                  nrow = 3,
##                  rel_heights = c(1, 0.1, 0.8), 
##                  labels = c("A", "")
##                  ## , scale = c(1, 0.85)
##                  )
## save_plot(file.path(saveDirectory,
##                     "fraction_CE_TD_similar_all_datasets_combined.pdf"),
##           gg9,
##           base_width = 12, base_height = 14
##           )




## datset number and data frame  -> figure
makePlot2 <- function(dataset_num, dataz) {
  # dataset
  x <- dataz[dataz$dataset_num==dataset_num,]

    ## RDU: 
    ## Number of columns for A4. 
    ## A4: 21 * 29 : 1.4
    ## Each approx square

    n_panels0 <- length(unique(x$sourceGenotypeTrimmed))
    n_col0 <- round(sqrt(n_panels0/1.4))
    
    
  p <- ggplot(x,aes(x=comp_type,y=js,color=comp_type)) +
      ## geom_boxplot(fill=NA) + ## one or two points. No need for boxplot
      ##    geom_point() +
      geom_jitter(width=0.1, size = 2, height = 0) +
    ylim(-0.01,1.01) +
    scale_color_brewer(type="seq",palette="Set1") +
    theme_light() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom") +
    labs(color="Comparison",
         y="Jensen-Shannon distance") +
    ggtitle(unique(x$dataset_id)) +
    facet_wrap(~sourceGenotypeTrimmed, ncol = n_col0)
  
  # scaling factor
  width_per_panel <- 2.5
  height_per_panel <- 1.8
  n_panels <- length(unique(x$sourceGenotypeTrimmed))
  n_rows <- p %>%
    ggplot2::ggplot_build() %>%
    magrittr::extract2('layout') %>% 
    magrittr::extract2('layout') %>%
    magrittr::extract2('ROW') %>%
    unique() %>%
    length()
    n_cols <- ceiling(n_panels/n_rows)

  
    
  title_lines <-  1 + max(str_count(x$sourceGenotypeTrimmed,"\n"))
  n_rows <- n_rows*(1 + 1/(10-1)*(title_lines-1)) # ~10 lines = one panel's height
  

    filename <- paste0("CE_TD_similar_plots_dataset_",
                      unique(x$dataset_id),".pdf",collapse="",sep="")
    filename <- file.path(saveDirectory,filename)
    ggsave(filename = filename,plot = p,device = "pdf",dpi = 600,
         width = n_cols*width_per_panel,height = n_rows*height_per_panel,
         limitsize = FALSE)
  
}

olddata2 <- olddata[olddata$sourceGenotype_nMutRel < 1, ]

individual_plots <- TRUE
if(individual_plots) {
pboptions(type="txt")
pblapply(unique(olddata2$dataset_num),
                function(i) makePlot2(i, dataz = olddata2),
                cl=detectCores())
}




