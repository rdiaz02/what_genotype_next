library(OncoSimulR) ## for sampledGenotypes

## This first part, until </here> is from run-all-biol-examples.R which is
## the same as  run-all-biol-examples-cores-cbn.R
## A few lines have been commented

date()

## arguments.script <- commandArgs(trailingOnly = TRUE)
## cores <- as.numeric(arguments.script[1])
## dataset <- as.numeric(arguments.script[2])


source("code-all-methods-2-trans-matrix.R")



## RhpcBLASctl::blas_set_num_threads(cores)
## RhpcBLASctl::omp_set_num_threads(cores)



###############################################

####  The data we use

## simple wrapper to show dimensions and turn into matrix
dim_mat <- function(x) {
    print(dim(x))
    return(as.matrix(x))
}

breast_schill <- dim_mat(readRDS(file="./MHN/data/BreastCancer.rds"))
colon_schill <- dim_mat(readRDS(file="./MHN/data/ColorectalCancer.rds"))
renal_schill <- dim_mat(readRDS(file="./MHN/data/RenalCellCarcinoma.rds"))
gliob_schill0 <- dim_mat(readRDS(file="./MHN/data/Glioblastoma.rds")) 
## Note  0 and 1 as floats. Change.
gliob_schill <- gliob_schill0
storage.mode(gliob_schill) <- "integer"
stopifnot(max(abs(gliob_schill - gliob_schill0)) == 0)
rm(gliob_schill0)
load("all_biol_example_data_from_Diaz_Uriarte_Vasallo_2019.RData")

## rm the five genes that are not in Figure 6
colnames(gliob_schill)[c(16, 19, 17, 14, 15)]
gliob_schill_15 <- gliob_schill[, -c(16, 19, 17, 14, 15)]

all_data <- c(
    all_biol_example_data
   ,
    list(
        breast_schill = breast_schill
    , colon_schill = colon_schill
    , renal_schill = renal_schill
    , gliob_schill = gliob_schill
    , gliob_schill_15 = gliob_schill_15
      )
)

## No data set has a constant column. Good.
any(unlist(lapply(all_data, any_constant_col)))

ncols <- lapply(all_data, ncol)
all_data <- all_data[order(unlist(ncols))]

lapply(all_data, ncol)
## Schill's gliobl occupies position 22
names(all_data)

## CAPRESE and CAPRI have problems with column names that other methods
## handle without trouble. Oh well. So change them.
change_column_names <- function(x) {
    ## None of the character combinations that are used to replace exist
    colnames(x) <- gsub("+", "_p_", colnames(x), fixed = TRUE)
    colnames(x) <- gsub("-", "_n_", colnames(x), fixed = TRUE)
    colnames(x) <- gsub("(", "_Pl_", colnames(x), fixed = TRUE)
    colnames(x) <- gsub(")", "_Pr_", colnames(x), fixed = TRUE)
    colnames(x) <- paste0("G_", colnames(x))
    return(x)
}

all_data <- lapply(all_data, change_column_names)



rm(all_biol_example_data)
rm(breast_schill, colon_schill, renal_schill, gliob_schill, gliob_schill_15)

## </here>


## Recall we keep all data sets to a max number of columns of 15
## x <- pre_process(x, remove.constant = FALSE,
##                      min.freq = 0, max.cols = 15)


limit_to_15 <- function(x) {
    x <- df_2_mat_integer(x)
    x2 <- pre_process(x, remove.constant = FALSE,
                      min.freq = 0, max.cols = 15)
    return(list(Genotypes = sampledGenotypes(x2),
                data = x2,
                original_dims = dim(x),
                reduced_dims = dim(x2)))
}


lapply(all_data, dim)



genots_data <- lapply(all_data, limit_to_15)
save(file = "genots_data.RData", genots_data, compress = FALSE)

genots <- lapply(genots_data, function(x) x$Genotypes)
save(file = "genots.RData", genots, compress = FALSE)



date()
