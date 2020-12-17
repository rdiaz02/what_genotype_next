##

date()

arguments.script <- commandArgs(trailingOnly = TRUE)
cores <- as.numeric(arguments.script[1])
dataset <- as.numeric(arguments.script[2])


source("code-all-methods-2-trans-matrix-max-genes-15.R")



RhpcBLASctl::blas_set_num_threads(cores)
RhpcBLASctl::omp_set_num_threads(cores)



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


## This is what we set in here. Recall CBN will also set it.
cat("\n What R thinks about number of threads/processes\n")

RhpcBLASctl::get_num_procs()
RhpcBLASctl::get_num_cores()
RhpcBLASctl::blas_get_num_procs()
RhpcBLASctl::omp_get_num_procs()
RhpcBLASctl::omp_get_max_threads()


this_data <- all_data[[dataset]]
dataname <- names(all_data)[dataset]
outname <- paste0("transition_matrices_num", dataset, "_dataname_",
                  dataname, "_cores_", cores, "_cores_cbn_version",
                  "_max_genes_15")
filename <- paste0(outname, ".RData")



cat("\n##############################################\n")
cat("\n Doing dataset  ", dataset, " named ", dataname, " with cores ", cores,
    " with max genes 15", "\n")


total_time <- system.time(
    assign(outname, all_methods_2_trans_mat(this_data, cores_cbn = cores))
)["elapsed"]


save(list = outname, file = filename)

cat("\n Total time = ", total_time)

cat("\n Done \n")


date()
