## Find out best JS and best method and rank methods
date()
rm(list = ls())
library(data.table)
setDTthreads(threads = 0)

## wide_* data files created in paired-wide-data.R
load("wide_genotype.RData")
load("wide_genotype_no_average.RData")
load("wide_array.RData")
load("wide_num_mut.RData")

colnames_genot <- grep("_js$", colnames(wide_genotype), value = TRUE)
colnames_genot2 <- grep("_js$", colnames(wide_genotype_no_average),
                        value = TRUE)
stopifnot(identical(colnames_genot, colnames_genot2))

colnames_arr <- grep("_js_w_sampl$",
                     colnames(wide_array), value = TRUE)

colnames_num_mut <- grep("_js_w_sampl$",
                         colnames(wide_num_mut), value = TRUE)

stopifnot(identical(colnames_arr, colnames_num_mut))

## Return the minimum, which is the minimum (best)
## and rank. Designed for data tables.
my_min_which_rank <- function(x) {
    ## if(is.data.table(x)) {
    ##     mnames <- colnames(x)
    ## } else if(is.vector(x)) {
    ##     mnames <- names(x)
    ## }
    mnames <- colnames(x)
    ## deal with NA and NaN
    ## if((length(na.omit(x)) == 0) 
    ##    ||
    ##    (nrow(na.omit(x)) == 0 )){ ## wrong: this is for apply-like or by row ops
    if(all(is.na(x)) ) {
        mx <- NA_real_
        ww <- NA_integer_
        single_best_method <- NA_character_
        number_best <- NA_integer_
        all_best_methods <- NA_character_
    } else {
        mx <- min(x, na.rm = TRUE)
        ww <- which(x == mx)
        stopifnot(length(ww) > 0)
        if(length(ww) == 1) {
            single_best_method <- mnames[ww]
            number_best <- 1L
            all_best_methods <- single_best_method
        } else {
            all_best_methods <- paste(mnames[ww], sep = ", ",
                                      collapse = ", ")
            number_best <- length(ww)
            ## Prevent biases from alphabetic ordering
            ww <- sample(ww, 1)
            single_best_method <- NA_character_
        }
    }

    rx <- rank(x, na.last = "keep")
    names(rx) <- paste0("rank_", mnames)

    return(
        c(list(min_js = mx,
               which_single_best_js = single_best_method,
               number_best_js = number_best,
               all_best_js = all_best_methods),
          do.call(list, as.list(rx))
          )
    )
}


## Create data.table with the ranks
rank_and_merge <- function(dt, cnames) {
    ## do this outside
    ## dt[, rownum := .I]
    out1 <- dt[, my_min_which_rank(.SD),
               .SDcols = cnames,
               by = rownum]
    nr1 <- nrow(dt)
    nr2 <- nrow(out1)
    finalout <- merge(dt, out1, by = "rownum")
    nr3 <- nrow(finalout)
    stopifnot(identical(nr1, nr2))
    stopifnot(identical(nr1, nr3))
    return(finalout)
}

date()

wide_num_mut[, rownum := .I]
wide_array[, rownum := .I]
wide_genotype[, rownum := .I]
wide_genotype_no_average[, rownum := .I]

wide_num_mut_rank <- rank_and_merge(wide_num_mut, colnames_num_mut)
date()

wide_array_rank <- rank_and_merge(wide_array, colnames_arr)
date()

wide_genotype_rank <- rank_and_merge(wide_genotype, colnames_genot)
date()
wide_genotype_no_average_rank <- rank_and_merge(wide_genotype_no_average,
                                                colnames_genot)
date()

save(file = "wide_array_rank.RData", wide_array_rank, compress = FALSE)
save(file = "wide_genotype_rank.RData", wide_genotype_rank, compress = FALSE)
save(file = "wide_genotype_no_average_rank.RData",
     wide_genotype_no_average_rank, compress = FALSE)
save(file = "wide_num_mut_rank.RData", wide_num_mut_rank, compress = FALSE)

gc()
date()

