- We run all the methods on the following data sets:

  - The RData is the same as file all_biol_example_data.RData, in
    Diaz-Uriarte & Vasallo, 2019.

  - The remaing data sets used are from Schill et al. Three of them were
    previously used by Gerstung et al., 2009. Read from the MHN directory
    (which you can obtain from the repo for the suppl mat of Schill et
    al., MHN paper).



- We limit all analyses to the 15 columns with largest frequency (some
  data sets have hundreds of columns) as the attempts to run data sets
  with more than 16 columns failed.

  - Launched from launch-run-rerun.sh which calls
    run-all-biol-examples-max-genes-15.R, that limits runs to the first 15
    genes. (But see below for Schill's gliobl.)

    The code in run-all-biol-examples-max-genes-15.R is basically the same
  as in run-all-biol-examples-cores-cbn.R, except for adding "max genes
  genes 15" in the cat's and _max_genes_15 in the RData. And they also
  differ in that -15.R renamed genes (see below).

  Both differ in which file they source (code-all-methods-2-trans-matrix.R
  and code-all-methods-2-trans-matrix-max-genes-15.R). The code in
  code-all-methods-2-trans-matrix.R and
  code-all-methods-2-trans-matrix-max-genes-15.R is basically the same
  except for using a max number of 20 or 15 genes.

  So datasets 20 to 27 are all trimmed to have 15 columns.


  - We had to rename gene names because of CAPRESE and CAPRI. This is done
  in run-all-biol-examples-max-genes-15.R: see function
  "change_column_names". When plotting, we of course removed the prefix
  and postfix to the genes.
  


- For the Glioblastoma data, we have two runs:
  - one with a data set called gliob_schill_15: this uses only 15 genes,
    the same 15 that are shown in Figure 6 of their paper.

  - another one with 15 genes, but with trimming based on frequency
    (gliob_schill); thus, this is redundant.

  - The Glioblastoma data was not used in the paper as it overlaps with
    other data.

- Output:
  - We use sparse matrices. Doing "as.matrix(output)" gives you the dense
   matrix (which could be a huge object, more than can be fit in RAM)

  - Not all methods make all genotypes accessible. Thus, some transition
    matrices are comparatively small (e.g., ACML_co: OT give a 325 by 325,
    CBN 85 by 85, and MHN the full 2048 x 2048)

- This directory: Note changes to in cbn-process.R (w.r.t. the one used in
  Diaz-Uriarte and Vasallo), in fuction call.external.cbn, so that we can
  set OMP cores.






