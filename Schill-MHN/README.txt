About the files under the Schill-MHN directory
==============================================

1. For the files of this directory to be useful, you should first obtain
   the original code from Schill et al., available from
   https://github.com/RudiSchill/MHN.

   We have not added it to this repository since there is no license file
   in that repo.

   The version we used was downloaded in January 2020; as of 2021-05-12,
   that repository was last updated on 2018-08-16.

   So clone that repository and leave it as MHN subirectory under this
   directory. (e.g, ~/what_genotype_next/Schill-MHN/MHN) as several of the
   files below will try to source the code from ./MHN. Alternatively,
   modify the files and the path to your liking.


2. The main files are schill-trans-mat.R and run-Schill-MHN-trans-mat.R.

   schill-trans-mat.R contains the code to obtain the transition matrices
   between genotypes from the rate matrices returned by MHN. Several
   implementations are provided. The defaults are used in "do_MHN" and
   "do_MHN2" that only differ in the second using spare matrices.

   This file is sourced from
   ../run-biol-examples/code-all-methods-2-trans-matrix-max-genes-15.R
   (which is sources from
   ../run-biol-examples/run-all-biol-examples-max-genes-15.R, the script
   that is called to run all the biological examples).


   run-Schill-MHN-trans-mat.R sources the former and contains a couple of
   functions that will process the simulated files. This is the script
   that actually processes the simulated data.


3. Additional files:

   "pre-process.R" contains utility functions for merging identical
   columns. This is useful mainly in the analysis of the simulated
   data. This file is the same as the pre-process file in
   ../run-biol-examples

   "local-maxima.R" contains a few examples to illustrate the use of the
   TD approaches.


