
# Code and data for "Conditional prediction of consecutive tumor evolution using cancer progression models: What genotype comes next?"

## Juan Diaz-Colunga, Ramon Diaz-Uriarte	
## Published in PLoS Computational Biology, https://doi.org/10.1371/journal.pcbi.1009055

----------------------------------------------------------------------------------


## Main scripts and directories

Here we provide an overview of the main scripts and directories. Other
files and directories will be mentioned, as needed, below. Where deemed of
used for readers, we also provide .Rout files. RData and rds files are
provided when their size permits.

Executing the following scripts in the indicated order and with default
options creates a ./data directory within the current one where transition
matrices extracted from both the OncoSimulR simulations and the CPMs are
stored. It also generates a table.RData object containing comparisons
across said matrices.

**oncoFunctions.R** - Contains functions called by many other scripts and some tests.

**next-genotype_transitionMatrix-sim.R** - Extracts transition matrices directly from simulations. By default, matrices are stored in the ./sim directory.

**next-genotype_transitionMatrix-cpm.R** - Extracts transition matrices from the outputs of the CPMs. By default, matrices are stored in the ./cpm directory.

**next-genotype_transitionMatrix-timeDiscretizedCPMs.R** - Works just like the previos script, but for Schill's MHN and the time-discretized versions of CBN_ot and MCCBN. Output matrices are also stored under the ./cpm directory by default.

**structData.R** - Takes the matrices under the ./sim and ./cpm directories (by default, but different directories can be provided) and groups them into nested lists stored in .RData objects and saved in the ./data directory (by default).

**makeTable.R** - Reads the content of the ./data directory (by default) and creates the table.rds file containing a data frame with the comparative statistics across transition matrices.

**oncoPlots.R** - Sourcing this file provides a set of functions that can
be used to produce plots taking the table generated by makeTable as
input. Not used.

**merge-additional-info.R** - Run after
makeTable.R. Merges additional info and outputs table-replicates.rds
(approx. 20 GB file).

**average-and-array-statistics.R** - Run after
merge-additional-info.R. Computes per-array
statistics (with different weighting schemes) and averages over
replicates, genotypes, and number of mutations, appropriately. 

**makeTableObservedFreqs.R** - Can be run at any time, but necessarily
before merge-additional-info.R. Returns observed frequencies of genotypes
in actual samples used. 

**makeTableRanksFitness.R** - Can be run at any time, but necessarily
before merge-additional-info.R. Returns fitness ranks of genotypes in
fitness landscapes.

**makeTableLocalMax.R** - Can be run at any time, but necessarily
before merge-additional-info.R. Returns frequency with which each genotype
is a local maximum in simulations.

**makeTableWeights.R** - Can be run at any time, but necessarily
before merge-additional-info.R. Returns the true frequency with which a
genotype is the most abundant genotype during the simulations, as well as
the frequency of each genotype in the complete 20000 samples. These two
are the two sets of weights used.

**makeTable-to-diff-wrt-null-2.R** - Computes difference in JS between the
null model and the method output. 

**paired-wide-data.R** - Created wide (or paired-wide) data frames, so finding
out best method or performing paired tests is easier.

**create_rank_data.R** - Rank methods and obtain the best method in each
case. Requires wide data. Comparisons made at the number mutations, array,
and single genotype (per replicate) case.

**create_wide_g_long_g.R** - Create minimal, concise and with factors,
data sets for the glmertree models.

**data_for_weighted_glmertree_plots.R** - (under plots-glms) Create data
set to produce plots from lmertree fits. This produces
""data_for_weighted_glmertree_plots.RData". 

**merge-compsfixNulls-wide_g.R** and **explore-reduce-method-comp.R**:
comparisons among methods: merging with wide\_g data for lmertree
fits. For some tables below. (E.g., in file when_good_conditions_0677.R)

**method-comp-biol-data.R** Comparison of predictions between methods on
the biological data, under /run-biol-examples.


**./plot-glms** Directory with code for fitting the lmertrees and plotting
the results. Additionally, also scripts for producing tables in
Supplementary material of performance when similar CE and TD
(when_good_conditions_0677.R)

**./run-biol-examples** Code and data to run the CPM analyses and obtain the
transition matrices on the cancer data sets.

**./biol-data-plots** Plots in manuscript and supplementary material for the
output from the analyses on the cancer data sets.

**./Schill-MHN** Code from Schill et al., for MHN plus additional code for
us to obtain transition matrices. Since there is no license in the
original repository (https://github.com/RudiSchill/MHN) we have not
provided it. The code should be left under Shcill-MHN/MHN (without
creating a new subdirectory). A README.txt further explains the files and
functions under this directory.


**CBN-with-bug-fix-and-likelihood** A directory that contains the
compressed directory of H-CBN with a couple of bug fixes and enhancements
by RDU, as explained in the Supplementary Material. The original code from
https://bsse.ethz.ch/cbg/software/ct-cbn.html is under the GNU GPL (v. 3),
the same used here. The authors of the code are Niko
Beerenwinkel, Moritz Gerstung, and Seth Sullivant.

**Fitnesslandscape-characteristics** Code for plots of fitness landscapes
characteristics. One of them is used in the Supplementary Material
(gamma-rsign-obs-peaks-FL.pdf)


## Pipeline

00. The time discretized CBN and MCCBN transition matrices as well as MHNs
    are needed for analyses below. (Before the next-genotype-... are
    run). Run Time-discretized-CBN-MCCBN/time-discretized-CBN.R and
	Schill-MHN/run-Schill-MHN-trans-mat.R

0. At any time, but before running merge-additional-info.R, run each of
   makeTableObservedFreqs.R, makeTableRanksFitness.R, makeTableLocalMax.R,
   and makeTableWeights.R.
    
	Output:
    
    - makeTableObservedFreqs.R: outObsFreqs.RData (~ 95 MB, but fast script)
    
    - makeTableRanksFitness.R: outRanksFitness.RData. (~ 43 MB, but
             fast script)
    
    - makeTableWeights.R: weightsAll.RData [plus intermediate RData:
           outfrerqsSampled and outFrerqsTrue]. (~ 17 MB; takes about 6 h)
    
	- makeTableLocalMax.R: outLocalMax.RData (~ 560 KB; takes about 7 h)


1. Run all three scripts next-genotype_transitionMatrix-*

2. Run structData

3. Run makeTable.R. The -scratch flag is needed the first time this script
   runs, then once the temporary file makeTable_tmp.rds file has been
   created, formatting changes can be made and implemented by running
   makeTable with no flags. If any upstream changes are made, a new
   -scratch run will be needed to incorporate them.

4. Run merge-additional-file-info.R. Input:
   pre-table-replicates.RData (plus fl_sampl.RData and the output from
   each of the files run in step 0). Output: table-replicates.rds.
   
5. Run average-and-array-statistics.R. Input:
   table-replicates.rds. Output: FIXME: zz, spell out.

6. At this point, the generated tables can be used for downstream analyses and plots.

7. Run makeTable-to-diff-wrt-null-2.R: Input:
   data_no_any_over_replicate.rds, array_statistics_over_replicate.rds,
   data_no_any.rds, num_mut_statistics_over_replicate.rds. Output:
   genotype_diff_null.RData, array_diff_null.RData,
   genotype_diff_null_no_average.RData, num_mut_diff_null.RData.
   
8. Run paired-data.R. Input: genotype_diff_null.RData,
   array_diff_null.RData, genotype_diff_null_no_average.RData,
   num_mut_diff_null.RData. Output: wide_array.RData, wide_genotype.RData,
   wide_genotype_no_average.RData, wide_num_mut.RData.
   
9. Run create_rank_data.R, then create wide_g_long_g.R.  Note that we
    eliminate cases with sampledProp == 0 (create_wide_g_long_g.R), then
    those with sampledProp < 1e-3 (fit-lmtree scripts and
    modelse-beta-regression.R), as well as the last genotype ---7
    mutations or 10 mutations, for landscapes with 7 and 10 loci,
    respectively--- (fit-lmtree scripts and modelse-beta-regression.R).

10. Run the lmertree fits.
    In subdirectory plots-glms/fits-Weighted-01-sqrt. A launch.sh will launch the
    lmertree fits. Once done, plotting with plotting-lmertree-fits-weighted-01-sqrt.R


11. Analyze the cancer data sets. Under ./run-biol-examples 

    - Launch analyses with CPMs. See details under README.txt. File
      launch-run-rerun.sh launches all R processes.
	  
	- Once analyses with CPMs are completed, compare output:
      method-comp-biol-data.R
	  

12. Run the two files under ./biol-data-plots to generate the plots for
the cancer data sets (main manuscript and supplementary material).


## Contact

For questions or comments, please contact 
[J. Diaz-Colunga](mailto:juan.diazcolunga@yale.edu) or [R. Diaz-Uriarte](mailto:r.diaz@uam.es).




