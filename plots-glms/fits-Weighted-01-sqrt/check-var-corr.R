ff <- dir(pattern = glob2rx("*.RData"))

load(ff[1])

VarCorr(lmertree_sqrt_W_Obs_minsize_01_fl_specific_10genes_gamma_epist_withFL_Inf_NM_allFL_Inf)
##  Groups            Name        Std.Dev.
##  id:sourceGenotype (Intercept) 0.148630
##  id:replicate      (Intercept) 0.012193
##  id                (Intercept) 0.066412
## Residual                      0.026711
rm(list = ls(pattern = "lmertree*"))

load(ff[2])

VarCorr(lmertree_sqrt_W_Obs_minsize_01_fl_specific_7genes_gamma_epist_withFL_Inf_NM_allFL_Inf)
 ## Groups            Name        Std.Dev.
 ## id:sourceGenotype (Intercept) 0.164961
 ## id:replicate      (Intercept) 0.022325
 ## id                (Intercept) 0.048077
 ## Residual                      0.028442

rm(list = ls(pattern = "lmertree*"))

load(ff[3])

VarCorr(get("lmertree_sqrt_W_Obs_minsize_01_fLandscape-rsign-peaks-recoded_10genes_Inf_NM_allFL_Inf"))
 ## Groups            Name        Std.Dev.
 ## id:sourceGenotype (Intercept) 0.133893
 ## id:replicate      (Intercept) 0.027864
 ## id                (Intercept) 0.076799
 ## Residual                      0.138810

load(ff[4])

VarCorr(get("lmertree_sqrt_W_Obs_minsize_01_fLandscape-rsign-peaks-recoded_7genes_Inf_NM_allFL_Inf"))

 ## Groups            Name        Std.Dev.
 ## id:sourceGenotype (Intercept) 0.163042
 ## id:replicate      (Intercept) 0.026274
 ## id                (Intercept) 0.081556
 ## Residual                      0.126074
