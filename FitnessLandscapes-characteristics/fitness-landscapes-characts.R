## Simplified and modified from
## /home/ramon/Proyectos/JDC-oncolab/FitnessLandscape-characteristics/data-plots-landscape-features.R
## used in Diaz-Uriarte and Vasallo, 2019.


## See that file for further details about variables


library(ggplot2)
library(latex2exp)
library(scales)

load("df_out_clonal_interf_stats_for_plot_anal.RData")
load("s_df.RData")




## "Local_peaks" as "Local maxima"
##   ltype_rename in plots-coefs.R takes a df and returns a df.
##   this here is more cumbersome, but safer.
ltype_rename <- function(x) {
    oldv <- x
    newv <- as.character(x)
    newv[newv == "Local_peaks"] <- "Local maxima"
    newv <- factor(newv,
                   levels = c("Represent.",
                              "Local maxima",
                              "RMF"))
    ## Paranoid checks
    tt <- table(oldv, newv)
    dtt <- diag(tt)
    stopifnot(sum(dtt) == length(x))
    stopifnot(colnames(tt) == c("Represent.", "Local maxima", "RMF"))
    stopifnot(rownames(tt) == c("Represent.", "Local_peaks", "RMF"))
    ## x$landscapeType <- newv
    ## x$landscapeType_before_rename <- oldv
    ## Nope, return only the variable: easier to catch assigning 
    ## to the wrong object
    return(newv)
}


df_out$type_Landscape <- ltype_rename(df_out$type_Landscape)
## FIXME: is this needed? sf?
s_df$type_Landscape <- ltype_rename(s_df$type_Landscape)




## Done inside local to modify a column name for nicer rendering
## and consitent naming with other figures
local({
    colnames(df_out)[143]
    colnames(df_out)[143] <- "typeLandscape"
    df_out$numObservedPeaks <- cut_interval(df_out$num_observed_peaks, 2)
    ## Understanding some splits because of the perfect separation of RMF
## and almost perfect of local max
pdf(file = "gamma-rsign-obs-peaks-FL.pdf", height = 5.8, width = 8.1)
ggplot(data = df_out,
      aes(color = typeLandscape,
          shape = typeLandscape,
          x = epist_rsign,
          y = gamma,
          )) +
    ## scale_y_continuous(trans = log1p_trans()) +
    scale_x_continuous(trans = log1p_trans()) +
    facet_wrap(. ~ numObservedPeaks, scales = "free", labeller = label_both) +
    geom_point(size = 1) + ##0.5) +
    labs(y = TeX("$\\gamma$ (gamma)"), 
         x = "Fraction reciprocal sign epistasis"
         ) ## +
    ## ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))
dev.off()
})


## Continue if you want a lot more plots of fitness landscape charact
stop()

pdf(file = "fitness-landscapes-characts.pdf", height = 5.8, width = 8.1)
## pdf(file = "fitness-landscapes-characts%03d.pdf",
##     onefile = FALSE,
##     height = 5.8, width = 8.1)

## Fig. D in S2 Fig of Diaz-Uriarte and Vasallo.
ggplot(data = df_out, 
       aes(color = Init_Size,
           y = freq_most_freq_mean_no450,
           x = Ngenes
           )) +    
    facet_grid(Mutation ~ type_Landscape) +
    geom_boxplot() +
    labs(x = "Number of genes",
         y = "Average frequency of most frequent genotype"
         ) +
    ggtitle("Clonal interference") +
    theme(plot.title = element_text(hjust = 0.5))


## Fig. E in S2 Fig of Diaz-Uriarte and Vasallo.
## Both variables are measuring the same phenomenon 
ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation, ## Ngenes,
          x = freq_most_freq_mean_no450,
          y = how_many_gt_5p_mean_no450
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 0.8) +
    labs(x = "Average frequency of most frequent genotype",
         y = "Average number of clones with frequency > 5%"
         ) +
    ggtitle("Clonal interference") +
    theme(plot.title = element_text(hjust = 0.5))

ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = pom_h,
          x = lod_h,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    labs(y = "POM diversity",
         x = "LOD diversity"
         ) +
    ggtitle("LOD and POM diversity") +
    theme(plot.title = element_text(hjust = 0.5))


ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = lod_h,
          x = how_many_gt_5p_mean_no450,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    labs(y = "LOD diversity",
         x = "Average number of clones with frequency > 5%"
         ) +
    ggtitle("Clonal interference and POM diversity") +
    theme(plot.title = element_text(hjust = 0.5))

ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = lod_h,
          x = freq_most_freq_mean_no450,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    labs(y = "LOD diversity",
         x = "Average frequency of most frequent genotype"
         ) +
    ggtitle("Clonal interference and POM diversity") +
    theme(plot.title = element_text(hjust = 0.5))

ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = pom_h,
          x = how_many_gt_5p_mean_no450,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    labs(y = "POM diversity",
         x = "Average number of clones with frequency > 5%"
         ) +
    ggtitle("Clonal interference and POM diversity") +
    theme(plot.title = element_text(hjust = 0.5))

## ## The above is clearer
## ggplot(data = df_out,
##       aes(color = Init_Size,
##           shape = Ngenes,
##           x = pom_h,
##           y = freq_most_freq_mean_no450,
##           )) +    
##     facet_grid(Mutation ~ type_Landscape) +
##     geom_point(size = 0.5) +
##     labs(x = "POM diversity",
##          y = "Average frequency of most frequent genotype"
##          ) +
##     ggtitle("Clonal interference and evolutionary predictability") +
##     theme(plot.title = element_text(hjust = 0.5))



## About the variables n_peaks, num_observed_peaks and num_local_peaks
## How they are obtained can be seen in sample-simuls-A-B.R, these lines:

##     df1$num_observed_peaks <- out0$number_last_lod
##     df1$num_accessible_genots <- length(out0$ruggified_dag$accessible_genots)
##     df1$num_local_peaks <- out0$ruggified_dag$num_peaks_noback
##     ## MAGELLAN stats
##     df1 <- cbind(df1, out0$ruggified_dag$magellan_stats)


## n_peaks comes directly from Magellan (see function Magellan_stats). It
## is not very useful for us, as it gives number of peaks in the
## landscape, regardless of whether or not they are accessible, and thus
## not account for the no-back-mutation assumption.

## num_local_peaks is computed by us under the no-back-mutation
## assumption. This is the variable peaks_noback in the
## ruggify-functions.R file, the "DAG_to_rugged" function. And actually
## calls and OncoSimulR function, "fast_peaks". That is an improvement,
## but ignores completely the dynamics: some local peaks in the fitness
## landscape might actually never be visited because of mutation rates or
## their fitness.

## num_observed_peaks: that gives the number of local peaks in the
## landscape that are actually visited in the evolutionary simulations.




## num_local_peaks, _num_observed_peaks, diversity_observed_peaks


ggplot(data = dplyr::filter(df_out, type_Landscape != "Represent."),
       aes(color = Init_Size,
           shape = Mutation,
           y = num_observed_peaks,
           x = num_local_peaks
           )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1, position = "jitter") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Number of local fitness maxima",
         y = "Number of observed local fitness maxima") +
    ggtitle("Static fitness landscape characteristics") +
    theme(plot.title = element_text(hjust = 0.5))


ggplot(data = dplyr::filter(df_out, type_Landscape != "Represent."),
       aes(color = Init_Size,
           shape = Mutation,
           y = diversity_observed_peaks,
           x = num_local_peaks
           )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Number of local fitness maxima",
         y = "Diversity of observed local fitness maxima") +
    ggtitle("Static fitness landscape characteristics") +
    theme(plot.title = element_text(hjust = 0.5))


## Fig. C in S2 Fig of Diaz-Uriarte and Vasallo.
ggplot(data = dplyr::filter(df_out, type_Landscape != "Represent."),
       aes(color = Init_Size,
           shape = Mutation,
           y = epist_rsign,
           x = num_local_peaks
           )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) +
    scale_x_log10() +
    scale_y_sqrt(breaks = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4)) +
    ## scale_y_log10() +
    ## scale_y_continuous(trans = log1p_trans()) +
    ## coord_trans(y = expm1_trans()) +
    labs(x = "Number of local fitness maxima",
         y = "Fraction reciprocal sign epistasis") +
    ggtitle("Static fitness landscape characteristics") +
    theme(plot.title = element_text(hjust = 0.5))


ggplot(data = dplyr::filter(df_out, type_Landscape != "Represent."),
       aes(color = Init_Size,
           shape = Mutation,
           y = epist_rsign,
           x = num_local_peaks
           )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) +
    scale_x_log10() +
    scale_y_sqrt(breaks = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4)) +
    ## scale_y_log10() +
    ## scale_y_continuous(trans = log1p_trans()) +
    ## coord_trans(y = expm1_trans()) +
    labs(x = "Number of local fitness maxima",
         y = "Fraction reciprocal sign epistasis") +
    ggtitle("Static fitness landscape characteristics") +
    theme(plot.title = element_text(hjust = 0.5))





## Nothing relevant
ggplot(data = dplyr::filter(df_out, type_Landscape != "Represent."),
      aes(color = Init_Size,
          shape = Mutation,
          y = pom_h,
          x = num_local_peaks,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
        scale_x_log10() +
    labs(y = "POM diversity",
         x = "Number of local fitness maxima"
         ) +
    ggtitle("Local peaks and POM diversity") +
    theme(plot.title = element_text(hjust = 0.5))


## Just a minor positive one
ggplot(data = dplyr::filter(df_out, type_Landscape != "Represent."),
      aes(color = Init_Size,
          shape = Mutation,
          y = pom_h,
          x = num_observed_peaks,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
        scale_x_log10() +
    labs(y = "POM diversity",
         x = "Number of observed local fitness maxima"
         ) +
    ggtitle("Observed local peaks and POM diversity") +
    theme(plot.title = element_text(hjust = 0.5))

ggplot(data = dplyr::filter(df_out, type_Landscape != "Represent."),
      aes(color = Init_Size,
          shape = Mutation,
          y = lod_h,
          x = num_observed_peaks,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
        scale_x_log10() +
    labs(y = "LOD diversity",
         x = "Number of observed local fitness maxima"
         ) +
    ggtitle("Observed local peaks and POM diversity") +
    theme(plot.title = element_text(hjust = 0.5))


## Much clearer
ggplot(data = dplyr::filter(df_out, type_Landscape != "Represent."),
      aes(color = Init_Size,
          shape = Mutation,
          y = pom_h,
          x = diversity_observed_peaks,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    labs(y = "POM diversity",
         x = "Diversity of observed local fitness maxima"
         ) +
    ggtitle("Diversity of observed local peaks and POM diversity") +
    theme(plot.title = element_text(hjust = 0.5))


ggplot(data = dplyr::filter(df_out, type_Landscape != "Represent."),
       ## data = df_out, 
       aes(color = Init_Size,
           y = diversity_observed_peaks,
           x = Ngenes
           )) +    
    facet_grid(Mutation ~ type_Landscape) +
    geom_boxplot() +
    labs(x = "Number of genes",
         y = "Diversity of observed local fitness maxima"
         ) +
    ggtitle("Diversity of observed local fitness maxima") +
    theme(plot.title = element_text(hjust = 0.5))


## Nothing remarkable
ggplot(data = dplyr::filter(df_out, type_Landscape != "Represent."),##df_out,
      aes(color = Init_Size,
          shape = Mutation,
          x = num_observed_peaks,
          y = how_many_gt_5p_mean_no450,
          )) +    
    facet_grid(Mutation ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    labs(x = "Number of observed local fitness maxima",
         y = "Average number of clones with frequency > 5%"
         ) +
    ggtitle("Clonal interference and observed local peaks") +
    theme(plot.title = element_text(hjust = 0.5))


ggplot(data = dplyr::filter(df_out, type_Landscape != "Represent."),##df_out,
      aes(color = Init_Size,
          shape = Mutation,
          x = num_local_peaks,
          y = how_many_gt_5p_mean_no450,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    labs(x = "Number of local fitness maxima",
         y = "Average number of clones with frequency > 5%"
         ) +
    ggtitle("Clonal interference and local peaks") +
    theme(plot.title = element_text(hjust = 0.5))

ggplot(data = dplyr::filter(df_out, type_Landscape != "Represent."),##df_out,
      aes(color = Init_Size,
          shape = Mutation,
          x = num_observed_peaks,
          y = freq_most_freq_mean_no450,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    labs(x = "Number of observed local fitness maxima",
         y = "Average frequency of most frequent genotype"
         ) +
    ggtitle("Clonal interference and observed local peaks") +
    theme(plot.title = element_text(hjust = 0.5))




## Weak but present, and stronger, if less steep, in size 2000
ggplot(data = dplyr::filter(df_out, type_Landscape != "Represent."),##df_out,
      aes(color = Init_Size,
          shape = Mutation,
          x = diversity_observed_peaks,
          y = how_many_gt_5p_mean_no450,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    labs(x = "Diversity of observed local fitness maxima",
         y = "Average number of clones with frequency > 5%"
         ) +
    ggtitle("Clonal interference and diversity observed local peaks") +
    theme(plot.title = element_text(hjust = 0.5))



ggplot( ## data = dplyr::filter(df_out, type_Landscape != "Represent."),
       data = df_out, 
       aes(color = Init_Size,
           y = gamma,
           x = Ngenes
           )) +    
    facet_grid(Mutation ~ type_Landscape) +
    geom_boxplot() +
    labs(x = "Number of genes",
         y = TeX("$\\gamma$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))



ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = gamma,
          x = epist_rsign,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    labs(x = "Fraction reciprocal sign epistasis",
         y = TeX("$\\gamma$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))

ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = gamma,
          x =  num_local_peaks,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    scale_x_log10() +
    labs(x = "Number of local fitness maxima",
         y = TeX("$\\gamma$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))


ggplot(data = dplyr::filter(df_out, type_Landscape != "Represent."),
       aes(color = Init_Size,
           shape = Mutation,
           y = gamma,
           x = num_observed_peaks
           )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) +
    scale_x_log10() +
##    scale_y_sqrt(breaks = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4)) +
    ## scale_y_log10() +
    ## scale_y_continuous(trans = log1p_trans()) +
    ## coord_trans(y = expm1_trans()) +
    labs(x = "Number of observed local fitness maxima",
         y = TeX("$\\gamma")) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))



ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = gamma,
          x = num_observed_peaks,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    labs(x = "Number of observed local fitness maxima",
         y = TeX("$\\gamma$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))


ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Ngenes,
          x = gamma,
          y = freq_most_freq_mean_no450,
          )) +    
    facet_grid(Mutation ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    labs(y = "Average frequency of most frequent genotype",
         x = TeX("$\\gamma$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))

ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          x = gamma,
          y = lod_h,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    labs(y = "LOD diversity",
         x = TeX("$\\gamma$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))


## yes, it makes sense
summary(rowSums(df_out[, c("w.1.", "w.2.", "w.3..")]))

## CBN (and OT) can encode epistasis of higher order

ggplot( ## data = dplyr::filter(df_out, type_Landscape != "Represent."),
    ## data = dplyr::filter(df_out, tree == TRUE),
    data = df_out,
       aes(color = Init_Size,
           y = w.1.,
           x = Ngenes
           )) +    
    facet_grid(Mutation ~ type_Landscape) +
    geom_boxplot() +
    labs(x = "Number of genes",
         y = TeX("Fourier expansion: $W_1$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))


ggplot( ## data = dplyr::filter(df_out, type_Landscape != "Represent."),
    ##data = dplyr::filter(df_out, tree == TRUE),
    data = df_out,
       aes(color = Init_Size,
           y = w.2.,
           x = Ngenes
           )) +    
    facet_grid(Mutation ~ type_Landscape) +
    geom_boxplot() +
    labs(x = "Number of genes",
         y = TeX("Fourier expansion: $W_2$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))


ggplot( ## data = dplyr::filter(df_out, type_Landscape != "Represent."),
    ##data = dplyr::filter(df_out, tree == TRUE),
    data= df_out,
       aes(color = Init_Size,
           y = w.3..,
           x = Ngenes
           )) +    
    facet_grid(Mutation ~ type_Landscape) +
    geom_boxplot() +
    labs(x = "Number of genes",
         y = TeX("Fourier expansion: $W_{3+}$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))



ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = w.1.,
          x = epist_rsign,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    labs(x = "Fraction reciprocal sign epistasis",
         y = TeX("Fourier expansion: $W_1$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))


ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = w.2.,
          x = epist_rsign,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    labs(x = "Fraction reciprocal sign epistasis",
         y = TeX("Fourier expansion: $W_2$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))


## Very interesting plot.
ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = w.3..,
          x = epist_rsign,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    ## scale_x_log10() +
    ## coord_trans(y = "log1p") + 
    labs(x = "Fraction reciprocal sign epistasis",
         y = TeX("Fourier expansion: $W_{3+}$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))

ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = w.3..,
          x = epist_rsign,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    ## scale_x_log10() +
    labs(x = "Fraction reciprocal sign epistasis",
         y = TeX("Fourier expansion: $W_{3+}$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))



ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = w.3..,
          x = w.2.,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    ## scale_x_log10() +
    labs(x = TeX("Fourier expansion: $W_{2}$"),
         y = TeX("Fourier expansion: $W_{3+}$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))



ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = w.1.,
          x = num_local_peaks,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    scale_x_log10() +
    labs(x = "Number of local fitness maxima",
         y = TeX("Fourier expansion: $W_1$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))

ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = w.2.,
          x = num_local_peaks,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    scale_x_log10() +
    labs(x = "Number of local fitness maxima",
         y = TeX("Fourier expansion: $W_2$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))

ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = w.3..,
          x = num_local_peaks,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    scale_x_log10() +
    labs(x = "Number of local fitness maxima",
         y = TeX("Fourier expansion: $W_3+$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))





ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = w.1.,
          x = num_observed_peaks,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    labs(x = "Number of observed local fitness maxima",
         y = TeX("Fourier expansion: $W_1$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))


ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = w.2.,
          x = num_observed_peaks,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    labs(x = "Number of observed local fitness maxima",
         y = TeX("Fourier expansion: $W_2$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))


ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = w.3..,
          x = num_observed_peaks,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    labs(x = "Number of observed local fitness maxima",
         y = TeX("Fourier expansion: $W_3+$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))




ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = gamma,
          x = w.1.,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    labs(y = TeX("$\\gamma"),
         x = TeX("Fourier expansion: $W_1$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))



ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = gamma,
          x = w.2.,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    labs(y = TeX("$\\gamma"),
         x = TeX("Fourier expansion: $W_2$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))



ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = gamma,
          x = w.3..,
          )) +    
    facet_grid(Ngenes ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    labs(y = TeX("$\\gamma"),
         x = TeX("Fourier expansion: $W_3+$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))






ggplot(data = s_df,
       aes(color = Detection, ## Sampling,
           y = sampledGenotypesNumber,
           x = Init_Size
           )) +    
    facet_grid(Ngenes + Mutation ~ type_Landscape) +
    geom_boxplot(lwd = 0.3, outlier.size = 0.2)  +
    labs(x = "Init Size",
         y = "Number of different sampled genotypes") +
    ggtitle("Samples' characteristics.") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = c("plum3", "lightgreen", "orange"))


ggplot(data = s_df,
       aes(color = Detection, ## Sampling,
           y = sampledGenotypesNumber,
           x = Init_Size
           )) +    
    facet_grid(Ngenes + Mutation ~ type_Landscape) +
    geom_boxplot(lwd = 0.3, outlier.size = 0.2)  +
    labs(x = "Init Size",
         y = "Number of different sampled genotypes") +
    ggtitle("Samples' characteristics.") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = c("plum3", "lightgreen", "orange"))

ggplot(data = s_df,
       aes(color = Detection, ## Sampling,
           y = Mean_muts,
           x = Init_Size
           )) +    
    facet_grid(Ngenes + Mutation ~ type_Landscape) +
    geom_boxplot(lwd = 0.3, outlier.size = 0.2)  +
    labs(x = "Init Size",
         y = "Mean number of mutations in sampled genotypes") +
    ggtitle("Samples' characteristics.") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = c("plum3", "lightgreen", "orange"))

ggplot(data = s_df,
       aes(color = Detection, ## Sampling,
           y = Skewness_muts,
           x = Init_Size
           )) +    
    facet_grid(Ngenes + Mutation ~ type_Landscape) +
    geom_boxplot(lwd = 0.3, outlier.size = 0.2)  +
    labs(x = "Init Size",
         y = "Skewness number of mutations in sampled genotypes") +
    ggtitle("Samples' characteristics.") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = c("plum3", "lightgreen", "orange"))

ggplot(data = s_df,
       aes(color = Detection, ## Sampling,
           y = Stdev_muts,
           x = Init_Size
           )) +    
    facet_grid(Ngenes + Mutation ~ type_Landscape) +
    geom_boxplot(lwd = 0.3, outlier.size = 0.2)  +
    labs(x = "Init Size",
         y = "Standard deviation number of mutations in sampled genotypes") +
    ggtitle("Samples' characteristics.") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = c("plum3", "lightgreen", "orange"))


ggplot(data = s_df,
       aes(color = Init_Size, ## Sampling,
           shape = Mutation,
           y = sampledGenotypesNumber,
           x = Mean_muts
           )) +    
    facet_grid(Ngenes + Detection ~ type_Landscape,
               scales = "free") +
    geom_point(size = 1)  +
    labs(x = "Mean number of mutations in sampled genotypes",
         y = "Number of different sampled genotypes") +
    ggtitle("Samples' characteristics.") +
    theme(plot.title = element_text(hjust = 0.5)) 


ggplot(data = s_df,
       aes(color = Init_Size, ## Sampling,
           shape = Mutation,
           y = Skewness_muts,
           x = Mean_muts
           )) +    
    facet_grid(Ngenes  ~ type_Landscape,
               scales = "free") +
    geom_point(size = 1)  +
    labs(x = "Mean number of mutations in sampled genotypes",
         y = "Skewness number of mutations in sampled genotypes") +
    ggtitle("Samples' characteristics.") +
    theme(plot.title = element_text(hjust = 0.5)) 


ggplot(data = s_df,
       aes(color = Init_Size, ## Sampling,
           shape = Mutation,
           y = Stdev_muts,
           x = Mean_muts
           )) +    
    facet_grid(Ngenes  ~ type_Landscape,
               scales = "free") +
    geom_point(size = 1)  +
    labs(x = "Mean number of mutations in sampled genotypes",
         y = "Standard deviation number of mutations in sampled genotypes") +
    ggtitle("Samples' characteristics.") +
    theme(plot.title = element_text(hjust = 0.5)) 

ggplot(data = s_df,
       aes(color = Detection,
           shape = Mutation,
           y = Stdev_muts,
           x = Mean_muts
           )) +    
    facet_grid(Ngenes + Detection ~ type_Landscape,
               scales = "free") +
    geom_point(size = 1)  +
    labs(x = "Mean number of mutations in sampled genotypes",
         y = "Standard deviation number of mutations in sampled genotypes") +
    ggtitle("Samples' characteristics.") +
    theme(plot.title = element_text(hjust = 0.5)) 

ggplot(data = s_df,
       aes(color = Init_Size, ## Sampling,
           shape = Mutation,
           x = sampledGenotypesNumber,
           y = Stdev_muts
           )) +    
    facet_grid(Ngenes + Detection  ~ type_Landscape,
               scales = "free") +
    geom_point(size = 1)  +
    labs(y = "Standard deviation number of mutations in sampled genotypes",
         x = "Number of different sampled genotypes") +
    ggtitle("Samples' characteristics.") +
    theme(plot.title = element_text(hjust = 0.5)) 

dev.off()




## For the plots for the glmertree ideas: models with fitness landscape,
## without fitness landscape, and correlated variables

## FIXME: add scatterplot of epistRSign and gamma, with different colors for the three
## landscapes, maybe split by something else

## Number of local fitness maxima: Clear difference by landscape
ggplot(data = df_out, 
       aes(color = Init_Size,
           y = num_observed_peaks,
           x = Ngenes
           )) +    
    facet_grid(Mutation ~ type_Landscape) +
    geom_boxplot() +
    labs(x = "Number of genes",
         y = "Number of observed local fitness maxima"
         ) +
    ggtitle("Fitness landscape and observed peaks") +
    theme(plot.title = element_text(hjust = 0.5))

## Reciprocal sign epistasis: Clear difference by landscape
ggplot(data = df_out, 
       aes(color = Init_Size,
           y = epist_rsign,
           x = Ngenes
           )) +    
    facet_grid(Mutation ~ type_Landscape) +
    geom_boxplot() +
    labs(x = "Number of genes",
         y = "Fraction reciprocal sign epistasis"
         ) +
    ggtitle("Fitness landscape and reciprocal sign epistasis peaks") +
    theme(plot.title = element_text(hjust = 0.5))


## Gamma: clear difference by landscape (much smaller in RMF)
ggplot(data = df_out, 
       aes(color = Init_Size,
           y = gamma,
           x = Ngenes
           )) +    
    facet_grid(Mutation ~ type_Landscape) +
    geom_boxplot() +
    labs(x = "Number of genes",
         y =  TeX("$\\gamma$")
         ) +
    ggtitle("Fitness landscape and gamma") +
    theme(plot.title = element_text(hjust = 0.5))



## SSWM: Landscape, per se, has little effect 
ggplot(data = df_out, 
       aes(color = Init_Size,
           y = freq_most_freq_mean_no450,
           x = Ngenes
           )) +    
    facet_grid(Mutation ~ type_Landscape) +
    geom_boxplot() +
    labs(x = "Number of genes",
         y = "Average frequency of most frequent genotype"
         ) +
    ggtitle("Fitness landscape and SSWM") +
    theme(plot.title = element_text(hjust = 0.5))

## SSWM: no relationship with gamma
ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Ngenes,
          x = gamma,
          y = freq_most_freq_mean_no450,
          )) +    
    facet_grid(Mutation ~ type_Landscape) +
    geom_point(size = 1) + ##0.5) +
    labs(y = "Average frequency of most frequent genotype",
         x = TeX("$\\gamma$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))

## SSWM: no relationship with reciprocal sign epist.
ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Ngenes,
          x = epist_rsign,
          y = freq_most_freq_mean_no450,
          )) +    
    facet_grid(Mutation ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    labs(y = "Average frequency of most frequent genotype",
         x = "Fraction reciprocal sign epistasis"
         ) +
    ## scale_x_log10() +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))



## gamma and reciprocal sign epist: related, and relation changes by FL
ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = gamma,
          x = epist_rsign,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    labs(x = "Fraction reciprocal sign epistasis",
         y = TeX("$\\gamma$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))

## gamma and number of observed peaks: related, and relation changes by FL
ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = gamma,
          x = num_observed_peaks,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    labs(x = "Number of observed local fitness maxima",
         y = TeX("$\\gamma$")
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))


## reciprocal sign epistasis and observed peaks related in Local Maxima
ggplot(data = df_out,
      aes(color = Init_Size,
          shape = Mutation,
          y = epist_rsign,
          x = num_observed_peaks,
          )) +    
    facet_grid(Ngenes ~ type_Landscape, scales = "free") +
    geom_point(size = 1) + ##0.5) +
    labs(x = "Number of observed local fitness maxima",
         y = "Fraction reciprocal sign epistasis"
         ) +
    ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))


df_out$numObservedPeaks <- cut_interval(df_out$num_observed_peaks, 2)


## Understanding some splits because of the perfect separation of RMF
## and almost perfect of local max
pdf(file = "gamma-rsign-obs-peaks-FL.pdf", height = 5.8, width = 8.1)
ggplot(data = df_out,
      aes(color = type_Landscape,
          shape = type_Landscape,
          x = epist_rsign,
          y = gamma,
          )) +
    ## scale_y_continuous(trans = log1p_trans()) +
    scale_x_continuous(trans = log1p_trans()) +
    facet_wrap(. ~ numObservedPeaks, scales = "free", labeller = label_both) +
    geom_point(size = 1) + ##0.5) +
    labs(y = TeX("$\\gamma$ (gamma)"), 
         x = "Fraction reciprocal sign epistasis"
         ) ## +
    ## ggtitle("Other epistasis statistics") +
    theme(plot.title = element_text(hjust = 0.5))
dev.off()





## For details about the Fourier coefficients and examples,
## see fourier-examples.R

## The 1st, 2nd, and 3 and above are in
## the short stats, under f[1], [2], f[3+]

## In the former magellan these were called
## w[1], w[2], w[3+]



## Candidate variables

## Average frequency of most frequent genotype [measure of SSWM]
## Number of observed local fitness maxima
## LOD diversity??

## gamma. Actually:
##      only gamma, only W3+, only W2, only fraction reciprocal sign epist.


## Measures of samples
## - sampledGenotypesNumber:
##     Number of different sampled genotypes
## - Mean_muts: Mean number of mutations in sampled genotypes
##       (of course, strong association with detection regime)




## Longer story

## Average frequency of most frequent genotype [measure of SSWM]

## LOD diversity (but not POM diversity): ??: high corr with SSWM with
## small init size and correlated with gamma too in representable and Local max.

## Number of observed local fitness maxima: corr with reciprocal sign
## epist in Local Max fl (not with SSWM, though). Strong to Mild with
## gamma in Local max. And shows no variation in the representable.

## gamma: defined for all fitness landscape not trivially, not correlated
## with SSWM. Moderately correlated with LOD diversity.

## Fraction reciprocal sign epist: only part of the story and strongly
## correlated, in funny changing ways with landscape type, with
## gamma. Also, by decree it is 0 in representable.

## W1, W2, W3: highly correlated among themselves, correlated with gamma,
## correlated with reciprocal sign epist. And each W only a part of the
## picture. W2, W3+ aols affected by number of genes in funny ways (which
## probably makes sense, as number of interactions of genes much larger in
## 10 genes)


## Skewness and kurtosis in number of mutations are far from intuitive, we
## need to refer to a normal, etc. And strong association with number of
## mutations in representable

## - Stdev_muts: standard deviation. Also associated with fitness
##       and sampling regime. But not a simple relationship??
##       Associated to mean muts, specially in RMF and local maxima


