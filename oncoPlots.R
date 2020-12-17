library(ggplot2)
library(dplyr)









### GENERAL UTILITY FUNCTIONS

### aggregate table (wrapper)
aggregateTable <- function(df,by,FUN) {
  
  # function to use when aggregating
  f <- function(x) {
    if(is.numeric(x)) return(FUN(x,na.rm=T))
    if(is.character(x)) return(paste(unique(x),collapse=" | "))
    if(is.factor(x)) return(f(as.character(x)))
  }
  
  # aggregate
  x <- aggregate(df,by=by,
                 FUN=f)
  x <- x[,colnames(x) %in% colnames(df)]
  x$replicate <- NULL
  
  return(x)
  
}

### paste sets of factors mantaining their order given a hierarchy
pasteFactors <- function(...,sep=" ",hierarchy=NULL) {
  
  f <- list(...)
  
  if(length(f)>1) {
      
    names(f) <- LETTERS[1:length(f)]
    
    # order of hierarchy
    if (length(hierarchy)) f <- f[hierarchy]
    f <- f[length(f):1]
    
    # generate levels of output
    l <- expand.grid(lapply(lapply(f,as.factor),levels))
    l <- l[,sort(names(f))]
    l <- apply(l,1,paste,collapse=sep)
    
    # paste input
    f <- f[sort(names(f))]
    fout <- as.data.frame(
      matrix(rep(NA,length(f)*max(unlist(lapply(f,length)))),
                   ncol=length(f))
      )
    for (i in 1:length(f)) {
      fout[,i] <- f[[i]]
    }
    fout <- apply(fout,1,paste,collapse=sep)
    fout <- factor(fout,levels=l)
    
    return(fout)
    
  } else return(unlist(...))
  
}









### PLOT GENERATORS

### JS across different methods
oncoPlot_stat_vs_cpm <- function(data=data,stat="js",
                                 ylab="Jensen-Shannon distance",
                                 lineType="solid",lineType_null="dashed",
                                 nullDiff=F) {
  
  # subset and aggregate table
  colsToKeep <- c("numGenes","typeLandscape","detect","size_split","method",stat)
  x <- data[data$sourceGenotype=="any",colsToKeep]
  colnames(x)[colnames(x)==stat] <- "y"
  
  xn <- unique(x[x$method=="null",])
  xn <- aggregateTable(xn,by=list(xn$numGenes,
                                  xn$typeLandscape),
                       FUN=mean)
  
  x <- x[x$method!="null",]
  x <- aggregateTable(x,by=list(x$numGenes,
                                x$typeLandscape,
                                x$detect,
                                x$size_split,
                                x$method),
                      FUN=mean)
  
  # order of factors
  x$size_split <- factor(x$size_split,levels=c(50,200,4000,NA))
  x$detect <- factor(x$detect,levels=c("small","uniform","large",NA))
  x$numGenes <- factor(x$numGenes,levels=c(10,7))
  x$typeLandscape <- factor(x$typeLandscape,levels=c("Represent.",
                                                     "Local maxima",
                                                     "RMF"))
  x$method <- factor(x$method,levels=c("MHN",
                                       "OT",
                                       "CBN",
                                       "MCCBN",
                                       "OT_uw",
                                       "CAPRESE",
                                       "CBN_uw",
                                       "MCCBN_uw",
                                       "CAPRI_AIC",
                                       "CAPRI_BIC",
                                       "MHN_td",
                                       "CBN_td",
                                       "MCCBN_td"))
  
  # grid
  x$xgrid <- pasteFactors(
    pasteFactors("Genes: ",x$numGenes,sep=""),
    pasteFactors(x$typeLandscape),
    sep="\n"
  )
  x$ygrid <- pasteFactors("Detect.: ",x$detect,sep="\n")
  
  # null model
  xn <- unique(merge(xn,
                     unique(x[,c("numGenes","typeLandscape","xgrid","ygrid")])))
  xn$numGenes <- factor(xn$numGenes,levels=c(10,7))
  xn$typeLandscape <- factor(xn$typeLandscape,levels=c("Represent.",
                                                       "Local maxima",
                                                       "RMF"))
  
  # remove unused columns
  x[,c("numGenes","typeLandscape","detect")] <- NULL
  xn[,c("detect","size_split","method","numGenes","typeLandscape")] <- NULL
  
  # calculate relative differences (if nullDiff is TRUE)
  if (nullDiff) {
    for (i in 1:length(levels(x$xgrid))) {
      for (j in 1:length(levels(x$ygrid))) {
        n <- x$xgrid==levels(x$xgrid)[i] & x$ygrid==levels(x$ygrid)[j]
        nn <- xn$xgrid==levels(x$xgrid)[i] & xn$ygrid==levels(x$ygrid)[j]
        
        x$y[n] <- xn$y[nn] - x$y[n]
      }
    }
    ylab <- paste(ylab,"(relative to null model)",sep=" ")
  }
  
  #plot
  p <- ggplot() +
    geom_point(data=x,
               aes(x=method,y=y,color=size_split,group=size_split),
               size=1) +
    geom_line(data=x,
              aes(x=method,y=y,color=size_split,group=size_split),
              linetype=lineType) +
    facet_grid(xgrid ~ ygrid) +
    theme(axis.text.x=element_text(angle=90,
                                   hjust=0.95,
                                   vjust=0.5)) +
    labs(x="",
         y=ylab,
         color="Sample size") +
    scale_color_manual(values=c("orangered",
                                "royalblue",
                                "limegreen"))
  
  if(!nullDiff) {
    p <- p +
      geom_hline(data=xn,
                 aes(yintercept=y),
                 linetype=lineType_null) +l <-
      ylim(0,1)
  } else {
    p <- p +
      ylim(floor(10*min(x$y))/10,
           ceiling(10*max(x$y))/10)
  }
  
  p
  return(p)
  
}

### JS as a function of the # of mutations of the source genotype
### data: num_mut_statistics_over_replicate
oncoPlot_stat_vs_nmut <- function(data,stat="js_w_sampl",
                                  ylab="Jensen-Shannon distance",
                                  lineType="solid",lineType_null="dashed",
                                  nullDiff=F,showLastPoint=T) {
  
  # subset and aggregate table
  colsToKeep <- c("numGenes","typeLandscape","detect","size_split",
                  "method","sourceGenotype_nMut",stat)
  x <- data[,colsToKeep]
  colnames(x)[colnames(x)==stat] <- "y"
  
  # remove last point if specified
  if (!showLastPoint) {
    x <- x[x$numGenes!=x$sourceGenotype_nMut,]
  }
  
  xn <- x[x$method=="null",]
  xn <- aggregateTable(xn,by=list(xn$numGenes,
                                  xn$typeLandscape,
                                  xn$sourceGenotype_nMut),
                      FUN=mean)
  
  x <- x[x$method!="null",]
  x <- aggregateTable(x,by=list(x$numGenes,
                                x$typeLandscape,
                                x$detect,
                                x$size_split,
                                x$method,
                                x$sourceGenotype_nMut),
                      FUN=mean)
  
  # order of factors
  x$size_split <- factor(x$size_split,levels=c(50,200,4000))
  x$numGenes <- factor(x$numGenes,levels=c(10,7))
  x$detect <- factor(x$detect,levels=c("small","uniform","large"))
  x$typeLandscape <- factor(x$typeLandscape,levels=c("Represent.",
                                                     "Local maxima",
                                                     "RMF"))
  x$method <- factor(x$method,levels=c("MHN",
                                       "OT",
                                       "CBN",
                                       "MCCBN",
                                       "OT_uw",
                                       "CAPRESE",
                                       "CBN_uw",
                                       "MCCBN_uw",
                                       "CAPRI_AIC",
                                       "CAPRI_BIC",
                                       "MHN_td",
                                       "CBN_td",
                                       "MCCBN_td"))
  
  # grid
  x$xgrid <- pasteFactors(
    pasteFactors("Genes: ",x$numGenes,sep=""),
    pasteFactors(x$typeLandscape),
    pasteFactors("Detect.: ",x$detect,sep=""),
    sep="\n"
  )
  x$ygrid <- x$method
  
  # null model
  xn <- unique(merge(xn,
                     unique(x[,c("numGenes","typeLandscape","xgrid","ygrid")])))
  
  # remove unused columns
  x[,c("numGenes","typeLandscape","detect","method")] <- NULL
  xn[,c("numGenes","typeLandscape","detect","size_split","method")] <- NULL
  
  # merge (for null model substraction if necessary)
  if (nullDiff) {
    xn <- merge(x,xn,
                by.x=c("xgrid","ygrid","sourceGenotype_nMut"),
                by.y=c("xgrid","ygrid","sourceGenotype_nMut"),
                suffixes=c("",".null"))
  }
  
  # plot
  p <- ggplot() +
    facet_grid(xgrid ~ ygrid) +
    scale_x_continuous(breaks=seq(0,10,by=2)) +
    scale_color_manual(values=c("orangered",
                                "royalblue",
                                "limegreen")) +
    labs(x="Number of mutations of source genotype",
         y=ylab,
         color="Sample size") +
    theme(panel.border=element_rect(colour="black",fill=NA))
  
  if(nullDiff) {
    
    p <- p +
      geom_line(data=xn,
                aes(x=sourceGenotype_nMut,y=(y.null-y),
                    color=size_split,group=size_split),
                linetype=lineType) +
      ylim(floor(10*min(xn$y.null-xn$y))/10,
           ceiling(10*max(xn$y.null-xn$y))/10) +
      geom_hline(yintercept=0,linetype=lineType_null,color="black")
    
  } else {
    
    p <- p + 
      geom_line(data=x,
                aes(x=sourceGenotype_nMut,y=y,
                    color=size_split,group=size_split),
                linetype=lineType) +
      geom_line(data=xn,
                aes(x=sourceGenotype_nMut,y=y),
                linetype=lineType_null,color="black") +
      ylim(0,1)
    
  }
  
  p
  return(p)

}
  
  
## Methods order:
## MHN_td, CBN_td, MCCBN_td, MHN, CBN, MCCBN, OT, ...
  

  






