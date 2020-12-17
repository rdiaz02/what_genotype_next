
list_party_packs <- c("partykit", "REEMtree", "evtree", "glmertree", "pmr", "party")
missing_party <- setdiff(list_party_packs, installed.packages())
if(length(missing_party)) install.packages(missing_party)
rm(missing_party, list_party_packs)

date()
library(data.table)
setDTthreads(threads = 0)

library(partykit) ## Don't use cluster, or you'll miss between cluster-unit factors!
## see "clustered-trees.R"
library(glmertree)
library(optimx) ## use different optimizers
library(minqa) ## use different optimizers


## For lme. See timings-lmertree.R
library(RhpcBLASctl)


## simple boxplot
plot_best <- function(fit, main = NULL) {
    partykit:::plot.lmtree(fit$tree, tp_args = list(id = FALSE),
                           main = paste0("\t\t\t", main)) 
}

## several boxplots
plot_method <- function(fit, main = NULL) {
    partykit:::plot.lmtree(fit$tree,
                           tp_args = list(id = FALSE),
                           gp = gpar(cex = 0.7, rot = 90, las = 3),
                           main = paste0("\t\t\t", main))
}

## wrapper for lmertree
##  localfit and localdata should be defined in the script that calls this
##    minsizef unlikely to have a large effect if we do post-pruning,
##    but should lead to faster fits and allow larger depths
tree_best <- function(landscape, depth,
                      fitf = localformula,
                      fitdata = localdata,
                      cores = 9, alpha = 0.01,
                      lmer.control = lmerControl(),
                      weights = FALSE,
                      minsizef = NULL) {

    if (landscape == "Represent.") {
        max_g_7 <- 5
        max_g_10 <- 8
    } else {
        max_g_7 <- 6
        max_g_10 <- 9
    }

    
    if(landscape == "all") {
        ## Does not work: typeLandscape_f added after a set of parentheses
        ## fitf <- update(fitf, ~ . + typeLandscape_f)
        fitf <- as.formula(paste(paste(deparse(fitf), collapse = " "),
                                 " + typeLandscape_f"))
        lscapev <- c("Represent.", "Local maxima", "RMF")
    } else if(landscape == "allFL") {
        fitf <- as.formula(paste(paste(deparse(fitf), collapse = " "),
                                 " + typeLandscape_f"))
        lscapev <- c("Rep", "LM", "RMF")
    } else if (landscape == "all_nofitFL") {
        lscapev <- c("Represent.", "Local maxima", "RMF")
    } else if (landscape == "all_nofitFL_rec") {
        lscapev <- c("Rep", "LM", "RMF")
    } else {
        lscapev <- landscape
    }

    responsev <- all.vars(fitf)[1]
    allvars <- all.vars(fitf)
    ## using sampledFreq lead to numerical problems 
    if(weights) allvars <- c(allvars, "sampledProp")

    datam <- fitdata[
        !is.na(get(responsev)) 
        & (
            ((numGenes == 7) & (nMut <= max_g_7))
            |
            ((numGenes == 10) & (nMut <= max_g_10))
        )
        & (typeLandscape_f %in% lscapev )
       , ..allvars
    ]

    ## cannot turn nMut into factor before that subsetting
    datam[, `:=` (nMut = as.factor(nMut))]
    if("sample_size"  %in% colnames(datam))
        datam[, `:=` (sample_size = as.factor(sample_size))]

    if(any(is.na(datam))) stop("missing values??")
    cat("\n nrow data set = ", nrow(datam), "\n")



    if(!is.null(minsizef)) {
        minsize <- as.integer(round(minsizef * nrow(datam)))
        cat("\n    Using minsize = ", minsize, "\n")
    }

    if(weights) {
        cat("\n   sum of weights = ", sum(datam$sampledProp), "\n")
        if(!all.equal(sum(datam$sampledProp),
                      nrow(datam)))
            message("\n\n sum of weights != nrow  !!!!! \n\n")
    }
    
    ## Yes, repetition, but other options are very cumbersome
    ## as I cannot pass an ifelse expression to weights
    ## nor a character string. Thus the logic above of passing a TRUE/FALSE,
    ## etc

    if(!weights) {
        return(lmertree(fitf,
                  , data = as.data.frame(datam)
                  , maxdepth = depth
                  , weights = NULL
                  , alpha = alpha
                  , cores = cores
                  , minsize = minsize
                  , lmer.control = lmer.control))
    } else {
        return(lmertree(fitf,
                  , data = as.data.frame(datam)
                  , maxdepth = depth
                  , weights = sampledProp
                  , alpha = alpha
                  , cores = cores
                  , minsize = minsize
                  , lmer.control = lmer.control))
    }
}





## fit, give time, print object, and save all objects with root name rootn
tree_best_time_print_save <- function(landscape, depth, rootn, ...) {

    name <- paste0(rootn, "_", gsub(" ", "_", landscape), "_", depth)
    cat(paste("\n\n Fitting landscape = ", landscape,
              "; depth = ", depth, ". Object name = ", name))
    cat("\n")
    print(system.time(assign(name,
                             tree_best(landscape, depth, ...),
                             envir = .GlobalEnv)))

    cat("\n")

    print(get(name, envir = .GlobalEnv)$tree)
    save(file = paste0(rootn, "_fits.RData"),
         list = ls(pattern = glob2rx(paste0(rootn, "*")),
                   envir = .GlobalEnv),
         compress = FALSE)
    cat("\n     done\n")
}







## vertical x axis labels
node_bivplot <- function (mobobj, which = NULL, id = TRUE, pop = TRUE, pointcol = "black", 
    pointcex = 0.5, boxcol = "black", boxwidth = 0.5, boxfill = "lightgray", 
    bg = "white", fitmean = TRUE, linecol = "red", cdplot = FALSE, 
    fivenum = TRUE, breaks = NULL, ylines = NULL, xlab = FALSE, 
    ylab = FALSE, margins = rep(1.5, 4), mainlab = NULL, ...) {
    mf <- model.frame(mobobj)
    y <- Formula::model.part(mobobj$info$Formula, mf, lhs = 1L, 
        rhs = 0L)
    if (isTRUE(ylab)) 
        ylab <- names(y)
    if (identical(ylab, FALSE)) 
        ylab <- ""
    if (is.null(ylines)) 
        ylines <- ifelse(identical(ylab, ""), 0, 2)
    y <- y[[1L]]
    X <- Formula::model.part(mobobj$info$Formula, mf, lhs = 0L, 
        rhs = 1L)
    fitted <- mobobj$fitted[["(fitted)"]]
    if (inherits(X, "try-error")) {
        rval <- switch(class(y)[1L], Surv = node_surv(mobobj, 
            id = id, mainlab = mainlab, ...), factor = node_barplot(mobobj, 
            id = id, mainlab = mainlab, ...), ordered = node_barplot(mobobj, 
            id = id, mainlab = mainlab, ...), node_boxplot(mobobj, 
            ...))
        return(rval)
    }
    if (is.factor(y)) 
        y <- factor(y, levels = rev(levels(y)))
    if (is.null(which)) 
        which <- 1L:NCOL(X)
    X <- X[, which, drop = FALSE]
    k <- NCOL(X)
    xlab <- if (!identical(xlab, FALSE)) {
        if (isTRUE(xlab)) 
            colnames(X)
        else rep(xlab, length.out = k)
    }
    else rep("", k)
    if (is.factor(y)) {
        if (!requireNamespace("vcd")) 
            stop(sprintf("Package %s is required for spine/CD plots", 
                sQuote("vcd")))
        if (cdplot) {
            num_fun <- function(x, y, yfit, i, name, ...) {
                vcd::cd_plot(x, y, xlab = xlab[i], ylab = ylab, 
                  name = name, newpage = FALSE, margins = margins, 
                  pop = FALSE, ...)
                if (fitmean) {
                  grid.lines(x, yfit, default.units = "native", 
                    gp = gpar(col = linecol))
                  if (pop) 
                    popViewport()
                  else upViewport()
                }
                else {
                  if (pop) 
                    popViewport()
                  else upViewport()
                }
            }
        }
        else {
            xscale <- if (is.null(breaks)) {
                if (fivenum) 
                  lapply(X, function(z) {
                    if (is.factor(z)) 
                      1
                    else fivenum(z)
                  })
                else lapply(X, function(z) {
                  if (is.factor(z)) 
                    1
                  else hist(z, plot = FALSE)$breaks
                })
            }
            else {
                if (is.list(breaks)) 
                  breaks
                else list(breaks)
            }
            num_fun <- function(x, y, yfit, i, name, ...) {
                vcd::spine(x, y, xlab = xlab[i], ylab = ylab, 
                  name = name, newpage = FALSE, margins = margins, 
                  pop = FALSE, breaks = xscale[[i]], ...)
                if (fitmean) {
                  xaux <- cut(x, breaks = xscale[[i]], include.lowest = TRUE)
                  yfit <- unlist(tapply(yfit, xaux, mean))
                  xaux <- prop.table(table(xaux))
                  xaux <- cumsum(xaux) - xaux/2
                  grid.lines(xaux, yfit, default.units = "native", 
                    gp = gpar(col = linecol))
                  grid.points(xaux, yfit, default.units = "native", 
                    gp = gpar(col = linecol, cex = pointcex), 
                    pch = 19)
                  if (pop) 
                    popViewport()
                  else upViewport()
                }
                else {
                  if (pop) 
                    popViewport()
                  else upViewport()
                }
            }
        }
        cat_fun <- function(x, y, yfit, i, name, ...) {
            vcd::spine(x, y, xlab = xlab[i], ylab = ylab, name = name, 
                newpage = FALSE, margins = margins, pop = FALSE, 
                ...)
            if (fitmean) {
                yfit <- unlist(tapply(yfit, x, mean))
                xaux <- prop.table(table(x))
                xaux <- cumsum(xaux + 0.02) - xaux/2 - 0.02
                grid.lines(xaux, yfit, default.units = "native", 
                  gp = gpar(col = linecol))
                grid.points(xaux, yfit, default.units = "native", 
                  gp = gpar(col = linecol, cex = pointcex), pch = 19)
                if (pop) 
                  popViewport()
                else upViewport()
            }
            else {
                if (pop) 
                  popViewport()
                else upViewport()
            }
        }
    }
    else {
        xscale <- sapply(X, function(z) {
            if (is.factor(z)) 
                c(1, length(levels(z)))
            else range(z)
        })
        yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
        num_fun <- function(x, y, yfit, i, name, ...) {
            xscale[, i] <- xscale[, i] + c(-0.1, 0.1) * diff(xscale[, 
                i])
            pushViewport(plotViewport(margins = margins, name = name, 
                yscale = yscale, xscale = xscale[, i]))
            grid.points(x, y, gp = gpar(col = pointcol, cex = pointcex))
            if (fitmean) {
                grid.lines(x, yfit, default.units = "native", 
                  gp = gpar(col = linecol))
            }
            grid.xaxis(at = c(ceiling(xscale[1L, i] * 10), floor(xscale[2L, 
                i] * 10))/10)
            grid.yaxis(at = c(ceiling(yscale[1L]), floor(yscale[2L])))
            grid.rect(gp = gpar(fill = "transparent"))
            if (ylab != "") 
                grid.text(ylab, y = unit(0.5, "npc"), x = unit(-2.5, 
                  "lines"), rot = 90)
            if (xlab[i] != "") 
                grid.text(xlab[i], x = unit(0.5, "npc"), y = unit(-2, 
                  "lines"))
            if (pop) 
                popViewport()
            else upViewport()
        }
        cat_fun <- function(x, y, yfit, i, name, ...) {
            xlev <- levels(x)
            pushViewport(plotViewport(margins = margins, name = name, 
                yscale = yscale, xscale = c(0.3, xscale[2L, i] + 
                  0.7)))
            for (j in seq_along(xlev)) {
                by <- boxplot(y[x == xlev[j]], plot = FALSE, las = 3)
                xl <- j - boxwidth/4
                xr <- j + boxwidth/4
                grid.lines(unit(c(xl, xr), "native"), unit(by$stats[1L], 
                  "native"), gp = gpar(col = boxcol))
                grid.lines(unit(j, "native"), unit(by$stats[1L:2L], 
                  "native"), gp = gpar(col = boxcol, lty = 2))
                grid.rect(unit(j, "native"), unit(by$stats[2L], 
                  "native"), width = unit(boxwidth, "native"), 
                  height = unit(diff(by$stats[2:3]), "native"), 
                  just = c("center", "bottom"), gp = gpar(col = boxcol, 
                    fill = boxfill))
                grid.rect(unit(j, "native"), unit(by$stats[3L], 
                  "native"), width = unit(boxwidth, "native"), 
                  height = unit(diff(by$stats[3L:4L]), "native"), 
                  just = c("center", "bottom"), gp = gpar(col = boxcol, 
                    fill = boxfill))
                grid.lines(unit(j, "native"), unit(by$stats[4L:5L], 
                  "native"), gp = gpar(col = boxcol, lty = 2))
                grid.lines(unit(c(xl, xr), "native"), unit(by$stats[5L], 
                  "native"), gp = gpar(col = boxcol))
                n <- length(by$out)
                if (n > 0L) {
                  grid.points(unit(rep.int(j, n), "native"), 
                    unit(by$out, "native"), size = unit(0.5, 
                      "char"), gp = gpar(col = boxcol))
                }
            }
            if (fitmean) {
                yfit <- unlist(tapply(yfit, x, mean))
                grid.lines(seq_along(xlev), yfit, default.units = "native", 
                  gp = gpar(col = linecol))
                grid.points(seq_along(xlev), yfit, default.units = "native", 
                  gp = gpar(col = linecol, cex = pointcex), pch = 19)
            }
            grid.rect(gp = gpar(fill = "transparent"))
            grid.xaxis(at = 1L:length(xlev), label = xlev, ## RDU: next line
                       edits = gEdit(gPath="labels", rot=90))
            grid.yaxis(at = c(ceiling(yscale[1L]), floor(yscale[2L])))
            if (ylab != "") 
                grid.text(ylab, y = unit(0.5, "npc"), x = unit(-3, 
                  "lines"), rot = 90)
            if (xlab[i] != "") 
                grid.text(xlab[i], x = unit(0.5, "npc"), y = unit(-2, 
                  "lines"))
            if (pop) 
                popViewport()
            else upViewport()
        }
    }
    rval <- function(node) {
        nid <- id_node(node)
        ix <- fitted %in% nodeids(mobobj, from = nid, terminal = TRUE)
        y <- y[ix]
        top_vp <- viewport(layout = grid.layout(nrow = k, ncol = 2, 
            widths = unit(c(ylines, 1), c("lines", "null")), 
            heights = unit(k, "null")), width = unit(1, "npc"), 
            height = unit(1, "npc") - unit(2, "lines"), name = paste("node_mob", 
                nid, sep = ""))
        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = bg, col = 0))
        top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
        pushViewport(top)
        if (is.null(mainlab)) {
            mainlab <- if (id) {
                function(id, nobs) sprintf("Node %s (n = %s)", 
                  id, nobs)
            }
            else {
                function(id, nobs) sprintf("n = %s", nobs)
            }
        }
        if (is.function(mainlab)) {
            mainlab <- mainlab(nid, info_node(node)$nobs)
        }
        grid.text(mainlab, y = unit(1, "npc") - unit(0.75, "lines"))
        popViewport()
        for (i in 1L:k) {
            xi <- X[ix, i]
            o <- order(xi)
            yi <- y[o]
            xi <- xi[o]
            yfit <- if (is.null(node$info$object)) {
                fitted(refit.modelparty(mobobj, node = nid))[o]
            }
            else {
                fitted(node$info$object)[o]
            }
            plot_vpi <- viewport(layout.pos.col = 2L, layout.pos.row = i)
            pushViewport(plot_vpi)
            if (is.factor(xi)) 
                cat_fun(xi, yi, yfit, i, paste("node_mob", nid, 
                  "-", i, sep = ""), ...)
            else num_fun(xi, yi, yfit, i, paste("node_mob", nid, 
                "-", i, sep = ""), ...)
            if (pop) 
                popViewport()
            else upViewport()
        }
        if (pop) 
            popViewport()
        else upViewport()
    }
    return(rval)
}
class(node_bivplot) <- "grapcon_generator"
environment(node_bivplot) <- asNamespace('partykit')
assignInNamespace("node_bivplot", node_bivplot, ns = "partykit")


multivar_ctree_plot <- function(x, main = NULL) {
    plot(x, terminal_panel = partykit::node_barplot(x, id = FALSE,
                                          just = "center",
                                          rot = 90,
                                          gp = gpar(cex = .7)),
         main = main)
}


library(codetools)
checkUsageEnv(env = .GlobalEnv)
