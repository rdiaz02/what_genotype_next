## grapcon generators. Most of these are modifications of partykit code.
## So that they are available to the code, these needs to be sourced
## (simpler than C-M-x from ESS)

rm(my_node_bp4)
my_node_bp4 <- function(obj, col = "black", fill = "lightgray",
                        bg = "white", 
                        width = 0.5,
                        ## yscale = c(0, 1), ## NULL,
                        ylines = 3,
                        cex = 0.5, ## 0.5
                        id = TRUE,
                        mainlab = NULL,
                        size_x_text = 11,
                        angle_x_text = 45,
                        size_y_text = 11,
                        gp = gpar()) {
    ## cat("\n aqui")
    ## node_id <- id_node(obj)
    ## cat("\n the node_i is ", node_id, "\n")
    ## if oldplot is TRUE, reproduce original behavior
    ## of node_boxplot without adding the barplot
    if(!is.null(terminal_plot_data[[item_in_data]]$weights))
        stop("my_node_bp4 does not deal with weights. Use my_node_bp5, or write code for it")
    cat("\n names of obj are")
    cat(names(obj))
    cat("\n done names")
    
    oldplot <- FALSE
    rval <- function(node){
        ## cat("\n inside rval")
        nid <- id_node(node)
        ## cat("\n\n inside rval nid", nid, "\n")

        ## We could be passing a part of a tree
        ## names has the original numbers
        item_in_data <- names(obj)[nid]
        ## cat("  item_in_data = ", item_in_data, "\n")

        yn <- terminal_plot_data[[item_in_data]]$min_js
        y <- yn
        ## cat(summary(yn))
        ## cat("\n done sum\n")
        ## cat("\n yscale is : ")
        ## yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
        ## cat(yscale)
        yscale <- c(0, 1) + c(-0.1, 0.1) * diff(range(c(0, 1)))
        cat("\n")
        is_best_data <- terminal_plot_data[[item_in_data]]$is_best
        df1 <- as.data.frame(colMeans(is_best_data))
        colnames(df1) <- "height"
        df1$names <- rownames(df1)

        ## ## FIXME: testing what data we have
        ## dat <- data_party(obj, nid)
        ## cat("\n dim dat", dim(dat))
        ## cat("\n colnames dat, ", colnames(dat))
        
        ## yn <- dat[["(response)"]]
        wn <- NULL
        ## wn <- dat[["(weights)"]]
        if (is.null(wn)) 
            wn <- rep(1, length(yn))
        x <- boxplot(rep.int(yn, wn), plot = FALSE)
        if(oldplot) {
            top_vp <- viewport(layout = grid.layout(nrow = 2,
                                                    ncol = 3, 
                                                    widths = unit(c(ylines, 1, 1),
                                                                  c("lines", "null", 
                                                                    "lines")),
                                                    heights = unit(c(1, 1),
                                                                   c("lines",
                                                                     "null"))),
                               width = unit(1, "npc"),
                               height = unit(1, "npc") - unit(2, "lines"),
                               name = paste("node_boxplot", nid, sep = ""),
                               gp = gp)
        } else {
            top_vp <- viewport(layout = grid.layout(nrow = 3, ## 2
                                                    ncol = 3, 
                                                    widths = unit(c(ylines, 1, 1),
                                                                  c("lines", "null", 
                                                                    "lines")),
                                                    ## heights = unit(c(1, 1),
                                                    ##                c("lines",
                                                    ##                  "null"))),
                                                    heights = unit(c(ifelse(id, 2, 1),
                                                                     0.5, 0.5),
                                                                   c("lines",
                                                                     "null",
                                                                     "null"))),
                               width = unit(1, "npc"),
                               height = unit(1, "npc") - unit(2, "lines"),
                               name = paste("node_boxplot", nid, sep = ""),
                               gp = gp)
        }
        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = bg, col = 0))
        top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
        pushViewport(top)
        ## Fixme: obtain nobs from length of y
        if (is.null(mainlab)) {
            mainlab <- if (id) {
                           function(id, nobs) sprintf("Node %s\n (n = %s)", 
                                                      id, nobs)
                       }
                       else {
                           function(id, nobs) sprintf("n = %s", nobs)
                       }
        }
        if (is.function(mainlab)) {
            mainlab <- mainlab(names(obj)[nid], sum(wn))
        }
        grid.text(mainlab)
        popViewport()
        ## cat("\n yscale")
        ## cat(summary(yscale))
        plot <- viewport(layout.pos.col = 2,
                         layout.pos.row = 2, 
                         xscale = c(0, 1), yscale = yscale,
                         ## name = paste0("node_boxplot", 
                         ##               nid, "plot"),
                         clip = FALSE)
        pushViewport(plot)
        grid.yaxis()
        grid.rect(gp = gpar(fill = "transparent"))
        grid.clip()
        xl <- 0.5 - width/4
        xr <- 0.5 + width/4
        grid.lines(unit(c(xl, xr), "npc"), unit(x$stats[1], "native"), 
                   gp = gpar(col = col))
        grid.lines(unit(0.5, "npc"), unit(x$stats[1:2], "native"), 
                   gp = gpar(col = col, lty = 2))
        grid.rect(unit(0.5, "npc"), unit(x$stats[2], "native"), 
                  width = unit(width, "npc"),
                  height = unit(diff(x$stats[c(2, 
                                               4)]), "native"),
                  just = c("center", "bottom"), 
                  gp = gpar(col = col, fill = fill))
        grid.lines(unit(c(0.5 - width/2, 0.5 + width/2), "npc"), 
                   unit(x$stats[3], "native"), gp = gpar(col = col, 
                                                         lwd = 2))
        grid.lines(unit(0.5, "npc"), unit(x$stats[4:5], "native"), 
                   gp = gpar(col = col, lty = 2))

        grid.lines(unit(c(xl, xr), "npc"), unit(x$stats[5], "native"), 
                   gp = gpar(col = col))
        n <- length(x$out)
        if (n > 0) {
            index <- 1:n
            if (length(index) > 0) 
                grid.points(unit(rep.int(0.5, length(index)), 
                  "npc"), unit(x$out[index], "native"), size = unit(cex, 
                  "char"), gp = gpar(col = col))
        }
        if(oldplot) {
            upViewport(2)            
        }
        if(!oldplot) {
            popViewport()
            bottomplot <- viewport(layout.pos.col = 1:2,
                                   layout.pos.row = 3)

            ggp1 <- ggplot(data = df1, aes(x = names, y = height)) + 
                geom_bar(stat = "identity") +
                theme(axis.title.x = element_blank()) +
                theme(axis.title.y = element_blank()) +
                ylim(0, 1) +
                theme(axis.text.x = element_text(
                          size = size_x_text,
                          angle= angle_x_text,
                          vjust = 1,
                          hjust = 1),
                      axis.text.y = element_text(
                          size = size_y_text)
                      )
            
            print(ggp1,
                  vp = bottomplot)
            upViewport(1)  
        }
        
    }
    return(rval)
}
class(my_node_bp4) <- "grapcon_generator"


## Like 4, but rotate the barplot
rm(my_node_bp5)
my_node_bp5 <- function(obj, col = "black", fill = "lightgray",
                        bg = "white", 
                        width = 0.5,
                        ## yscale = NULL,
                        ylines = 3,
                        cex = 0.5, ## 0.5
                        id = TRUE,
                        mainlab = NULL,
                        size_x_text = 11,
                        size_y_text = 11,
                        angle_x_text = 45,
                        Ntot = nrow(obj$fitted),
                        sqrt_transf = FALSE,
                        show_mean = FALSE,
                        show_bottom = "is_best", ## is_best, rank, js
                        gp = gpar()) {
    ## cat("\n aqui")
    ## node_id <- id_node(obj)
    ## cat("\n the node_i is ", node_id, "\n")
    ## if oldplot is TRUE, reproduce original behavior
    ## of node_boxplot without adding the barplot

    cat("\n names of obj are")
    cat(names(obj))
    cat("\n done names")
    
    oldplot <- FALSE
    rval <- function(node){
        ## cat("\n inside rval")
        nid <- id_node(node)
        ## cat("\n\n inside rval nid", nid, "\n")

        ## We could be passing a part of a tree
        ## names has the original numbers
        item_in_data <- names(obj)[nid]
        ## cat("  item_in_data = ", item_in_data, "\n")

        ## Will not work with the split trees
        ## Ntot <- nrow(obj$fitted)
        
        ## Here
        ## if using transformation, backtransform here
        ## And I think I can always rm y completely: FIXME or not create it
        yn <- terminal_plot_data[[item_in_data]]$min_js
        if(sqrt_transf) yn <- yn^2
        ## y <- yn
        ## cat(summary(yn))
        ## cat("\n done sum\n")
        ## yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
        yscale <- c(0, 1) + c(-0.1, 0.1) * diff(range(c(0, 1)))
        is_best_data <- terminal_plot_data[[item_in_data]]$is_best

        
        if(!is.null(terminal_plot_data[[item_in_data]]$weights)) {
            ## cat("\n   doing weighted df for barplot")
            tmpb <- as.matrix(is_best_data)
            df1 <-
                as.data.frame(colWeightedMeans(
                    tmpb,
                    w = terminal_plot_data[[item_in_data]]$weights))
        } else {
            df1 <- as.data.frame(colMeans(is_best_data))
        }
        colnames(df1) <- "height"
        df1$names <- rownames(df1)
        ## cat("  item_in_data = ", item_in_data, ". Printing df1 \n")
        ## print(df1)
        ## testing what data we have
        ## dat <- data_party(obj, nid)
        ## cat("\n dim dat", dim(dat))
        ## cat("\n colnames dat, ", colnames(dat))

        ## Weights
        ## find the smallest value so int mult is OK
        ## ceiling(1/min(weights))
        ## and then, round
        
        if(is.null(terminal_plot_data[[item_in_data]]$weights)) {
            wn <- NULL
            wn <- rep(1, length(yn))
            weighted <- FALSE
        } else {
            ## Ensure all weights a minimum of 1
            wn0 <- terminal_plot_data[[item_in_data]]$weights
            ## cat("\n       sum of weights is ", sum(wn0))
            wn <- round(wn0/min(wn0))
            weighted <- TRUE
        }
        ## yn <- dat[["(response)"]]
        ## wn <- dat[["(weights)"]]
        ## if (is.null(wn)) 
        ##    wn <- rep(1, length(yn))
        x <- boxplot(rep.int(yn, wn), plot = FALSE)
        if(show_mean) meany <- weighted.mean(yn, wn)
        ## as wn used below
        if(weighted) wn <- wn0
        ## The outliers: only the unique values
        ## specially important with weights, o.w., this is slow as hell
        ## (reps of the values)
        x$out <- unique(x$out)
        if(oldplot) {
            top_vp <- viewport(layout = grid.layout(nrow = 2,
                                                    ncol = 3, 
                                                    widths = unit(c(ylines, 1, 1),
                                                                  c("lines", "null", 
                                                                    "lines")),
                                                    heights = unit(c(1, 1),
                                                                   c("lines",
                                                                     "null"))),
                               width = unit(1, "npc"),
                               height = unit(1, "npc") - unit(2, "lines"),
                               name = paste("node_boxplot", nid, sep = ""),
                               gp = gp)
        } else {
            top_vp <- viewport(layout = grid.layout(nrow = 3, ## 2
                                                    ncol = 3, 
                                                    widths = unit(c(ylines, 1, 1),
                                                                  c("lines", "null", 
                                                                    "lines")),
                                                    ## heights = unit(c(1, 1),
                                                    ##                c("lines",
                                                    ##                  "null"))),
                                                    heights = unit(c(ifelse(id, 2, 1),
                                                                     0.5, 0.5),
                                                                   c("lines",
                                                                     "null",
                                                                     "null"))),
                               width = unit(1, "npc"),
                               height = unit(1, "npc") - unit(2, "lines"),
                               name = paste("node_boxplot", nid, sep = ""),
                               gp = gp)
        }
        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = bg, col = 0))
        top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
        pushViewport(top)
        ## Fixme: obtain nobs from length of y
        if (is.null(mainlab)) {
            mainlab <- if (id) {
                           ## function(id, nobs) sprintf("Node %s\n (n = %s)\n (%s%%)", 
                           ##                            id, nobs, round(100 * nobs/Ntot))
                           function(id, nobs) sprintf("Node %s\n (n = %s; %s%%)", 
                                                      id, nobs, round(100 * nobs/Ntot))

                       }
                       else {
                           function(id, nobs) sprintf("n = %s", nobs)
                       }
        }
        if (is.function(mainlab)) {
            ## cat("\n    nid =", nid, " In this mainlab we have: Ntot = ", Ntot,
            ##     ".  round(sum(wn)) ", round(sum(wn)), "\n")
            mainlab <- mainlab(names(obj)[nid], round(sum(wn)))
            ## mainlab <- mainlab(names(obj)[nid], round(sum(wn)),  round(100 * sum(wn)/Ntot))
        }
        grid.text(mainlab)
        popViewport()
        ## cat("\n yscale")
        ## cat(summary(yscale))
        plot <- viewport(layout.pos.col = 2,
                         layout.pos.row = 2, 
                         xscale = c(0, 1), yscale = yscale,
                         ## name = paste0("node_boxplot", 
                         ##               nid, "plot"),
                         clip = FALSE)
        pushViewport(plot)
        grid.yaxis()
        grid.rect(gp = gpar(fill = "transparent"))
        grid.clip()
        xl <- 0.5 - width/4
        xr <- 0.5 + width/4
        grid.lines(unit(c(xl, xr), "npc"), unit(x$stats[1], "native"), 
                   gp = gpar(col = col))
        grid.lines(unit(0.5, "npc"), unit(x$stats[1:2], "native"), 
                   gp = gpar(col = col, lty = 2))
        grid.rect(unit(0.5, "npc"), unit(x$stats[2], "native"), 
                  width = unit(width, "npc"),
                  height = unit(diff(x$stats[c(2, 
                                               4)]), "native"),
                  just = c("center", "bottom"), 
                  gp = gpar(col = col, fill = fill))
        grid.lines(unit(c(0.5 - width/2, 0.5 + width/2), "npc"), 
                   unit(x$stats[3], "native"), gp = gpar(col = col, 
                                                         lwd = 2))
        grid.lines(unit(0.5, "npc"), unit(x$stats[4:5], "native"), 
                   gp = gpar(col = col, lty = 2))

        grid.lines(unit(c(xl, xr), "npc"), unit(x$stats[5], "native"), 
                   gp = gpar(col = col))
        if(show_mean) {
            grid.lines(unit(c(.5 - width, .5 + width), "npc"), unit(meany, "native"), 
                       gp = gpar(col = "blue", lty = 2))
            grid.lines(unit(c(xl, xr), "npc"), unit(x$stats[3], "native"), 
                       gp = gpar(col = "black", lty = 3))
        }
        n <- length(x$out)
        ## showing points misleading if weighted?
        if ((n > 0) ) { # & !weighted) {
            index <- 1:n
            if (length(index) > 0) 
                grid.points(unit(rep.int(0.5, length(index)), 
                  "npc"), unit(x$out[index], "native"), size = unit(cex, 
                  "char"), gp = gpar(col = col))
        }
        if(oldplot) {
            upViewport(2)            
        }
        if(!oldplot) {
            popViewport()
            bottomplot <- viewport(layout.pos.col = 1:2,
                                   layout.pos.row = 3)
            if(show_bottom == "is_best") {
                ggp1 <- ggplot(data = df1, aes(x = names, y = height)) + 
                    geom_bar(stat = "identity")  +
                    ylim(0, 1) +
                    ## scale_y_reverse(limits = c(1, 0)) +
                    coord_flip() 
            } else {
                ggp1 <- ggplot(data = df1, aes(x = names, y = height)) + 
                    geom_point(size = 1) +
                    coord_flip()
                if(show_bottom == "rank")
                    ggp1 <- ggp1 + ylim(1, 13) ## + scale_y_reverse(limits = c(13, 1))
                if(show_bottom == "js")
                    ggp1 <- ggp1 + ylim(0, 1) ## + scale_y_reverse(limits = c(1, 0))
            }

            ggp1 <- ggp1 +
                    theme(axis.title.x = element_blank()) +
                    theme(axis.title.y = element_blank()) +
                    theme(axis.text.x = element_text(
                              size = size_x_text,
                              angle = angle_x_text,
                              vjust = 1,
                              hjust = 1),
                          axis.text.y = element_text(
                              size = size_y_text)
                          ) +
                theme(panel.grid.major =
                          element_line(size = (0.2), colour="grey"))
                
            print(ggp1,
                  vp = bottomplot)
            upViewport(1)  
        }
        
    }
    return(rval)
}
class(my_node_bp5) <- "grapcon_generator"




## using .list.rules.party in
## From: https://stackoverflow.com/a/41976697/3670562
## this is a slightly modified version, with parenteses
my.list.rules.party <- function(x, i = NULL, ...) {
    if (is.null(i)) i <- nodeids(x, terminal = TRUE)
    if (length(i) > 1) {
        ret <- sapply(i, my.list.rules.party, x = x)
        names(ret) <- if (is.character(i)) i else names(x)[i]
        return(ret)
    }
    if (is.character(i) && !is.null(names(x)))
        i <- which(names(x) %in% i)
    stopifnot(length(i) == 1 & is.numeric(i))
    stopifnot(i <= length(x) & i >= 1)
    i <- as.integer(i)
    dat <- data_party(x, i)  
    if (!is.null(x$fitted)) {
        findx <- which("(fitted)" == names(dat))[1]  
        fit <- dat[,findx:ncol(dat), drop = FALSE]   
        dat <- dat[,-(findx:ncol(dat)), drop = FALSE]
        if (ncol(dat) == 0)
            dat <- x$data
    } else {
        fit <- NULL  
        dat <- x$data
    }

    rule <- c()

    recFun <- function(node) {
        if (id_node(node) == i) return(NULL)   
        kid <- sapply(kids_node(node), id_node)
        whichkid <- max(which(kid <= i))
        split <- split_node(node)
        ivar <- varid_split(split)
        svar <- names(dat)[ivar]
        index <- index_split(split)
        if (is.factor(dat[, svar])) {
            if (is.null(index)) 
                index <- ((1:nlevels(dat[, svar])) > breaks_split(split)) + 1
            slevels <- levels(dat[, svar])[index == whichkid]
            srule <- paste(svar, " %in% c(\"", 
                paste(slevels, collapse = "\", \"", sep = ""), "\")",
                sep = "")
        } else {
            if (is.null(index)) index <- 1:length(kid)
            breaks <- cbind(c(-Inf, breaks_split(split)), 
                            c(breaks_split(split), Inf))
            sbreak <- breaks[index == whichkid,]
            right <- right_split(split)
            srule <- c()
            if (is.finite(sbreak[1]))
                srule <- c(srule, 
                    paste(svar, ifelse(right, ">", ">="), sbreak[1]))
            if (is.finite(sbreak[2]))
                srule <- c(srule, 
                    paste(svar, ifelse(right, "<=", "<"), sbreak[2]))
            srule <- paste(srule, collapse = " & ")
        }
        rule <<- c(rule, srule)
        return(recFun(node[[whichkid]]))
    }
    node <- recFun(node_party(x))
    rule <- paste("(", rule, ")")
    paste(rule, collapse = " & ")
}







rm(my_node_empty)
my_node_empty <- function(obj, col = "black", fill = "lightgray",
                        bg = "white", 
                        width = 0.5, yscale = NULL, ylines = 3,
                        cex = 0.5, id = FALSE, ## TRUE
                        mainlab = NULL, gp = gpar()) {
    cat("\n aqui\n")
    ## node_id <- id_node(obj)
    ## cat("\n the node_i is ", node_id, "\n")
    oldplot <- FALSE
    rval <- function(node){
        ## cat("\n inside rval")
        nid <- id_node(node)
        ## cat("\n inside rval nid", nid, "\n")
    }
    
    ##     cat("\n inside rval")
    ##     nid <- id_node(node)
    ##     cat("\n inside rval nid", nid, "\n")
    ##     y <- terminal_plot_data[[nid]]$min_js
    ##     cat(summary(y))
    ##     cat("\n done sum\n")
    ##     is_best_data <- terminal_plot_data[[nid]]$is_best

    ##     df1 <- as.data.frame(colMeans(is_best_data))
    ##     colnames(df1) <- "height"
    ##     df1$names <- rownames(df1)

        
    ##     yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
    ##     ## dat <- data_party(obj, nid)
    ##     yn <- y
    ##     ## yn <- dat[["(response)"]]
    ##     wn <- NULL
    ##     ## wn <- dat[["(weights)"]]
    ##     if (is.null(wn)) 
    ##         wn <- rep(1, length(yn))
    ##     x <- boxplot(rep.int(yn, wn), plot = FALSE)
    ##     if(oldplot) {
    ##         top_vp <- viewport(layout = grid.layout(nrow = 2,
    ##                                                 ncol = 3, 
    ##                                                 widths = unit(c(ylines, 1, 1),
    ##                                                               c("lines", "null", 
    ##                                                                 "lines")),
    ##                                                 heights = unit(c(1, 1),
    ##                                                                c("lines",
    ##                                                                  "null"))),
    ##                            width = unit(1, "npc"),
    ##                            height = unit(1, "npc") - unit(2, "lines"),
    ##                            name = paste("node_boxplot", nid, sep = ""),
    ##                            gp = gp)
    ##     } else {
    ##         top_vp <- viewport(layout = grid.layout(nrow = 3, ## 2
    ##                                                 ncol = 3, 
    ##                                                 widths = unit(c(ylines, 1, 1),
    ##                                                               c("lines", "null", 
    ##                                                                 "lines")),
    ##                                                 ## heights = unit(c(1, 1),
    ##                                                 ##                c("lines",
    ##                                                 ##                  "null"))),
    ##                                                 heights = unit(c(1, 0.5, 0.5),
    ##                                                                c("lines",
    ##                                                                  "null",
    ##                                                                  "null"))),
    ##                            width = unit(1, "npc"),
    ##                            height = unit(1, "npc") - unit(2, "lines"),
    ##                            name = paste("node_boxplot", nid, sep = ""),
    ##                            gp = gp)
    ##     }
    ##     pushViewport(top_vp)
    ##     grid.rect(gp = gpar(fill = bg, col = 0))
    ##     top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
    ##     pushViewport(top)
    ##     ## Fixme: obtain nobs from length of y
    ##     if (is.null(mainlab)) {
    ##         mainlab <- if (id) {
    ##                        function(id, nobs) sprintf("Node %s (n = %s)", 
    ##                                                   id, nobs)
    ##                    }
    ##                    else {
    ##                        function(id, nobs) sprintf("n = %s", nobs)
    ##                    }
    ##     }
    ##     if (is.function(mainlab)) {
    ##         mainlab <- mainlab(names(obj)[nid], sum(wn))
    ##     }
    ##     grid.text(mainlab)
    ##     popViewport()
    ##     plot <- viewport(layout.pos.col = 2,
    ##                      layout.pos.row = 2, 
    ##                      xscale = c(0, 1), yscale = yscale,
    ##                      name = paste0("node_boxplot", 
    ##                                    nid, "plot"), clip = FALSE)
    ##     pushViewport(plot)
    ##     grid.yaxis()
    ##     grid.rect(gp = gpar(fill = "transparent"))
    ##     grid.clip()
    ##     xl <- 0.5 - width/4
    ##     xr <- 0.5 + width/4
    ##     grid.lines(unit(c(xl, xr), "npc"), unit(x$stats[1], "native"), 
    ##                gp = gpar(col = col))
    ##     grid.lines(unit(0.5, "npc"), unit(x$stats[1:2], "native"), 
    ##                gp = gpar(col = col, lty = 2))
    ##     grid.rect(unit(0.5, "npc"), unit(x$stats[2], "native"), 
    ##               width = unit(width, "npc"),
    ##               height = unit(diff(x$stats[c(2, 
    ##                                            4)]), "native"),
    ##               just = c("center", "bottom"), 
    ##               gp = gpar(col = col, fill = fill))
    ##     grid.lines(unit(c(0.5 - width/2, 0.5 + width/2), "npc"), 
    ##                unit(x$stats[3], "native"), gp = gpar(col = col, 
    ##                                                      lwd = 2))
    ##     grid.lines(unit(0.5, "npc"), unit(x$stats[4:5], "native"), 
    ##                gp = gpar(col = col, lty = 2))

    ##     grid.lines(unit(c(xl, xr), "npc"), unit(x$stats[5], "native"), 
    ##                gp = gpar(col = col))

    ##     if(oldplot) {
    ##         upViewport(2)            
    ##     }
    ##     if(!oldplot) {
    ##         popViewport()
    ##         bottomplot <- viewport(layout.pos.col = 1:2,
    ##                                layout.pos.row = 3)

    ##         ggp1 <- ggplot(data = df1, aes(x = names, y = height)) + 
    ##             geom_bar(stat = "identity") +
    ##             theme(axis.title.x = element_blank()) +
    ##             theme(axis.title.y = element_blank()) +
    ##             ylim(0, 1) +
    ##             theme(axis.text.x = element_text(angle= 45,
    ##                                              vjust = 1,
    ##                                              hjust = 1))
            
    ##         print(ggp1,
    ##               vp = bottomplot)
    ##         upViewport(1)  
    ##     }
        
    ## }
    return(rval)
}
class(my_node_empty) <- "grapcon_generator"



rm(my_node_bp_original)
## straight copy of partykit::node_boxplot
## commented out on top
my_node_bp_original <- function (obj, col = "black",
                                 fill = "lightgray", bg = "white", 
    width = 0.5, yscale = NULL, ylines = 3, cex = 0.5, id = TRUE, 
    mainlab = NULL, gp = gpar()) 
{
    ## y <- obj$fitted[["(response)"]]
    ## stopifnot(is.numeric(y))
    ## if (is.null(yscale)) 
    ##     yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
    rval <- function(node) {
        nid <- id_node(node)
        cat("\n I am at node ",  nid , "\n")
        dat <- data_party(obj, nid)
        yn <- dat[["(response)"]]
        cat(summary(yn))
        y <- yn
        yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
        wn <- dat[["(weights)"]]
        if (is.null(wn)) 
            wn <- rep(1, length(yn))
        x <- boxplot(rep.int(yn, wn), plot = FALSE)
        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3, 
            widths = unit(c(ylines, 1, 1), c("lines", "null", 
                "lines")), heights = unit(c(1, 1), c("lines", 
                "null"))), width = unit(1, "npc"), height = unit(1, 
            "npc") - unit(2, "lines"), name = paste("node_boxplot", 
            nid, sep = ""), gp = gp)
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
            mainlab <- mainlab(names(obj)[nid], sum(wn))
        }
        grid.text(mainlab)
        popViewport()
        ## cat("\n yscale")
        ## cat(summary(yscale))
        plot <- viewport(layout.pos.col = 2, layout.pos.row = 2, 
                         xscale = c(0, 1),
                         yscale = yscale,
                         ## name = paste0("node_boxplot", 
                         ##               nid, "plot"),
                         clip = FALSE)
        pushViewport(plot)
        grid.yaxis()
        grid.rect(gp = gpar(fill = "transparent"))
        grid.clip()
        xl <- 0.5 - width/4
        xr <- 0.5 + width/4
        grid.lines(unit(c(xl, xr), "npc"),
                   unit(x$stats[1], "native"), 
            gp = gpar(col = col))
        grid.lines(unit(0.5, "npc"), unit(x$stats[1:2], "native"), 
            gp = gpar(col = col, lty = 2))
        grid.rect(unit(0.5, "npc"), unit(x$stats[2], "native"), 
            width = unit(width, "npc"), height = unit(diff(x$stats[c(2, 
                4)]), "native"), just = c("center", "bottom"), 
            gp = gpar(col = col, fill = fill))
        grid.lines(unit(c(0.5 - width/2, 0.5 + width/2), "npc"), 
            unit(x$stats[3], "native"), gp = gpar(col = col, 
                lwd = 2))
        grid.lines(unit(0.5, "npc"), unit(x$stats[4:5], "native"), 
            gp = gpar(col = col, lty = 2))
        grid.lines(unit(c(xl, xr), "npc"), unit(x$stats[5], "native"), 
            gp = gpar(col = col))
        n <- length(x$out)
        if (n > 0) {
            index <- 1:n
            if (length(index) > 0) 
                grid.points(unit(rep.int(0.5, length(index)), 
                  "npc"), unit(x$out[index], "native"), size = unit(cex, 
                  "char"), gp = gpar(col = col))
        }
        upViewport(2)
    }
    return(rval)
}
class(my_node_bp_original) <- "grapcon_generator"





## rm(my_node_bp4)
## my_node_number <- function(obj, col = "black", fill = "lightgray",
##                            bg = "white", 
##                            width = 0.5, yscale = NULL, ylines = 3,
##                            cex = 0.5, id = TRUE,
##                            mainlab = NULL, gp = gpar()) {
##     cat("\n aqui")
##     ## node_id <- id_node(obj)
##     ## cat("\n the node_i is ", node_id, "\n")
##     ## if oldplot is TRUE, reproduce original behavior
##     ## of node_boxplot without adding the barplot

##     cat("\n names of obj are")
##     cat(names(obj))
##     cat("\n done names")
    
##     oldplot <- FALSE
##     rval <- function(node){
##         ## cat("\n inside rval")
##         nid <- id_node(node)
##         cat("\n\n inside rval nid", nid, "\n")

##         ## We could be passing a part of a tree
##         ## names has the original numbers
##         item_in_data <- names(obj)[nid]
##         cat("  item_in_data = ", item_in_data, "\n")

##         yn <- terminal_plot_data[[item_in_data]]$min_js
##         y <- yn
##         ## cat(summary(yn))
##         ## cat("\n done sum\n")
##         yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))

##         is_best_data <- terminal_plot_data[[item_in_data]]$is_best
##         df1 <- as.data.frame(colMeans(is_best_data))
##         colnames(df1) <- "height"
##         df1$names <- rownames(df1)

        
##         ## dat <- data_party(obj, nid)
##         ## yn <- dat[["(response)"]]
##         wn <- NULL
##         ## wn <- dat[["(weights)"]]
##         if (is.null(wn)) 
##             wn <- rep(1, length(yn))
##         x <- boxplot(rep.int(yn, wn), plot = FALSE)
##         if(oldplot) {
##             top_vp <- viewport(layout = grid.layout(nrow = 2,
##                                                     ncol = 3, 
##                                                     widths = unit(c(ylines, 1, 1),
##                                                                   c("lines", "null", 
##                                                                     "lines")),
##                                                     heights = unit(c(1, 1),
##                                                                    c("lines",
##                                                                      "null"))),
##                                width = unit(1, "npc"),
##                                height = unit(1, "npc") - unit(2, "lines"),
##                                name = paste("node_boxplot", nid, sep = ""),
##                                gp = gp)
##         } else {
##             top_vp <- viewport(layout = grid.layout(nrow = 3, ## 2
##                                                     ncol = 3, 
##                                                     widths = unit(c(ylines, 1, 1),
##                                                                   c("lines", "null", 
##                                                                     "lines")),
##                                                     ## heights = unit(c(1, 1),
##                                                     ##                c("lines",
##                                                     ##                  "null"))),
##                                                     heights = unit(c(ifelse(id, 2, 1),
##                                                                      0.5, 0.5),
##                                                                    c("lines",
##                                                                      "null",
##                                                                      "null"))),
##                                width = unit(1, "npc"),
##                                height = unit(1, "npc") - unit(2, "lines"),
##                                name = paste("node_boxplot", nid, sep = ""),
##                                gp = gp)
##         }
##         pushViewport(top_vp)
##         grid.rect(gp = gpar(fill = bg, col = 0))
##         top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
##         pushViewport(top)
##         ## Fixme: obtain nobs from length of y
##         if (is.null(mainlab)) {
##             mainlab <- if (id) {
##                            function(id, nobs) sprintf("Node %s\n (n = %s)", 
##                                                       id, nobs)
##                        }
##                        else {
##                            function(id, nobs) sprintf("n = %s", nobs)
##                        }
##         }
##         if (is.function(mainlab)) {
##             mainlab <- mainlab(names(obj)[nid], sum(wn))
##         }
##         grid.text(mainlab)
##         popViewport()
##             upViewport(1)  
##         }
##     }
##     return(rval)
## }
## class(my_node_number) <- "grapcon_generator"


## from node_inner
my_node_number2 <- function (obj, id = TRUE,
                             pval = FALSE,
                             abbreviate = FALSE, fill = "white", 
                             gp = gpar(),
                             yunit = 1) {
    meta <- obj$data
    nam <- names(obj)
    extract_label <- function(node) {
        if (is.terminal(node)) 
            return(rep.int("", 2L))
        varlab <- character_split(split_node(node), meta)$name
        if (abbreviate > 0L) 
            varlab <- abbreviate(varlab, as.integer(abbreviate))
        if (pval) {
            nullna <- function(x) is.null(x) || is.na(x)
            pval <- suppressWarnings(try(!nullna(info_node(node)$p.value), 
                silent = TRUE))
            pval <- if (inherits(pval, "try-error")) 
                FALSE
            else pval
        }
        if (pval) {
            pvalue <- node$info$p.value
            plab <- ifelse(pvalue < 10^(-3L), paste("p <", 10^(-3L)), 
                paste("p =", round(pvalue, digits = 3L)))
        }
        else {
            plab <- ""
        }
        return(c(varlab, plab))
    }
    maxstr <- function(node) {
        lab <- extract_label(node)
        klab <- if (is.terminal(node)) 
            ""
        else unlist(lapply(kids_node(node), maxstr))
        lab <- c(lab, klab)
        lab <- unlist(lapply(lab, function(x) strsplit(x, "\n")))
        lab <- lab[which.max(nchar(lab))]
        if (length(lab) < 1L) 
            lab <- ""
        return(lab)
    }
    nstr <- maxstr(node_party(obj))
    if (nchar(nstr) < 6) 
        nstr <- "aAAAAa"
    rval <- function(node) {
        pushViewport(viewport(gp = gp, name = paste("node_inner", 
            id_node(node), "_gpar", sep = "")))
        node_vp <- viewport(x = unit(0.5, "npc"), y = unit(0.5, 
            "npc"), width = unit(1, "strwidth", nstr) * 1.3, 
            height = unit(3, "lines"), name = paste("node_inner", 
                id_node(node), sep = ""), gp = gp)
        pushViewport(node_vp)
        xell <- c(seq(0, 0.2, by = 0.01), seq(0.2, 0.8, by = 0.05), 
            seq(0.8, 1, by = 0.01))
        yell <- sqrt(xell * (1 - xell))
        lab <- extract_label(node)
        fill <- rep(fill, length.out = 2L)
        
        ## grid.polygon(x = unit(c(xell, rev(xell)), "npc"), y = unit(c(yell, 
        ##     -yell) + 0.5, "npc"), gp = gpar(fill = fill[1]))
        ## grid.text(lab[1L], y = unit(1.5 + 0.5 * (lab[2L] != ""), 
        ##     "lines"))
        if (lab[2L] != "") 
            grid.text(lab[2L], y = unit(1, "lines"))
        if (id) {
            nodeIDvp <- viewport(x = unit(0.5, "npc"),
                                 y = unit(yunit, #1, ## FIXME: when they leave blans
                "npc"), width = max(unit(1, "lines"), unit(1.3, 
                "strwidth", nam[id_node(node)])), height = max(unit(1, 
                "lines"), unit(1.3, "strheight", nam[id_node(node)])))
            pushViewport(nodeIDvp)
            grid.rect(gp = gpar(fill = fill[2]))
            grid.text(nam[id_node(node)])
            popViewport()
        }
        upViewport(2)
    }
    return(rval)
}
class(my_node_number2) <- "grapcon_generator"




my_node_inner1 <- function (obj, id = TRUE, pval = TRUE, abbreviate = FALSE, fill = "white", 
    gp = gpar()) 
{
    meta <- obj$data
    nam <- names(obj)
    extract_label <- function(node) {
        if (is.terminal(node)) 
            return(rep.int("", 2L))
        varlab <- character_split(split_node(node), meta)$name
        if (abbreviate > 0L) 
            varlab <- abbreviate(varlab, as.integer(abbreviate))
        if (pval) {
            nullna <- function(x) is.null(x) || is.na(x)
            pval <- suppressWarnings(try(!nullna(info_node(node)$p.value), 
                silent = TRUE))
            pval <- if (inherits(pval, "try-error")) 
                FALSE
            else pval
        }
        if (pval) {
            pvalue <- node$info$p.value
            plab <- ifelse(pvalue < 10^(-3L), paste("p <", 10^(-3L)), 
                paste("p =", round(pvalue, digits = 3L)))
        }
        else {
            plab <- ""
        }
        return(c(varlab, plab))
    }
    maxstr <- function(node) {
        lab <- extract_label(node)
        klab <- if (is.terminal(node)) 
            ""
        else unlist(lapply(kids_node(node), maxstr))
        lab <- c(lab, klab)
        lab <- unlist(lapply(lab, function(x) strsplit(x, "\n")))
        lab <- lab[which.max(nchar(lab))]
        if (length(lab) < 1L) 
            lab <- ""
        return(lab)
    }
    nstr <- maxstr(node_party(obj))
    if (nchar(nstr) < 6) 
        nstr <- "aAAAAa"
    rval <- function(node) {
        ## FIXME: testing what data we have
        nid <- id_node(node)
        cat("\n\n inside internal node nid", nid, "\n")
        ## We could be passing a part of a tree
        ## names has the original numbers
        item_in_data <- names(obj)[nid]
        cat("  ... internal node: item_in_data = ", item_in_data, "\n")
        dat <- data_party(obj, nid)
        cat("\n ...  internal node: dim dat", dim(dat))
        cat("\n ...  internal node: colnames dat, ", colnames(dat))
        ## end testing

        
        pushViewport(viewport(gp = gp, name = paste("node_inner", 
            id_node(node), "_gpar", sep = "")))
        node_vp <- viewport(x = unit(0.5, "npc"), y = unit(0.5, 
            "npc"), width = unit(1, "strwidth", nstr) * 1.3, 
            height = unit(3, "lines"), name = paste("node_inner", 
                id_node(node), sep = ""), gp = gp)
        pushViewport(node_vp)
        xell <- c(seq(0, 0.2, by = 0.01), seq(0.2, 0.8, by = 0.05), 
            seq(0.8, 1, by = 0.01))
        yell <- sqrt(xell * (1 - xell))
        lab <- extract_label(node)
        fill <- rep(fill, length.out = 2L)
        grid.polygon(x = unit(c(xell, rev(xell)), "npc"), y = unit(c(yell, 
            -yell) + 0.5, "npc"), gp = gpar(fill = fill[1]))
        grid.text(lab[1L], y = unit(1.5 + 0.5 * (lab[2L] != ""), 
            "lines"))
        if (lab[2L] != "") 
            grid.text(lab[2L], y = unit(1, "lines"))
        if (id) {
            nodeIDvp <- viewport(x = unit(0.5, "npc"), y = unit(1, 
                "npc"), width = max(unit(1, "lines"), unit(1.3, 
                "strwidth", nam[id_node(node)])), height = max(unit(1, 
                "lines"), unit(1.3, "strheight", nam[id_node(node)])))
            pushViewport(nodeIDvp)
            grid.rect(gp = gpar(fill = fill[2]))
            grid.text(nam[id_node(node)])
            popViewport()
        }
        upViewport(2)
    }
    return(rval)
}
class(my_node_inner1) <- "grapcon_generator"





my_node_inner2 <- function (obj, id = TRUE, pval = TRUE,
                            abbreviate = FALSE, fill = "white",
                            width_inner = 14, height_inner = 11,
                            ## size_x_text = , angle_x_text = 
                            gp = gpar()) 
{
    meta <- obj$data
    nam <- names(obj)
    extract_label <- function(node) {
        if (is.terminal(node)) 
            return(rep.int("", 2L))
        varlab <- character_split(split_node(node), meta)$name
        if (abbreviate > 0L) 
            varlab <- abbreviate(varlab, as.integer(abbreviate))
        if (pval) {
            nullna <- function(x) is.null(x) || is.na(x)
            pval <- suppressWarnings(try(!nullna(info_node(node)$p.value), 
                silent = TRUE))
            pval <- if (inherits(pval, "try-error")) 
                FALSE
            else pval
        }
        if (pval) {
            pvalue <- node$info$p.value
            plab <- ifelse(pvalue < 10^(-3L), paste("p <", 10^(-3L)), 
                paste("p =", round(pvalue, digits = 3L)))
        }
        else {
            plab <- ""
        }
        return(c(varlab, plab))
    }
    maxstr <- function(node) {
        lab <- extract_label(node)
        klab <- if (is.terminal(node)) 
            ""
        else unlist(lapply(kids_node(node), maxstr))
        lab <- c(lab, klab)
        lab <- unlist(lapply(lab, function(x) strsplit(x, "\n")))
        lab <- lab[which.max(nchar(lab))]
        if (length(lab) < 1L) 
            lab <- ""
        return(lab)
    }
    nstr <- maxstr(node_party(obj))
    if (nchar(nstr) < 6) 
        nstr <- "aAAAAa"

    rval <- function(node) { ##, obj) {
        nid <- id_node(node)
        item_in_data <- names(obj)[nid]
        cat("\n node ", nid)
        spn <- split_node(node)
        ## addinfo <- info_node(node)
        ## print(addinfo)
        
        terminal <- numerical <- categorical <- FALSE
        if(!is.null(spn$breaks)) {
            numerical <- TRUE
            if(length(spn$breaks) > 1) stop("length breaks > 1")
        } 
        if(!is.null(spn$index)) {
            categorical <- TRUE
        }
        if(is.null(spn$index) & is.null(spn$breaks)) terminal <- TRUE
        
        ## cat("\n      numerical = ", numerical)
        ## cat( "  terminal = ", terminal)
        ## FIXME_zu 
        ## FIXME: change assignment of this_data
        ## by that in my_node_bp5
        ##  coco <- terminal_plot_data[[item_in_data]][[]]
        ## Nope, terminal_plot_data as of now only contains terminal data
        
        this_data <- data_party(obj, nid)
        this_var <- this_data[, spn$varid]
        if(is.factor(this_var) & !categorical) stop("factor and not categorical")

        ## ## why? from their code
        ## meta <- obj$data
        ## varlab <- character_split(split_node(node), meta)$name

        varlab2 <- colnames(this_data)[spn$varid]
        ## cat("\n     length this var is = ", length(this_var))
        ## cat("\n     name of the var is ", varlab)
        ## cat("\n     name2 of the var is ", varlab2)

        pushViewport(viewport(gp = gp,
                              name = paste("node_inner", 
                                           id_node(node), "_gpar", sep = "")))
        node_vp <- viewport(x = unit(0.5, "npc"),
                            y = unit(0.5, "npc"),
                            width = unit(width_inner, "lines"), ## 7                            
                            ## width = unit(1, "strwidth", nstr) * 1.3, 
                            height = unit(height_inner, "lines"), ## unit(3, "lines"),
                            name = paste("node_inner", 
                                         id_node(node), sep = ""),
                            gp = gp)
        ## pushViewport(node_vp)   
        if(numerical) {
            ## print(summary(this_var))
                print(density_color_3(this_var, spn$breaks, varname = varlab2,
                                      nodenum = item_in_data),
                      vp = node_vp)
        }
        if(categorical) {
            print(categorical_plot(this_var, spn$index, varname = varlab2,
                                   nodenum = item_in_data),
                  vp = node_vp)
        }
        upViewport(1)
    }

    
    ## rval <- function(node) {
    ##     ## FIXME: testing what data we have
    ##     nid <- id_node(node)
    ##     cat("\n\n inside internal node nid", nid, "\n")
    ##     ## We could be passing a part of a tree
    ##     ## names has the original numbers
    ##     item_in_data <- names(obj)[nid]
    ##     cat("  ... internal node: item_in_data = ", item_in_data, "\n")
    ##     dat <- data_party(obj, nid)
    ##     cat("\n ...  internal node: dim dat", dim(dat))
    ##     cat("\n ...  internal node: colnames dat, ", colnames(dat))
    ##     ## end testing

        
    ##     pushViewport(viewport(gp = gp, name = paste("node_inner", 
    ##         id_node(node), "_gpar", sep = "")))
    ##     node_vp <- viewport(x = unit(0.5, "npc"), y = unit(0.5, 
    ##         "npc"), width = unit(1, "strwidth", nstr) * 1.3, 
    ##         height = unit(3, "lines"), name = paste("node_inner", 
    ##             id_node(node), sep = ""), gp = gp)
    ##     pushViewport(node_vp)
    ##     xell <- c(seq(0, 0.2, by = 0.01), seq(0.2, 0.8, by = 0.05), 
    ##         seq(0.8, 1, by = 0.01))
    ##     yell <- sqrt(xell * (1 - xell))
    ##     lab <- extract_label(node)
    ##     fill <- rep(fill, length.out = 2L)
    ##     grid.polygon(x = unit(c(xell, rev(xell)), "npc"), y = unit(c(yell, 
    ##         -yell) + 0.5, "npc"), gp = gpar(fill = fill[1]))
    ##     grid.text(lab[1L], y = unit(1.5 + 0.5 * (lab[2L] != ""), 
    ##         "lines"))
    ##     if (lab[2L] != "") 
    ##         grid.text(lab[2L], y = unit(1, "lines"))
    ##     if (id) {
    ##         nodeIDvp <- viewport(x = unit(0.5, "npc"), y = unit(1, 
    ##             "npc"), width = max(unit(1, "lines"), unit(1.3, 
    ##             "strwidth", nam[id_node(node)])), height = max(unit(1, 
    ##             "lines"), unit(1.3, "strheight", nam[id_node(node)])))
    ##         pushViewport(nodeIDvp)
    ##         grid.rect(gp = gpar(fill = fill[2]))
    ##         grid.text(nam[id_node(node)])
    ##         popViewport()
    ##     }
    ##     upViewport(2)
    ## }
    return(rval)
}
class(my_node_inner2) <- "grapcon_generator"


## wrapper
my_tree_plot <- function(obj, terminal_panel = my_node_bp4,
                         inner_panel = my_node_inner2,
                         pop = TRUE,
                         tp_args = NULL,
                         ip_args = NULL,
                         main = NULL) {
    if(is.null(ip_args)) {
        ip_args <- list(pval = FALSE)
    } else {
        ip_args <- c(list(pval = FALSE), ip_args)
    }
        
    partykit:::plot.constparty(obj,
                               terminal_panel = terminal_panel,
                               inner_panel = inner_panel,
                               ip_args = ip_args,
                               tp_args = tp_args,
                               main = main,
                               pop = pop)  
}


## fitnessRank is numerical, but we want a histogram
## because a smoothed density makes little sense
## But the same could be argued about numObservedPeaks
density_color_3 <- function(data, cut,
                            varname = NULL,
                            nodenum = NULL,
                            rightc = "lightsalmon",
                            leftc = "#92C5DE") {
    if(varname %in% c("fitnessRank", "numObservedPeaks"))
       histogram_color_2(data, cut, varname,
                         nodenum)
    else 
       density_color_2(data, cut, varname, nodenum)

}

## https://stackoverflow.com/a/31216266/3670562
density_color_2 <- function(data, cut,
                            varname = NULL,
                            nodenum = NULL,
                            rightc = "lightsalmon",
                            leftc = "#92C5DE") {
    df.plot <- data.frame(dat = data)
    p <- ggplot(df.plot, aes(dat)) +
        geom_density(fill = leftc) ## +
        ## geom_vline(xintercept = cut)
    d <- ggplot_build(p)$data[[1]]
    p <- p + geom_area(data = subset(d, x > cut),
                       aes(x = x, y = y),
                       fill = rightc) ## +
        ## geom_segment(x = cut, xend = cut,
        ##              y = 0, yend = approx(x = d$x,
        ##                                   y = d$y,
        ##                                   xout = cut)$y,
        ##              colour = "black", size = 1)
    p <- p + theme(plot.title = element_text(hjust = 0.5)) +
        xlab(varname) +
        theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
              axis.title.x = element_text(color="black", size=20, face="bold")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste("Node ", nodenum)) +
        theme(plot.background = element_rect(fill = "white"))  +
        theme(axis.text.x = element_text(
                  size = 12,
                  angle= 45,
                  vjust = 1,
                  hjust = 1))
    ## Breaks the coloring logic and doesn't really help
    ## if(varname %in% c("propLocalMax", "observedProp", "fitnessRank") )
    ##     p <- p + scale_x_continuous(trans = "log1p")

    ## if(varname %in% c("propLocalMax", "observedProp", "fitnessRank") )
    ##      p <- p + scale_y_continuous(trans = "log1p")
    
    p
}


## Base on
## https://stackoverflow.com/questions/10016316/how-to-manually-fill-colors-in-a-ggplot2-histogram
histogram_color_2 <- function(data, cut,
                            varname = NULL,
                            nodenum = NULL,
                            rightc = "lightsalmon",
                            leftc = "#92C5DE") {

    dd <- sort(data)
    ## an ugly rightc colored line at the bottom through the range
    ## color <- factor(ifelse(dd < cut, 1, 2))
    ## df.plot <- data.frame(x = dd, color = color)
    ## p <- ggplot(df.plot, aes(x = x, fill = color, color = color)) +
    ##     geom_histogram(position = "identity")
    
    unique_dd <- unique(dd)
    color_code <- ifelse(unique_dd <= cut, 1, 2)
    colors <- c(leftc, rightc)[color_code]
    df.plot2 <- data.frame(x = dd)
    
    p <- ggplot(df.plot2, aes(x = x)) +
        geom_bar(fill = colors)
    
    p <- p + theme(plot.title = element_text(hjust = 0.5)) +
        xlab(varname) +
        theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
              axis.title.x = element_text(color="black", size=20, face="bold")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste("Node ", nodenum)) +
        theme(plot.background = element_rect(fill = "white"))  +
        theme(axis.text.x = element_text(
                  size = 12,
                  angle= 45,
                  vjust = 1,
                  hjust = 1)) 
    p
}



categorical_plot <- function(data, index, varname = NULL,
                             nodenum = NULL,
                             rightc = "lightsalmon",
                             leftc = "#92C5DE") {

    cat("\n inside categorical_plot")
    cat("\ n length data = ", length(data))
    ## rm the non-existent values: do below
    tdata <- table(data)
    ## tdata <- tdata[tdata > 0]
    df <- as.data.frame(tdata)
    ## cat("\n print df")
    ## print(df)
    ## cat("\n")
    colnames(df)[1] <- "Var1"
    df$index <- index
    ## Reorder levels of factors, just in case
    df <- df[order(df$index), ]
    df <- df[df$Freq > 0, ]
    df$Var1 <- factor(df$Var1, levels = df$Var1)
    if(varname == "typeLandscape_f") varname <- "typeLandscape"
    p <- ggplot(data = df, aes(x = Var1, y = Freq, fill = as.factor(index))) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values=c(leftc, rightc)) +
        theme(legend.position="none") +
        xlab(varname) +
        ggtitle(paste("Node ", nodenum)) +
        theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
              axis.title.x = element_text(color="black", size=20, face="bold")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(plot.background = element_rect(fill = "white")) +
        theme(axis.text.x = element_text(
                  size = 12,
                  angle= 45,
                  vjust = 1,
                  hjust = 1))
        ## theme(axis.title.x = element_blank())
    p
}



