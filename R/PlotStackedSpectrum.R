#' Plot a spectrum as a stacked bar chart
#' 
#' @param bg.cat.1 XXXX
#' 
#' @param target.cat XXXX
#' 
#' @param background.tite XXXX
#' 
#' @param target.title XXXXX
#' 
#' @param set.neg.zero XXXX
#' 
#' @export

Plot1StackedSpectrum <- function(bg.cat.1, 
                                   target.cat.1,
                                   background.title = "Background",
                                   target.title = "Target", 
                                   set.neg.zero = TRUE) {
  
  zeros <- which(target.cat.1 < 0)
  neg.catalog <- target.cat.1
  
  if (set.neg.zero) {
    target.cat.1[zeros] <- 0
  }
  neg.catalog[-zeros] <- 0
  neg.catalog <- neg.catalog * -1
  
  to.plot <- cbind(bg.cat.1, target.cat.1)
  
  bp <- ICAMS::PlotCatalog(to.plot)$plot.object
  
  bpz <- bp[zeros]
  # message(paste(bpz, collapse = " "))
  
  y0 <- bg.cat.1[zeros]
  y1 <- bg.cat.1[zeros] - neg.catalog[zeros]
  if (attr(target.cat.1, "catalog.type") == "density") {
    # message("adjusting for density")
    y0 <- y0 * 1e6
    y1 <- y1 * 1e6
  }
  
  # Legend
  minus.col <- "purple"
  segments(x0 = bpz, 
           y0 = y0, 
           y1 = y1, 
           col = minus.col,
           lwd = 3,
           lty = 1)
  
  legend.txt <- c(target.title, background.title)
  legend.col <- c("grey35", "red")
  if (length(bpz) > 0) {
    legend.txt <- c(legend.txt, "Background minus target")
    legend.col <- c(legend.col, minus.col)
  }
  legend("topright", legend = legend.txt, fill = legend.col, cex = 0.8)

  return(invisible(bp))
}




Not.used.PlotStackedSpectrum.R <- function(bg.catalog, target.catalog, set.neg.zero = TRUE) {
  
  xlabels <- FALSE
  
  for (i in 1:ncol(bg.catalog)) {
    if (FALSE && i == ncol(bg.catalog)) { # FALSE an experiment
      xlabels <- TRUE
      par(mar = c(last.bottom.mar, 4, 3, 2))
    }
    
    target.cat.1 <- target.catalog[ , i, drop = FALSE]
    bg.cat.1     <- bg.catalog[ , i, drop = FALSE]    
    
    zeros <- which(target.cat.1 < 0)
    neg.catalog <- target.cat.1
    
    if (set.neg.zero) {
      target.cat.1[zeros] <- 0
    }
    neg.catalog[-zeros] <- 0
    neg.catalog <- neg.catalog * -1
    
    to.plot <- cbind(bg.cat.1, target.cat.1)
    
    bp <- ICAMS::PlotCatalog(
      to.plot,
      upper   = (i == 1),
      xlabels = xlabels)$plot.object
    
    bpz <- bp[zeros]
    # message(paste(bpz, collapse = " "))
    
    y0 <- bg.cat.1[zeros]
    y1 <- bg.cat.1[zeros] - neg.catalog[zeros]
    if (attr(target.catalog, "catalog.type") == "density") {
      # message("adjusting for density")
      y0 <- y0 * 1e6
      y1 <- y1 * 1e6
    }
    
    # message(paste(y0, collapse = "y"))
    # message(paste(y1, collapse = "z"))
    minus.col <- "purple"
    segments(x0 = bpz, 
             y0 = y0, 
             y1 = y1, 
             col = minus.col,
             lwd = 3,
             lty = 1)
    
    legend.txt <- c(params$whichnitro, "Background")
    legend.col <- c("grey35", "red")
    if (length(bpz) > 0) {
      legend.txt <- c(legend.txt, paste("Background -", params$whichnitro))
      legend.col <- c(legend.col, minus.col)
    }
    legend("topright", legend = legend.txt, fill = legend.col, cex = 0.8)
  }
}

