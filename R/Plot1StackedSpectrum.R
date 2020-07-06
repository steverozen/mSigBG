#' Plot a spectrum as a stacked bar chart
#' 
#' @param background.spectrum Partial spectra due to a background signature.
#' 
#' @param target.spectrum Partial spectra due to a target signature.
#' 
#' @param background.title A title for the legend for the partial
#'   spectra due to a background signature.
#' 
#' @param target.title A title for the legend for the partial
#'   spectra due to a target signature.
#' 
#' @param set.neg.zero Sometimes after subtraction a part of the
#'   spectra due to a target signatures is negative. If this
#'   argument is true, set these to 0.
#' 
#' @export
#' 
#' @import graphics
#' 

Plot1StackedSpectrum <- function(background.spectrum, 
                                 target.spectrum,
                                 background.title = "Background",
                                 target.title     = "Target", 
                                 set.neg.zero     = TRUE) {
  
  zeros <- which(target.spectrum < 0)
  neg.catalog <- target.spectrum
  
  if (set.neg.zero) {
    target.spectrum[zeros] <- 0
  }
  neg.catalog[-zeros] <- 0
  neg.catalog <- neg.catalog * -1
  
  to.plot <- cbind(background.spectrum, target.spectrum)
  
  bp <- ICAMS::PlotCatalog(to.plot)$plot.object
  
  bpz <- bp[zeros]
  # message(paste(bpz, collapse = " "))
  
  y0 <- background.spectrum[zeros]
  y1 <- background.spectrum[zeros] - neg.catalog[zeros]
  if (attr(target.spectrum, "catalog.type") == "density") {
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

