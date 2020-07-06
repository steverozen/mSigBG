#' Convert spectra to signatures and then plot mean with "error" bars
#' 
#' @inheritParams  MeanOfSpectraAsSig
#' 
#' @export
#'
#' @return The mean of the spectra as a signature, the
#'   constituent spectra as signatures, and the y positions of the
#'   arrowheads.


PlotSpectraAsSigsWithUncertainty <- function(spectra, title = "Mean.as.signature") {
  xx <- MeanOfSpectraAsSig(spectra = spectra, title = title)
  arrow.tops <- apply(xx$constituent.sigs, 1, max)
  arrow.bottoms <- apply(xx$constituent.sigs, 1, min)
  
  options(digits = 2)
  
  bp <- ICAMS::PlotCatalog(
    catalog = xx$mean.sig[ , 1, drop = FALSE],
    upper   = TRUE,
    xlabels = TRUE,
    cex     = 0.8,
    ylim    = c(min(arrow.bottoms), max(arrow.tops) + 0.005)
  )
  AddArrows(bp$plot.object, arrow.tops, arrow.bottoms)
  xx$arrow.tops    <- arrow.tops
  xx$arrow.bottoms <- arrow.bottoms
  return(invisible(xx))
}

AddArrows <- function(bp, tops, bottoms) {
  oldopt <- getOption("warn")
  on.exit(options(warn = oldopt))
  options(warn = -1) # Does not seem to turn off warnings
  which0 <- which((tops - bottoms) == 0)
  tops[which0] <- tops[which0] + 1e-5 
  suppressWarnings(
    # Necessary because generates warnings for 0-length arrows
    arrows(
      x0     = bp,
      y0     = tops,    # location of up arrow tips
      y1     = bottoms, # location of down arrow tips
      angle  = 90,      # use "T" arrow heads
      code   = 3,       # use arrow heads at both ends
      length = 0.025    # width of arrow heads
    ))
  # Try this to get rid of the warnings:
  assign("last.warning", NULL, envir = baseenv())
  
}
