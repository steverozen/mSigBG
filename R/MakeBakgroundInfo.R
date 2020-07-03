#' Estimate a background signature.
#' 
#' @param bg.specta The spectra from which to compute the background information.
#' 
#' @param title The name of the single column in the output.
#' 
#' @return An single column \code{\link[ICAMS]{ICAMS}} catalog.
#' 
#' @export
#' 
MakeBackgroundInfo <- function(bg.spectra, title = "Background.sig") {
  
  return(list(background.sig    = MeanOfSpectraAsSig(bg.spectra, title),
              sig.nbinom.size   = 10,
              count.nbinom.mu   = mean(colSums(bg.spectra)), 
              count.nbinom.size = 20))
}
