#' Estimate a background signature and its characteristics.
#' 
#' @param bg.specta The spectra from which to compute the background information.
#' 
#' @param title The name of the single column of the signature in the
#'   \code{background.sig} element of the output.
#' 
#' @return A list with the elements \describe{
#'    \item{\code{background.sig}}{The background signatures as an \code{\link[ICAMS]{ICAMS}} catalog.}
#'    \item{\code{sig.nbinom.size}}{The negative binomial \code{size} dispersion parameter for the 
#'    background signature profile. See \code{\link[stats]{NegBinomial}}.}
#'    \item{\code{count.nbinom.mu}}{XXXXXX}
#'    \item{\code{count.nbinom.size}}{XXXXXX}
#' }
#' 
#' 
#' An single column .
#' 
#' @export
#' 
MakeBackgroundInfo <- function(bg.spectra, title = "Background.sig") {
  
  return(list(background.sig    = MeanOfSpectraAsSig(bg.spectra, title),
              sig.nbinom.size   = 10,
              count.nbinom.mu   = mean(colSums(bg.spectra)), 
              count.nbinom.size = 20))
}
