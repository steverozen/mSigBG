#' Estimate a background signature and its characteristics.
#' 
#' @param bg.specta The spectra from which to compute the background information.
#' 
#' @param title The name of the single column of the signature in the
#'   \code{background.sig} element of the output.
#' 
#' @return A list with the elements \describe{
#'    \item{\code{background.sig}}{The background signatures as a
#'    single-column \code{\link[ICAMS]{ICAMS}} catalog.}
#'    \item{\code{sig.nbinom.size}}{The negative binomial \code{size} dispersion parameter for the 
#'    background signature profile. See \code{\link[stats]{NegBinomial}}.}
#'    \item{\code{count.nbinom.mu}}{The mean of the numbers of mutations in \code{bg.spectra}.}
#'    \item{\code{count.nbinom.size}}{The negative binomial \code{size} dispersion parameter for the 
#'    numbers of mutations caused by the background signature (i.e. \code{count.nbinom.mu}).
#'    See \code{\link[stats]{NegBinomial}}.}
#'    \item{\code{input.spectra}}{The \code{bg.spectra} used to estimate the background.}
#' }
#' 
#' @export
#' 
MakeBackgroundInfo <- function(bg.spectra, title = "Background.sig") {
  return(list(background.sig    = MeanOfSpectraAsSig(bg.spectra, title),
              sig.nbinom.size   = 10,
              count.nbinom.mu   = mean(colSums(bg.spectra)), 
              count.nbinom.size = 10,
              input.spectra     = bg.spectra))
}
