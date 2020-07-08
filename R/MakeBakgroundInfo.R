#' Estimate a background signature and its characteristics.
#' 
#' @param bg.spectra The spectra from which to compute the background information.
#' 
#' @param title The name of the single column of the signature in the
#'   \code{background.sig} element of the output.
#'   
#' @param sig.nbinom.size The negative binomial \code{size} dispersion parameter for the 
#'    background signature profile. See \code{\link[stats]{NegBinomial}}.
#'    Smaller is more dispersed.
#' 
#' @param count.nbinom.size The negative binomial \code{size} dispersion parameter for the 
#'    numbers of mutations caused by the background signature 
#'    (i.e. \code{count.nbinom.mu}).
#'    See \code{\link[stats]{NegBinomial}}.
#'    Smaller is more dispersed.
#' 
#' @return A list with the elements \describe{
#'    \item{\code{background.sig}}{The background signatures as a
#'    single-column \code{\link[ICAMS]{ICAMS}} catalog.}
#'    \item{\code{sig.nbinom.size}}{See input argument \code{sig.nbinom.size}.}
#'    \item{\code{count.nbinom.mu}}{The mean of the numbers of mutations in \code{bg.spectra}.}
#'    \item{\code{count.nbinom.size}}{See input argument \code{count.nbinom.size}.}
#'    \item{\code{input.spectra}}{The \code{bg.spectra} used to estimate the background.}
#' }
#' 
#' @export
#' 
MakeBackgroundInfo <- function(bg.spectra, 
                               title             = "Background.sig",
                               sig.nbinom.size   = 10,
                               count.nbinom.size = 10) {
  return(list(background.sig    = MeanOfSpectraAsSig(bg.spectra, title),
              sig.nbinom.size   = sig.nbinom.size,
              count.nbinom.mu   = mean(colSums(bg.spectra)), 
              count.nbinom.size = count.nbinom.size,
              input.spectra     = bg.spectra))
}
