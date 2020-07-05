#' Calculate the inferred target signature and inferred background and target components of input spectra
#' 
#' 
#' @inheritParams SeparateSignatureFromBackground
#' 
#' @param sig.name A name for the inferred signature
#' 
#' @export
#' 
#' @return A list with the elements\describe{
#' 
#' \item{\code{inferred.target.spectra}}{The mutations counts inferred
#'   to be from the target signature.}
#' 
#' }

SeparateSignatureAndSpectra <- function(
  spectra,
  bg.sig.info,
  m.opts = NULL,
  start.b.fraction = 0.1,
  sig.name = "Inferred.sig") {
  
  ret <- SeparateSignatureFromBackground(
    spectra          = spectra,
    bg.sig.info      = bg.sig.info,
    m.opts           = m.opts,
    start.b.fraction = start.b.fraction)
  
  inferred.sig <- 
    ICAMS::as.catalog(
      object         = ret$inferred.target.sig, 
      catalog.type   = "counts.signature",
      region         = "genome",
      abundance      = attr(spectra, "abundance"),
      infer.rownames = TRUE)
  
  colnames(inferred.sig) <- sig.name
  
  ret$inferred.target.sig.as.catalog <- inferred.sig
  
  inferred.bg.spectra <- 
    round(bg.sig.info$background.sig %*% 
            matrix(ret$exposures.to.bg.sig, nrow = 1))
  
  inferred.bg.spectra <- ICAMS::as.catalog(inferred.bg.spectra,
                                           catalog.type   = "counts",
                                           region         = "genome",
                                           abundance      = attr(spectra, "abundance"),
                                           infer.rownames = TRUE)
  colnames(inferred.bg.spectra) <- paste0(colnames(spectra), "-inf.bg.spect")
  
  ret$inferred.bg.spectra <- inferred.bg.spectra
  
  ret$inferred.target.spectra <- spectra - inferred.bg.spectra
  
  # TODO: try an alternative estimation, in which each peak, p, in the specturm is
  # allocated to according to the ratio
  # (total.background.count) * p_(in background sig) : (total.target.count) * p_(in target sig)
  
  colnames(ret$inferred.target.spectra) <- paste0(colnames(spectra),".inf.tg.spect")
  
  return(ret)
  
}