#' Return the mean of multiple spectra as a signature.
#' 
#' @param spectra An \code{\link[ICAMS]{ICAMS}} spectrum catalog.
#'   Convert each spectrum to a signature and then compute the
#'   mean.
#'   
#' @param title The name of the output signature.
#'   
#' @export
#' 
MeanOfSpectraAsSig <- function(spectra, title = "sig.from.spectra.mean") {
  
  if (is.null(spectra)) {
    stop("Argument spectra is NULL")
  }
  
  ctype <- attr(spectra, "catalog.type", exact = TRUE)
  
  if (is.null(ctype)) {
    stop("Cannot call MeanOfSpectraAsSig when catalog.type is NULL")
  } else if (ctype == "counts") {
    tctype <- "counts.signature"
  } else if (ctype == "density") {
    tctype <- "density.signature"
  } else {
    stop("Cannot call MeanOfSpectraAsSig when catalog.type is ", ctype)
  }
  
  target.abundance <- attr(spectra, "abundance", exact = TRUE)
  target.region    <- attr(spectra, "region",    exact = TRUE)

  # stopifnot(!is.null(target.abundance))
  if (is.null(target.abundance)) {
    warning("Using NULL target.abundance in MeanOfSpectraAsSig")
  }

  sigs <- ICAMS::TransformCatalog(spectra, 
                                  target.catalog.type = tctype,
                                  target.abundance = target.abundance)
  
  mean.sig <- apply(X = sigs, MARGIN = 1, mean)

  mean.sig <- matrix(mean.sig, ncol = 1)
  
  mean.sig <- 
    ICAMS::as.catalog(mean.sig, 
                      catalog.type   = tctype, 
                      region         = target.region, 
                      abundance      = target.abundance,
                      infer.rownames = TRUE)
  
  colnames(mean.sig) <- title
  
  return(list(mean.sig = mean.sig, constituent.sigs = sigs))
}
