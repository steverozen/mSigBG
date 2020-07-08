#' Specifications of background signatures
#'
#'
#' @format A list with one element for each cell line. 
#' Each element of the list
#' is in turn a list with the elements \describe{
#' 
#'    \item{\code{background.sig}}{The background signatures as a
#'    single-column \code{\link[ICAMS]{ICAMS}} catalog.}
#'  
#'    \item{\code{sig.nbinom.size}}{The negative binomial \code{size} dispersion parameter for the 
#'    background signature profile. See \code{\link[stats]{NegBinomial}}.}
#'    
#'    \item{\code{count.nbinom.mu}}{The mean of the numbers of mutations in \code{input.spectra}.}
#'   
#'    \item{\code{count.nbinom.size}}{The negative binomial \code{size} dispersion parameter for the 
#'    numbers of mutations caused by the background signature (i.e. \code{count.nbinom.mu}).
#'    See \code{\link[stats]{NegBinomial}}.}
#'    
#'    \item{\code{input.spectra}}{The \code{bg.spectra} used to estimate the background.}
#' }
#' 
#' @source \code{background.info} 
#' was estimated from \code{\link{HepG2.background.spectra}} and
#' \code{\link{MCF10A.background.spectra}}.
#'
#' @name background.info
#' 
#' @examples 
#' background.info[["HepG2"]]$count.nbinom.mu
#' background.info[["HepG2"]]$count.nbinom.size
#' background.info[["HepG2"]]$sig.nbinom.size
#' background.info[["HepG2"]]$background.sig[1:3, ]
#' 
#' 
#' 
"background.info"


#' Example spectra of cell lines exposed to cisplatin.
#' 
#' @format A list of \code{\link[ICAMS]{ICAMS}} \code{counts} spectra catalogs.
#'
#' @name example.spectra
#' 
#' @examples 
#' rowSums(example.spectra[["MCF10A.cisplatin"]])[1:3]
"example.spectra"


#' Background spectra for the HepG2 and MCF-10A cell lines.
#' 
#' @format An \code{\link[ICAMS]{ICAMS}} \code{counts} catalog.
#' 
#' @name background.spectra
#' 
"HepG2.background.spectra"

 
#' @name background.spectra
#' 
"MCF10A.background.spectra"
