#' Subtract the maximum attributable mutations of a known signature from a spectrum.
#'
#' @description
#' Given a single observed mutational spectrum and a known signature to subtract,
#' find the maximum number of mutations (\code{n_subtract}) attributable to the
#' known signature such that a Poisson-approximated multinomial draw of
#' \code{n_subtract} mutations from \code{sig_to_subtract} has at least
#' \code{target_prob} probability of not exceeding the observed spectrum count
#' in any channel. In other words, for every channel \eqn{j}, the draw from
#' \code{sig_to_subtract} should plausibly fit within the observed \code{spectrum[j]}.
#'
#' The residual signature profile is computed by subtracting
#' \code{n_subtract * sig_to_subtract} from the spectrum, clamping negative
#' values to zero, and normalizing.
#'
#' @param spectrum A single observed spectrum as a numeric vector or a
#'   single-column matrix (e.g., an \code{\link[ICAMS]{ICAMS}} catalog).
#'
#' @param sig_to_subtract A signature profile (non-negative numeric vector
#'   summing to 1) representing the signature to subtract. Can be a
#'   single-column matrix or a plain numeric vector of length
#'   matching \code{spectrum}.
#'
#' @param target_prob The target probability that a multinomial draw of
#'   \code{n_subtract} mutations from \code{sig_to_subtract} does not
#'   exceed the observed count in any channel. If \code{NULL}, defaults
#'   to 0.05 for spectra with <= 100 channels or 0.01 for longer spectra.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{n_subtract}}{Estimated maximum number of mutations
#'     attributable to \code{sig_to_subtract}.}
#'   \item{\code{n_residual}}{Remaining mutations (\code{N - n_subtract}).}
#'   \item{\code{residual_sig}}{Residual signature profile (non-negative
#'     numeric vector summing to 1).}
#'   \item{\code{prob_no_exceed}}{The actual P(no exceed) at the
#'     estimated \code{n_subtract}.}
#'   \item{\code{target_prob}}{The target probability used.}
#' }
#'
#' @export
max_subtract_signature <- function(spectrum,
                                   sig_to_subtract,
                                   target_prob = NULL) {

  # Extract numeric vectors from matrix inputs if needed
  if (is.matrix(spectrum)) {
    stopifnot(ncol(spectrum) == 1)
    channel_names <- rownames(spectrum)
    spectrum_vec <- as.numeric(spectrum[, 1])
  } else {
    channel_names <- names(spectrum)
    spectrum_vec <- as.numeric(spectrum)
  }

  if (is.matrix(sig_to_subtract)) {
    stopifnot(ncol(sig_to_subtract) == 1)
    sig_vec <- as.numeric(sig_to_subtract[, 1])
  } else {
    sig_vec <- as.numeric(sig_to_subtract)
  }

  # Validate inputs
  stopifnot(length(spectrum_vec) == length(sig_vec))
  stopifnot(all(spectrum_vec >= 0))
  stopifnot(all(sig_vec >= 0))
  stopifnot(abs(sum(sig_vec) - 1) < 1e-6)

  N <- sum(spectrum_vec)
  n_channels <- length(spectrum_vec)

  if (is.null(target_prob)) {
    target_prob <- if (n_channels <= 100) 0.05 else 0.01
  }

  # Bisect to find n_subtract where P(no exceed) = target_prob
  n_subtract <- bisect_n_subtract(sig_vec, spectrum_vec, target_prob)

  # Compute residual signature
  residual_counts <- spectrum_vec - n_subtract * sig_vec
  residual_counts <- pmax(residual_counts, 0)
  sum_residual <- sum(residual_counts)

  if (sum_residual > 0) {
    residual_sig <- residual_counts / sum_residual
  } else {
    residual_sig <- rep(1 / n_channels, n_channels)
  }

  if (!is.null(channel_names)) {
    names(residual_sig) <- channel_names
  }

  prob_actual <- prob_no_exceed_poisson(n_subtract, sig_vec, spectrum_vec)

  list(
    n_subtract     = n_subtract,
    n_residual     = N - n_subtract,
    residual_sig   = residual_sig,
    prob_no_exceed = prob_actual,
    target_prob    = target_prob
  )
}


#' Probability that a Poisson-approximated multinomial draw does not exceed
#' observed counts in any channel.
#'
#' Approximates each channel of a \code{Multinomial(n, sig)} draw as
#' independent \code{Poisson(n * sig_j)}. This is conservative because
#' multinomial channels are negatively correlated.
#'
#' @param n_subtract Number of mutations drawn.
#' @param sig Signature profile (numeric vector summing to 1).
#' @param spectrum Observed counts per channel (numeric vector).
#'
#' @return Probability (scalar between 0 and 1).
#'
#' @keywords internal
prob_no_exceed_poisson <- function(n_subtract, sig, spectrum) {
  lambdas <- n_subtract * sig
  exp(sum(stats::ppois(spectrum, lambda = lambdas, log.p = TRUE)))
}


#' Bisect to find n_subtract where P(no exceed) equals target_prob.
#'
#' @param sig Signature profile (numeric vector summing to 1).
#' @param spectrum Observed counts per channel (numeric vector).
#' @param target_prob Target probability.
#' @param tol Bisection tolerance.
#'
#' @return Estimated n_subtract (numeric scalar).
#'
#' @keywords internal
bisect_n_subtract <- function(sig, spectrum, target_prob, tol = 0.5) {
  lo <- 0
  hi <- sum(spectrum)

  if (prob_no_exceed_poisson(0, sig, spectrum) <= target_prob) return(0)
  if (prob_no_exceed_poisson(hi, sig, spectrum) >= target_prob) return(hi)

  while (hi - lo > tol) {
    mid <- (lo + hi) / 2
    if (prob_no_exceed_poisson(mid, sig, spectrum) > target_prob) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  (lo + hi) / 2
}
