#' Subtract the maximum attributable mutations of a known signature from a spectrum.
#'
#' @description
#' Given a single observed mutational spectrum and a known signature to subtract,
#' find the maximum number of mutations (\code{n_subtract}) attributable to the
#' known signature such that a Poisson-approximated multinomial draw of
#' \code{n_subtract} mutations from \code{sig_to_subtract} plausibly fits
#' within the observed spectrum. Specifically, for every channel \eqn{j},
#' \code{P(Poisson(n_subtract * sig_j) <= spectrum_j)} must meet a
#' per-channel threshold derived from \code{target_prob} using the
#' \enc{Šidák}{Sidak} correction, which adjusts for the number of channels
#' so that the result is approximately invariant to spectrum dimensionality
#' (e.g., SBS96 vs SBS1536).
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
#' @param target_prob The target joint probability that a Poisson-approximated
#'   multinomial draw of \code{n_subtract} mutations from
#'   \code{sig_to_subtract} does not exceed the observed count in any channel.
#'   A per-channel threshold is derived via the \enc{Šidák}{Sidak} correction:
#'   \code{target_prob^(1/n_channels)}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{n_subtract}}{Estimated maximum number of mutations
#'     attributable to \code{sig_to_subtract}.}
#'   \item{\code{n_residual}}{Remaining mutations (\code{N - n_subtract}).}
#'   \item{\code{residual_sig}}{Residual signature profile (non-negative
#'     numeric vector summing to 1).}
#'   \item{\code{min_channel_prob}}{The minimum per-channel P(not exceed) at
#'     the estimated \code{n_subtract}.}
#'   \item{\code{per_channel_threshold}}{The \enc{Šidák}{Sidak}-corrected
#'     per-channel threshold used.}
#'   \item{\code{target_prob}}{The target joint probability used.}
#' }
#'
#' @export
old_max_subtract_signature <- function(spectrum,
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
    target_prob <- 0.05
  }

  # Šidák correction: per-channel threshold from desired joint probability
  per_channel_threshold <- target_prob^(1 / n_channels)

  # Bisect to find n_subtract where min per-channel P(not exceed) = threshold
  n_subtract <- bisect_n_subtract(sig_vec, spectrum_vec, per_channel_threshold)

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

  min_ch_prob <- min_channel_prob(n_subtract, sig_vec, spectrum_vec)

  list(
    n_subtract            = n_subtract,
    n_residual            = N - n_subtract,
    residual_sig          = residual_sig,
    min_channel_prob      = min_ch_prob,
    per_channel_threshold = per_channel_threshold,
    target_prob           = target_prob
  )
}


#' Minimum per-channel probability that a Poisson draw does not exceed the
#' observed count.
#'
#' For each channel \eqn{j}, computes
#' \code{P(Poisson(n_subtract * sig_j) <= spectrum_j)} and returns the
#' minimum across channels. Used with a \enc{Šidák}{Sidak}-corrected
#' per-channel threshold so the criterion is approximately invariant to the
#' number of channels.
#'
#' @param n_subtract Number of mutations drawn.
#' @param sig Signature profile (numeric vector summing to 1).
#' @param spectrum Observed counts per channel (numeric vector).
#'
#' @return Minimum per-channel probability (scalar between 0 and 1).
#'
#' @keywords internal
min_channel_prob <- function(n_subtract, sig, spectrum) {
  lambdas <- n_subtract * sig
  min(stats::ppois(spectrum, lambda = lambdas))
}


#' Bisect to find n_subtract where the minimum per-channel probability
#' equals the per-channel threshold.
#'
#' @param sig Signature profile (numeric vector summing to 1).
#' @param spectrum Observed counts per channel (numeric vector).
#' @param per_channel_threshold \enc{Šidák}{Sidak}-corrected per-channel
#'   probability threshold.
#' @param tol Bisection tolerance.
#'
#' @return Estimated n_subtract (numeric scalar).
#'
#' @keywords internal
bisect_n_subtract <- function(sig, spectrum, per_channel_threshold, tol = 0.5) {
  lo <- 0
  hi <- sum(spectrum)

  if (min_channel_prob(0, sig, spectrum) <= per_channel_threshold) return(0)
  if (min_channel_prob(hi, sig, spectrum) >= per_channel_threshold) return(hi)

  while (hi - lo > tol) {
    mid <- (lo + hi) / 2
    if (min_channel_prob(mid, sig, spectrum) > per_channel_threshold) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  (lo + hi) / 2
}
