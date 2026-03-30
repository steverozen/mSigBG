#' Subtract the maximum attributable mutations with a limit on total negative
#' residual counts.
#'
#' @description
#' Given a single observed mutational spectrum and a known signature to subtract,
#' find the maximum number of mutations (\code{n_subtract}) such that the
#' sum of absolute values of negative residuals does not exceed
#' \code{max_neg_fraction * sum(spectrum)}.
#'
#' The residual signature profile is computed by clamping negative values to
#' zero and normalizing.
#'
#' @param spectrum A single observed spectrum as a numeric vector or a
#'   single-column matrix (e.g., an \code{\link[ICAMS]{ICAMS}} catalog).
#'
#' @param sig_to_subtract A signature profile (non-negative numeric vector
#'   summing to 1) representing the signature to subtract. Can be a
#'   single-column matrix or a plain numeric vector of length
#'   matching \code{spectrum}.
#'
#' @param max_neg_fraction Maximum total negative residual count as a fraction
#'   of \code{sum(spectrum)}.
#'
#' @param nbinom.size Dispersion parameter for the negative binomial
#'   distribution used in the Monte Carlo p-value simulation; smaller
#'   values mean more overdispersion.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{n_subtract}}{Maximum number of mutations subtractable
#'     within the negative-count limit.}
#'   \item{\code{n_residual}}{Remaining mutations (\code{N - n_subtract}).}
#'   \item{\code{residual_sig}}{Residual signature profile (non-negative
#'     numeric vector summing to 1).}
#'   \item{\code{total_negative}}{Sum of absolute values of negative residuals
#'     at the estimated \code{n_subtract}.}
#'   \item{\code{n_negative_channels}}{Number of channels with negative
#'     residuals at the estimated \code{n_subtract}.}
#'   \item{\code{prob_ge_total_negative}}{Monte Carlo estimated probability
#'     that negative-binomial-sampled spectra (using the estimated
#'     \code{n_subtract}, \code{residual_sig}, and \code{sig_to_subtract})
#'     would produce a total negative count >= the observed
#'     \code{total_negative}.}
#' }
#'
#' @export
max_subtract_signature <- function(spectrum, sig_to_subtract,
                          max_neg_fraction = 0.05,
                          nbinom.size = 10) {

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
  max_neg_count <- max_neg_fraction * N

  # Bisect to find max n_subtract where total negative <= max_neg_count
  lo <- 0
  hi <- N
  tol <- 0.5

  total_neg <- function(n) {
    residual <- spectrum_vec - n * sig_vec
    sum(abs(residual[residual < 0]))
  }

  if (total_neg(hi) <= max_neg_count) {
    n_subtract <- hi
  } else {
    while (hi - lo > tol) {
      mid <- (lo + hi) / 2
      if (total_neg(mid) <= max_neg_count) {
        lo <- mid
      } else {
        hi <- mid
      }
    }
    n_subtract <- lo
  }

  # Compute residual signature
  residual_counts <- spectrum_vec - n_subtract * sig_vec
  neg_mask <- residual_counts < 0
  total_negative <- sum(abs(residual_counts[neg_mask]))
  n_negative_channels <- sum(neg_mask)
  residual_counts <- pmax(residual_counts, 0)
  sum_residual <- sum(residual_counts)

  n_channels <- length(spectrum_vec)
  if (sum_residual > 0) {
    residual_sig <- residual_counts / sum_residual
  } else {
    residual_sig <- rep(1 / n_channels, n_channels)
  }

  if (!is.null(channel_names)) {
    names(residual_sig) <- channel_names
  }

  # Estimate probability of getting >= total_negative negative mutations
  # by Monte Carlo simulation. Under the model, the observed spectrum is
  # n_subtract draws from sig_to_subtract + n_residual draws from residual_sig.
  # Simulate negative binomial draws per channel for each component, subtract
  # the expected sig_to_subtract contribution, and measure total negative.
  n_residual <- N - n_subtract
  n_sim <- 1000
  expected_subtract <- n_subtract * sig_vec
  sim_neg_totals <- vapply(seq_len(n_sim), function(i) {
    sim_subtract <- stats::rnbinom(n_channels, mu = expected_subtract,
                                   size = nbinom.size)
    sim_residual <- stats::rnbinom(n_channels, mu = n_residual * residual_sig,
                                   size = nbinom.size)
    sim_spectrum <- sim_subtract + sim_residual
    sim_resid <- sim_spectrum - expected_subtract
    sum(abs(sim_resid[sim_resid < 0]))
  }, numeric(1))
  prob_ge_total_negative <- mean(sim_neg_totals >= total_negative)

  list(
    n_subtract             = n_subtract,
    n_residual             = n_residual,
    residual_sig           = residual_sig,
    total_negative         = total_negative,
    n_negative_channels    = n_negative_channels,
    prob_ge_total_negative = prob_ge_total_negative
  )
}
