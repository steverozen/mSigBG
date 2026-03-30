# Accuracy tests for max_subtract_signature
# Runs 20 random trials each for 89, 96, and 476 features.
# Reports how well n_subtract recovers the true value.

devtools::load_all()

run_accuracy_trials <- function(n_features, n_trials = 20) {
  results <- data.frame(
    trial      = integer(n_trials),
    true_n_a   = numeric(n_trials),
    true_n_b   = numeric(n_trials),
    n_subtract = numeric(n_trials),
    ratio      = numeric(n_trials),
    cos_sim    = numeric(n_trials)
  )

  for (i in seq_len(n_trials)) {
    sig_a <- runif(n_features); sig_a <- sig_a / sum(sig_a)
    sig_b <- runif(n_features); sig_b <- sig_b / sum(sig_b)

    # Random mixture: n_a is 10-40% of total, total is 1000-10000
    total <- sample(1000:10000, 1)
    frac_a <- runif(1, 0.1, 0.4)
    n_a <- round(total * frac_a)
    n_b <- total - n_a

    spec_a <- as.numeric(rmultinom(1, size = n_a, prob = sig_a))
    spec_b <- as.numeric(rmultinom(1, size = n_b, prob = sig_b))
    spectrum <- spec_a + spec_b

    result <- max_subtract_signature(spectrum, sig_a, target_prob = 0.01)

    cos_sim <- lsa::cosine(result$residual_sig, sig_b)[1, 1]

    results$trial[i]      <- i
    results$true_n_a[i]   <- n_a
    results$true_n_b[i]   <- n_b
    results$n_subtract[i] <- round(result$n_subtract, 1)
    results$ratio[i]      <- round(result$n_subtract / n_a, 3)
    results$cos_sim[i]    <- round(cos_sim, 4)
  }

  results
}

set.seed(2024)

for (nf in c(89, 96, 476)) {
  cat(strrep("=", 70), "\n")
  cat(sprintf("Features: %d  (20 trials)\n", nf))
  cat(strrep("=", 70), "\n")

  res <- run_accuracy_trials(nf)
  print(res, row.names = FALSE)

  cat("\nSummary:\n")
  cat(sprintf("  n_subtract / true_n_a  —  mean: %.3f  median: %.3f  range: [%.3f, %.3f]\n",
              mean(res$ratio), median(res$ratio),
              min(res$ratio), max(res$ratio)))
  cat(sprintf("  cosine similarity      —  mean: %.4f  median: %.4f  range: [%.4f, %.4f]\n",
              mean(res$cos_sim), median(res$cos_sim),
              min(res$cos_sim), max(res$cos_sim)))
  cat("\n")
}
