
test_that("max_subtract_signature basic properties with synthetic data", {
  set.seed(42)

  sig_a <- runif(89); sig_a <- sig_a / sum(sig_a)
  sig_b <- runif(89); sig_b <- sig_b / sum(sig_b)

  n_a <- 200
  n_b <- 800
  spec_a <- as.numeric(rmultinom(1, size = n_a, prob = sig_a))
  spec_b <- as.numeric(rmultinom(1, size = n_b, prob = sig_b))
  spectrum <- spec_a + spec_b

  result <- max_subtract_signature(spectrum, sig_a, max_neg_fraction = 0.02)

  # residual_sig is a valid probability distribution
  expect_equal(sum(result$residual_sig), 1, tolerance = 1e-10)
  expect_true(all(result$residual_sig >= 0))

  # n_subtract + n_residual = N
  expect_equal(result$n_subtract + result$n_residual, sum(spectrum))

  # total_negative within the limit
  expect_lte(result$total_negative, 0.02 * sum(spectrum) + 1)

  # n_subtract is non-negative and bounded
  expect_gte(result$n_subtract, 0)
  expect_lte(result$n_subtract, sum(spectrum))
})


test_that("max_subtract_signature works with matrix input", {
  set.seed(123)

  sig_a <- runif(89); sig_a <- sig_a / sum(sig_a)
  sig_b <- runif(89); sig_b <- sig_b / sum(sig_b)
  spectrum <- matrix(
    as.numeric(rmultinom(1, size = 500, prob = sig_a)) +
    as.numeric(rmultinom(1, size = 500, prob = sig_b)),
    ncol = 1
  )
  sig_mat <- matrix(sig_a, ncol = 1)

  result <- max_subtract_signature(spectrum, sig_mat)

  expect_equal(sum(result$residual_sig), 1, tolerance = 1e-10)
  expect_true(all(result$residual_sig >= 0))
  expect_gte(result$n_subtract, 0)
  expect_lte(result$n_subtract, sum(spectrum))
})


test_that("max_subtract_signature returns near-zero when spectrum has no contribution", {
  set.seed(77)

  # sig_a has all weight in channels 1-44, sig_b in channels 45-89
  sig_a <- c(runif(44), rep(0, 45)); sig_a <- sig_a / sum(sig_a)
  sig_b <- c(rep(0, 44), runif(45)); sig_b <- sig_b / sum(sig_b)

  # Spectrum is purely from sig_b
  spectrum <- as.numeric(rmultinom(1, size = 1000, prob = sig_b))

  result <- max_subtract_signature(spectrum, sig_a)

  # Should subtract very few mutations since sig_a channels are empty
  expect_lt(result$n_subtract, 50)
})


test_that("max_subtract_signature realistic 89-channel at max_neg_fraction=0.02", {
  fixture <- readRDS(test_path("fixtures", "test_89_realistic.rds"))

  set.seed(999)
  result <- max_subtract_signature(fixture$spectrum, fixture$sig_subtract,
                                   max_neg_fraction = 0.02)

  # n_subtract should be close to true (ratio ~1.01 in our trials)
  ratio <- result$n_subtract / fixture$n_target
  expect_gt(ratio, 0.85)
  expect_lt(ratio, 1.15)

  # Residual should resemble true background
  cos_sim <- lsa::cosine(result$residual_sig, fixture$true_bg_sig)[1, 1]
  expect_gt(cos_sim, 0.9)

  # total_negative within limit
  expect_lte(result$total_negative, 0.02 * sum(fixture$spectrum) + 1)

  # Valid probability distribution
  expect_equal(sum(result$residual_sig), 1, tolerance = 1e-10)
  expect_true(all(result$residual_sig >= 0))

  # Monte Carlo p-value should not be extreme
  expect_gt(result$prob_ge_total_negative, 0.01)
})


test_that("max_subtract_signature realistic 89-channel at max_neg_fraction=0.01", {
  fixture <- readRDS(test_path("fixtures", "test_89_realistic.rds"))

  set.seed(999)
  result <- max_subtract_signature(fixture$spectrum, fixture$sig_subtract,
                                   max_neg_fraction = 0.01)

  # At 0.01, method is more conservative — ratio ~0.96 in our trials
  ratio <- result$n_subtract / fixture$n_target
  expect_gt(ratio, 0.80)
  expect_lt(ratio, 1.10)

  # Residual should still resemble true background
  cos_sim <- lsa::cosine(result$residual_sig, fixture$true_bg_sig)[1, 1]
  expect_gt(cos_sim, 0.85)

  # total_negative within limit
  expect_lte(result$total_negative, 0.01 * sum(fixture$spectrum) + 1)

  # Valid probability distribution
  expect_equal(sum(result$residual_sig), 1, tolerance = 1e-10)
  expect_true(all(result$residual_sig >= 0))

  # Monte Carlo p-value should not be extreme
  expect_gt(result$prob_ge_total_negative, 0.01)
})


test_that("stricter max_neg_fraction subtracts fewer mutations", {
  fixture <- readRDS(test_path("fixtures", "test_89_realistic.rds"))

  set.seed(999)
  r_loose <- max_subtract_signature(fixture$spectrum, fixture$sig_subtract,
                                    max_neg_fraction = 0.02)
  set.seed(999)
  r_strict <- max_subtract_signature(fixture$spectrum, fixture$sig_subtract,
                                     max_neg_fraction = 0.01)

  expect_gt(r_loose$n_subtract, r_strict$n_subtract)
})


# --- SBS96 tests using COSMIC v3.3 signatures ---

test_that("max_subtract_signature realistic SBS96 at max_neg_fraction=0.02", {
  fixture <- readRDS(test_path("fixtures", "test_96_realistic.rds"))

  set.seed(999)
  result <- max_subtract_signature(fixture$spectrum, fixture$sig_subtract,
                                   max_neg_fraction = 0.02)

  # n_subtract should be close to true (ratio ~1.05 in our trials)
  ratio <- result$n_subtract / fixture$n_target
  expect_gt(ratio, 0.85)
  expect_lt(ratio, 1.15)

  # Residual should resemble true background
  cos_sim <- lsa::cosine(result$residual_sig, fixture$true_bg_sig)[1, 1]
  expect_gt(cos_sim, 0.9)

  # total_negative within limit
  expect_lte(result$total_negative, 0.02 * sum(fixture$spectrum) + 1)

  # Valid probability distribution
  expect_equal(sum(result$residual_sig), 1, tolerance = 1e-10)
  expect_true(all(result$residual_sig >= 0))

  # Monte Carlo p-value should not be extreme
  expect_gt(result$prob_ge_total_negative, 0.01)
})


test_that("max_subtract_signature realistic SBS96 at max_neg_fraction=0.01", {
  fixture <- readRDS(test_path("fixtures", "test_96_realistic.rds"))

  set.seed(999)
  result <- max_subtract_signature(fixture$spectrum, fixture$sig_subtract,
                                   max_neg_fraction = 0.01)

  # At 0.01, more conservative — ratio ~0.99 in our trials
  ratio <- result$n_subtract / fixture$n_target
  expect_gt(ratio, 0.80)
  expect_lt(ratio, 1.10)

  # Residual should resemble true background
  cos_sim <- lsa::cosine(result$residual_sig, fixture$true_bg_sig)[1, 1]
  expect_gt(cos_sim, 0.85)

  # total_negative within limit
  expect_lte(result$total_negative, 0.01 * sum(fixture$spectrum) + 1)

  # Valid probability distribution
  expect_equal(sum(result$residual_sig), 1, tolerance = 1e-10)
  expect_true(all(result$residual_sig >= 0))

  # Monte Carlo p-value should not be extreme
  expect_gt(result$prob_ge_total_negative, 0.01)
})


test_that("stricter max_neg_fraction subtracts fewer mutations for SBS96", {
  fixture <- readRDS(test_path("fixtures", "test_96_realistic.rds"))

  set.seed(999)
  r_loose <- max_subtract_signature(fixture$spectrum, fixture$sig_subtract,
                                    max_neg_fraction = 0.02)
  set.seed(999)
  r_strict <- max_subtract_signature(fixture$spectrum, fixture$sig_subtract,
                                     max_neg_fraction = 0.01)

  expect_gt(r_loose$n_subtract, r_strict$n_subtract)
})
