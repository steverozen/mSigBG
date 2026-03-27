
test_that("max_subtract_signature recovers minority signature subtraction", {
  set.seed(42)

  sig_a <- runif(89); sig_a <- sig_a / sum(sig_a)
  sig_b <- runif(89); sig_b <- sig_b / sum(sig_b)

  n_a <- 200
  n_b <- 800
  spec_a <- as.numeric(rmultinom(1, size = n_a, prob = sig_a))
  spec_b <- as.numeric(rmultinom(1, size = n_b, prob = sig_b))
  spectrum <- spec_a + spec_b

  result <- max_subtract_signature(spectrum, sig_a)

  # n_subtract should be in the right ballpark (within 2x of true)
  expect_gt(result$n_subtract, 50)
  expect_lt(result$n_subtract, 400)

  # Residual signature should be similar to sig_b
  cos_sim <- lsa::cosine(result$residual_sig, sig_b)[1, 1]
  expect_gt(cos_sim, 0.9)

  # residual_sig is a valid probability distribution
  expect_equal(sum(result$residual_sig), 1, tolerance = 1e-10)
  expect_true(all(result$residual_sig >= 0))

  # n_subtract + n_residual = N
  expect_equal(result$n_subtract + result$n_residual, sum(spectrum))
})


test_that("max_subtract_signature respects target_prob parameter", {
  set.seed(99)

  sig_a <- runif(89); sig_a <- sig_a / sum(sig_a)
  sig_b <- runif(89); sig_b <- sig_b / sum(sig_b)
  spectrum <- as.numeric(rmultinom(1, size = 300, prob = sig_a)) +
              as.numeric(rmultinom(1, size = 700, prob = sig_b))

  r_low  <- max_subtract_signature(spectrum, sig_a, target_prob = 0.01)
  r_high <- max_subtract_signature(spectrum, sig_a, target_prob = 0.9)

  # Lower target_prob allows subtracting more mutations
  expect_gt(r_low$n_subtract, r_high$n_subtract)
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
  expect_true(result$n_subtract >= 0)
  expect_true(result$n_subtract <= sum(spectrum))
})


test_that("max_subtract_signature returns zero when spectrum has no contribution", {
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
