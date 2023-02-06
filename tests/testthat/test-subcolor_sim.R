
test_that("subcolor_sim", {
  tmp <- mSigBG::example.spectra$HepG2.cisplatin
  retval <- subcolor_sim(tmp[ , 1, drop = FALSE],
                         tmp[ , 2, drop = FALSE])
  
  testthat::expect_equal(
    retval, 
    list(
      by.color = c("C>A" = 0.993039146389516,
                   "C>G" = 0.988993496603531,
                   "C>T" = 0.999088643094178,
                   "T>A" = 0.997937136924211,
                   "T>C" = 0.981165079520172,
                   "T>G" = 0.97535717152383), 
      weighted.avg = 0.993937745433499))
  
})

