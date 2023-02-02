
test_that("subcolor_sim", {
  tmp <- mSigBG::example.spectra$HepG2.cisplatin
  retval <- subcolor_sim(tmp[ , 1, drop = FALSE],
                         tmp[ , 2, drop = FALSE])
  
  testthat::expect_equal(retval, background.info[["MCF10A"]])
  
})
