
test_that("MakeBackgroundInfo", {

  retval <- MakeBackgroundInfo(bg.spectra = MCF10A.background.spectra,
                               title = "MCF10A.background")

  testthat::expect_equal(retval, background.info[["MCF10A"]])

})
