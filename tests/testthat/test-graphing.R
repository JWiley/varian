test_that("vmp_plot returns a ggplot2 graph", {
  set.seed(1234) # make reproducible
  x <- matrix(rnorm(1000), ncol = 2)

  pm <- vmp_plot(x, plot=FALSE)
  expect_s3_class(pm$Individual[[1]], "ggplot")
})
