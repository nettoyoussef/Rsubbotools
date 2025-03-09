skip_on_cran()
test_that("Subbofit Matrix of Covariance:", {
  # Generate tests for the covariance matrices
  data_test <- generate_datasets(.5)
  fit <- subbofit(data_test$x)
  covar_mt <- subbo_covar_r(fit, 10^6)
  colnames(covar_mt) <- c("b", "a", "m")
  expect_equal(fit$matrix, covar_mt)
})
