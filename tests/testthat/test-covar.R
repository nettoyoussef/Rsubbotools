test_that("Subbofit Matrix of Covariance:", {
  # Generate tests for the covariance matrices
  data_test <- generate_datasets(.5)
  fit <- subbofit(data_test$x, verb = 3)
  covar_mt <- subbo_covar_r(subbo_test, 10^6)

  expect_equal(fit$matrix, covar_mt)
})
