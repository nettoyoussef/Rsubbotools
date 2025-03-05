##############################################################################

# Sepfit
paste0("Sepfit")
skip_on_cran()

test_that("SubLaplace:", {
  # sepfit -V 1 < sublaplace.txt
  #
  #--- FINAL RESULT --------------------------------------------------
  #                           | correlation matrix
  #     value     std.err     |  mu      si      la      al
  # mu= -0.0002411 0.001      |    --   0.0000  0.0000 -0.0000
  # si=  4.004     0.001      | -0.0000    --   0.0000 -0.0000
  # la=  0.0002655 0.001      | -0.0000 -0.0000    --  -0.0000
  # al=  0.5001    0.001      | -0.0000 -0.0000 -0.0000    --
  #
  #                           Upper triangle: covariances
  #                           Lower triangle: correlation coefficients
  #-------------------------------------------------------------------
  #
  # mu           sigma         lambda        alpha         log-like
  #-2.4108e-04   4.0035e+00   2.6547e-04   5.0010e-01   3.3857e+00

  orig_value <-
    generate_orig_dt(
      coef = c(-2.4108e-04, 4.0035e+00, 2.6547e-04, 5.0010e-01),
      log_likelihood = 3.3857e+00,
      std_error = c(0.001, 0.001, 0.001, 0.001)
      # we pass the transposed matrix and the code corrects it
      , matrix =
        c(
          NA, 0.0000, 0.0000, -0.0000,
          -0.0000, NA, 0.0000, -0.0000,
          -0.0000, -0.0000, NA, -0.0000,
          -0.0000, -0.0000, -0.0000, NA
        ),
      distribution = "sepfit"
    )
  check_fits(orig_value, .5, sepfit)
  data_test <- generate_datasets(.5)
  subbo_test <- sepfit(data_test$x, verb = 3)
})


test_that("Laplace:", {
  # sepfit -V 1 < laplace.txt
  #
  #--- FINAL RESULT --------------------------------------------------
  #                           | correlation matrix
  #     value     std.err     |  mu      si      la      al
  # mu=  0.002005  0.001      |    --   0.0000  0.0000 -0.0000
  # si=  0.9995    0.001      | -0.0000    --   0.0000 -0.0000
  # la= -0.001532  0.001      | -0.0000 -0.0000    --  -0.0000
  # al=  1.001     0.001      | -0.0000 -0.0000 -0.0000    --
  #
  #                           Upper triangle: covariances
  #                           Lower triangle: correlation coefficients
  #-------------------------------------------------------------------
  #
  # mu           sigma         lambda        alpha         log-like
  # 2.0055e-03   9.9954e-01  -1.5315e-03   1.0013e+00   1.6923e+00

  orig_value <-
    generate_orig_dt(
      coef = c(2.0055e-03, 9.9954e-01, -1.5315e-03, 1.0013e+00),
      log_likelihood = 1.6923e+00,
      std_error = c(0.001, 0.001, 0.001, 0.001),
      matrix =
        c(
          NA, 0.0000, 0.0000, -0.0000,
          -0.0000, NA, 0.0000, -0.0000,
          -0.0000, -0.0000, NA, -0.0000,
          -0.0000, -0.0000, -0.0000, NA
        ),
      distribution = "sepfit"
    )
  check_fits(orig_value, 1, sepfit)
})

test_that("Subnormal:", {
  # sepfit -V 1 < subnormal.txt
  #
  #--- FINAL RESULT --------------------------------------------------
  #                           | correlation matrix
  #     value     std.err     |  mu      si      la      al
  # mu= -0.0003492 0.001      |    --   0.0000  0.0000 -0.0000
  # si=  0.7632    0.001      | -0.0000    --   0.0000 -0.0000
  # la= -0.0003853 0.001      | -0.0000 -0.0000    --  -0.0000
  # al=  1.501     0.001      | -0.0000 -0.0000 -0.0000    --
  #
  #                           Upper triangle: covariances
  #                           Lower triangle: correlation coefficients
  #-------------------------------------------------------------------
  #
  # mu           sigma         lambda        alpha         log-like
  #-3.4922e-04   7.6317e-01  -3.8533e-04   1.5010e+00   1.2572e+00

  orig_value <-
    generate_orig_dt(
      coef = c(-3.4922e-04, 7.6317e-01, -3.8533e-04, 1.5010e+00),
      log_likelihood = 1.2572e+00,
      std_error = c(0.001, 0.001, 0.001, 0.001),
      matrix =
        c(
          NA, 0.0000, 0.0000, -0.0000,
          -0.0000, NA, 0.0000, -0.0000,
          -0.0000, -0.0000, NA, -0.0000,
          -0.0000, -0.0000, -0.0000, NA
        ),
      distribution = "sepfit"
    )
  check_fits(orig_value, 1.5, sepfit)
})


test_that("Normal:", {
  # sepfit -V 1 < normal.txt
  #
  #--- FINAL RESULT --------------------------------------------------
  #                           | correlation matrix
  #     value     std.err     |  mu      si      la      al
  # mu=  1.697e-05 0.001      |    --   0.0000  0.0000 -0.0000
  # si=  0.7077    0.001      | -0.0000    --   0.0000 -0.0000
  # la=  7.18e-05  0.001      | -0.0000 -0.0000    --  -0.0000
  # al=  2.003     0.001      | -0.0000 -0.0000 -0.0000    --
  #
  #                           Upper triangle: covariances
  #                           Lower triangle: correlation coefficients
  #-------------------------------------------------------------------
  #
  # mu           sigma         lambda        alpha         log-like
  # 1.6975e-05   7.0771e-01   7.1803e-05   2.0032e+00   1.0725e+00

  orig_value <-
    generate_orig_dt(
      coef = c(1.6975e-05, 7.0771e-01, 7.1803e-05, 2.0032e+00),
      log_likelihood = 1.0725e+00,
      std_error = c(0.001, 0.001, 0.001, 0.001),
      matrix =
        c(
          NA, 0.0000, 0.0000, -0.0000,
          -0.0000, NA, 0.0000, -0.0000,
          -0.0000, -0.0000, NA, -0.0000,
          -0.0000, -0.0000, -0.0000, NA
        ),
      distribution = "sepfit"
    )
  check_fits(orig_value, 2, sepfit)
})

test_that("SuperNormal:", {
  # sepfit -V 1 < supernormal.txt
  #
  #--- FINAL RESULT --------------------------------------------------
  #                           | correlation matrix
  #     value     std.err     |  mu      si      la      al
  # mu= -0.0005465 0.001      |    --   0.0000  0.0000 -0.0000
  # si=  0.6928    0.001      | -0.0000    --   0.0000 -0.0000
  # la=  7.536e-05 0.001      | -0.0000 -0.0000    --  -0.0000
  # al=  2.498     0.001      | -0.0000 -0.0000 -0.0000    --
  #
  #                           Upper triangle: covariances
  #                           Lower triangle: correlation coefficients
  #-------------------------------------------------------------------
  #
  # mu           sigma         lambda        alpha         log-like
  #-5.4652e-04   6.9284e-01   7.5359e-05   2.4982e+00   9.7328e-01

  orig_value <-
    generate_orig_dt(
      coef = c(-5.4652e-04, 6.9284e-01, 7.5359e-05, 2.4982e+00),
      log_likelihood = 9.7328e-01,
      std_error = c(0.001, 0.001, 0.001, 0.001),
      matrix =
        c(
          NA, 0.0000, 0.0000, -0.0000,
          -0.0000, NA, 0.0000, -0.0000,
          -0.0000, -0.0000, NA, -0.0000,
          -0.0000, -0.0000, -0.0000, NA
        ),
      distribution = "sepfit"
    )
  check_fits(orig_value, 2.5, sepfit)
})

##############################################################################
