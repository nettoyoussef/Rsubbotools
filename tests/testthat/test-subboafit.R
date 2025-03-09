##############################################################################

# Subboafit
# checks against original implementation of subbotools
# Depends on using the same sample to fit as the subbotools package
# since R changes the rng in each version, to work this test must be manually
# updated each time or use a constant sample
skip()
paste0("Subboafit")

test_that("SubLaplace:", {
  # subboafit -V 1 < sublaplace.txt
  #
  #--- FINAL RESULT --------------------------------------------------
  #                           | correlation matrix
  #     value     std.err     |  bl      br      al      ar      m
  # bl=  0.4988    0.0009042  |    --   0.0000  0.0000  0.0000 -0.0000
  # br=  0.4999    0.0009073  |  0.6438    --   0.0000  0.0000  0.0000
  # al=  3.998     0.007684   |  0.2077  0.5035    --   0.0000  0.0000
  # ar=  3.992     0.007665   |  0.5044  0.2104  0.3945    --  -0.0000
  # m = -0.0005344 -nan       | -nan -nan -nan -nan    --
  #
  #                           Upper triangle: covariances
  #                           Lower triangle: correlation coefficients
  #-------------------------------------------------------------------
  #
  # bl           br           al           ar           m            log-like
  #  4.9876e-01   4.9988e-01   3.9980e+00   3.9921e+00  -5.3440e-04   3.3857e+00


  orig_value <-
    generate_orig_dt(
      coef =
        c(4.9876e-01, 4.9988e-01, 3.9980e+00, 3.9921e+00, -5.3440e-04),
      log_likelihood = 3.3857e+00,
      std_error = c(0.0009042, 0.0009073, 0.007684, 0.007665, NaN)
      # we pass the transposed matrix and the code corrects it
      , matrix =
        c(
          NA, 0.0000, 0.0000, 0.0000, -0.0000,
          0.6438, NA, 0.0000, 0.0000, 0.0000,
          0.2077, 0.5035, NA, 0.0001, 0.0000,
          0.5044, 0.2104, 0.3945, NA, -0.0000,
          NaN, NaN, NaN, NaN, NA
        ),
      distribution = "subboafit"
    )
  check_fits(orig_value, .5, subboafit)


  # Generate tests for the covariance matrices
  data_test <- generate_datasets(.5)
  subbo_test <- subboafit(data_test$x, verb = 3)
})


test_that("Laplace:", {
  # subboafit -V 1 < laplace.txt
  #
  #--- FINAL RESULT --------------------------------------------------
  #                           | correlation matrix
  #     value     std.err     |  bl      br      al      ar      m
  # bl=  1.003     0.002636   |    --  -0.0000  0.0000 -0.0000  0.0000
  # br=  1.001     0.002634   | -0.0021    --  -0.0000  0.0000 -0.0000
  # al=  1.002     0.001803   |  0.6199 -0.0034    --  -0.0000  0.0000
  # ar=  0.9977    0.001798   | -0.0025  0.6182 -0.0039    --  -0.0000
  # m =  0.003416  0.002343   |  0.5632 -0.5641  0.5881 -0.5877    --
  #
  #                           Upper triangle: covariances
  #                           Lower triangle: correlation coefficients
  #-------------------------------------------------------------------
  #
  # bl           br           al           ar           m            log-like
  # 1.0027e+00   1.0007e+00   1.0020e+00   9.9769e-01   3.4161e-03   1.6923e+00

  orig_value <-
    generate_orig_dt(
      coef = c(1.0027e+00, 1.0007e+00, 1.0020e+00, 9.9769e-01, 3.4161e-03),
      log_likelihood = 1.6923e+00,
      std_error = c(0.002636, 0.002634, 0.001803, 0.001798, 0.002343),
      matrix =
        c(
          NA, -0.0000, 0.0000, -0.0000, 0.0000,
          -0.0021, NA, -0.0000, 0.0000, -0.0000,
          0.6199, -0.0034, NA, -0.0000, 0.0000,
          -0.0025, 0.6182, -0.0039, NA, -0.0000,
          0.5632, -0.5641, 0.5881, -0.5877, NA
        ),
      distribution = "subboafit"
    )
  check_fits(orig_value, 1, subboafit)
})

test_that("Subnormal:", {
  # subboafit -V 1 < subnormal.txt
  #
  #--- FINAL RESULT --------------------------------------------------
  #                           | correlation matrix
  #     value     std.err     |  bl      br      al      ar      m
  # bl=  1.505     0.005966   |    --  -0.0000  0.0000 -0.0000  0.0000
  # br=  1.503     0.005967   | -0.4643    --  -0.0000  0.0000 -0.0000
  # al=  0.765     0.002499   |  0.8585 -0.6045    --  -0.0000  0.0000
  # ar=  0.7625    0.002495   | -0.6040  0.8585 -0.7549    --  -0.0000
  # m =  0.001722  0.003918   |  0.7780 -0.7786  0.9149 -0.9150    --
  #
  #                           Upper triangle: covariances
  #                           Lower triangle: correlation coefficients
  #-------------------------------------------------------------------
  #
  # bl           br           al           ar           m            log-like
  # 1.5053e+00   1.5029e+00   7.6504e-01   7.6245e-01   1.7217e-03   1.2572e+00

  orig_value <-
    generate_orig_dt(
      coef = c(1.5053e+00, 1.5029e+00, 7.6504e-01, 7.6245e-01, 1.7217e-03),
      log_likelihood = 1.2572e+00,
      std_error = c(0.005966, 0.005967, 0.002499, 0.002495, 0.003918),
      matrix =
        c(
          NA, -0.0000, 0.0000, -0.0000, 0.0000,
          -0.4643, NA, -0.0000, 0.0000, -0.0000,
          0.8585, -0.6045, NA, -0.0000, 0.0000,
          -0.6040, 0.8585, -0.7549, NA, -0.0000,
          0.7780, -0.7786, 0.9149, -0.9150, NA
        ),
      distribution = "subboafit"
    )
  check_fits(orig_value, 1.5, subboafit)
})


test_that("Normal:", {
  # subboafit -V 1 < normal.txt
  #
  #--- FINAL RESULT --------------------------------------------------
  #                           | correlation matrix
  #     value     std.err     |  bl      br      al      ar      m
  # bl=  2.012     0.01143    |    --  -0.0001  0.0000 -0.0000  0.0001
  # br=  1.997     0.01138    | -0.6925    --  -0.0000  0.0000 -0.0001
  # al=  0.7107    0.003887   |  0.9248 -0.8091    --  -0.0000  0.0000
  # ar=  0.7049    0.00387    | -0.8080  0.9251 -0.9227    --  -0.0000
  # m =  0.004351  0.005779   |  0.8724 -0.8735  0.9729 -0.9730    --
  #
  #                           Upper triangle: covariances
  #                           Lower triangle: correlation coefficients
  #-------------------------------------------------------------------
  #
  # bl           br           al           ar           m            log-like
  # 2.0121e+00   1.9968e+00   7.1068e-01   7.0494e-01   4.3510e-03   1.0725e+00

  orig_value <-
    generate_orig_dt(
      coef = c(2.0121e+00, 1.9968e+00, 7.1068e-01, 7.0494e-01, 4.3510e-03),
      log_likelihood = 1.0725e+00,
      std_error = c(0.01143, 0.01138, 0.003887, 0.00387, 0.005779),
      matrix =
        c(
          NA, -0.0001, 0.0000, -0.0000, 0.0001,
          -0.6925, NA, -0.0000, 0.0000, -0.0001,
          0.9248, -0.8091, NA, -0.0000, 0.0000,
          -0.8079, 0.9251, -0.9227, NA, -0.0000,
          0.8724, -0.8735, 0.9729, -0.9730, NA
        ),
      distribution = "subboafit"
    )
  check_fits(orig_value, 2, subboafit)
})

test_that("SuperNormal:", {
  # subboafit -V 1 < supernormal.txt
  #
  #--- FINAL RESULT --------------------------------------------------
  #                           | correlation matrix
  #     value     std.err     |  bl      br      al      ar      m
  # bl=  2.523     0.01928    |    --  -0.0003  0.0001 -0.0001  0.0001
  # br=  2.474     0.01919    | -0.8063    --  -0.0001  0.0001 -0.0001
  # al=  0.7015    0.005593   |  0.9520 -0.8903    --  -0.0000  0.0000
  # ar=  0.6842    0.005548   | -0.8879  0.9531 -0.9668    --  -0.0000
  # m =  0.01212   0.007865   |  0.9178 -0.9203  0.9879 -0.9882    --
  #
  #                           Upper triangle: covariances
  #                           Lower triangle: correlation coefficients
  #-------------------------------------------------------------------
  #
  # bl           br           al           ar           m            log-like
  # 2.5234e+00   2.4742e+00   7.0150e-01   6.8424e-01   1.2117e-02   9.7328e-01


  orig_value <-
    generate_orig_dt(
      coef = c(2.5234e+00, 2.4742e+00, 7.0150e-01, 6.8424e-01, 1.2117e-02),
      log_likelihood = 9.7328e-01,
      std_error = c(0.01928, 0.01919, 0.005593, 0.005548, 0.007865),
      matrix =
        c(
          NA, -0.0003, 0.0001, -0.0001, 0.0001,
          -0.8063, NA, -0.0001, 0.0001, -0.0001,
          0.9520, -0.8903, NA, -0.0000, 0.0000,
          -0.8879, 0.9531, -0.9668, NA, -0.0000,
          0.9178, -0.9203, 0.9879, -0.9882, NA
        ),
      distribution = "subboafit"
    )
  check_fits(orig_value, 2.5, subboafit)
})

##############################################################################
