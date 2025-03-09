##############################################################################

# Laplafit
# checks against original implementation of subbotools
# Depends on using the same sample to fit as the subbotools package
# since R changes the rng in each version, to work this test must be manually
# updated each time or use a constant sample
skip()
paste0("Laplafit")

test_that("SubLaplace:", {
  # laplafit -V 3 < sublaplace.txt
  #
  #--- FINAL RESULT -----------------------------------
  #                           | correlation matrix
  #      value     std.err    |  m       a
  # m  =  0.001586  0.07193    |  1.0000 -0.0834
  # al =  5.997     8.337e-05  | -0.0834  1.0000
  #----------------------------------------------------
  # m            a            log-like
  # 1.5861e-03   5.9971e+00   3.4844e+00

  orig_value <-
    generate_orig_dt(
      coef = c(1.5861e-03, 5.9971e+00),
      log_likelihood = 3.4844e+00,
      std_error = c(0.07193, 8.337e-05)
      # we pass the transposed matrix and the code corrects it
      , matrix =
        c(
          NA, -0.0834,
          -0.0834, NA
        ),
      distribution = "laplafit"
    )
  check_fits(orig_value, .5, laplafit)
})


test_that("Laplace:", {
  # laplafit -V 3 < laplace.txt
  #
  #--- FINAL RESULT -----------------------------------
  #                           | correlation matrix
  #      value     std.err    |  m       a
  # m  =  0.0007161 0.001996   |  1.0000 -0.5004
  # al =  0.9991    0.0005004  | -0.5004  1.0000
  #----------------------------------------------------
  # m            a            log-like
  # 7.1608e-04   9.9911e-01   1.6923e+00

  orig_value <-
    generate_orig_dt(
      coef = c(7.1608e-04, 9.9911e-01),
      log_likelihood = 1.6923e+00,
      std_error = c(0.001996, 0.0005004),
      matrix =
        c(
          NA, -0.5004,
          -0.5004, NA
        ),
      distribution = "laplafit"
    )
  check_fits(orig_value, 1, laplafit)
})

test_that("Subnormal:", {
  # laplafit -V 3 < subnormal.txt
  #
  #--- FINAL RESULT -----------------------------------
  #                           | correlation matrix
  #      value     std.err    |  m       a
  # m  =  0.0001101 0.0008696  |  1.0000 -0.7583
  # al =  0.6594    0.0007583  | -0.7583  1.0000
  #----------------------------------------------------
  # m            a            log-like
  # 1.1009e-04   6.5941e-01   1.2767e+00

  orig_value <-
    generate_orig_dt(
      coef = c(1.1009e-04, 6.5941e-01),
      log_likelihood = 1.2767e+00,
      std_error = c(0.0008696, 0.0007583),
      matrix =
        c(
          NA, -0.7583,
          -0.7583, NA
        ),
      distribution = "laplafit"
    )
  check_fits(orig_value, 1.5, laplafit)
})


test_that("Normal:", {
  # laplafit -V 3 < normal.txt
  #
  #--- FINAL RESULT -----------------------------------
  #                           | correlation matrix
  #      value     std.err    |  m       a
  # m  =  0.0003519 0.000637   |  1.0000 -0.8860
  # al =  0.5644    0.000886   | -0.8860  1.0000
  #----------------------------------------------------
  # m            a            log-like
  # 3.5190e-04   5.6436e-01   1.1211e+00

  orig_value <-
    generate_orig_dt(
      coef = c(3.5190e-04, 5.6436e-01),
      log_likelihood = 1.1211e+00,
      std_error = c(0.000637, 0.000886),
      matrix =
        c(
          NA, -0.8860,
          -0.8860, NA
        ),
      distribution = "laplafit"
    )
  check_fits(orig_value, 2, laplafit)
})

test_that("SuperNormal:", {
  # laplafit -V 3 < supernormal.txt
  #
  #--- FINAL RESULT -----------------------------------
  #                           | correlation matrix
  #      value     std.err    |  m       a
  # m  =  0.0003685 0.0005506  |  1.0000 -0.9529
  # al =  0.5247    0.0009529  | -0.9529  1.0000
  #----------------------------------------------------
  # m            a            log-like
  # 3.6851e-04   5.2471e-01   1.0482e+00

  orig_value <-
    generate_orig_dt(
      coef = c(3.6851e-04, 5.2471e-01),
      log_likelihood = 1.0482e+00,
      std_error = c(0.0005506, 0.0009529),
      matrix =
        c(
          NA, -0.9529,
          -0.9529, NA
        ),
      distribution = "laplafit"
    )
  check_fits(orig_value, 2.5, laplafit)
})

##############################################################################
