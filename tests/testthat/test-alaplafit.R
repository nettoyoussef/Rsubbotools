##############################################################################

# Alaplafit
# checks against original implementation of subbotools
# Depends on using the same sample to fit as the subbotools package
# since R changes the rng in each version, to work this test must be manually
# updated each time or use a constant sample
skip()
paste0("Alaplafit")


test_that("SubLaplace:", {
  # laplaafit -V 3 < sublaplace.txt
  #
  #--- FINAL RESULT -----------------------------------
  #                           | correlation matrix
  #      value     std.err    |  m       al     ar
  # m  =  0.003138  0.07193    |  1.0000  0.0834 -0.0833
  # al =  6.001     0.07198    |  0.0834  1.0000  0.0000
  # ar =  5.993     0.07188    | -0.0833  0.0000  1.0000
  #----------------------------------------------------
  # m            al            ar            log-like
  # 3.1378e-03   6.0011e+00   5.9930e+00   3.4844e+00

  orig_value <-
    generate_orig_dt(
      coef = c(3.1378e-03, 6.0011e+00, 5.9930e+00),
      log_likelihood = 3.4844e+00,
      std_error = c(0.07193, 0.07198, 0.07188)
      # we pass the transposed matrix and the code corrects it
      , matrix =
        c(
          NA, 0.0834, -0.0833,
          0.0834, NA, 0.0000,
          -0.0833, NA, NA
        ),
      distribution = "alaplafit"
    )
  check_fits(orig_value, .5, alaplafit)
})


test_that("Laplace:", {
  # laplaafit -V 3 < laplace.txt
  #
  #--- FINAL RESULT -----------------------------------
  #                           | correlation matrix
  #      value     std.err    |  m       al     ar
  # m  =  0.002405  0.001996   |  1.0000  0.5009 -0.5000
  # al =  1.001     0.002      |  0.5009  1.0000  0.0000
  # ar =  0.9974    0.001993   | -0.5000  0.0000  1.0000
  #----------------------------------------------------
  # m            al            ar            log-like
  # 2.4055e-03   1.0008e+00   9.9740e-01   1.6923e+00

  orig_value <-
    generate_orig_dt(
      coef = c(2.4055e-03, 1.0008e+00, 9.9740e-01),
      log_likelihood = 1.6923e+00,
      std_error = c(0.001996, 0.002, 0.001993),
      matrix =
        c(
          1.0000, 0.5009, -0.5000,
          0.5009, 1.0000, 0.0000,
          -0.5000, 0.0000, 1.0000
        ),
      distribution = "alaplafit"
    )
  check_fits(orig_value, 1, alaplafit)
})

test_that("Subnormal:", {
  # laplaafit -V 3 < subnormal.txt
  #
  #--- FINAL RESULT -----------------------------------
  #                           | correlation matrix
  #      value     std.err    |  m       al     ar
  # m  =  0.001456  0.0008696  |  1.0000  0.7588 -0.7577
  # al =  0.6604    0.0008709  |  0.7588  1.0000  0.0000
  # ar =  0.6584    0.0008683  | -0.7577  0.0000  1.0000
  #----------------------------------------------------
  # m            al            ar            log-like
  # 1.4560e-03   6.6038e-01   6.5843e-01   1.2767e+00


  orig_value <-
    generate_orig_dt(
      coef = c(1.4560e-03, 6.6038e-01, 6.5843e-01),
      log_likelihood = 1.2767e+00,
      std_error = c(0.0008696, 0.0008709, 0.0008683),
      matrix =
        c(
          1.0000, 0.7588, -0.7577,
          0.7588, 1.0000, 0.0000,
          -0.7577, 0.0000, 1.0000
        ),
      distribution = "alaplafit"
    )
  check_fits(orig_value, 1.5, alaplafit)
})


test_that("Normal:", {
  # laplaafit -V 3 < normal.txt
  #
  #--- FINAL RESULT ------------------------------------
  #                            | correlation matrix
  #       value     std.err    |  m       al     ar
  # m  =  0.001671  0.000637   |  1.0000  0.8866 -0.8853
  # al =  0.5652    0.0006379  |  0.8866  1.0000  0.0000
  # ar =  0.5635    0.0006361  | -0.8853  0.0000  1.0000
  #-----------------------------------------------------
  # m            al            ar            log-like
  # 1.6709e-03   5.6517e-01   5.6354e-01   1.1211e+00


  orig_value <-
    generate_orig_dt(
      coef = c(1.6709e-03, 5.6517e-01, 5.6354e-01),
      log_likelihood = 1.1211e+00,
      std_error = c(0.000637, 0.0006379, 0.0006361),
      matrix =
        c(
          1.0000, 0.8866, -0.8853,
          0.8866, 1.0000, 0.0000,
          -0.8853, 0.0000, 1.0000
        ),
      distribution = "alaplafit"
    )
  check_fits(orig_value, 2, alaplafit)
})

test_that("SuperNormal:", {
  # laplaafit -V 3 < supernormal.txt
  #
  #--- FINAL RESULT ------------------------------------
  #                            | correlation matrix
  #       value     std.err    |  m       al     ar
  # m  =  0.002735  0.0005506  |  1.0000  0.9542 -0.9516
  # al =  0.5261    0.0005521  |  0.9542  1.0000  0.0000
  # ar =  0.5233    0.0005491  | -0.9516  0.0000  1.0000
  #-----------------------------------------------------
  # m            al            ar            log-like
  # 2.7352e-03   5.2614e-01   5.2327e-01   1.0482e+00

  orig_value <-
    generate_orig_dt(
      coef = c(2.7352e-03, 5.2614e-01, 5.2327e-01),
      log_likelihood = 1.0482e+00,
      std_error = c(0.0005506, 0.0005521, 0.0005491),
      matrix =
        c(
          1.0000, 0.9542, -0.9516,
          0.9542, 1.0000, 0.0000,
          -0.9516, 0.0000, 1.0000
        ),
      distribution = "alaplafit"
    )
  check_fits(orig_value, 2.5, alaplafit)
})

##############################################################################
