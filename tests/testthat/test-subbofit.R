##############################################################################

# Subbofit
# checks against original implementation of subbotools
# Depends on using the same sample to fit as the subbotools package
# since R changes the rng in each version, to work this test must be manually
# updated each time or use a constant sample

skip()
paste0("Subbofit")

test_that("Subbofit Method of moments:", {
  set.seed(1)
  x <- rsubbo(1e6, 0, 1, 1.5)
  para_mm <- sd(x) / (sum(abs(x - mean(x))) / length(x))

  mm_objf <- function(x, l) {
    # lgamma(x) == log(abs(gamma(x)))
    dtmp1 <- exp(x)
    y <- .5 * lgamma(3. * dtmp1) + .5 * lgamma(dtmp1) -
      lgamma(2. * dtmp1) - (l)

    return(y)
  }

  mm_objdf <- function(x) {
    dtmp1 <- exp(x)
    dy <- dtmp1 * (1.5 * digamma(3. * dtmp1) + .5 * digamma(dtmp1) -
      2. * digamma(2. * dtmp1))

    return(dy)
  }


  mm_fl <- function(l) {
    function(x) {
      mm_objf(x, l)
    }
  }

  mm_f <- mm_fl(log(para_mm))

  test_steffensen <- steffensen(0, mm_f, mm_objdf, max_iter = 3)

  expect_equal(mm_subbotin(para_mm), exp(-test_steffensen$x[3]))
})

##############################################################################

paste0("Test fit against subfamilies of the AEP\n")

paste0(
  " This test check if the following fits work as expected:
- sepfit
- subboafish
- subboafit
- subbofish
- subbofit
- subbolafit
- laplafit
- alaplafit

The tests are conducted in the following way.
First, we generate 5 datasets contemplating each of the special cases of the
rpower generating function, which are related to the b parameter:
- b < 1 or b >4
- b == 1 (Laplacian special case)
- 1 < b < 2
- b == 2 (Gaussian special case)
- 2 < b < 4

We then run each of the original subbotools functions on these datasets with
default parameters (we hard code the results on the test) and compare the
answers against the Rsubbotools implementation.

"
)

test_that("SubLaplace:", {
  # subbofit -V 1 < sublaplace.txt
  #
  #
  #--- FINAL RESULT ----------------------------------
  #
  #     value     std.err     |  b       a       m
  # b =  0.4993    0.0008212  |    --   0.0000  0.0000
  # a =  3.995     0.006409   |  0.4709    --   0.0000
  # m =  0.0006382 -nan       | -nan -nan    --
  #
  #           Upper triangle: covariances
  #           Lower triangle: correlation coefficients
  #---------------------------------------------------
  #
  # b            a            m            log-like
  # 4.9932e-01   3.9950e+00   6.3821e-04   3.3857e+00

  orig_value <-
    generate_orig_dt(
      coef = c(4.9932e-01, 3.9950e+00, 6.3821e-04),
      log_likelihood = 3.3857e+00,
      std_error = c(0.0008212, 0.006409, NaN)
      # we pass the transposed matrix and the code corrects it
      , matrix =
        c(
          NA, 0.0000, 0.0000,
          0.4709, NA, 0.0000,
          NaN, NaN, NA
        ),
      distribution = "subbofit"
    )

  check_fits(orig_value, .5, subbofit)
})

test_that("Laplace:", {
  # subbofit -V 1 < laplace.txt
  #
  #--- FINAL RESULT ----------------------------------
  #
  #     value     std.err     |  b       a       m
  # b =  1.002     0.001861   |    --   0.0000 -0.0000
  # a =  0.9998    0.001271   |  0.6180    --  -0.0000
  # m =  0.0007125 0.001001   | -0.0000 -0.0000    --
  #
  #           Upper triangle: covariances
  #           Lower triangle: correlation coefficients
  #---------------------------------------------------
  #
  # b            a            m            log-like
  # 1.0017e+00   9.9983e-01   7.1252e-04   1.6923e+00

  orig_value <-
    generate_orig_dt(
      coef = c(1.0017e+00, 9.9983e-01, 7.1252e-04),
      log_likelihood = 1.6923e+00,
      std_error = c(0.001861, 0.001271, 0.001001),
      matrix =
        c(
          NA, 0.0000, -0.0000,
          0.6180, NA, -0.0000,
          -0.0000, -0.0000, NA
        ),
      distribution = "subbofit"
    )
  check_fits(orig_value, 1, subbofit)
})

test_that("Subnormal:", {
  # subbofit -V 1 < subnormal.txt
  #
  #--- FINAL RESULT ----------------------------------
  #
  #     value     std.err     |  b       a       m
  # b =  1.504     0.003088   |    --   0.0000 -0.0000
  # a =  0.7637    0.0008741  |  0.7018    --  -0.0000
  # m = -0.0002804 0.0008213  | -0.0000 -0.0000    --
  #
  #           Upper triangle: covariances
  #           Lower triangle: correlation coefficients
  #---------------------------------------------------
  #
  # b            a            m            log-like
  # 1.5041e+00   7.6375e-01  -2.8042e-04   1.2572e+00

  orig_value <-
    generate_orig_dt(
      coef = c(1.5041e+00, 7.6375e-01, -2.8042e-04),
      log_likelihood = 1.2572e+00,
      std_error = c(0.003088, 0.0008741, 0.0008213),
      matrix =
        c(
          NA, 0.0000, -0.0000,
          0.7018, NA, -0.0000,
          -0.0000, -0.0000, NA
        ),
      distribution = "subbofit"
    )
  check_fits(orig_value, 1.5, subbofit)
})

test_that("Normal:", {
  # subbofit -V 1 < normal.txt
  #
  #--- FINAL RESULT ----------------------------------
  #
  #     value     std.err     |  b       a       m
  # b =  2.004     0.004473   |    --   0.0000 -0.0000
  # a =  0.7078    0.0007626  |  0.7551    --  -0.0000
  # m =  3.197e-05 0.0007072  | -0.0000 -0.0000    --
  #
  #           Upper triangle: covariances
  #           Lower triangle: correlation coefficients
  #---------------------------------------------------
  #
  # b            a            m            log-like
  #  2.0044e+00   7.0781e-01   3.1969e-05   1.0725e+00

  orig_value <-
    generate_orig_dt(
      coef = c(2.0044e+00, 7.0781e-01, 3.1969e-05),
      log_likelihood = 1.0725e+00,
      std_error = c(0.004473, 0.0007626, 0.0007072),
      matrix =
        c(
          NA, 0.0000, -0.0000,
          0.7551, NA, -0.0000,
          -0.0000, -0.0000, NA
        ),
      distribution = "subbofit"
    )
  check_fits(orig_value, 2, subbofit)
})

test_that("SuperNormal:", {
  # subbofit -V 1 < supernormal.txt
  #
  #--- FINAL RESULT ----------------------------------
  #
  #     value     std.err     |  b       a       m
  # b =  2.499     0.005986   |    --   0.0000 -0.0000
  # a =  0.6929    0.0007172  |  0.7915    --  -0.0000
  # m = -0.0003205 0.0006301  | -0.0000 -0.0000    --
  #
  #           Upper triangle: covariances
  #           Lower triangle: correlation coefficients
  #---------------------------------------------------
  #
  # b            a            m            log-like
  #  2.4988e+00   6.9286e-01  -3.2053e-04   9.7328e-01

  orig_value <-
    generate_orig_dt(
      coef = c(2.4988e+00, 6.9286e-01, -3.2053e-04),
      log_likelihood = 9.7328e-01,
      std_error = c(0.005986, 0.0007172, 0.0006301),
      matrix =
        c(
          NA, 0.0000, -0.0000,
          0.7915, NA, -0.0000,
          -0.0000, -0.0000, NA
        ),
      distribution = "subbofit"
    )
  check_fits(orig_value, 2.5, subbofit)
})

##############################################################################
