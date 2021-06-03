
##############################################################################

# Laplafit

paste0("Laplafit")
library(testthat)
source(paste0(getwd(), "/tests/testthat/test-functions.R"))

test_that("SubLaplace:", {

    # laplafit -V 3 < sublaplace.txt
    #--- START LOADING DATA
    #<<< [stdin]
    #
    #--- FINAL RESULT -----------------------------------
    #                           | correlation matrix
    #      value     std.err    |  m       a
    #m  =  0.003172  0.2877     |  1.0000 -0.0417
    #al =  11.99     4.169e-05  | -0.0417  1.0000
    #----------------------------------------------------
    # m            a            log-like
    # 3.1721e-03   1.1994e+01   4.1776e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(3.1721e-03, 1.1994e+01)
           ,log_likelihood = 4.1776e+00
          , std_error      = c(0.2877, 4.169e-05)
            # we pass the transposed matrix and the code corrects it
          , matrix         =
                c(
                    NA, -0.0417
                  , -0.0417, NA
                )
          , distribution   = "laplafit"
        )
    check_fits(orig_value, .5, laplafit)
}
)


test_that("Laplace:", {

    # laplafit -V 3 < laplace.txt
    #--- START LOADING DATA
    #<<< [stdin]
    #
    #--- FINAL RESULT -----------------------------------
    #                           | correlation matrix
    #      value     std.err    |  m       a
    #m  =  0.001432  0.007986   |  1.0000 -0.2502
    #al =  1.998     0.0002502  | -0.2502  1.0000
    #----------------------------------------------------
    # m            a            log-like
    # 1.4322e-03   1.9982e+00   2.3854e+00



    orig_value <-
        generate_orig_dt(
            coef           = c(1.4322e-03, 1.9982e+00)
           ,log_likelihood = 2.3854e+00
          , std_error      = c(0.007986, 0.0002502)
          , matrix         =
                c(
                    NA, -0.2502
                  , -0.2502, NA
                )
          , distribution   = "laplafit"
        )
    check_fits(orig_value, 1, laplafit)
}
)

test_that("Subnormal:", {

    # laplafit -V 3 < subnormal.txt
    #--- START LOADING DATA
    #<<< [stdin]
    #
    #--- FINAL RESULT -----------------------------------
    #                           | correlation matrix
    #      value     std.err    |  m       a
    #m  =  0.0002202 0.003479   |  1.0000 -0.3791
    #al =  1.319     0.0003791  | -0.3791  1.0000
    #----------------------------------------------------
    # m            a            log-like
    # 2.2018e-04   1.3188e+00   1.9699e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(2.2018e-04, 1.3188e+00)
           ,log_likelihood = 1.9699e+00
          , std_error      = c(0.003479, 0.0003791)
          , matrix         =
                c(
                    NA, -0.3791
                  , -0.3791, NA
                )
          , distribution   = "laplafit"
        )
    check_fits(orig_value, 1.5, laplafit)
}
)


test_that("Normal:", {

    # laplafit -V 3 < normal.txt
    #--- START LOADING DATA
    #<<< [stdin]
    #
    #--- FINAL RESULT -----------------------------------
    #                           | correlation matrix
    #      value     std.err    |  m       a
    #m  =  0.0007038 0.002548   |  1.0000 -0.4430
    #al =  1.129     0.000443   | -0.4430  1.0000
    #----------------------------------------------------
    # m            a            log-like
    # 7.0381e-04   1.1287e+00   1.8142e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(7.0381e-04, 1.1287e+00)
           ,log_likelihood = 1.8142e+00
          , std_error      = c(0.002548, 0.000443)
          , matrix         =
                c(
                    NA, -0.4430
                  , -0.4430, NA
                )
          , distribution   = "laplafit"
        )
    check_fits(orig_value, 2, laplafit)
}
)

test_that("SuperNormal:", {

    # laplafit -V 3 < supernormal.txt
    #--- START LOADING DATA
    #<<< [stdin]
    #
    #--- FINAL RESULT -----------------------------------
    #                           | correlation matrix
    #      value     std.err    |  m       a
    #m  =  0.000737  0.002203   |  1.0000 -0.4765
    #al =  1.049     0.0004765  | -0.4765  1.0000
    #----------------------------------------------------
    # m            a            log-like
    # 7.3702e-04   1.0494e+00   1.7414e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(7.3702e-04, 1.0494e+00)
           ,log_likelihood = 1.7414e+00
          , std_error      = c(0.002203, 0.0004765)
          , matrix         =
                c(
                    NA, -0.4765
                  , -0.4765, NA
                )
          , distribution   = "laplafit"
        )
    check_fits(orig_value, 2.5, laplafit)
}
)

##############################################################################
