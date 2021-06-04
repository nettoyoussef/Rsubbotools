
##############################################################################

# Sepfit
paste0("Sepfit")
library(testthat)
source(paste0(getwd(), "/tests/testthat/test-functions.R"))

test_that("SubLaplace:", {

    # sepfit -V 1 < sublaplace.txt
    #
    #--- FINAL RESULT --------------------------------------------------
    #                           | correlation matrix
    #     value     std.err     |  mu      si      la      al
    # mu=  0.01002   0.001      |    --   0.0000  0.0000 -0.0000
    # si=  8.001     0.001      | -0.0000    --   0.0000 -0.0000
    # la= -0.0009501 0.001      | -0.0000 -0.0000    --  -0.0000
    # al=  0.5       0.001      | -0.0000 -0.0000 -0.0000    --
    #
    #                           Upper triangle: covariances
    #                           Lower triangle: correlation coefficients
    #-------------------------------------------------------------------
    #
    # mu           sigma       lambda        alpha        log-like
    # 1.0025e-02   8.0006e+00  -9.5011e-04   4.9998e-01   4.0788e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(1.0025e-02, 8.0006e+00, -9.5011e-04, 4.9998e-01)
           ,log_likelihood = 4.0788e+00
          , std_error      = c(0.001, 0.001, 0.001, 0.001)
            # we pass the transposed matrix and the code corrects it
          , matrix         =
                c(
                    NA, 0.0000, 0.0000, -0.0000
                   ,-0.0000, NA, 0.0000, -0.0000
                   ,-0.0000, -0.0000, NA, -0.0000
                   ,-0.0000, -0.0000, -0.0000, NA
                )
          , distribution   = "sepfit"
        )
    check_fits(orig_value, .5, sepfit)
    data_test   <- generate_datasets(.5)
    subbo_test <- sepfit(data_test$x, verb = 3)
}
)


test_that("Laplace:", {

    #sepfit -V 1 < laplace.txt
    #
    #--- FINAL RESULT --------------------------------------------------
    #                           | correlation matrix
    #     value     std.err     |  mu      si      la      al
    # mu=  0.00243   0.001      |    --   0.0000  0.0000 -0.0000
    # si=  1.998     0.001      | -0.0000    --   0.0000 -0.0000
    # la= -0.001605  0.001      | -0.0000 -0.0000    --  -0.0000
    # al=  1.001     0.001      | -0.0000 -0.0000 -0.0000    --
    #
    #                           Upper triangle: covariances
    #                           Lower triangle: correlation coefficients
    #-------------------------------------------------------------------
    #
    # mu           sigma         lambda        alpha         log-like
    # 2.4301e-03   1.9983e+00  -1.6047e-03   1.0011e+00   2.3854e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(2.4301e-03, 1.9983e+00, -1.6047e-03, 1.0011e+00)
           ,log_likelihood = 2.3854e+00
          , std_error      = c(0.001, 0.001, 0.001, 0.001)
          , matrix         =
                c(
                    NA, 0.0000, 0.0000, -0.0000
                   ,-0.0000, NA, 0.0000, -0.0000
                   ,-0.0000, -0.0000, NA, -0.0000
                   ,-0.0000, -0.0000, -0.0000, NA
                )
          , distribution   = "sepfit"
        )
    check_fits(orig_value, 1, sepfit)
}
)

test_that("Subnormal:", {

    # sepfit -V 1 < subnormal.txt
    #
    #--- FINAL RESULT --------------------------------------------------
    #                           | correlation matrix
    #     value     std.err     |  mu      si      la      al
    # mu= -4.194e-05 0.001      |    --   0.0000  0.0000 -0.0000
    # si=  1.527     0.001      | -0.0000    --   0.0000 -0.0000
    # la= -0.0006582 0.001      | -0.0000 -0.0000    --  -0.0000
    # al=  1.504     0.001      | -0.0000 -0.0000 -0.0000    --
    #
    #                           Upper triangle: covariances
    #                           Lower triangle: correlation coefficients
    #-------------------------------------------------------------------
    #
    # mu           sigma         lambda        alpha         log-like
    # -4.1938e-05   1.5274e+00  -6.5823e-04   1.5040e+00   1.9504e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(-4.1938e-05, 1.5274e+00, -6.5823e-04, 1.5040e+00)
           ,log_likelihood = 1.9504e+00
          , std_error      = c(0.001, 0.001, 0.001, 0.001)
          , matrix         =
                c(
                    NA, 0.0000, 0.0000, -0.0000
                   ,-0.0000, NA, 0.0000, -0.0000
                   ,-0.0000, -0.0000, NA, -0.0000
                   ,-0.0000, -0.0000, -0.0000, NA
                )
          , distribution   = "sepfit"
        )
    check_fits(orig_value, 1.5, sepfit)
}
)


test_that("Normal:", {

    # sepfit -V 1 < normal.txt
    #
    #--- FINAL RESULT --------------------------------------------------
    #                           | correlation matrix
    #     value     std.err     |  mu      si      la      al
    # mu=  0.0002461 0.001      |    --   0.0000  0.0000 -0.0000
    # si=  1.414     0.001      | -0.0000    --   0.0000 -0.0000
    # la=  0.0001059 0.001      | -0.0000 -0.0000    --  -0.0000
    # al=  2         0.001      | -0.0000 -0.0000 -0.0000    --
    #
    #                           Upper triangle: covariances
    #                           Lower triangle: correlation coefficients
    #-------------------------------------------------------------------
    #
    # mu           sigma         lambda        alpha         log-like
    # 2.4610e-04   1.4140e+00   1.0590e-04   2.0003e+00   1.7657e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(2.4610e-04, 1.4140e+00, 1.0590e-04, 2.0003e+00)
           ,log_likelihood = 1.7657e+00
          , std_error      = c(0.001, 0.001, 0.001, 0.001)
          , matrix         =
                c(
                    NA, 0.0000, 0.0000, -0.0000
                   ,-0.0000, NA, 0.0000, -0.0000
                   ,-0.0000, -0.0000, NA, -0.0000
                   ,-0.0000, -0.0000, -0.0000, NA
                )
          , distribution   = "sepfit"
        )
    check_fits(orig_value, 2, sepfit)
}
)

test_that("SuperNormal:", {

    # sepfit -V 1 < supernormal.txt
    #
    #--- FINAL RESULT --------------------------------------------------
    #                           | correlation matrix
    #     value     std.err     |  mu      si      la      al
    # mu= -0.0005149 0.001      |    --   0.0000  0.0000 -0.0000
    # si=  1.386     0.001      | -0.0000    --   0.0000 -0.0000
    # la= -0.0004806 0.001      | -0.0000 -0.0000    --  -0.0000
    # al=  2.502     0.001      | -0.0000 -0.0000 -0.0000    --
    #
    #                           Upper triangle: covariances
    #                           Lower triangle: correlation coefficients
    #-------------------------------------------------------------------
    #
    # mu           sigma         lambda        alpha         log-like
    # -5.1489e-04   1.3864e+00  -4.8063e-04   2.5015e+00   1.6664e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(-5.1489e-04, 1.3864e+00, -4.8063e-04, 2.5015e+00)
           ,log_likelihood = 1.6664e+00
          , std_error      = c(0.001, 0.001, 0.001, 0.001)
          , matrix         =
                c(
                    NA, 0.0000, 0.0000, -0.0000
                   ,-0.0000, NA, 0.0000, -0.0000
                   ,-0.0000, -0.0000, NA, -0.0000
                   ,-0.0000, -0.0000, -0.0000, NA
                )
          , distribution   = "sepfit"
        )
    check_fits(orig_value, 2.5, sepfit)
}
)

##############################################################################
