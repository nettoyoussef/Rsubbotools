
##############################################################################

# Subbolafit

paste0("Subbolafit")
library(testthat)
source(paste0(getwd(), "/tests/testthat/test-functions.R"))

# this routine does not output the std error for the parameters, neither does
# it returns the covariance matrix

test_that("SubLaplace:", {

    # subbolafit -V 1 < sublaplace.txt
    # bl           br           a            m            log-likelihood
    # 4.989885e-01 4.996534e-01 7.990020e+00 -4.306577e-05 4.078818e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(4.989885e-01, 4.996534e-01, 7.990020e+00, -4.306577e-05)
          , log_likelihood = 4.078818e+00
          , distribution   = "subbolafit"
        )
    check_fits(orig_value, .5, subbolafit)
}
)

test_that("Laplace:", {

    # subbolafit -V 1 < laplace.txt
    # bl           br           a            m            log-likelihood
    # 1.000817e+00 1.002628e+00 1.999673e+00 3.943432e-04 2.385405e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(1.000817e+00, 1.002628e+00, 1.999673e+00, 3.943432e-04)
          , log_likelihood = 2.385405e+00
          , distribution   = "subbolafit"
        )
    check_fits(orig_value, 1, subbolafit)
}
)

test_that("Subnormal:", {

    # subbolafit -V 1 < subnormal.txt
    # bl           br           a             m            log-likelihood
    # 1.502682e+00 1.505532e+00 1.527493e+00 -7.968545e-04 1.950368e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(1.502682e+00, 1.505532e+00, 1.527493e+00, -7.968545e-04)
          , log_likelihood = 1.950368e+00
          , distribution   = "subbolafit"
        )

    check_fits(orig_value, 1.5, subbolafit)
}
)


test_that("Normal:", {

    # subbolafit -V 1 < normal.txt
    # bl           br           a            m            log-likelihood
    # 2.004459e+00 2.004388e+00 1.415614e+00 6.393844e-05 1.765696e+00

    orig_value <-
        generate_orig_dt(
            coef = c(2.004459e+00, 2.004388e+00, 1.415614e+00, 6.393844e-05)
           ,log_likelihood = 1.765696e+00
        )

    check_fits(orig_value, 2, subbolafit)
}
)

test_that("SuperNormal:", {

    # subbolafit -V 1 < supernormal.txt
    # bl           br           a            m            log-likelihood
    # 2.495521e+00 2.502094e+00 1.385714e+00 -2.454156e-04 1.666428e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(2.495521e+00, 2.502094e+00, 1.385714e+00, -2.454156e-04)
          , log_likelihood = 1.666428e+00
          , distribution   = "subbolafit"
        )
    check_fits(orig_value, 2.5, subbolafit)
}
)

##############################################################################
