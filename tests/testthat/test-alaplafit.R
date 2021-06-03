
##############################################################################

# Alaplafit

paste0("Alaplafit")
library(testthat)
source(paste0(getwd(), "/tests/testthat/test-functions.R"))


test_that("SubLaplace:", {

    # laplaafit -V 3 < sublaplace.txt
    #--- START LOADING DATA
    #<<< [stdin]
    #
    #--- FINAL RESULT -----------------------------------
    #                           | correlation matrix
    #      value     std.err    |  m       al     ar
    #m  =  0.006276  0.2877     |  1.0000  0.0417 -0.0417
    #al =  12        0.2879     |  0.0417  1.0000  0.0000
    #ar =  11.99     0.2875     | -0.0417  0.0000  1.0000
    #----------------------------------------------------
    # m            al            ar            log-like
    # 6.2755e-03   1.2002e+01   1.1986e+01   4.1776e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(6.2755e-03, 1.2002e+01, 1.1986e+01)
          , log_likelihood = 4.1776e+00
          , std_error      = c(0.2877, 0.2879, 0.2875)
          # we pass the transposed matrix and the code corrects it
          , matrix         =
                c(
                    NA, 0.0417, -0.0417
                  , 0.0417, NA, 0.0000
                  , -0.0417, NA, 1.0000
                )
          , distribution   = "alaplafit"
        )
    check_fits(orig_value, .5, alaplafit)
}
)


test_that("Laplace:", {

    # laplaafit -V 3 < laplace.txt
    #--- START LOADING DATA
    #<<< [stdin]
    #
    #--- FINAL RESULT -----------------------------------
    #                           | correlation matrix
    #      value     std.err    |  m       al     ar
    #m  =  0.004811  0.007986   |  1.0000  0.2504 -0.2500
    #al =  2.002     0.007999   |  0.2504  1.0000  0.0000
    #ar =  1.995     0.007972   | -0.2500  0.0000  1.0000
    #----------------------------------------------------
    # m            al            ar            log-like
    # 4.8109e-03   2.0016e+00   1.9948e+00   2.3854e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(4.8109e-03, 2.0016e+00, 1.9948e+00)
           ,log_likelihood = 2.3854e+00
          , std_error      = c(0.007986, 0.007999, 0.007972)
          , matrix         =
                c(
                    1.0000, 0.2504, -0.2500
                  , 0.2504, 1.0000, 0.0000
                  , -0.2500, 0.0000, 1.0000
                )
          , distribution   = "alaplafit"
        )
    check_fits(orig_value, 1, alaplafit)
}
)

test_that("Subnormal:", {

    # laplaafit -V 3 < subnormal.txt
    # --- START LOADING DATA
    # <<< [stdin]
    #
    # --- FINAL RESULT -----------------------------------
    #                            | correlation matrix
    #       value     std.err    |  m       al     ar
    # m  =  0.002912  0.003479   |  1.0000  0.3794 -0.3788
    # al =  1.321     0.003484   |  0.3794  1.0000  0.0000
    # ar =  1.317     0.003473   | -0.3788  0.0000  1.0000
    # ----------------------------------------------------
    # m            al            ar            log-like
    # 2.9120e-03   1.3208e+00   1.3169e+00   1.9699e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(2.9120e-03, 1.3208e+00, 1.3169e+00)
          , log_likelihood = 1.9699e+00
          , std_error      = c(0.003479, 0.003484, 0.003473)
          , matrix         =
                c(
                    1.0000, 0.3794, -0.3788
                  , 0.3794, 1.0000, 0.0000
                  , -0.3788, 0.0000, 1.0000
                )
          , distribution   = "alaplafit"
        )
    check_fits(orig_value, 1.5, alaplafit)

}
)


test_that("Normal:", {

    # laplaafit -V 3 < normal.txt
    # --- START LOADING DATA
    # <<< [stdin]
    #
    # --- FINAL RESULT -----------------------------------
    #                            | correlation matrix
    #       value     std.err    |  m       al     ar
    # m  =  0.003342  0.002548   |  1.0000  0.4433 -0.4427
    # al =  1.13      0.002552   |  0.4433  1.0000  0.0000
    # ar =  1.127     0.002544   | -0.4427  0.0000  1.0000
    # ----------------------------------------------------
    # m            al            ar            log-like
    # 3.3419e-03   1.1303e+00   1.1271e+00   1.8142e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(3.3419e-03, 1.1303e+00, 1.1271e+00)
          , log_likelihood = 1.8142e+00
          , std_error      = c(0.002548, 0.002552, 0.002544)
          , matrix         =
                c(
                    1.0000, 0.4433, -0.4427
                   ,0.4433, 1.0000, 0.0000
                   ,-0.4427, 0.0000, 1.0000
                )
          , distribution   = "alaplafit"
        )
    check_fits(orig_value, 2, alaplafit)
}
)

test_that("SuperNormal:", {


    # laplaafit -V 3 < supernormal.txt
    #--- START LOADING DATA
    #<<< [stdin]
    #
    #--- FINAL RESULT -----------------------------------
    #                           | correlation matrix
    #      value     std.err    |  m       al     ar
    #m  =  0.00547   0.002203   |  1.0000  0.4771 -0.4758
    #al =  1.052     0.002209   |  0.4771  1.0000  0.0000
    #ar =  1.047     0.002197   | -0.4758  0.0000  1.0000
    #----------------------------------------------------
    # m            al            ar            log-like
    # 5.4703e-03   1.0523e+00   1.0465e+00   1.7414e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(5.4703e-03, 1.0523e+00, 1.0465e+00)
          , log_likelihood = 1.7414e+00
          , std_error      = c(0.002203, 0.002209, 0.002197)
          , matrix         =
                c(
                    1.0000, 0.4771, -0.4758
                  , 0.4771, 1.0000, 0.0000
                  , -0.4758, 0.0000, 1.0000
                )
          , distribution   = "alaplafit"
        )
    check_fits(orig_value, 2.5, alaplafit)
}
)

##############################################################################
