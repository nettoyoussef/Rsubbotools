
##############################################################################

# Subboafit

paste0("Subboafit")
library(testthat)
source(paste0(getwd(), "/tests/testthat/test-functions.R"))

test_that("SubLaplace:", {


    # subboafit -V 1 < sublaplace.txt
    #
    #subboafit -V 1 < sublaplace.txt
    #
    #--- FINAL RESULT --------------------------------------------------
    #                           | correlation matrix
    #     value     std.err     |  bl      br      al      ar      m
    # bl=  0.4993    0.0009056  |    --   0.0000  0.0000  0.0000 -0.0000
    # br=  0.4993    0.0009059  |  0.6440    --   0.0000  0.0000  0.0000
    # al=  7.994     0.01536    |  0.2093  0.5040    --   0.0001  0.0000
    # ar=  7.986     0.01534    |  0.5039  0.2089  0.3945    --  -0.0000
    # m =  0.007982  -nan       | -nan -nan -nan -nan    --
    #
    #                           Upper triangle: covariances
    #                           Lower triangle: correlation coefficients
    #-------------------------------------------------------------------
    #
    # bl           br           al           ar           m            log-like
    #4.9932e-01   4.9934e-01   7.9943e+00   7.9859e+00   7.9821e-03   4.0788e+00

    orig_value <-
        generate_orig_dt(
            coef           =
                c(4.9932e-01, 4.9934e-01, 7.9943e+00, 7.9859e+00, 7.9821e-03)
          , log_likelihood = 4.0788e+00
          , std_error      = c(0.0009056, 0.0009059, 0.01536, 0.01534, NaN)
          # we pass the transposed matrix and the code corrects it
          , matrix         =
                c(
                   NA, 0.0000, 0.0000, 0.0000, -0.0000
                 , 0.6440, NA, 0.0000, 0.0000, 0.0000
                 , 0.2093, 0.5040, NA, 0.0001, 0.0000
                 , 0.5039, 0.2089, 0.3945, NA, -0.0000
                 , NA, NA, NA, NA, NA
                )
          , distribution   = "subboafit"
        )
    check_fits(orig_value, .5, subboafit)


    # Generate tests for the covariance matrices
    data_test   <- generate_datasets(.5)
    subbo_test <- subboafit(data_test$x, verb = 3)

}
)


test_that("Laplace:", {

    # subboafit -V 1 < laplace.txt
    #
    #--- FINAL RESULT --------------------------------------------------
    #                           | correlation matrix
    #     value     std.err     |  bl      br      al      ar      m
    # bl=  1.003     0.002636   |    --  -0.0000  0.0000 -0.0000  0.0000
    # br=  1.001     0.002634   | -0.0021    --  -0.0000  0.0000 -0.0000
    # al=  2.004     0.003607   |  0.6199 -0.0034    --  -0.0000  0.0000
    # ar=  1.995     0.003595   | -0.0025  0.6182 -0.0039    --  -0.0000
    # m =  0.006832  0.004686   |  0.5632 -0.5641  0.5881 -0.5877    --
    #
    #                           Upper triangle: covariances
    #                           Lower triangle: correlation coefficients
    #-------------------------------------------------------------------
    #
    # bl           br           al           ar           m            log-like
    # 1.0027e+00   1.0007e+00   2.0039e+00   1.9954e+00   6.8321e-03   2.3854e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(1.0027e+00, 1.0007e+00, 2.0039e+00, 1.9954e+00, 6.8321e-03)
           ,log_likelihood = 2.3854e+00
          , std_error      = c(0.002636, 0.002634, 0.003607, 0.003595, 0.004686)
          , matrix         =
                c(
                    NA, -0.0000, 0.0000, -0.0000, 0.0000
                   ,-0.0021, NA, -0.0000, 0.0000, -0.0000
                   ,0.6199, -0.0034, NA, -0.0000, 0.0000
                   ,-0.0025, 0.6182, -0.0039, NA, -0.0000
                   ,0.5632, -0.5641, 0.5881, -0.5877, NA
                )
          , distribution   = "subboafit"
        )
    check_fits(orig_value, 1, subboafit)

}
)

test_that("Subnormal:", {

    # subboafit -V 1 < subnormal.txt
    #
    #--- FINAL RESULT --------------------------------------------------
    #                           | correlation matrix
    #     value     std.err     |  bl      br      al      ar      m
    # bl=  1.505     0.005966   |    --  -0.0000  0.0000 -0.0000  0.0000
    # br=  1.503     0.005967   | -0.4643    --  -0.0000  0.0000 -0.0000
    # al=  1.53      0.004998   |  0.8585 -0.6045    --  -0.0000  0.0000
    # ar=  1.525     0.00499    | -0.6040  0.8585 -0.7549    --  -0.0000
    # m =  0.003443  0.007837   |  0.7780 -0.7786  0.9149 -0.9150    --
    #
    #                           Upper triangle: covariances
    #                           Lower triangle: correlation coefficients
    #-------------------------------------------------------------------
    #
    # bl           br           al           ar           m            log-like
    # 1.5053e+00   1.5029e+00   1.5301e+00   1.5249e+00   3.4434e-03   1.9504e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(1.5053e+00, 1.5029e+00, 1.5301e+00, 1.5249e+00, 3.4434e-03)
          , log_likelihood = 1.9504e+00
          , std_error      = c(0.005966, 0.005967, 0.004998, 0.00499, 0.007837)
          , matrix         =
                c(
                    NA, -0.0000, 0.0000, -0.0000, 0.0000
                   ,-0.4643, NA, -0.0000, 0.0000, -0.0000
                   ,0.8585, -0.6045, NA, -0.0000, 0.0000
                   ,-0.6040, 0.8585, -0.7549, NA, -0.0000
                   ,0.7780, -0.7786, 0.9149, -0.9150, NA
                )
          , distribution   = "subboafit"
        )
    check_fits(orig_value, 1.5, subboafit)

}
)


test_that("Normal:", {

    # subboafit -V 1 < normal.txt
    #
    #--- FINAL RESULT --------------------------------------------------
    #                           | correlation matrix
    #     value     std.err     |  bl      br      al      ar      m
    # bl=  2.012     0.01143    |    --  -0.0001  0.0001 -0.0001  0.0001
    # br=  1.997     0.01138    | -0.6925    --  -0.0001  0.0001 -0.0001
    # al=  1.421     0.007773   |  0.9248 -0.8091    --  -0.0001  0.0001
    # ar=  1.41      0.007741   | -0.8079  0.9251 -0.9227    --  -0.0001
    # m =  0.008723  0.01156    |  0.8724 -0.8735  0.9729 -0.9730    --
    #
    #                           Upper triangle: covariances
    #                           Lower triangle: correlation coefficients
    #-------------------------------------------------------------------
    #
    # bl           br           al           ar           m            log-like
    # 2.0121e+00   1.9968e+00   1.4214e+00   1.4099e+00   8.7232e-03   1.7657e+00


    orig_value <-
        generate_orig_dt(
            coef           = c(2.0121e+00, 1.9968e+00, 1.4214e+00, 1.4099e+00, 8.7232e-0)
          , log_likelihood = 1.7657e+00
          , std_error      = c(0.01143, 0.01138, 0.007773, 0.007741, 0.01156)
          , matrix         =
                c(
                    NA, -0.0001, 0.0001, -0.0001, 0.0001
                   ,-0.6925, NA, -0.0001, 0.0001, -0.0001
                   ,0.9248, -0.8091, NA, -0.0001, 0.0001
                   ,-0.8079, 0.9251, -0.9227, NA, -0.0001
                   ,0.8724, -0.8735, 0.9729, -0.9730, NA
                )
          , distribution   = "subboafit"
        )
    check_fits(orig_value, 2, subboafit)

}
)

test_that("SuperNormal:", {

    # subboafit -V 1 < supernormal.txt
    #
    #--- FINAL RESULT --------------------------------------------------
    #                           | correlation matrix
    #     value     std.err     |  bl      br      al      ar      m
    # bl=  2.523     0.01928    |    --  -0.0003  0.0002 -0.0002  0.0003
    # br=  2.474     0.01919    | -0.8063    --  -0.0002  0.0002 -0.0003
    # al=  1.403     0.01119    |  0.9520 -0.8903    --  -0.0001  0.0002
    # ar=  1.368     0.0111     | -0.8879  0.9531 -0.9668    --  -0.0002
    # m =  0.02425   0.01573    |  0.9178 -0.9203  0.9879 -0.9882    --
    #
    #                           Upper triangle: covariances
    #                           Lower triangle: correlation coefficients
    #-------------------------------------------------------------------
    #
    # bl           br           al           ar           m            log-like
    # 2.5235e+00   2.4742e+00   1.4030e+00   1.3685e+00   2.4252e-02   1.6664e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(2.5235e+00, 2.4742e+00, 1.4030e+00, 1.3685e+00, 2.4252e-02)
          , log_likelihood = 1.6664e+00
          , std_error      = c(0.01928, 0.01919, 0.01119, 0.0111, 0.01573)
          , matrix         =
                c(
                    NA, -0.0003, 0.0002, -0.0002, 0.0003
                   ,-0.8063, NA, -0.0002, 0.0002, -0.0003
                   ,0.9520, -0.8903, NA, -0.0001, 0.0002
                   ,-0.8879, 0.9531, -0.9668, NA, -0.0002
                   ,0.9178, -0.9203, 0.9879, -0.9882, NA
                )
          , distribution   = "subboafit"
        )
    check_fits(orig_value, 2.5, subboafit)

}
)

##############################################################################
