##############################################################################

# Subbofit

paste0("Subbofit")
library(testthat)
source(paste0(getwd(), "/tests/testthat/test-functions.R"))

test_that("Subbofit Method of moments:", {
    set.seed(1)
    x = rsubbo(1e6, 0, 1, 1.5)
    para_mm = sd(x)/(sum(abs(x- mean(x)))/length(x))

    mm_objf <- function(x, l){

        # lgamma(x) == log(abs(gamma(x)))
        dtmp1 = exp(x);
        y = .5*lgamma(3.*dtmp1)+.5*lgamma(dtmp1)-
            lgamma(2.*dtmp1)-(l);

        return(y);
    }

    mm_objdf <- function(x){

        dtmp1 = exp(x);
        dy = dtmp1*(1.5*digamma(3.*dtmp1)+.5*digamma(dtmp1)-
                    2.*digamma(2.*dtmp1));

        return(dy);
    }


    mm_fl <- function(l){
        function(x){
            mm_obj(x, l)
        }
    }

    mm_f <- mm_fl(log(para_mm))

    test_steffensen <- steffensen(0,mm_f, mm_objdf, max_iter = 3)

    expect_equal(mm(para_mm), exp(-test_steffensen$x[3]))
}
)

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
    # a =  7.99      0.01282    |  0.4709    --   0.0000
    # m = -0.001069  -nan       | -nan -nan    --
    #
    #           Upper triangle: covariances
    #           Lower triangle: correlation coefficients
    #---------------------------------------------------
    #
    # b            a            m            log-like
    # 4.9932e-01   7.9900e+00  -1.0688e-03   4.0788e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(4.9932e-01, 7.9900e+00, -1.0688e-03)
           ,log_likelihood = 4.0788e+00
           ,std_error      = c(0.0008212, 0.01282, NaN)
          # we pass the transposed matrix and the code corrects it
          , matrix         =
                c(
                    NA     , 0.0000, 0.0000
                   , 0.4709, NA    , 0.0000
                   , NaN    , NaN    , NA
                )
          , distribution   = "subbofit"
        )

    check_fits(orig_value, .5, subbofit)
}
)

test_that("Laplace:", {

    # subbofit -V 1 < laplace.txt
    #
    #--- FINAL RESULT ----------------------------------
    #
    #     value     std.err     |  b       a       m
    # b =  1.002     0.001861   |    --   0.0000 -0.0000
    # a =  2         0.002541   |  0.6180    --  -0.0000
    # m =  0.001425  0.002002   | -0.0000 -0.0000    --
    #
    #           Upper triangle: covariances
    #           Lower triangle: correlation coefficients
    #---------------------------------------------------
    #
    # b            a            m            log-like
    #  1.0017e+00   1.9997e+00   1.4250e-03   2.3854e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(1.002, 2, 0.001425)
           ,log_likelihood = 2.3854e+00
           ,std_error      = c(0.001861, 0.002541, 0.002002)
           ,matrix         =
                c(
                    NA, 0.0000, -0.0000
                   ,0.6180, NA, -0.0000
                   ,-0.0000, -0.0000, NA
                )
          , distribution   = "subbofit"
        )
    check_fits(orig_value, 1, subbofit)

}
)

test_that("Subnormal:", {


    # subbofit -V 1 < subnormal.txt
    #
    #--- FINAL RESULT ----------------------------------
    #
    #     value     std.err     |  b       a       m
    # b =  1.504     0.003088   |    --   0.0000 -0.0000
    # a =  1.527     0.001748   |  0.7018    --  -0.0000
    # m = -0.0005608 0.001643   | -0.0000 -0.0000    --
    #
    #           Upper triangle: covariances
    #           Lower triangle: correlation coefficients
    #---------------------------------------------------
    #
    # b            a            m            log-like
    #  1.5041e+00   1.5275e+00  -5.6084e-04   1.9504e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(1.504, 1.527, -0.0005608)
           ,log_likelihood = 1.9504e+00
           ,std_error      = c(0.003088, 0.001748, 0.001643)
          , matrix         =
                c(
                    NA, 0.0000, -0.0000
                   ,0.7018, NA, -0.0000
                   ,-0.0000, -0.0000, NA
                )
          , distribution   = "subbofit"
        )
    check_fits(orig_value, 1.5, subbofit)
}
)

test_that("Normal:", {


    # subbofit -V 1 < normal.txt
    #
    #--- FINAL RESULT ----------------------------------
    #
    #     value     std.err     |  b       a       m
    # b =  2.004     0.004473   |    --   0.0000 -0.0000
    # a =  1.416     0.001525   |  0.7551    --  -0.0000
    # m =  6.394e-05 0.001414   | -0.0000 -0.0000    --
    #
    #           Upper triangle: covariances
    #           Lower triangle: correlation coefficients
    #---------------------------------------------------
    #
    # b            a            m            log-like
    # 2.0044e+00   1.4156e+00   6.3938e-05   1.7657e+00


    orig_value <-
        generate_orig_dt(
            coef           = c(2.0044e+00, 1.4156e+00, 6.3938e-05)
           ,log_likelihood = 1.7657e+00
           ,std_error      = c(0.004473, 0.001525, 0.001414)
           ,matrix         =
                c(
                    NA, 0.0000, -0.0000
                   ,0.7551, NA, -0.0000
                   ,-0.0000, -0.0000, NA
                )
          , distribution   = "subbofit"
        )
    check_fits(orig_value, 2, subbofit)

}
)

test_that("SuperNormal:", {

    # subbofit -V 1 < supernormal.txt
    #
    #--- Final RESULT ----------------------------------
    #
    #     value     std.err     |  b       a       m
    # b =  2.499     0.005986   |    --   0.0000 -0.0000
    # a =  1.386     0.001434   |  0.7915    --  -0.0000
    # m = -0.0006411 0.00126    | -0.0000 -0.0000    --
    #
    #           Upper triangle: covariances
    #           Lower triangle: correlation coefficients
    #---------------------------------------------------
    #
    # b            a            m            log-like
    # 2.4988e+00   1.3857e+00  -6.4107e-04   1.6664e+00

    orig_value <-
        generate_orig_dt(
            coef           = c(2.499,  1.386, -0.0006411)
           ,log_likelihood = 1.6664e+00
           ,std_error      = c(0.005986,  0.001434, 0.00126)
           ,matrix         =
                c(
                    NA, 0.0000, -0.0000
                   ,0.7915, NA, -0.0000
                   ,-0.0000, -0.0000, NA
                )
          , distribution   = "subbofit"
        )
    check_fits(orig_value, 2.5, subbofit)
}
)

##############################################################################
