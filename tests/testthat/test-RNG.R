
##############################################################################

# RNGs

paste0("Random Number Generators")

library(testthat)
source(paste0(getwd(), "/tests/testthat/test-functions.R"))

paste0("Check if setting seed on RNGs work")

test_that("rpower: ", {

    # values for parameters
    # tests each conditional statement
    # SubLaplace, Laplace, mid Norm-Lap,
    # Normal and superNormal, on this order
    b <- c(0.7, 1, 1.5, 2, 2.5)

    lapply(1:length(b), function(i){
        x = compare_RNG(rpower, list(100, 2, b[i]))
        expect_identical(x[["x"]], x[["y"]], label = paste0("parameter b=", b[i]))
    })
})


test_that("rsubbo: ", {

    x = compare_RNG(rsubbo, list())
    expect_identical(x[["x"]], x[["y"]])
})

test_that("rasubbo: ", {

    x = compare_RNG(rasubbo, list())
    expect_identical(x[["x"]], x[["y"]])
})


test_that("rlaplace: ", {

    x = compare_RNG(rlaplace, list(n = 100, a =1))
    expect_identical(x[["x"]], x[["y"]])
})

test_that("rgamma_mt: ", {

    # test argument for all if conditions
    a <- c(0.5, 1, 2)

    lapply(1:length(b), function(i){
        x = compare_RNG(rgamma_mt, list(100, a[i], 1))
        expect_identical(x[["x"]], x[["y"]], label = paste0("parameter b=", b[i]))
    })


})


paste0("Test arguments for functions")

test_that("rgamma_mt: ", {

    # the 'a' parameter must be positive
    expect_error(rgamma_mt(100, -1, 1) )

})

##############################################################################
