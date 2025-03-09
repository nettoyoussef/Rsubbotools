##############################################################################

# RNGs
library(testthat)
paste0("Random Number Generators")


##############################################################################

paste0("Check if setting seed on RNGs work")

test_that("rpower: ", {
  test_function_rng_seed(rpower, list(n = 100))
})
test_that("rsubbo: ", {
  test_function_rng_seed(rsubbo, list(n = 100))
})
test_that("rasubbo: ", {
  test_function_rng_seed(rasubbo, list(n = 100))
})
test_that("rlaplace: ", {
  test_function_rng_seed(rlaplace, list(n = 100))
})
test_that("ralaplace: ", {
  test_function_rng_seed(ralaplace, list(n = 100))
})
test_that("rgamma_c: ", {
  test_function_rng_seed(rgamma_c, list(n = 100))
})

##############################################################################

paste0("Check if RNG results are consistent with R's rnorm")

b <- c(10^-3, 10^-2, 10^-1, 0, 10, 10^2, 10^3)

# Basically compares the value of the second ar with the third for the vector b
## Functions may differ in the number of calls to seeds, so I cannot vectorize
## I have to run a loop and set a seed for a unitary sample.
test_that("rsubbo: ", {
  test_function_rng(100, rnorm(1), rsubbo(1))
})
test_that("rasubbo: ", {
  test_function_rng(100, rnorm(1), rasubbo(1))
})
test_that("rpower: ", {
  test_function_rng(100, rnorm(1), rpower(1, 0, sqrt(2), 2))
})
test_that("rpower - m!=0:", {
  test_function_rng(100, rnorm(1, 10), rpower(1, 10, sqrt(2), 2))
})

##############################################################################

paste0("Check if RNG results are consistent with rlaplace")

paste0("Laplacian Density")

test_that("rsubbo: ", {
  test_function_rng(100, rlaplace(1), rsubbo(1, 0, 1, 1))
})
test_that("rasubbo: ", {
  test_function_rng(100, rlaplace(1), rasubbo(1, 0, 1, 1, 1, 1))
})
test_that("ralaplace: ", {
  test_function_rng(100, rlaplace(1), ralaplace(1, 0, 1, 1))
})
test_that("rpower: ", {
  test_function_rng(100, rlaplace(1), rpower(1, 0, 1, 1))
})
test_that("rpower - m!=0: ", {
  test_function_rng(100, rlaplace(1, 10), rpower(1, 10, 1, 1))
})

##############################################################################

paste0("Check if RNG results are consistent with rsubbo")

# rsubbo(n, 0, a, b) == rpower(n, 0, a * pow(b, 1./b), b)
# rsubbo(n, 0, a, b) == rasubbo(n, 0, a, a, b, b)
# rpower(n, )
# E = a*G(1/b)^{1/b}

test_that("rpower - m!=0: ", {
  test_function_rng(100, rlaplace(1, 10), rpower(1, 10, 1, 1))
})

##############################################################################

paste0("Check if RNG results are consistent with ralaplace")

# failling
test_that("rasubbo: ", {
  test_function_rng(100, ralaplace(1, 0, 1, 2), rasubbo(1, 0, 1, 2, 1, 1))
})

##############################################################################

paste0("Test arguments for functions")

test_that("rgamma_c: ", {
  # the 'a' parameter must be positive
  expect_error(rgamma_c(100, -1, 1))
})
