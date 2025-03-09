##############################################################################

# Quantiles

paste0("Quantile Functions")

library(testthat)

paste0("Check if quantile generators work")

##############################################################################

paste0("Gaussian Quantiles")

b <- seq(0, 1, by = 0.1)

paste0("Gaussian CDF")

# Basically compares the value of the second ar with the third for the vector b
test_that("qsubbo: ", {
  test_function(b, qnorm(b), qsubbo(b), "parameter b=")
})
test_that("qsubbo - m!= 0: ", {
  test_function(b, qnorm(b, 10, 1), qsubbo(b, 10), "parameter b=")
})
test_that("qasubbo: ", {
  test_function(b, qnorm(b), qasubbo(b), "parameter b=")
})
test_that("qasubbo - m!= 0: ", {
  test_function(b, qnorm(b, 10, 1), qasubbo(b, 10), "parameter b=")
})
test_that("qpower: ", {
  test_function(b, qnorm(b), qpower(b, 0, sqrt(2), 2), "parameter b=")
})
test_that("qpower - m!= 0: ", {
  test_function(b, qnorm(b, 10, 1), qpower(b, 10, sqrt(2), 2), "parameter b=")
})
test_that("qsep: ", {
  test_function(b, qnorm(b), qsep(b, 0, 1, 2, 0), "parameter b=")
})
test_that("qsep - m!= 0: ", {
  test_function(b, qnorm(b, 10, 1), qsep(b, 10, 1, 2, 0), "parameter b=")
})

##############################################################################

paste0("Laplace Quantiles")

test_that("qsubbo: ", {
  test_function(b, qlaplace(b), qsubbo(b, 0, 1, 1), "parameter b=")
})
test_that("qsubbo - m!= 0: ", {
  test_function(b, qlaplace(b, 10), qsubbo(b, 10, 1, 1), "parameter b=")
})
test_that("qasubbo: ", {
  test_function(b, qlaplace(b), qasubbo(b, 0, 1, 1, 1, 1), "parameter b=")
})
test_that("qasubbo - m!= 0: ", {
  test_function(b, qlaplace(b, 10), qasubbo(b, 10, 1, 1, 1, 1), "parameter b=")
})
test_that("qpower: ", {
  test_function(b, qlaplace(b), qpower(b, 0, 1, 1), "parameter b=")
})
test_that("qpower - m!= 0: ", {
  test_function(b, qlaplace(b, 10), qpower(b, 10, 1, 1), "parameter b=")
})
test_that("qalaplace: ", {
  test_function(b, qlaplace(b), qalaplace(b, 0, 1, 1), "parameter b=")
})
test_that("qalaplace - m!= 0: ", {
  test_function(b, qlaplace(b, 10), qalaplace(b, 10, 1, 1), "parameter b=")
})
test_that("qsep: ", {
  test_function(b, qlaplace(b), qsep(b, 0, 1, 1, 0), "parameter b=")
})
test_that("qsep - m!= 0: ", {
  test_function(b, qlaplace(b, 10), qsep(b, 10, 1, 1, 0), "parameter b=")
})

##############################################################################

paste0("Asymmetric Laplace Quantiles")

test_that("qasubbo: ", {
  test_function(b, qalaplace(b, 0, 1, 2), qasubbo(b, 0, 1, 2, 1, 1), "parameter b=")
})
test_that("qasubbo - m!= 0: ", {
  test_function(b, qalaplace(b, 10, 1, 2), qasubbo(b, 10, 1, 2, 1, 1), "parameter b=")
})
