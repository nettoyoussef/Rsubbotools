##############################################################################

# CDFs

paste0("Cumulative Distribution Functions")

library(testthat)

paste0("Check if CDF generators work")

paste0("Gaussian CDF")

b <- c(10^-3, 10^-2, 10^-1, 0, 10, 10^2, 10^3)
# Basically compares the value of the second ar with the third for the vector b
test_that("psubbo: ", {
  test_function(b, pnorm(b), psubbo(b), "parameter b=")
})
test_that("psubbo - m!= 0: ", {
  test_function(b, pnorm(b, 10, 1), psubbo(b, 10), "parameter b=")
})
test_that("pasubbo: ", {
  test_function(b, pnorm(b), pasubbo(b), "parameter b=")
})
test_that("pasubbo - m!= 0: ", {
  test_function(b, pnorm(b, 10, 1), pasubbo(b, 10), "parameter b=")
})
test_that("psep: ", {
  test_function(b, pnorm(b), psep(b, 0, 1, 2, 0), "parameter b=")
})
test_that("psep - m!= 0: ", {
  test_function(b, pnorm(b, 10, 1), psep(b, 10, 1, 2, 0), "parameter b=")
})
test_that("ppower: ", {
  test_function(b, pnorm(b), ppower(b, 0, sqrt(2), 2), "parameter b=")
})
test_that("ppower - m!= 0: ", {
  test_function(b, pnorm(b, 10, 1), ppower(b, 10, sqrt(2), 2), "parameter b=")
})

##############################################################################

paste0("Laplace CDF")

test_that("psubbo: ", {
  test_function(b, plaplace(b), psubbo(b, 0, 1, 1), "parameter b=")
})
test_that("psubbo - m!= 0: ", {
  test_function(b, plaplace(b, 10), psubbo(b, 10, 1, 1), "parameter b=")
})
test_that("pasubbo: ", {
  test_function(b, plaplace(b), pasubbo(b, 0, 1, 1, 1, 1), "parameter b=")
})
test_that("pasubbo - m!= 0: ", {
  test_function(b, plaplace(b, 10), pasubbo(b, 10, 1, 1, 1, 1), "parameter b=")
})
test_that("palaplace: ", {
  test_function(b, plaplace(b), palaplace(b, 0, 1, 1), "parameter b=")
})
test_that("palaplace - m!= 0: ", {
  test_function(b, plaplace(b, 10), palaplace(b, 10, 1, 1), "parameter b=")
})

# for psep, the laplace in the current parametrization is
# not exact for numbers above the location parameter
# we approximate the laplace distribution
# 1-(exp(-|x-m|)/2) by
# exp(x-m)/2
# It is close up to 6 decimal places, but not the same
test_that("psep: ", {
  test_function(b, round(plaplace(b), 6), round(psep(b, 0, 1, 1, 0), 6), "parameter b=")
})
test_that("psep - m!= 0: ", {
  test_function(b, plaplace(b, 10), psep(b, 10, 1, 1, 0), "parameter b=")
})
test_that("ppower: ", {
  test_function(b, plaplace(b), ppower(b, 0, 1, 1), "parameter b=")
})
test_that("ppower - m!= 0: ", {
  test_function(b, plaplace(b, 10, 1), ppower(b, 10, 1, 1), "parameter b=")
})

##############################################################################

paste0("Asymmetric Laplace Quantiles")

test_that("qasubbo: ", {
  test_function(b, palaplace(b, 0, 1, 2), pasubbo(b, 0, 1, 2, 1, 1), "parameter b=")
})
test_that("qasubbo - m!= 0: ", {
  test_function(b, palaplace(b, 10, 1, 2), pasubbo(b, 10, 1, 2, 1, 1), "parameter b=")
})
