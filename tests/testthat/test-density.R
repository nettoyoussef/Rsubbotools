##############################################################################

# Densities

paste0("Density Generators")

library(testthat)

paste0("Check if density generators work")

paste0("Gaussian ")

b <- c(10^-3, 10^-2, 10^-1, 0, 10, 10^2, 10^3)

# Basically compares the value of the second ar with the third for the vector b
test_that("dsubbo: ", {
  test_function(b, dnorm(b), dsubbo(b), "parameter b=")
})
test_that("dsubbo - m!= 0: ", {
  test_function(b, dnorm(b, 10, 1), dsubbo(b, 10), "parameter b=")
})
test_that("dasubbo: ", {
  test_function(b, dnorm(b), dasubbo(b), "parameter b=")
})
test_that("dasubbo - m!= 0: ", {
  test_function(b, dnorm(b, 10, 1), dasubbo(b, 10), "parameter b=")
})
test_that("dsep: ", {
  test_function(b, dnorm(b), dsep(b, 0, 1, 2, 0), "parameter b=")
})
test_that("dsep - m!= 0: ", {
  test_function(b, dnorm(b, 10, 1), dsep(b, 10, 1, 2, 0), "parameter b=")
})
test_that("dpower: ", {
  test_function(b, dnorm(b), dpower(b, 0, sqrt(2), 2), "parameter b=")
})
test_that("dpower - m!= 0: ", {
  test_function(b, dnorm(b, 10, 1), dpower(b, 10, sqrt(2), 2), "parameter b=")
})

##############################################################################

paste0("Laplacian Density")

test_that("dsubbo: ", {
  test_function(b, dlaplace(b), dsubbo(b, 0, 1, 1), "parameter b=")
})
test_that("dsubbo - m!= 0: ", {
  test_function(b, dlaplace(b, 10), dsubbo(b, 10, 1, 1), "parameter b=")
})
test_that("dasubbo: ", {
  test_function(b, dlaplace(b), dasubbo(b, 0, 1, 1, 1, 1), "parameter b=")
})
test_that("dasubbo - m!= 0: ", {
  test_function(b, dlaplace(b, 10), dasubbo(b, 10, 1, 1, 1, 1), "parameter b=")
})
test_that("dalaplace: ", {
  test_function(b, dlaplace(b), dalaplace(b, 0, 1, 1), "parameter b=")
})
test_that("dalaplace - m!= 0: ", {
  test_function(b, dlaplace(b, 10), dalaplace(b, 10, 1, 1), "parameter b=")
})
test_that("dsep: ", {
  test_function(b, dlaplace(b), dsep(b, 0, 1, 1, 0), "parameter b=")
})
test_that("dsep - m!= 0: ", {
  test_function(b, dlaplace(b, 10), dsep(b, 10, 1, 1, 0), "parameter b=")
})
test_that("dpower: ", {
  test_function(b, dlaplace(b), dpower(b, 0, 1, 1), "parameter b=")
})
test_that("dpower - m!= 0: ", {
  test_function(b, dlaplace(b, 10, 1), dpower(b, 10, 1, 1), "parameter b=")
})

##############################################################################

paste0("Asymmetric Laplace Quantiles")

test_that("dasubbo: ", {
  test_function(b, dalaplace(b, 0, 1, 2), dasubbo(b, 0, 1, 2, 1, 1), "parameter b=")
})
test_that("dasubbo - m!= 0: ", {
  test_function(b, dalaplace(b, 10, 1, 2), dasubbo(b, 10, 1, 2, 1, 1), "parameter b=")
})

##############################################################################
