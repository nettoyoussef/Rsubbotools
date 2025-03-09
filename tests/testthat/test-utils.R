library(testthat)

test_that("sortRcpp is equal to sort in R: ", {
  x <- round(rnorm(100) * 100, 0)
  y <- sapply(x, function(i) i) # to make deep copy
  sortRcpp(x)
  expect_equal(x, sort(y))
})



test_that("Newton/Steffensen Method works:", {
  # test utils
  g <- function(x) {
    return(x^3 - 27)
  }

  # success with initial value == 0.3
  optim <- (newton(0.3, g))$x
  expect_equal(round(optim[length(optim)], 4), 3)

  expect_error(steffensen(0.3, g))
  # success with initial value == 1
  optim <- (steffensen(1., g))$x
  expect_equal(round(optim[length(optim)], 4), 3)
})
