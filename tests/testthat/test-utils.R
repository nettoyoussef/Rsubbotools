library(testthat)

test_that("sortRcpp is equal to sort in R: ", {
  x <- round(rnorm(100) * 100, 0)
  y <- sapply(x, function(i) i) # to make deep copy
  sortRcpp(x)
  expect_equal(x, sort(y))
})



test_that("Newton Method works:", {


    # test utils
    g <- function(x){

        return(x^3-27)    
    }

    newton(0.3,g)
    expect_error(steffensen(0.3,g))
    expect_success(steffensen(1.,g))


})
