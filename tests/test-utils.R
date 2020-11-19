


library(testthat)
library("pryr")
#pryr::address(x)

test_that("sortRcpp is equal to sort in R: ", {

    x = round(rnorm(100)*100,0)
    y = sapply(x, function(i) i) # to make deep copy
    sortRcpp(x)
    expect_equal(x, sort(y))
})

