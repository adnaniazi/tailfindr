context("Rectify")

test_that("Rectification(flipping) is working correctly", {
    data <- c(1, 2, 3, 4, 5, 0, -1, -2, -3, -4, -5)
    expected <- c(1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5)
    expect_equal(rectify(data), expected)
})
