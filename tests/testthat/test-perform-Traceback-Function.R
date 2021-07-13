test_that("returns error if performTraceback receives a non string sequenceA", {
    matrix <- matrix(0, nrow = 3, ncol = 3)
    sequenceA <- 1
    expect_error(performTraceback(matrix,sequenceA,"GCATT"))
})

test_that("returns error if performTraceback receives a non string sequenceA", {
    matrix <- matrix(0, nrow = 3, ncol = 3)
    sequenceB <- 1
    expect_error(performTraceback(matrix,"GTT",sequenceB))
})

test_that("returns error if performTraceback receives a non numeric matrix", {
    matrix <- matrix("-", nrow = 3, ncol = 3)
    expect_error(performTraceback(matrix,"GTT","GCATT"))
})

test_that("returns error if performTraceback receives a NULL matrix", {
    expect_error(performTraceback(NULL,"GTT","GCATT"))
})

test_that("returns error if performTraceback receives a NULL string as a sequenceA", {
    matrix <- matrix(0, nrow = 3, ncol = 3)
    sequenceA <- NULL
    expect_error(performTraceback(matrix,sequenceA,"GCATT"))
})

test_that("returns error if performTraceback receives an empty string as a sequenceB", {
    matrix <- matrix(0, nrow = 3, ncol = 3)
    sequenceB <- ""
    expect_error(performTraceback(matrix,"GTT",sequenceB))
})

test_that("returns error if performTraceback receives a matrix that is not corresponding to the length of sequenceA", {
    matrix <- matrix(0, nrow = 3, ncol = 6)
    sequenceA <- "GTT"
    sequenceB <- "GCATT"
    expect_error(performTraceback(matrix,sequenceA,sequenceB))
})

test_that("returns error if performTraceback receives a matrix that is not corresponding to the length of sequenceB", {
    matrix <- matrix(0, nrow = 4, ncol = 5)
    sequenceA <- "GTT"
    sequenceB <- "GCATT"
    expect_error(performTraceback(matrix,sequenceA,sequenceB))
})
