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
