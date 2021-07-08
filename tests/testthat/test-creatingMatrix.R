test_that("returns error if creatingMatrix receives a non numeric matchCost", {
    matchCost <- "a malevolent string"
    expect_error(creatingMatrix(matchCost,-3,-4,"GTT","GCATT"))
})

test_that("returns error if creatingMatrix receives a non numeric mismatchCost", {
    mismatchCost <- "a malevolent string"
    expect_error(creatingMatrix(7,mismatchCost,-4,"GTT","GCATT"))
})

test_that("returns error if creatingMatrix receives a non numeric gapPenalty", {
    gapPenalty <- "a malevolent string"
    expect_error(creatingMatrix(7,-3,gapPenalty,"GTT","GCATT"))
})

test_that("returns error if creatingMatrix receives a non string sequenceA", {
    sequenceA <- 1
    expect_error(creatingMatrix(7,-3,-4,sequenceA,"GCATT"))
})

test_that("returns error if creatingMatrix receives a non string sequenceB", {
    sequenceB <- 1
    expect_error(creatingMatrix(7,-3,-4,"GTT",sequenceB))
})
