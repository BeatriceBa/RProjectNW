test_that("returns error if isItAMatch receives a non string sequenceA", {
    expect_error(isItAMatch(1, "sequenceB", 1, 1))
})

test_that("returns error if isItAMatch receives a non string sequenceB", {
    expect_error(isItAMatch("sequenceA", 2, 2, 2))
})

test_that("returns error if isItAMatch receives a non numeric i", {
    expect_error(isItAMatch("sequenceA", "sequenceB", "-", 2))
})

test_that("returns error if isItAMatch receives a non numeric j", {
    expect_error(isItAMatch("sequenceA", "sequenceB", 2, "-"))
})

test_that("returns error if isItAMatch receives a NULL string as a sequenceA", {
    sequenceA <- NULL
    expect_error(isItAMatch(sequenceA,"GCATT", 3, 3))
})

test_that("returns error if isItAMatch receives an empty string as a sequenceB", {
    sequenceB <- ""
    expect_error(isItAMatch("GTT",sequenceB, 3, 3))
})


