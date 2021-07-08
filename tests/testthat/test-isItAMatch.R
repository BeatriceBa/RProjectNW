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
