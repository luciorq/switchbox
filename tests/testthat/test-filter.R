test_that("swap_filter_wilcoxon returns character vector of feature names", {
  data(trainingData, package = "switchbox")
  result <- swap_filter_wilcoxon(trainingGroup, matTraining, n_features = 20)
  expect_type(result, "character")
  expect_true(length(result) <= 20)
  expect_true(all(result %in% rownames(matTraining)))
})

test_that("swap_filter_wilcoxon with up_down = FALSE returns features", {
  data(trainingData, package = "switchbox")
  result <- swap_filter_wilcoxon(
    trainingGroup,
    matTraining,
    n_features = 10,
    up_down = FALSE
  )
  expect_type(result, "character")
  expect_true(length(result) <= 10)
})

test_that("swap_filter_wilcoxon rejects invalid input", {
  x <- matrix(
    1:6,
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2", "s3"))
  )
  expect_error(
    swap_filter_wilcoxon(c("a", "b", "a"), x),
    "factor"
  )
})
