test_that("matTraining dataset loads correctly", {
  data(trainingData, package = "switchbox")
  expect_true(is.matrix(matTraining))
  expect_true(is.numeric(matTraining))
  expect_equal(nrow(matTraining), 70)
  expect_equal(ncol(matTraining), 78)
  expect_true(!is.null(rownames(matTraining)))
})

test_that("trainingGroup dataset loads correctly", {
  data(trainingData, package = "switchbox")
  expect_true(is.factor(trainingGroup))
  expect_equal(length(trainingGroup), 78)
  expect_equal(length(levels(trainingGroup)), 2)
})

test_that("matTesting dataset loads correctly", {
  data(testingData, package = "switchbox")
  expect_true(is.matrix(matTesting))
  expect_true(is.numeric(matTesting))
  expect_equal(ncol(matTesting), 307)
})

test_that("testingGroup dataset loads correctly", {
  data(testingData, package = "switchbox")
  expect_true(is.factor(testingGroup))
  expect_equal(length(testingGroup), 307)
})
