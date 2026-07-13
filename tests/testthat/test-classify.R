test_that("swap_ktsp_statistics returns expected structure", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:5)
  stats <- swap_ktsp_statistics(matTraining, classifier)

  expect_type(stats, "list")
  expect_true("statistics" %in% names(stats))
  expect_true("comparisons" %in% names(stats))
  expect_length(stats$statistics, ncol(matTraining))
  expect_equal(nrow(stats$comparisons), ncol(matTraining))
  expect_equal(ncol(stats$comparisons), nrow(classifier$TSPs))
})

test_that("swap_ktsp_statistics works with single sample", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:5)
  stats <- swap_ktsp_statistics(matTraining[, 1, drop = FALSE], classifier)
  expect_length(stats$statistics, 1)
})

test_that("swap_ktsp_classify returns factor with correct levels", {
  data(trainingData, package = "switchbox")
  data(testingData, package = "switchbox")
  classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:5)
  preds <- swap_ktsp_classify(matTesting, classifier)

  expect_s3_class(preds, "factor")
  expect_length(preds, ncol(matTesting))
  expect_true(all(levels(preds) %in% classifier$labels))
})

test_that("swap_ktsp_classify predictions are deterministic", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:5)
  preds1 <- swap_ktsp_classify(matTraining, classifier)
  preds2 <- swap_ktsp_classify(matTraining, classifier)
  expect_identical(preds1, preds2)
})

test_that("swap_ktsp_statistics works with tie-aware classifier", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_ktsp(
    matTraining,
    trainingGroup,
    handle_ties = TRUE,
    k_range = 3:5
  )
  stats <- swap_ktsp_statistics(matTraining, classifier)
  expect_length(stats$statistics, ncol(matTraining))
})

test_that("swap_ktsp_statistics rejects non-numeric matrix", {
  classifier <- list(TSPs = matrix(c("g1", "g2"), nrow = 1))
  x <- matrix(
    letters[1:4],
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2"))
  )
  expect_error(
    swap_ktsp_statistics(x, classifier),
    "numeric matrix"
  )
})
