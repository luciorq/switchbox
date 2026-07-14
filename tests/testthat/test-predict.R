test_that("swap_prediction_stats returns correct metrics", {
  preds <- factor(c("a", "b", "a", "b", "a"), levels = c("a", "b"))
  truth <- factor(c("a", "b", "a", "a", "b"), levels = c("a", "b"))

  stats <- swap_prediction_stats(preds, truth, classes = c("a", "b"))
  expect_named(
    stats,
    c("accuracy", "sensitivity", "specificity", "balanced_accuracy")
  )
  expect_equal(stats["accuracy"], c(accuracy = 3 / 5))
})

test_that("swap_prediction_stats with perfect predictions", {
  preds <- factor(c("a", "b", "a", "b"))
  truth <- factor(c("a", "b", "a", "b"))

  stats <- swap_prediction_stats(preds, truth)
  expect_equal(stats["accuracy"], c(accuracy = 1.0))
  expect_equal(stats["sensitivity"], c(sensitivity = 1.0))
  expect_equal(stats["specificity"], c(specificity = 1.0))
})

test_that("swap_prediction_stats includes AUC with decision values", {
  preds <- factor(c("a", "b", "a", "b"))
  truth <- factor(c("a", "b", "a", "b"))
  dv <- c(1, -1, 1, -1)

  stats <- swap_prediction_stats(preds, truth, decision_values = dv)
  expect_true("auc" %in% names(stats))
})

test_that("swap_prediction_stats rejects mismatched lengths", {
  expect_error(
    swap_prediction_stats(factor(c("a", "b")), factor("a")),
    "same length"
  )
})

test_that("swap_ktsp_result returns expected structure", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:5)
  result <- swap_ktsp_result(classifier, matTraining, trainingGroup)

  expect_type(result, "list")
  expect_true("stats" %in% names(result))
  expect_true("roc" %in% names(result))
  expect_true("accuracy" %in% names(result$stats))
})

test_that("swap_ktsp_result includes predictions when requested", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:5)
  result <- swap_ktsp_result(
    classifier,
    matTraining,
    trainingGroup,
    predictions = TRUE,
    decision_values = TRUE
  )
  expect_true("predictions" %in% names(result))
  expect_true("decision_values" %in% names(result))
})

test_that("swap_ktsp_result's roc component matches its own reported AUC", {
  # The AUC implied by result$roc (positive class = higher decision value)
  # must agree with result$stats["auc"], which is computed independently by
  # swap_prediction_stats() from the same decision values. If compute_roc()
  # is called with the positive/negative classes swapped, the two diverge
  # (auc_from_roc == 1 - stats_auc), which also corrupts swap_plot_roc().
  data(trainingData, package = "switchbox")
  data(testingData, package = "switchbox")
  classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:8)
  result <- swap_ktsp_result(classifier, matTesting, testingGroup)

  auc_from_roc <- compute_auc(roc_data = result$roc)
  expect_equal(unname(auc_from_roc), unname(result$stats[["auc"]]))
})

test_that("swap_train_test_results evaluates train and test", {
  data(trainingData, package = "switchbox")
  data(testingData, package = "switchbox")
  result <- swap_train_test_results(
    matTraining,
    trainingGroup,
    matTesting,
    testingGroup,
    k_range = 3:5
  )

  expect_type(result, "list")
  expect_true("classifier" %in% names(result))
  expect_true("train" %in% names(result))
  expect_true("test" %in% names(result))
  expect_true("accuracy" %in% names(result$train))
  expect_true("accuracy" %in% names(result$test))
})
