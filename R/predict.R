#' Compute Prediction Statistics
#'
#' Calculates accuracy, sensitivity, specificity, balanced accuracy, and
#' optionally AUC from predictions and true labels.
#'
#' @param predictions A factor or vector of predicted class labels.
#' @param truth A factor or vector of true class labels.
#' @param classes Character vector of length 2: `c(positive, negative)`. If
#'   `NULL`, derived from `levels(truth)`.
#' @param decision_values An optional numeric vector of decision/score values.
#'   If provided, AUC is also computed.
#'
#' @return A named numeric vector with elements `accuracy`, `sensitivity`,
#'   `specificity`, `balanced_accuracy`, and optionally `auc`.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:8)
#' preds <- swap_ktsp_classify(matTraining, classifier)
#' swap_prediction_stats(preds, trainingGroup)
swap_prediction_stats <- function(
  predictions,
  truth,
  classes = NULL,
  decision_values = NULL
) {
  if (length(predictions) != length(truth)) {
    stop(
      "Predictions and truth vectors must have the same length.",
      call. = FALSE
    )
  }

  if (is.null(classes)) {
    classes <- rev(levels(factor(truth)))
  }

  accuracy <- sum(predictions == truth) / length(predictions)
  sensitivity <- sum(
    predictions[truth == classes[1]] == truth[truth == classes[1]]
  ) /
    sum(truth == classes[1])
  specificity <- sum(
    predictions[truth == classes[2]] == truth[truth == classes[2]]
  ) /
    sum(truth == classes[2])
  balanced_accuracy <- (sensitivity + specificity) / 2

  if (is.null(decision_values)) {
    c(
      accuracy = accuracy,
      sensitivity = sensitivity,
      specificity = specificity,
      balanced_accuracy = balanced_accuracy
    )
  } else {
    auc <- tryCatch(
      compute_auc(truth, decision_values, classes),
      error = function(e) NA_real_
    )
    c(
      accuracy = accuracy,
      sensitivity = sensitivity,
      specificity = specificity,
      balanced_accuracy = balanced_accuracy,
      auc = auc
    )
  }
}


#' Get KTSP Classification Result
#'
#' Classifies samples and computes prediction statistics in a single call.
#'
#' @param classifier A KTSP classifier.
#' @param x Numeric matrix (features x samples).
#' @param y Factor of true class labels.
#' @param classes Character vector: `c(case, control)`. If `NULL`, derived
#'   from `levels(y)`.
#' @param predictions Logical; include predictions in the output?
#' @param decision_values Logical; include decision values in the output?
#'
#' @return A list with:
#'   - `stats`: named numeric vector of performance statistics
#'   - `roc`: ROC curve data frame (from internal `compute_roc()`)
#'   - `predictions` (if requested): factor of predictions
#'   - `decision_values` (if requested): numeric vector
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:8)
#' result <- swap_ktsp_result(classifier, matTraining, trainingGroup)
#' result$stats
swap_ktsp_result <- function(
  classifier,
  x,
  y,
  classes = NULL,
  predictions = FALSE,
  decision_values = FALSE
) {
  if (is.null(classes)) {
    classes <- rev(levels(factor(y)))
  }

  fit_predictions <- swap_ktsp_classify(x, classifier)
  fit_dv <- swap_ktsp_statistics(x, classifier)$statistics

  stats_result <- swap_prediction_stats(
    predictions = fit_predictions,
    truth = y,
    classes = classes,
    decision_values = fit_dv
  )

  roc_data <- tryCatch(
    compute_roc(y, fit_dv, classes = rev(classes)),
    error = function(e) NULL
  )

  result <- list(stats = stats_result, roc = roc_data)
  if (predictions) {
    result$predictions <- fit_predictions
  }
  if (decision_values) {
    result$decision_values <- fit_dv
  }
  result
}


#' Train and Test a KTSP Classifier
#'
#' Convenience function that trains a KTSP classifier on training data and
#' evaluates it on both training and test data.
#'
#' @param train_x Training data matrix (features x samples).
#' @param train_y Training labels factor.
#' @param test_x Test data matrix.
#' @param test_y Test labels factor.
#' @param classes Character vector: `c(case, control)`.
#' @param predictions Logical; include predictions?
#' @param decision_values Logical; include decision values?
#' @param ... Additional arguments passed to [swap_train_ktsp()].
#'
#' @return A list with components:
#'   - `classifier`: the trained KTSP classifier
#'   - `train`: training set statistics
#'   - `test`: test set statistics
#'   - `train_roc`, `test_roc`: ROC curve data frames
#'   - `train_predictions`, `test_predictions` (if requested)
#'   - `train_decision_values`, `test_decision_values` (if requested)
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' data(testingData)
#' result <- swap_train_test_results(
#'   matTraining, trainingGroup,
#'   matTesting, testingGroup
#' )
#' result$test
swap_train_test_results <- function(
  train_x,
  train_y,
  test_x,
  test_y,
  classes = NULL,
  predictions = FALSE,
  decision_values = FALSE,
  ...
) {
  if (is.null(classes)) {
    classes <- rev(levels(factor(train_y)))
  }

  fit <- swap_train_ktsp(train_x, train_y, ...)

  train_result <- swap_ktsp_result(
    fit,
    train_x,
    train_y,
    classes,
    predictions,
    decision_values
  )
  test_result <- swap_ktsp_result(
    fit,
    test_x,
    test_y,
    classes,
    predictions,
    decision_values
  )

  result <- list(
    classifier = fit,
    train = train_result$stats,
    test = test_result$stats,
    train_roc = train_result$roc,
    test_roc = test_result$roc
  )

  if (predictions) {
    result$train_predictions <- train_result$predictions
    result$test_predictions <- test_result$predictions
  }
  if (decision_values) {
    result$train_decision_values <- train_result$decision_values
    result$test_decision_values <- test_result$decision_values
  }

  result
}

# Internal helper used by k_by_measurement
get_ktsp_result <- function(
  fit,
  x,
  y,
  classes = NULL,
  predictions = FALSE,
  decision_values = FALSE
) {
  swap_ktsp_result(fit, x, y, classes, predictions, decision_values)
}
