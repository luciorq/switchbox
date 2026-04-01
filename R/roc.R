#' Compute ROC Curve
#'
#' Computes the Receiver Operating Characteristic (ROC) curve from predicted
#' decision values and true binary labels.
#'
#' @param truth A factor or vector of true class labels with exactly two levels.
#' @param decision_values A numeric vector of decision/score values. Higher
#'   values should correspond to the positive class (second level of `truth`).
#' @param classes Character vector of length 2: `c(positive, negative)`. If
#'   `NULL`, levels of `truth` are used (second level = positive).
#'
#' @returns A data frame with columns:
#'   - `threshold`: the decision threshold
#'   - `tpr`: true positive rate (sensitivity)
#'   - `fpr`: false positive rate (1 - specificity)
#'
#' @noRd
compute_roc <- function(truth, decision_values, classes = NULL) {
  if (is.null(classes)) {
    lvls <- levels(factor(truth))
    classes <- c(lvls[2], lvls[1])
  }

  positive <- classes[1]
  negative <- classes[2]

  is_positive <- truth == positive
  is_negative <- truth == negative

  n_pos <- sum(is_positive)
  n_neg <- sum(is_negative)

  if (n_pos == 0L || n_neg == 0L) {
    return(data.frame(threshold = NA_real_, tpr = NA_real_, fpr = NA_real_))
  }

  # Sort by decreasing decision value
  ord <- order(decision_values, decreasing = TRUE)
  sorted_truth <- is_positive[ord]
  sorted_values <- decision_values[ord]

  # Compute cumulative TP and FP
  tp_cumsum <- cumsum(sorted_truth)
  fp_cumsum <- cumsum(!sorted_truth)

  tpr <- tp_cumsum / n_pos
  fpr <- fp_cumsum / n_neg

  # Add origin point (0, 0)
  tpr <- c(0, tpr)
  fpr <- c(0, fpr)
  thresholds <- c(Inf, sorted_values)

  data.frame(
    threshold = thresholds,
    tpr = tpr,
    fpr = fpr
  )
}


#' Compute Area Under the ROC Curve (AUC)
#'
#' Uses the trapezoidal rule to compute the AUC from ROC data or directly
#' from predictions and truth labels.
#'
#' @param truth A factor or vector of true class labels.
#' @param decision_values A numeric vector of decision values.
#' @param classes Character vector of length 2: `c(positive, negative)`.
#' @param roc_data Precomputed ROC data frame (from `compute_roc()`). If
#'   provided, `truth`, `decision_values`, and `classes` are ignored.
#'
#' @returns A single numeric value representing the AUC (between 0 and 1).
#'
#' @noRd
compute_auc <- function(
  truth = NULL,
  decision_values = NULL,
  classes = NULL,
  roc_data = NULL
) {
  if (is.null(roc_data)) {
    roc_data <- compute_roc(truth, decision_values, classes)
  }

  if (nrow(roc_data) < 2L || any(is.na(roc_data$tpr))) {
    return(NA_real_)
  }

  # Trapezoidal rule
  n <- nrow(roc_data)
  dx <- diff(roc_data$fpr)
  avg_y <- (roc_data$tpr[-n] + roc_data$tpr[-1]) / 2
  abs(sum(dx * avg_y))
}
