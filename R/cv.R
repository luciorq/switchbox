#' K-Fold Cross-Validation for KTSP
#'
#' Performs K-fold cross-validation by training a KTSP classifier on K-1 folds
#' and evaluating on the held-out fold, then combining results.
#'
#' @param x Numeric matrix (features x samples).
#' @param y Factor with two levels.
#' @param classes Character vector: `c(case, control)`. If `NULL`, derived
#'   from `levels(y)`.
#' @param k Integer; number of folds. Default is `4`.
#' @param folds An optional pre-specified list of index vectors (one per fold).
#'   If provided, `k` is ignored.
#' @param randomize Logical; randomly permute samples before splitting?
#'   Default is `TRUE`.
#' @param ... Additional arguments passed to [swap_train_ktsp()].
#'
#' @returns A list with:
#'   - `cv`: list of per-fold results
#'   - `folds`: the fold index lists used
#'   - `predictions`: combined predictions across all folds
#'   - `decision_values`: combined decision values
#'   - `stats`: overall prediction statistics
#'   - `truth`: the true labels in fold order
#'   - `randomized_indices`: permutation indices used (if randomized)
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' set.seed(42)
#' cv_result <- swap_ktsp_cv(matTraining, trainingGroup, k = 3)
#' cv_result$stats
swap_ktsp_cv <- function(
  x,
  y,
  classes = NULL,
  k = 4,
  folds = NULL,
  randomize = TRUE,
  ...
) {
  ix <- seq_len(ncol(x))

  if (is.null(folds)) {
    if (randomize) {
      ix <- sample(ix)
    }
    x <- x[, ix, drop = FALSE]
    y <- y[ix]
    folds <- swap_get_kfold_indices(y, k)
  } else {
    k <- length(folds)
  }

  if (is.null(classes)) {
    classes <- rev(levels(factor(y)))
  }

  cv <- vector("list", k)
  for (i in seq_len(k)) {
    cv[[i]] <- swap_train_test_results(
      train_x = x[, -folds[[i]], drop = FALSE],
      train_y = y[-folds[[i]]],
      test_x = x[, folds[[i]], drop = FALSE],
      test_y = y[folds[[i]]],
      classes = classes,
      predictions = TRUE,
      decision_values = TRUE,
      ...
    )
  }

  truth <- y[unlist(folds)]
  preds <- unlist(lapply(cv, function(z) z$test_predictions))
  dvs <- unlist(lapply(cv, function(z) {
    z$test_decision_values / nrow(z$classifier$TSPs)
  }))

  stats_result <- swap_prediction_stats(
    predictions = preds,
    truth = truth,
    classes = classes,
    decision_values = dvs
  )

  list(
    cv = cv,
    folds = folds,
    predictions = preds,
    decision_values = dvs,
    stats = stats_result,
    truth = truth,
    randomized_indices = ix
  )
}


#' Leave-One-Out Cross-Validation for KTSP
#'
#' Performs leave-one-out cross-validation: each sample is held out once while
#' a KTSP classifier is trained on the remaining samples.
#'
#' @param x Numeric matrix (features x samples).
#' @param y Factor with two levels.
#' @param classes Character vector: `c(case, control)`. If `NULL`, derived
#'   from `levels(y)`.
#' @param ... Additional arguments passed to [swap_train_ktsp()].
#'
#' @returns A list with:
#'   - `loo`: list of per-sample results
#'   - `predictions`: combined predictions
#'   - `decision_values`: combined decision values
#'   - `stats`: overall prediction statistics
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' # LOO is slow; using a subset for the example
#' small_x <- matTraining[1:20, 1:20]
#' small_y <- trainingGroup[1:20]
#' loo_result <- swap_ktsp_loo(small_x, small_y)
#' loo_result$stats
swap_ktsp_loo <- function(x, y, classes = NULL, ...) {
  if (is.null(classes)) {
    classes <- rev(levels(factor(y)))
  }

  n <- ncol(x)
  loo <- vector("list", n)

  for (i in seq_len(n)) {
    loo[[i]] <- swap_train_test_results(
      train_x = x[, -i, drop = FALSE],
      train_y = y[-i],
      test_x = x[, i, drop = FALSE],
      test_y = y[i],
      classes = classes,
      predictions = TRUE,
      decision_values = TRUE,
      ...
    )
  }

  preds <- unlist(lapply(loo, function(z) z$test_predictions))
  dvs <- unlist(lapply(loo, function(z) {
    z$test_decision_values / nrow(z$classifier$TSPs)
  }))

  stats_result <- swap_prediction_stats(
    predictions = preds,
    truth = y,
    classes = classes,
    decision_values = dvs
  )

  list(
    loo = loo,
    predictions = preds,
    decision_values = dvs,
    stats = stats_result
  )
}


#' Generate Stratified K-Fold Indices
#'
#' Creates a list of index vectors for K-fold cross-validation, maintaining
#' approximate class balance across folds.
#'
#' @param y A factor or vector of class labels.
#' @param k Integer; number of folds. Default is `4`.
#'
#' @returns A list of length `k`, where each element is an integer vector of
#'   sample indices for that fold.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' folds <- swap_get_kfold_indices(trainingGroup, k = 5)
#' lengths(folds)
swap_get_kfold_indices <- function(y, k = 4) {
  folds <- vector("list", k)

  for (label_i in lapply(unique(y), function(yi) sample(which(y == yi)))) {
    fold_assignment <- ((seq_along(label_i) - 1L) %% k) + 1L
    for (j in seq_len(k)) {
      folds[[j]] <- c(folds[[j]], label_i[fold_assignment == j])
    }
  }

  folds
}
