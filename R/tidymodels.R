#' Convert a Data Frame to a switchbox Matrix
#'
#' Converts a data frame or tibble with samples in rows and features in
#' columns (the tidymodels convention) to a numeric matrix with features in
#' rows and samples in columns (the switchbox convention).
#'
#' @param data A data frame, tibble, or matrix with samples in rows and
#'   numeric features in columns.
#' @param outcome Optional; a character string naming the column that contains
#'   the class labels. If provided, that column is excluded from the matrix and
#'   returned separately as a factor.
#'
#' @returns If `outcome` is `NULL`, a numeric matrix with features in rows and
#'   samples in columns. If `outcome` is provided, a list with components:
#'   - `x`: the numeric matrix (features x samples)
#'   - `y`: a factor of class labels
#'
#' @export
#'
#' @examples
#' df <- data.frame(
#'   class = factor(c("a", "b", "a", "b")),
#'   gene1 = c(1.2, 3.4, 1.5, 3.1),
#'   gene2 = c(5.6, 2.1, 5.9, 2.3)
#' )
#' result <- swap_from_tidy(df, outcome = "class")
#' result$x
#' result$y
swap_from_tidy <- function(data, outcome = NULL) {
  if (is.matrix(data)) {
    if (!is.numeric(data)) {
      stop("Matrix must be numeric.", call. = FALSE)
    }
    return(t(data))
  }

  if (!is.null(outcome)) {
    if (!outcome %in% names(data)) {
      stop(
        sprintf("Column '%s' not found in data.", outcome),
        call. = FALSE
      )
    }
    y <- factor(data[[outcome]])
    data <- data[, setdiff(names(data), outcome), drop = FALSE]
    x <- t(as.matrix(data))
    return(list(x = x, y = y))
  }

  t(as.matrix(data))
}


#' Convert a switchbox Matrix to a Tibble
#'
#' Converts a numeric matrix with features in rows and samples in columns
#' (the switchbox convention) to a data frame with samples in rows and
#' features in columns (the tidymodels convention).
#'
#' @param x A numeric matrix with features in rows and samples in columns.
#' @param y Optional; a factor or vector of class labels. If provided, it is
#'   added as a column named by `outcome`.
#' @param outcome Character; column name for the class labels. Default is
#'   `"class"`.
#'
#' @returns A data frame with samples in rows. If `y` is provided, it includes
#'   a factor column for the outcome.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' tidy_df <- swap_to_tidy(matTraining, trainingGroup, outcome = "prognosis")
#' head(tidy_df[, 1:5])
swap_to_tidy <- function(x, y = NULL, outcome = "class") {
  df <- as.data.frame(t(x))

  if (!is.null(y)) {
    df[[outcome]] <- factor(y)
    # Move outcome to first column
    df <- df[, c(outcome, setdiff(names(df), outcome)), drop = FALSE]
  }

  df
}


#' Evaluate a KTSP Classifier and Return a Tidy Results Tibble
#'
#' Classifies samples using a KTSP classifier and returns results in a format
#' compatible with yardstick metric functions. The returned data frame contains
#' columns for the true labels, predicted class, and optionally the normalized
#' decision value.
#'
#' @param classifier A KTSP classifier as returned by [swap_train_ktsp()] or
#'   [swap_train_1tsp()].
#' @param x A numeric matrix (features x samples).
#' @param y A factor of true class labels.
#'
#' @returns A data frame with columns:
#'   - `truth`: factor of true class labels
#'   - `.pred_class`: factor of predicted class labels
#'   - `.decision_value`: numeric decision value (voting statistic normalized
#'     by the number of pairs)
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:5)
#' results <- swap_tidy_result(classifier, matTraining, trainingGroup)
#' head(results)
swap_tidy_result <- function(classifier, x, y) {
  preds <- swap_ktsp_classify(x, classifier)
  stats <- swap_ktsp_statistics(x, classifier)$statistics
  n_pairs <- nrow(classifier$TSPs)

  data.frame(
    truth = factor(y, levels = levels(preds)),
    .pred_class = preds,
    .decision_value = stats / n_pairs,
    row.names = NULL
  )
}
