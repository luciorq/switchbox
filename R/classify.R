#' Compute KTSP Statistics (Votes)
#'
#' Computes the per-sample KTSP voting statistics by comparing the feature
#' values within each top-scoring pair. The default statistic is the sum
#' of votes minus K/2 (majority-wins principle).
#'
#' @param x A numeric matrix (features x samples). Must contain the features
#'   referenced in `classifier$TSPs`.
#' @param classifier A KTSP classifier as returned by [swap_train_ktsp()] or
#'   [swap_train_1tsp()].
#' @param combine_fn An optional function to aggregate the per-TSP votes. If
#'   not provided, uses the default signed voting statistic.
#'
#' @returns A list with:
#'   - `statistics`: named numeric vector of aggregate per-sample statistics
#'   - `comparisons`: matrix (samples x pairs) of individual TSP comparisons
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:8)
#' stats <- swap_ktsp_statistics(matTraining, classifier)
#' head(stats$statistics)
swap_ktsp_statistics <- function(x, classifier, combine_fn) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("'x' must be a numeric matrix.", call. = FALSE)
  }

  # Handle single-sample case
  if (ncol(x) == 1L) {
    x <- as.matrix(x)
  }

  tsps <- classifier$TSPs
  if (!is.matrix(tsps)) {
    tsps <- matrix(tsps, nrow = 1)
  }
  n_pairs <- nrow(tsps)

  has_tie_vote <- !is.null(classifier$tie_vote) &&
    !all(classifier$tie_vote == 0) &&
    !all(classifier$tie_vote == "both" | classifier$tie_vote == 0)

  if (has_tie_vote) {
    # Tie-aware comparisons
    comparisons <- matrix(nrow = ncol(x), ncol = n_pairs)
    comparison_names <- character(n_pairs)
    rownames(comparisons) <- colnames(x)

    for (i in seq_len(n_pairs)) {
      g1 <- x[tsps[i, 1], ]
      g2 <- x[tsps[i, 2], ]

      tv <- as.character(classifier$tie_vote[i])
      if (tv == "both") {
        # Equality counted as 0.5 (FIXED: operator precedence bug)
        comparisons[, i] <- (g1 > g2) + 0.5 * (g1 == g2)
        comparison_names[i] <- sprintf("%s>%s", tsps[i, 1], tsps[i, 2])
      } else if (tv == classifier$labels[2]) {
        # Equality counted in favour of class 1
        comparisons[, i] <- as.numeric(g1 > g2)
        comparison_names[i] <- sprintf("%s>=%s", tsps[i, 1], tsps[i, 2])
      } else {
        # Equality counted in favour of class 0
        comparisons[, i] <- as.numeric(g1 >= g2)
        comparison_names[i] <- sprintf("%s>>%s", tsps[i, 1], tsps[i, 2])
      }
    }
    colnames(comparisons) <- comparison_names

    if (missing(combine_fn)) {
      ktsp_stat <- rowSums(comparisons) - n_pairs / 2
      names(ktsp_stat) <- colnames(x)
    } else {
      ktsp_stat <- apply(comparisons, 1, combine_fn)
    }
  } else {
    # Standard (no-tie) comparisons
    comparisons <- t(
      x[tsps[, 1], , drop = FALSE] > x[tsps[, 2], , drop = FALSE]
    )

    if (n_pairs == 1L) {
      comparisons <- matrix(comparisons, ncol = 1)
    }

    colnames(comparisons) <- vapply(
      seq_len(n_pairs),
      function(i) sprintf("%s>%s", tsps[i, 1], tsps[i, 2]),
      character(1)
    )
    rownames(comparisons) <- colnames(x)

    stats1 <- x[tsps[, 1], , drop = FALSE] > x[tsps[, 2], , drop = FALSE]
    stats2 <- x[tsps[, 1], , drop = FALSE] < x[tsps[, 2], , drop = FALSE]

    if (n_pairs == 1L) {
      stats1 <- t(matrix(stats1, ncol = 1))
      stats2 <- t(matrix(stats2, ncol = 1))
    }

    if (missing(combine_fn)) {
      s1 <- colSums(as.matrix(stats1))
      s2 <- colSums(as.matrix(stats2))
      ktsp_stat <- s1 - s2
      names(ktsp_stat) <- colnames(x)
    } else {
      ktsp_stat <- apply(stats1, 2, combine_fn)
    }
  }

  list(statistics = ktsp_stat, comparisons = comparisons)
}


#' Classify Samples Using a KTSP Classifier
#'
#' Predicts the class of each sample by applying the KTSP voting rule. By
#' default, uses the majority-wins principle (statistic > 0).
#'
#' @param x A numeric matrix (features x samples).
#' @param classifier A KTSP classifier from [swap_train_ktsp()] or
#'   [swap_train_1tsp()].
#' @param decision_fn An optional function that takes the aggregate voting
#'   statistic and returns a logical vector (TRUE = second class). If not
#'   provided, uses the default `statistic > 0`.
#'
#' @returns A factor of predicted class labels, with levels matching
#'   `classifier$labels`.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' data(testingData)
#' classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:8)
#' predictions <- swap_ktsp_classify(matTesting, classifier)
#' table(predictions, testingGroup)
swap_ktsp_classify <- function(x, classifier, decision_fn) {
  labels <- classifier$labels

  if (missing(decision_fn)) {
    ktsp_stat <- swap_ktsp_statistics(x, classifier)$statistics > 0
  } else {
    ktsp_stat <- swap_ktsp_statistics(x, classifier, decision_fn)$statistics
  }

  factor(
    ifelse(ktsp_stat, labels[[2]], labels[[1]]),
    levels = labels
  )
}
