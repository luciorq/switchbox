#' Train a K-Top-Scoring-Pair (KTSP) Classifier
#'
#' Trains a binary classifier based on K disjoint top-scoring pairs of
#' features. The number of pairs K is selected automatically from `k_range`
#' using a selection method (by default, a t-test-based criterion).
#'
#' @param x A numeric matrix with features in rows and samples in columns.
#'   Must have rownames.
#' @param y A factor with exactly two levels indicating the class of each
#'   sample.
#' @param classes Character vector of length 2: `c(case, control)`. If `NULL`,
#'   derived from `levels(y)`.
#' @param k_range Integer vector specifying the range of K values to consider.
#'   Default is `2:10`.
#' @param filter_fn Feature filtering function. Default is
#'   [swap_filter_wilcoxon()]. Set to `NULL` to use all features.
#' @param restricted_pairs An optional two-column character matrix specifying
#'   candidate feature pairs.
#' @param handle_ties Logical; use tie-handling scoring? Default is `FALSE`.
#' @param disjoint Logical; require feature pairs to be disjoint? Default is
#'   `TRUE`.
#' @param k_selection_fn Function to select K. Default is [swap_k_by_ttest()].
#'   Alternatives include [swap_k_by_measurement()].
#' @param k_opts List of additional options passed to `k_selection_fn`.
#' @param score_fn Scoring function. Default is
#'   [swap_calculate_signed_tsp_scores()].
#' @param score_opts Additional options passed to `score_fn`.
#' @param verbose Logical; print progress messages?
#' @param ... Additional arguments passed to `filter_fn`.
#'
#' @returns A list (the classifier) with components:
#'   - `name`: character string like `"5TSPS"`
#'   - `TSPs`: character matrix with K rows and 2 columns (gene1, gene2)
#'   - `score`: numeric vector of pair scores
#'   - `tie_vote`: factor indicating tie handling per pair
#'   - `labels`: character vector of class labels
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:8)
#' classifier$TSPs
swap_train_ktsp <- function(
  x,
  y,
  classes = NULL,
  k_range = 2:10,
  filter_fn = swap_filter_wilcoxon,
  restricted_pairs = NULL,
  handle_ties = FALSE,
  disjoint = TRUE,
  k_selection_fn = swap_k_by_ttest,
  k_opts = list(),
  score_fn = swap_calculate_signed_tsp_scores,
  score_opts = NULL,
  verbose = FALSE,
  ...
) {
  if (is.null(classes)) {
    classes <- rev(levels(y))
  }

  if (verbose) {
    message(sprintf(
      "Selecting %s (case) and %s (control) as class labels.",
      classes[1],
      classes[2]
    ))
  }

  # Calculate pair scores
  s <- swap_calculate_scores(
    x,
    y,
    classes,
    filter_fn,
    restricted_pairs,
    handle_ties,
    verbose,
    score_fn,
    score_opts,
    ...
  )

  max_k <- max(k_range)
  tsp_table <- make_tsp_table(s, max_k, disjoint)

  # Select final TSPs
  sel <- k_selection_fn(x, y, tsp_table, classes, k_range, k_opts)

  # Prepare output
  classifier <- list(name = sprintf("%dTSPS", length(sel)))
  classifier$TSPs <- as.matrix(
    tsp_table[sel, c("gene1", "gene2")],
    length(sel),
    2
  )
  classifier$score <- tsp_table$score[sel]
  classifier$tie_vote <- factor(
    tsp_table$tie_vote[sel],
    levels = 0:2,
    labels = c("both", rev(classes))
  )
  classifier$labels <- rev(classes)

  names(classifier$score) <- rownames(classifier$TSPs)
  names(classifier$tie_vote) <- rownames(classifier$TSPs)

  classifier
}


#' Train a Single Top-Scoring-Pair (1-TSP) Classifier
#'
#' Trains a binary classifier using the single best top-scoring pair of
#' features.
#'
#' @inheritParams swap_train_ktsp
#'
#' @returns A list (the classifier) with the same structure as
#'   [swap_train_ktsp()] but with exactly one pair.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' classifier <- swap_train_1tsp(matTraining, trainingGroup)
#' classifier$TSPs
swap_train_1tsp <- function(
  x,
  y,
  classes = NULL,
  filter_fn = swap_filter_wilcoxon,
  restricted_pairs = NULL,
  handle_ties = FALSE,
  disjoint = TRUE,
  score_fn = swap_calculate_signed_tsp_scores,
  score_opts = NULL,
  verbose = FALSE,
  ...
) {
  # Calculate pair scores
  s <- swap_calculate_scores(
    x,
    y,
    classes,
    filter_fn,
    restricted_pairs,
    handle_ties,
    verbose,
    score_fn,
    score_opts,
    ...
  )

  if (is.null(classes)) {
    classes <- rev(levels(y))
    if (verbose) {
      message(sprintf(
        "Selecting %s (case) and %s (control) as class labels.",
        classes[1],
        classes[2]
      ))
    }
  }

  j <- which.max(abs(s$score))[1]
  pair <- unlist(strsplit(names(s$score)[j], ","))
  score <- s$score[j]
  is_reversed <- score > 0 && isTRUE(s$signed)
  if (is_reversed) {
    pair <- rev(pair)
  }

  # Flip tie vote when reversing pair order
  tie_vote_change <- c(0L, 2L, 1L)
  ties <- if (is_reversed) {
    tie_vote_change[s$tie_vote[j] + 1L]
  } else {
    s$tie_vote[j]
  }
  score <- abs(score)

  tie_vote <- factor(ties, levels = 0:2, labels = c("both", s$labels))

  list(
    name = "1TSP",
    TSPs = matrix(pair, nrow = 1),
    score = score,
    labels = s$labels,
    tie_vote = tie_vote
  )
}


#' Select K by T-Test
#'
#' Selects the optimal number of top-scoring pairs (K) by maximizing the
#' t-test statistic between the KTSP vote distributions of the two classes.
#'
#' @param x Numeric matrix (features x samples).
#' @param y Factor with two levels.
#' @param score_table Data frame of scored pairs (from [swap_make_tsp_table()]).
#' @param classes Character vector of length 2: `c(case, control)`.
#' @param k_range Integer vector of candidate K values.
#' @param k_opts List of additional options (currently unused).
#'
#' @returns Integer vector of selected row indices into `score_table`.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' scores <- swap_calculate_scores(matTraining, trainingGroup)
#' tbl <- swap_make_tsp_table(scores, maxk = 10)
#' sel <- swap_k_by_ttest(
#'   matTraining, trainingGroup, tbl,
#'   rev(levels(trainingGroup)), 2:8, list()
#' )
swap_k_by_ttest <- function(
  x,
  y,
  score_table,
  classes,
  k_range,
  k_opts = list()
) {
  max_tsps <- nrow(score_table)
  s0 <- which(y == classes[2])
  s1 <- which(y == classes[1])

  best_tt <- -Inf
  best_k <- 1L
  classifier_test <- list(labels = rev(classes))

  for (k in seq_along(k_range)) {
    nk <- min(k_range[k], max_tsps)
    classifier_test$TSPs <- as.matrix(
      score_table[seq_len(nk), c("gene1", "gene2")],
      nrow = nk,
      ncol = 2
    )
    classifier_test$tie_vote <- score_table$tie_vote[seq_len(nk)]

    stat0 <- swap_ktsp_statistics(
      x[, s0, drop = FALSE],
      classifier_test
    )$statistics
    stat1 <- swap_ktsp_statistics(
      x[, s1, drop = FALSE],
      classifier_test
    )$statistics

    tt <- abs(mean(stat0) - mean(stat1)) /
      sqrt(stats::var(stat1) + stats::var(stat0) + 1e-9)

    if (abs(best_tt - tt) > 1e-7 && best_tt < tt) {
      best_tt <- tt
      best_k <- k
    }
  }

  seq_len(k_range[best_k])
}


#' Select K by Performance Measurement
#'
#' Selects the optimal K by maximizing a classification performance metric
#' (e.g. AUC, accuracy) on the training data.
#'
#' @param x Numeric matrix (features x samples).
#' @param y Factor with two levels.
#' @param score_table Data frame of scored pairs.
#' @param classes Character vector: `c(case, control)`.
#' @param k_range Integer vector of candidate K values.
#' @param k_opts List with optional elements:
#'   - `disjoint`: logical (default `TRUE`)
#'   - `measurement`: character, one of `"accuracy"`, `"sensitivity"`,
#'     `"specificity"`, `"balanced_accuracy"`, `"auc"` (default `"auc"`)
#'
#' @returns Integer vector of selected row indices.
#'
#' @export
swap_k_by_measurement <- function(
  x,
  y,
  score_table,
  classes,
  k_range,
  k_opts = list(
    disjoint = TRUE,
    measurement = "auc"
  )
) {
  disjoint <- k_opts$disjoint %||% TRUE
  measurement <- k_opts$measurement %||% "auc"

  valid <- c(
    "accuracy",
    "sensitivity",
    "specificity",
    "balanced_accuracy",
    "auc"
  )
  if (!measurement %in% valid) {
    stop(
      sprintf(
        "'measurement' must be one of: %s",
        paste(valid, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  max_k <- max(k_range)
  fit <- list(
    TSPs = as.matrix(score_table[1, 1:2], 1, 2),
    tie_vote = score_table[1, 4],
    labels = rev(classes)
  )

  ix <- 1L
  stat_values <- numeric(0)

  for (i in seq_len(min(max_k, nrow(score_table)))) {
    if (i > 1L) {
      if (
        disjoint &&
          (sum(
            as.matrix(score_table[i, 1:2]) %in%
              unique(as.vector(fit$TSPs))
          ) >
            0)
      ) {
        next
      }
      fit$TSPs <- as.matrix(
        rbind(fit$TSPs, score_table[i, 1:2]),
        length(ix) + 1,
        2
      )
      fit$tie_vote <- c(fit$tie_vote, score_table[i, 4])
      ix <- c(ix, i)
    }

    result <- get_ktsp_result(
      fit = fit,
      x = x,
      y = y,
      classes = classes
    )
    stat_values <- c(stat_values, result$stats[measurement])
  }

  ix[seq_len(which.max(stat_values))]
}


# Build sorted TSP table from scores
#
# @param scores List from swap_calculate_scores
# @param maxk Maximum number of pairs to include
# @param disjoint Require disjoint pairs?
# @return data.frame with columns gene1, gene2, score, tie_vote
#
# @noRd
make_tsp_table <- function(scores, maxk, disjoint = TRUE) {
  tie_vote_change <- c(0L, 2L, 1L)
  abs_score <- abs(scores$score)

  ix <- order(abs_score, decreasing = TRUE)

  gene1 <- character(0)
  gene2 <- character(0)
  score_vals <- numeric(0)
  tie_votes <- integer(0)

  for (i in seq_along(ix)) {
    j <- ix[i]
    pair <- unlist(strsplit(names(scores$score)[j], ","))

    # Skip if reverse pair exists
    if (
      (pair[1] %in% gene1 && pair[2] %in% gene2) ||
        (pair[1] %in% gene2 && pair[2] %in% gene1)
    ) {
      next
    }

    if (disjoint && (sum(pair %in% unique(c(gene1, gene2))) > 0)) {
      next
    }

    # If score positive and signed, reverse pair order
    is_reversed <- scores$score[j] > 0 && isTRUE(scores$signed)
    score_vals <- c(score_vals, abs_score[j])
    gene1 <- c(gene1, if (is_reversed) pair[2] else pair[1])
    gene2 <- c(gene2, if (is_reversed) pair[1] else pair[2])

    tv <- if (is_reversed) {
      tie_vote_change[scores$tie_vote[j] + 1L]
    } else {
      scores$tie_vote[j]
    }
    tie_votes <- c(tie_votes, tv)

    if (length(gene1) == maxk) break
  }

  data.frame(
    gene1 = gene1,
    gene2 = gene2,
    score = score_vals,
    tie_vote = tie_votes,
    stringsAsFactors = FALSE
  )
}


#' Create a TSP Score Table
#'
#' Assembles a sorted data frame of top-scoring feature pairs from computed
#' scores.
#'
#' @param scores A list as returned by [swap_calculate_scores()], containing
#'   `score`, `tie_vote`, and `signed` components.
#' @param maxk Integer; maximum number of pairs to include.
#' @param disjoint Logical; require pairs to be disjoint (no shared features)?
#'   Default is `TRUE`.
#'
#' @returns A data frame with columns `gene1`, `gene2`, `score`, and
#'   `tie_vote`, sorted by descending absolute score.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' scores <- swap_calculate_scores(matTraining, trainingGroup)
#' tbl <- swap_make_tsp_table(scores, maxk = 10)
#' head(tbl)
swap_make_tsp_table <- function(scores, maxk, disjoint = TRUE) {
  make_tsp_table(scores, maxk, disjoint)
}

# Null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x
