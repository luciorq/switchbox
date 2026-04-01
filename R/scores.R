#' Calculate Pairwise TSP Scores
#'
#' Computes pairwise Top-Scoring Pair (TSP) scores for all feature pairs or a
#' restricted set of pairs. This is the main entry point for score computation.
#'
#' @param x A numeric matrix with features (e.g. genes) in rows and samples in
#'   columns. Must have rownames.
#' @param y A factor with exactly two levels indicating the class of each
#'   sample.
#' @param classes Character vector of length 2 specifying `c(case, control)`.
#'   If `NULL`, derived from `levels(y)`.
#' @param filter_fn A function for feature filtering, or `NULL` to use all
#'   features. The function should accept `(y, x, ...)` and return a character
#'   vector of feature names. Default is [swap_filter_wilcoxon()].
#' @param restricted_pairs An optional two-column character matrix specifying
#'   candidate feature pairs to evaluate. Both columns must contain feature
#'   names present in `rownames(x)`.
#' @param handle_ties Logical; if `TRUE`, uses an extended scoring method that
#'   explicitly handles tied feature values.
#' @param verbose Logical; if `TRUE`, prints progress messages.
#' @param score_fn A function that computes pairwise scores. Must have the
#'   signature `(y, x1, x2, classes, restricted_pairs, handle_ties, verbose,
#'   score_opts)`. Default is [swap_calculate_signed_tsp_scores()].
#' @param score_opts A list of additional options passed to `score_fn`.
#' @param ... Additional arguments passed to `filter_fn`.
#'
#' @returns A list with components:
#'   - `score`: named numeric vector of pair scores (names are `"gene1,gene2"`)
#'   - `labels`: character vector of length 2 with class labels
#'   - `tie_vote`: integer vector of tie-vote types (0 = both, 1 = class 0
#'     favoured, 2 = class 1 favoured)
#'   - `signed`: logical indicating whether scores are signed
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' scores <- swap_calculate_scores(matTraining, trainingGroup)
#' head(scores$score)
swap_calculate_scores <- function(
  x,
  y,
  classes = NULL,
  filter_fn = swap_filter_wilcoxon,
  restricted_pairs = NULL,
  handle_ties = FALSE,
  verbose = FALSE,
  score_fn = swap_calculate_signed_tsp_scores,
  score_opts = list(),
  ...
) {
  check_input(y = y, x = x, restricted_pairs = restricted_pairs)

  formatted <- format_input(
    y = y,
    x = x,
    filter_fn = filter_fn,
    restricted_pairs = restricted_pairs,
    verbose = verbose,
    ...
  )

  scores <- score_fn(
    y = y,
    x1 = formatted$x,
    x2 = formatted$x,
    classes = classes,
    restricted_pairs = formatted$filtered_pairs,
    handle_ties = handle_ties,
    verbose = verbose,
    score_opts = score_opts
  )

  if (is.null(scores$signed)) {
    scores$signed <- FALSE
  }

  scores
}


#' Calculate Signed TSP Scores
#'
#' Computes the signed Top-Scoring Pair scores using an efficient C
#' implementation. The signed score for a pair (i, j) is:
#' `P(Xi < Xj | class1) - P(Xi < Xj | class0) - P(Xi > Xj | class1) +
#' P(Xi > Xj | class0)`.
#'
#' @param y A factor with exactly two levels.
#' @param x1 A numeric matrix (features x samples).
#' @param x2 A numeric matrix (features x samples). Defaults to `x1`.
#' @param classes Character vector of length 2: `c(case, control)`.
#' @param restricted_pairs Optional two-column character matrix.
#' @param handle_ties Logical; use tie-handling scoring method?
#' @param verbose Logical; print progress messages?
#' @param score_opts Additional options (unused, for API compatibility).
#'
#' @returns A list with `score`, `labels`, `tie_vote`, and `signed` components.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' scores <- swap_calculate_signed_tsp_scores(
#'   trainingGroup, matTraining
#' )
#' head(scores$score)
swap_calculate_signed_tsp_scores <- function(
  y,
  x1,
  x2 = NULL,
  classes = NULL,
  restricted_pairs = NULL,
  handle_ties = FALSE,
  verbose = FALSE,
  score_opts = list()
) {
  if (is.null(x2)) {
    x2 <- x1
  }

  scores <- calculate_signed_score(
    y = y,
    x1 = x1,
    x2 = x2,
    restricted_pairs = restricted_pairs,
    handle_ties = handle_ties,
    verbose = verbose,
    classes = classes
  )

  # Convert matrix scores to vector for unrestricted case
  if (is.matrix(scores$score)) {
    scores$score <- score_matrix_to_vector(scores$score)
    if (!is.null(scores$tie_vote)) {
      scores$tie_vote <- score_matrix_to_vector(scores$tie_vote)
    } else {
      scores$tie_vote <- rep(0L, length(scores$score))
    }
  } else {
    if (is.null(scores$tie_vote)) {
      scores$tie_vote <- rep(0L, length(scores$score))
    }
  }

  scores$signed <- TRUE
  scores
}


#' Calculate Basic TSP Scores
#'
#' Computes basic (unsigned) TSP scores using a pure R implementation.
#' The basic score for pair (i, j) is:
#' `P(Xi < Xj | class1) + P(Xi > Xj | class0)`.
#'
#' @inheritParams swap_calculate_signed_tsp_scores
#'
#' @returns A list with `score`, `labels`, `tie_vote`, and `signed` components.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' scores <- swap_calculate_basic_tsp_scores(
#'   trainingGroup, matTraining
#' )
#' head(scores$score)
swap_calculate_basic_tsp_scores <- function(
  y,
  x1,
  x2 = NULL,
  classes = NULL,
  restricted_pairs = NULL,
  handle_ties = FALSE,
  verbose = FALSE,
  score_opts = list()
) {
  if (is.null(x2)) {
    x2 <- x1
  }

  if (is.null(classes)) {
    lvls <- levels(y)
  } else {
    lvls <- rev(classes)
  }

  scores <- list(labels = lvls)

  if (is.null(restricted_pairs)) {
    score <- t(vapply(
      rownames(x1),
      function(xi) {
        vapply(
          rownames(x2),
          function(xj) {
            calculate_tsp_instance(x1[xi, ], x2[xj, ], y, lvls)
          },
          numeric(1)
        )
      },
      numeric(nrow(x2))
    ))

    rownames(score) <- rownames(x1)
    colnames(score) <- rownames(x2)
    scores$score <- score_matrix_to_vector(score)
  } else {
    score <- apply(restricted_pairs, 1, function(p) {
      calculate_tsp_instance(x1[p[1], ], x2[p[2], ], y, lvls)
    })
    names(score) <- apply(restricted_pairs, 1, function(p) {
      sprintf("%s,%s", p[1], p[2])
    })
    scores$score <- score
  }

  if (verbose && handle_ties) {
    message("Tie handling not available for basic TSP score.")
  }

  scores$tie_vote <- rep(0L, length(scores$score))
  scores$signed <- FALSE
  scores
}


# Internal: Core C-backed signed score computation
#
# This is the workhorse that calls the C routines.
#
# @noRd
calculate_signed_score <- function(
  y,
  x1,
  x2,
  restricted_pairs = NULL,
  handle_ties = FALSE,
  verbose = FALSE,
  classes = NULL
) {
  n <- length(y)
  m1 <- nrow(x1)
  m2 <- nrow(x2)

  # Rank the combined data
  data_tied <- apply(rbind(x1, x2), 2, rank)
  data1_tied <- data_tied[seq_len(m1), , drop = FALSE]
  data2_tied <- data_tied[(m1 + 1):(m1 + m2), , drop = FALSE]

  # Set up class labels
  if (is.null(classes)) {
    labels <- levels(y)
  } else {
    labels <- rev(classes)
  }

  situation_v <- as.vector(y)

  if (handle_ties) {
    if (is.null(restricted_pairs)) {
      # Unrestricted + ties
      if (verbose) {
        message(sprintf(
          "Computing scores for %d features (%s pairs).",
          m1,
          formatC(m1 * (m1 - 1) / 2)
        ))
      }

      d <- .C(
        "CalculateSignedScoreCoreTieHandler",
        as.integer(as.numeric(situation_v == labels[1])),
        as.integer(n),
        as.double(data1_tied),
        as.integer(m1),
        as.double(data2_tied),
        as.integer(m2),
        as.double(matrix(0, m1, m2)),
        as.double(matrix(0, m1, m2)),
        as.double(matrix(0, m1, m2)),
        as.double(matrix(0, m1, m2))
      )

      score <- matrix(d[[7]], nrow = m1)
      tie_vote <- matrix(d[[10]], nrow = m1)

      names1 <- rownames(x1)
      names2 <- rownames(x2)
      rownames(score) <- names1
      colnames(score) <- names2
      rownames(tie_vote) <- names1
      colnames(tie_vote) <- names2

      list(score = score, labels = labels, tie_vote = tie_vote)
    } else {
      # Restricted + ties -- FIXED: was calling wrong C function
      pairs_no <- nrow(restricted_pairs)
      pairs_ind <- matrix(0L, pairs_no, 2)
      pairs_ind[, 1] <- match(restricted_pairs[, 1], rownames(x1))
      pairs_ind[, 2] <- match(restricted_pairs[, 2], rownames(x2))

      if (verbose) {
        message(sprintf(
          "Computing scores for %d restricted pairs.",
          pairs_no
        ))
      }

      d <- .C(
        "CalculateSignedScoreRestrictedPairsCoreTieHandler",
        as.integer(as.numeric(situation_v == labels[1])),
        as.integer(n),
        as.double(data1_tied),
        as.integer(m1),
        as.double(data2_tied),
        as.integer(m2),
        as.integer(pairs_ind[, 1] - 1L),
        as.integer(pairs_ind[, 2] - 1L),
        as.integer(pairs_no),
        as.double(numeric(pairs_no)),
        as.double(numeric(pairs_no)),
        as.double(numeric(pairs_no)),
        as.double(numeric(pairs_no))
      )

      score <- d[[10]]
      names(score) <- sprintf(
        "%s,%s",
        restricted_pairs[, 1],
        restricted_pairs[, 2]
      )
      tie_vote <- d[[13]]

      list(score = score / 2, labels = labels, tie_vote = tie_vote)
    }
  } else {
    if (is.null(restricted_pairs)) {
      # Unrestricted, no ties
      if (verbose) {
        message(sprintf(
          "Computing scores for %d features (%s pairs).",
          m1,
          formatC(m1 * (m1 - 1) / 2)
        ))
      }

      d <- .C(
        "CalculateSignedScoreCore",
        as.integer(as.numeric(situation_v == labels[1])),
        as.integer(n),
        as.double(data1_tied),
        as.integer(m1),
        as.double(data2_tied),
        as.integer(m2),
        as.double(matrix(0, m1, m2)),
        as.double(matrix(0, m1, m2)),
        as.double(1),
        as.double(matrix(0, m1, m2)),
        as.double(1),
        as.double(matrix(0, m1, m2)),
        as.double(matrix(0, m1, m2))
      )

      score <- matrix(d[[7]], nrow = m1)

      names1 <- rownames(x1)
      names2 <- rownames(x2)
      rownames(score) <- names1
      colnames(score) <- names2

      list(score = score / 2, labels = labels)
    } else {
      # Restricted, no ties
      pairs_no <- nrow(restricted_pairs)
      pairs_ind <- matrix(0L, pairs_no, 2)
      pairs_ind[, 1] <- match(restricted_pairs[, 1], rownames(x1))
      pairs_ind[, 2] <- match(restricted_pairs[, 2], rownames(x2))

      if (verbose) {
        message(sprintf(
          "Computing scores for %d restricted pairs.",
          pairs_no
        ))
      }

      d <- .C(
        "CalculateSignedScoreRestrictedPairsCore",
        as.integer(as.numeric(situation_v == labels[1])),
        as.integer(n),
        as.double(data1_tied),
        as.integer(m1),
        as.double(data2_tied),
        as.integer(m2),
        as.integer(pairs_ind[, 1] - 1L),
        as.integer(pairs_ind[, 2] - 1L),
        as.integer(pairs_no),
        as.double(numeric(pairs_no)),
        as.double(numeric(pairs_no)),
        as.double(numeric(pairs_no))
      )

      score <- d[[10]]
      names(score) <- sprintf(
        "%s,%s",
        restricted_pairs[, 1],
        restricted_pairs[, 2]
      )

      list(score = score / 2, labels = labels)
    }
  }
}


# Compute basic TSP score for a single feature pair
# @noRd
calculate_tsp_instance <- function(xi, xj, y, levels) {
  (sum(xi[y == levels[1]] < xj[y == levels[1]]) /
    sum(y == levels[1])) +
    (sum(xi[y == levels[2]] > xj[y == levels[2]]) /
      sum(y == levels[2]))
}
