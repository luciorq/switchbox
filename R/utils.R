#' Convert Score Matrix to Named Vector
#'
#' Converts a matrix of pairwise scores to a named vector where names are
#' in `"gene1,gene2"` format.
#'
#' @param m A numeric matrix with rownames and colnames.
#'
#' @returns A named numeric vector.
#'
#' @export
#'
#' @examples
#' m <- matrix(1:9, 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
#' swap_score_matrix_to_vector(m)
swap_score_matrix_to_vector <- function(m) {
  score_matrix_to_vector(m)
}


#' Convert Named Score Vector to Matrix
#'
#' Converts a named vector of scores (with `"gene1,gene2"` names) back to
#' a matrix.
#'
#' @param v A named numeric vector.
#'
#' @returns A numeric matrix with appropriate rownames and colnames.
#'
#' @export
#'
#' @examples
#' v <- c("a,b" = 1, "a,c" = 2, "b,a" = 3, "b,c" = 4)
#' swap_score_vector_to_matrix(v)
swap_score_vector_to_matrix <- function(v) {
  score_vector_to_matrix(v)
}


#' Create Random Train/Test Split
#'
#' Randomly splits a dataset into training and testing sets while maintaining
#' class proportions.
#'
#' @param x Numeric matrix (features x samples).
#' @param y Factor of class labels.
#' @param classes Character vector: `c(case, control)`. If `NULL`, derived
#'   from `levels(y)`.
#' @param p Numeric; proportion of samples for training. Default is `0.5`.
#'
#' @returns A list with components:
#'   - `train_x`, `test_x`: data matrices
#'   - `train_y`, `test_y`: label factors
#'   - `train_ids`, `test_ids`: sample indices
#'   - `classes`: the class labels used
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' set.seed(42)
#' split <- swap_make_train_test_data(matTraining, trainingGroup)
#' dim(split$train_x)
#' dim(split$test_x)
swap_make_train_test_data <- function(x, y, classes = NULL, p = 0.5) {
  if (is.null(classes)) {
    classes <- rev(levels(y))
  }

  n0 <- sum(y == classes[2])
  n1 <- sum(y == classes[1])

  sel0 <- sample(seq_len(n0), floor(n0 * p))
  sel1 <- sample(seq_len(n1), floor(n1 * p))

  idx0 <- which(y == classes[2])
  idx1 <- which(y == classes[1])

  train_ids <- c(idx0[sel0], idx1[sel1])
  test_ids <- c(idx0[-sel0], idx1[-sel1])

  list(
    train_x = x[, train_ids, drop = FALSE],
    test_x = x[, test_ids, drop = FALSE],
    train_y = y[train_ids],
    test_y = y[test_ids],
    train_ids = train_ids,
    test_ids = test_ids,
    classes = classes
  )
}


# Internal: convert score matrix to named vector
# @noRd
score_matrix_to_vector <- function(m) {
  v <- as.vector(m)
  names(v) <- as.vector(t(vapply(
    rownames(m),
    function(rn) {
      vapply(
        colnames(m),
        function(cn) sprintf("%s,%s", rn, cn),
        character(1)
      )
    },
    character(ncol(m))
  )))
  v
}


# Internal: convert named vector to score matrix
# @noRd
score_vector_to_matrix <- function(v) {
  pairs <- strsplit(names(v), ",")
  rows <- sort(unique(vapply(pairs, `[`, character(1), 1)))
  cols <- sort(unique(vapply(pairs, `[`, character(1), 2)))

  m <- matrix(NA_real_, length(rows), length(cols), dimnames = list(rows, cols))
  for (i in seq_along(v)) {
    p <- pairs[[i]]
    m[p[1], p[2]] <- v[i]
  }
  m
}
