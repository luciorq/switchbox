#' Validate Inputs for switchbox Functions
#'
#' Checks that `y` is a two-level factor, `x` is a numeric matrix with
#' rownames and matching dimensions, and `restricted_pairs` (if provided)
#' is a two-column character matrix with names found in `x`.
#'
#' @param y A factor with exactly two levels representing class labels.
#' @param x A numeric matrix with features in rows and samples in columns.
#'   Must have rownames.
#' @param restricted_pairs An optional two-column character matrix specifying
#'   feature pairs. Both columns must contain strings matching `rownames(x)`.
#'
#' @return Invisible `NULL`. Called for its side effect of raising errors
#'   on invalid input.
#'
#' @noRd
check_input <- function(y, x, restricted_pairs) {
  # Check y

  if (!is.factor(y) || length(levels(y)) != 2L) {
    stop(
      "'y' must be a factor with exactly two levels.",
      call. = FALSE
    )
  }

  # Check x
  if (!missing(x) && !is.null(x)) {
    if (!is.matrix(x) || !is.numeric(x)) {
      stop("'x' must be a numeric matrix.", call. = FALSE)
    }
    if (is.null(rownames(x))) {
      stop("'x' must have rownames (feature names).", call. = FALSE)
    }
    if (length(y) != ncol(x)) {
      stop(
        "'x' must have as many columns as the length of 'y'.",
        call. = FALSE
      )
    }
  }

  # Check restricted_pairs
  if (!missing(restricted_pairs) && !is.null(restricted_pairs)) {
    if (!is.matrix(restricted_pairs) || !is.character(restricted_pairs)) {
      stop(
        "'restricted_pairs' must be a two-column character matrix.",
        call. = FALSE
      )
    }
    if (ncol(restricted_pairs) != 2L) {
      stop(
        "'restricted_pairs' must have exactly 2 columns.",
        call. = FALSE
      )
    }
    if (
      all(!restricted_pairs[, 1] %in% rownames(x)) ||
        all(!restricted_pairs[, 2] %in% rownames(x))
    ) {
      stop(
        "None of the 'restricted_pairs' match 'rownames(x)'.",
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}


#' Format and Filter Input Data
#'
#' Applies a feature filtering function and intersects with restricted pairs.
#'
#' @param y Factor with two levels.
#' @param x Numeric matrix (features x samples).
#' @param filter_fn Function for feature filtering, or `NULL` to skip.
#' @param restricted_pairs Optional two-column character matrix.
#' @param verbose Logical; print progress messages?
#' @param ... Additional arguments passed to `filter_fn`.
#'
#' @return A list with components:
#'   - `x`: the filtered numeric matrix
#'   - `filtered_pairs`: the filtered restricted pairs (or `NULL`)
#'
#' @noRd
format_input <- function(
  y,
  x,
  filter_fn,
  restricted_pairs = NULL,
  verbose = FALSE,
  ...
) {
  # Apply filter
  if (is.null(filter_fn)) {
    filtered <- rownames(x)
    if (verbose) message("No feature filtering applied.")
  } else {
    if (verbose) {
      message("Applying feature filter to 'x'...")
    }
    filtered <- filter_fn(y, x, ...)
  }

  x <- x[filtered, , drop = FALSE]

  # Unrestricted case
  if (is.null(restricted_pairs)) {
    if (nrow(x) < 4L) {
      stop("Not enough features left after filtering!", call. = FALSE)
    }
    return(list(x = x, filtered_pairs = NULL))
  }

  # Restricted case: intersect filtered features with pairs
  available <- intersect(as.vector(restricted_pairs), rownames(x))
  if (length(available) < 4L) {
    stop(
      "Not enough features left after filtering and pair restriction!",
      call. = FALSE
    )
  }

  if (verbose) {
    message("Restricting analysis to provided candidate TSPs.")
  }
  x <- x[available, , drop = FALSE]
  features <- rownames(x)

  keep <- which(
    (restricted_pairs[, 1] %in% features) &
      (restricted_pairs[, 2] %in% features)
  )
  if (length(keep) < 2L) {
    stop(
      "Not enough pairs left after filtering and restriction!",
      call. = FALSE
    )
  }

  list(x = x, filtered_pairs = restricted_pairs[keep, , drop = FALSE])
}
