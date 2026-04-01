#' Filter Features Using the Wilcoxon Rank-Sum Test
#'
#' Selects the top differentially expressed features based on the Wilcoxon
#' rank-sum test statistic. By default, returns equal numbers of features
#' positively and negatively associated with the phenotype.
#'
#' @param y A factor with exactly two levels indicating sample classes.
#' @param x A numeric matrix with features in rows and samples in columns.
#' @param n_features Integer; number of features to return. Default is `100`.
#' @param up_down Logical; if `TRUE` (default), return equal numbers of
#'   up-regulated and down-regulated features. If `FALSE`, return the top
#'   features by absolute test statistic regardless of direction.
#'
#' @returns A character vector of selected feature names.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' top_genes <- swap_filter_wilcoxon(trainingGroup, matTraining, n_features = 10)
#' top_genes
swap_filter_wilcoxon <- function(y, x, n_features = 100, up_down = TRUE) {
  check_input(y, x)

  # Compute Wilcoxon test statistic
  tied_data <- apply(x, 2, rank)
  tied_data_p <- t(apply(tied_data, 1, rank))

  n <- sum(y == levels(y)[1])
  m <- sum(y == levels(y)[2])

  sum_zeros <- apply(
    tied_data_p[, which(y == levels(y)[1]), drop = FALSE],
    1,
    sum
  )
  w_index <- (sum_zeros - n * (n + m + 1) / 2) /
    sqrt(n * m * (n + m + 1) / 12)

  if (up_down) {
    s <- order(w_index, decreasing = TRUE)
    len_s <- length(s)
    half <- round(n_features / 2)

    features_up <- s[seq_len(min(half, len_s))]
    features_down <- s[seq(max(len_s - half + 1L, 1L), len_s)]
    features_index <- unique(c(features_up, features_down))
  } else {
    s <- order(abs(w_index), decreasing = TRUE)
    features_index <- s[seq_len(min(n_features, length(s)))]
  }

  rownames(x)[features_index]
}
