#' Plot Gene Expression Values Across Samples
#'
#' Creates a line plot showing the expression values of one or more genes
#' across samples, separated by class.
#'
#' @param x Numeric matrix (features x samples).
#' @param y Factor of class labels.
#' @param classes Character vector: `c(case, control)`. If missing, derived
#'   from `levels(y)`.
#' @param genes Character vector of gene names to plot.
#' @param colors Character vector of colors for each gene. If not provided,
#'   uses the default ggplot2 color palette.
#'
#' @returns A [ggplot2::ggplot] object.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' swap_plot_genes(
#'   matTraining, trainingGroup,
#'   genes = rownames(matTraining)[1:3]
#' )
swap_plot_genes <- function(x, y, classes, genes, colors = NULL) {
  if (length(genes) < 1L) {
    stop("At least one gene name required.", call. = FALSE)
  }
  if (!all(genes %in% rownames(x))) {
    stop("Gene names not found in matrix rownames.", call. = FALSE)
  }

  if (missing(classes)) {
    classes <- rev(levels(y))
  }

  # Order samples: control first, then case
  sample_order <- c(
    which(y == classes[2]),
    which(y == classes[1])
  )
  x_ordered <- x[, sample_order, drop = FALSE]
  y_ordered <- y[sample_order]

  # Build data frame
  n_samples <- ncol(x_ordered)
  df_list <- lapply(genes, function(g) {
    data.frame(
      sample_index = seq_len(n_samples),
      expression = x_ordered[g, ],
      gene = g,
      class = as.character(y_ordered),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, df_list)
  df$gene <- factor(df$gene, levels = genes)

  # Class boundary
  n_control <- sum(y == classes[2])

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$sample_index,
      y = .data$expression,
      color = .data$gene,
      shape = .data$class
    )
  ) +
    ggplot2::geom_line(ggplot2::aes(group = .data$gene), linewidth = 0.4) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_vline(
      xintercept = n_control + 0.5,
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    ggplot2::scale_shape_manual(values = c(17, 19)) +
    ggplot2::labs(
      x = "Sample",
      y = "Expression",
      color = "Gene",
      shape = "Class"
    ) +
    ggplot2::theme_minimal()

  if (!is.null(colors) && length(colors) == length(genes)) {
    p <- p + ggplot2::scale_color_manual(values = colors)
  }

  p
}


#' Scatter Plot of a Gene Pair
#'
#' Creates a scatter plot of two genes, colored by class, with a diagonal
#' reference line at y = x.
#'
#' @param x Numeric matrix (features x samples).
#' @param y Factor of class labels.
#' @param classes Character vector: `c(case, control)`.
#' @param genes Character vector of exactly 2 gene names.
#' @param colors Character vector of 2 colors for the two classes.
#'
#' @returns A [ggplot2::ggplot] object.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' swap_plot_gene_pair_scatter(
#'   matTraining, trainingGroup,
#'   genes = rownames(matTraining)[1:2]
#' )
swap_plot_gene_pair_scatter <- function(x, y, classes, genes, colors = NULL) {
  if (length(genes) < 2L) {
    stop("Requires exactly two gene names.", call. = FALSE)
  }
  if (!all(genes[1:2] %in% rownames(x))) {
    stop("Gene names not found in matrix rownames.", call. = FALSE)
  }
  genes <- genes[1:2]

  if (missing(classes)) {
    classes <- rev(levels(y))
  }
  if (is.null(colors)) {
    colors <- c("#2b8cbe", "#e34a33")
  }

  df <- data.frame(
    gene1 = x[genes[1], ],
    gene2 = x[genes[2], ],
    class = as.character(y),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$gene1,
      y = .data$gene2,
      color = .data$class,
      shape = .data$class
    )
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::scale_color_manual(values = stats::setNames(colors, classes)) +
    ggplot2::scale_shape_manual(values = c(17, 19)) +
    ggplot2::labs(
      x = genes[1],
      y = genes[2],
      color = "Class",
      shape = "Class"
    ) +
    ggplot2::theme_minimal()
}


#' Plot ROC Curves for Train and Test Sets
#'
#' Creates overlaid ROC curves from the result of
#' [swap_train_test_results()].
#'
#' @param result A list with `train_roc` and `test_roc` components (data
#'   frames with `fpr` and `tpr` columns).
#' @param colors Character vector of 2 colors for train and test curves.
#'
#' @returns A [ggplot2::ggplot] object.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' data(testingData)
#' result <- swap_train_test_results(
#'   matTraining, trainingGroup,
#'   matTesting, testingGroup,
#'   decision_values = TRUE
#' )
#' swap_plot_roc(result)
swap_plot_roc <- function(result, colors = NULL) {
  if (!all(c("train_roc", "test_roc") %in% names(result))) {
    stop(
      "'result' must have 'train_roc' and 'test_roc' components.",
      call. = FALSE
    )
  }

  if (is.null(colors)) {
    colors <- c("#2b8cbe", "#e34a33")
  }

  train_auc <- compute_auc(roc_data = result$train_roc)
  test_auc <- compute_auc(roc_data = result$test_roc)

  df_train <- result$train_roc
  df_train$set <- sprintf("Train (AUC = %.3f)", train_auc)
  df_test <- result$test_roc
  df_test$set <- sprintf("Test (AUC = %.3f)", test_auc)

  df <- rbind(df_train, df_test)

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$fpr,
      y = .data$tpr,
      color = .data$set
    )
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "grey50"
    ) +
    ggplot2::scale_color_manual(
      values = stats::setNames(
        colors,
        c(
          sprintf("Train (AUC = %.3f)", train_auc),
          sprintf("Test (AUC = %.3f)", test_auc)
        )
      )
    ) +
    ggplot2::labs(
      x = "False Positive Rate",
      y = "True Positive Rate",
      color = NULL,
      title = "ROC Curve"
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal()
}


#' Heatmap of TSP Votes
#'
#' Creates a heatmap of the individual TSP comparison votes per sample.
#'
#' @param classifier A KTSP classifier.
#' @param x Numeric matrix (features x samples).
#' @param y Optional factor of class labels (used for column annotation).
#' @param combine_fn Optional combine function passed to
#'   [swap_ktsp_statistics()].
#'
#' @returns A [ggplot2::ggplot] object.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:8)
#' swap_plot_votes(classifier, matTraining, trainingGroup)
swap_plot_votes <- function(classifier, x, y = NULL, combine_fn) {
  if (missing(combine_fn)) {
    votes <- swap_ktsp_statistics(x, classifier)$comparisons
  } else {
    votes <- swap_ktsp_statistics(x, classifier, combine_fn)$comparisons
  }

  votes <- votes * 1 # Ensure numeric

  # Build data frame for ggplot
  # Use unique sample IDs for x-axis (rownames or indices)
  sample_ids <- rownames(votes)
  if (is.null(sample_ids)) {
    sample_ids <- as.character(seq_len(nrow(votes)))
  }

  pair_names <- colnames(votes)
  if (is.null(pair_names)) {
    pair_names <- as.character(seq_len(ncol(votes)))
  }

  # Expand to long format
  df <- expand.grid(
    sample = seq_len(nrow(votes)),
    pair = seq_len(ncol(votes))
  )
  df$vote <- as.vector(votes)
  df$sample_id <- sample_ids[df$sample]
  df$pair_label <- pair_names[df$pair]

  if (!is.null(y)) {
    df$class <- as.character(y)[df$sample]
  }

  # Order samples by hierarchical clustering if enough samples
  if (nrow(votes) > 2L) {
    hc <- stats::hclust(stats::dist(votes))
    df$sample_id <- factor(df$sample_id, levels = sample_ids[hc$order])
  } else {
    df$sample_id <- factor(df$sample_id)
  }
  df$pair_label <- factor(df$pair_label, levels = pair_names)

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$sample_id,
      y = .data$pair_label,
      fill = .data$vote
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
    ggplot2::labs(x = "Sample", y = "TSP", fill = "Vote") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        size = 6
      )
    )
}


#' Boxplot of a Gene Pair
#'
#' Creates side-by-side boxplots for two genes, optionally overlaying
#' individual data points.
#'
#' @param genes Character vector of exactly 2 gene names.
#' @param x Numeric matrix (features x samples).
#' @param y Optional factor of class labels (for point coloring).
#' @param classes Character vector: `c(case, control)`.
#' @param points Logical; overlay individual data points? Default is `FALSE`.
#' @param point_coloring Character; one of `"byGene"`, `"byDirection"`, or
#'   `"byClass"`. Controls how points are colored when `points = TRUE`.
#' @param colors Character vector of 2 colors for the gene boxplots.
#'
#' @returns A [ggplot2::ggplot] object.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' swap_plot_gene_pair_boxplot(
#'   rownames(matTraining)[1:2],
#'   matTraining, trainingGroup,
#'   points = TRUE
#' )
swap_plot_gene_pair_boxplot <- function(
  genes,
  x,
  y = NULL,
  classes = NULL,
  points = FALSE,
  point_coloring = "byGene",
  colors = NULL
) {
  if (length(genes) != 2L) {
    stop("Requires exactly two gene names.", call. = FALSE)
  }
  if (!all(genes %in% rownames(x))) {
    stop("Gene names not found in matrix rownames.", call. = FALSE)
  }
  if (is.null(colors)) {
    colors <- c("#2b8cbe", "#e34a33")
  }

  df <- data.frame(
    expression = c(x[genes[1], ], x[genes[2], ]),
    gene = factor(rep(genes, each = ncol(x)), levels = genes),
    stringsAsFactors = FALSE
  )

  if (!is.null(y)) {
    df$class <- rep(as.character(y), 2)
  }

  # Direction: gene1 < gene2
  direction <- ifelse(
    x[genes[1], ] < x[genes[2], ],
    paste0(genes[1], " < ", genes[2]),
    paste0(genes[1], " >= ", genes[2])
  )
  df$direction <- rep(direction, 2)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$gene,
      y = .data$expression
    )
  ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(color = .data$gene),
      outlier.shape = if (points) NA else 19,
      linewidth = 0.8
    ) +
    ggplot2::scale_color_manual(values = stats::setNames(colors, genes)) +
    ggplot2::labs(x = "Gene", y = "Expression", color = "Gene") +
    ggplot2::theme_minimal()

  if (points) {
    if (point_coloring == "byGene") {
      p <- p +
        ggplot2::geom_jitter(
          ggplot2::aes(color = .data$gene),
          width = 0.15,
          size = 1,
          alpha = 0.6
        )
    } else if (point_coloring == "byDirection") {
      p <- p +
        ggplot2::geom_jitter(
          ggplot2::aes(fill = .data$direction),
          width = 0.15,
          size = 1.5,
          alpha = 0.6,
          shape = 21,
          color = "grey30"
        ) +
        ggplot2::scale_fill_brewer(palette = "Set2")
    } else if (point_coloring == "byClass" && !is.null(y)) {
      p <- p +
        ggplot2::geom_jitter(
          ggplot2::aes(fill = .data$class),
          width = 0.15,
          size = 1.5,
          alpha = 0.6,
          shape = 21,
          color = "grey30"
        ) +
        ggplot2::scale_fill_brewer(palette = "Set1")
    }
  }

  p
}


#' Boxplot of a Gene Pair Split by Class
#'
#' Creates boxplots showing the expression of two genes, split by phenotype
#' class.
#'
#' @param genes Character vector of exactly 2 gene names.
#' @param x Numeric matrix (features x samples).
#' @param y Factor of class labels.
#' @param classes Character vector: `c(case, control)`.
#' @param points Logical; overlay data points?
#' @param ordering Character; `"byGene"` groups boxplots by gene, `"byClass"`
#'   groups by class.
#' @param colors Character vector of 2 colors for the two classes.
#'
#' @returns A [ggplot2::ggplot] object.
#'
#' @export
#'
#' @examples
#' data(trainingData)
#' swap_plot_gene_pair_classes_boxplot(
#'   rownames(matTraining)[1:2],
#'   matTraining, trainingGroup,
#'   points = TRUE
#' )
swap_plot_gene_pair_classes_boxplot <- function(
  genes,
  x,
  y,
  classes = NULL,
  points = FALSE,
  ordering = "byGene",
  colors = NULL
) {
  if (length(genes) != 2L) {
    stop("Requires exactly two gene names.", call. = FALSE)
  }
  if (!all(genes %in% rownames(x))) {
    stop("Gene names not found in matrix rownames.", call. = FALSE)
  }
  if (is.null(classes)) {
    classes <- rev(levels(y))
  }
  if (is.null(colors)) {
    colors <- c("#2b8cbe", "#e34a33")
  }

  # Build long-format data frame
  df <- data.frame(
    expression = c(
      x[genes[1], y == classes[2]],
      x[genes[1], y == classes[1]],
      x[genes[2], y == classes[2]],
      x[genes[2], y == classes[1]]
    ),
    gene = factor(rep(genes, each = length(y)), levels = genes),
    class = factor(
      c(
        rep(classes[2], sum(y == classes[2])),
        rep(classes[1], sum(y == classes[1])),
        rep(classes[2], sum(y == classes[2])),
        rep(classes[1], sum(y == classes[1]))
      ),
      levels = classes
    ),
    stringsAsFactors = FALSE
  )

  if (ordering == "byGene") {
    df$group <- interaction(df$gene, df$class, sep = " - ")
    level_order <- as.vector(t(outer(genes, classes, paste, sep = " - ")))
    df$group <- factor(df$group, levels = level_order)
  } else {
    df$group <- interaction(df$class, df$gene, sep = " - ")
    level_order <- as.vector(t(outer(classes, genes, paste, sep = " - ")))
    df$group <- factor(df$group, levels = level_order)
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$group,
      y = .data$expression,
      color = .data$class
    )
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = if (points) NA else 19,
      linewidth = 0.8
    ) +
    ggplot2::scale_color_manual(values = stats::setNames(colors, classes)) +
    ggplot2::labs(x = NULL, y = "Expression", color = "Class") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )

  if (points) {
    p <- p + ggplot2::geom_jitter(width = 0.15, size = 1, alpha = 0.5)
  }

  p
}

# Required for ggplot2::aes() with .data pronoun
utils::globalVariables(".data")
