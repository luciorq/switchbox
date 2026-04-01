test_that("swap_plot_genes returns a ggplot", {
  data(trainingData, package = "switchbox")
  p <- swap_plot_genes(
    matTraining,
    trainingGroup,
    genes = rownames(matTraining)[1:2]
  )
  expect_s3_class(p, "ggplot")
})

test_that("swap_plot_genes errors with missing genes", {
  data(trainingData, package = "switchbox")
  expect_error(
    swap_plot_genes(matTraining, trainingGroup, genes = "NOT_A_GENE"),
    "not found"
  )
})

test_that("swap_plot_genes errors with zero genes", {
  data(trainingData, package = "switchbox")
  expect_error(
    swap_plot_genes(matTraining, trainingGroup, genes = character(0)),
    "At least one"
  )
})

test_that("swap_plot_gene_pair_scatter returns a ggplot", {
  data(trainingData, package = "switchbox")
  p <- swap_plot_gene_pair_scatter(
    matTraining,
    trainingGroup,
    genes = rownames(matTraining)[1:2]
  )
  expect_s3_class(p, "ggplot")
})

test_that("swap_plot_gene_pair_scatter errors with < 2 genes", {
  data(trainingData, package = "switchbox")
  expect_error(
    swap_plot_gene_pair_scatter(
      matTraining,
      trainingGroup,
      genes = rownames(matTraining)[1]
    ),
    "two gene names"
  )
})

test_that("swap_plot_roc returns a ggplot", {
  data(trainingData, package = "switchbox")
  data(testingData, package = "switchbox")
  result <- swap_train_test_results(
    matTraining,
    trainingGroup,
    matTesting,
    testingGroup,
    k_range = 3:5,
    decision_values = TRUE
  )
  p <- swap_plot_roc(result)
  expect_s3_class(p, "ggplot")
})

test_that("swap_plot_roc errors without roc data", {
  expect_error(
    swap_plot_roc(list(train = 1, test = 2)),
    "train_roc"
  )
})

test_that("swap_plot_votes returns a ggplot", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:5)
  p <- swap_plot_votes(classifier, matTraining, trainingGroup)
  expect_s3_class(p, "ggplot")
})

test_that("swap_plot_gene_pair_boxplot returns a ggplot", {
  data(trainingData, package = "switchbox")
  p <- swap_plot_gene_pair_boxplot(
    rownames(matTraining)[1:2],
    matTraining,
    trainingGroup
  )
  expect_s3_class(p, "ggplot")
})

test_that("swap_plot_gene_pair_boxplot with points", {
  data(trainingData, package = "switchbox")
  for (col_mode in c("byGene", "byDirection", "byClass")) {
    p <- swap_plot_gene_pair_boxplot(
      rownames(matTraining)[1:2],
      matTraining,
      trainingGroup,
      points = TRUE,
      point_coloring = col_mode
    )
    expect_s3_class(p, "ggplot")
  }
})

test_that("swap_plot_gene_pair_classes_boxplot returns a ggplot", {
  data(trainingData, package = "switchbox")
  p <- swap_plot_gene_pair_classes_boxplot(
    rownames(matTraining)[1:2],
    matTraining,
    trainingGroup
  )
  expect_s3_class(p, "ggplot")
})

test_that("swap_plot_gene_pair_classes_boxplot with byClass ordering", {
  data(trainingData, package = "switchbox")
  p <- swap_plot_gene_pair_classes_boxplot(
    rownames(matTraining)[1:2],
    matTraining,
    trainingGroup,
    ordering = "byClass",
    points = TRUE
  )
  expect_s3_class(p, "ggplot")
})
