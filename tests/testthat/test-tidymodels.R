# -- swap_from_tidy -----------------------------------------------------------

test_that("swap_from_tidy converts data frame to features-in-rows matrix", {
  df <- data.frame(
    gene1 = c(1.2, 3.4, 1.5),
    gene2 = c(5.6, 2.1, 5.9)
  )
  mat <- swap_from_tidy(df)

  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 2L) # features

  expect_equal(ncol(mat), 3L) # samples
  expect_equal(rownames(mat), c("gene1", "gene2"))
  expect_equal(mat["gene1", ], c(1.2, 3.4, 1.5))
})

test_that("swap_from_tidy extracts outcome column when specified", {
  df <- data.frame(
    class = factor(c("a", "b", "a", "b")),
    gene1 = c(1.2, 3.4, 1.5, 3.1),
    gene2 = c(5.6, 2.1, 5.9, 2.3)
  )
  result <- swap_from_tidy(df, outcome = "class")

  expect_type(result, "list")
  expect_true("x" %in% names(result))
  expect_true("y" %in% names(result))
  expect_true(is.matrix(result$x))
  expect_s3_class(result$y, "factor")
  expect_equal(nrow(result$x), 2L)
  expect_equal(ncol(result$x), 4L)
  expect_equal(levels(result$y), c("a", "b"))
  # outcome column should not be in the matrix

  expect_false("class" %in% rownames(result$x))
})

test_that("swap_from_tidy transposes a numeric matrix", {
  mat <- matrix(1:6, nrow = 3, ncol = 2, dimnames = list(NULL, c("g1", "g2")))
  result <- swap_from_tidy(mat)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2L)
  expect_equal(ncol(result), 3L)
  expect_equal(rownames(result), c("g1", "g2"))
})

test_that("swap_from_tidy errors on non-numeric matrix", {
  mat <- matrix(letters[1:6], nrow = 3)
  expect_error(swap_from_tidy(mat), "numeric")
})

test_that("swap_from_tidy errors when outcome column is missing", {
  df <- data.frame(gene1 = c(1, 2), gene2 = c(3, 4))
  expect_error(swap_from_tidy(df, outcome = "class"), "not found")
})

test_that("swap_from_tidy converts character outcome to factor", {
  df <- data.frame(
    group = c("ctrl", "treat", "ctrl"),
    g1 = c(1, 2, 3),
    g2 = c(4, 5, 6)
  )
  result <- swap_from_tidy(df, outcome = "group")
  expect_s3_class(result$y, "factor")
  expect_equal(as.character(result$y), c("ctrl", "treat", "ctrl"))
})

# -- swap_to_tidy -------------------------------------------------------------

test_that("swap_to_tidy converts features-in-rows matrix to data frame", {
  mat <- matrix(
    c(1.2, 5.6, 3.4, 2.1, 1.5, 5.9),
    nrow = 2,
    dimnames = list(c("gene1", "gene2"), c("s1", "s2", "s3"))
  )
  df <- swap_to_tidy(mat)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 3L)
  expect_equal(ncol(df), 2L)
  expect_true("gene1" %in% names(df))
  expect_true("gene2" %in% names(df))
})

test_that("swap_to_tidy prepends outcome column when y is provided", {
  mat <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2"))
  )
  y <- factor(c("a", "b"))
  df <- swap_to_tidy(mat, y, outcome = "group")

  expect_equal(names(df)[1], "group")
  expect_s3_class(df$group, "factor")
  expect_equal(as.character(df$group), c("a", "b"))
  expect_equal(ncol(df), 3L) # group + g1 + g2
})

test_that("swap_to_tidy uses default outcome name 'class'", {
  mat <- matrix(1:4, nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
  y <- factor(c("x", "y"))
  df <- swap_to_tidy(mat, y)

  expect_equal(names(df)[1], "class")
})

test_that("swap_to_tidy without y returns only features", {
  mat <- matrix(1:6, nrow = 3, dimnames = list(c("a", "b", "c"), c("s1", "s2")))
  df <- swap_to_tidy(mat)

  expect_equal(ncol(df), 3L)
  expect_equal(nrow(df), 2L)
  expect_equal(names(df), c("a", "b", "c"))
})

# -- round-trip ---------------------------------------------------------------

test_that("swap_from_tidy and swap_to_tidy are inverse operations", {
  data(trainingData, package = "switchbox")
  tidy_df <- swap_to_tidy(matTraining, trainingGroup, outcome = "group")
  result <- swap_from_tidy(tidy_df, outcome = "group")

  expect_equal(result$x, matTraining)
  expect_equal(result$y, trainingGroup)
})

test_that("round-trip without outcome preserves matrix", {
  mat <- matrix(
    c(1.1, 2.2, 3.3, 4.4, 5.5, 6.6),
    nrow = 2,
    dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
  )
  result <- swap_from_tidy(swap_to_tidy(mat))
  expect_equal(result, mat)
})

# -- swap_tidy_result ---------------------------------------------------------

test_that("swap_tidy_result returns correct structure", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:5)
  results <- swap_tidy_result(classifier, matTraining, trainingGroup)

  expect_s3_class(results, "data.frame")
  expect_true("truth" %in% names(results))
  expect_true(".pred_class" %in% names(results))
  expect_true(".decision_value" %in% names(results))
  expect_equal(nrow(results), ncol(matTraining))
})

test_that("swap_tidy_result truth and pred_class are factors", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:5)
  results <- swap_tidy_result(classifier, matTraining, trainingGroup)

  expect_s3_class(results$truth, "factor")
  expect_s3_class(results$.pred_class, "factor")
  expect_equal(levels(results$truth), levels(results$.pred_class))
})

test_that("swap_tidy_result decision values are normalized", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:5)
  results <- swap_tidy_result(classifier, matTraining, trainingGroup)

  # Decision values are normalized by n_pairs, so bounded by [-1, 1]
  expect_true(all(results$.decision_value >= -1))
  expect_true(all(results$.decision_value <= 1))
})

test_that("swap_tidy_result works on test data", {
  data(trainingData, package = "switchbox")
  data(testingData, package = "switchbox")
  classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:5)
  results <- swap_tidy_result(classifier, matTesting, testingGroup)

  expect_equal(nrow(results), ncol(matTesting))
  expect_equal(as.character(results$truth), as.character(testingGroup))
})

test_that("swap_tidy_result works with 1TSP classifier", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_1tsp(matTraining, trainingGroup)
  results <- swap_tidy_result(classifier, matTraining, trainingGroup)

  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), ncol(matTraining))
  expect_true(all(
    c("truth", ".pred_class", ".decision_value") %in% names(results)
  ))
})
