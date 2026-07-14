test_that("swap_train_ktsp trains a valid classifier", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:5)

  expect_type(classifier, "list")
  expect_true("name" %in% names(classifier))
  expect_true("TSPs" %in% names(classifier))
  expect_true("score" %in% names(classifier))
  expect_true("tie_vote" %in% names(classifier))
  expect_true("labels" %in% names(classifier))

  expect_true(is.matrix(classifier$TSPs))
  expect_equal(ncol(classifier$TSPs), 2)
  expect_true(nrow(classifier$TSPs) >= 1)
  expect_match(classifier$name, "TSPS$")
})

test_that("swap_train_1tsp trains a 1-TSP classifier", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_1tsp(matTraining, trainingGroup)

  expect_equal(classifier$name, "1TSP")
  expect_equal(nrow(classifier$TSPs), 1)
  expect_equal(ncol(classifier$TSPs), 2)
  expect_true(classifier$score > 0)
})

test_that("swap_train_ktsp respects k_range", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 2:4)
  n_pairs <- nrow(classifier$TSPs)
  expect_true(n_pairs >= 2 && n_pairs <= 4)
})

test_that("swap_train_ktsp with filter_fn = NULL", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_ktsp(
    matTraining,
    trainingGroup,
    filter_fn = NULL,
    k_range = 2:3
  )
  expect_true(nrow(classifier$TSPs) >= 1)
})

test_that("swap_train_ktsp with handle_ties = TRUE", {
  data(trainingData, package = "switchbox")
  classifier <- swap_train_ktsp(
    matTraining,
    trainingGroup,
    handle_ties = TRUE,
    k_range = 2:4
  )
  expect_true(is.factor(classifier$tie_vote))
  expect_true(nrow(classifier$TSPs) >= 1)
})

test_that("swap_train_ktsp with restricted_pairs", {
  data(trainingData, package = "switchbox")
  genes <- rownames(matTraining)
  rp <- cbind(
    rep(genes[1:10], each = 5),
    rep(genes[11:15], times = 10)
  )
  classifier <- swap_train_ktsp(
    matTraining,
    trainingGroup,
    restricted_pairs = rp,
    k_range = 2:3
  )
  # All TSP genes should come from restricted pairs
  tsp_genes <- as.vector(classifier$TSPs)
  expect_true(all(tsp_genes %in% as.vector(rp)))
})

test_that("swap_k_by_ttest returns valid indices", {
  data(trainingData, package = "switchbox")
  scores <- swap_calculate_scores(matTraining, trainingGroup)
  tbl <- swap_make_tsp_table(scores, maxk = 10)
  classes <- rev(levels(trainingGroup))

  sel <- swap_k_by_ttest(matTraining, trainingGroup, tbl, classes, 2:6, list())
  expect_type(sel, "integer")
  expect_true(length(sel) >= 2 && length(sel) <= 6)
  expect_true(max(sel) <= nrow(tbl))
})

test_that("swap_k_by_measurement returns valid indices", {
  data(trainingData, package = "switchbox")
  scores <- swap_calculate_scores(matTraining, trainingGroup)
  tbl <- swap_make_tsp_table(scores, maxk = 10)
  classes <- rev(levels(trainingGroup))

  sel <- swap_k_by_measurement(
    matTraining,
    trainingGroup,
    tbl,
    classes,
    2:6,
    k_opts = list(measurement = "accuracy")
  )
  expect_type(sel, "integer")
  expect_true(length(sel) >= 1)
})

test_that("swap_make_tsp_table returns a data frame", {
  data(trainingData, package = "switchbox")
  scores <- swap_calculate_scores(matTraining, trainingGroup)
  tbl <- swap_make_tsp_table(scores, maxk = 5)

  expect_s3_class(tbl, "data.frame")
  expect_named(tbl, c("gene1", "gene2", "score", "tie_vote"))
  expect_true(nrow(tbl) <= 5)
  # Scores should be decreasing
  expect_true(all(diff(tbl$score) <= 1e-10))
})

test_that("swap_train_1tsp preserves gene order for unsigned scores", {
  # swap_calculate_basic_tsp_scores() produces unsigned (non-negative,
  # asymmetric) scores. Unlike the signed score, a positive score must NOT
  # trigger a pair reversal, since the basic score has no "flip if positive"
  # convention: score(gene1, gene2) != score(gene2, gene1) in general.
  data(trainingData, package = "switchbox")
  scores <- swap_calculate_basic_tsp_scores(trainingGroup, matTraining)
  expect_false(scores$signed)

  best <- names(scores$score)[which.max(abs(scores$score))]
  best_pair <- strsplit(best, ",", fixed = TRUE)[[1]]

  classifier <- swap_train_1tsp(
    matTraining,
    trainingGroup,
    score_fn = swap_calculate_basic_tsp_scores
  )

  expect_equal(classifier$TSPs[1, 1], best_pair[1])
  expect_equal(classifier$TSPs[1, 2], best_pair[2])
  expect_equal(unname(classifier$score), max(scores$score))
})

test_that("swap_train_1tsp still reverses pairs for signed scores", {
  data(trainingData, package = "switchbox")
  scores <- swap_calculate_signed_tsp_scores(trainingGroup, matTraining)
  expect_true(scores$signed)

  j <- which.max(abs(scores$score))
  raw_pair <- strsplit(names(scores$score)[j], ",", fixed = TRUE)[[1]]
  raw_score <- scores$score[[j]]

  classifier <- swap_train_1tsp(matTraining, trainingGroup)

  if (raw_score > 0) {
    expect_equal(classifier$TSPs[1, 1], raw_pair[2])
    expect_equal(classifier$TSPs[1, 2], raw_pair[1])
  } else {
    expect_equal(classifier$TSPs[1, 1], raw_pair[1])
    expect_equal(classifier$TSPs[1, 2], raw_pair[2])
  }
})
