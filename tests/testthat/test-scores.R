test_that("swap_calculate_scores returns expected structure", {
  data(trainingData, package = "switchbox")
  scores <- swap_calculate_scores(matTraining, trainingGroup)

  expect_type(scores, "list")
  expect_true("score" %in% names(scores))
  expect_true("labels" %in% names(scores))
  expect_true("tie_vote" %in% names(scores))
  expect_true("signed" %in% names(scores))
  expect_true(scores$signed)
  expect_type(scores$score, "double")
  expect_true(length(scores$score) > 0)
})

test_that("swap_calculate_signed_tsp_scores works", {
  data(trainingData, package = "switchbox")
  # Use a small subset for speed
  x_small <- matTraining[1:10, ]
  scores <- swap_calculate_signed_tsp_scores(trainingGroup, x_small)

  expect_type(scores, "list")
  expect_true(scores$signed)
  expect_true(length(scores$score) > 0)
  expect_true(all(grepl(",", names(scores$score))))
})

test_that("swap_calculate_basic_tsp_scores works", {
  data(trainingData, package = "switchbox")
  x_small <- matTraining[1:10, ]
  scores <- swap_calculate_basic_tsp_scores(trainingGroup, x_small)

  expect_type(scores, "list")
  expect_false(scores$signed)
  expect_true(length(scores$score) > 0)
})

test_that("swap_calculate_scores with restricted_pairs", {
  data(trainingData, package = "switchbox")
  genes <- rownames(matTraining)
  rp <- cbind(genes[1:5], genes[6:10])

  scores <- swap_calculate_scores(
    matTraining,
    trainingGroup,
    filter_fn = NULL,
    restricted_pairs = rp
  )
  expect_true(length(scores$score) == 5)
})

test_that("swap_calculate_scores with handle_ties = TRUE", {
  data(trainingData, package = "switchbox")
  scores <- swap_calculate_scores(
    matTraining,
    trainingGroup,
    handle_ties = TRUE
  )
  expect_type(scores, "list")
  expect_true(length(scores$score) > 0)
})

test_that("score names have comma-separated format", {
  data(trainingData, package = "switchbox")
  x_small <- matTraining[1:10, ]
  scores <- swap_calculate_signed_tsp_scores(trainingGroup, x_small)
  expect_true(all(grepl("^[^,]+,[^,]+$", names(scores$score))))
})
