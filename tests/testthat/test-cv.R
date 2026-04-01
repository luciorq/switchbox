test_that("swap_ktsp_cv runs k-fold cross-validation", {
  data(trainingData, package = "switchbox")
  set.seed(42)
  result <- swap_ktsp_cv(
    matTraining,
    trainingGroup,
    k = 3,
    k_range = 2:4
  )

  expect_type(result, "list")
  expect_true("cv" %in% names(result))
  expect_true("stats" %in% names(result))
  expect_true("predictions" %in% names(result))
  expect_true("folds" %in% names(result))
  expect_equal(length(result$cv), 3)
  expect_true("accuracy" %in% names(result$stats))
})

test_that("swap_ktsp_cv with pre-specified folds", {
  data(trainingData, package = "switchbox")
  folds <- list(1:26, 27:52, 53:78)
  result <- swap_ktsp_cv(
    matTraining,
    trainingGroup,
    folds = folds,
    k_range = 2:4
  )
  expect_equal(length(result$cv), 3)
})

test_that("swap_get_kfold_indices creates balanced folds", {
  data(trainingData, package = "switchbox")
  set.seed(42)
  folds <- swap_get_kfold_indices(trainingGroup, k = 4)

  expect_equal(length(folds), 4)
  # All indices are covered exactly once
  all_indices <- sort(unlist(folds))
  expect_equal(all_indices, seq_along(trainingGroup))
})

test_that("swap_get_kfold_indices maintains stratification", {
  y <- factor(rep(c("a", "b"), each = 20))
  set.seed(1)
  folds <- swap_get_kfold_indices(y, k = 4)

  # Each fold should have approximately equal class proportions
  for (fold in folds) {
    prop_a <- sum(y[fold] == "a") / length(fold)
    expect_true(abs(prop_a - 0.5) <= 0.15)
  }
})

test_that("swap_ktsp_loo runs leave-one-out", {
  data(trainingData, package = "switchbox")
  # Use a very small subset for speed
  small_x <- matTraining[1:20, 1:10]
  small_y <- trainingGroup[1:10]

  result <- swap_ktsp_loo(small_x, small_y, k_range = 2:3)

  expect_type(result, "list")
  expect_true("loo" %in% names(result))
  expect_true("stats" %in% names(result))
  expect_equal(length(result$loo), 10)
  expect_equal(length(result$predictions), 10)
})
