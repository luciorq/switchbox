test_that("compute_roc returns valid ROC curve", {
  truth <- factor(c(rep("pos", 5), rep("neg", 5)))
  # Perfect separation

  dv <- c(rep(1, 5), rep(-1, 5))
  roc <- compute_roc(truth, dv, classes = c("pos", "neg"))

  expect_s3_class(roc, "data.frame")
  expect_named(roc, c("threshold", "tpr", "fpr"))
  # Starts at (0, 0)
  expect_equal(roc$fpr[1], 0)
  expect_equal(roc$tpr[1], 0)
  # Ends at (1, 1)
  expect_equal(roc$fpr[nrow(roc)], 1)
  expect_equal(roc$tpr[nrow(roc)], 1)
})

test_that("compute_roc returns NA for degenerate input", {
  truth <- factor(c("pos", "pos", "pos"))
  dv <- c(1, 2, 3)
  roc <- compute_roc(truth, dv, classes = c("pos", "neg"))
  expect_true(is.na(roc$tpr[1]))
})

test_that("compute_auc returns 1 for perfect separation", {
  truth <- factor(c(rep("pos", 5), rep("neg", 5)))
  dv <- c(rep(1, 5), rep(-1, 5))
  auc <- compute_auc(truth, dv, classes = c("pos", "neg"))
  expect_equal(auc, 1.0)
})

test_that("compute_auc returns ~0.5 for random predictions", {
  set.seed(42)
  truth <- factor(rep(c("pos", "neg"), each = 500))
  dv <- rnorm(1000)
  auc <- compute_auc(truth, dv, classes = c("pos", "neg"))
  # Should be close to 0.5 with some tolerance

  expect_true(abs(auc - 0.5) < 0.1)
})

test_that("compute_auc works with precomputed roc_data", {
  truth <- factor(c(rep("pos", 5), rep("neg", 5)))
  dv <- c(rep(1, 5), rep(-1, 5))
  roc_data <- compute_roc(truth, dv, classes = c("pos", "neg"))
  auc <- compute_auc(roc_data = roc_data)
  expect_equal(auc, 1.0)
})

test_that("compute_auc returns NA for degenerate input", {
  auc <- compute_auc(
    roc_data = data.frame(
      threshold = NA_real_,
      tpr = NA_real_,
      fpr = NA_real_
    )
  )
  expect_true(is.na(auc))
})
