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
  expect_true(all(grepl(",", names(scores$score), fixed = TRUE)))
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

# -- Fixture tests for the C scoring kernels ----------------------------------
#
# The src/CalculateSignedScoreCore*.cc kernels are the original, peer-reviewed
# scoring routines from the switchBox paper (Afsari et al. 2015) and are
# treated as a trusted reference implementation, not re-derived here. These
# tests instead pin down (a) exact, hand-derivable values for the pure-R
# unsigned/basic score, which has no C-level ambiguity, and (b) invariants
# that must hold regardless of the exact sign convention used internally:
# the restricted-pairs C routine must agree with the unrestricted one for the
# same pair, and doubling/halving mistakes (like the one fixed in this
# release, see NEWS.md) must not silently reappear.

test_that("basic (unsigned) score is exactly hand-derivable for perfect separation", {
  # 2 genes, 4 samples, no ties: g1 > g2 in every "control" sample and
  # g1 < g2 in every "case" sample.
  x <- matrix(
    c(
      5, 1, # control s1: g1 > g2
      6, 2, # control s2: g1 > g2
      1, 5, # case s1:    g1 < g2
      2, 7 # case s2:    g1 < g2
    ),
    nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("s", 1:4))
  )
  y <- factor(c("control", "control", "case", "case"), levels = c("control", "case"))

  scores <- swap_calculate_basic_tsp_scores(y, x)

  # calculate_tsp_instance(xi, xj) = P(xi < xj | levels[1]) + P(xi > xj | levels[2])
  # (g2, g1): P(g2 < g1 | control) = 1, P(g2 > g1 | case) = 1 -> 2 (maximum)
  # (g1, g2): P(g1 < g2 | control) = 0, P(g1 > g2 | case) = 0 -> 0 (minimum)
  expect_equal(unname(scores$score[["g2,g1"]]), 2)
  expect_equal(unname(scores$score[["g1,g2"]]), 0)
  expect_false(scores$signed)
})

test_that("signed score (no ties) is stable and matches the restricted-pairs path", {
  x <- matrix(
    c(
      5, 1,
      6, 2,
      1, 5,
      2, 7
    ),
    nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("s", 1:4))
  )
  y <- factor(c("control", "control", "case", "case"), levels = c("control", "case"))

  unrestricted <- swap_calculate_signed_tsp_scores(y, x)
  restricted <- swap_calculate_signed_tsp_scores(
    y,
    x,
    restricted_pairs = matrix(c("g1", "g2"), nrow = 1)
  )

  # Pinned reference value from the current, trusted C kernel (deterministic
  # rank-based computation, no RNG involved -- reproducible exactly).
  expect_equal(unname(unrestricted$score[["g1,g2"]]), 1.000001, tolerance = 1e-6)

  # The two code paths (unrestricted matrix routine vs. restricted-pairs
  # routine) must agree exactly on the same pair.
  expect_equal(
    unname(restricted$score[["g1,g2"]]),
    unname(unrestricted$score[["g1,g2"]])
  )
})

test_that("signed score with ties agrees between restricted and unrestricted paths", {
  # Regression test for a bug where the unrestricted tie-handling branch of
  # calculate_signed_score() was missing the `/ 2` normalization applied by
  # every other branch (no-tie unrestricted, no-tie restricted, and the
  # restricted tie-handling branch), making handle_ties = TRUE scores exactly
  # double their correct magnitude whenever restricted_pairs was not used.
  x <- matrix(
    c(
      1, 2, # control s1: g1 < g2
      2, 2, # control s2: tie
      3, 1, # case s1:    g1 > g2
      1, 1 # case s2:    tie
    ),
    nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("s", 1:4))
  )
  y <- factor(c("control", "control", "case", "case"), levels = c("control", "case"))

  unrestricted <- swap_calculate_signed_tsp_scores(y, x, handle_ties = TRUE)
  restricted <- swap_calculate_signed_tsp_scores(
    y,
    x,
    restricted_pairs = matrix(c("g1", "g2"), nrow = 1),
    handle_ties = TRUE
  )

  expect_equal(unname(unrestricted$score[["g1,g2"]]), -0.2500005, tolerance = 1e-6)
  expect_equal(
    unname(restricted$score[[1]]),
    unname(unrestricted$score[["g1,g2"]])
  )
})
