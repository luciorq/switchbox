test_that("swap_score_matrix_to_vector converts correctly", {
  m <- matrix(1:9, 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  v <- swap_score_matrix_to_vector(m)

  expect_type(v, "integer")
  expect_true(all(grepl(",", names(v), fixed = TRUE)))
  expect_length(v, 9)
  expect_equal(v[["a,a"]], 1L)
})

test_that("swap_score_vector_to_matrix converts correctly", {
  v <- c("a,b" = 1, "a,c" = 2, "b,a" = 3, "b,c" = 4)
  m <- swap_score_vector_to_matrix(v)

  expect_true(is.matrix(m))
  expect_equal(m["a", "b"], 1)
  expect_equal(m["b", "c"], 4)
})

test_that("matrix-to-vector-to-matrix roundtrip", {
  m <- matrix(
    c(0.1, 0.2, 0.3, 0.4),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("g1", "g2"))
  )
  v <- swap_score_matrix_to_vector(m)
  m2 <- swap_score_vector_to_matrix(v)

  expect_equal(dim(m2), dim(m))
  expect_equal(m2["g1", "g2"], m["g1", "g2"])
  expect_equal(m2["g2", "g1"], m["g2", "g1"])
})

test_that("swap_make_train_test_data creates valid split", {
  data(trainingData, package = "switchbox")
  set.seed(42)
  split <- swap_make_train_test_data(matTraining, trainingGroup, p = 0.6)

  expect_type(split, "list")
  expect_true("train_x" %in% names(split))
  expect_true("test_x" %in% names(split))
  expect_true("train_y" %in% names(split))
  expect_true("test_y" %in% names(split))

  # No overlap between train and test
  expect_length(intersect(split$train_ids, split$test_ids), 0)
  # All samples accounted for
  expect_equal(
    sort(c(split$train_ids, split$test_ids)),
    seq_len(ncol(matTraining))
  )
})

test_that("swap_make_train_test_data respects proportion", {
  data(trainingData, package = "switchbox")
  set.seed(1)
  split <- swap_make_train_test_data(matTraining, trainingGroup, p = 0.7)

  n_train <- ncol(split$train_x)
  n_total <- ncol(matTraining)
  ratio <- n_train / n_total
  expect_true(ratio > 0.5 && ratio < 0.85)
})
