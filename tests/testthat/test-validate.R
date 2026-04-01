test_that("check_input rejects non-factor y", {
  x <- matrix(
    1:6,
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2", "s3"))
  )
  expect_error(
    check_input(y = c("a", "b", "a"), x = x),
    "factor with exactly two levels"
  )
})

test_that("check_input rejects y with != 2 levels", {
  x <- matrix(
    1:6,
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2", "s3"))
  )
  expect_error(
    check_input(y = factor(c("a", "b", "c")), x = x),
    "factor with exactly two levels"
  )
})

test_that("check_input rejects non-numeric x", {
  x <- matrix(letters[1:6], nrow = 2, dimnames = list(c("g1", "g2"), NULL))
  y <- factor(c("a", "b", "a"))
  expect_error(
    check_input(y = y, x = x),
    "numeric matrix"
  )
})

test_that("check_input rejects x without rownames", {
  x <- matrix(1:6, nrow = 2)
  y <- factor(c("a", "b", "a"))
  expect_error(
    check_input(y = y, x = x),
    "rownames"
  )
})

test_that("check_input rejects mismatched x/y dimensions", {
  x <- matrix(
    1:6,
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2", "s3"))
  )
  y <- factor(c("a", "b"))
  expect_error(
    check_input(y = y, x = x),
    "as many columns"
  )
})

test_that("check_input accepts valid inputs", {
  x <- matrix(
    1:6,
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2", "s3"))
  )
  y <- factor(c("a", "b", "a"))
  expect_invisible(check_input(y = y, x = x))
})

test_that("check_input rejects bad restricted_pairs", {
  x <- matrix(
    1:6,
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2", "s3"))
  )
  y <- factor(c("a", "b", "a"))

  # Not a matrix
  expect_error(
    check_input(y = y, x = x, restricted_pairs = c("g1", "g2")),
    "two-column character matrix"
  )

  # Wrong number of columns
  rp <- matrix(c("g1", "g2", "g1"), ncol = 3)
  expect_error(
    check_input(y = y, x = x, restricted_pairs = rp),
    "exactly 2 columns"
  )
})

test_that("format_input applies filter", {
  set.seed(1)
  x <- matrix(
    rnorm(200),
    nrow = 20,
    dimnames = list(paste0("g", 1:20), paste0("s", 1:10))
  )
  y <- factor(rep(c("a", "b"), each = 5))

  result <- format_input(
    y = y,
    x = x,
    filter_fn = swap_filter_wilcoxon,
    n_features = 10
  )
  expect_true(nrow(result$x) <= 10)
  expect_null(result$filtered_pairs)
})

test_that("format_input with filter_fn = NULL keeps all features", {
  x <- matrix(
    1:40,
    nrow = 10,
    dimnames = list(paste0("g", 1:10), paste0("s", 1:4))
  )
  y <- factor(c("a", "a", "b", "b"))

  result <- format_input(y = y, x = x, filter_fn = NULL)
  expect_equal(nrow(result$x), 10)
})
