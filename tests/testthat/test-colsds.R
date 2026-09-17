## .col_sds() decides which features clear the x.sds / y.sds variance floors.
## It replaced MatrixGenerics::colSds / sparseMatrixStats::colSds, so these tests
## lock it to stats::sd, to the dense/sparse agreement the filter assumes, and to
## the blocking being invisible.

test_that(".col_sds matches stats::sd on a dense matrix", {
  set.seed(11)
  m <- matrix(rnorm(40 * 7, mean = 3, sd = 2), nrow = 40)
  expect_equal(POLARIS:::.col_sds(m), apply(m, 2, stats::sd))
})

test_that(".col_sds matches stats::sd on a sparse matrix", {
  skip_if_not_installed("Matrix")
  set.seed(12)
  d <- matrix(0, nrow = 50, ncol = 6)
  d[sample(length(d), 90)] <- rnorm(90, mean = 1, sd = 3)
  s <- methods::as(d, "dgCMatrix")
  expect_true(methods::is(s, "sparseMatrix"))
  expect_equal(POLARIS:::.col_sds(s), apply(d, 2, stats::sd))
})

test_that("the dense and sparse paths agree on the same data", {
  set.seed(13)
  d <- matrix(0, nrow = 60, ncol = 8)
  d[sample(length(d), 150)] <- rnorm(150)
  expect_equal(POLARIS:::.col_sds(d),
               POLARIS:::.col_sds(methods::as(d, "dgCMatrix")))
})

test_that("column blocking does not change the dense result", {
  set.seed(14)
  m <- matrix(rnorm(30 * 20), nrow = 30)
  full <- POLARIS:::.col_sds(m, block = 1000L)
  expect_identical(POLARIS:::.col_sds(m, block = 1L), full)
  expect_identical(POLARIS:::.col_sds(m, block = 3L), full)
  expect_identical(POLARIS:::.col_sds(m, block = 7L), full)
})

test_that("a constant column gives exactly zero, not a negative root", {
  m <- cbind(rep(2.5, 25), rnorm(25))
  expect_identical(POLARIS:::.col_sds(m)[1], 0)
  s <- methods::as(cbind(rep(0, 25), c(rep(0, 24), 1)), "dgCMatrix")
  expect_identical(POLARIS:::.col_sds(s)[1], 0)
  expect_false(any(is.na(POLARIS:::.col_sds(s))))
})

test_that(".col_sds needs at least two rows", {
  expect_error(POLARIS:::.col_sds(matrix(1:3, nrow = 1)), "two cells")
})

test_that("feature names are carried through, dense and sparse alike", {
  set.seed(16)
  d <- matrix(rnorm(30 * 4), nrow = 30,
              dimnames = list(NULL, c("g1", "g2", "g3", "g4")))
  expect_named(POLARIS:::.col_sds(d), c("g1", "g2", "g3", "g4"))
  expect_named(POLARIS:::.col_sds(d, block = 2L), c("g1", "g2", "g3", "g4"))
  expect_named(POLARIS:::.col_sds(methods::as(d, "dgCMatrix")),
               c("g1", "g2", "g3", "g4"))
  ## and stay absent when the input has none
  expect_null(names(POLARIS:::.col_sds(matrix(rnorm(20), nrow = 5))))
})

test_that("the variance floor keeps the features it should", {
  set.seed(15)
  n <- 200
  m <- cbind(low = rnorm(n, sd = 0.10), high = rnorm(n, sd = 4), flat = rep(1, n))
  keep <- POLARIS:::.col_sds(m) > 0.95
  expect_equal(unname(keep), c(FALSE, TRUE, FALSE))
})
