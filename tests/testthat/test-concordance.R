fit_for_score <- function(n = 200) {
  set.seed(21)
  sim <- polarisSimulate(n = n, g = 60, p = 80, r = 5)
  list(sim = sim, fit = Polaris(sim$X, sim$Y, x.sds = 0, verbose = FALSE))
}

test_that("the default metric is cosine alone", {
  ## The POLARIS Methods state the score is computed "with its default cosine
  ## metric", but match.arg(several.ok = TRUE) against a full-choices default
  ## meant a bare call computed all four, including the slow per-cell cor loop.
  z <- fit_for_score()
  out <- SkymapSimScore(z$fit)
  expect_false(is.null(out$simscore.cosine))
  expect_null(out$simscore.euclidean)
  expect_null(out$simscore.inner.prod)
  expect_null(out$simscore.cor)
})

test_that("concordance is the rescaled cosine and lies in [0,1]", {
  z <- fit_for_score()
  out <- SkymapSimScore(z$fit)
  expect_equal(unname(out$concordance), (out$simscore.cosine + 1) / 2)
  expect_true(all(out$concordance >= 0 & out$concordance <= 1))
  expect_length(out$concordance, nrow(z$fit$skymap.U))
  expect_identical(names(out$concordance), rownames(z$fit$skymap.U))
})

test_that("metric = 'cor' returns one finite value per cell", {
  ## As data.frames, U[i, ] stayed a one-row data.frame, so cor() returned a
  ## k x k all-NA matrix and the stored result had length k^2 * n_cells.
  z <- fit_for_score(n = 100)
  out <- SkymapSimScore(z$fit, metric = "cor")
  expect_length(out$simscore.cor, 100)
  expect_false(anyNA(out$simscore.cor))
  expect_true(all(out$simscore.cor >= -1 & out$simscore.cor <= 1))
})

test_that("all four metrics can be requested together", {
  z <- fit_for_score(n = 80)
  out <- SkymapSimScore(z$fit, metric = c("cosine", "euclidean", "inner.prod", "cor"))
  for (m in c("cosine", "euclidean", "inner.prod", "cor"))
    expect_length(out[[paste0("simscore.", m)]], 80)
})

test_that("n.pc = 1 does not collapse to a vector", {
  ## U[, 1:1] on a one-column data.frame dropped to a numeric and rowSums()
  ## then errored with "must be an array of at least two dimensions".
  z <- fit_for_score(n = 80)
  expect_no_error(SkymapSimScore(z$fit, n.pc = 1))
  out <- SkymapSimScore(z$fit, n.pc = 1)
  expect_length(out$simscore.cosine, 80)
  ## cor needs 2 components and should say so rather than return nonsense
  expect_warning(SkymapSimScore(z$fit, metric = "cor", n.pc = 1),
                 "at least 2 components")
})

test_that("an oversized n.pc warns and is clamped", {
  z <- fit_for_score(n = 80)
  expect_warning(SkymapSimScore(z$fit, n.pc = 999), "exceeds")
})

test_that("discordant cells score lower, which is the point of the score", {
  set.seed(31)
  sim <- polarisSimulate(n = 400, g = 80, p = 100, r = 5, discordant = 0.3)
  fit <- SkymapSimScore(Polaris(sim$X, sim$Y, x.sds = 0, verbose = FALSE))
  lo <- mean(fit$concordance[sim$shared < 0.5])
  hi <- mean(fit$concordance[sim$shared >= 0.5])
  expect_lt(lo, hi)
})
