prov_fit <- function(...) {
  set.seed(7)
  sim <- polarisSimulate(n = 180, g = 55, p = 70, r = 4)
  list(sim = sim, fit = Polaris(sim$X, sim$Y, x.sds = 0, verbose = FALSE, ...))
}

test_that("polarisVersions reports POLARIS, R and the Imports", {
  v <- polarisVersions()
  expect_s3_class(v, "data.frame")
  expect_named(v, c("package", "version"))
  expect_identical(v$package[1:2], c("POLARIS", "R"))
  expect_false(anyNA(v$version))
  ## the packages the core actually needs must appear
  for (p in c("RSpectra", "irlba", "Matrix", "MatrixGenerics"))
    expect_true(p %in% v$package, info = p)
  ## base packages carry R's version and would be noise
  expect_false(any(c("methods", "stats", "utils", "parallel") %in% v$package))
  expect_gt(nrow(polarisVersions("all")), nrow(v))
})

test_that("a fit records the software that produced it", {
  ## Supplementary Table 2 needs versions, and a whole-archive search of the
  ## pre-package code returned none. Every fit now carries its own.
  fit <- prov_fit()$fit
  vs <- fit$input.param$versions
  expect_true(is.list(vs))
  expect_true("POLARIS" %in% names(vs))
  expect_true("R" %in% names(vs))
  expect_identical(vs$POLARIS, as.character(utils::packageVersion("POLARIS")))
})

test_that("polarisFitSummary tabulates dimensions and selected ranks", {
  z <- prov_fit()
  s <- polarisFitSummary(demo = z$fit)
  expect_equal(nrow(s), 1)
  expect_equal(s$dataset, "demo")
  expect_equal(s$n.cells, 180)
  expect_equal(s$r3, ncol(z$fit$skymap.U))
  expect_equal(s$r4, ncol(z$fit$skymap.V))
  expect_equal(s$r5, ncol(z$fit$skymap.P))
  expect_equal(s$T1, 1.01)
  expect_equal(s$x.sds, 0)
  expect_false(is.na(s$polaris.version))
})

test_that("polarisFitSummary takes several fits, named or in a list", {
  a <- prov_fit()$fit
  b <- prov_fit(T3 = 1.2)$fit
  s <- polarisFitSummary(one = a, two = b)
  expect_equal(nrow(s), 2)
  expect_equal(s$dataset, c("one", "two"))
  expect_equal(polarisFitSummary(list(one = a, two = b)), s)
  ## unnamed fits get numbered rather than erroring
  expect_equal(nrow(polarisFitSummary(a, b)), 2)
  expect_error(polarisFitSummary(), "at least one skymap")
})

test_that("polarisFitSummary recovers ranks from a pre-package fit", {
  ## The old code captured input.param before defaults resolved, so the ranks
  ## were NULL. They have to fall back to the embedding widths.
  fit <- prov_fit()$fit
  fit$input.param <- list(T1 = 1.01, T2 = 1.01, T3 = 1.1)   # as the old code saved it
  s <- polarisFitSummary(old = fit)
  expect_equal(s$r3, ncol(fit$skymap.U))
  expect_equal(s$r5, ncol(fit$skymap.P))
  expect_true(is.na(s$r1))
})

test_that("polarisSpectrum exposes the ratios the denoiser thresholds on", {
  z <- prov_fit()
  sp <- polarisSpectrum(z$fit, "X")
  expect_named(sp, c("modality", "component", "singular_value", "ratio",
                     "threshold", "retained"))
  expect_true(all(sp$modality == "X"))
  expect_equal(sp$component, seq_len(nrow(sp)))
  ## singular values are non-increasing, so every ratio is at least 1
  expect_true(all(sp$ratio[-nrow(sp)] >= 1 - 1e-12))
  expect_true(is.na(sp$ratio[nrow(sp)]))
  ## the retained block is a prefix, and its width is r1
  expect_equal(sum(sp$retained), z$fit$input.param$r1)
  expect_true(all(diff(as.integer(sp$retained)) <= 0))
  ## and the rule itself: the last retained component clears the threshold
  k <- sum(sp$retained)
  expect_gt(sp$ratio[k], sp$threshold[1])
})

test_that("polarisSpectrum covers all three modalities and the full feature spectrum", {
  fit <- prov_fit()$fit
  all3 <- polarisSpectrum(fit)
  expect_setequal(unique(all3$modality), c("X", "Y", "feature"))
  ## the feature spectrum must extend BEYOND r5, else the denoiser figure
  ## cannot be redrawn from a saved fit
  fe <- all3[all3$modality == "feature", ]
  expect_gt(nrow(fe), fit$input.param$r5)
  expect_equal(nrow(fe), length(fit$stdev.feature.full))
  expect_equal(sum(fe$retained), fit$input.param$r5)
  ## the retained values match the truncated slot the linkage score uses
  expect_equal(fe$singular_value[fe$retained], as.numeric(fit$stdev.feature))
})

test_that("polarisSpectrum degrades gracefully on a fit without the full spectrum", {
  fit <- prov_fit()$fit
  fit$stdev.feature.full <- NULL          # a pre-package fit
  fe <- polarisSpectrum(fit, "feature")
  expect_equal(nrow(fe), length(fit$stdev.feature))
  expect_true(all(fe$retained))
})
