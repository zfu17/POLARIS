## Each test below that references a bug is pinning behaviour that was wrong in
## the pre-package scripts, so a regression would be caught rather than
## rediscovered.

sim_fit <- function(...) {
  set.seed(11)
  sim <- polarisSimulate(n = 150, g = 50, p = 70, r = 4)
  list(sim = sim, fit = Polaris(sim$X, sim$Y, x.sds = 0, verbose = FALSE, ...))
}

test_that("Polaris returns a well-formed skymap", {
  z <- sim_fit()
  fit <- z$fit
  expect_s3_class(fit, "skymap")
  expect_equal(nrow(fit$skymap.cell), 150)
  expect_equal(ncol(fit$skymap.cell), ncol(fit$skymap.U) + ncol(fit$skymap.V))
  expect_identical(rownames(fit$skymap.cell), rownames(z$sim$X))
  expect_identical(rownames(fit$skymap.P), colnames(z$sim$X))
  expect_identical(rownames(fit$skymap.Q), colnames(z$sim$Y))
  ## U, V, P, Q must be matrices, not data.frames: as data.frames,
  ## SkymapSimScore(metric = "cor") silently produced an all-NA result of the
  ## wrong length.
  for (s in c("skymap.U", "skymap.V", "skymap.P", "skymap.Q"))
    expect_true(is.matrix(fit[[s]]), info = s)
})

test_that("the fit records the ranks it actually used", {
  fit <- sim_fit()$fit
  ip <- fit$input.param
  ## The pre-package version captured input.param before any default resolved,
  ## so the saved object recorded only the thresholds and the ranks were NULL.
  for (nm in c("r1", "r2", "r3", "r4", "r5"))
    expect_true(is.numeric(ip[[nm]]) && length(ip[[nm]]) == 1L, info = nm)
  expect_equal(ip$r3, ncol(fit$skymap.U))
  expect_equal(ip$r5, ncol(fit$skymap.P))
})

test_that("Polaris is deterministic and leaves the caller's RNG alone", {
  set.seed(11); sim <- polarisSimulate(n = 120, g = 40, p = 50, r = 3)
  a <- Polaris(sim$X, sim$Y, x.sds = 0, verbose = FALSE)
  b <- Polaris(sim$X, sim$Y, x.sds = 0, verbose = FALSE)
  expect_equal(a$stdev.cell, b$stdev.cell, tolerance = 0)

  ## irlba draws a random start, so without the internal seed handling a fit
  ## would advance the user's stream and break their reproducibility.
  set.seed(99); before <- stats::runif(1)
  set.seed(99); invisible(Polaris(sim$X, sim$Y, x.sds = 0, verbose = FALSE))
  after <- stats::runif(1)
  expect_equal(before, after)
})

test_that("mismatched or unnamed cells are caught", {
  set.seed(11); sim <- polarisSimulate(n = 80, g = 30, p = 40, r = 3)
  expect_error(Polaris(sim$X, sim$Y[1:79, ], verbose = FALSE),
               "same number of cells")
  shuffled <- sim$Y[sample(nrow(sim$Y)), ]
  expect_error(Polaris(sim$X, shuffled, x.sds = 0, verbose = FALSE),
               "differently ordered cell names")
  bare.x <- sim$X; bare.y <- sim$Y
  rownames(bare.x) <- rownames(bare.y) <- NULL
  expect_warning(Polaris(bare.x, bare.y, x.sds = 0, verbose = FALSE),
                 "no cell names")
})

test_that("an over-strict variance floor is reported, not left to fail later", {
  set.seed(11); sim <- polarisSimulate(n = 80, g = 30, p = 40, r = 3)
  ## x.sds = 0.95 on non-z-scored input silently emptied X in the pre-package
  ## code, and the failure surfaced much further downstream.
  expect_error(Polaris(sim$X, sim$Y, x.sds = 100, verbose = FALSE),
               "standard deviation above x.sds")
})

test_that("a flat spectrum gives a named error, not '1:-Inf'", {
  ## max(which(...)) returned -Inf when no ratio cleared the threshold, and the
  ## fit then died on `1:r` with "result would be too long a vector".
  set.seed(11); sim <- polarisSimulate(n = 100, g = 40, p = 50, r = 3)
  expect_error(Polaris(sim$X, sim$Y, T1 = 1e6, x.sds = 0, verbose = FALSE),
               "No spectral gap")
})

test_that("modalities with 50 or fewer features still fit", {
  ## min(50, ncol(X)) equalled min(dim) here, which made irlba hard-error.
  set.seed(11); sim <- polarisSimulate(n = 200, g = 40, p = 45, r = 4)
  expect_no_error(Polaris(sim$X, sim$Y, x.sds = 0, verbose = FALSE))
})

test_that("gene.chr.ref is optional and its absence is explained", {
  set.seed(11); sim <- polarisSimulate(n = 120, g = 40, p = 50, r = 3)
  fit <- Polaris(sim$X, sim$Y, x.sds = 0, verbose = FALSE)
  expect_null(fit$skymap.feature.chr)
  ## The default gene.chr.ref = NULL used to fail inside mclapply and return a
  ## list of try-error objects that looked like a successful fit.
  expect_error(SkymapLinkageTable(fit, gene.anno = sim$anno, verbose = FALSE),
               "no per-chromosome feature embedding")
})

test_that("the per-chromosome embedding is built when an annotation is given", {
  set.seed(11); sim <- polarisSimulate(n = 200, g = 60, p = 90, r = 5)
  fit <- Polaris(sim$X, sim$Y, gene.chr.ref = sim$anno, x.sds = 0, verbose = FALSE)
  expect_true(is.list(fit$skymap.feature.chr))
  expect_true("chr1" %in% names(fit$skymap.feature.chr))
  cc <- fit$skymap.feature.chr[["chr1"]]
  expect_named(cc, c("P", "Q", "stdev"))
  expect_equal(nrow(cc$P) , 60)
})
