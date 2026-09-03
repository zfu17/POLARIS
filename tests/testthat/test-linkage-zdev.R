linkage_fit <- function() {
  set.seed(41)
  sim <- polarisSimulate(n = 250, g = 60, p = 90, r = 5)
  fit <- Polaris(sim$X, sim$Y, gene.chr.ref = sim$anno, x.sds = 0, verbose = FALSE)
  list(sim = sim, fit = fit)
}

test_that("the linkage table has the documented schema", {
  z <- linkage_fit()
  tab <- SkymapLinkageTable(z$fit, gene.anno = z$sim$anno, verbose = FALSE)
  expect_true(all(c("gene", "peak", "chr", "dist", "score_raw",
                    "SE", "score_adj", "p_lfsr", "q_lfsr") %in% names(tab)))
  expect_gt(nrow(tab), 0)
  expect_true(all(tab$dist >= 0))
  expect_true(all(tab$SE >= 0, na.rm = TRUE))
  ## sorted by decreasing score_adj
  expect_equal(tab$score_adj, sort(tab$score_adj, decreasing = TRUE),
               tolerance = 1e-12)
})

test_that("se = FALSE returns only the raw score", {
  ## The roxygen used to claim SkymapLinkageTable(skymap) gave "score_raw + dist
  ## only", but the stored SVDs are always present so SE was always computed.
  z <- linkage_fit()
  tab <- SkymapLinkageTable(z$fit, gene.anno = z$sim$anno, se = FALSE,
                            verbose = FALSE)
  expect_false("SE" %in% names(tab))
  expect_true("score_raw" %in% names(tab))
})

test_that("a missing annotation is an explicit error", {
  z <- linkage_fit()
  expect_error(SkymapLinkageTable(z$fit, verbose = FALSE), "is required")
  ## the old default was a hard-coded absolute path on the author's cluster
  expect_error(SkymapLinkageTable(z$fit, gene.anno = "/some/path.rds",
                                  verbose = FALSE), "must be a GRanges")
})

test_that("supplying X and Y warns that it is not the published estimator", {
  z <- linkage_fit()
  expect_warning(
    SkymapLinkageTable(z$fit, gene.anno = z$sim$anno, X = z$sim$X, Y = z$sim$Y,
                       verbose = FALSE),
    "NOT the estimator")
})

test_that("the window bounds the candidate set", {
  z <- linkage_fit()
  narrow <- SkymapLinkageTable(z$fit, gene.anno = z$sim$anno, window = 5e4,
                               se = FALSE, verbose = FALSE)
  wide   <- SkymapLinkageTable(z$fit, gene.anno = z$sim$anno, window = 1e6,
                               se = FALSE, verbose = FALSE)
  expect_lt(nrow(narrow), nrow(wide))
  expect_true(all(narrow$dist <= 5e4))
})

test_that("the jackknife SE matches a brute-force leave-one-cell-out", {
  ## Validates the closed form in se.gene() against the definition, on a small
  ## fit where n refits are affordable.
  set.seed(51)
  sim <- polarisSimulate(n = 60, g = 25, p = 30, r = 3)
  fit <- Polaris(sim$X, sim$Y, gene.chr.ref = sim$anno, x.sds = 0, verbose = FALSE)
  tab <- SkymapLinkageTable(fit, gene.anno = sim$anno, verbose = FALSE)
  tab <- tab[!is.na(tab$SE), ]
  skip_if(nrow(tab) < 1)

  C  <- as.matrix(fit$skymap.cell)
  Xs <- fit$X.svd; Ys <- fit$Y.svd
  n  <- nrow(C)
  K  <- C %*% t(C)

  row <- tab[1, ]
  gi <- match(row$gene, rownames(fit$skymap.P))
  pj <- match(row$peak, rownames(fit$skymap.Q))
  x <- as.numeric(Xs$u %*% (Xs$d * Xs$v[gi, ]))
  y <- as.numeric(Ys$u %*% (Ys$d * Ys$v[pj, ]))

  s <- as.numeric(t(x) %*% K %*% y)
  loo <- vapply(seq_len(n), function(i) {
    keep <- setdiff(seq_len(n), i)
    as.numeric(t(x[keep]) %*% K[keep, keep] %*% y[keep])
  }, numeric(1))
  brute <- sqrt((n - 1) / n * sum((loo - mean(loo))^2))
  expect_equal(row$SE, brute, tolerance = 1e-8)
})

test_that("SkymapZdev standardizes within the documented distance strata", {
  tab <- data.frame(
    gene = rep(c("A", "B"), each = 4),
    peak = sprintf("chr1-%d-%d", seq(1, 71, by = 10), seq(2, 72, by = 10)),
    chr = "chr1",
    dist = c(5e3, 2e4, 3e5, 9e5, 6e3, 2.2e4, 3.2e5, 9.5e5),
    score_raw = c(10, 8, 6, 4, 8, 6, 4, 2),
    SE = rep(2, 8))
  out <- SkymapZdev(tab)
  expect_true(all(c("distbin", "mu0", "z_dev", "z_dev_fl") %in% names(out)))
  expect_equal(nrow(out), nrow(tab))
  ## within each stratum the two pairs average to mu0
  agg <- tapply(out$score_raw, out$distbin, mean)
  expect_equal(out$mu0, as.numeric(agg[as.character(out$distbin)]))
  expect_equal(out$z_dev, (out$score_raw - out$mu0) / out$SE)
  ## the strata are the published ones
  expect_equal(POLARIS_DIST_BREAKS,
               c(0, 1e4, 2.5e4, 5e4, 1e5, 2.5e5, 5e5, 1e6, Inf))
  expect_equal(POLARIS_MIN_DIST, 1e4)
})

test_that("SkymapZdev leaves rows alone and reports missing SE", {
  tab <- data.frame(gene = "A", peak = "chr1-1-2", chr = "chr1",
                    dist = c(2e4, 3e4), score_raw = c(1, 2),
                    SE = c(1, NA))
  out <- SkymapZdev(tab)
  expect_equal(nrow(out), 2)
  expect_true(is.na(out$z_dev[2]))
  expect_error(SkymapZdev(tab[, c("gene", "dist")]), "missing the column")
})

test_that("distances beyond the largest breakpoint are refused", {
  tab <- data.frame(gene = "A", peak = "chr1-1-2", chr = "chr1",
                    dist = 5e6, score_raw = 1, SE = 1)
  expect_no_error(SkymapZdev(tab))   # Inf is the last break, so this is fine
  expect_error(SkymapZdev(tab, breaks = c(0, 1e4, 1e6)), "exceed")
})
