## Gen 1 (2026-06-08) versus Gen 3 (current, = the package) on the same input.
##
## Why this test exists. The provenance audit in
## manuscript/PACKAGE_OPEN_QUESTIONS.md establishes two things:
##   * every manuscript figure and artifact dates from 2026-06-06..06-10, i.e.
##     Gen 1 (polaris_function.R.bak_2026-06-21, mtime Jun 8);
##   * the package is bit-identical to Gen 3 (regression_vs_script.R).
## The missing link is Gen 1 == Gen 3. If they agree, the published figures are
## reproducible by the released package and Code availability is unproblematic.
## If they do not, that has to be known before submission.
##
## Expected differences, by construction rather than by error:
##   * Gen 1 denoises X with RSpectra::svds; Gen 3 uses irlba. Both are exact to
##     solver tolerance, so agreement should be ~1e-8, not bit-level.
##   * Gen 1 takes svds(dense n x n product, k = 50), so prod.svd has 50
##     components; Gen 3's low-rank identity yields only min(r1, r2). Only the
##     leading max(r3, r4) are ever used, so stdev.cell should still agree.
##   * r5 WILL differ: Gen 1's dense genes x peaks path can return up to 50
##     components, Gen 3's QR path caps at r3 + r4. This is the Jun 21 change.
.libPaths(c("/n/holystore01/LABS/xlin/Lab/ziqifu/Rlib_polaris", .libPaths()))

SCRIPTS <- "/n/holystore01/LABS/xlin/Lab/ziqifu/POLAR/scripts"
BASE    <- "/n/holystore01/LABS/xlin/Lab/ziqifu/POLAR"

suppressMessages({
  library(dplyr); library(glue); library(stringr); library(tibble)
  library(Matrix); library(RSpectra); library(irlba); library(parallel)
  library(matrixStats); library(sparseMatrixStats); library(SparseArray)
  library(GenomicRanges); library(uwot); library(BiocNeighbors); library(Seurat)
})
Rcpp::sourceCpp(file.path(SCRIPTS, "RCPP/mat_eigen.cpp"))

gen1 <- new.env()
source(file.path(SCRIPTS, "polaris_function.R.bak_2026-06-21"), local = gen1)

obj  <- readRDS(file.path(BASE, "neighborhood_evaluation/prediction_validation/seurat_obj_celltype/pbmc_sorted_10k_B_seurat_filtered.rds"))
anno <- readRDS("/n/holystore01/LABS/xlin/Lab/ziqifu/tools/gtf/cellranger_hg38_genetype.rds")
X <- t(obj[["RNA"]]$scale.data)
Y <- t(obj[["ATAC"]]$data)
cat(sprintf("input: X %d x %d, Y %d x %d\n\n", nrow(X), ncol(X), nrow(Y), ncol(Y)))

t0 <- Sys.time()
g1 <- get("Polaris", envir = gen1)(X, Y, gene.chr.ref = as.data.frame(anno), mc.cores = 1)
g1 <- get("SkymapSimScore", envir = gen1)(g1, metric = "cosine")
cat(sprintf("\nGen 1: %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

t0 <- Sys.time()
g3 <- POLARIS::Polaris(X, Y, gene.chr.ref = anno, mc.cores = 1, verbose = FALSE)
g3 <- POLARIS::SkymapSimScore(g3)
cat(sprintf("Gen 3 (package): %.1f s\n\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

pass <- TRUE
chk <- function(label, ok, detail = "") {
  pass <<- pass && isTRUE(ok)
  cat(sprintf("[%s] %-48s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label, detail))
}
note <- function(label, detail) cat(sprintf("[NOTE] %-48s %s\n", label, detail))
colcor <- function(A, B) {
  A <- as.matrix(A); B <- as.matrix(B)
  k <- min(ncol(A), ncol(B))
  vapply(seq_len(k), function(j) abs(stats::cor(A[, j], B[, j])), numeric(1))
}
maxad <- function(a, b) {
  a <- as.numeric(as.matrix(a)); b <- as.numeric(as.matrix(b))
  k <- min(length(a), length(b)); max(abs(a[seq_len(k)] - b[seq_len(k)]))
}

## ---- the feature space must be identical: same filter, same input ----------
chk("gene set identical", identical(rownames(g3$skymap.P), rownames(g1$skymap.P)),
    sprintf("%d vs %d", nrow(g3$skymap.P), nrow(g1$skymap.P)))
chk("peak set identical", identical(rownames(g3$skymap.Q), rownames(g1$skymap.Q)),
    sprintf("%d vs %d", nrow(g3$skymap.Q), nrow(g1$skymap.Q)))

## ---- cell-side ranks and embedding -----------------------------------------
chk("r3 identical", ncol(g3$skymap.U) == ncol(g1$skymap.U),
    sprintf("%d vs %d", ncol(g3$skymap.U), ncol(g1$skymap.U)))
chk("r4 identical", ncol(g3$skymap.V) == ncol(g1$skymap.V),
    sprintf("%d vs %d", ncol(g3$skymap.V), ncol(g1$skymap.V)))
note("r5 (expected to differ: Jun 21 QR change)",
     sprintf("Gen3 %d vs Gen1 %d", ncol(g3$skymap.P), ncol(g1$skymap.P)))

d <- maxad(g3$stdev.cell, g1$stdev.cell)
chk("stdev.cell agrees < 1e-6", d < 1e-6, sprintf("%.3e", d))

for (s in c("skymap.U", "skymap.V")) {
  cc <- colcor(g3[[s]], g1[[s]])
  chk(paste(s, "min |col cor| > 0.9999"), min(cc) > 0.9999, sprintf("min %.10f", min(cc)))
}

## ---- the headline statistic -------------------------------------------------
cs <- stats::cor(g3$simscore.cosine, g1$simscore.cosine)
chk("concordance cor > 0.999999", cs > 0.999999, sprintf("%.12f", cs))
d <- max(abs(g3$simscore.cosine - g1$simscore.cosine))
chk("concordance max abs diff < 1e-6", d < 1e-6, sprintf("%.3e", d))

## ---- feature loadings, over the components both generations retain ---------
kP <- min(ncol(g3$skymap.P), ncol(g1$skymap.P))
cc <- colcor(g3$skymap.P, g1$skymap.P)
chk(sprintf("skymap.P first %d cols |cor| > 0.9999", kP), min(cc) > 0.9999,
    sprintf("min %.10f", min(cc)))
cc <- colcor(g3$skymap.Q, g1$skymap.Q)
chk(sprintf("skymap.Q first %d cols |cor| > 0.9999", kP), min(cc) > 0.9999,
    sprintf("min %.10f", min(cc)))

## ---- THE decisive comparison: the published linkage numbers ----------------
t1 <- get("SkymapLinkageTable", envir = gen1)(
  g1, gene.anno = anno, window = 1e6, mc.cores = 1, verbose = FALSE)
t3 <- POLARIS::SkymapLinkageTable(
  g3, gene.anno = anno, window = 1e6, mc.cores = 1, verbose = FALSE)

chk("linkage table same n rows", nrow(t3) == nrow(t1),
    sprintf("%d vs %d", nrow(t3), nrow(t1)))
if (nrow(t3) == nrow(t1)) {
  o <- match(paste(t3$gene, t3$peak), paste(t1$gene, t1$peak))
  chk("linkage pairs identical", !anyNA(o))
  if (!anyNA(o)) {
    for (cc in c("dist", "score_raw", "SE", "score_adj")) {
      a <- t3[[cc]]; b <- t1[[cc]][o]
      ok <- is.finite(a) & is.finite(b)
      ad <- max(abs(a[ok] - b[ok]))
      rr <- stats::cor(a[ok], b[ok])
      chk(sprintf("linkage %-9s cor > 0.9999", cc), rr > 0.9999,
          sprintf("cor %.10f, max abs diff %.3e", rr, ad))
    }
    ## Ranking is what the benchmarks consume, so compare the ORDER too.
    z3 <- POLARIS::SkymapZdev(t3); z1 <- t1[o, ]
    z1$distbin <- cut(z1$dist, POLARIS::POLARIS_DIST_BREAKS,
                      include.lowest = TRUE, right = FALSE)
    mu0 <- tapply(z1$score_raw, z1$distbin, mean, na.rm = TRUE)
    z1$z_dev <- (z1$score_raw - as.numeric(mu0[as.character(z1$distbin)])) / z1$SE
    ok <- is.finite(z3$z_dev) & is.finite(z1$z_dev)
    rr <- stats::cor(z3$z_dev[ok], z1$z_dev[ok], method = "spearman")
    chk("z_dev Spearman rank cor > 0.999", rr > 0.999, sprintf("rho %.10f", rr))
  }
}

cat("\n", if (pass) "GEN 1 AND GEN 3 AGREE" else "GEN 1 AND GEN 3 DIFFER", "\n", sep = "")
