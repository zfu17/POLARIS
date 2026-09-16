## Does the package's reconcile.ranks=FALSE path reproduce a PUBLISHED unreconciled fit?
## Technique from the 2026-09-10 audit: a fit can be classified from its OWN stored
## X.svd/Y.svd by rebuilding the product SVD both ways and comparing to stdev.cell.
## No Seurat object and no refit required.
.libPaths(c("/n/holystore01/LABS/xlin/Lab/ziqifu/Rlib_polaris", .libPaths()))
B <- "/n/holystore01/LABS/xlin/Lab/ziqifu/POLAR"

## Exactly the joint-embedding block of R/polaris.R, parameterised by convention.
rebuild <- function(Xs, Ys, r1, r2) {
  U1 <- Xs$u[, seq_len(r1), drop = FALSE]
  U2 <- Ys$u[, seq_len(r2), drop = FALSE]
  core <- svd((Xs$d[seq_len(r1)] * crossprod(U1, U2)) * rep(Ys$d[seq_len(r2)], each = r1))
  core$d
}
relerr <- function(a, b) { k <- min(length(a), length(b))
  max(abs(a[seq_len(k)] - b[seq_len(k)])) / max(abs(b[seq_len(k)])) }

fits <- c(
  `pbmc_sorted_10k (Fig 2/3)`   = file.path(B, "data_proc/10XPBMC/pbmc_sorted_10k_polaris_output.rds"),
  `pbmc_unsorted_10k (Fig 2/3)` = file.path(B, "data_proc/10XPBMC/pbmc_unsorted_10k_polaris_output.rds"),
  `pbmc_sorted_3k (Supp 2/3)`   = file.path(B, "data_proc/10XPBMC/pbmc_sorted_3k_polaris_output.rds"),
  `B cell (per-cell-type)`      = file.path(B, "neighborhood_evaluation/prediction_validation/skymap/pbmc_sorted_10k_B_polaris_output.rds"),
  `CD8 T (per-cell-type)`       = file.path(B, "neighborhood_evaluation/prediction_validation/skymap/pbmc_sorted_10k_CD8_T_polaris_output.rds"))

cat(sprintf("%-30s %5s %5s %12s %12s  %s\n",
            "fit", "r1", "r2", "unreconciled", "collapsed", "verdict"))
for (nm in names(fits)) {
  f <- fits[[nm]]; if (!file.exists(f)) { cat(nm, "MISSING\n"); next }
  sk <- readRDS(f); Xs <- sk$X.svd; Ys <- sk$Y.svd
  target <- as.numeric(sk$stdev.cell)
  ## recover the per-modality ranks the ratio rule would have chosen
  rr <- function(d, T) { r <- d[-length(d)] / d[-1L]; max(which(r > T)) }
  r1 <- rr(Xs$d, 1.01); r2 <- rr(Ys$d, 1.01); m <- min(r1, r2)
  e_un <- relerr(rebuild(Xs, Ys, r1, r2), target)
  e_co <- relerr(rebuild(Xs, Ys, m,  m ), target)
  verdict <- if (e_un < 1e-8 && e_co > 1e-3) "UNRECONCILED" else
             if (e_co < 1e-8 && e_un > 1e-3) "COLLAPSED" else
             if (r1 == r2) "identical (r1==r2)" else "AMBIGUOUS"
  cat(sprintf("%-30s %5d %5d %12.2e %12.2e  %s\n", nm, r1, r2, e_un, e_co, verdict))
  rm(sk); gc(verbose = FALSE)
}
