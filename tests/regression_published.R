## Regression test: does the packaged POLARIS reproduce the published fit?
## Compares against the stored pbmc_sorted_10k_B_polaris_output.rds, which was
## produced by runPolaris.R with the pre-package polaris_function.R.
.libPaths(c("/n/holystore01/LABS/xlin/Lab/ziqifu/Rlib_polaris", .libPaths()))
suppressMessages({library(POLARIS); library(Seurat)})

BASE <- "/n/holystore01/LABS/xlin/Lab/ziqifu/POLAR"
obj  <- readRDS(file.path(BASE, "neighborhood_evaluation/prediction_validation/seurat_obj_celltype/pbmc_sorted_10k_B_seurat_filtered.rds"))
old  <- readRDS(file.path(BASE, "neighborhood_evaluation/prediction_validation/skymap/pbmc_sorted_10k_B_polaris_output.rds"))
anno <- readRDS("/n/holystore01/LABS/xlin/Lab/ziqifu/tools/gtf/cellranger_hg38_genetype.rds")

## Exactly as runPolaris.R builds them
X <- t(obj[["RNA"]]$scale.data)
Y <- t(obj[["ATAC"]]$data)
cat(sprintf("input: X %d x %d, Y %d x %d\n", nrow(X), ncol(X), nrow(Y), ncol(Y)))

t0 <- Sys.time()
new <- Polaris(X, Y, gene.chr.ref = anno, mc.cores = 1L, verbose = TRUE)
new <- SkymapSimScore(new)
cat(sprintf("\nfit took %.1f s\n\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

pass <- TRUE
chk <- function(label, ok, detail = "") {
  pass <<- pass && isTRUE(ok)
  cat(sprintf("[%s] %-46s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label, detail))
}

## ---- shapes and feature sets -------------------------------------------------
chk("skymap.cell dim", identical(dim(new$skymap.cell), dim(old$skymap.cell)),
    sprintf("%s vs %s", paste(dim(new$skymap.cell), collapse="x"),
            paste(dim(old$skymap.cell), collapse="x")))
chk("gene set identical", identical(rownames(new$skymap.P), rownames(old$skymap.P)),
    sprintf("%d vs %d genes", nrow(new$skymap.P), nrow(old$skymap.P)))
chk("peak set identical", identical(rownames(new$skymap.Q), rownames(old$skymap.Q)),
    sprintf("%d vs %d peaks", nrow(new$skymap.Q), nrow(old$skymap.Q)))
chk("cell names identical", identical(rownames(new$skymap.cell), rownames(old$skymap.cell)))
chk("r5 (feature rank)", identical(ncol(new$skymap.P), ncol(old$skymap.P)),
    sprintf("%d vs %d", ncol(new$skymap.P), ncol(old$skymap.P)))

## ---- singular values ---------------------------------------------------------
reldiff <- function(a, b) max(abs(a - b) / pmax(abs(b), .Machine$double.eps))
d1 <- reldiff(new$stdev.cell, old$stdev.cell)
chk("stdev.cell rel. diff < 1e-6", d1 < 1e-6, sprintf("%.3e", d1))
d2 <- reldiff(new$stdev.feature, old$stdev.feature)
chk("stdev.feature rel. diff < 1e-6", d2 < 1e-6, sprintf("%.3e", d2))

## ---- embeddings, up to the per-column sign the SVD leaves free ---------------
colcor <- function(A, B) {
  A <- as.matrix(A); B <- as.matrix(B)
  vapply(seq_len(ncol(A)), function(j) abs(stats::cor(A[, j], B[, j])), numeric(1))
}
for (s in c("skymap.U", "skymap.V", "skymap.P", "skymap.Q")) {
  cc <- colcor(new[[s]], old[[s]])
  chk(paste(s, "min |col cor| > 0.9999"), min(cc) > 0.9999,
      sprintf("min %.8f", min(cc)))
}

## ---- the headline statistic --------------------------------------------------
cs <- stats::cor(new$simscore.cosine, old$simscore.cosine)
chk("concordance cor > 0.99999", cs > 0.99999, sprintf("%.10f", cs))
md <- max(abs(new$simscore.cosine - old$simscore.cosine))
chk("concordance max abs diff < 1e-6", md < 1e-6, sprintf("%.3e", md))

## ---- per-chromosome embedding ------------------------------------------------
onull <- vapply(old$skymap.feature.chr, is.null, logical(1))
nnull <- vapply(new$skymap.feature.chr, is.null, logical(1))
chk("chromosomes decomposed",
    identical(sort(names(which(!nnull))), sort(names(which(!onull)))),
    sprintf("new %d, old %d", sum(!nnull), sum(!onull)))

ch <- intersect(names(which(!nnull)), names(which(!onull)))
if (length(ch)) {
  worst <- 0; worstch <- ""
  for (c0 in ch) {
    a <- new$skymap.feature.chr[[c0]]; b <- old$skymap.feature.chr[[c0]]
    if (!identical(dim(a$P), dim(b$P))) { worst <- Inf; worstch <- c0; break }
    r <- reldiff(a$stdev, b$stdev)
    if (r > worst) { worst <- r; worstch <- c0 }
  }
  chk("per-chr stdev rel. diff < 1e-6", worst < 1e-6,
      sprintf("worst %.3e on %s", worst, worstch))
}

cat("\n", if (pass) "ALL CHECKS PASSED" else "SOME CHECKS FAILED", "\n", sep = "")
saveRDS(new, "/n/holystore01/LABS/xlin/Lab/ziqifu/Rlib_polaris/new_B_skymap.rds")
