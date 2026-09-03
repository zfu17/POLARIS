## Definitive reproduction test: the packaged POLARIS versus the CURRENT
## polaris_function.R, on the same input with the same RNG start.
##
## This is the test that matters. The stored *_polaris_output.rds artifacts date
## from Jun 6 2026, before the Jun 21 QR optimization of the feature embedding,
## so they carry a larger r5 (19 vs 15) whose extra components are null-space
## noise. Comparing against them measures that change, not this port.
.libPaths(c("/n/holystore01/LABS/xlin/Lab/ziqifu/Rlib_polaris", .libPaths()))

SCRIPTS <- "/n/holystore01/LABS/xlin/Lab/ziqifu/POLAR/scripts"
BASE    <- "/n/holystore01/LABS/xlin/Lab/ziqifu/POLAR"

## The pre-package core calls many symbols bare, so its driver's preamble has to
## be reproduced before it can be sourced at all.
suppressMessages({
  library(dplyr); library(glue); library(stringr); library(tibble)
  library(Matrix); library(RSpectra); library(irlba); library(parallel)
  library(matrixStats); library(sparseMatrixStats); library(SparseArray)
  library(GenomicRanges); library(uwot); library(BiocNeighbors)
  library(Seurat)
})
Rcpp::sourceCpp(file.path(SCRIPTS, "RCPP/mat_eigen.cpp"))

## Source the old core into its own environment so its Polaris() does not
## collide with the package's.
old.env <- new.env()
source(file.path(SCRIPTS, "polaris_function.R"), local = old.env)
OldPolaris        <- get("Polaris", envir = old.env)
OldSkymapSimScore <- get("SkymapSimScore", envir = old.env)

obj  <- readRDS(file.path(BASE, "neighborhood_evaluation/prediction_validation/seurat_obj_celltype/pbmc_sorted_10k_B_seurat_filtered.rds"))
anno <- readRDS("/n/holystore01/LABS/xlin/Lab/ziqifu/tools/gtf/cellranger_hg38_genetype.rds")
X <- t(obj[["RNA"]]$scale.data)
Y <- t(obj[["ATAC"]]$data)
cat(sprintf("input: X %d x %d, Y %d x %d\n\n", nrow(X), ncol(X), nrow(Y), ncol(Y)))

## The package seeds irlba internally with seed = 1; give the old code the same
## RNG start so any remaining difference is attributable to the port itself.
set.seed(1)
t0 <- Sys.time()
old <- OldPolaris(X, Y, gene.chr.ref = as.data.frame(anno), mc.cores = 1)
old <- OldSkymapSimScore(old, metric = "cosine")
cat(sprintf("\nold script: %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

t0 <- Sys.time()
new <- POLARIS::Polaris(X, Y, gene.chr.ref = anno, mc.cores = 1, verbose = FALSE)
new <- POLARIS::SkymapSimScore(new)
cat(sprintf("package:    %.1f s\n\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

pass <- TRUE
chk <- function(label, ok, detail = "") {
  pass <<- pass && isTRUE(ok)
  cat(sprintf("[%s] %-44s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label, detail))
}
maxad <- function(a, b) max(abs(as.numeric(as.matrix(a)) - as.numeric(as.matrix(b))))

chk("gene set identical", identical(rownames(new$skymap.P), rownames(old$skymap.P)),
    sprintf("%d genes", nrow(new$skymap.P)))
chk("peak set identical", identical(rownames(new$skymap.Q), rownames(old$skymap.Q)),
    sprintf("%d peaks", nrow(new$skymap.Q)))
chk("ranks identical (r3, r4, r5)",
    ncol(new$skymap.U) == ncol(old$skymap.U) &&
    ncol(new$skymap.V) == ncol(old$skymap.V) &&
    ncol(new$skymap.P) == ncol(old$skymap.P),
    sprintf("r3=%d r4=%d r5=%d", ncol(new$skymap.U), ncol(new$skymap.V),
            ncol(new$skymap.P)))

for (s in c("skymap.cell", "skymap.U", "skymap.V", "skymap.P", "skymap.Q")) {
  if (!identical(dim(as.matrix(new[[s]])), dim(as.matrix(old[[s]])))) {
    chk(paste(s, "dims"), FALSE, "dimension mismatch"); next
  }
  d <- maxad(new[[s]], old[[s]])
  chk(paste(s, "max abs diff < 1e-10"), d < 1e-10, sprintf("%.3e", d))
}
for (s in c("stdev.cell", "stdev.feature")) {
  d <- maxad(new[[s]], old[[s]])
  chk(paste(s, "max abs diff < 1e-10"), d < 1e-10, sprintf("%.3e", d))
}

d <- maxad(new$simscore.cosine, old$simscore.cosine)
chk("concordance (raw cosine) identical", d < 1e-12, sprintf("%.3e", d))
chk("concordance rescaled to [0,1]",
    all(new$concordance >= 0 & new$concordance <= 1) &&
    maxad(new$concordance, (old$simscore.cosine + 1) / 2) < 1e-12)

on <- names(which(!vapply(old$skymap.feature.chr, is.null, logical(1))))
nn <- names(which(!vapply(new$skymap.feature.chr, is.null, logical(1))))
chk("same chromosomes decomposed", identical(sort(on), sort(nn)),
    sprintf("%d vs %d", length(nn), length(on)))
worst <- 0; worstch <- ""
for (c0 in intersect(on, nn)) {
  a <- new$skymap.feature.chr[[c0]]; b <- old$skymap.feature.chr[[c0]]
  if (!identical(dim(a$P), dim(b$P)) || !identical(dim(a$Q), dim(b$Q))) {
    worst <- Inf; worstch <- c0; break
  }
  d <- max(maxad(a$P, b$P), maxad(a$Q, b$Q), maxad(a$stdev, b$stdev))
  if (d > worst) { worst <- d; worstch <- c0 }
}
chk("per-chr P/Q/stdev max abs diff < 1e-10", worst < 1e-10,
    sprintf("worst %.3e on %s", worst, worstch))

## Linkage table: the published estimator (no X/Y passed)
told <- get("SkymapLinkageTable", envir = old.env)(
  old, gene.anno = anno, window = 1e6, mc.cores = 1, verbose = FALSE)
tnew <- POLARIS::SkymapLinkageTable(
  new, gene.anno = anno, window = 1e6, mc.cores = 1, verbose = FALSE)
chk("linkage table same n rows", nrow(tnew) == nrow(told),
    sprintf("%d vs %d", nrow(tnew), nrow(told)))
if (nrow(tnew) == nrow(told)) {
  key <- function(t) paste(t$gene, t$peak)
  o <- match(key(tnew), key(told))
  chk("linkage pairs identical", !anyNA(o))
  if (!anyNA(o)) {
    for (cc in c("score_raw", "SE", "score_adj", "dist")) {
      d <- max(abs(tnew[[cc]] - told[[cc]][o]), na.rm = TRUE)
      chk(paste("linkage", cc, "max abs diff < 1e-9"), d < 1e-9, sprintf("%.3e", d))
    }
  }
}

cat("\n", if (pass) "REPRODUCTION VERIFIED" else "DIFFERENCES FOUND", "\n", sep = "")
