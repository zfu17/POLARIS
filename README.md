# POLARIS

Concordance-aware integration of paired single-cell multiomic data.

POLARIS takes two matrices measured on the same cells and returns three things
from one closed-form spectral decomposition:

- a **per-cell concordance score**, saying how far the two modalities agree
  about each individual cell;
- a **joint cell embedding** that preserves structure supported by *either*
  modality, rather than only what they share;
- a **cross-modal feature embedding** that scores candidate feature pairs
  (regulatory element to gene, protein to gene) with a standard error attached
  to every pair.

It runs on CPU, needs no pretraining or accelerator, and scales to tens of
thousands of cells in minutes.

## Installation

POLARIS imports several Bioconductor packages, so install it through
BiocManager, which resolves CRAN and Bioconductor dependencies together:

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install("zfu17/POLARIS")
```

The package contains compiled code, so a C++ toolchain is required: Rtools on
Windows, the Xcode command line tools on macOS. Linux systems generally have one
already. Add `build_vignettes = TRUE` to also build the vignette locally, which
needs pandoc 2.8 or newer; RStudio bundles a suitable one.

## Quick start

```r
library(POLARIS)

# X and Y are cells x features, with the same cells in the same order.
# For scMultiome: X = t(RNA scale.data), Y = t(TF-IDF ATAC data).
anno <- polarisTSSRef("hg38")
fit  <- Polaris(X, Y, gene.chr.ref = anno, mc.cores = 4)

# per-cell concordance, rescaled to [0,1]
fit <- SkymapSimScore(fit)
summary(fit$concordance)

# cis gene-peak linkage with per-link standard errors
tab <- SkymapLinkageTable(fit, gene.anno = anno, window = 1e6, mc.cores = 4)

# rank by the distance-aware z-statistic, not by raw magnitude
tab <- SkymapZdev(tab)
tab <- tab[tab$dist >= POLARIS_MIN_DIST, ]
head(tab[order(-tab$z_dev), c("gene", "peak", "dist", "z_dev")])
```

From a Seurat object, `RunPolaris()` assembles the arguments and returns the
joint embedding as a `DimReduc`:

```r
res <- RunPolaris(obj, assay.x = "RNA", assay.y = "ATAC", mc.cores = 4)
```

See `vignette("POLARIS")` for the full workflow, the preprocessing conventions,
and the CITE-seq and cell-type-specific paths. The vignette is installed when
you pass `build_vignettes = TRUE`; its source is `vignettes/POLARIS.Rmd`.

## Two things worth reading before your first run

**Orientation.** `X` and `Y` are **cells in rows**. This is the transpose of how
Seurat stores assays.

**The `x.sds` floor.** `Polaris()` keeps features of `X` whose standard
deviation exceeds `x.sds` (default `0.95`). This deliberately selects a broader
feature set than a highly-variable-gene list, but it assumes `X` is z-scored.
If you pass a log-normalized matrix instead, almost every feature falls below
the floor and the fit stops with an empty `X`. Pass `x.sds = 0` in that case.

## Citation

See `citation("POLARIS")`.

## License

MIT
