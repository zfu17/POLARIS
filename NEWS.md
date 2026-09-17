# POLARIS 0.99.1

## Compatibility

- POLARIS now installs on R 4.2 and later, down from R 4.4. Nothing in the
  package ever used a feature newer than R 4.0; the old floor simply recorded
  the R the package was built on.

- The two Bioconductor packages of the matrixStats family, `MatrixGenerics` and
  `sparseMatrixStats`, are no longer imported. They supplied one thing, the
  per-feature standard deviations behind the `x.sds` and `y.sds` variance
  floors, and their `colSds` forwards `useNames = NA` on Bioconductor 3.16 and
  earlier, which `matrixStats` 1.0 and later reject. Any fit on R 4.2 or older
  therefore died with `Argument 'useNames' must be either TRUE or FALSE`.
  POLARIS computes the standard deviation directly, so the retained feature set
  no longer depends on which Bioconductor release a user happens to have, and
  the install pulls two fewer Bioconductor packages.

  **Results are unchanged.** Those standard deviations feed only the comparison
  `sd > floor`; the values never enter the decomposition. On the published
  granulocyte-sorted PBMC 10k matrices (10,280 cells; dense X of 36,601 genes,
  sparse Y of 143,883 peaks) not one feature changes its keep-or-drop decision.
  The standard deviations agree with the previous implementation to 1.5e-13
  (dense) and 4.6e-15 (sparse), while the nearest feature to its floor sits
  3.4e-05 and 1.2e-01 away, margins larger by factors of 2e8 and 3e13.

# POLARIS 0.99.0

First packaged release. Previously distributed as a sourced script
(`polaris_function.R`).

Verified **bit-identical** to that script on the PBMC B-cell dataset: the cell
embedding, both paired embeddings, the feature loadings, all singular values,
the concordance score, all 23 per-chromosome decompositions, and every column of
a 1,216,164-row linkage table agree to `0.000e+00`.

## New

- `SkymapZdev()` computes the genomic-distance-aware z-statistic that the
  POLARIS benchmarks rank by. It previously existed only in analysis scripts.
  `POLARIS_DIST_BREAKS` and `POLARIS_MIN_DIST` expose the published strata.
- `SkymapSimScore()` returns `concordance`, the cosine rescaled to `[0,1]`,
  alongside the raw `simscore.cosine`. Rescaling was previously done by hand at
  plotting time.
- `polarisTSSRef()` and `polarisMakeTSSRef()` supply gene annotation, either
  from a bundled hg38 table or built from any TxDb, EnsDb or GTF.
- `RunPolaris()` fits directly from a Seurat object and returns the joint
  embedding as a `DimReduc`.
- `polarisSimulate()` generates paired modalities with a known per-cell
  concordance, for examples and tests.
- `skymap` is now an S3 class with `print` and `summary` methods. It remains a
  plain named list, so existing code that indexes a fit by name is unaffected.
- Fits are reproducible. The truncated SVD of `X` draws a random start, so
  POLARIS seeds it locally and restores the caller's RNG state on exit.
- `input.param` now records the ranks actually selected, not only the
  thresholds.

## Fixes

- `SkymapUMAP()` and `SkymapFindNeighbors()` no longer fail on a bare call. Their
  default `slot` listed `skymap.feature`, which `Polaris()` had stopped
  returning.
- `SkymapUMAP.Chr()` and `SkymapFindNeighbors.Chr()` work again. They assumed a
  combined per-chromosome matrix; `Polaris()` stores `list(P, Q, stdev)`.
- `SkymapSimScore()` defaults to cosine alone, as documented. It previously
  computed all four metrics on a bare call.
- `SkymapSimScore(metric = "cor")` returns one value per cell. Because the
  embeddings were data.frames, it previously returned an all-`NA` result of the
  wrong length.
- `SkymapSimScore()` works with a single component; it previously failed in
  `rowSums()`.
- A flat singular-value spectrum now reports that no spectral gap was found,
  instead of failing with `Error in 1:r : result would be too long a vector`.
- Modalities with 50 or fewer features now fit; the rank request equalled
  `min(dim)`, which the iterative solver rejects.
- A convergence failure in the truncated SVD is surfaced rather than suppressed.
- A failed worker in a parallel step now raises. Errors previously became
  `try-error` objects inside the returned fit, which travelled downstream
  disguised as a successful result.
- `X` and `Y` are checked for matching cell names, not just matching row counts.
- Features off the primary assembly are reported rather than silently dropped.
- No absolute paths remain in any function's defaults.
- The C++ kernels ship in `src/` instead of being compiled by each caller from a
  hard-coded path, and no longer use the removed `Eigen::MappedSparseMatrix`.

## Documentation

- `SkymapLinkageTable()` documents which vectors the standard error actually
  uses: a rank-50 reconstruction, the denoiser's ceiling, not the selected rank
  `r`. Passing `X` and `Y` switches to the raw matrices, a different and
  unpublished estimator, and now warns.
- `x.sds` is documented as a method parameter with its rationale.
