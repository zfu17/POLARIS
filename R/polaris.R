#' Fit POLARIS to two paired modalities
#'
#' Performs the POLARIS spectral integration of two matrices measured on the
#' same cells and returns a [skymap] object carrying the joint cell embedding,
#' the paired cell embeddings used for the concordance score, and the cross-modal
#' feature embedding.
#'
#' @param X Numeric matrix for modality 1, **cells in rows, features in
#'   columns**, centered and scaled. For scMultiome this is the RNA
#'   `scale.data`, transposed.
#' @param Y Numeric matrix for modality 2, cells in rows, features in columns,
#'   with the same cells in the same order as `X`. For scMultiome this is the
#'   TF-IDF normalized ATAC `data`, transposed; for CITE-seq the CLR-normalized
#'   ADT `data`, transposed.
#' @param T1,T2 Singular-value ratio thresholds of the denoiser for `X` and `Y`.
#'   Components are retained while the spectrum still decays and truncated once
#'   it flattens. Defaults `1.01`.
#' @param T3 Ratio threshold for the feature-embedding rank `r5`. Default `1.1`.
#' @param r1,r2 Ranks retained for `X` and `Y`. Derived from `T1` and `T2` unless
#'   given. Note that the two modality ranks are set to their common minimum, so
#'   supplying unequal values has no effect beyond the smaller of the two.
#' @param r3,r4 Dimensions of \eqn{\hat{U}} and \eqn{\hat{V}}. Default
#'   `min(r1, r2)`.
#' @param r5 Dimension of `P` and `Q`. Derived from `T3` unless given.
#' @param gene.chr.ref A `GRanges` with a `gene_name` metadata column and
#'   strand, giving the genomic location of the features of `X`. Required for the
#'   per-chromosome feature decomposition that [SkymapLinkageTable()] scores. If
#'   `NULL` the per-chromosome step is skipped with a message and the fit still
#'   returns the genome-wide `skymap.P` and `skymap.Q`, which is the right
#'   behaviour for modality pairs without genomic coordinates such as CITE-seq.
#'   See [polarisTSSRef()].
#' @param chrs Chromosomes to decompose. Default `NULL` uses those actually
#'   present in both `gene.chr.ref` and the feature names of `Y`.
#' @param x.sds,y.sds Standard-deviation floors for retaining features of `X`
#'   and `Y`. A feature is kept when its standard deviation exceeds the floor.
#'   The defaults are deliberately asymmetric: `x.sds = 0.95` keeps genes that
#'   retain appreciable variance, which for globally scaled data subset to one
#'   cell population is a broader and more permissive selection than a
#'   highly-variable-gene list, while `y.sds = 0` removes only constant features
#'   of the second modality. Set `x.sds = 0` for modality pairs where `X` is not
#'   z-scored, as the CITE-seq analyses do.
#' @param mc.cores Cores for the per-chromosome decomposition.
#' @param seed Seed for the starting vector of the `X` decomposition. The
#'   truncated SVD used for `X` draws a random start, so without a fixed seed
#'   two fits of the same data differ slightly (about `1e-9` in the singular
#'   values) and the caller's random stream is advanced as a side effect.
#'   POLARIS therefore fixes the seed locally and restores the previous RNG state
#'   afterwards, so a fit is bit-reproducible and leaves your random stream
#'   untouched. Pass `NULL` for the older non-deterministic behaviour.
#' @param verbose Report progress.
#'
#' @return A [skymap] object.
#'
#' @details
#' The optimizer of the POLARIS objective is obtained from the SVD of
#' \eqn{(XX^\top)^{1/2}(YY^\top)^{1/2}}. Rather than forming those `n` x `n`
#' operators, this implementation uses the exact identity
#' \deqn{(XX^\top)^{1/2}(YY^\top)^{1/2} = U_1 \left[D_1 (U_1^\top U_2) D_2\right] U_2^\top,}
#' where \eqn{X = U_1 D_1 V_1^\top} and \eqn{Y = U_2 D_2 V_2^\top} are the thin
#' SVDs of the two denoised modalities. Since \eqn{U_1} and \eqn{U_2} have
#' orthonormal columns, the SVD of the small bracketed `r` x `r` matrix yields
#' the SVD of the full operator, reducing an `n` x `n` problem to an `r` x `r`
#' one. The feature-loading product is likewise
#' decomposed through its thin QR factors instead of materializing a dense
#' genes-by-features matrix. Both give the identical decomposition.
#'
#' Ranks are selected by the ratio rule
#' \eqn{k = \max\{j : \sigma_j / \sigma_{j+1} > T\}}.
#'
#' @examples
#' set.seed(1)
#' sim <- polarisSimulate(n = 200, g = 60, p = 80, r = 5)
#' fit <- Polaris(sim$X, sim$Y, x.sds = 0, verbose = FALSE)
#' fit
#'
#' @seealso [SkymapSimScore()], [SkymapLinkageTable()]
#' @export
Polaris <- function(X, Y,
                    T1 = 1.01, T2 = 1.01, T3 = 1.1,
                    r1 = NULL, r2 = NULL, r3 = NULL, r4 = NULL, r5 = NULL,
                    gene.chr.ref = NULL, chrs = NULL,
                    x.sds = 0.95, y.sds = 0,
                    mc.cores = 1L, seed = 1L, verbose = TRUE) {

  say <- function(...) if (isTRUE(verbose)) message(...)

  ## ---- cell conformability -------------------------------------------------
  if (nrow(X) != nrow(Y))
    stop(sprintf(paste0("X and Y must have the same number of cells (rows): ",
                        "X has %d, Y has %d. Both must be cells x features."),
                 nrow(X), nrow(Y)), call. = FALSE)

  ## Equal row COUNTS are not enough: the whole method assumes row i of X and
  ## row i of Y are the same cell. A reordered Y would silently mispair the
  ## modalities and produce a plausible but meaningless fit.
  if (!is.null(rownames(X)) && !is.null(rownames(Y))) {
    if (!identical(rownames(X), rownames(Y)))
      stop(paste0("X and Y have the same number of cells but different or ",
                  "differently ordered cell names. Reorder Y to match X, ",
                  "e.g. Y <- Y[rownames(X), , drop = FALSE]."), call. = FALSE)
  } else {
    warning(paste0("X and/or Y have no cell names, so the pairing of rows ",
                   "cannot be verified. POLARIS assumes row i of X and row i ",
                   "of Y are the same cell."), call. = FALSE)
  }
  if (anyDuplicated(rownames(X)))
    stop("X has duplicated cell names; deduplicate before fitting.", call. = FALSE)

  ## ---- feature filtering ---------------------------------------------------
  ## The pre-package code round-tripped both modalities through
  ## SVT_SparseMatrix before computing column standard deviations. That was
  ## unnecessary: MatrixGenerics::colSds dispatches on base matrices and on
  ## dgCMatrix alike and returns identical values, so the conversion only cost
  ## memory. It also relied on `colSds` resolving to the right namespace by
  ## library() ordering, since matrixStats::colSds is not generic and errors on
  ## any S4 matrix.
  keep.x <- .col_sds(X) > x.sds
  keep.y <- .col_sds(Y) > y.sds
  if (!any(keep.x))
    stop(sprintf(paste0("No feature of X has standard deviation above x.sds = %g, ",
                        "so X is empty. If X is not z-scored (for example ",
                        "log-normalized rather than scaled) pass x.sds = 0."),
                 x.sds), call. = FALSE)
  if (!any(keep.y))
    stop(sprintf("No feature of Y has standard deviation above y.sds = %g.",
                 y.sds), call. = FALSE)
  say(sprintf("Retaining %d of %d features of X (sd > %g) and %d of %d of Y (sd > %g).",
              sum(keep.x), length(keep.x), x.sds,
              sum(keep.y), length(keep.y), y.sds))
  X <- X[, keep.x, drop = FALSE]
  Y <- Y[, keep.y, drop = FALSE]

  ## dgCMatrix is preferred, but its non-zero count must fit in an integer, so a
  ## large dense modality falls back to a base matrix.
  tmp <- try(methods::as(X, "dgCMatrix"), silent = TRUE)
  X <- if (inherits(tmp, "try-error")) as.matrix(X) else tmp
  tmp <- try(methods::as(Y, "dgCMatrix"), silent = TRUE)
  Y <- if (inherits(tmp, "try-error")) as.matrix(Y) else tmp

  ## ---- denoise each modality ----------------------------------------------
  ## Solver choice is measured, not assumed. X is effectively dense (scale.data),
  ## which suits irlba; Y is genuinely sparse TF-IDF, where RSpectra wins. Both
  ## return the same subspace (minCos 1.000000, d1 to 4e-16 at n = 35,000).
  say("Denoising X and Y.")
  kx <- .svd_k(X)
  ky <- .svd_k(Y)
  X.svd <- .thin_svd(X, kx, solver = "irlba", seed = seed)
  Y.svd <- .thin_svd(Y, ky, solver = "svds")

  r1.auto <- .rank_from_ratio(X.svd$d, T1, "X")
  r2.auto <- .rank_from_ratio(Y.svd$d, T2, "Y")
  if (is.null(r1)) r1 <- r1.auto
  if (is.null(r2)) r2 <- r2.auto

  ## The two modality ranks are set to their common minimum (POLARIS Methods,
  ## Rank selection). Warn rather than silently overriding a user's choice.
  if (r1 != r2) {
    say(sprintf("Modality ranks r1 = %d and r2 = %d set to their common minimum, %d.",
                r1, r2, min(r1, r2)))
    r1 <- r2 <- min(r1, r2)
  }
  if (r1 > length(X.svd$d) || r2 > length(Y.svd$d))
    stop(sprintf("Requested rank %d exceeds the %d components computed.",
                 max(r1, r2), min(length(X.svd$d), length(Y.svd$d))), call. = FALSE)

  ## ---- joint cell embedding -----------------------------------------------
  ## The operator being decomposed is (XX')^(1/2) (YY')^(1/2), which is the
  ## objective in the POLARIS Methods. With X = U1 D1 V1' and Y = U2 D2 V2',
  ##   (XX')^(1/2) (YY')^(1/2) = U1 D1 U1' U2 D2 U2' = U1 [D1 (U1'U2) D2] U2'.
  ## U1 and U2 have orthonormal columns, so svd(core) = W S Z' gives the SVD of
  ## the full operator as (U1 W) S (U2 Z)', from an r x r problem instead of
  ## three dense n x n matrices.
  ##
  ## NOTE the exponents. The comment in the pre-package polaris_function.R wrote
  ## this identity as "XX'YY' = U1 D1 (U1'U2) D2 U2'", which is not true: XX'YY'
  ## would carry D1^2 and D2^2 and its singular values are ~80x larger. The CODE
  ## was always right (it uses D1 and D2, i.e. the square roots) and matches the
  ## Methods; only the comment was wrong. Verified numerically: the reconstruction
  ## matches (XX')^(1/2)(YY')^(1/2) to 3.6e-15 and differs from XX'YY' by O(600).
  say("Computing U and V.")
  U1 <- X.svd$u[, seq_len(r1), drop = FALSE]
  U2 <- Y.svd$u[, seq_len(r2), drop = FALSE]
  core <- svd((X.svd$d[seq_len(r1)] * crossprod(U1, U2)) *
                rep(Y.svd$d[seq_len(r2)], each = r1))
  prod.svd <- list(u = U1 %*% core$u, v = U2 %*% core$v, d = core$d)
  rm(U1, U2, core)

  if (is.null(r3)) r3 <- min(r1, r2)
  if (is.null(r4)) r4 <- min(r1, r2)
  ## r3 and r4 index prod.svd, which has only min(r1, r2) components. Unchecked,
  ## an over-large r3 produced NA columns rather than an error.
  avail <- length(prod.svd$d)
  if (r3 > avail || r4 > avail)
    stop(sprintf(paste0("r3 = %d and r4 = %d cannot exceed the %d components ",
                        "available from the joint decomposition (min(r1, r2))."),
                 r3, r4, avail), call. = FALSE)
  if (r3 < 1L || r4 < 1L)
    stop("r3 and r4 must be at least 1.", call. = FALSE)

  skymap.cell <- cbind(prod.svd$u[, seq_len(r3), drop = FALSE],
                       prod.svd$v[, seq_len(r4), drop = FALSE])
  dimnames(skymap.cell) <- list(rownames(X), paste0("skyPC_", seq_len(r3 + r4)))

  ## Kept as matrices, not data.frames. The original stored U and V as
  ## data.frames, which made SkymapSimScore(metric = "cor") return an all-NA
  ## matrix of the wrong length, because U[i, ] on a data.frame stays a one-row
  ## data.frame instead of collapsing to a vector. Every consumer in the
  ## analysis scripts already wrapped these in as.matrix().
  skymap.U <- prod.svd$u[, seq_len(r3), drop = FALSE]
  dimnames(skymap.U) <- list(rownames(X), paste0("uPC_", seq_len(r3)))
  skymap.V <- prod.svd$v[, seq_len(r4), drop = FALSE]
  dimnames(skymap.V) <- list(rownames(X), paste0("vPC_", seq_len(r4)))

  stdev.cell <- prod.svd$d[seq_len(max(r3, r4))]

  ## ---- feature loadings ---------------------------------------------------
  inf.emb.X <- if (methods::is(X, "sparseMatrix"))
    matrixMultSparseDense(Matrix::t(X), skymap.cell)
  else matrixMultiplyEigen(t(X), skymap.cell)
  rownames(inf.emb.X) <- colnames(X)

  inf.emb.Y <- if (methods::is(Y, "sparseMatrix"))
    matrixMultSparseDense(Matrix::t(Y), skymap.cell)
  else matrixMultiplyEigen(t(Y), skymap.cell)
  rownames(inf.emb.Y) <- colnames(Y)

  ## ---- joint feature embedding --------------------------------------------
  ## XwY = inf.emb.X %*% t(inf.emb.Y) has rank <= ncol(inf.emb.X); take its exact
  ## SVD from the thin QR factors rather than materializing the dense
  ## genes x features matrix.
  say("Computing P and Q.")
  qa <- qr(inf.emb.X)
  qb <- qr(inf.emb.Y)
  core.fe <- svd(qr.R(qa) %*% t(qr.R(qb)))
  XwY.svd <- list(u = qr.Q(qa) %*% core.fe$u,
                  v = qr.Q(qb) %*% core.fe$v,
                  d = core.fe$d)
  if (is.null(r5)) r5 <- .rank_from_ratio(XwY.svd$d, T3, "the feature embedding")
  r5 <- min(r5, length(XwY.svd$d))

  skymap.P <- XwY.svd$u[, seq_len(r5), drop = FALSE]
  dimnames(skymap.P) <- list(colnames(X), paste0("pPC_", seq_len(r5)))
  skymap.Q <- XwY.svd$v[, seq_len(r5), drop = FALSE]
  dimnames(skymap.Q) <- list(colnames(Y), paste0("qPC_", seq_len(r5)))
  stdev.feature <- XwY.svd$d[seq_len(r5)]
  ## The full spectrum is kept as well (at most r3 + r4 values). The pre-package
  ## fit discarded everything past r5, which made the denoiser spectrum figure
  ## impossible to redraw from a saved object. See polarisSpectrum().
  stdev.feature.full <- XwY.svd$d

  ## ---- per-chromosome feature embedding -----------------------------------
  skymap.feature.chr <- .feature_chr_embedding(
    inf.emb.X = inf.emb.X, inf.emb.Y = inf.emb.Y,
    gene.names = colnames(X), feature.names = colnames(Y),
    gene.chr.ref = gene.chr.ref, chrs = chrs, r5 = r5,
    mc.cores = mc.cores, verbose = verbose)

  .new_skymap(list(
    skymap.cell        = skymap.cell,
    skymap.U           = skymap.U,
    skymap.V           = skymap.V,
    stdev.cell         = stdev.cell,
    skymap.P           = skymap.P,
    skymap.Q           = skymap.Q,
    stdev.feature      = stdev.feature,
    stdev.feature.full = stdev.feature.full,
    skymap.feature.chr = skymap.feature.chr,
    X.svd              = X.svd,
    Y.svd              = Y.svd,
    ## Records the ranks that were actually used, not only the thresholds. The
    ## POLARIS manuscript tabulates the selected r and r' per dataset, and the
    ## earlier version captured this list before any default was resolved, so
    ## the saved provenance did not contain them.
    input.param        = list(T1 = T1, T2 = T2, T3 = T3,
                              r1 = r1, r2 = r2, r3 = r3, r4 = r4, r5 = r5,
                              r1.auto = r1.auto, r2.auto = r2.auto,
                              x.sds = x.sds, y.sds = y.sds, seed = seed,
                              n.cells = nrow(skymap.cell),
                              n.features.X = ncol(X), n.features.Y = ncol(Y),
                              polaris.version = as.character(
                                utils::packageVersion("POLARIS")),
                              ## Records the software that produced this fit, so
                              ## a saved skymap stays self-documenting even if
                              ## the environment later changes. Surfaced by
                              ## polarisFitSummary() and polarisVersions().
                              versions = .capture_versions())
  ))
}


#' Per-chromosome decomposition of the cross-modal loading product
#' @noRd
.feature_chr_embedding <- function(inf.emb.X, inf.emb.Y,
                                   gene.names, feature.names,
                                   gene.chr.ref, chrs, r5,
                                   mc.cores, verbose) {
  say <- function(...) if (isTRUE(verbose)) message(...)

  if (is.null(gene.chr.ref)) {
    say(paste0("gene.chr.ref not supplied, so the per-chromosome feature ",
               "embedding is skipped. skymap.P and skymap.Q are still returned. ",
               "SkymapLinkageTable() needs the per-chromosome embedding; supply ",
               "gene.chr.ref to compute it."))
    return(NULL)
  }
  if (!methods::is(gene.chr.ref, "GRanges"))
    stop("`gene.chr.ref` must be a GRanges. See polarisTSSRef().", call. = FALSE)
  if (!"gene_name" %in% names(S4Vectors::mcols(gene.chr.ref)))
    stop("`gene.chr.ref` must have a `gene_name` metadata column.", call. = FALSE)

  ref.chr  <- as.character(GenomicRanges::seqnames(gene.chr.ref))
  ref.gene <- gene.chr.ref$gene_name
  feat.chr <- .feature_chr(feature.names)

  n.drop <- sum(is.na(feat.chr))
  if (n.drop)
    say(sprintf(paste0("%d of %d features of Y are not on the primary assembly ",
                       "(chrM, scaffolds or patch contigs) and are excluded from ",
                       "the per-chromosome embedding."),
                n.drop, length(feat.chr)))

  if (is.null(chrs)) {
    present <- intersect(unique(stats::na.omit(feat.chr)), unique(ref.chr))
    chrs <- .default_chrs()[.default_chrs() %in% present]
    if (!length(chrs)) chrs <- present
  }
  if (!length(chrs))
    stop(paste0("No chromosome is shared between the feature names of Y and ",
                "gene.chr.ref. Check that peaks are named chr-start-end and that ",
                "gene.chr.ref uses the same chromosome naming."), call. = FALSE)

  say(sprintf("Computing the feature embedding for %d chromosomes.", length(chrs)))

  res <- .maybe_mclapply(chrs, function(chr.idx) {
    gi <- which(gene.names %in% ref.gene[ref.chr == chr.idx])
    pj <- which(!is.na(feat.chr) & feat.chr == chr.idx)
    k <- min(r5, length(gi), length(pj))
    if (k < 3L) return(NULL)
    A <- matrixMultiplyEigen(inf.emb.X[gi, , drop = FALSE],
                             t(inf.emb.Y[pj, , drop = FALSE]))
    s <- suppressWarnings(RSpectra::svds(A, k = k))
    P <- s$u[, seq_len(k), drop = FALSE]
    Q <- s$v[, seq_len(k), drop = FALSE]
    rownames(P) <- gene.names[gi]
    rownames(Q) <- feature.names[pj]
    list(P = P, Q = Q, stdev = s$d[seq_len(k)])
  }, mc.cores = mc.cores, what = "chromosomes")

  names(res) <- chrs
  empty <- vapply(res, is.null, logical(1))
  if (all(empty)) {
    say(paste0("No chromosome had at least 3 genes and 3 features in common. ",
               "This is expected when the features of Y carry no genomic ",
               "coordinates, as for CITE-seq; use skymap.P and skymap.Q."))
  } else if (any(empty)) {
    say(sprintf("Skipped %s (fewer than 3 genes or features).",
                paste(chrs[empty], collapse = ", ")))
  }
  res
}
