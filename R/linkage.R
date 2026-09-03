#' Cis gene-peak linkage table with per-link standard errors
#'
#' Enumerates every within-chromosome gene-peak pair whose peak midpoint lies
#' within `window` of the gene's strand-aware transcription start site, scores
#' each pair by the singular-value-weighted inner product of its loadings, and
#' attaches a closed-form leave-one-cell-out standard error.
#'
#' @param skymap A [skymap] from [Polaris()], fitted with `gene.chr.ref` so that
#'   it carries the per-chromosome feature embedding.
#' @param gene.anno A `GRanges` with a `gene_name` metadata column and strand,
#'   used to place transcription start sites. See [polarisTSSRef()]. Must be the
#'   annotation for the genome build the peaks were called against.
#' @param window Half-width in bp of the cis window around the TSS. Peaks whose
#'   **midpoint** falls within `window` of the TSS are candidates. Default
#'   `1e6`, the value used throughout the POLARIS manuscript.
#' @param se Compute per-link standard errors, and from them the
#'   empirical-Bayes shrunken score and local false sign rate. Default `TRUE`.
#'   Setting `FALSE` returns `score_raw` and `dist` only and is much faster.
#' @param X,Y Optional cells-by-features matrices. **Supplying these changes the
#'   estimator.** See Details. Leave them `NULL` to reproduce the published
#'   standard errors.
#' @param tau2 Empirical-Bayes prior variance. Estimated by moments from the
#'   data when `NULL`.
#' @param chrs Restrict to these chromosomes. Default all present.
#' @param mc.cores Cores for the per-gene standard-error computation.
#' @param verbose Report progress.
#'
#' @return A `data.frame` sorted by decreasing `score_adj` (or `score_raw` when
#'   `se = FALSE`), with columns `gene`, `peak`, `chr`, `dist`, `score_raw`, and
#'   when `se = TRUE` also `SE`, `score_adj`, `p_lfsr`, `q_lfsr`.
#'
#' @details
#' The raw linkage score is the singular-value-weighted inner product
#' \deqn{\mathrm{IPw}(g,p) = P_g^\top \mathrm{diag}(d) Q_p,}
#' computed within each chromosome, so no gene is ever scored against a peak on
#' another chromosome.
#'
#' The standard error is a fixed-embedding leave-one-cell-out jackknife that
#' does not refit POLARIS. Writing the score as a bilinear form
#' \eqn{\mathrm{IPw} = x^\top K y} in the cell kernel \eqn{K = WW^\top}, all `n`
#' omitted-cell scores follow in closed form from two matrix-vector products,
#' and \eqn{\mathrm{SE} = \sqrt{\tfrac{n-1}{n}\sum_i (s_{(-i)} - \bar{s})^2}}.
#'
#' `z_dev`, the genomic-distance-aware statistic that the POLARIS benchmarks
#' rank by, is not a column here. Compute it with [SkymapZdev()], which needs
#' the whole table in order to form the within-distance-bin baselines.
#'
#' @section Which vectors the standard error uses:
#' The POLARIS standard error is defined on the **rank-`r` denoised** gene and
#' peak vectors. That is what you get by leaving `X` and `Y` as `NULL`: they are
#' reconstructed from the SVDs stored in the fit, and this is the estimator
#' behind every published POLARIS result.
#'
#' If you pass `X` and `Y`, the standard error is instead computed on those raw
#' normalized matrices. That is a different estimator, it is not the one in the
#' manuscript, and it will change `SE`, `score_adj` and any `z_dev` derived from
#' them. The argument is retained for comparison only, and warns when used.
#'
#' @examples
#' \dontrun{
#' anno <- polarisTSSRef("hg38")
#' fit  <- Polaris(X, Y, gene.chr.ref = anno, mc.cores = 4)
#' tab  <- SkymapLinkageTable(fit, gene.anno = anno, mc.cores = 4)
#' tab  <- SkymapZdev(tab)
#' head(tab[order(-tab$z_dev), ])
#' }
#'
#' @seealso [SkymapZdev()], [polarisTSSRef()]
#' @export
SkymapLinkageTable <- function(skymap,
                               gene.anno,
                               window = 1e6,
                               se = TRUE,
                               X = NULL, Y = NULL,
                               tau2 = NULL,
                               chrs = NULL,
                               mc.cores = 1L,
                               verbose = TRUE) {

  say <- function(...) if (isTRUE(verbose)) message(...)
  .check_skymap(skymap)

  if (is.null(skymap$skymap.feature.chr))
    stop(paste0("This skymap has no per-chromosome feature embedding. Re-run ",
                "Polaris() with gene.chr.ref supplied."), call. = FALSE)

  if (missing(gene.anno) || is.null(gene.anno))
    stop(paste0("`gene.anno` is required: a GRanges with a gene_name column and ",
                "strand, for the same genome build the peaks were called on. ",
                "Use polarisTSSRef('hg38'), polarisTSSRef('mm10'), or ",
                "polarisMakeTSSRef() for another build."), call. = FALSE)
  if (is.character(gene.anno))
    stop(paste0("`gene.anno` must be a GRanges, not a file path. Read it first, ",
                "or use polarisTSSRef()."), call. = FALSE)
  if (!methods::is(gene.anno, "GRanges"))
    stop("`gene.anno` must be a GRanges.", call. = FALSE)
  if (!"gene_name" %in% names(S4Vectors::mcols(gene.anno)))
    stop("`gene.anno` must have a `gene_name` metadata column.", call. = FALSE)

  ## ---- strand-aware TSS ----------------------------------------------------
  ## resize(fix = "start") is interpreted relative to strand, so a minus-strand
  ## gene is anchored at its higher coordinate and start() then returns the true
  ## TSS. Requires real strand information: with strand "*" every gene is
  ## treated as plus.
  if (any(as.character(GenomicRanges::strand(gene.anno)) == "*"))
    warning(paste0("Some ranges in `gene.anno` are unstranded ('*'). Those genes ",
                   "will be treated as plus-strand, so their TSS may be the ",
                   "wrong end of the gene."), call. = FALSE)

  tss.gr <- GenomicRanges::resize(gene.anno, width = 1L, fix = "start")
  tss.df <- data.frame(gene = gene.anno$gene_name,
                       chr  = as.character(GenomicRanges::seqnames(gene.anno)),
                       tss  = GenomicRanges::start(tss.gr),
                       stringsAsFactors = FALSE)
  dup <- duplicated(tss.df$gene)
  if (any(dup))
    say(sprintf(paste0("%d gene symbols appear more than once in gene.anno; ",
                       "keeping the first occurrence of each."), sum(dup)))
  tss.df <- tss.df[!dup, ]
  rownames(tss.df) <- tss.df$gene

  feat <- skymap$skymap.feature.chr
  if (is.null(chrs)) chrs <- names(feat)
  chrs <- chrs[vapply(chrs, function(ch) !is.null(feat[[ch]]), logical(1))]
  if (!length(chrs))
    stop("No chromosome in this skymap has a feature embedding.", call. = FALSE)

  ## ---- enumerate cis pairs and score them ---------------------------------
  build.chr <- function(ch) {
    cc <- feat[[ch]]
    P <- as.matrix(cc$P); Q <- as.matrix(cc$Q); d <- cc$stdev
    genes <- rownames(P); peaks <- rownames(Q)
    pk <- .parse_peaks(peaks)
    g.tss <- tss.df[genes, "tss"]
    g.chr <- tss.df[genes, "chr"]
    keep.g <- which(!is.na(g.tss) & g.chr == ch)
    if (!length(keep.g)) return(NULL)
    Pd <- P * rep(d, each = nrow(P))            # scales column k by d[k]
    out <- vector("list", length(keep.g))
    for (ii in seq_along(keep.g)) {
      gi <- keep.g[ii]
      t0 <- g.tss[gi]
      cand <- which(abs(pk$mid - t0) <= window)
      if (!length(cand)) next
      sc <- as.vector(Q[cand, , drop = FALSE] %*% Pd[gi, ])
      out[[ii]] <- data.frame(gene = genes[gi], peak = peaks[cand], chr = ch,
                              dist = abs(pk$mid[cand] - t0), score_raw = sc,
                              stringsAsFactors = FALSE)
    }
    do.call(rbind, out)
  }

  say(sprintf("Enumerating cis pairs over %d chromosomes.", length(chrs)))
  tab <- do.call(rbind, .maybe_mclapply(chrs, build.chr, mc.cores, "chromosomes"))
  if (is.null(tab) || !nrow(tab))
    stop(paste0("No cis pair found. Check that gene.anno matches the genome the ",
                "peaks were called on, and that `window` is not too small."),
         call. = FALSE)
  rownames(tab) <- NULL
  say(sprintf("%s cis pairs over %s genes.",
              format(nrow(tab), big.mark = ","),
              format(length(unique(tab$gene)), big.mark = ",")))

  if (!se) return(tab[order(-tab$score_raw), ])

  ## ---- per-link jackknife standard error ----------------------------------
  if (is.null(skymap$skymap.cell))
    stop("skymap.cell is needed for the standard error.", call. = FALSE)
  C <- as.matrix(skymap$skymap.cell)
  cells <- rownames(C)
  if (is.null(cells) || !length(cells))
    stop(paste0("skymap.cell has no cell names, so cells cannot be matched. ",
                "Re-run Polaris() with named rows on X and Y."), call. = FALSE)

  genes.u <- unique(tab$gene)
  peaks.u <- unique(tab$peak)
  have.raw <- !is.null(X) && !is.null(Y)

  if (have.raw) {
    warning(paste0("X and Y were supplied, so the standard error is computed on ",
                   "the raw normalized matrices. This is NOT the estimator used ",
                   "in the POLARIS manuscript, which is defined on the rank-r ",
                   "denoised vectors. Omit X and Y to reproduce published ",
                   "values."), call. = FALSE)
    if (is.null(rownames(X)) || is.null(rownames(Y)))
      stop("X and Y must have cell barcodes as rownames.", call. = FALSE)
    keep <- intersect(intersect(cells, rownames(X)), rownames(Y))
    if (length(keep) < 0.5 * nrow(C))
      stop(sprintf(paste0("Only %d of %d skymap cells were found in X and Y; ",
                          "check that the barcodes match."),
                   length(keep), nrow(C)), call. = FALSE)
    cells <- keep
    C <- C[cells, , drop = FALSE]
    X <- X[cells, genes.u, drop = FALSE]
    Y <- Y[cells, peaks.u, drop = FALSE]
  } else {
    if (is.null(skymap$X.svd) || is.null(skymap$Y.svd))
      stop(paste0("This skymap has no stored SVDs, so the denoised vectors ",
                  "cannot be reconstructed. Re-fit with Polaris(), or pass X ",
                  "and Y (a different estimator)."), call. = FALSE)
    if (is.null(rownames(skymap$skymap.P)) || is.null(rownames(skymap$skymap.Q)))
      stop("skymap.P/Q need rownames to map features onto the stored SVDs.",
           call. = FALSE)
    Xs <- skymap$X.svd; Ys <- skymap$Y.svd
    if (nrow(Xs$u) != nrow(C) || nrow(Ys$u) != nrow(C))
      stop("The stored SVDs do not match skymap.cell in cell count.", call. = FALSE)
    gi <- match(genes.u, rownames(skymap$skymap.P))
    pj <- match(peaks.u, rownames(skymap$skymap.Q))
    genes.u <- genes.u[!is.na(gi)]; gi <- gi[!is.na(gi)]
    peaks.u <- peaks.u[!is.na(pj)]; pj <- pj[!is.na(pj)]
    X <- Xs$u %*% (Xs$d * t(Xs$v[gi, , drop = FALSE]))
    Y <- Ys$u %*% (Ys$d * t(Ys$v[pj, , drop = FALSE]))
    dimnames(X) <- list(cells, genes.u)
    dimnames(Y) <- list(cells, peaks.u)
  }

  n <- length(cells)
  infX <- as.matrix(Matrix::crossprod(X, C))   # genes x D
  infY <- as.matrix(Matrix::crossprod(Y, C))   # peaks x D
  mcell <- rowSums(C^2)                        # diag(K)

  se.gene <- function(g, pk) {
    u  <- as.numeric(X[, g])
    tg <- as.numeric(C %*% infX[g, ])          # K x
    v  <- tg - u * mcell                       # K x - x * diag(K)
    s  <- C %*% t(infY[pk, , drop = FALSE])    # K y,  n x |pk|
    Yp <- as.matrix(Y[, pk, drop = FALSE])
    w1 <- u^2; w2 <- u * v; w3 <- v^2
    sumD  <- as.numeric(crossprod(s, u)) + as.numeric(crossprod(Yp, v))
    sumD2 <- as.numeric(crossprod(s^2, w1)) +
      2 * as.numeric(crossprod(s * Yp, w2)) +
      as.numeric(crossprod(Yp^2, w3))
    sqrt(pmax((n - 1) / n * (sumD2 - sumD^2 / n), 0))
  }

  say(sprintf("Computing leave-one-cell-out standard errors for %s pairs.",
              format(nrow(tab), big.mark = ",")))
  scorable <- tab$gene %in% genes.u & tab$peak %in% peaks.u
  sp <- split(which(scorable), tab$gene[scorable])
  se.list <- .maybe_mclapply(names(sp), function(g) {
    idx <- sp[[g]]
    data.frame(idx = idx, SE = se.gene(g, tab$peak[idx]))
  }, mc.cores, "gene standard-error chunks")
  se.df <- do.call(rbind, se.list)

  tab$SE <- NA_real_
  tab$SE[se.df$idx] <- se.df$SE
  ## The pre-package version could leave silent NAs here when an mclapply worker
  ## was killed, and order() then demoted those genes to the bottom of the table
  ## without any diagnostic. .maybe_mclapply now raises on worker failure; this
  ## checks the remaining path, features absent from the stored SVDs.
  if (anyNA(tab$SE))
    say(sprintf(paste0("%d of %d pairs have no standard error because the ",
                       "feature is absent from the stored SVDs; their score_adj ",
                       "and q_lfsr are NA."),
                sum(is.na(tab$SE)), nrow(tab)))

  ## ---- empirical-Bayes shrinkage and local false sign rate ----------------
  if (is.null(tau2))
    tau2 <- max(0, stats::var(tab$score_raw, na.rm = TRUE) -
                  mean(tab$SE^2, na.rm = TRUE))
  if (tau2 <= 0)
    warning(paste0("The empirical-Bayes prior variance estimated as 0, meaning ",
                   "the average squared standard error exceeds the variance of ",
                   "the scores. score_adj will be all zero and q_lfsr all 1. ",
                   "Rank by score_raw or by SkymapZdev() instead."),
            call. = FALSE)

  w <- tau2 / (tau2 + tab$SE^2)
  tab$score_adj <- w * tab$score_raw
  postsd <- sqrt(w) * tab$SE
  ## postsd == 0 is set to p = 1 rather than NA. This is deliberately the
  ## published behaviour: switching it to NA would change the number of
  ## non-missing p values and therefore every BH-adjusted q_lfsr in the table.
  ## It is conservative (such a pair can never be called) and only arises in
  ## degenerate cases, tau2 == 0 or an exactly-zero SE.
  tab$p_lfsr <- ifelse(!is.na(postsd) & postsd > 0,
                       stats::pnorm(-abs(tab$score_adj) / postsd), 1)
  tab$q_lfsr <- stats::p.adjust(tab$p_lfsr, method = "BH")
  attr(tab, "tau2") <- tau2

  tab[order(-tab$score_adj), ]
}
