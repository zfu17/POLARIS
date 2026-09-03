#' Simulate a paired two-modality dataset
#'
#' Generates two matrices on the same cells that share a low-rank latent
#' structure, with a per-cell weight controlling how strongly each cell's
#' second modality reflects the shared factors. Cells with a low weight are
#' cross-modality discordant by construction, so the simulation gives a ground
#' truth for the concordance score.
#'
#' Intended for examples, tests and teaching rather than for benchmarking.
#'
#' @param n Number of cells.
#' @param g Number of features in modality 1.
#' @param p Number of features in modality 2.
#' @param r Number of shared latent factors.
#' @param noise Standard deviation of the independent noise added to each
#'   modality.
#' @param discordant Fraction of cells whose second modality is driven by
#'   private rather than shared factors. These are the low-concordance cells.
#' @param feature.names Whether to name modality-2 features as genomic
#'   intervals (`chr1-100-600`, ...) so the simulation can be passed through
#'   [SkymapLinkageTable()]. Default `TRUE`.
#'
#' @return A list with
#'   \describe{
#'     \item{`X`, `Y`}{cells-by-features matrices, centered and scaled.}
#'     \item{`Z`}{the shared latent factors, `n` x `r`.}
#'     \item{`shared`}{per-cell weight on the shared factors, in `[0,1]`. Cells
#'       with a low value should receive a low concordance score.}
#'     \item{`anno`}{a `GRanges` placing the modality-1 features on `chr1`, so
#'       the simulation also exercises the linkage path.}
#'   }
#'
#' @examples
#' set.seed(1)
#' sim <- polarisSimulate(n = 200, g = 60, p = 80, r = 5)
#' fit <- SkymapSimScore(Polaris(sim$X, sim$Y, x.sds = 0, verbose = FALSE))
#' # discordant cells should score lower
#' tapply(fit$concordance, sim$shared < 0.5, mean)
#'
#' @export
polarisSimulate <- function(n = 300, g = 80, p = 120, r = 5,
                            noise = 0.5, discordant = 0.2,
                            feature.names = TRUE) {

  if (r >= min(g, p)) stop("`r` must be smaller than both `g` and `p`.", call. = FALSE)
  if (discordant < 0 || discordant > 1)
    stop("`discordant` must be between 0 and 1.", call. = FALSE)

  ## Shared latent factors, with decaying strength so the spectrum has a gap for
  ## the rank-selection rule to find.
  Z <- matrix(stats::rnorm(n * r), n, r)
  Z <- Z * rep(sqrt(seq(r, 1)), each = n)

  ## Private factors for modality 2, used by the discordant cells.
  Zp <- matrix(stats::rnorm(n * r), n, r) * rep(sqrt(seq(r, 1)), each = n)

  shared <- rep(1, n)
  n.dis <- round(discordant * n)
  if (n.dis > 0) shared[sample.int(n, n.dis)] <- stats::runif(n.dis, 0, 0.4)

  Bx <- matrix(stats::rnorm(r * g), r, g)
  By <- matrix(stats::rnorm(r * p), r, p)

  X <- Z %*% Bx + matrix(stats::rnorm(n * g, sd = noise), n, g)
  Zy <- shared * Z + sqrt(pmax(1 - shared^2, 0)) * Zp
  Y <- Zy %*% By + matrix(stats::rnorm(n * p, sd = noise), n, p)

  cells <- sprintf("cell_%04d", seq_len(n))
  genes <- sprintf("Gene%03d", seq_len(g))
  if (isTRUE(feature.names)) {
    st <- seq(1000, by = 5000, length.out = p)
    feats <- sprintf("chr1-%d-%d", st, st + 500)
  } else {
    feats <- sprintf("Feature%03d", seq_len(p))
  }

  X <- scale(X); Y <- scale(Y)
  dimnames(X) <- list(cells, genes)
  dimnames(Y) <- list(cells, feats)

  ## Place the modality-1 features on chr1, interleaved with the modality-2
  ## intervals, alternating strand so the strand-aware TSS logic is exercised.
  gs <- seq(1000, by = 5000 * max(1, floor(p / g)), length.out = g)
  anno <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = gs, end = gs + 2000),
    strand   = rep(c("+", "-"), length.out = g),
    gene_name = genes,
    gene_type = "protein_coding")

  list(X = X, Y = Y, Z = Z, shared = shared, anno = anno)
}
