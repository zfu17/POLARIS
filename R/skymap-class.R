#' The skymap object
#'
#' A `skymap` is the fitted POLARIS model returned by [Polaris()]. It is a plain
#' named list carrying an S3 class attribute, so every element remains reachable
#' with `$` exactly as before this code became a package. Existing analysis
#' scripts that index a POLARIS fit by name continue to work unchanged; the class
#' only adds `print` and `summary` methods and lets other functions check that
#' they were handed a real fit.
#'
#' @section Elements:
#' \describe{
#'   \item{`skymap.cell`}{`n` x `(r3 + r4)` matrix. The joint cell embedding
#'     \eqn{W = [\hat{U} | \hat{V}]}. Rows are cells, named by barcode.}
#'   \item{`skymap.U`, `skymap.V`}{`n` x `r3` and `n` x `r4` matrices. The paired
#'     cell embeddings. Their columns are matched one to one through the SVD,
#'     which is what makes the per-cell concordance score meaningful.}
#'   \item{`stdev.cell`}{singular values of the cell embedding.}
#'   \item{`skymap.P`, `skymap.Q`}{`g` x `r5` and `p` x `r5` matrices. Genome-wide
#'     gene and second-modality feature loadings.}
#'   \item{`stdev.feature`}{singular-value weights \eqn{d} of the feature
#'     embedding. These enter the linkage score.}
#'   \item{`skymap.feature.chr`}{named list, one entry per chromosome, each
#'     `list(P, Q, stdev)`. This is the per-chromosome decomposition that
#'     [SkymapLinkageTable()] scores, so every linkage is within-chromosome.}
#'   \item{`X.svd`, `Y.svd`}{the truncated SVDs of the two denoised modalities.
#'     Retained because [SkymapLinkageTable()] reconstructs the denoised feature
#'     vectors from them to compute the jackknife standard error.}
#'   \item{`input.param`}{the thresholds supplied and the ranks actually
#'     selected. Unlike earlier versions this records the resolved ranks, not
#'     just the thresholds.}
#' }
#'
#' Optional elements are added in place by downstream functions:
#' `concordance` and `simscore.*` by [SkymapSimScore()], `umap.*` by
#' [SkymapUMAP()], and `nb.*` by [SkymapFindNeighbors()].
#'
#' @name skymap
#' @aliases skymap-class
NULL

#' @noRd
.new_skymap <- function(x) {
  structure(x, class = c("skymap", "list"))
}

#' @noRd
.check_skymap <- function(skymap, need = character()) {
  if (!inherits(skymap, "skymap")) {
    ## Accept a bare list so that fits saved by the pre-package scripts still work.
    if (!is.list(skymap) || is.null(skymap$skymap.cell))
      stop("`skymap` must be a skymap object returned by Polaris().", call. = FALSE)
  }
  missing <- need[!vapply(need, function(n) !is.null(skymap[[n]]), logical(1))]
  if (length(missing))
    stop(sprintf("This skymap has no %s. %s",
                 paste(missing, collapse = ", "),
                 "Re-run Polaris(), or the function that produces it."),
         call. = FALSE)
  invisible(TRUE)
}

#' @export
print.skymap <- function(x, ...) {
  n <- nrow(x$skymap.cell)
  cat(sprintf("A skymap: POLARIS fit on %s cells\n", format(n, big.mark = ",")))
  cat(sprintf("  cell embedding   W  %s x %d   (U: %d, V: %d)\n",
              format(n, big.mark = ","), ncol(x$skymap.cell),
              ncol(x$skymap.U), ncol(x$skymap.V)))
  if (!is.null(x$skymap.P))
    cat(sprintf("  feature loadings P  %s x %d   Q  %s x %d\n",
                format(nrow(x$skymap.P), big.mark = ","), ncol(x$skymap.P),
                format(nrow(x$skymap.Q), big.mark = ","), ncol(x$skymap.Q)))
  if (!is.null(x$skymap.feature.chr)) {
    ok <- sum(!vapply(x$skymap.feature.chr, is.null, logical(1)))
    cat(sprintf("  per-chromosome   %d of %d chromosomes decomposed\n",
                ok, length(x$skymap.feature.chr)))
  }
  ip <- x$input.param
  if (!is.null(ip))
    cat(sprintf("  thresholds       T1=%g T2=%g T3=%g\n",
                ip$T1, ip$T2, ip$T3))
  extra <- setdiff(names(x), c("skymap.cell", "skymap.U", "skymap.V", "stdev.cell",
                               "skymap.P", "skymap.Q", "stdev.feature",
                               "skymap.feature.chr", "X.svd", "Y.svd",
                               "input.param"))
  if (length(extra))
    cat(sprintf("  also present     %s\n", paste(extra, collapse = ", ")))
  invisible(x)
}

#' @export
summary.skymap <- function(object, ...) {
  print(object)
  if (!is.null(object$concordance)) {
    q <- stats::quantile(object$concordance, c(0, .25, .5, .75, 1), na.rm = TRUE)
    cat("\nconcordance (rescaled to [0,1]):\n")
    print(round(q, 4))
  }
  invisible(object)
}
