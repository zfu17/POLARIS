#' Per-cell cross-modality concordance score
#'
#' Computes the agreement between a cell's coordinates in the two paired cell
#' embeddings \eqn{\hat{U}(r)} and \eqn{\hat{V}(r)}. Because their columns are
#' matched one to one through the SVD, a cell's coordinates in the two
#' embeddings are directly comparable.
#'
#' @param skymap A [skymap] from [Polaris()].
#' @param metric One or more of `"cosine"`, `"euclidean"`, `"inner.prod"`,
#'   `"cor"`. Defaults to `"cosine"` alone, which is the metric the POLARIS
#'   concordance score is defined on.
#' @param n.pc Number of paired components used. Defaults to all of them.
#'
#' @return The `skymap` with elements added in place:
#'   \describe{
#'     \item{`concordance`}{the rescaled score in \eqn{[0,1]}, present whenever
#'       `"cosine"` is requested. This is the quantity the POLARIS manuscript
#'       calls the concordance score and the one to use for display and
#'       thresholding.}
#'     \item{`simscore.cosine`}{the raw cosine similarity in \eqn{[-1,1]}.}
#'     \item{`simscore.euclidean`, `simscore.inner.prod`, `simscore.cor`}{the
#'       other metrics, when requested.}
#'   }
#'
#' @details
#' The concordance score of cell \eqn{i} is the cosine similarity of its paired
#' coordinates, rescaled to the unit interval:
#' \deqn{\mathrm{concordance}_i = \tfrac{1}{2}\left(\cos\left(\hat{U}(r)_{[i,]},
#'   \hat{V}(r)_{[i,]}\right) + 1\right) \in [0,1].}
#' Concordance near 1 indicates a cell whose local structure is similarly
#' supported across modalities; lower values flag uneven cross-modality support.
#'
#' Both forms are returned because they differ by a monotone transform and so
#' rank cells identically, but only the rescaled form is bounded in
#' \eqn{[0,1]}. Returning `concordance` directly removes the need to rescale by
#' hand at plotting time.
#'
#' All four metrics are invariant to the joint rotation that the POLARIS
#' optimizer is defined up to, so each is well defined.
#'
#' @examples
#' set.seed(1)
#' sim <- polarisSimulate(n = 200, g = 60, p = 80, r = 5)
#' fit <- Polaris(sim$X, sim$Y, x.sds = 0, verbose = FALSE)
#' fit <- SkymapSimScore(fit)
#' summary(fit$concordance)
#'
#' @export
SkymapSimScore <- function(skymap,
                           metric = c("cosine", "euclidean", "inner.prod", "cor"),
                           n.pc = NULL) {

  .check_skymap(skymap, need = c("skymap.U", "skymap.V"))

  ## match.arg with several.ok picks up the FULL default vector when the caller
  ## supplies nothing, which in the pre-package version meant a bare call
  ## silently computed all four metrics rather than the documented cosine. Take
  ## the first choice as the default instead.
  if (missing(metric)) metric <- "cosine"
  metric <- match.arg(metric, c("cosine", "euclidean", "inner.prod", "cor"),
                      several.ok = TRUE)

  U <- as.matrix(skymap$skymap.U)
  V <- as.matrix(skymap$skymap.V)

  ## Clamp against BOTH factors: r3 and r4 are separately settable in Polaris(),
  ## so ncol(U) and ncol(V) can differ and slicing V by ncol(U) would go out of
  ## bounds.
  avail <- min(ncol(U), ncol(V))
  if (is.null(n.pc)) n.pc <- avail
  if (n.pc > avail) {
    warning(sprintf("n.pc = %d exceeds the %d paired components available; using %d.",
                    n.pc, avail, avail), call. = FALSE)
    n.pc <- avail
  }
  if (n.pc < 1L) stop("n.pc must be at least 1.", call. = FALSE)
  ## drop = FALSE matters: with n.pc = 1 a plain [ , 1:1 ] collapses to a vector
  ## and rowSums() then fails.
  U <- U[, seq_len(n.pc), drop = FALSE]
  V <- V[, seq_len(n.pc), drop = FALSE]

  for (m in metric) {
    skymap[[paste0("simscore.", m)]] <- switch(
      m,
      cosine     = as.vector(rowSums(U * V) /
                               (sqrt(rowSums(U * U)) * sqrt(rowSums(V * V)))),
      euclidean  = as.vector(1 - sqrt(rowSums((U - V)^2))),
      inner.prod = as.vector(rowSums(U * V)),
      ## Needs at least 2 components, and needs U/V to be matrices so that
      ## U[i, ] collapses to a vector. As data.frames this returned a k x k
      ## all-NA matrix per cell.
      cor        = if (n.pc < 2L) {
        warning("metric = 'cor' needs at least 2 components; skipping.",
                call. = FALSE)
        NULL
      } else {
        vapply(seq_len(nrow(U)), function(i) stats::cor(U[i, ], V[i, ]),
               numeric(1))
      })
  }

  if ("cosine" %in% metric) {
    skymap$concordance <- (skymap$simscore.cosine + 1) / 2
    names(skymap$concordance) <- rownames(U)
  }

  skymap
}
