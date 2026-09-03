#' Genomic-distance breakpoints used by POLARIS
#'
#' The fixed distance strata, in bp, within which linkage scores are
#' standardized by [SkymapZdev()]. These are the breakpoints used throughout the
#' POLARIS manuscript.
#'
#' @format A numeric vector of 9 breakpoints, from 0 to `Inf`.
#' @export
POLARIS_DIST_BREAKS <- c(0, 1e4, 2.5e4, 5e4, 1e5, 2.5e5, 5e5, 1e6, Inf)

#' Minimum peak-to-TSS distance for distance-stratified analyses
#'
#' Pairs closer than this were excluded from every distance-stratified analysis
#' in the POLARIS manuscript, because the Hi-C loop callers used as ground truth
#' have an effective resolution floor near this distance.
#'
#' @format A single number, 10000.
#' @export
POLARIS_MIN_DIST <- 1e4

#' Distance-aware z-statistic for gene-peak linkage
#'
#' Standardizes each raw linkage score within its genomic-distance stratum,
#' giving the statistic the POLARIS gene-peak benchmarks rank by. Raw scores
#' decay with genomic distance, so ranking by magnitude recovers proximity
#' rather than regulatory signal.
#'
#' @param tab A linkage table from [SkymapLinkageTable()], with `dist`,
#'   `score_raw` and `SE` columns.
#' @param breaks Distance breakpoints in bp. Default [POLARIS_DIST_BREAKS].
#' @param floor.se Also return `z_dev_fl`, a variant with `SE` floored at its
#'   `floor.quantile` quantile. This tames inflation from very small standard
#'   errors and was used for the threshold-level permutation FDR. Default `TRUE`.
#' @param floor.quantile Quantile at which to floor `SE`. Default `0.10`.
#'
#' @return `tab` with added columns:
#'   \describe{
#'     \item{`distbin`}{the distance stratum, a factor.}
#'     \item{`mu0`}{the mean `score_raw` of that stratum.}
#'     \item{`z_dev`}{the distance-aware z-statistic.}
#'     \item{`z_dev_fl`}{the SE-floored variant, when `floor.se = TRUE`.}
#'   }
#'   Rows are not removed and the row order is unchanged.
#'
#' @details
#' \deqn{z_{\mathrm{dev}}(g,p) = \frac{\mathrm{IPw}(g,p) -
#'   \mu_0[\mathrm{bin}(\delta)]}{\mathrm{SE}(g,p)}}
#' where \eqn{\mu_0} is the mean raw score in the pair's distance stratum and
#' \eqn{\delta} is the peak-midpoint-to-TSS distance.
#'
#' `z_dev` is a variance-standardized **ranking** statistic. It is not a
#' calibrated probability or a false-discovery rate: the standard error is a
#' fixed-embedding jackknife conditional on the fitted model and computed on the
#' denoised modalities.
#'
#' Pairs whose `SE` is missing or not positive get `z_dev = NA`; `mu0` is still
#' computed from all pairs in the stratum.
#'
#' @section Excluding proximal pairs:
#' Distance-stratified analyses in the POLARIS manuscript exclude pairs closer
#' than [POLARIS_MIN_DIST] (10 kb). This function does **not** drop them, so
#' that it never silently discards rows. Filter explicitly:
#'
#' ```
#' tab <- SkymapZdev(tab)
#' tab <- tab[tab$dist >= POLARIS_MIN_DIST, ]
#' ```
#'
#' Because `mu0` is computed per stratum, filtering before or after makes no
#' difference to the `z_dev` of the retained pairs.
#'
#' @examples
#' tab <- data.frame(
#'   gene = rep(c("A", "B"), each = 5),
#'   peak = paste0("chr1-", seq(1, 91, by = 10), "-", seq(2, 92, by = 10)),
#'   chr = "chr1",
#'   dist = c(5e3, 2e4, 6e4, 3e5, 8e5, 8e3, 3e4, 7e4, 4e5, 9e5),
#'   score_raw = c(5, 4, 3, 2, 1, 4.5, 3.5, 2.5, 1.5, 0.5),
#'   SE = rep(0.5, 10))
#' SkymapZdev(tab)[, c("dist", "distbin", "mu0", "z_dev")]
#'
#' @seealso [SkymapLinkageTable()]
#' @export
SkymapZdev <- function(tab,
                       breaks = POLARIS_DIST_BREAKS,
                       floor.se = TRUE,
                       floor.quantile = 0.10) {

  if (!is.data.frame(tab))
    stop("`tab` must be a data.frame from SkymapLinkageTable().", call. = FALSE)
  need <- c("dist", "score_raw", "SE")
  miss <- setdiff(need, names(tab))
  if (length(miss))
    stop(sprintf(paste0("`tab` is missing the column(s) %s. Run ",
                        "SkymapLinkageTable() with se = TRUE."),
                 paste(miss, collapse = ", ")), call. = FALSE)
  if (max(tab$dist, na.rm = TRUE) > max(breaks))
    stop(sprintf(paste0("Distances up to %g exceed the largest breakpoint %g, ",
                        "so some pairs would fall outside every stratum."),
                 max(tab$dist, na.rm = TRUE), max(breaks)), call. = FALSE)

  tab$distbin <- cut(tab$dist, breaks, include.lowest = TRUE, right = FALSE)

  ## mu0 is the stratum mean of the raw score, over all pairs in the stratum.
  mu0 <- tapply(tab$score_raw, tab$distbin, mean, na.rm = TRUE)
  tab$mu0 <- as.numeric(mu0[as.character(tab$distbin)])

  ok <- is.finite(tab$SE) & tab$SE > 0
  tab$z_dev <- NA_real_
  tab$z_dev[ok] <- (tab$score_raw[ok] - tab$mu0[ok]) / tab$SE[ok]

  if (isTRUE(floor.se)) {
    fl <- as.numeric(stats::quantile(tab$SE[ok], floor.quantile, na.rm = TRUE))
    tab$z_dev_fl <- NA_real_
    tab$z_dev_fl[ok] <- (tab$score_raw[ok] - tab$mu0[ok]) / pmax(tab$SE[ok], fl)
    attr(tab, "se.floor") <- fl
  }
  tab
}
