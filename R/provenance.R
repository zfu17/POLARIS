## Provenance and audit helpers.
##
## These exist because the POLARIS manuscript needs three things that the
## pre-package code could not supply:
##   * Supplementary Table 2 wants a version for every tool, and a whole-archive
##     search returned zero version strings.
##   * Supplementary Table 3 wants the selected r and r' per dataset, but the old
##     `input.param` was captured before any default resolved, so it recorded
##     only the thresholds and left the ranks NULL.
##   * Supplementary Fig. 1 shows the denoiser spectrum, but the fit discarded
##     everything past the selected rank, so the figure could not be redrawn.

#' Versions of POLARIS and its dependencies
#'
#' Reports the version of POLARIS, of R, and of the packages POLARIS depends on,
#' in a form that can be pasted straight into a manuscript software table.
#'
#' @param which `"imports"` (default) reports POLARIS, R and everything in the
#'   `Imports` field. `"all"` adds the `Suggests` packages that are installed.
#' @param installed.only Drop rows for packages that are not installed. Default
#'   `TRUE`.
#'
#' @return A `data.frame` with columns `package` and `version`, POLARIS first,
#'   then R, then dependencies alphabetically.
#'
#' @details
#' Every fit also records this at the time it ran, in
#' `skymap$input.param$versions`, so a saved `skymap` documents the software that
#' produced it even if the environment later changes. [polarisFitSummary()]
#' surfaces that per fit.
#'
#' @examples
#' polarisVersions()
#'
#' @seealso [polarisFitSummary()]
#' @export
polarisVersions <- function(which = c("imports", "all"), installed.only = TRUE) {
  which <- match.arg(which)
  fields <- if (identical(which, "all")) c("Imports", "Suggests") else "Imports"
  pkgs <- unlist(lapply(fields, function(f) {
    d <- utils::packageDescription("POLARIS", fields = f)
    if (is.na(d) || !nzchar(d)) return(character())
    p <- trimws(strsplit(d, ",")[[1]])
    p <- sub("\\s*\\(.*\\)$", "", p)          # drop version constraints
    p[nzchar(p)]
  }), use.names = FALSE)
  ## base packages carry R's own version, so reporting them separately is noise
  pkgs <- setdiff(unique(pkgs), c("methods", "stats", "utils", "parallel",
                                  "graphics", "grDevices", "tools", "grid"))

  ver <- vapply(pkgs, function(p)
    tryCatch(as.character(utils::packageVersion(p)),
             error = function(e) NA_character_), character(1))

  out <- rbind(
    data.frame(package = "POLARIS",
               version = as.character(utils::packageVersion("POLARIS")),
               stringsAsFactors = FALSE),
    data.frame(package = "R",
               version = paste(R.version$major, R.version$minor, sep = "."),
               stringsAsFactors = FALSE),
    data.frame(package = pkgs[order(tolower(pkgs))],
               version = ver[order(tolower(pkgs))],
               stringsAsFactors = FALSE))
  rownames(out) <- NULL
  if (isTRUE(installed.only)) out <- out[!is.na(out$version), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @noRd
.capture_versions <- function() {
  v <- tryCatch(polarisVersions("imports"), error = function(e) NULL)
  if (is.null(v)) return(NULL)
  stats::setNames(as.list(v$version), v$package)
}

#' Tabulate the settings and selected ranks of one or more fits
#'
#' One row per fit, giving the data dimensions, the thresholds supplied and the
#' ranks the denoiser actually selected. This is the table the POLARIS
#' manuscript reports per dataset.
#'
#' @param ... [skymap] objects, or a single named `list` of them. Names become
#'   the `dataset` column; unnamed fits are numbered.
#'
#' @return A `data.frame`, one row per fit: `dataset`, `n.cells`,
#'   `n.features.X`, `n.features.Y`, `T1`, `T2`, `T3`, `r1`, `r2`, `r3`, `r4`,
#'   `r5`, `x.sds`, `y.sds`, `reconcile.ranks`, `seed`, `polaris.version`.
#'
#' @details
#' `r1` and `r2` are the modality ranks after being set to their common minimum;
#' `r1.auto` and `r2.auto`, the values the ratio rule chose before that, are kept
#' in `input.param` but not tabulated here.
#'
#' Fits produced by the pre-package script will show `NA` for the ranks, because
#' that version captured `input.param` before the defaults were resolved. Recover
#' those with `ncol(fit$skymap.U)` and `ncol(fit$skymap.P)`.
#'
#' @examples
#' set.seed(1)
#' a <- Polaris(polarisSimulate(n = 150, g = 50, p = 60, r = 4)$X,
#'              polarisSimulate(n = 150, g = 50, p = 60, r = 4)$Y,
#'              x.sds = 0, verbose = FALSE)
#' polarisFitSummary(pbmc = a)
#'
#' @seealso [polarisVersions()], [polarisSpectrum()]
#' @export
polarisFitSummary <- function(...) {
  fits <- list(...)
  if (length(fits) == 1L && is.list(fits[[1]]) && !inherits(fits[[1]], "skymap"))
    fits <- fits[[1]]
  if (!length(fits)) stop("Supply at least one skymap.", call. = FALSE)
  nm <- names(fits)
  if (is.null(nm)) nm <- rep("", length(fits))
  nm[!nzchar(nm)] <- paste0("fit", seq_along(fits))[!nzchar(nm)]

  grab <- function(ip, f, default = NA) {
    v <- ip[[f]]
    if (is.null(v) || !length(v)) default else v[[1]]
  }

  do.call(rbind, lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    .check_skymap(fit)
    ip <- fit$input.param
    if (is.null(ip)) ip <- list()
    data.frame(
      dataset      = nm[i],
      n.cells      = nrow(fit$skymap.cell),
      n.features.X = grab(ip, "n.features.X", if (!is.null(fit$skymap.P)) nrow(fit$skymap.P) else NA),
      n.features.Y = grab(ip, "n.features.Y", if (!is.null(fit$skymap.Q)) nrow(fit$skymap.Q) else NA),
      T1 = grab(ip, "T1"), T2 = grab(ip, "T2"), T3 = grab(ip, "T3"),
      r1 = grab(ip, "r1"), r2 = grab(ip, "r2"),
      ## fall back to the embedding widths, which is how these have to be
      ## recovered from a pre-package fit
      r3 = grab(ip, "r3", ncol(fit$skymap.U)),
      r4 = grab(ip, "r4", ncol(fit$skymap.V)),
      r5 = grab(ip, "r5", if (!is.null(fit$skymap.P)) ncol(fit$skymap.P) else NA),
      x.sds = grab(ip, "x.sds"), y.sds = grab(ip, "y.sds"),
      ## Which rank convention the fit used. Decides which published analyses it
      ## reproduces, and is invisible from ncol(skymap.cell).
      reconcile.ranks = grab(ip, "reconcile.ranks", NA),
      seed  = grab(ip, "seed"),
      polaris.version = grab(ip, "polaris.version", NA_character_),
      stringsAsFactors = FALSE, row.names = NULL)
  }))
}

#' The denoiser spectrum and where the rank was cut
#'
#' Returns the singular values a fit retained, the consecutive ratios the
#' denoiser thresholds on, and which components were kept. Use it to audit rank
#' selection, or to draw the spectrum figure.
#'
#' @param skymap A [skymap] from [Polaris()].
#' @param modality `"X"`, `"Y"`, `"feature"`, or `"all"` (default) to stack all
#'   three.
#'
#' @return A `data.frame` with `modality`, `component`, `singular_value`,
#'   `ratio` (\eqn{\sigma_j / \sigma_{j+1}}, `NA` for the last component),
#'   `threshold`, and `retained`.
#'
#' @details
#' The rank rule is
#' \deqn{k = \max\{j : 1 \le j \le K-1,\ \sigma_j / \sigma_{j+1} > T\},}
#' where \eqn{K} is the number of components **computed**, not the number of
#' cells. That bound matters: `Polaris()` computes at most
#' `min(50, min(dim) - 1)` components, so the selected rank can never exceed
#' that cap however slowly the spectrum decays. It is also why a threshold of
#' `T = 1.01` can be reported as selecting anywhere from ten to thirty-odd
#' components on different datasets — the answer depends on the data and on the
#' cap, not on the threshold alone.
#'
#' For `"feature"` the full spectrum is available only for fits made by this
#' package, which stores it as `stdev.feature.full`. Earlier fits kept only the
#' retained values, so `ratio` runs out at the selected rank.
#'
#' @examples
#' set.seed(1)
#' sim <- polarisSimulate(n = 200, g = 60, p = 80, r = 5)
#' fit <- Polaris(sim$X, sim$Y, x.sds = 0, verbose = FALSE)
#' head(polarisSpectrum(fit, "X"))
#'
#' @seealso [polarisFitSummary()]
#' @export
polarisSpectrum <- function(skymap, modality = c("all", "X", "Y", "feature")) {
  .check_skymap(skymap)
  modality <- match.arg(modality)
  ip <- skymap$input.param
  if (is.null(ip)) ip <- list()

  one <- function(d, k, thr, lab) {
    if (is.null(d) || !length(d)) return(NULL)
    d <- as.numeric(d)
    ratio <- c(d[-length(d)] / d[-1L], NA_real_)
    data.frame(modality = lab,
               component = seq_along(d),
               singular_value = d,
               ratio = ratio,
               threshold = if (is.null(thr) || !length(thr)) NA_real_ else thr[[1]],
               retained = if (is.null(k) || !length(k)) NA else seq_along(d) <= k[[1]],
               stringsAsFactors = FALSE, row.names = NULL)
  }

  parts <- list()
  if (modality %in% c("all", "X"))
    parts$X <- one(skymap$X.svd$d, ip$r1, ip$T1, "X")
  if (modality %in% c("all", "Y"))
    parts$Y <- one(skymap$Y.svd$d, ip$r2, ip$T2, "Y")
  if (modality %in% c("all", "feature")) {
    fd <- if (!is.null(skymap$stdev.feature.full)) skymap$stdev.feature.full
          else skymap$stdev.feature
    parts$feature <- one(fd, ip$r5, ip$T3, "feature")
  }
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts))
    stop("This skymap carries no spectrum for the requested modality.",
         call. = FALSE)
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}
