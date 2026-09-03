## Internal helpers. Nothing here is exported.

#' Default primary-assembly chromosome names
#'
#' The human primary assembly, in the order POLARIS has always used.
#' @noRd
.default_chrs <- function() paste0("chr", c(seq_len(22), "X", "Y"))

#' Extract a chromosome label from feature (peak) names
#'
#' Peak names are expected to encode the chromosome first, as `chr1-100-200`
#' or `chr1:100-200`. The pattern deliberately matches only the primary
#' assembly, so features on chrM, on unplaced scaffolds (`GL000194.1`) and on
#' patch contigs (`chr1_KI270706v1_random`) return `NA` and are excluded from
#' the per-chromosome feature embedding. That is the historical behaviour; what
#' is new is that [Polaris()] now reports how many were dropped instead of
#' doing it silently.
#'
#' The trailing word boundary is what keeps `chr2` from matching `chr20`: `\\d+`
#' is greedy and the boundary is satisfied by the following separator.
#' @noRd
.feature_chr <- function(x) {
  pat <- "\\bchr\\s*(?:\\d+|[XY])\\b"
  out <- rep(NA_character_, length(x))
  mm <- regexpr(pat, x, perl = TRUE, ignore.case = TRUE)
  ok <- mm != -1L
  out[ok] <- regmatches(x, mm)
  out
}

#' Parse peak names into chromosome, start, end and midpoint
#'
#' Accepts `chr-start-end` and `chr:start-end`. Unlike the original
#' `do.call(rbind, strsplit(peaks, '[:-]'))`, a name that does not split into
#' exactly three fields is reported rather than silently recycled into the
#' matrix.
#' @noRd
.parse_peaks <- function(peaks) {
  parts <- strsplit(peaks, "[:-]")
  nf <- lengths(parts)
  if (any(nf != 3L)) {
    bad <- utils::head(peaks[nf != 3L], 3L)
    stop(sprintf(paste0("%d feature name(s) do not parse as chr-start-end or ",
                        "chr:start-end, for example: %s. Rename the features ",
                        "of the second modality, or drop them before fitting."),
                 sum(nf != 3L), paste(bad, collapse = ", ")), call. = FALSE)
  }
  m <- matrix(unlist(parts, use.names = FALSE), ncol = 3L, byrow = TRUE)
  start <- suppressWarnings(as.numeric(m[, 2L]))
  end   <- suppressWarnings(as.numeric(m[, 3L]))
  if (anyNA(start) || anyNA(end)) {
    bad <- utils::head(peaks[is.na(start) | is.na(end)], 3L)
    stop(sprintf("Non-numeric coordinates in feature name(s): %s.",
                 paste(bad, collapse = ", ")), call. = FALSE)
  }
  data.frame(peak = peaks, chr = m[, 1L], start = start, end = end,
             mid = (start + end) / 2, stringsAsFactors = FALSE)
}

#' Rank from the singular-value ratio rule
#'
#' Implements the denoiser of the POLARIS Methods,
#' \eqn{k = \max\{j : \sigma_j / \sigma_{j+1} > T\}}.
#'
#' The original code wrote this as `which(...) %>% max()`, which returns `-Inf`
#' with a warning when no ratio clears `T`, and then failed several lines later
#' with `Error in 1:r : result would be too long a vector`. Here the empty case
#' is named.
#' @noRd
.rank_from_ratio <- function(d, threshold, what) {
  d <- d[is.finite(d) & d > 0]
  if (length(d) < 2L)
    stop(sprintf("Cannot select a rank for %s: fewer than two positive singular values.",
                 what), call. = FALSE)
  ratios <- d[-length(d)] / d[-1L]
  hit <- which(ratios > threshold)
  if (!length(hit)) {
    stop(sprintf(paste0("No spectral gap above T = %g in %s: the singular-value ",
                        "spectrum is flat, so the automatic denoiser cannot pick a ",
                        "rank. Lower the threshold, or set the rank directly."),
                 threshold, what), call. = FALSE)
  }
  max(hit)
}

#' Number of components a thin SVD solver can return
#'
#' Both `irlba::irlba` and `RSpectra::svds` require strictly fewer components
#' than `min(dim)`. The original `min(50, ncol(X))` equals `min(dim)` whenever a
#' modality has 50 or fewer features, which makes `irlba` hard-error. Capping at
#' `min(dim) - 1` keeps the historical value of 50 for every real dataset while
#' letting small inputs through.
#' @noRd
.svd_k <- function(x, cap = 50L) {
  lim <- min(dim(x)) - 1L
  if (lim < 2L)
    stop(sprintf(paste0("Input is too small to decompose: %d x %d. At least 3 ",
                        "cells and 3 features are needed."),
                 nrow(x), ncol(x)), call. = FALSE)
  min(cap, lim)
}

#' Stop if any element of an mclapply result is an error
#'
#' `parallel::mclapply` does not raise; a worker that fails leaves a `try-error`
#' object in the returned list. The original code assigned names over the result
#' and returned it, so a failed chromosome travelled downstream disguised as a
#' successful fit. This turns that into an error at the point of failure.
#' @noRd
.check_mclapply <- function(res, what) {
  bad <- vapply(res, function(z) inherits(z, "try-error"), logical(1))
  if (any(bad)) {
    msg <- conditionMessage(attr(res[[which(bad)[1L]]], "condition"))
    stop(sprintf(paste0("%d of %d %s failed. First error: %s\n",
                        "If this is a memory failure, lower mc.cores."),
                 sum(bad), length(res), what, msg), call. = FALSE)
  }
  res
}

#' lapply or mclapply depending on mc.cores, with error checking
#' @noRd
.maybe_mclapply <- function(X, FUN, mc.cores, what) {
  if (mc.cores > 1L && .Platform$OS.type != "windows") {
    .check_mclapply(parallel::mclapply(X, FUN, mc.cores = mc.cores), what)
  } else {
    lapply(X, FUN)
  }
}

#' Run a function under a fixed seed without disturbing the caller's RNG
#'
#' `irlba::irlba` draws a random starting vector, which has two consequences we
#' do not want in a package: consecutive fits of the same data differ (by about
#' 1e-9 in the singular values, so scientifically irrelevant but not bit
#' reproducible), and the caller's random stream is advanced as a side effect, so
#' anything the user seeded around a POLARIS fit stops being reproducible.
#'
#' This saves `.Random.seed`, sets a fixed seed, evaluates, and restores the
#' previous state on exit. `RSpectra::svds` needs none of this; it is already
#' deterministic.
#' @noRd
.deterministic <- function(seed, fun) {
  if (is.null(seed)) return(fun())
  has <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (has) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (has) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  fun()
}

#' Truncated SVD, choosing a solver that suits the shape
#'
#' Solver choice for the two denoising SVDs is measured, not assumed. On real
#' scMultiome input, `X` is effectively dense (scale.data) which suits irlba,
#' while `Y` is genuinely sparse TF-IDF where RSpectra wins; using each where it
#' wins takes the block down 10-13%. Both return the same subspace (minCos
#' 1.000000, leading singular value to 4e-16 at n = 35,000).
#'
#' Neither iterative solver is appropriate when the requested rank is a large
#' fraction of `min(dim)`: irlba then warns that it did not converge and its
#' trailing directions are unreliable, which matters because POLARIS uses all of
#' them. Below that ratio an exact full SVD is both faster and correct, so small
#' inputs take that path. `RSpectra::svds` already does this internally; irlba
#' does not.
#' @noRd
.thin_svd <- function(x, k, solver = c("irlba", "svds"), seed = NULL) {
  solver <- match.arg(solver)
  lim <- min(dim(x))
  ## Iterative methods are only worth it when k is a small fraction of min(dim).
  ## For the published datasets (k = 50 of 816 cells) this is comfortably true,
  ## so the historical solver is used and results are unchanged.
  if (k >= lim / 3) {
    s <- svd(as.matrix(x), nu = k, nv = k)
    return(list(u = s$u, v = s$v, d = s$d[seq_len(k)]))
  }
  if (solver == "svds") {
    return(suppressWarnings(RSpectra::svds(x, k = k)))
  }
  .deterministic(seed, function() withCallingHandlers(
    irlba::irlba(x, nv = k, nu = k),
    warning = function(w) {
      ## Never swallow a convergence failure: the retained subspace may be
      ## wrong, and Polaris uses all of it.
      if (grepl("did not converge", conditionMessage(w), fixed = TRUE))
        warning(sprintf(paste0("The truncated SVD did not converge (%d ",
                               "components of a %d x %d matrix). The denoised ",
                               "subspace may be unreliable."),
                        k, nrow(x), ncol(x)), call. = FALSE)
      invokeRestart("muffleWarning")
    }))
}

#' Column standard deviations, dispatching explicitly on sparseness
#'
#' `MatrixGenerics::colSds` handles both classes, but the dgCMatrix method lives
#' in sparseMatrixStats, which MatrixGenerics only loads lazily from its own
#' Suggests. Calling it directly makes the dependency real rather than
#' incidental. Both routes return identical values.
#'
#' The pre-package code called `colSds` unqualified, which resolved by
#' library() attach order: matrixStats::colSds is not generic and errors on any
#' S4 matrix, so the fit worked only because the drivers happened to attach
#' sparseMatrixStats and SparseArray afterwards.
#' @noRd
.col_sds <- function(x) {
  if (methods::is(x, "sparseMatrix")) sparseMatrixStats::colSds(x)
  else MatrixGenerics::colSds(x)
}
