#' Nearest-neighbor graphs on skymap embeddings
#'
#' @param skymap A [skymap] from [Polaris()].
#' @param slot Which embedding(s) to build a graph on. One or more of
#'   `"skymap.cell"`, `"skymap.U"`, `"skymap.V"`, `"skymap.P"`, `"skymap.Q"`, or
#'   `"all"`. Defaults to `"skymap.cell"`.
#' @param n.neighbors Neighbors per point.
#' @param metric Distance for the Annoy backend. Note the capitalized spelling
#'   BiocNeighbors expects, e.g. `"Cosine"`, which differs from uwot's
#'   `"cosine"` in [SkymapUMAP()].
#' @param mc.cores Cores, one slot per core.
#' @param verbose Report progress.
#' @param ... Further arguments for [BiocNeighbors::findKNN()].
#'
#' @return The `skymap` with results added under `skymap$NN.graph`, one entry
#'   per slot named `<slot>.nb`, each carrying `index` and `distance` matrices
#'   with rownames taken from the embedding. Repeated calls accumulate, with a
#'   same-named entry replaced, matching the storage contract the analysis
#'   scripts read (for example `skymap$NN.graph$skymap.cell.nb`).
#'
#' @examples
#' set.seed(1)
#' sim <- polarisSimulate(n = 150, g = 50, p = 60, r = 4)
#' fit <- Polaris(sim$X, sim$Y, x.sds = 0, verbose = FALSE)
#' fit <- SkymapFindNeighbors(fit, n.neighbors = 10, verbose = FALSE)
#' dim(fit$NN.graph$skymap.cell.nb$index)
#'
#' @export
SkymapFindNeighbors <- function(skymap,
                                slot = "skymap.cell",
                                n.neighbors = 300L,
                                metric = "Cosine",
                                mc.cores = 1L,
                                verbose = TRUE, ...) {

  .check_skymap(skymap)
  slot <- .resolve_slots(slot, skymap)

  res <- .maybe_mclapply(slot, function(s) {
    if (isTRUE(verbose)) message("Neighbors on ", s, ".")
    .knn(as.matrix(skymap[[s]]), n.neighbors, metric, ...)
  }, mc.cores, "neighbor graphs")

  names(res) <- paste0(slot, ".nb")
  skymap$NN.graph <- .merge_slot(skymap$NN.graph, res)
  skymap
}

#' Cross-modality nearest neighbors within each chromosome
#'
#' For every chromosome, finds each gene's nearest second-modality features and
#' each feature's nearest genes in the shared per-chromosome loading space.
#'
#' @param skymap A [skymap] fitted with `gene.chr.ref`.
#' @param chr Chromosomes to process, or `"all"`.
#' @param n.neighbors Neighbors per query.
#' @param metric Distance for the Annoy backend, e.g. `"Cosine"`.
#' @param mc.cores Cores, one chromosome per core.
#' @param verbose Report progress.
#'
#' @return The `skymap` with `skymap$NN.graph$skymap.feature.chr.nb`, a list
#'   keyed by chromosome, each holding `gene.to.feature` and `feature.to.gene`.
#'
#' @details
#' Rewritten for the `list(P, Q, stdev)` layout that `Polaris()` stores. The
#' pre-package version indexed `skymap.feature.chr[[chr]]` as a single matrix
#' and called `rownames()` on it, which is `NULL` for a list, so it could not
#' run on any current fit.
#'
#' @export
SkymapFindNeighbors.Chr <- function(skymap, chr = "all",
                                    n.neighbors = 200L,
                                    metric = "Cosine",
                                    mc.cores = 1L,
                                    verbose = TRUE) {

  .check_skymap(skymap)
  feat <- .resolve_chrs(skymap, chr)

  res <- .maybe_mclapply(names(feat), function(ch) {
    if (isTRUE(verbose)) message("Cross-modality neighbors on ", ch, ".")
    cc <- feat[[ch]]
    P <- as.matrix(cc$P); Q <- as.matrix(cc$Q)
    kg <- min(n.neighbors, nrow(P))
    kq <- min(n.neighbors, nrow(Q))
    list(
      gene.to.feature = .knn_query(Q, P, kq, metric, rownames(P), rownames(Q)),
      feature.to.gene = .knn_query(P, Q, kg, metric, rownames(Q), rownames(P)))
  }, mc.cores, "chromosomes")

  names(res) <- names(feat)
  skymap$NN.graph[["skymap.feature.chr.nb"]] <-
    .merge_slot(skymap$NN.graph[["skymap.feature.chr.nb"]], res)
  skymap
}

#' @noRd
.knn <- function(m, k, metric, ...) {
  k <- min(k, nrow(m) - 1L)
  if (k < 1L)
    stop(sprintf("Need at least 2 rows to build a neighbor graph; got %d.",
                 nrow(m)), call. = FALSE)
  r <- BiocNeighbors::findKNN(m, k = k,
                              BNPARAM = BiocNeighbors::AnnoyParam(distance = metric),
                              ...)
  rownames(r$index) <- rownames(m)
  rownames(r$distance) <- rownames(m)
  r
}

#' @noRd
.knn_query <- function(reference, query, k, metric, query.names, ref.names) {
  k <- min(k, nrow(reference))
  if (k < 1L) return(NULL)
  r <- BiocNeighbors::queryKNN(
    X = reference, query = query, k = k,
    BNPARAM = BiocNeighbors::AnnoyParam(distance = metric))
  rownames(r$index) <- query.names
  rownames(r$distance) <- query.names
  ## Translate reference row numbers into feature names, which is what a caller
  ## actually wants from a cross-modality graph.
  r$name <- matrix(ref.names[r$index], nrow = nrow(r$index),
                   dimnames = list(query.names, NULL))
  r
}
