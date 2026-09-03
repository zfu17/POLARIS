#' Fit POLARIS to a Seurat object
#'
#' Convenience wrapper that pulls the two modality matrices out of a Seurat
#' object, orients them as POLARIS expects, and calls [Polaris()]. It also
#' returns the joint cell embedding as a Seurat `DimReduc`, so the fit can be
#' used directly with `FindNeighbors`, `FindClusters` and `RunUMAP`.
#'
#' `Polaris()` itself takes plain matrices and is the reproducible entry point.
#' This function only assembles its arguments.
#'
#' @param object A `Seurat` object with both modalities as assays.
#' @param assay.x,assay.y Assay names. Defaults `"RNA"` and `"ATAC"`.
#' @param layer.x,layer.y Layer (slot) to take from each assay. The defaults are
#'   the conventions used throughout the POLARIS manuscript: `"scale.data"` for
#'   RNA, so \eqn{X} is centered and scaled, and `"data"` for ATAC, which is
#'   expected to hold TF-IDF normalized values from `Signac::RunTFIDF`. For
#'   CITE-seq use `assay.y = "ADT"` with CLR-normalized `"data"`.
#' @param gene.chr.ref Gene annotation. If `NULL` and `assay.y` carries a Signac
#'   annotation, that is used, which is the safest source because it is
#'   guaranteed to match the genome the peaks were called on. Pass a `GRanges`
#'   to override, or `FALSE` to skip the per-chromosome step.
#' @param reduction.name Name under which to store the joint cell embedding.
#' @param ... Passed to [Polaris()], including `x.sds`, `T1`, `T2`, `T3`,
#'   `mc.cores` and `seed`.
#'
#' @return A list with `skymap` (the [skymap] fit) and `object` (the input
#'   Seurat object with the joint embedding added as `reduction.name`).
#'
#' @section Preprocessing:
#' POLARIS does not normalize for you, because the normalization is a
#' methodological choice that belongs to your analysis. The conventions used in
#' the manuscript are:
#'
#' ```
#' obj <- NormalizeData(obj, assay = "RNA",
#'                      normalization.method = "LogNormalize",
#'                      scale.factor = 1e6)
#' obj <- FindVariableFeatures(obj, nfeatures = 3000)
#' # NOTE: ScaleData is deliberately run on ALL genes, not on the variable
#' # features. POLARIS selects its own feature space by the x.sds variance
#' # floor, which is broader than a highly-variable-gene list.
#' obj <- ScaleData(obj, assay = "RNA", features = rownames(obj))
#' DefaultAssay(obj) <- "ATAC"
#' obj <- Signac::RunTFIDF(obj)
#' ```
#'
#' Note the interaction with `x.sds`: its default of `0.95` assumes \eqn{X} is
#' z-scored. If you pass a log-normalized layer instead, nearly every gene falls
#' below the floor and the fit will stop with an empty \eqn{X}. Pass
#' `x.sds = 0` in that case.
#'
#' @examples
#' \dontrun{
#' res <- RunPolaris(obj, assay.x = "RNA", assay.y = "ATAC", mc.cores = 4)
#' res$object <- Seurat::FindNeighbors(res$object, reduction = "polaris",
#'                                     dims = 1:ncol(res$skymap$skymap.cell))
#' res$skymap <- SkymapSimScore(res$skymap)
#' res$object$concordance <- res$skymap$concordance[colnames(res$object)]
#' }
#'
#' @seealso [Polaris()]
#' @export
RunPolaris <- function(object,
                       assay.x = "RNA", assay.y = "ATAC",
                       layer.x = "scale.data", layer.y = "data",
                       gene.chr.ref = NULL,
                       reduction.name = "polaris",
                       ...) {

  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("RunPolaris() needs the Seurat package.", call. = FALSE)
  if (!methods::is(object, "Seurat"))
    stop("`object` must be a Seurat object. For matrices use Polaris().",
         call. = FALSE)

  for (a in c(assay.x, assay.y))
    if (!a %in% names(object@assays))
      stop(sprintf("No assay named '%s'. Available: %s.", a,
                   paste(names(object@assays), collapse = ", ")), call. = FALSE)

  grab <- function(assay, layer) {
    m <- SeuratObject::LayerData(object, assay = assay, layer = layer)
    if (is.null(m) || !length(m))
      stop(sprintf(paste0("Layer '%s' of assay '%s' is empty. For scale.data, ",
                          "run ScaleData(obj, assay = '%s', features = ",
                          "rownames(obj)) first."), layer, assay, assay),
           call. = FALSE)
    Matrix::t(m)                      # features x cells -> cells x features
  }

  X <- grab(assay.x, layer.x)
  Y <- grab(assay.y, layer.y)
  ## Seurat can hold different cell sets per assay; align on the intersection in
  ## the object's own column order.
  common <- intersect(rownames(X), rownames(Y))
  if (!length(common))
    stop("The two assays share no cells.", call. = FALSE)
  if (length(common) < nrow(X) || length(common) < nrow(Y))
    message(sprintf("Using the %d cells present in both assays.", length(common)))
  X <- X[common, , drop = FALSE]
  Y <- Y[common, , drop = FALSE]

  if (isFALSE(gene.chr.ref)) {
    gene.chr.ref <- NULL
  } else if (is.null(gene.chr.ref)) {
    gene.chr.ref <- tryCatch({
      a <- object[[assay.y]]
      if (!is.null(methods::slot(a, "annotation"))) methods::slot(a, "annotation") else NULL
    }, error = function(e) NULL)
    if (!is.null(gene.chr.ref)) {
      message("Using the annotation attached to assay '", assay.y, "'.")
      if (!"gene_name" %in% names(S4Vectors::mcols(gene.chr.ref)))
        stop(paste0("The attached annotation has no `gene_name` column. Pass ",
                    "gene.chr.ref explicitly, e.g. polarisTSSRef('hg38')."),
             call. = FALSE)
    } else {
      message(paste0("No annotation attached to assay '", assay.y, "', so the ",
                     "per-chromosome step is skipped. Pass gene.chr.ref to ",
                     "enable SkymapLinkageTable()."))
    }
  }

  sk <- Polaris(X, Y, gene.chr.ref = gene.chr.ref, ...)

  object[[reduction.name]] <- SeuratObject::CreateDimReducObject(
    embeddings = sk$skymap.cell,
    stdev      = as.numeric(sk$stdev.cell),
    key        = "skyPC_",
    assay      = assay.x)

  list(skymap = sk, object = object)
}
