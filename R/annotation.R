#' Gene annotation for POLARIS linkage
#'
#' Returns a `GRanges` of gene bodies with a `gene_name` column and strand,
#' suitable for the `gene.chr.ref` argument of [Polaris()] and the `gene.anno`
#' argument of [SkymapLinkageTable()].
#'
#' @param genome Name of a bundled reference table. Currently only `"hg38"` is
#'   bundled; call with any other value to see what is available. For every other
#'   genome and build use [polarisMakeTSSRef()].
#'
#' @return A stranded `GRanges` with `gene_name` and `gene_type` columns.
#'
#' @details
#' **Prefer the annotation you already have.** The annotation must match the
#' genome build the peaks were called against, and a Signac ATAC assay already
#' carries exactly the right object:
#'
#' ```
#' anno <- Signac::Annotation(obj[["ATAC"]])
#' ```
#'
#' That is guaranteed consistent with your peaks, which a bundled table is not.
#' Use `polarisTSSRef()` when no annotation is attached, and
#' [polarisMakeTSSRef()] for any build the package does not bundle.
#'
#' The bundled tables carry only the columns POLARIS reads: gene name, sequence
#' name, TSS-bearing coordinate, strand and gene type. They are derived by
#' `data-raw/make_tss_ref.R`, which records their provenance.
#'
#' @examples
#' anno <- polarisTSSRef("hg38")
#' head(anno)
#'
#' @seealso [polarisMakeTSSRef()]
#' @export
polarisTSSRef <- function(genome = "hg38") {
  dir <- system.file("extdata", package = "POLARIS")
  have <- sub("^tss_", "", sub("\\.rds$", "",
              list.files(dir, pattern = "^tss_.*\\.rds$")))
  f <- file.path(dir, paste0("tss_", genome, ".rds"))
  if (!file.exists(f))
    stop(sprintf(paste0("No bundled annotation for '%s'. Bundled: %s. For any ",
                        "other genome or build, construct one from a TxDb, EnsDb ",
                        "or GTF with polarisMakeTSSRef()."),
                 genome, paste(have, collapse = ", ")), call. = FALSE)
  d <- readRDS(f)
  GenomicRanges::GRanges(
    seqnames = d$chr,
    ranges   = IRanges::IRanges(start = d$start, end = d$end),
    strand   = d$strand,
    gene_name = d$gene,
    gene_type = d$gene_type)
}

#' Build a gene annotation from a TxDb, EnsDb or GTF
#'
#' Turns any gene-level annotation source into the `GRanges` that POLARIS
#' expects, so the package is not limited to the builds it bundles.
#'
#' @param x One of: a `TxDb` or `EnsDb` object; a `GRanges` that already has
#'   gene ranges; or a path to a GTF or GFF file.
#' @param gene.name.column For a `GRanges` or GTF input, the metadata column
#'   holding gene symbols. Default `"gene_name"`, falling back to `"gene_id"`.
#'
#' @return A stranded `GRanges` with a `gene_name` column.
#'
#' @details
#' Strand is essential. POLARIS places transcription start sites with
#' `GenomicRanges::resize(fix = "start")`, which is strand-aware, so an
#' unstranded annotation silently assigns minus-strand genes the wrong end of
#' the gene. This function errors if strand information is absent.
#'
#' @examples
#' \dontrun{
#' # From a Bioconductor annotation package
#' anno <- polarisMakeTSSRef(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
#'
#' # From the GTF that the reads were aligned against, which is the safest source
#' anno <- polarisMakeTSSRef("refdata-cellranger-arc-GRCh38/genes/genes.gtf.gz")
#' }
#'
#' @seealso [polarisTSSRef()]
#' @export
polarisMakeTSSRef <- function(x, gene.name.column = "gene_name") {

  gr <- if (is.character(x)) {
    if (!requireNamespace("rtracklayer", quietly = TRUE))
      stop(paste0("Reading a GTF needs the rtracklayer package. Install it, or ",
                  "pass a TxDb, EnsDb or GRanges instead."), call. = FALSE)
    g <- rtracklayer::import(x)
    if ("type" %in% names(S4Vectors::mcols(g))) g <- g[g$type == "gene"]
    g
  } else if (methods::is(x, "GRanges")) {
    x
  } else if (methods::is(x, "TxDb") || methods::is(x, "EnsDb")) {
    if (!requireNamespace("GenomicFeatures", quietly = TRUE))
      stop("Extracting genes from a TxDb/EnsDb needs the GenomicFeatures package.",
           call. = FALSE)
    GenomicFeatures::genes(x)
  } else {
    stop(paste0("`x` must be a TxDb, an EnsDb, a GRanges, or a path to a ",
                "GTF/GFF file."), call. = FALSE)
  }

  mc <- names(S4Vectors::mcols(gr))
  col <- if (gene.name.column %in% mc) gene.name.column
  else if ("gene_name" %in% mc) "gene_name"
  else if ("symbol" %in% mc) "symbol"
  else if ("gene_id" %in% mc) "gene_id"
  else NA_character_
  if (is.na(col))
    stop(sprintf(paste0("No gene-name column found. Available: %s. Set ",
                        "gene.name.column."), paste(mc, collapse = ", ")),
         call. = FALSE)

  gr$gene_name <- as.character(S4Vectors::mcols(gr)[[col]])

  st <- as.character(GenomicRanges::strand(gr))
  if (all(st == "*"))
    stop(paste0("`x` carries no strand information. POLARIS needs strand to ",
                "place transcription start sites; without it every minus-strand ",
                "gene would get the wrong end of the gene as its TSS."),
         call. = FALSE)
  if (any(st == "*"))
    warning(sprintf(paste0("%d ranges are unstranded ('*') and will be treated ",
                           "as plus-strand."), sum(st == "*")), call. = FALSE)

  gr <- gr[!is.na(gr$gene_name) & nzchar(gr$gene_name)]
  if (!length(gr))
    stop("No usable gene ranges after filtering.", call. = FALSE)
  gr
}
