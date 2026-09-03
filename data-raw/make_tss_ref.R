## Build the slim gene-annotation tables bundled in inst/extdata.
##
## Provenance
## ----------
## hg38: derived from the Cell Ranger ARC GRCh38 reference GTF (the annotation
## the POLARIS manuscript's 10x multiome data were aligned and quantified
## against), reduced to gene-level ranges. The source object used here is the
## project's cached GRanges:
##   /n/holystore01/LABS/xlin/Lab/ziqifu/tools/gtf/cellranger_hg38_genetype.rds
## To regenerate from the GTF directly instead:
##   gr <- rtracklayer::import("refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz")
##   gr <- gr[gr$type == "gene"]
##
## Only the columns POLARIS reads are kept: gene symbol, sequence name, the gene
## body's start and end, strand, and gene type. Strand is essential, because
## GenomicRanges::resize(fix = "start") is strand-aware and is how POLARIS
## places transcription start sites. All 39 sequence levels are retained,
## including scaffolds; Polaris() itself restricts to the primary assembly.
##
## Run from the package root with R >= 4.4.

suppressMessages(library(GenomicRanges))

slim <- function(gr) {
  stopifnot(!any(as.character(strand(gr)) == "*"))
  data.frame(
    gene      = as.character(gr$gene_name),
    chr       = as.character(seqnames(gr)),
    start     = start(gr),
    end       = end(gr),
    strand    = as.character(strand(gr)),
    gene_type = as.character(gr$gene_type),
    stringsAsFactors = FALSE)
}

src <- "/n/holystore01/LABS/xlin/Lab/ziqifu/tools/gtf/cellranger_hg38_genetype.rds"
if (!file.exists(src))
  stop("Source annotation not found: ", src)

hg38 <- slim(readRDS(src))
message(sprintf("hg38: %d genes, %d sequence levels, strand: %s",
                nrow(hg38), length(unique(hg38$chr)),
                paste(names(table(hg38$strand)), table(hg38$strand),
                      sep = "=", collapse = " ")))

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
saveRDS(hg38, "inst/extdata/tss_hg38.rds", compress = "xz")
message(sprintf("wrote inst/extdata/tss_hg38.rds (%s KB)",
                round(file.size("inst/extdata/tss_hg38.rds") / 1024)))

## mm10 is not bundled yet. To add it:
##   ens  <- EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
##   gr   <- GenomicFeatures::genes(ens)
##   gr$gene_name <- gr$symbol
##   gr$gene_type <- gr$gene_biotype
##   seqlevelsStyle(gr) <- "UCSC"
##   saveRDS(slim(gr), "inst/extdata/tss_mm10.rds", compress = "xz")
