#' POLARIS: concordance-aware integration of paired single-cell multiomic data
#'
#' POLARIS takes two matrices measured on the same cells and returns, from one
#' closed-form spectral decomposition, a per-cell cross-modality concordance
#' score, a joint cell embedding, and a cross-modal feature embedding that
#' scores candidate feature pairs with an uncertainty.
#'
#' @section Entry point:
#' [Polaris()] performs the fit and returns a [skymap] object. Every other
#' function in the package consumes that object:
#'
#' \describe{
#'   \item{[SkymapSimScore()]}{per-cell concordance score}
#'   \item{[SkymapUMAP()], [SkymapFindNeighbors()]}{embeddings and graphs}
#'   \item{[SkymapLinkageTable()]}{cis gene-peak linkage with jackknife SE}
#'   \item{[SkymapZdev()]}{genomic-distance-aware z-statistic}
#' }
#'
#' @section Notation:
#' The code follows the notation of the POLARIS manuscript. \eqn{X} (`n` cells
#' by `g` genes) and \eqn{Y} (`n` cells by `p` features) are the centered and
#' scaled modality matrices. \eqn{\hat{U}(r)} and \eqn{\hat{V}(r)} are the
#' paired cell embeddings, \eqn{W = [\hat{U} | \hat{V}]} the joint cell
#' embedding, and \eqn{P}, \eqn{Q}, \eqn{d} the gene loadings, second-modality
#' loadings and singular-value weights of the joint feature embedding.
#'
#' @keywords internal
#' @useDynLib POLARIS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
