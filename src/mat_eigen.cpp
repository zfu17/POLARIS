// Dense and sparse-dense matrix products used by Polaris() for the feature-loading
// projections. Ported from POLAR/scripts/RCPP/mat_eigen.cpp with three changes:
//   1. Eigen::MappedSparseMatrix -> Eigen::Map<Eigen::SparseMatrix<double> >
//      (MappedSparseMatrix is deprecated and slated for removal in Eigen 4).
//   2. An explicit conformability check on the sparse-dense product. Eigen's own
//      assertion is compiled out because R defines -DNDEBUG, so a mismatched
//      product was undefined behaviour rather than an error.
//   3. matrixMultSparseSparse, svdEigen and covEigen are not carried over; they
//      were never called from any analysis script. svdEigen in particular returned
//      V transposed in a field named "V", which would silently mislead any future
//      caller using the base svd() idiom.
// Numerics are unchanged: both kernels are plain Eigen products, as before.

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// Types in the exported signatures are fully qualified on purpose. Rcpp
// generates src/RcppExports.cpp as a separate translation unit that reproduces
// these signatures verbatim without any using-declarations from this file, so
// bare `Map<MatrixXd>` fails to compile there. The pre-package version was
// sourceCpp'd as a single file, where the using-declarations were in scope.

// [[Rcpp::export]]
Eigen::MatrixXd matrixMultiplyEigen(const Eigen::Map<Eigen::MatrixXd> &A,
                                    const Eigen::Map<Eigen::MatrixXd> &B) {
  if (A.cols() != B.rows()) {
    Rcpp::stop("Non-conformable arguments: A is %dx%d, B is %dx%d.",
               A.rows(), A.cols(), B.rows(), B.cols());
  }
  return A * B;
}

// [[Rcpp::export]]
Eigen::MatrixXd matrixMultSparseDense(
    const Eigen::Map<Eigen::SparseMatrix<double> > &A,
    const Eigen::Map<Eigen::MatrixXd> &B) {
  if (A.cols() != B.rows()) {
    Rcpp::stop("Non-conformable arguments: A is %dx%d, B is %dx%d.",
               A.rows(), A.cols(), B.rows(), B.cols());
  }
  return A * B;
}
