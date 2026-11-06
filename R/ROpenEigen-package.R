#' ROpenEigen: R Bindings for Eigen 3
#'
#' A thin, object-oriented interface to the Eigen 3 C++ template library.
#'
#' @section Two kinds of object:
#'
#' `eigenMat` holds a dense matrix on the C++ side. It behaves like a base R
#' matrix for the idioms that matter: `A %*% B`, `t(A)`, `solve(A, b)`, `A[i, j]`,
#' `A + 1`, `abs(A)`, `sum(A)`, `Re(A)`, `dim(A)`. Use [as.matrix()] to bring the
#' data back to R.
#'
#' Decomposition objects (`eigenLLT`, `eigenLU`, `eigenSVD`, ...) hold a
#' factorisation that is computed once and reused. Their methods are reached
#' with `$`: `llt$solve(b)`, `svd$singular_values()`. `names(obj)` lists them.
#'
#' @section Cost model:
#'
#' Crossing the R/C++ boundary copies once, in the direction you asked for.
#' Every operation after that runs entirely inside Eigen. This is the opposite
#' trade-off from a zero-copy view: it makes R's copy-on-modify semantics hold,
#' at the price of one copy at the edge.
#'
#' @section What is deliberately missing:
#'
#' Integer matrices. Eigen's `MatrixXi` and `MatrixXd` are distinct types that do
#' not interoperate, so keeping integers would mean writing every operation
#' twice. Integer input is promoted to double, as `RcppEigen::as<MatrixXd>`
#' does.
#'
#' Sparse matrices are not wrapped yet, although the headers are installed.
#'
#' @keywords internal
#' @useDynLib ROpenEigen, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
