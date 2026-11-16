\name{eigen_matrix_theory}
\alias{eigen_matrix_theory}
\alias{eigen_rref}
\alias{eigen_nullspace}
\alias{eigen_colspace}
\alias{eigen_rowspace}
\alias{eigen_gram}
\alias{eigen_gram_schmidt}
\alias{eigen_projection}
\alias{eigen_leading_minors}
\alias{eigen_adjugate}
\alias{eigen_charpoly}
\alias{eigen_companion}
\alias{eigen_doolittle}
\alias{eigen_ldlt_simple}
\alias{eigen_jordan}
\alias{eigen_kronecker}
\alias{eigen_hadamard}
\alias{eigen_vec}
\alias{eigen_inertia}
\alias{eigen_spectral_radius}
\alias{eigen_cond}
\alias{eigen_polar}
\alias{eigen_full_rank}
\alias{eigen_ginverse}
\title{Matrix-theory helpers built on ROpenEigen}
\description{
  R-level helpers for the algorithms that appear in matrix-theory textbooks but
  are deliberately not provided by Eigen: row reduction, null spaces, Gram
  matrices, characteristic polynomials, Jordan block shapes, matrix products,
  inertia, condition numbers, full-rank factorisations and generalized inverses.
}
\usage{
eigen_rref(x, tol = 1e-10)
eigen_nullspace(x, tol = 1e-10)
eigen_colspace(x, tol = 1e-10)
eigen_rowspace(x, tol = 1e-10)
eigen_gram(x)
eigen_gram_schmidt(x, tol = 1e-12)
eigen_projection(x, tol = 1e-10)
eigen_leading_minors(x)
eigen_adjugate(x)
eigen_charpoly(x)
eigen_companion(coeffs)
eigen_doolittle(x, tol = 1e-12)
eigen_ldlt_simple(x, tol = 1e-12)
eigen_jordan(x, tol = 1e-8)
eigen_kronecker(a, b)
eigen_hadamard(a, b)
eigen_vec(x, byrow = FALSE)
eigen_inertia(x, tol = 1e-10)
eigen_spectral_radius(x)
eigen_cond(x, p = 2)
eigen_polar(x)
eigen_full_rank(x, tol = 1e-10)
eigen_ginverse(x, tol = 1e-10)
}
\arguments{
  \item{x, a, b}{A matrix or an \code{eigenMat} object.}
  \item{coeffs}{Characteristic-polynomial coefficients in descending order.}
  \item{tol}{Numerical tolerance.}
  \item{p}{Norm type for \code{eigen_cond}; 2 uses the reciprocal condition number.}
  \item{byrow}{Vectorise by rows instead of by columns.}
}
\value{
  Lists or matrices as described in the package README and the matrix-theory
  notebooks that use them.
}
\examples{
A <- matrix(c(3, 2, -1, -1, 0, 0, -1, 3, 0), 3, byrow = TRUE)
eigen_rref(A)$rank
eigen_charpoly(A)
eigen_doolittle(A)
}
