# 矩阵分解
#
# 每个分解既是「一次性函数」（eigen_llt(A)$solve(b)），也是「可复用对象」
# （llt <- eigen_llt(A); llt$solve(b1); llt$solve(b2)）。这是 Eigen 相对
# base R 的 solve() 的核心优势：solve(A, b) 每次都重新分解，而分解对象
# 只算一次。
#
# 方法表注册在 .reigen_cls 里，用 $ 取用。方法一律以 self 为第一形参。

# ---------------------------------------------------------------- 构造

#' Matrix decompositions
#'
#' Each function factorises a matrix once and returns a decomposition object.
#' The object keeps the factorisation on the C++ side, so repeated
#' `$solve()` calls reuse it instead of refactorising.
#'
#' `eigen_llt()` requires a symmetric positive-definite matrix and stops
#' otherwise, as [base::chol()] does; `eigen_ldlt()` only requires symmetry.
#' Both read **only the lower triangle** -- passing an asymmetric matrix will
#' not raise an error, it will silently be treated as symmetric. Use
#' [eigen_svd()] or [eigen_full_lu()] when you do not know the structure of the
#' matrix.
#'
#' The `eigen_llt()` check matters more than it looks: Eigen reports a failed
#' Cholesky only through `info()`, while `$solve()` and `$reconstructed()` keep
#' returning values -- factorisations of a *different* matrix. Letting that
#' through would turn a definiteness problem into a silently wrong answer.
#'
#' `eigen_eigen()` is the self-adjoint (symmetric / Hermitian) eigensolver and
#' also reads only the lower triangle. For a general square matrix use
#' `eigen_geigen()`, which always returns complex eigenvalues.
#'
#' `eigen_svd()` is the only decomposition that accepts a rectangular matrix
#' *and* a rank-deficient one without complaint. It has two backends: Eigen's
#' `JacobiSVD`, which is accurate but has a very large constant factor, and
#' `BDCSVD`, a divide-and-conquer variant that is roughly 4x faster at
#' 100x100, 8x at 200x200 and 16x at 400x400. Both inherit the same
#' `threshold()`/`rank()` implementation, so switching backends does not move
#' the rank cut-off or the pseudo-inverse truncation point. `method = "auto"`
#' (the default) picks `"bdc"` when `min(dim(x)) >= 16` and `"jacobi"`
#' below that, where the difference is not worth a second code path.
#' `$backend()` reports which one actually ran.
#'
#' `eigen_full_lu()` accepts a rectangular matrix too. `solve()`, `rank()`,
#' `matrixL()`, `matrixU()`, `rcond()` and the
#' permutation matrices all work there; `det()` and `inverse()` are square-only
#' and raise an error otherwise -- Eigen defines them only for square input,
#' and calling them on a rectangular factorisation reads out of bounds.
#'
#' @param x A matrix, `eigenMat`, or anything [eigenMat()] accepts. It must be
#'   finite: decompositions reject `NaN` and `Inf`, as [base::svd()] and
#'   [base::qr()] do. Use [eigen_all_finite()] to test first.
#' @param method For [eigen_svd()] only: `"auto"` (default), `"jacobi"` or
#'   `"bdc"`. All three produce the same factorisation; they differ in speed.
#' @return A decomposition object. Its methods are listed by `names()`.
#' @name eigen-decompositions
NULL

#' @rdname eigen-decompositions
#' @export
eigen_llt <- function(x) reigen_llt_new(.as_eigen(x))

#' @rdname eigen-decompositions
#' @export
eigen_ldlt <- function(x) reigen_ldlt_new(.as_eigen(x))

#' @rdname eigen-decompositions
#' @export
eigen_lu <- function(x) reigen_lu_new(.as_eigen(x))

#' @rdname eigen-decompositions
#' @export
eigen_full_lu <- function(x) reigen_fulllu_new(.as_eigen(x))

#' @rdname eigen-decompositions
#' @export
eigen_qr <- function(x) reigen_qr_new(.as_eigen(x))

#' @rdname eigen-decompositions
#' @export
eigen_colpiv_qr <- function(x) reigen_colpivqr_new(.as_eigen(x))

#' @rdname eigen-decompositions
#' @export
eigen_svd <- function(x, method = c("auto", "jacobi", "bdc")) {
  method <- match.arg(method)
  reigen_svd_new(.as_eigen(x), match(method, c("auto", "jacobi", "bdc")) - 1L)
}

#' @rdname eigen-decompositions
#' @export
eigen_eigen <- function(x) reigen_eigen_new(.as_eigen(x))

#' @rdname eigen-decompositions
#' @export
eigen_geigen <- function(x) reigen_ceigen_new(.as_eigen(x))

# ---------------------------------------------------------------- 缺省右端项

# 分解对象上没有 dim，形状挂在 C++ 构造时写的 eigenDim 属性上。
.eigen_dim <- function(x) {
  d <- attr(x, "eigenDim")
  if (is.null(d)) stop("object carries no dimension attribute", call. = FALSE)
  d
}

.eigen_rhs <- function(self, b) {
  if (!is.null(b)) return(.as_eigen(b))
  n <- .eigen_dim(self)[1]
  reigen_mat_identity(n)
}

# ---------------------------------------------------------------- 方法注册

.reigen_reg("eigenLLT", m = list(
  solve         = function(self, b = NULL) reigen_llt_solve(self, .eigen_rhs(self, b)),
  matrixL       = function(self) reigen_llt_matrixL(self),
  matrixU       = function(self) reigen_llt_matrixU(self),
  reconstructed = function(self) reigen_llt_reconstructed(self),
  logdet        = function(self) reigen_llt_logdet(self),
  rcond         = function(self) reigen_llt_rcond(self),
  is_ok         = function(self) reigen_llt_ok(self)
))

.reigen_reg("eigenLDLT", m = list(
  solve         = function(self, b = NULL) reigen_ldlt_solve(self, .eigen_rhs(self, b)),
  matrixL       = function(self) reigen_ldlt_matrixL(self),
  matrixU       = function(self) reigen_ldlt_matrixU(self),
  vectorD       = function(self) reigen_ldlt_vectorD(self),
  reconstructed = function(self) reigen_ldlt_reconstructed(self),
  logdet        = function(self) reigen_ldlt_logdet(self),
  rcond         = function(self) reigen_ldlt_rcond(self),
  is_ok         = function(self) reigen_ldlt_ok(self)
))

.reigen_reg("eigenLU", m = list(
  solve         = function(self, b = NULL) reigen_lu_solve(self, .eigen_rhs(self, b)),
  matrixL       = function(self) reigen_lu_matrixL(self),
  matrixU       = function(self) reigen_lu_matrixU(self),
  permutationP  = function(self) reigen_lu_permutationP(self),
  inverse       = function(self) reigen_lu_inverse(self),
  reconstructed = function(self) reigen_lu_reconstructed(self),
  det           = function(self) reigen_lu_det(self)
))

.reigen_reg("eigenFullPivLU", m = list(
  solve         = function(self, b = NULL) reigen_fulllu_solve(self, .eigen_rhs(self, b)),
  matrixL       = function(self) reigen_fulllu_matrixL(self),
  matrixU       = function(self) reigen_fulllu_matrixU(self),
  permutationP  = function(self) reigen_fulllu_permutationP(self),
  permutationQ  = function(self) reigen_fulllu_permutationQ(self),
  inverse       = function(self) reigen_fulllu_inverse(self),
  reconstructed = function(self) reigen_fulllu_reconstructed(self),
  det           = function(self) reigen_fulllu_det(self),
  rank          = function(self) reigen_fulllu_rank(self),
  rcond         = function(self) reigen_fulllu_rcond(self),
  is_invertible = function(self) reigen_fulllu_invertible(self)
))

.reigen_reg("eigenQR", m = list(
  solve         = function(self, b = NULL) reigen_qr_solve(self, .eigen_rhs(self, b)),
  matrixQ       = function(self) reigen_qr_matrixQ(self),
  matrixR       = function(self) reigen_qr_matrixR(self),
  reconstructed = function(self) reigen_qr_reconstructed(self),
  log_abs_det   = function(self) reigen_qr_log_abs_det(self)
))

.reigen_reg("eigenColPivQR", m = list(
  solve         = function(self, b = NULL) reigen_colpivqr_solve(self, .eigen_rhs(self, b)),
  matrixQ       = function(self) reigen_colpivqr_matrixQ(self),
  matrixR       = function(self) reigen_colpivqr_matrixR(self),
  permutationP  = function(self) reigen_colpivqr_permutationP(self),
  reconstructed = function(self) reigen_colpivqr_reconstructed(self),
  rank          = function(self) reigen_colpivqr_rank(self),
  max_pivot     = function(self) reigen_colpivqr_max_pivot(self),
  kernel_dim    = function(self) reigen_colpivqr_kernel_dim(self),
  is_invertible = function(self) reigen_colpivqr_invertible(self),
  log_abs_det   = function(self) reigen_colpivqr_log_abs_det(self)
))

.reigen_reg("eigenSVD", m = list(
  solve          = function(self, b = NULL) reigen_svd_solve(self, .eigen_rhs(self, b)),
  singular_values = function(self) reigen_svd_singular_values(self),
  u              = function(self) reigen_svd_matrixU(self),
  v              = function(self) reigen_svd_matrixV(self),
  pseudo_inverse = function(self) reigen_svd_pseudo_inverse(self),
  reconstructed  = function(self) reigen_svd_reconstructed(self),
  rank           = function(self) reigen_svd_rank(self),
  nonzero        = function(self) reigen_svd_nonzero(self),
  threshold      = function(self) reigen_svd_threshold(self),
  backend        = function(self) reigen_svd_backend(self)
))

.reigen_reg("eigenSelfAdjointEigen", m = list(
  values  = function(self) reigen_eigen_values(self),
  vectors = function(self) reigen_eigen_vectors(self),
  solve   = function(self, b = NULL) reigen_eigen_solve(self, .eigen_rhs(self, b)),
  is_ok   = function(self) reigen_eigen_ok(self)
))

.reigen_reg("eigenComplexEigen", m = list(
  values  = function(self) reigen_ceigen_values(self),
  vectors = function(self) reigen_ceigen_vectors(self),
  is_ok   = function(self) reigen_ceigen_ok(self)
))

# ---------------------------------------------------------------- 泛型接入

#' @export
dim.eigenObject <- function(x) .eigen_dim(x)

# solve() 作用在分解对象上时，直接走它的 solve 方法：
# solve(llt, b) 与 llt$solve(b) 等价。base R 的 solve(A, b) 每次都重分解，
# 这里不是——这是把分解对象传进来的全部意义。
# 与 solve.eigenMat 一致：b 是无 dim 的向量时返回向量。
#' @export
solve.eigenObject <- function(a, b, ...) {
  if (missing(b)) return(a$solve())
  drop <- is.null(dim(b)) && !inherits(b, "eigenMat")
  r <- a$solve(b)
  if (drop) as.vector(r) else r
}

# ---------------------------------------------------------------- 一次性快捷方式

#' Eigenvalues of a symmetric matrix
#'
#' Shorthand for `eigen_eigen(x)$values()`. Eigenvalues come back in ascending
#' order, as Eigen returns them.
#'
#' @param x A symmetric matrix or `eigenMat`.
#' @return An `eigenMat` column vector.
#' @export
eigen_values <- function(x) reigen_eigen_values(eigen_eigen(x))

#' @rdname eigen_values
#' @export
eigen_vectors <- function(x) reigen_eigen_vectors(eigen_eigen(x))

#' Eigenvalues of a general square matrix
#'
#' Always complex, even when the input is real: the eigenvalues of a real matrix
#' need not be real. Use [eigen_values()] when the matrix is known symmetric --
#' it is faster and returns real numbers.
#'
#' @param x A square matrix or `eigenMat`.
#' @return An `eigenMat` complex column vector.
#' @export
eigen_gvalues <- function(x) reigen_ceigen_values(eigen_geigen(x))

#' @rdname eigen_gvalues
#' @export
eigen_gvectors <- function(x) reigen_ceigen_vectors(eigen_geigen(x))

#' Singular values
#' @param x A matrix or `eigenMat`.
#' @return An `eigenMat` column vector, in descending order.
#' @export
eigen_singular_values <- function(x) reigen_svd_singular_values(eigen_svd(x))

#' Moore-Penrose pseudo-inverse
#' @param x A matrix or `eigenMat`.
#' @return An `eigenMat`.
#' @export
eigen_pinv <- function(x) reigen_svd_pseudo_inverse(eigen_svd(x))

#' Cholesky factor
#'
#' @param x A symmetric positive-definite matrix.
#' @param upper Return the upper factor instead of the lower one.
#' @return An `eigenMat`.
#' @export
eigen_chol <- function(x, upper = FALSE) {
  # 非正定会在这里直接报错（见 src/reigen_decomp.cpp 的 reigen_llt_new），
  # 所以不需要再查 is_ok()。
  llt <- eigen_llt(x)
  if (upper) llt$matrixU() else llt$matrixL()
}
