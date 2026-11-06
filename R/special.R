# 特殊函数与矩阵函数
#
# 这两组都来自 Eigen 的 unsupported/ 模块，但性质完全不同：
#
#   - 特殊函数是逐元素的。eigen_erf(x) 对矩阵的每个元素取 erf，
#     与 base R 的 sin()/gamma() 一样，形状不变。
#   - 矩阵函数不是逐元素的。eigen_expm(A) 是矩阵指数 sum A^k/k!，
#     不是每个元素取 exp；所以它们名字里带 m，且只接受方阵。
#
# 全部返回 eigenMat，可以直接接 %*%、solve() 等后续运算。

# ---------------------------------------------------------------- 特殊函数

#' Special functions from Eigen's SpecialFunctions module
#'
#' Elementwise special functions.  Each one takes a matrix (or anything
#' [eigenMat()] accepts) and returns an `eigenMat` of the same shape; scalars
#' broadcast the way they do everywhere else in the package.
#'
#' The functions are Eigen's, not base R's, but the definitions are the standard
#' ones: `eigen_erf()` is the error function, `eigen_erfc()` its complement,
#' `eigen_gammaln()` the log-gamma, `eigen_digamma()` the logarithmic
#' derivative, `eigen_ndtri()` the inverse of the standard normal CDF (so
#' `eigen_ndtri(0.975)` is about 1.96), `eigen_igamma(a, x)` and
#' `eigen_igammac(a, x)` the regularised (upper) incomplete gamma functions,
#' `eigen_betainc(a, b, x)` the regularised incomplete beta, `eigen_zeta(x, q)`
#' the Hurwitz zeta and `eigen_polygamma(n, x)` the n-th polygamma.
#'
#' `eigen_bessel()` covers the twelve Bessel functions Eigen implements; the
#' trailing `e` in a type name means the exponentially scaled variant (see the
#' Eigen documentation for the exact definition), which is what you want when
#' the unscaled value overflows.
#'
#' Complex input is rejected: Eigen implements these for real scalars only.
#'
#' @param x,a,b,n,q Numeric input.  In the two-argument functions the first
#'   argument is the parameter and the second the value, matching Eigen:
#'   `eigen_igamma(a, x)`, `eigen_zeta(x, q)`, `eigen_polygamma(n, x)`.
#' @param type For `eigen_bessel()`: one of `"i0"`, `"i0e"`, `"i1"`, `"i1e"`,
#'   `"k0"`, `"k0e"`, `"k1"`, `"k1e"`, `"j0"`, `"j1"`, `"y0"`, `"y1"`.
#' @return An `eigenMat` with the same dimensions as the (recycled) input.
#' @name eigen-special
NULL

#' @rdname eigen-special
#' @export
eigen_erf <- function(x) reigen_special(.as_eigen(x), 0L)

#' @rdname eigen-special
#' @export
eigen_erfc <- function(x) reigen_special(.as_eigen(x), 1L)

#' @rdname eigen-special
#' @export
eigen_gammaln <- function(x) reigen_special(.as_eigen(x), 2L)

#' @rdname eigen-special
#' @export
eigen_digamma <- function(x) reigen_special(.as_eigen(x), 3L)

#' @rdname eigen-special
#' @export
eigen_ndtri <- function(x) reigen_special(.as_eigen(x), 4L)

#' @rdname eigen-special
#' @export
eigen_igamma <- function(a, x) reigen_special2(.as_eigen(a), .as_eigen(x), 0L)

#' @rdname eigen-special
#' @export
eigen_igammac <- function(a, x) reigen_special2(.as_eigen(a), .as_eigen(x), 1L)

#' @rdname eigen-special
#' @export
eigen_zeta <- function(x, q = 1) reigen_special2(.as_eigen(x), .as_eigen(q), 2L)

#' @rdname eigen-special
#' @export
eigen_polygamma <- function(n, x) reigen_special2(.as_eigen(n), .as_eigen(x), 3L)

#' @rdname eigen-special
#' @export
eigen_betainc <- function(a, b, x) reigen_betainc(a, b, .as_eigen(x))

#' @rdname eigen-special
#' @export
eigen_bessel <- function(x, type = "i0") {
  codes <- c(i0 = 10L, i0e = 11L, i1 = 12L, i1e = 13L,
             k0 = 14L, k0e = 15L, k1 = 16L, k1e = 17L,
             j0 = 18L, j1 = 19L, y0 = 20L, y1 = 21L)
  if (length(type) != 1L || !type %in% names(codes)) {
    stop("type must be one of: ", paste(names(codes), collapse = ", "), call. = FALSE)
  }
  reigen_special(.as_eigen(x), codes[[type]])
}

# ---------------------------------------------------------------- 矩阵函数

#' Matrix functions from Eigen's MatrixFunctions module
#'
#' These are **not** elementwise: `eigen_expm(A)` is the matrix exponential
#' `sum_k A^k / k!`, `eigen_logm(A)` the inverse of it, `eigen_sqrtm(A)` the
#' matrix square root and `eigen_powm(A, p)` the matrix power.  For the
#' elementwise versions use the `Math` group (`exp(A)` on an `eigenMat` already
#' works and is elementwise, like base R).
#'
#' All four require a square real matrix.  A negative real eigenvalue makes the
#' real logarithm/square root undefined; Eigen returns `NaN` silently in that
#' case, so the wrappers stop instead and say why.
#'
#' @param x A square real matrix.
#' @param p A single number: the exponent for `eigen_powm()`.
#' @return An `eigenMat`.
#' @name eigen-matrix-functions
NULL

#' @rdname eigen-matrix-functions
#' @export
eigen_expm <- function(x) reigen_matrix_function(.as_eigen(x), 0L, 0)

#' @rdname eigen-matrix-functions
#' @export
eigen_logm <- function(x) reigen_matrix_function(.as_eigen(x), 1L, 0)

#' @rdname eigen-matrix-functions
#' @export
eigen_sqrtm <- function(x) reigen_matrix_function(.as_eigen(x), 2L, 0)

#' @rdname eigen-matrix-functions
#' @export
eigen_powm <- function(x, p) reigen_matrix_function(.as_eigen(x), 3L, as.numeric(p)[1])
