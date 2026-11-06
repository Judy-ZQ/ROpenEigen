# 不在 base R 泛型覆盖范围内的操作
#
# 能挂到泛型上的都已经挂在 R/mat.R 里了（Ops / Math / Summary / %*% / t /
# solve / determinant / Re / Im / [ ）。这里只收泛型够不到的那些：迹、秩、
# 条件数、逆、拼接、块、叉积、近似比较。
#
# det / norm / diag / trace 都不是 S3 泛型（isS3stdGeneric 对内部泛型会返回
# FALSE，不能拿它判断），所以只能以 eigen_ 前缀的独立函数提供。

#' Trace of a matrix
#'
#' `eigen_logdet()` returns `log|det|` computed by summing `log|U_ii|` over the
#' LU factorisation, so it stays finite where `eigen_det()` would overflow.
#' The sign is not included — use `determinant()` if you need it.
#'
#' @param x An `eigenMat`, or an `eigenSparse`.
#' @return A length-1 numeric (or complex) vector.
#' @export
eigen_trace <- function(x) {
  # 迹只碰主对角线，稀疏侧在统计的那一遍里顺带累加出来，不需要稠密化。
  # 缺失值也要单独判：非对角线上的 NA 对迹没有影响，所以不能用全矩阵的计数，
  # 得用 trace_na / trace_nan（只数对角线上的）。
  if (inherits(x, "eigenSparse")) {
    st <- reigen_sparse_stats(x)
    if (st[["trace_na"]]  > 0) return(NA_real_)
    if (st[["trace_nan"]] > 0) return(NaN)
    return(st[["trace"]])
  }
  reigen_mat_trace(.as_eigen(x))
}

#' Numerical rank
#'
#' Uses a complete-pivot LU and counts pivots that exceed a threshold relative
#' to the largest pivot.
#'
#' The default threshold is `eps * min(nrow(x), ncol(x))`, i.e. about `1.3e-15`
#' for a 6x6 matrix. That is **far tighter than base R's** [base::qr()], which
#' uses `tol = 1e-7`. The two therefore disagree on ill-conditioned input, and
#' both are right for their own definition: a matrix whose smallest singular
#' value is `1e-12` times the largest is full rank under Eigen's threshold and
#' rank-deficient under R's. Pass `tol = 1e-7` to reproduce `qr(x)$rank`.
#'
#' @param x An `eigenMat`.
#' @param tol Threshold as a multiple of the largest pivot; `0` keeps Eigen's
#'   default (`eps * min(nrow, ncol)`).
#' @return A length-1 numeric vector.
#' @export
eigen_rank <- function(x, tol = 0) reigen_rank(.as_eigen(x), as.numeric(tol))

#' Reciprocal condition number
#'
#' `sigma_min / sigma_max` from the SVD -- the reciprocal of the 2-norm
#' condition number, so `1 / kappa(x, exact = TRUE)`. Defined for square and
#' rectangular matrices alike, and `0` for a singular one.
#'
#' This is deliberately not Eigen's `FullPivLU::rcond()`: that estimate is
#' only defined for square input (it asserts `rows() == cols()`, and R compiles
#' with `-DNDEBUG`, so the assert vanishes and a rectangular matrix reads out
#' of bounds), and on an exactly singular matrix it returns 1 where the answer
#' is 0.
#'
#' Every decomposition object's `$rcond()` method returns this same quantity --
#' `eigen_llt(x)$rcond()`, `eigen_ldlt(x)$rcond()` and
#' `eigen_full_lu(x)$rcond()` all agree with `eigen_rcond(x)`. `LLT` and `LDLT`
#' ship an `rcond()` in Eigen that computes something else entirely, so the
#' wrapper recomputes from the reconstructed matrix instead; one matrix must not
#' have two condition numbers depending on which door you came in through.
#'
#' @param x An `eigenMat`.
#' @return A length-1 numeric vector.
#' @export
eigen_rcond <- function(x) reigen_rcond(.as_eigen(x))

#' Matrix inverse
#' @param x An `eigenMat`.
#' @return An `eigenMat`.
#' @export
eigen_inverse <- function(x) reigen_inverse(.as_eigen(x))

#' Solve a linear system
#'
#' Square systems go through partial-pivot LU; rectangular ones go through
#' column-pivoting Householder QR and give the least-squares solution, matching
#' [qr.solve()]. Passing a decomposition object as `a` reuses its factorisation
#' instead of refactorising.
#'
#' @param a Coefficient matrix, or a decomposition object.
#' @param b Right-hand side.
#' @return An `eigenMat`.
#' @export
eigen_solve <- function(a, b) {
  # eigenMat 的类向量是 c("eigenMat", "eigenObject")——共用 $ 和 print 的代价。
  # 但它没有 $solve 方法，只判 "eigenObject" 会让 eigen_solve(E, b) 掉进
  # $.eigenObject，报一句 "eigenMat has no method 'solve'"。分解对象才是
  # c("eigenLLT", "eigenObject") 这种形态。
  if (inherits(a, "eigenObject") && !inherits(a, "eigenMat")) return(a$solve(b))
  reigen_solve(.as_eigen(a), .as_eigen(b))
}

#' Determinant
#'
#' `det()` is not an S3 generic in base R, so it cannot be hooked directly.
#' It does still work, because its definition is
#' `sign * exp(determinant(x)$modulus)` and `determinant()` *is* a generic that
#' dispatches on `eigenMat` — but it has to exponentiate, so it overflows to
#' `Inf` or underflows to `0` once the determinant leaves the double range.
#' `eigen_det()` returns the value without that detour; `eigen_logdet()` returns
#' `log|det|`, which is the safe choice for anything large.
#'
#' @param x An `eigenMat`.
#' @return A length-1 numeric, or complex for a complex matrix.
#' @export
eigen_det <- function(x) reigen_det(.as_eigen(x))

#' Matrix norms
#'
#' Mirrors [base::norm()], which is also not an S3 generic. `type` accepts
#' `"O"`/`"1"` (max absolute column sum), `"I"` (max absolute row sum),
#' `"F"` (Frobenius), `"M"` (max modulus) and `"2"` (spectral). `"p"` gives the
#' entrywise p-norm of the flattened matrix.
#'
#' For a sparse matrix the spectral norm is the one type that cannot be read off
#' a single pass: the largest singular value is the square root of the largest
#' eigenvalue of `A'A`, and a sparse matrix has no eigendecomposition (Eigen's
#' `SelfAdjointEigenSolver` takes dense input only). `eigen_norm()` therefore
#' runs a power iteration `z <- A' (A v)`, `v <- z/|z|`, using `|z|` as the
#' running estimate of the largest eigenvalue and stopping on a relative change
#' below `1e-12`. The iteration converges linearly at the ratio of the two
#' largest eigenvalues, so it is fast unless the two leading singular values are
#' nearly equal; if it has not converged after 10000 steps it stops rather than
#' returning a value that looks converged and is not.
#'
#' @param x An `eigenMat`, or an `eigenSparse`. Sparse input never densifies:
#'   `"F"` and `"M"` are one `O(nnz)` scan, `"O"` and `"I"` are `O(nnz)` bin sums
#'   over the stored entries, and `"2"` is the power iteration above. `"p"` has
#'   no sparse meaning and stops with an explicit error.
#' @param type Norm type.
#' @param p Exponent for `type = "p"`.
#' @return A length-1 numeric.
#' @export
eigen_norm <- function(x, type = c("O", "I", "F", "M", "2"), p = 2) {
  if (length(type) > 1L) type <- type[1L]
  if (identical(type, "1")) type <- "O"
  if (inherits(x, "eigenSparse")) return(.sparse_norm(x, type))
  code <- switch(type, O = 1L, I = 2L, F = 0L, M = 3L, `2` = 4L, p = 5L,
                 stop(sprintf("unsupported norm type '%s'", type), call. = FALSE))
  reigen_norm(.as_eigen(x), code, as.numeric(p))
}

# 稀疏矩阵的五种范数。全都是 O(nnz) 或 O(nnz × 迭代数)，没有一条会稠密化。
#
#   F  —— Frobenius，统计那一遍顺带算出来的 **缩放后** 平方和（st[["frob"]]）。
#         不能用 sqrt(st[["sumsq"]])：sumsq 是朴素的，对应 base 的 sum(x^2)，
#         在 1e-200 上下溢成 0、在 1e300 上溢成 Inf。范数要的是缩放版。
#   M  —— max|a_ij|。隐含零是 0，所以要把 0 拉进候选：全正矩阵的 M 是
#         max，不是 max(|min|)；全零矩阵的 M 是 0。
#   O  —— 列绝对值之和的最大值；I 是行的。两者共用 reigen_sparse_norm_oi()
#         的一遍扫描（base 的 tabulate() 没有 weights 参数，别试）。
#   "2" —— 谱范数，幂迭代，见 reigen_sparse_norm2()。
.sparse_norm <- function(x, type) {
  st <- reigen_sparse_stats(x)
  # base 的四个非谱范数对 NA 都是"传上去"：sum/max 遇 NA 得 NA。谱范数走
  # svd，base 在那里报错。分开处理，别让一个 NA 变成"没收敛"。
  bad <- .sparse_na(st)
  if (!is.null(bad)) {
    if (!identical(type, "2")) return(bad)
    stop("eigen_norm(type = \"2\") on an eigenSparse needs finite values; the matrix has missing entries",
         call. = FALSE)
  }
  if (identical(type, "F")) return(st[["frob"]])
  if (identical(type, "M")) return(max(abs(st[["min"]]), abs(st[["max"]])))
  if (type %in% c("O", "I")) {
    v <- reigen_sparse_norm_oi(x)
    return(if (type == "O") v[[1L]] else v[[2L]])
  }
  if (identical(type, "2")) {
    v <- reigen_sparse_norm2(x, 1e-12, 10000L)
    if (is.na(v))
      stop(paste0("eigen_norm(type = \"2\") did not converge in 10000 power ",
                  "iterations; the two largest singular values are too close. ",
                  "Use type = \"F\", or densify with as.matrix() and call ",
                  "eigen_norm() on the dense matrix."), call. = FALSE)
    return(v)
  }
  stop(sprintf(paste0("eigen_norm() on an eigenSparse does not support type = \"%s\". ",
                      "Supported: \"O\" (or \"1\"), \"I\", \"F\", \"M\", \"2\". ",
                      "type = \"p\" has no sparse meaning here."), type),
       call. = FALSE)
}

#' @rdname eigen_det
#' @export
eigen_is_complex <- function(x) {
  # 稀疏层只支持 double，所以答案恒为 FALSE——但仍然要分流，否则 .as_eigen()
  # 会直接把它拒掉，一个恒为 FALSE 的问题变成一次报错。
  if (inherits(x, "eigenSparse")) return(FALSE)
  reigen_mat_is_complex(.as_eigen(x))
}

#' @rdname eigen_trace
#' @export
eigen_logdet <- function(x) reigen_logdet(.as_eigen(x))

#' Inner and cross products
#'
#' `eigen_dot()` is `sum(a * b)` -- a plain elementwise product, summed. It
#' follows R's convention and deliberately does **not** conjugate the first
#' argument, which means it does **not** agree with Eigen's
#' `MatrixBase::dot()`: Eigen's own documentation says that "if the scalar type
#' is complex numbers, then this function returns the hermitian
#' (sesquilinear) dot product, conjugate-linear in the first variable".
#'
#' For a real vector the two conventions coincide, so the difference only shows
#' up in complex arithmetic. Use `sum(Conj(a) * b)` or `eigen_adjoint(a) %*% b`
#' when you want the Hermitian inner product.
#'
#' Consequently `eigen_dot(x, x)` is *not* the squared norm of a complex `x`:
#' it is `sum(x * x)`, which can be complex. Use [eigen_sumsq()] (or
#' [eigen_norm()]`(x, "F")^2`) for the sum of squared moduli, `sum(abs(x)^2)`.
#'
#' @param a,b `eigenMat` or numeric.
#' @return `eigen_dot()` returns a length-1 numeric or complex vector;
#'   `eigen_cross()` returns an `eigenMat` for two length-3 vectors, and a plain
#'   length-1 numeric for two length-2 vectors -- the 2-D cross product is the
#'   scalar z component, as `Eigen::Vector2d::cross()` returns it.
#' @name eigen-products
NULL

#' @rdname eigen-products
#' @export
eigen_dot <- function(a, b) as.vector(reigen_dot(.as_eigen(a), .as_eigen(b)))

#' @rdname eigen-products
#' @export
eigen_cross <- function(a, b) reigen_cross(.as_eigen(a), .as_eigen(b))

#' Outer product
#' @param a,b `eigenMat` or numeric.
#' @return An `eigenMat` of size `length(a)` by `length(b)`.
#' @export
eigen_outer <- function(a, b) reigen_matmul(.as_eigen(a), t(.as_eigen(b)))

#' Conjugate transpose
#' @param x An `eigenMat`, or an `eigenSparse`.
#' @return An `eigenMat`, or an `eigenSparse`.
#' @export
eigen_adjoint <- function(x) {
  # 稀疏层只有 double，共轭转置就是转置，而且结果仍然稀疏。
  if (inherits(x, "eigenSparse")) return(t(x))
  reigen_mat_adjoint(.as_eigen(x))
}

#' Bind matrices
#' @param ... Matrices or vectors.
#' @param along `"row"` for `rbind`, `"col"` for `cbind`.
#' @return An `eigenMat`.
#' @export
eigen_bind <- function(..., along = c("row", "col")) {
  along <- match.arg(along)
  xs <- lapply(list(...), .as_eigen)
  if (length(xs) == 0L) stop("nothing to bind", call. = FALSE)
  side <- if (along == "row") 1L else 2L
  Reduce(function(p, q) reigen_mat_bind(p, q, side), xs)
}

#' @rdname eigen_bind
#' @export
eigen_rbind <- function(...) eigen_bind(..., along = "row")

#' @rdname eigen_bind
#' @export
eigen_cbind <- function(...) eigen_bind(..., along = "col")

#' Tile a matrix
#' @param x An `eigenMat`.
#' @param nrow,ncol Repetition counts.
#' @return An `eigenMat`.
#' @export
eigen_rep <- function(x, nrow = 1L, ncol = 1L)
  reigen_mat_rep(.as_eigen(x), as.integer(nrow), as.integer(ncol))

#' Reshape a matrix
#' @param x An `eigenMat`.
#' @param nrow,ncol New dimensions; the product must match the element count.
#' @return An `eigenMat`.
#' @export
eigen_reshape <- function(x, nrow, ncol)
  reigen_mat_reshape(.as_eigen(x), as.integer(nrow), as.integer(ncol))

#' Extract a submatrix
#' @param x An `eigenMat`.
#' @param row,col Integer ranges (1-based, inclusive). `NULL` means "to the end".
#' @return An `eigenMat`.
#' @export
eigen_block <- function(x, row = NULL, col = NULL) {
  d <- dim(x)
  r0 <- if (is.null(row)) 0L else as.integer(min(row))
  r1 <- if (is.null(row)) 0L else as.integer(max(row))
  c0 <- if (is.null(col)) 0L else as.integer(min(col))
  c1 <- if (is.null(col)) 0L else as.integer(max(col))
  reigen_mat_block(x, r0, r1, c0, c1)
}

#' Approximate equality
#' @param a,b `eigenMat` or numeric.
#' @param tol Relative tolerance.
#' @return A length-1 logical.
#' @export
eigen_approx <- function(a, b, tol = 1e-10)
  reigen_approx(.as_eigen(a), .as_eigen(b), as.numeric(tol))

#' Finiteness checks
#'
#' `is.finite()` works elementwise and returns a logical matrix; these three
#' return a single logical and stop at the first offending entry, so they remain
#' the right tool for a precondition guard.
#' @param x An `eigenMat`.
#' @return A length-1 logical.
#' @name eigen-finite
NULL

#' @rdname eigen-finite
#' @export
eigen_all_finite <- function(x) .eigen_finite(x, 0L)

#' @rdname eigen-finite
#' @export
eigen_has_nan <- function(x) .eigen_finite(x, 1L)

#' @rdname eigen-finite
#' @export
eigen_has_inf <- function(x) .eigen_finite(x, 2L)

# 稀疏对象不能走 .as_eigen()——它不是 eigenMat。这里显式分流，而不是把
# eigen_all_finite 变成泛型：只有两个分支，泛型带来的间接层不划算。
.eigen_finite <- function(x, code) {
  if (inherits(x, "eigenSparse")) return(reigen_sparse_finite(x, code))
  reigen_finite(.as_eigen(x), code)
}

#' Position of the minimum or maximum
#'
#' Returns the position of the extreme value in **column-major** order, which is
#' the order [base::which.min()] uses, so the two agree — including on ties
#' (the first occurrence wins) and on missing values (`NA` and `NaN` are
#' skipped; `Inf` is a value, not a missing one).
#'
#' For an `eigenSparse` the implicit zeros are part of the candidate set: the
#' minimum of a matrix that is mostly zeros is `0`, and its position is the
#' first column-major cell whose value is `0`, whether that cell is stored or
#' not. That is exactly what `which.min(as.matrix(x))` returns. If every element
#' is missing, `value`, `row` and `col` all come back `NA`.
#'
#' @param x An `eigenMat`, or an `eigenSparse`.
#' @return A list with `value`, `row` and `col` (1-based).
#' @export
eigen_which_min <- function(x) {
  if (inherits(x, "eigenSparse")) return(reigen_sparse_which(x, 0L))
  reigen_which(.as_eigen(x), 0L)
}

#' @rdname eigen_which_min
#' @export
eigen_which_max <- function(x) {
  if (inherits(x, "eigenSparse")) return(reigen_sparse_which(x, 1L))
  reigen_which(.as_eigen(x), 1L)
}

# reigen_sparse_stats() 是给 na.rm = TRUE 那条路用的：它**跳过**缺失值累加，
# 同时把 NA 和 NaN 分开数。这个设计本身是对的，但它意味着每个调用方都必须
# 自己决定"要不要把缺失值传上去"。
#
# 漏掉这一步不会报错，只会给出一个看起来合理的偏小值。实测过：
#   A = matrix(c(1, 2, NA, 4, 5, 6), 2, 3)
#   eigen_sum(S)   -> 18   （sum(S) 是 NA）
#   eigen_prod(S)  -> 240  （prod(S) 是 NA）
#   mean(S)        -> 3    （mean(A) 是 NA）
# 同一个量两扇门给出不同答案，正是 consistency.R 存在的理由。
#
# base 区分 NA 和 NaN：sum(c(1, NA)) 是 NA，sum(c(1, NaN)) 是 NaN。两者 is.na()
# 都为真但 identical() 为假，所以不能一律返回 NA_real_。NA 优先——和 base 一致，
# 混合时 sum(c(NA, NaN)) 是 NA。
.sparse_na <- function(st) {
  if (st[["na"]]  > 0) return(NA_real_)
  if (st[["nan"]] > 0) return(NaN)
  NULL
}

# 取一个统计量，缺失值存在且没被要求跳过时先把它传上去。
.sparse_stat <- function(x, field) {
  st <- reigen_sparse_stats(x)
  bad <- .sparse_na(st)
  if (!is.null(bad)) return(bad)
  st[[field]]
}

#' Reductions that have no base R generic
#'
#' These have no `na.rm` argument, matching the behaviour of the plain
#' `sum()`/`prod()`/`mean()` they mirror: a missing value propagates, so
#' `eigen_sum()` of a matrix containing `NA` is `NA` and one containing `NaN` is
#' `NaN`. `NA` wins when both are present, as it does in base R.
#'
#' Use `sum(x, na.rm = TRUE)` on an `eigenSparse` when you want the missing
#' values dropped; that path is what the underlying single pass is built for,
#' and it never densifies.
#'
#' @param x An `eigenMat`, or an `eigenSparse` (computed from the stored
#'   entries in `O(nnz)`, with implicit zeros taking part).
#' @return A length-1 numeric.
#' @export
eigen_sum <- function(x) {
  if (inherits(x, "eigenSparse")) return(.sparse_stat(x, "sum"))
  reigen_reduce(.as_eigen(x), 0L)
}

#' @rdname eigen_sum
#' @export
eigen_prod <- function(x) {
  if (inherits(x, "eigenSparse")) return(.sparse_stat(x, "prod"))
  reigen_reduce(.as_eigen(x), 1L)
}

#' @rdname eigen_sum
#' @export
eigen_mean <- function(x) {
  if (inherits(x, "eigenSparse")) {
    st <- reigen_sparse_stats(x)
    bad <- .sparse_na(st)
    if (!is.null(bad)) return(bad)
    return(st[["sum"]] / st[["total"]])
  }
  reigen_reduce(.as_eigen(x), 2L)
}

#' @rdname eigen_sum
#' @export
eigen_sumsq <- function(x) {
  if (inherits(x, "eigenSparse")) return(.sparse_stat(x, "sumsq"))
  reigen_reduce(.as_eigen(x), 5L)
}
