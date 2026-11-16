# ===========================================================================
# matrix_theory.R —— 矩阵论教材常用的 R 层辅助算法
#
# 这些函数刻意不放进 C++：教材算法（行最简形、零空间、特征多项式、Jordan 块、
# 减号逆等）依赖主元策略或浮点阈值，数值上并没有唯一答案；Eigen 本身也刻意
# 不提供它们。这里用已导出的 C++ 能力（秩、SVD、特征值、分解）作地基，
# 在 R 层补出教学需要的入口，便于教材笔记本直接调用。
# ===========================================================================

.eigen_as_matrix <- function(x) as.matrix(x)

## ---------------------------------------------------------------- 行最简形
#' @export
eigen_rref <- function(x, tol = 1e-10) {
  A <- .eigen_as_matrix(x) * 1
  m <- nrow(A); n <- ncol(A); r <- 0L; piv <- integer(0)
  for (j in seq_len(n)) {
    if (r >= m) break
    col <- abs(A[(r + 1L):m, j])
    if (!length(col) || max(col) <= tol) next
    i <- (r + 1L):m; i <- i[which.max(col)]
    if (i != r + 1L) A[c(r + 1L, i), ] <- A[c(i, r + 1L), ]
    A[r + 1L, ] <- A[r + 1L, ] / A[r + 1L, j]
    for (k in seq_len(m)) if (k != r + 1L && abs(A[k, j]) > tol)
      A[k, ] <- A[k, ] - A[k, j] * A[r + 1L, ]
    r <- r + 1L; piv <- c(piv, j)
  }
  A[abs(A) < tol] <- 0
  list(rref = A, rank = r, pivots = piv)
}

## ---------------------------------------------------------------- 零空间 / 列空间
#' @export
eigen_nullspace <- function(x, tol = 1e-10) {
  R <- eigen_rref(x, tol)
  n <- ncol(.eigen_as_matrix(x))
  free <- setdiff(seq_len(n), R$pivots)
  if (!length(free)) return(matrix(0, n, 0L))
  N <- matrix(0, n, length(free))
  for (k in seq_along(free)) {
    N[free[k], k] <- 1
    for (i in seq_along(R$pivots)) N[R$pivots[i], k] <- -R$rref[i, free[k]]
  }
  N
}

#' @export
eigen_colspace <- function(x, tol = 1e-10) {
  A <- .eigen_as_matrix(x)
  piv <- eigen_rref(A, tol)$pivots
  if (!length(piv)) return(matrix(0, nrow(A), 0L))
  A[, piv, drop = FALSE]
}

#' @export
eigen_rowspace <- function(x, tol = 1e-10) {
  A <- .eigen_as_matrix(x)
  R <- eigen_rref(t(A), tol)
  piv <- R$pivots
  if (!length(piv)) return(matrix(0, ncol(A), 0L))
  A[, piv, drop = FALSE]
}

## ---------------------------------------------------------------- Gram / 正交化 / 投影
#' @export
eigen_gram <- function(x) {
  A <- .eigen_as_matrix(x)
  if (is.complex(A)) Conj(t(A)) %*% A else t(A) %*% A
}

#' @export
eigen_gram_schmidt <- function(x, tol = 1e-12) {
  A <- .eigen_as_matrix(x) * 1
  m <- nrow(A); n <- ncol(A)
  Q <- matrix(0, m, n); R <- matrix(0, n, n)
  for (j in seq_len(n)) {
    v <- A[, j]
    for (i in seq_len(j - 1L)) {
      R[i, j] <- if (is.complex(A)) sum(Conj(Q[, i]) * A[, j]) else sum(Q[, i] * A[, j])
      v <- v - R[i, j] * Q[, i]
    }
    R[j, j] <- sqrt(sum(Mod(v) ^ 2))
    if (!is.finite(R[j, j]) || R[j, j] < tol)
      stop("eigen_gram_schmidt: 列向量线性相关")
    Q[, j] <- v / R[j, j]
  }
  list(Q = Q, R = R)
}

#' @export
eigen_projection <- function(x, tol = 1e-10) {
  Q <- eigen_gram_schmidt(eigen_colspace(x, tol), tol)$Q
  Q %*% if (is.complex(Q)) Conj(t(Q)) else t(Q)
}

## ---------------------------------------------------------------- 顺序主子式 / 伴随 / 特征多项式
#' @export
eigen_leading_minors <- function(x) {
  A <- .eigen_as_matrix(x)
  vapply(seq_len(min(nrow(A), ncol(A))), function(k)
    det(A[seq_len(k), seq_len(k), drop = FALSE]), numeric(1))
}

#' @export
eigen_adjugate <- function(x) {
  A <- .eigen_as_matrix(x)
  if (nrow(A) != ncol(A)) stop("eigen_adjugate: 需要方阵")
  n <- nrow(A)
  C <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n))
    C[i, j] <- (-1) ^ (i + j) * det(A[-i, -j, drop = FALSE])
  Conj(t(C))
}

#' @export
eigen_charpoly <- function(x) {
  A <- .eigen_as_matrix(x) * 1
  if (nrow(A) != ncol(A)) stop("eigen_charpoly: 需要方阵")
  n <- nrow(A); c0 <- 1; B <- matrix(0, n, n)
  for (k in seq_len(n)) {
    B <- A %*% B + c0[k] * diag(n)
    c0[k + 1L] <- -sum(diag(A %*% B)) / k
  }
  c0
}

#' @export
eigen_companion <- function(coeffs) {
  a <- as.vector(coeffs)
  n <- length(a) - 1L
  if (n < 1L) stop("eigen_companion: 至少需要两个系数")
  C <- matrix(0, n, n)
  if (n > 1L) C[2:n, 1:(n - 1L)] <- diag(n - 1L)
  C[, n] <- -a[-1L] / a[1L]
  C
}

## ---------------------------------------------------------------- 不换行分解
#' @export
eigen_doolittle <- function(x, tol = 1e-12) {
  A <- .eigen_as_matrix(x) * 1; n <- nrow(A)
  if (n != ncol(A)) stop("eigen_doolittle: 需要方阵")
  L <- diag(1, n); U <- matrix(0, n, n)
  for (k in seq_len(n)) {
    prev <- seq_len(k - 1L)
    pivot <- A[k, k] - if (k > 1L) sum(L[k, prev] * U[prev, k]) else 0
    if (abs(pivot) < tol)
      stop(sprintf("eigen_doolittle: 第 %d 个主元为 0，不换行分解失败", k))
    for (j in seq.int(k, n))
      U[k, j] <- A[k, j] - if (k > 1L) sum(L[k, prev] * U[prev, j]) else 0
    if (k < n) for (i in seq.int(k + 1L, n))
      L[i, k] <- (A[i, k] - if (k > 1L) sum(L[i, prev] * U[prev, k]) else 0) / U[k, k]
  }
  list(L = L, U = U)
}

#' @export
eigen_ldlt_simple <- function(x, tol = 1e-12) {
  A <- .eigen_as_matrix(x)
  d <- eigen_doolittle(A, tol)
  list(L = d$L, D = diag(diag(d$U)))
}

## ---------------------------------------------------------------- Jordan 块型
#' @export
eigen_jordan <- function(x, tol = 1e-8) {
  A <- .eigen_as_matrix(x) * 1
  n <- nrow(A)
  if (n != ncol(A)) stop("eigen_jordan: 需要方阵")
  vals <- as.vector(eigen_gvalues(eigenMat(A)))
  ## Eigen 的复特征值在重根处会有 1e-8 级抖动，分组容度放宽到 1e-6
  gtol <- max(tol, 1e-6)
  key <- paste0(round(Re(vals) / gtol), ":", round(Im(vals) / gtol))
  grp <- split(seq_along(vals), key)
  J <- matrix(0, n, n); blocks <- list(); pos <- 1L
  for (g in grp) {
    lam <- vals[g[1L]]
    ## 数值特征值带微小虚部/误差；构造幂零矩阵前先按容差取代表元
    lam0 <- as.complex(complex(real = round(Re(lam) / gtol) * gtol,
                               imaginary = round(Im(lam) / gtol) * gtol))
    m <- length(g)
    M <- A - lam0 * diag(n)
    d <- 0L                      # d[k] = dim ker(M^k)，先放 d_0 = 0
    P <- diag(n)
    repeat {
      P <- P %*% M
      d <- c(d, n - eigen_rank(eigenMat(P), tol = max(tol, 1e-7)))
      k <- length(d) - 1L
      if (d[k + 1L] >= m || k >= n) break
    }
    sizes <- integer(0)
    for (k in seq_len(length(d) - 1L)) {
      inc1 <- d[k + 1L] - d[k]
      inc2 <- if (k < length(d) - 1L) d[k + 2L] - d[k + 1L] else 0L
      if (inc1 > inc2) sizes <- c(sizes, rep.int(k, inc1 - inc2))
    }
    if (!length(sizes) || sum(sizes) != m) sizes <- m   # 兜底：单块
    sizes <- sort(sizes, decreasing = TRUE)
    for (sz in sizes) {
      J[pos:(pos + sz - 1L), pos:(pos + sz - 1L)] <- lam
      if (sz > 1L) for (t in seq_len(sz - 1L)) J[pos + t - 1L, pos + t] <- 1
      blocks[[length(blocks) + 1L]] <- list(value = lam, size = sz)
      pos <- pos + sz
    }
  }
  list(values = vals, blocks = blocks, J = J)
}

## ---------------------------------------------------------------- 特殊积 / 向量化
#' @export
eigen_kronecker <- function(a, b) eigenMat(kronecker(.eigen_as_matrix(a), .eigen_as_matrix(b)))

#' @export
eigen_hadamard <- function(a, b) eigenMat(.eigen_as_matrix(a) * .eigen_as_matrix(b))

#' @export
eigen_vec <- function(x, byrow = FALSE) {
  A <- .eigen_as_matrix(x)
  eigenVec(as.vector(if (byrow) t(A) else A))
}

## ---------------------------------------------------------------- 二次型 / 谱量
#' @export
eigen_inertia <- function(x, tol = 1e-10) {
  A <- .eigen_as_matrix(x)
  lam <- if (is.complex(A)) as.vector(eigen_gvalues(eigenMat(A))) else
    as.vector(eigen_eigen(eigenMat(A))$values())
  lam <- Re(lam)
  c(positive = sum(lam > tol), negative = sum(lam < -tol),
    zero = sum(abs(lam) <= tol))
}

#' @export
eigen_spectral_radius <- function(x) {
  max(Mod(as.vector(eigen_gvalues(eigenMat(x)))))
}

#' @export
eigen_cond <- function(x, p = 2) {
  if (p == 2) return(1 / eigen_rcond(eigenMat(x)))
  A <- .eigen_as_matrix(x)
  as.numeric(eigen_norm(eigenMat(A), p) * eigen_norm(eigenMat(solve(A)), p))
}

#' @export
eigen_polar <- function(x) {
  A <- .eigen_as_matrix(x)
  s <- eigen_svd(eigenMat(A))
  U <- as.matrix(s$u()); V <- as.matrix(s$v())
  eigenMat(U %*% if (is.complex(V)) Conj(t(V)) else t(V))
}

## ---------------------------------------------------------------- 满秩分解 / 减号逆
#' @export
eigen_full_rank <- function(x, tol = 1e-10) {
  A <- .eigen_as_matrix(x)
  R <- eigen_rref(A, tol)
  if (!R$rank) return(list(B = matrix(0, nrow(A), 0L), C = matrix(0, 0L, ncol(A))))
  B <- A[, R$pivots, drop = FALSE]
  C <- R$rref[seq_len(R$rank), , drop = FALSE]
  list(B = B, C = C)
}

#' @export
eigen_ginverse <- function(x, tol = 1e-10) {
  A <- .eigen_as_matrix(x)
  fr <- eigen_full_rank(A, tol)
  B <- fr$B; C <- fr$C
  Cpinv <- Conj(t(C)) %*% solve(C %*% Conj(t(C)))
  Bpinv <- solve(Conj(t(B)) %*% B) %*% Conj(t(B))
  Cpinv %*% Bpinv
}
