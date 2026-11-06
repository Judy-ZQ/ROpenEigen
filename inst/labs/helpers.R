# ===========================================================================
# helpers.R —— 实验教程公共工具
#
# 这里放的是「教材要求手算、但 ROpenEigen 出于数值原因刻意不提供」的那几个
# 算法：行最简形、零空间、代数余子式、特征多项式、Gram-Schmidt、逆序数。
# 每个都写成十几行的纯 R，配 ROpenEigen 的 det / rank / dot 做交叉验证。
#
# 为什么不放进包里：
#   * 行最简形依赖主元选取策略，浮点下没有唯一答案，工程上应当用 QR/SVD；
#   * 特征多项式展开在数值上是病态的（Wilkinson 多项式），工程上直接求特征值；
#   * 零空间用 SVD 的右奇异向量比 rref 稳。
# 实验里保留手写版本，是为了让「教材算法」和「数值算法」的差别看得见。
#
# 用法：source("helpers.R")，或 source(system.file("labs","helpers.R",package="ROpenEigen"))
# ===========================================================================

# ---------------------------------------------------------------- 1. 排列与逆序数
# 教材 §1.1 用逆序数定义 n 阶行列式。逆序数 = 前面比它大的数的个数之和。
inversions <- function(p) {
  p <- as.integer(p)
  sum(vapply(seq_along(p), function(i) sum(p[seq_len(i - 1)] > p[i]), integer(1)))
}

# n! 个排列全枚举（n <= 8 才现实），每行一个排列
permutations <- function(n) {
  if (n == 1L) return(matrix(1L, 1, 1))
  sub <- permutations(n - 1L)
  out <- matrix(0L, nrow(sub) * n, n)
  k <- 0L
  for (pos in seq_len(n)) {                 # 把 n 插到第 pos 个位置
    for (r in seq_len(nrow(sub))) {
      k <- k + 1L
      out[k, ] <- append(sub[r, ], n, after = pos - 1L)
    }
  }
  out
}

# 按定义算行列式：Σ (-1)^τ(p) a_{1,p1} a_{2,p2} ... a_{n,pn}
# 只用 base R，用来和 eigen_det() 对账（n <= 7）
det_by_definition <- function(A) {
  n <- nrow(A); P <- permutations(n)
  s <- ifelse(apply(P, 1, inversions) %% 2 == 0, 1, -1)
  sum(s * apply(P, 1, function(p) prod(A[cbind(seq_len(n), p)])))
}

# ---------------------------------------------------------------- 2. 余子式与伴随矩阵
# 教材 §1.4：M_ij 是划去第 i 行第 j 列剩下的 n-1 阶行列式，
# A_ij = (-1)^(i+j) M_ij 是代数余子式。
cofactor_matrix <- function(A) {
  n <- nrow(A)
  C <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    C[i, j] <- (-1)^(i + j) * det(A[-i, -j, drop = FALSE])
  }
  C
}

# 伴随矩阵 adj(A) = (A_ij)^T，满足 A %*% adj(A) = det(A) * I
adjugate <- function(A) t(cofactor_matrix(A))

# ---------------------------------------------------------------- 3. 行最简形
# 教材 §2.5 的初等行变换。返回 list(rref=, rank=, pivots=, ops=)。
# 主元选取用「当前列绝对值最大」（部分主元），这是数值上唯一站得住的做法；
# 教材里手算时用「取第一个非零」，两者在整数矩阵上通常同解。
rref <- function(A, tol = 1e-10) {
  A <- as.matrix(A); A <- matrix(as.numeric(A), nrow(A), ncol(A))
  m <- nrow(A); n <- ncol(A); r <- 0L; piv <- integer(0)
  for (j in seq_len(n)) {
    if (r >= m) break
    col <- abs(A[(r + 1L):m, j])
    if (max(col) <= tol) next                    # 整列为零，跳过
    i <- (r + 1L):m; i <- i[which.max(col)]
    if (i != r + 1L) A[c(r + 1L, i), ] <- A[c(i, r + 1L), ]   # 换行
    A[r + 1L, ] <- A[r + 1L, ] / A[r + 1L, j]                  # 主元化为 1
    for (k in seq_len(m)) {
      if (k != r + 1L && abs(A[k, j]) > tol) {
        A[k, ] <- A[k, ] - A[k, j] * A[r + 1L, ]               # 其余行消成 0
      }
    }
    r <- r + 1L; piv <- c(piv, j)
  }
  A[abs(A) < tol] <- 0
  list(rref = A, rank = r, pivots = piv)
}

# ---------------------------------------------------------------- 4. 零空间（基础解系）
# 教材 §4.2：Ax = 0 的基础解系含 n - r 个向量。这里直接由 rref 读出来：
# 每个自由变量取 1、其余自由变量取 0，回代得到主元变量的值。
nullspace <- function(A, tol = 1e-10) {
  R <- rref(A, tol); n <- ncol(as.matrix(A)); piv <- R$pivots
  free <- setdiff(seq_len(n), piv)
  if (!length(free)) return(matrix(0, n, 0))
  N <- matrix(0, n, length(free))
  for (k in seq_along(free)) {
    N[free[k], k] <- 1
    for (i in seq_along(piv)) N[piv[i], k] <- -R$rref[i, free[k]]
  }
  N
}

# 非齐次 Ax = b 的一个特解（无解返回 NULL）
particular_solution <- function(A, b, tol = 1e-10) {
  Ab <- cbind(A, b)
  R <- rref(Ab, tol)
  n <- ncol(A)
  if (R$rank < rref(A, tol)$rank) return(NULL)      # rank(A) < rank(A|b)：无解
  piv <- R$pivots[R$pivots <= n]
  x <- numeric(n)
  x[piv] <- R$rref[seq_along(piv), n + 1]
  x
}

# ---------------------------------------------------------------- 5. 特征多项式
# 教材 §5.1：|λI - A| 展开成 λ 的多项式。Faddeev-LeVerrier 递推：
#   c_0 = 1,  B_0 = 0
#   c_k = -tr(A B_{k-1}) / k,  B_k = A B_{k-1} + c_k I
# 返回降幂系数 c(1, c1, ..., cn)，即 λ^n + c1 λ^(n-1) + ... + cn。
# 数值上极不稳定（教科书里的经典反例见 lab05），只适合小矩阵与教学演示。
charpoly_faddeev <- function(A) {
  A <- as.matrix(A); n <- nrow(A)
  ck <- numeric(n + 1); ck[1] <- 1
  B <- matrix(0, n, n)
  for (k in seq_len(n)) {
    B <- A %*% B + ck[k] * diag(n)
    ck[k + 1] <- -sum(diag(A %*% B)) / k
  }
  ck
}

# 由根构造多项式系数（升幂）：prod (x - r_i)
poly_from_roots <- function(r) {
  p <- 1
  for (ri in r) p <- c(0, p) - c(ri * p, 0)
  p
}

# ---------------------------------------------------------------- 6. Gram-Schmidt
# 教材 §3.4：把线性无关组化为正交组再单位化。返回 list(Q=正交单位列, R=上三角)。
# 与 QR 分解的关系：A = Q R，Q 的列就是 Schmidt 的结果（可能差一个符号）。
gram_schmidt <- function(A, tol = 1e-12) {
  A <- as.matrix(A); n <- ncol(A); m <- nrow(A)
  Q <- matrix(0, m, n); R <- matrix(0, n, n)
  for (j in seq_len(n)) {
    v <- A[, j]
    for (i in seq_len(j - 1L)) {
      R[i, j] <- sum(Q[, i] * A[, j])
      v <- v - R[i, j] * Q[, i]
    }
    R[j, j] <- sqrt(sum(v^2))
    if (R[j, j] < tol) stop("列向量线性相关，Gram-Schmidt 失效")
    Q[, j] <- v / R[j, j]
  }
  list(Q = Q, R = R)
}

# 带度量的 Schmidt 正交化：内积由 Gram 矩阵给出，<u,v> = u' G v。
# 处理「多项式空间」这类抽象内积空间时用这个：把基向量写成坐标，
# 内积不再是坐标的点积，而是坐标乘 Gram 矩阵。
gram_schmidt_metric <- function(A, G, tol = 1e-12) {
  A <- as.matrix(A); n <- ncol(A)
  Q <- matrix(0, nrow(A), n)
  ip <- function(u, v) as.numeric(sum(u * (G %*% v)))
  for (j in seq_len(n)) {
    v <- A[, j]
    for (i in seq_len(j - 1L)) {
      v <- v - ip(v, Q[, i]) / ip(Q[, i], Q[, i]) * Q[, i]
    }
    nv <- sqrt(ip(v, v))
    if (nv < tol) stop("向量组在这组内积下线性相关")
    Q[, j] <- v / nv
  }
  Q
}

# 去掉零空间的投影矩阵（到 A 的列空间的正交投影）
projection_matrix <- function(A) {
  Q <- qr.Q(qr(A)); Q %*% t(Q)
}

# ---------------------------------------------------------------- 7. 二次型
# 由二次型 f = x'Ax 的「上三角写法」得到对称矩阵：
# 教材 §6.1 的例题给的是 a_ij + a_ji 的合并形式，这里收 (i,j,系数) 三元组。
quad_matrix <- function(n, terms) {
  A <- matrix(0, n, n)
  for (t in terms) {                       # t = c(i, j, 系数)
    i <- t[1]; j <- t[2]; v <- t[3]
    if (i == j) A[i, j] <- A[i, j] + v
    else { A[i, j] <- A[i, j] + v / 2; A[j, i] <- A[j, i] + v / 2 }
  }
  A
}

# 代入求值：f(x) = x' A x
quad_eval <- function(A, x) as.numeric(t(x) %*% A %*% x)

# 顺序主子式（教材 §6.3 定理 6.3.2 的判别法）
leading_minors <- function(A) {
  n <- nrow(A)
  vapply(seq_len(n), function(k) det(A[seq_len(k), seq_len(k), drop = FALSE]), numeric(1))
}

# 惯性指数：正、负、零特征值个数（教材 §6.3 定理 6.3.3）
inertia <- function(A, tol = 1e-10) {
  lam <- Re(eigen(A, symmetric = TRUE, only.values = TRUE)$values)
  c(positive = sum(lam > tol), negative = sum(lam < -tol), zero = sum(abs(lam) <= tol))
}

# ---------------------------------------------------------------- 8. 先判条件，再求解
# 这一条是本教程反复强调的：ROpenEigen 的 solve() 把工作交给 Eigen 的 LU 分解，
# 而 LU 分解不做奇异检测——对奇异矩阵它照样算下去，返回一张 NaN/Inf 的矩阵，
# 而 base R 的 solve() 会直接报 "system is computationally singular"。
# 所以教学代码里求逆/解方程统一走下面两个入口：先看 rcond，再决定要不要算。
safe_solve <- function(A, b = NULL, tol = 1e-12) {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) stop("safe_solve: 需要方阵", call. = FALSE)
  rc <- eigen_rcond(eigenMat(A))
  if (!is.finite(rc) || rc < tol) {
    stop(sprintf("矩阵接近奇异（rcond = %.3e < %.0e），解不可信；请改用 SVD 或最小二乘",
                 rc, tol), call. = FALSE)
  }
  if (is.null(b)) solve(eigenMat(A)) else solve(eigenMat(A), eigenVec(b))
}

safe_inverse <- function(A, tol = 1e-12) safe_solve(A, NULL, tol)

# ---------------------------------------------------------------- 9. 输出小工具
# 让实验输出整齐：sprintf 风格的 cat，末尾换行
say <- function(...) cat(..., "\n", sep = "")

# 误差报告：相对误差，避免 1e-16 这种绝对量级看不出好坏
relerr <- function(a, b) {
  a <- as.matrix(a); b <- as.matrix(b)
  max(abs(a - b)) / max(1, max(abs(b)))
}

# 断言：不满足就停，用来把「实验结论」变成可执行的自检
expect <- function(ok, what) {
  if (!isTRUE(ok)) stop("实验失败：", what, call. = FALSE)
  invisible(TRUE)
}

expect_near <- function(a, b, tol = 1e-10, what = "数值不符") {
  e <- relerr(a, b)
  if (!is.finite(e) || e > tol) {
    stop(sprintf("%s（相对误差 %.3e > %.0e）", what, e, tol), call. = FALSE)
  }
  invisible(TRUE)
}
