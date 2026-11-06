# ===========================================================================
# 08 · 稀疏矩阵
#
# 对应 EigenExample：
#   special_examples/Tutorial_sparse_example.cpp
#   special_examples/Tutorial_sparse_example_details.cpp
#   snippets/SparseMatrix_coeffs.cpp
#   snippets/BiCGSTAB_simple.cpp
#   snippets/BiCGSTAB_step_by_step.cpp
#
# 运行：Rscript 08-sparse.R
#
# 对照要点：Eigen 的稀疏矩阵是**列主序压缩**结构，装配靠 Triplet 列表；
# ROpenEigen 的 eigen_sparse() 直接收 (i, j, x) 三元组，所以 cpp 里
# "buildProblem 往 vector<Triplet> 里塞、再 setFromTriplets" 的那一套，
# 在 R 里就是"攒好三个向量、一次调用"。
# ===========================================================================

suppressPackageStartupMessages(library(ROpenEigen))

# ---- SparseMatrix_coeffs.cpp ------------------------------------------------
# cpp 用 insert() 逐项插入，再 makeCompressed()，然后看 nonZeros() 与 coeffs()。
# R 里一次性给三元组即可；重复下标会相加（Eigen 的 setFromTriplets 语义）。
# 注意 cpp 里 A.insert(1,2)=0 那个 0 是**显式存下来**的，所以 nonZeros 是 3；
# R 侧同样保留显式零。
i <- c(2, 1, 3)          # cpp 是 0-based：(1,2) (0,1) (2,0)
j <- c(3, 2, 1)
x <- c(0, 1, 2)
A <- eigen_sparse(i, j, x, nrow = 3, ncol = 3)
print(as.matrix(A))
cat("stored nonzeros:", eigen_sparse_nnz(A), " 值:", eigen_triplets(A)$x, "\n")

# ---- Tutorial_sparse_example_details.cpp ------------------------------------
# cpp 里 buildProblem() 的细节：每个未知量连到自己和上下左右四个邻居，
# 边界上是 -1、中心是 4，右端项在边界中间取 1。下面把它照搬到 R。
# （cpp 用 n=300 的图，这里用 n=6 便于打印；改 n 即可放大。）
build_problem <- function(n) {
  m <- n * n
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  b <- numeric(m)
  idx <- function(r, c) (r - 1) * n + c          # 1-based 的 (r, c) → 线性下标
  for (c in 1:n) for (r in 1:n) {
    k <- idx(r, c)
    ii <- c(ii, k); jj <- c(jj, k); xx <- c(xx, 4)
    if (r > 1) { ii <- c(ii, k); jj <- c(jj, idx(r - 1, c)); xx <- c(xx, -1) }
    if (r < n) { ii <- c(ii, k); jj <- c(jj, idx(r + 1, c)); xx <- c(xx, -1) }
    if (c > 1) { ii <- c(ii, k); jj <- c(jj, idx(r, c - 1)); xx <- c(xx, -1) }
    if (c < n) { ii <- c(ii, k); jj <- c(jj, idx(r, c + 1)); xx <- c(xx, -1) }
    if (c == 1 || c == n || r == 1 || r == n) b[k] <- 0
  }
  # 与 cpp 相同的边界条件：图像左右两列中点附近取 1
  b[idx(round(n / 2), 1)] <- 1
  b[idx(round(n / 2), n)] <- 1
  list(A = eigen_sparse(ii, jj, xx, nrow = m, ncol = m), b = b, n = n)
}

n <- 6
P <- build_problem(n)
cat("未知量:", n * n, " 非零元:", eigen_sparse_nnz(P$A), "\n")

# ---- Tutorial_sparse_example.cpp --------------------------------------------
# cpp:  SimplicialCholesky<SpMat> chol(A);  x = chol.solve(b);
# R:    eigen_sparse_solve(A, b, "llt") —— 矩阵对称正定，LLT 最快。
x <- eigen_sparse_solve(P$A, P$b, "llt")
cat("残差:", max(abs(as.vector(as.matrix(P$A %*% x)) - P$b)), "\n")
cat("解的前 6 个:", round(as.vector(x)[1:6], 4), "\n")

# 换成迭代解法，接口一致，只是多两个属性看收敛情况：
xcg <- eigen_sparse_solve(P$A, P$b, "cg", tol = 1e-12, maxiter = 2000)
cat("CG 迭代次数:", attr(xcg, "iterations"),
    " 相对残差:", attr(xcg, "residual"), "\n")

# 把解铺回 n x n 的网格，就是 cpp 里 saveAsBitmap() 画的那张图：
print(matrix(round(as.vector(x), 3), n, n, byrow = TRUE))

# ---- BiCGSTAB_simple.cpp ----------------------------------------------------
# cpp 是伪代码（A、b 用注释省略了），重点是 solver.compute(A) 之后的
# iterations()/error()。R 里这两个数就是返回值的属性。
# 下面给一个真实的**非对称**稀疏系统：一维对流扩散，下对角 -1、上对角 -3、
# 对角 4。对称的 Poisson 问题用 CG 更合适，BiCGSTAB 的舞台正是这种非对称阵。
n <- 200
Acd <- eigen_sparse(
  i = c(1:n, 2:n, 1:(n - 1)),
  j = c(1:n, 1:(n - 1), 2:n),
  x = c(rep(4, n), rep(-1, n - 1), rep(-3, n - 1)),
  nrow = n, ncol = n)
cat("对称?", isSymmetric(as.matrix(Acd)), "\n")
bcd <- as.numeric(1:n)
x1 <- eigen_sparse_solve(Acd, bcd, "bicgstab")
cat("BiCGSTAB 迭代:", attr(x1, "iterations"),
    " 残差:", max(abs(as.vector(as.matrix(Acd %*% x1)) - bcd)), "\n")

# 同一个分解再解一次（cpp 注释里的 "update b; solve again"）：
bcd2 <- as.numeric(n:1)
x2 <- eigen_sparse_solve(Acd, bcd2, "bicgstab")
cat("第二次解残差:", max(abs(as.vector(as.matrix(Acd %*% x2)) - bcd2)), "\n")

# 上面 BiCGSTAB 只迭代 1 步就收敛，是因为对角占优 + 对角预条件把它变成了
# 接近单位矩阵的问题。要看清"迭代"这件事，换成 2-D Poisson（对称正定）用 CG：
cat("--- CG 在 2-D Poisson 上的收敛 ---", "\n")
for (mi in c(1, 2, 4, 10, 30)) {
  r <- tryCatch(eigen_sparse_solve(P$A, P$b, "cg", tol = 1e-10, maxiter = mi),
                error = function(e) NULL)
  if (is.null(r)) cat(sprintf("maxiter=%-4d 未收敛（按设计报错）\n", mi))
  else cat(sprintf("maxiter=%-4d 迭代 %d 次后收敛，残差 %.3e\n", mi,
                   attr(r, "iterations"),
                   max(abs(as.vector(as.matrix(P$A %*% r)) - P$b))))
}

# ---- 稀疏与稠密的互操作 -----------------------------------------------------
# as.matrix() 会把稀疏矩阵摊成稠密（小心大矩阵）；%*%、t()、solve() 都派发。
S <- eigen_sparse(c(1, 1, 2, 3), c(1, 2, 2, 3), c(2, 1, 3, 4), nrow = 3, ncol = 3)
cat("nnz:", eigen_sparse_nnz(S), "\n")
print(as.matrix(t(S)))
print(as.matrix(S %*% eigenVec(c(1, 1, 1))))
cat("单个元素 S(2,2):", eigen_sparse_get(S, 2, 2), " S(1,3):",
    eigen_sparse_get(S, 1, 3), "\n")

cat("\n08-sparse.R 结束\n")
