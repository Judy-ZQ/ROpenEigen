# ===========================================================================
# 06 · 解线性方程组
#
# 对应 EigenExample：
#   examples/TutorialLinAlgExSolveColPivHouseholderQR.cpp
#   examples/TutorialLinAlgExSolveLDLT.cpp
#   examples/TutorialLinAlgExComputeSolveError.cpp
#   examples/TutorialLinAlgComputeTwice.cpp
#   examples/TutorialLinAlgInverseDeterminant.cpp
#   examples/TutorialLinAlgRankRevealing.cpp
#   examples/TutorialLinAlgSVDSolve.cpp
#   examples/TutorialLinAlgSelfAdjointEigenSolver.cpp
#   examples/TutorialLinAlgSetThreshold.cpp
#   examples/Tutorial_PartialLU_solve.cpp
#   examples/TutorialInplaceLU.cpp
#   examples/matrixfree_cg.cpp
#   snippets/Tutorial_solve_{matrix_inverse,multiple_rhs,reuse_decomposition,
#                            singular,triangular,triangular_inplace}.cpp
#   snippets/MatrixBase_{inverse,computeInverseWithCheck,
#                        computeInverseAndDetWithCheck,eigenvalues}.cpp
#   snippets/class_FullPivLU.cpp
#
# 运行：Rscript 06-linear-solve.R
#
# 对照要点：Eigen 里"A.colPivHouseholderQr().solve(b)"这种链式写法，在 R 里有
# 两种对应：
#   一次性  —— eigen_solve(A, b)，内部选一种分解、算完就丢；
#   可复用  —— 先 eigen_colpiv_qr(A) 得到分解对象，再 $solve(b)。分解只做一次，
#              这正是 ROpenEigen 存在的理由（R 里 solve() 每次都重新分解）。
# ===========================================================================

suppressPackageStartupMessages(library(ROpenEigen))

# ---- TutorialLinAlgExSolveColPivHouseholderQR.cpp ---------------------------
# cpp:  A << 1,2,3, 4,5,6, 7,8,10;  b << 3,3,4;  x = A.colPivHouseholderQr().solve(b)
A <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 10), 3, 3, byrow = TRUE)
b <- c(3, 3, 4)
print(eigen_solve(eigenMat(A), b))
# 指定分解器时用对象形式：
print(eigen_colpiv_qr(eigenMat(A))$solve(b))

# ---- TutorialLinAlgExSolveLDLT.cpp ------------------------------------------
# cpp 用 A.ldlt().solve(b)，A 是对称正定矩阵。
A <- matrix(c(2, -1, -1, 2), 2, 2, byrow = TRUE)
b <- c(1, 0)
print(eigen_ldlt(eigenMat(A))$solve(b))

# ---- TutorialLinAlgExComputeSolveError.cpp ----------------------------------
# cpp 算 ‖Ax - b‖ 检查解的残差。R 里一行就够。
x <- as.vector(eigen_solve(eigenMat(A), b))
cat("residual:", sqrt(sum((A %*% x - b)^2)), "\n")

# ---- TutorialLinAlgComputeTwice.cpp -----------------------------------------
# cpp 演示"同一次分解算两次结构不同的右端项"。R 里复用分解对象即可，
# 分解只发生一次——这是与 base R 的 solve() 最大的区别。
A <- matrix(c(4, 1, 1, 3), 2, 2, byrow = TRUE)
lu <- eigen_lu(eigenMat(A))
print(lu$solve(c(1, 2)))
print(lu$solve(c(3, 4)))
# 对照：solve(A, ...) 每次都重新分解
print(solve(A, c(1, 2)))

# ---- TutorialLinAlgInverseDeterminant.cpp -----------------------------------
# cpp:  A.determinant(), A.inverse()
A <- matrix(c(1, 2, 1, 2, 1, 0, -1, 1, 2), 3, 3, byrow = TRUE)
cat("det:", eigen_det(eigenMat(A)), "\n")
print(eigen_inverse(eigenMat(A)))

# ---- MatrixBase_inverse.cpp -------------------------------------------------
print(eigen_inverse(eigenMat(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))))

# ---- MatrixBase_computeInverseWithCheck.cpp / ...AndDetWithCheck.cpp --------
# cpp 在求逆前先看是否可逆（FullPivLU 的 isInvertible()），避免对奇异矩阵
# 求逆得到一堆 inf。R 里同样要先检查：eigen_inverse() 走的是 Eigen 的
# 「假定可逆」的快路径，对奇异矩阵**不报错**，只会返回毫无意义的值。
S <- matrix(c(1, 2, 2, 4), 2, 2, byrow = TRUE)     # 奇异
cat("rcond:", eigen_rcond(eigenMat(S)), " 可逆?", eigen_rcond(eigenMat(S)) > 1e-12,
    " rank:", eigen_rank(eigenMat(S)), "\n")
cat("行列式:", eigen_det(eigenMat(S)), "\n")
cat("对奇异矩阵求逆的结果（无意义，只用来演示为什么必须先检查）:\n")
print(eigen_inverse(eigenMat(S)))
cat("有意义的做法是伪逆:\n")
print(eigen_pinv(eigenMat(S)))
cat("（base R 的 solve() 在这种情况下会直接报 'computationally singular'；\n")
cat("  eigen_inverse() 不会，因为它刻意不做那次 rcond 检查以省一遍 SVD。）\n")

# ---- TutorialLinAlgRankRevealing.cpp ----------------------------------------
# cpp 用 FullPivLU 取 rank、null-space 基、column-space 基。
A <- matrix(c(1, 2, 5, 2, 1, 4, 3, 0, 3), 3, 3, byrow = TRUE)
lu <- eigen_full_lu(eigenMat(A))
cat("rank:", lu$rank(), "\n")
# ROpenEigen 没有 $kernel()/$image()，用 SVD 得到同一组基：
sv <- eigen_svd(eigenMat(A))
V <- as.matrix(sv$v()); s <- as.vector(as.matrix(sv$singular_values()))
cat("奇异值:", round(s, 6), "\n")
print(round(V[, s < 1e-10 * max(s), drop = FALSE], 6))    # 零空间基
print(round(as.matrix(sv$u())[, s > 1e-10 * max(s), drop = FALSE], 6))  # 列空间基

# ---- TutorialLinAlgSVDSolve.cpp ---------------------------------------------
# cpp 用 bdcSvd().solve(b) 求最小二乘解（A 是 3x2，超定）。
set.seed(42)
A <- matrix(runif(6, -1, 1), 3, 2)
b <- runif(3)
print(eigen_svd(eigenMat(A))$solve(b))          # 最小二乘
cat("正规方程对照:", solve(t(A) %*% A, t(A) %*% b), "\n")

# ---- TutorialLinAlgSelfAdjointEigenSolver.cpp -------------------------------
# cpp:  SelfAdjointEigenSolver<Matrix2f> eigensolver(A); eigensolver.eigenvalues()
A <- matrix(c(1, 2, 2, 3), 2, 2, byrow = TRUE)
es <- eigen_eigen(eigenMat(A))
print(as.vector(as.matrix(es$values())))
print(es$vectors())

# ---- MatrixBase_eigenvalues.cpp ---------------------------------------------
# 同上的简写形式：eigen_values() 直接给特征值向量。
print(as.vector(as.matrix(eigen_values(eigenMat(A)))))

# ---- TutorialLinAlgSetThreshold.cpp -----------------------------------------
# cpp 调 rank-revealing 分解的 setThreshold()，控制"多小算零"。
# ROpenEigen 里对应的是 SVD 的 threshold 与 full LU 的 rank 判定：
sv <- eigen_svd(eigenMat(A))
cat("默认阈值:", sv$threshold(), " 秩:", sv$rank(), "\n")

# ---- Tutorial_PartialLU_solve.cpp -------------------------------------------
# cpp 用 PartialPivLU 解方程，并说明它比一般 inverse 稳定。
A <- matrix(c(3, 1, 1, 2), 2, 2, byrow = TRUE)
lu <- eigen_lu(eigenMat(A))
print(lu$solve(c(9, 8)))
print(as.matrix(lu$inverse() %*% c(9, 8)))

# ---- TutorialInplaceLU.cpp --------------------------------------------------
# cpp 在分解对象上原地解、看 L/U/P 三个因子。
A <- matrix(c(4, 3, 6, 3), 2, 2, byrow = TRUE)
lu <- eigen_lu(eigenMat(A))
print(lu$matrixL()); print(lu$matrixU()); print(lu$permutationP())
cat("重构误差:", max(abs(as.matrix(lu$reconstructed()) - A)), "\n")

# ---- Tutorial_solve_matrix_inverse.cpp --------------------------------------
# cpp:  x = A.inverse() * b  与  x = A.lu().solve(b) 的对比
A <- eigenMat(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
b <- c(5, 11)
print(eigen_inverse(A) %*% b)
print(eigen_lu(A)$solve(b))

# ---- Tutorial_solve_multiple_rhs.cpp ----------------------------------------
# cpp:  一次解多个右端项（矩阵 B），比逐个解快。
A <- eigenMat(matrix(c(2, 1, 1, 3), 2, 2, byrow = TRUE))
B <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3)
print(eigen_lu(A)$solve(B))

# ---- Tutorial_solve_reuse_decomposition.cpp ---------------------------------
# cpp 强调"分解一次、解多次"。R 里就是复用同一个分解对象。
lu <- eigen_lu(A)
for (k in 1:3) cat("rhs", k, ":", as.vector(lu$solve(c(k, 2 * k))), "\n")

# ---- Tutorial_solve_singular.cpp --------------------------------------------
# cpp 对奇异矩阵求解，靠 rank-revealing 分解给出"最小二乘意义上的解"。
S <- matrix(c(1, 1, 1, 1), 2, 2, byrow = TRUE)
cat("rank:", eigen_rank(eigenMat(S)), "\n")
print(eigen_full_lu(eigenMat(S))$solve(c(1, 1)))    # 无解时给出一个解
print(eigen_pinv(eigenMat(S)) %*% c(1, 1))          # 伪逆解

# ---- Tutorial_solve_triangular.cpp / _inplace.cpp ---------------------------
# cpp 用 triangularView().solve() 解三角系统；R 里是 backsolve()/forwardsolve()。
U <- matrix(c(2, 1, 0, 3), 2, 2, byrow = TRUE)      # 上三角
L <- t(U)                                            # 下三角
b <- c(4, 9)
print(backsolve(U, b))
print(forwardsolve(L, b))

# ---- class_FullPivLU.cpp ----------------------------------------------------
# cpp 逐个展示 FullPivLU 的接口。R 侧同名方法都在 $ 里，列一遍：
A <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 10), 3, 3, byrow = TRUE)
lu <- eigen_full_lu(eigenMat(A))
cat("rank:", lu$rank(), " 可逆:", lu$is_invertible(), " rcond:", lu$rcond(), "\n")
cat("det:", lu$det(), "\n")
cat("重构误差:", max(abs(as.matrix(lu$reconstructed()) - A)), "\n")
print(lu$permutationP()); print(lu$permutationQ())

# ---- matrixfree_cg.cpp ------------------------------------------------------
# cpp 演示"不给矩阵、只给算子"的共轭梯度：用户提供 A*v 的计算方式。
# R 里同样的想法是把矩阵乘写成一个函数，自己迭代——这里给一个最小实现，
# 与 eigen_sparse_solve(S, b, "cg") 的结果对照，说明封装里跑的是同一个算法。

A <- matrix(c(4, 1, 1, 3), 2, 2, byrow = TRUE)
b <- c(1, 2)

cg <- function(Afun, b, x = rep(0, length(b)), tol = 1e-12, maxit = 100) {
  r <- b - Afun(x); p <- r; rs <- sum(r * r)
  for (k in seq_len(maxit)) {
    Ap <- Afun(p)
    alpha <- rs / sum(p * Ap)
    x <- x + alpha * p
    r <- r - alpha * Ap
    rs_new <- sum(r * r)
    if (sqrt(rs_new) < tol) break
    p <- r + (rs_new / rs) * p
    rs <- rs_new
  }
  list(x = x, iterations = k)
}

res <- cg(function(v) A %*% v, b)
cat("matrix-free CG:", res$x, " 迭代", res$iterations, "次\n")
cat("直接解:", solve(A, b), "\n")

cat("\n06-linear-solve.R 结束\n")
