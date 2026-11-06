# ===========================================================================
# 07 · 分解本身：LLT / LU / QR / SVD / 特征值 / Schur / Hessenberg
#
# 对应 EigenExample：
#   snippets/LLT_example.cpp / LLT_solve.cpp
#   snippets/PartialPivLU_solve.cpp / FullPivLU_{solve,kernel,image}.cpp
#   snippets/HouseholderQR_{solve,householderQ}.cpp
#   snippets/ColPivHouseholderQR_solve.cpp / FullPivHouseholderQR_solve.cpp
#   snippets/JacobiSVD_basic.cpp / LeastSquaresQR.cpp / LeastSquaresNormalEquations.cpp
#   snippets/SelfAdjointEigenSolver_*.cpp                     （9 个）
#   snippets/EigenSolver_*.cpp                                （5 个）
#   snippets/ComplexEigenSolver_*.cpp                         （3 个）
#   snippets/ComplexSchur_*.cpp / RealSchur_*.cpp             （5 个）
#   snippets/HessenbergDecomposition_*.cpp                    （3 个）
#   snippets/Tridiagonalization_*.cpp                         （6 个）
#   snippets/RealQZ_compute.cpp / GeneralizedEigenSolver.cpp
#   snippets/HouseholderSequence_HouseholderSequence.cpp
#   snippets/Jacobi_makeGivens.cpp / Jacobi_makeJacobi.cpp
#   snippets/Triangular_solve.cpp
#
# 运行：Rscript 07-decompositions.R
#
# 本章与上一章的分工：上一章是"拿分解去解方程"，本章是"把分解对象本身
# 摊开看" —— 因子、置换、奇异值、特征向量、秩、行列式。
# 凡是 ROpenEigen 没有对应物的（Schur、Hessenberg、Tridiagonalization），
# 下面会明确写"没有"，并给出最接近的做法，而不是硬凑一个。
# ===========================================================================

suppressPackageStartupMessages(library(ROpenEigen))

A <- matrix(c(4, 1, 1, 3), 2, 2, byrow = TRUE)        # 对称正定

# ---- LLT_example.cpp --------------------------------------------------------
# cpp 构造 LLT，检查 info()，再取 matrixL()。R 里 $is_ok() 就是那个检查。
llt <- eigen_llt(eigenMat(A))
cat("ok:", llt$is_ok(), "\n")
print(llt$matrixL())
cat("L L^T 重构误差:", max(abs(as.matrix(llt$reconstructed()) - A)), "\n")
cat("logdet:", llt$logdet(), " 对照:", determinant(A, logarithm = TRUE)$modulus, "\n")

# ---- LLT_solve.cpp ----------------------------------------------------------
print(llt$solve(c(1, 2)))

# ---- eigen_chol：R 风格的 Cholesky -----------------------------------------
# 包另外提供 eigen_chol()，返回的是**矩阵**（上三角），与 base R 的 chol() 对齐。
print(eigen_chol(eigenMat(A)))
print(chol(A))

# ---- PartialPivLU_solve.cpp -------------------------------------------------
lu <- eigen_lu(eigenMat(matrix(c(3, 1, 1, 2), 2, 2, byrow = TRUE)))
print(lu$solve(c(9, 8)))
print(lu$permutationP())

# ---- FullPivLU_solve.cpp ----------------------------------------------------
lu <- eigen_full_lu(eigenMat(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 10), 3, 3, byrow = TRUE)))
print(lu$solve(c(3, 3, 4)))

# ---- FullPivLU_kernel.cpp / FullPivLU_image.cpp -----------------------------
# cpp 取零空间基与列空间基。ROpenEigen 的 full LU 没有 kernel()/image()，
# 用 SVD 得到同一组基：右奇异向量中标量值≈0 的那些张成零空间，非零的张成值域。
A3 <- matrix(c(1, 2, 5, 2, 1, 4, 3, 0, 3), 3, 3, byrow = TRUE)
sv <- eigen_svd(eigenMat(A3))
s <- as.vector(as.matrix(sv$singular_values()))
V <- as.matrix(sv$v()); U <- as.matrix(sv$u())
cat("奇异值:", round(s, 6), "\n")
cat("零空间基（dim =", sum(s < 1e-10 * max(s)), "）:\n")
print(round(V[, s < 1e-10 * max(s), drop = FALSE], 6))
cat("列空间基:\n")
print(round(U[, s > 1e-10 * max(s), drop = FALSE], 6))

# ---- HouseholderQR_solve.cpp ------------------------------------------------
A <- matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = TRUE)   # 超定
b <- c(1, 2, 3)
print(eigen_qr(eigenMat(A))$solve(b))                  # 最小二乘
print(qr.solve(A, b))                                  # base R 对照

# ---- HouseholderQR_householderQ.cpp -----------------------------------------
# cpp 取 Q（Householder 序列本身）。R 里 $matrixQ() 直接给稠密 Q。
qr <- eigen_qr(eigenMat(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)))
print(qr$matrixQ())
print(qr$matrixR())

# ---- ColPivHouseholderQR_solve.cpp ------------------------------------------
A <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 10), 3, 3, byrow = TRUE)
q <- eigen_colpiv_qr(eigenMat(A))
print(q$solve(c(3, 3, 4)))
cat("rank:", q$rank(), " 可逆:", q$is_invertible(), "\n")

# ---- FullPivHouseholderQR_solve.cpp -----------------------------------------
# ROpenEigen 只封装了列主元 QR（colPiv），全主元 QR 没有单独暴露：
# Eigen 的 FullPivHouseholderQR 主要用来做"最稳的秩判定"，R 里用 full LU 或 SVD
# 得到同样信息。下面用 SVD 给出秩与数值秩阈值。
sv <- eigen_svd(eigenMat(A))
cat("rank(SVD):", sv$rank(), " threshold:", sv$threshold(), "\n")

# ---- JacobiSVD_basic.cpp ----------------------------------------------------
# cpp 用 JacobiSVD<MatrixXf> 算 2x2 的奇异值。R 里 eigen_svd(method="jacobi")。
A <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
sj <- eigen_svd(eigenMat(A), method = "jacobi")
sb <- eigen_svd(eigenMat(A), method = "bdc")
cat("jacobi:", as.vector(as.matrix(sj$singular_values())), "\n")
cat("bdc   :", as.vector(as.matrix(sb$singular_values())), "\n")
cat("base R:", svd(A)$d, "\n")

# ---- LeastSquaresQR.cpp -----------------------------------------------------
# cpp 比较 QR 与正规方程两种最小二乘。QR 更稳，正规方程会放大条件数。
A <- matrix(c(1, 1, 1, 1, 1, 2, 1, 3, 1, 4), 5, 2)
b <- c(2, 3, 5, 6, 8)
print(eigen_qr(eigenMat(A))$solve(b))
print(solve(t(A) %*% A, t(A) %*% b))

# ---- LeastSquaresNormalEquations.cpp ----------------------------------------
# 正规方程写法在 R 里就是上面第二行；条件数平方的关系可以看一眼：
cat("cond(A) =", kappa(A), " cond(A'A) =", kappa(t(A) %*% A), "\n")

# ---- SelfAdjointEigenSolver_eigenvalues.cpp ---------------------------------
A <- matrix(c(1, 2, 2, 3), 2, 2, byrow = TRUE)
es <- eigen_eigen(eigenMat(A))
print(as.vector(as.matrix(es$values())))

# ---- SelfAdjointEigenSolver_eigenvectors.cpp --------------------------------
print(es$vectors())

# ---- SelfAdjointEigenSolver_compute_MatrixType*.cpp ------------------------
# cpp 演示"构造时算"与"构造后 compute()"两种用法。R 里构造即计算，
# 没有单独的 compute() 步骤（要重算就再调一次 eigen_eigen()）。
print(as.vector(as.matrix(eigen_eigen(eigenMat(A))$values())))

# ---- SelfAdjointEigenSolver_operatorSqrt.cpp -------------------------------
# cpp:  A^{1/2} = V diag(sqrt(λ)) V^T。ROpenEigen 直接提供 eigen_sqrtm()。
# 注意 A 必须是对称正定的，否则没有实的平方根 —— 换一个特征值都为正的矩阵：
Aspd <- matrix(c(4, 1, 1, 3), 2, 2, byrow = TRUE)
es <- eigen_eigen(eigenMat(Aspd))
V <- as.matrix(es$vectors()); lam <- as.vector(as.matrix(es$values()))
print(eigen_sqrtm(eigenMat(Aspd)))
cat("手算对照误差:",
    max(abs(as.matrix(eigen_sqrtm(eigenMat(Aspd))) - V %*% diag(sqrt(lam)) %*% t(V))), "\n")

# ---- SelfAdjointEigenSolver_operatorInverseSqrt.cpp ------------------------
# cpp:  A^{-1/2}。R 里用同样的谱分解，或者 powm(A, -0.5)。
cat("powm(A, -0.5) 与谱分解的差:",
    max(abs(as.matrix(eigen_powm(eigenMat(Aspd), -0.5)) -
            V %*% diag(lam^(-0.5)) %*% t(V))), "\n")

# ---- EigenSolver_eigenvalues.cpp / _eigenvectors.cpp -----------------------
# 一般（非对称）实矩阵：cpp 的 EigenSolver 给复数特征值。R 里 eigen_eigen()
# 对非对称输入会自动改用 ComplexEigenSolver，返回的 $values() 是复数。
N <- matrix(c(0, -1, 1, 0), 2, 2, byrow = TRUE)     # 旋转 90°，特征值 ±i
ge <- eigen_geigen(eigenMat(N))
print(as.vector(as.matrix(ge$values())))
print(ge$vectors())

# ---- EigenSolver_compute.cpp / EigenSolver_EigenSolver_MatrixType.cpp ------
# 构造即计算，同上。这里补一个"实矩阵有实特征值"的例子：
print(as.vector(as.matrix(eigen_geigen(eigenMat(matrix(c(1, 2, 0, 3), 2, 2,
                                                     byrow = TRUE)))$values())))

# ---- EigenSolver_pseudoEigenvectors.cpp ------------------------------------
# cpp 的 pseudoEigenvectors() 是给"实 Schur 形式"用的实向量。R 里取实部：
cat("实部:\n"); print(round(Re(as.matrix(ge$vectors())), 6))

# ---- ComplexEigenSolver_compute.cpp / _eigenvalues.cpp / _eigenvectors.cpp -
# cpp 用 ComplexEigenSolver 直接处理复矩阵。R 里同样构造复矩阵喂给 eigen_geigen()。
C <- matrix(c(1 + 1i, 0, 0, 2 - 1i), 2, 2, byrow = TRUE)
ce <- eigen_geigen(eigenMat(C))
print(as.vector(as.matrix(ce$values())))

# ---- ComplexSchur_compute.cpp / _matrixT.cpp / _matrixU.cpp ----------------
# 复 Schur 分解 A = U T U^*，ROpenEigen 没有封装 ComplexSchur。
# 最接近的做法是直接用特征分解（对可对角化矩阵 T 就是特征值对角阵）：
cat("Schur T 的对角（= 特征值）:", as.vector(as.matrix(ce$values())), "\n")

# ---- RealSchur_compute.cpp / RealSchur_RealSchur_MatrixType.cpp ------------
# 实 Schur 分解 A = Q T Q^T。base R 没有，Matrix 包（ROpenEigen 的 Suggests）
# 里的 Schur() 有，返回 S4 对象，槽位用 $T / $Q 取。
S <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
if (requireNamespace("Matrix", quietly = TRUE)) {
  sc <- Matrix::Schur(S)
  print(round(sc$T, 6))
  print(round(sc$Q, 6))
  print(max(abs(sc$Q %*% sc$T %*% t(sc$Q) - S)))   # 重构误差
} else {
  print("未安装 Matrix 包，跳过实 Schur 示例")
}

# ---- HessenbergDecomposition_compute.cpp / _matrixH.cpp / _packedMatrix.cpp
# ROpenEigen 没有 Hessenberg 分解，base R 也没有。这里给出定义式的做法：
# 用 Householder 把矩阵逐步化为上 Hessenberg（教学用，不追求效率）。
hessenberg <- function(M) {
  H <- as.matrix(M); n <- nrow(H)
  if (n > 2) for (k in 1:(n - 2)) {
    x <- H[(k + 1):n, k]
    v <- x; v[1] <- v[1] + sign(x[1]) * sqrt(sum(x^2))
    if (sqrt(sum(v^2)) < 1e-14) next
    v <- v / sqrt(sum(v^2))
    H[(k + 1):n, ] <- H[(k + 1):n, ] - 2 * v %*% (t(v) %*% H[(k + 1):n, ])
    H[, (k + 1):n] <- H[, (k + 1):n] - 2 * (H[, (k + 1):n] %*% v) %*% t(v)
  }
  H
}
H <- hessenberg(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 17, 18), 4, 4))
below <- outer(1:4, 1:4, "-") > 1          # 第一条次对角线以下的那些位置
cat("下 Hessenberg 以下的元素应≈0:", max(abs(H[below])), "\n")
print(round(H, 6))

# ---- Tridiagonalization_*.cpp（6 个）---------------------------------------
# 对称矩阵的三对角化：eigen_tridiagonal 没有封装。base R 也没有直接函数，
# 但可以对对称矩阵做 Hessenberg（此时结果就是三对角）：
S <- matrix(c(4, 1, 0, 1, 3, 1, 0, 1, 2), 3, 3, byrow = TRUE)
T3 <- hessenberg(S)
print(round(T3, 6))
cat("副对角以外应≈0:", max(abs(T3[row(T3) - col(T3) > 1])), "\n")

# ---- RealQZ_compute.cpp / GeneralizedEigenSolver.cpp -----------------------
# 广义特征值问题 A x = λ B x，以及它的 Schur 形式（QZ）。
# R 侧都没有直接接口：base R 的 eigen() 不接受第二个矩阵，Matrix 包也没导出
# qz()。能算的是**同一个谱**的等价问题 B^{-1} A x = λ x —— 数学上等价，
# 数值上不如 QZ 稳（B 病态时会放大误差），所以这里只演示，不推荐当生产写法。
Aq <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
Bq <- matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE)
print(round(eigen(solve(Bq, Aq))$values, 6))
print(paste("cond(B) =", round(kappa(Bq), 3), "条件数越大这个替代越不可信"))

# ---- HouseholderSequence_HouseholderSequence.cpp ---------------------------
# HouseholderSequence 是 HouseholderQR 内部的那个序列对象，用来对多个向量
# 批量施加/转置施加。R 里没有暴露，等价物是 $matrixQ() 后显式相乘。
q <- eigen_qr(eigenMat(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = TRUE)))
Q <- as.matrix(q$matrixQ())
cat("Q 正交性误差:", max(abs(t(Q) %*% Q - diag(ncol(Q)))), "\n")

# ---- Jacobi_makeGivens.cpp / Jacobi_makeJacobi.cpp ------------------------
# 两个片段构造 2x2 旋转把某个元素清零。ROpenEigen 里可以直接用
# eigen_rotation2d(angle) 造同一个旋转。
a <- 3; b <- 4
# 要把 (a, b) 转到 x 轴上（第二个分量变 0），旋转角是 -atan2(b, a)：
# eigen_rotation2d(θ) 是 [[cos,-sin],[sin,cos]]，顺时针转 θ 才把 (a,b) 压到轴上。
G <- as.matrix(eigen_rotation2d(-atan2(b, a)))
print(round(G, 6))
cat("把 (a,b) 转到 x 轴上:", round(G %*% c(a, b), 6), "（第二个分量应为 0）\n")

# ---- Triangular_solve.cpp --------------------------------------------------
# cpp 的 triangularView<Upper>().solve(b)，R 里是 backsolve()。
U <- matrix(c(2, 1, 0, 3), 2, 2, byrow = TRUE)
print(backsolve(U, c(4, 9)))

cat("\n07-decompositions.R 结束\n")
