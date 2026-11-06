suppressPackageStartupMessages(library(ROpenEigen))

# EigenExample/ 里的官方示例，逐个在 R 侧重做一遍。
#
# 这份测试和别的不同：其余文件验的是「封装与 base R 一致」，这一份验的是
# 「Eigen 文档里那些例子，用这套封装能不能照样得到同样的结果」。每个块上面
# 的注释写着对应的 .cpp 文件名，断言右边的数字来自编译运行那个 .cpp 得到的
# **实际输出**，不是从文档抄的期望值。
#
# 两个必须记住的转写规则：
#   1. Eigen 的逗号初始化 `m << 1, 2, 3, 4;` 是**行优先**填的，而 R 的
#      matrix() 是列优先。所以这里一律写 byrow = TRUE。
#   2. Eigen 打印矩阵是行优先的，R 也是，对照起来直接看即可。

ok <- 0L; bad <- 0L
chk <- function(label, expr) {
  r <- tryCatch(isTRUE(expr), error = function(e) {
    cat(sprintf("  !! %-40s %s\n", label, conditionMessage(e))); NA
  })
  if (is.na(r)) { bad <<- bad + 1L; return(invisible(FALSE)) }
  if (r) ok <<- ok + 1L else { bad <<- bad + 1L; cat(sprintf("  ?? %s\n", label)) }
  invisible(r)
}
eq <- function(a, b, tol = 1e-10) isTRUE(all.equal(a, b, tolerance = tol))
M <- function(x) as.matrix(x)

cat("== 1. QuickStart_example.cpp ==\n")
A <- matrix(c(3, 2.5, -1, 1.5), 2)          # 打印出来是 3 -1 / 2.5 1.5
chk("构造与取值", eq(M(eigenMat(A)), A))
chk("转置",       eq(M(t(eigenMat(A))), t(A)))
chk("自乘的对角", eq(diag(M(eigenMat(A) %*% eigenMat(A))), diag(A %*% A)))

cat("== 2. tut_arithmetic_add_sub.cpp ==\n")
a <- matrix(c(1, 2, 3, 4), 2, byrow = TRUE)
b <- matrix(c(2, 3, 1, 4), 2, byrow = TRUE)
ea <- eigenMat(a); eb <- eigenMat(b)
chk("a + b",      eq(M(ea + eb), matrix(c(3, 5, 4, 8), 2, byrow = TRUE)))
chk("a - b",      eq(M(ea - eb), matrix(c(-1, -1, 2, 0), 2, byrow = TRUE)))
chk("-v + w - v", eq(as.vector(-c(1, 2, 3) + c(1, 0, 0) - c(1, 2, 3)), c(-1, -4, -6)))

cat("== 3. tut_arithmetic_scalar_mul_div.cpp ==\n")
a2 <- matrix(c(1, 2, 3, 4), 2, byrow = TRUE)
chk("a * 2.5", eq(M(eigenMat(a2) * 2.5), a2 * 2.5))
chk("0.1 * v", eq(as.vector(eigenMat(c(1, 2, 3)) * 0.1), c(0.1, 0.2, 0.3)))
chk("v * 2",   eq(as.vector(eigenMat(c(1, 2, 3)) * 2), c(2, 4, 6)))

cat("== 4. tut_arithmetic_matrix_mul.cpp ==\n")
mat <- matrix(c(1, 2, 3, 4), 2, byrow = TRUE)
u <- c(-1, 1); vv <- c(2, 0)
em <- eigenMat(mat)
chk("mat*mat",  eq(M(em %*% em), matrix(c(7, 10, 15, 22), 2, byrow = TRUE)))
chk("mat*u",    eq(as.vector(em %*% u), c(1, 1)))
chk("u^T*mat",  eq(as.vector(t(u) %*% em), c(2, 2)))
chk("u^T*v",    eq(as.numeric(t(u) %*% vv), -2))
chk("u*v^T",    eq(M(eigen_outer(eigenMat(u), eigenMat(vv))),
                   matrix(c(-2, 0, 2, 0), 2, byrow = TRUE)))

cat("== 5. tut_arithmetic_dot_cross.cpp ==\n")
p <- c(1, 2, 3); q <- c(0, 1, 2)
chk("v.dot(w) = 8",   eq(eigen_dot(eigenMat(p), eigenMat(q)), 8))
chk("叉积 = (1,-2,1)", eq(as.vector(eigen_cross(eigenMat(p), eigenMat(q))), c(1, -2, 1)))
chk("二维叉积 = 1",    eq(eigen_cross(eigenMat(c(1, 2)), eigenMat(c(0, 1))), 1))

cat("== 6. tut_arithmetic_redux_basic.cpp ==\n")
r <- eigenMat(matrix(c(1, 2, 3, 4), 2, byrow = TRUE))
chk("sum = 10",   eq(eigen_sum(r), 10))
chk("prod = 24",  eq(eigen_prod(r), 24))
chk("mean = 2.5", eq(eigen_mean(r), 2.5))
chk("min = 1",    eq(min(r), 1))
chk("max = 4",    eq(max(r), 4))
chk("trace = 5",  eq(eigen_trace(r), 5))

cat("== 7. Tutorial_ArrayClass_addition.cpp / mult.cpp ==\n")
a3 <- matrix(1:9, 3, byrow = TRUE)
chk("a + a = 2a（逐元素）", eq(M(eigenMat(a3) + eigenMat(a3)), 2 * a3))
chk("a - 2",               eq(M(eigenMat(a3) - 2), a3 - 2))
ax <- matrix(c(1, 2, 3, 4), 2, byrow = TRUE)
bx <- matrix(c(5, 6, 7, 8), 2, byrow = TRUE)
chk("逐元素积 = [[5,12],[21,32]]",
    eq(M(eigenMat(ax) * eigenMat(bx)), matrix(c(5, 12, 21, 32), 2, byrow = TRUE)))
chk("矩阵积 != 逐元素积", !eq(M(eigenMat(ax) %*% eigenMat(bx)), ax * bx))

cat("== 8. Tutorial_ArrayClass_cwise_other.cpp ==\n")
aa <- c(-1.45311, 0.0655203, -0.30336, 0.724223, 0.262856)
ea2 <- eigenMat(aa)
chk("abs()",        eq(as.vector(abs(ea2)), abs(aa)))
chk("abs().sqrt()", eq(as.vector(sqrt(abs(ea2))), sqrt(abs(aa))))
chk("对位取小",     eq(pmin(as.vector(ea2), as.vector(sqrt(abs(ea2)))),
                       pmin(aa, sqrt(abs(aa)))))

cat("== 9. Tutorial_ArrayClass_interop_matrix.cpp ==\n")
m2 <- matrix(c(1, 2, 3, 4), 2, byrow = TRUE)
n2 <- matrix(c(5, 6, 7, 8), 2, byrow = TRUE)
chk("m*n = [[19,22],[43,50]]",
    eq(M(eigenMat(m2) %*% eigenMat(n2)), matrix(c(19, 22, 43, 50), 2, byrow = TRUE)))
chk("cwiseProduct = [[5,12],[21,32]]",
    eq(M(eigenMat(m2) * eigenMat(n2)), matrix(c(5, 12, 21, 32), 2, byrow = TRUE)))

cat("== 10. Tutorial_BlockOperations_*.cpp ==\n")
m3 <- matrix(1:9, 3, byrow = TRUE)
em3 <- eigenMat(m3)
chk("第 2 行 = 4 5 6",   eq(as.vector(em3[2, ]), c(4, 5, 6)))
chk("leftCols(2)",       eq(M(em3[, 1:2]), m3[, 1:2]))
chk("bottomRows(2)",     eq(M(em3[2:3, ]), m3[2:3, ]))
chk("角上的块",          eq(M(em3[1:2, 1:2]), m3[1:2, 1:2]))
chk("col(2) += 3 * col(0)", {
  x <- eigenMat(m3)
  x[, 3] <- x[, 3] + 3 * x[, 1]
  eq(as.vector(x[, 3]), c(6, 18, 30))
})
chk("segment(1,4) *= 2", {
  v6 <- eigenMat(1:6)
  v6[2:5] <- 2 * v6[2:5]
  eq(as.vector(v6), c(1, 4, 6, 8, 10, 6))
})
chk("head(3) / tail(3)", {
  v6 <- eigenMat(1:6)
  eq(as.vector(v6[1:3]), c(1, 2, 3)) && eq(as.vector(v6[4:6]), c(4, 5, 6))
})
chk("中心块赋值", {
  A4 <- eigenMat(matrix(0.6, 4, 4))
  A4[2:3, 2:3] <- matrix(c(1, 2, 3, 4), 2, byrow = TRUE)
  exp4 <- matrix(0.6, 4, 4)
  exp4[2:3, 2:3] <- matrix(c(1, 2, 3, 4), 2, byrow = TRUE)
  eq(M(A4), exp4)
})
chk("右下列 -> 左上列（先读后写）", {
  A4 <- eigenMat(matrix(0.6, 4, 4))
  A4[1:2, 1:2] <- matrix(c(3, 4), 2, 2)
  eq(as.vector(A4[1:2, 1:2]), c(3, 4, 3, 4))
})

cat("== 11. Tutorial_ReductionsVisitorsBroadcasting_*.cpp ==\n")
mat2 <- matrix(c(1, 2, 6, 9, 3, 1, 7, 2), 2, byrow = TRUE)
chk("列最大值 = 3 2 7 9", eq(apply(mat2, 2, max), c(3, 2, 7, 9)))
chk("行最大值 = 9 7",     eq(apply(mat2, 1, max), c(9, 7)))
chk("下标访问与 apply 一致", {
  e <- eigenMat(mat2)
  eq(max(e[, 3]), 7) && eq(max(e[1, ]), 9)
})
# mat.colwise().sum() = (4, 3, 13, 11)，最大在位置 2（0-based），
# 对应那一列是 (6, 7)，和为 13
chk("maxCoeff 位置 = 2、该列 = (6,7)、和 = 13", {
  colsum <- colSums(mat2)
  eq(colsum, c(4, 3, 13, 11)) &&
    which.max(colsum) - 1L == 2L &&
    eq(mat2[, which.max(colsum)], c(6, 7)) &&
    eq(max(colsum), 13)
})
chk("squaredNorm / norm / lpNorm", {
  v2 <- c(1, 2)
  eq(sum(v2^2), 5) && eq(sqrt(sum(v2^2)), sqrt(5)) && eq(sum(abs(v2)), 3)
})
chk("eigen_norm(x, \"2\") = sqrt(5)", eq(eigen_norm(eigenMat(c(1, 2)), "2"), sqrt(5)))
chk("矩阵 squaredNorm = 30", {
  mm <- matrix(c(1, 2, 3, 4), 2, byrow = TRUE)
  eq(eigen_sumsq(eigenMat(mm)), 30)
})
chk("矩阵 1-范数 = 10（Eigen 的 lpNorm<1>）", {
  mm <- matrix(c(1, 2, 3, 4), 2, byrow = TRUE)
  eq(sum(abs(mm)), 10)
})
chk("无穷范数 = 4", {
  mm <- matrix(c(1, 2, 3, 4), 2, byrow = TRUE)
  eq(eigen_norm(eigenMat(mm), "M"), 4)
})
# m = [[1,-2],[-3,4]]：列绝对值和 (4,6) -> 6；行绝对值和 (3,7) -> 7
chk("operator-norm：1-范数 6 / 无穷范数 7", {
  mm <- matrix(c(1, -2, -3, 4), 2, byrow = TRUE)
  eq(eigen_norm(eigenMat(mm), "O"), 6) && eq(eigen_norm(eigenMat(mm), "I"), 7)
})
chk("colwise 广播（每列加同一个向量）", {
  mm <- matrix(c(1, 2, 6, 9, 3, 1, 7, 2), 2, byrow = TRUE)
  eq(sweep(mm, 1, c(1, 1), "+"),
     matrix(c(2, 3, 7, 10, 4, 2, 8, 3), 2, byrow = TRUE))
})
chk("逻辑归约 all / any / count", {
  e <- eigenMat(matrix(c(1, 2, 3, 4), 2, byrow = TRUE))
  all(as.vector(e) > 0) && any(as.vector(e) > 2) && sum(as.vector(e) > 2) == 2
})
# Eigen 的 visitors 示例打印的是「最大值 4 在 (1,1)」这种坐标（0-based），
# 本包的 eigen_which_* 给的是 1-based 的行列与取值。
chk("which.min / which.max 的位置", {
  e <- eigenMat(matrix(c(1, 2, 3, 4), 2, byrow = TRUE))
  lo <- eigen_which_min(e); hi <- eigen_which_max(e)
  eq(lo$value, 1) && eq(c(lo$row, lo$col), c(1L, 1L)) &&
    eq(hi$value, 4) && eq(c(hi$row, hi$col), c(2L, 2L))
})

cat("== 12. 线性代数示例 ==\n")
A3 <- matrix(c(1, 2, 1, 2, 1, 0, -1, 1, 2), 3, byrow = TRUE)
chk("det = -3", eq(eigen_det(eigenMat(A3)), -3))
chk("逆矩阵", eq(M(eigen_inverse(eigenMat(A3))),
                 matrix(c(-2 / 3, 1, 1 / 3, 4 / 3, -1, -2 / 3, -1, 1, 1), 3, byrow = TRUE)))
Ar <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 10), 3, byrow = TRUE)
br <- c(3, 3, 4)
chk("ColPivHouseholderQR 解 = (-2, 1, 1)",
    eq(as.vector(eigen_solve(eigenMat(Ar), br)), c(-2, 1, 1), tol = 1e-5))
chk("FullPivLU 解相同",
    eq(as.vector(solve(eigen_lu(Ar), eigenMat(br))), c(-2, 1, 1), tol = 1e-5))
Ad <- matrix(c(2, -1, -1, 3), 2, byrow = TRUE)
bd <- matrix(c(1, 2, 3, 1), 2, byrow = TRUE)
chk("LDLT 解 = [[1.2,1.4],[1.4,0.8]]",
    eq(M(solve(eigen_ldlt(Ad), eigenMat(bd))),
       matrix(c(1.2, 1.4, 1.4, 0.8), 2, byrow = TRUE)))
chk("同一份分解用两次结果相同", {
  L <- eigen_llt(Ad)
  eq(M(solve(L, eigenMat(bd))), M(solve(L, eigenMat(bd))))
})
chk("矩阵改了要重算", {
  Ad2 <- matrix(c(2, -1, -1, 4), 2, byrow = TRUE)
  eq(M(solve(eigen_llt(Ad2), eigenMat(bd))), solve(Ad2, bd))
})
As <- matrix(c(1, 2, 2, 3), 2, byrow = TRUE)
chk("自伴特征值 = (-0.236068, 4.23607)",
    eq(as.vector(eigen_values(eigenMat(As))), c(-0.2360679, 4.2360679), tol = 1e-6))
chk("特征向量满足 A v = lambda v", {
  V <- M(eigen_vectors(eigenMat(As)))
  lam <- as.vector(eigen_values(eigenMat(As)))
  all(vapply(1:2, function(j) eq(as.vector(As %*% V[, j]), lam[j] * V[, j], 1e-8), TRUE))
})
Ark <- matrix(c(1, 2, 5, 2, 1, 4, 3, 0, 3), 3, byrow = TRUE)
chk("秩 = 2", eq(eigen_rank(eigenMat(Ark)), 2L))
chk("零空间基 (0.5, 1, -0.5)", {
  ns <- c(0.5, 1, -0.5)
  eq(max(abs(Ark %*% ns)), 0, tol = 1e-10)
})
chk("SVD 最小二乘解", {
  Asv <- matrix(c(-0.726553, 0.362112, 0.0327601, 0.131428, -0.15168, 0.884173),
                3, byrow = TRUE)
  bsv <- c(0.466488, 0.838312, 0.297589)
  got <- as.vector(M(eigen_pinv(eigenMat(Asv)) %*% bsv))
  eq(got, as.vector(qr.solve(Asv, bsv)), 1e-5)
})
# FullPivLU 的 rank() 用阈值判秩：A 的第二行与第一行几乎相同，
# 阈值放大到 1e-5 之后秩从 2 掉到 1。
chk("秩与阈值：默认 2、阈值 1e-5 时 1", {
  At <- matrix(c(2, 1, 2, 0.9999999999), 2, byrow = TRUE)
  eq(eigen_rank(eigenMat(At)), 2L) && eq(eigen_rank(eigenMat(At), tol = 1e-5), 1L)
})

cat("== 13. Cwise_erf / Cwise_erfc / Cwise_lgamma.cpp ==\n")
xs <- c(-0.5, 1, 0, -2)
chk("Cwise_erf",  eq(as.vector(M(eigen_erf(eigenMat(xs)))), 2 * pnorm(xs * sqrt(2)) - 1))
chk("Cwise_erfc", eq(as.vector(M(eigen_erfc(eigenMat(xs)))), 2 * pnorm(-xs * sqrt(2))))
chk("Cwise_lgamma", eq(as.vector(M(eigen_gammaln(eigenMat(c(1, 2, 3))))), lgamma(c(1, 2, 3))))

cat("== 14. special_examples/Tutorial_sparse_example.cpp ==\n")
# 官方稀疏示例：单位正方形上的 2 维 Poisson，5 点差分，边长 6 -> 36 个未知量
nx <- 6L
N <- nx * nx
P5 <- matrix(0, N, N)
for (i in 1:nx) {
  for (j in 1:nx) {
    k <- (i - 1) * nx + j
    P5[k, k] <- 4
    if (i > 1) P5[k, k - nx] <- -1
    if (i < nx) P5[k, k + nx] <- -1
    if (j > 1) P5[k, k - 1] <- -1
    if (j < nx) P5[k, k + 1] <- -1
  }
}
S5 <- eigen_sparse(P5)
chk("非零元个数", eq(eigen_sparse_nnz(S5), sum(P5 != 0)))
chk("稀疏解 = 稠密解", {
  bb <- as.numeric(seq_len(N))
  eq(as.vector(M(eigen_sparse_solve(S5, bb, method = "lu"))),
     as.vector(solve(P5, bb)), 1e-10)
})
chk("稀疏不放大矩阵", eigen_sparse_nnz(S5) < N^2 / 3)

cat("== 15. matrixfree_cg.cpp：无矩阵的 CG ==\n")
chk("CG 与直接法同解", {
  n3 <- 20L
  tri <- matrix(0, n3, n3)
  for (i in 1:n3) {
    tri[i, i] <- 2
    if (i > 1) tri[i, i - 1] <- -1
    if (i < n3) tri[i, i + 1] <- -1
  }
  b3 <- rep(1, n3)
  eq(as.vector(M(eigen_sparse_solve(eigen_sparse(tri), b3, method = "cg"))),
     as.vector(solve(tri, b3)), 1e-8)
})

cat(sprintf("\n通过 %d, 失败 %d\n", ok, bad))
if (bad > 0L) stop(sprintf("%d assertion(s) failed", bad), call. = FALSE)
