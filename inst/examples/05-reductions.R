# ===========================================================================
# 05 · 归约、访问者与广播
#
# 对应 EigenExample：
#   examples/Tutorial_ReductionsVisitorsBroadcasting_{colwise,rowwise,visitors,
#              reductions_bool,reductions_norm,maxnorm,reductions_operatornorm,
#              broadcast_simple,broadcast_1nn,broadcast_simple_rowwise}.cpp
#   snippets/PartialRedux_{sum,prod,count,maxCoeff,minCoeff,norm,squaredNorm}.cpp
#   snippets/VectorwiseOp_homogeneous.cpp / Vectorwise_reverse.cpp
#   snippets/DirectionWise_{replicate,replicate_int,hnormalized}.cpp
#   snippets/MatrixBase_{all,prod,operatorNorm,isDiagonal,isIdentity,isOnes,
#                        isZero,isOrthogonal,isUnitary}.cpp
#   snippets/MatrixBase_colwise.cpp / MatrixBase_rowwise.cpp
#   snippets/tut_arithmetic_redux_minmax.cpp
#   snippets/SelfAdjointView_{eigenvalues,operatorNorm}.cpp
#
# 运行：Rscript 05-reductions.R
#
# 广播这一节值得单独说明：R 的向量回收（recycling）是**按列**进行的，
# 恰好等于 Eigen 的 colwise() 广播；rowwise() 广播则要显式转置。
# ===========================================================================

suppressPackageStartupMessages(library(ROpenEigen))

mat <- matrix(c(1, 2, 6, 9, 3, 1, 7, 2), 2, 4, byrow = TRUE)   # 2x4

# ---- Tutorial_ReductionsVisitorsBroadcasting_colwise.cpp --------------------
# cpp:  mat.colwise().maxCoeff()      R: apply(mat, 2, max)
print(apply(mat, 2, max))

# ---- Tutorial_ReductionsVisitorsBroadcasting_rowwise.cpp --------------------
# cpp:  mat.rowwise().maxCoeff()      R: apply(mat, 1, max)
print(apply(mat, 1, max))

# ---- PartialRedux_sum.cpp / _prod.cpp / _count.cpp --------------------------
# colwise().sum() / .prod() / .count(谓词)
print(apply(mat, 2, sum))
print(apply(mat, 2, prod))
print(apply(mat, 2, function(x) sum(x > 2)))       # count(x > 2)

# ---- PartialRedux_maxCoeff.cpp / _minCoeff.cpp ------------------------------
print(apply(mat, 2, max))
print(apply(mat, 2, min))

# ---- PartialRedux_norm.cpp / _squaredNorm.cpp -------------------------------
# colwise().norm() 是每列的欧氏范数；R 里 sqrt(colSums(mat^2))。
print(sqrt(colSums(mat^2)))
print(colSums(mat^2))

# ---- Tutorial_ReductionsVisitorsBroadcasting_visitors.cpp -------------------
# cpp 用 maxCoeff(&row, &col) 同时拿到值和位置。R 里 which.max + arrayInd。
m <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
k <- which.max(m)
cat("Max:", m[k], "at:", arrayInd(k, dim(m))[1] - 1, ",",
    arrayInd(k, dim(m))[2] - 1, "\n")     # 减 1 换成 Eigen 的 0-based

k <- which.min(m)
cat("Min:", m[k], "at:", arrayInd(k, dim(m))[1] - 1, ",",
    arrayInd(k, dim(m))[2] - 1, "\n")

# 包里的 eigen_which_max()/eigen_which_min() 直接给 1-based 的线性下标：
print(eigen_which_max(eigenMat(m)))

# ---- Tutorial_ReductionsVisitorsBroadcasting_reductions_bool.cpp -----------
# cpp:  (m.array() > 0).all() / .any() / .count()
m <- matrix(c(1, -2, 3, 4), 2, 2, byrow = TRUE)
print(all(m > 0))
print(any(m > 0))
print(sum(m > 0))

# ---- Tutorial_ReductionsVisitorsBroadcasting_reductions_norm.cpp ------------
# cpp 比较 v.norm() 与 v.squaredNorm()；R 里 norm(v, "2") 与 sum(v^2)。
v <- c(3, 4)
cat("norm:", sqrt(sum(v^2)), " squared:", sum(v^2), "\n")

# ---- Tutorial_ReductionsVisitorsBroadcasting_maxnorm.cpp --------------------
# cpp:  m.maxCoeff(), m.minCoeff() 以及按行/列的 maxCoeff。
print(max(m)); print(min(m))

# ---- Tutorial_ReductionsVisitorsBroadcasting_reductions_operatornorm.cpp ----
# cpp:  m.operatorNorm() 即谱范数（最大奇异值）。R 里可以直接用 eigen_norm()。
m <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
print(eigen_norm(eigenMat(m), "2"))

# ---- MatrixBase_operatorNorm.cpp --------------------------------------------
print(eigen_norm(eigenMat(m), "2"))

# ---- MatrixBase_prod.cpp / MatrixBase_all.cpp -------------------------------
# cpp:  m.prod() 连乘；m.all() 是否全为真
print(eigen_prod(eigenMat(m)))
print(all(m != 0))

# ---- MatrixBase_isDiagonal.cpp / isIdentity.cpp / isOnes.cpp / isZero.cpp ---
# cpp 的这几个断言返回 bool。R 里没有对应函数，直接检验非对角元/差值。
is_diag <- function(x, tol = 0) {
  x <- as.matrix(x); all(abs(x[row(x) != col(x)]) <= tol)
}
is_ident <- function(x, tol = 0) {
  x <- as.matrix(x)
  nrow(x) == ncol(x) && all(abs(x - diag(nrow(x))) <= tol)
}
cat("是对角阵:", is_diag(diag(3)), " 是单位阵:", is_ident(diag(3)), "\n")
cat("一般矩阵:", is_diag(matrix(1:4, 2)), " 是单位阵:", is_ident(matrix(1:4, 2)), "\n")
cat("全 1:", all(matrix(1, 2, 2) == 1), " 全 0:", all(matrix(0, 2, 2) == 0), "\n")

# ---- MatrixBase_isOrthogonal.cpp / isUnitary.cpp ----------------------------
# cpp:  m.isOrthogonal(1e-4) —— 判 A^T A ≈ I。R 里一步算残差。
q <- matrix(c(0, -1, 1, 0), 2, 2, byrow = TRUE)      # 旋转 90°，正交
cat("正交残差:", max(abs(t(q) %*% q - diag(2))), "\n")

qc <- matrix(c(1, 1i, 1i, 1) / sqrt(2), 2, 2)        # 酉矩阵
cat("酉残差:", max(Mod(Conj(t(qc)) %*% qc - diag(2))), "\n")

# ---- tut_arithmetic_redux_minmax.cpp ----------------------------------------
# cpp 同时取 minCoeff/maxCoeff 与它们的位置。
m <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
cat("min:", min(m), "at", which.min(m), " max:", max(m), "at", which.max(m), "\n")

# ---- Tutorial_ReductionsVisitorsBroadcasting_broadcast_simple.cpp -----------
# cpp:  mat.colwise() += v     把 2 维向量加到每一列
# R 的回收按列进行，写 mat + v 就是同一个意思。
v <- c(0, 1)
print(eigenMat(mat + v))

# ---- Tutorial_ReductionsVisitorsBroadcasting_broadcast_simple_rowwise.cpp ---
# cpp:  mat.rowwise() += v.transpose()     把 4 维行向量加到每一行
# R 里按行广播要显式转置：t(t(mat) + v)
v <- c(0, 1, 2, 3)
print(eigenMat(t(t(mat) + v)))

# ---- Tutorial_ReductionsVisitorsBroadcasting_broadcast_1nn.cpp --------------
# cpp 求"离 v 最近的一列"：先 colwise() - v，再 colwise().squaredNorm()，
# 再取最小值的下标。R 里一行 apply 就能写完。
m <- matrix(c(1, 23, 6, 9, 3, 11, 7, 2), 2, 4, byrow = TRUE)
v <- c(2, 3)
d <- colSums((m - v)^2)
k <- which.min(d)
cat("Nearest neighbour is column", k - 1, ":\n")     # cpp 输出 0-based
print(m[, k])

# ---- Vectorwise_reverse.cpp -------------------------------------------------
# cpp:  m.rowwise().reverse() / m.colwise().reverse()
m <- matrix(1:12, 3, 4, byrow = TRUE)
print(eigenMat(t(apply(m, 1, rev))))     # 每行反转
print(eigenMat(apply(m, 2, rev)))        # 每列反转

# ---- VectorwiseOp_homogeneous.cpp / DirectionWise_hnormalized.cpp -----------
# 按列齐次化 = 末尾补一行 1；按行齐次化 = 末尾补一列 1。
print(eigenMat(rbind(m, 1)))
print(eigenMat(cbind(m, 1)))

# ---- DirectionWise_replicate.cpp / DirectionWise_replicate_int.cpp ----------
# cpp:  v.replicate(2, 3)：纵向复制 2 份、横向 3 份，2x1 的向量变成 4x3。
v <- c(1, 2)
print(eigen_rep(eigenVec(v), 2, 3))

# ---- SelfAdjointView_eigenvalues.cpp ---------------------------------------
# cpp 用 selfadjointView<Lower>().eigenvalues() 只读下三角求特征值。
# R 里 eigendecomposition 都吃完整矩阵，先把下三角补成对称矩阵。
lower <- matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE)
sym <- lower + t(lower) - diag(diag(lower))
print(sort(as.vector(as.matrix(eigen_eigen(eigenMat(sym))$values())), decreasing = TRUE))

# ---- SelfAdjointView_operatorNorm.cpp --------------------------------------
# 对称矩阵的谱范数 = 特征值绝对值的最大值。
cat("operatorNorm:", max(abs(eigen(as.matrix(sym), symmetric = TRUE)$values)), "\n")

# ===========================================================================
# 尾声：一个把逐元素运算与归约都用上的完整例子
# ---- demos/mandelbrot/mandelbrot.cpp --------------------------------------
# cpp 的数学内核只有一句话：对复平面上每个格点 c，迭代 z <- z^2 + c，
# 直到 |z| > 2 或到达最大迭代次数，逃逸所需的步数决定该点颜色。
# Qt 那一层只负责把结果刷成像素；内核本身是"整块 Array 做逐元素运算"，
# 在 R 里就是向量化的矩阵运算 —— 不用写任何内层循环。
# ===========================================================================
mandel_iter <- function(xr, yr, maxit = 50) {
  nr <- length(yr); nc <- length(xr)
  cr <- matrix(rep(xr, each = nr), nr, nc)     # broadcast：每列同一个实部
  ci <- matrix(rep(yr, times = nc), nr, nc)    # 每行同一个虚部
  zr <- zi <- matrix(0, nr, nc)
  cnt <- matrix(maxit, nr, nc)
  live <- matrix(TRUE, nr, nc)
  for (k in seq_len(maxit)) {
    zr2 <- zr * zr - zi * zi + cr              # (a+bi)^2 + c 的实部
    zi  <- 2 * zr * zi + ci                    #               虚部
    zr  <- zr2
    esc <- live & (zr * zr + zi * zi > 4)
    cnt[esc] <- k
    live[esc] <- FALSE
    zr[!live] <- 0; zi[!live] <- 0             # 已逃逸的点固定在 0，别溢出
    if (!any(live)) break
  }
  cnt
}

cnt <- mandel_iter(seq(-2, 0.5, length.out = 48), seq(-1.2, 1.2, length.out = 24))
cat("最高迭代次数:", max(cnt), "\n")
cat("逃逸（< 50 步）的格点占比:", round(mean(cnt < 50), 4), "\n")

# 用字符画看一眼形状：'@' 表示一直没逃逸（在集合内），空格表示 1 步就逃逸。
ramp <- c(" ", ".", "+", "*", "#", "@")   # 迭代次数越多越"深"
idx <- matrix(cut(as.vector(cnt), breaks = c(-Inf, 1, 3, 6, 12, 25, Inf),
                  labels = FALSE), nrow(cnt))
for (i in seq_len(nrow(idx))) cat(paste(ramp[idx[i, ]], collapse = ""), "\n")

cat("\n05-reductions.R 结束\n")
