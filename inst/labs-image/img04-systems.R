# ===========================================================================
# 实验四 · 线性方程组：把模糊的图还原回来
#
# 对应教材：第四章 §4.1–§4.3（同解、齐次、非齐次、最小二乘）
# 建议学时：3
#
# 实验目的
#   1. 把「图像模糊」写成一个线性方程组 Ax = b —— 未知数是原图，b 是模糊图；
#   2. 看到病态：条件数一大，直接解出来的图全是噪声；
#   3. 用正则化 (A'A + λI)x = A'b 把它救回来，并理解 λ 在做什么；
#   4. 用超高维的最小二乘拟合一条直线（320 个方程、2 个未知数）。
#
# 运行：Rscript img04-systems.R   或在 RStudio 里 source()
# ===========================================================================

# 定位 imghelpers.R：脚本自己所在目录 -> 当前目录（及往上的包源码树）-> 已安装的包。
# 三种入口都能跑：Rscript imgXX.R、在 RStudio 里 source()、包内回归测试。
.lab_dir <- local({
  cand <- character(0)
  # (1) 脚本自己在哪里：Rscript 把路径放进 --file=，source() / RStudio 把它放进 ofile
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) cand <- dirname(sub("^--file=", "", f[1]))
  for (fr in rev(sys.frames())) {
    of <- fr$ofile
    if (is.character(of) && length(of) == 1L && file.exists(of)) {
      cand <- c(cand, dirname(of)); break
    }
  }
  # (2) 当前目录，以及从当前目录往上 4 层找包源码树（<root>/inst/labs-image）
  d <- getwd()
  for (i in 0:4) {
    cand <- c(cand, d, file.path(d, "inst", "labs-image"), file.path(d, "labs-image"))
    d <- dirname(d)
  }
  # (3) 已安装的包（安装副本可能比源码树旧，所以放在最后）
  cand <- c(cand, system.file("labs-image", package = "ROpenEigen"))
  hit <- cand[file.exists(file.path(cand, "imghelpers.R"))]
  if (!length(hit))
    stop("找不到 imghelpers.R。三种跑法都支持：在 inst/labs-image 下 Rscript、",
         "在 RStudio 里 source()、或先安装 ROpenEigen。", call. = FALSE)
  hit[1]
})
source(file.path(.lab_dir, "imghelpers.R"))
suppressPackageStartupMessages({
  library(ROpenCV)
  library(ROpenEigen)
})

A <- il_gray(il_read())
h <- nrow(A); w <- ncol(A)
say("测试图：", h, " × ", w, "，图会存到：", il_outdir())

# ===========================================================================
# 一、模糊 = 矩阵乘
# ===========================================================================
say("\n== 4.1 一行像素的模糊，就是一个 320×320 的方程组 ==")

# 高斯模糊核（15 个抽头，σ = 3），归一化后和为 1：模糊不改变总亮度
tap <- 7
kk <- dnorm(-tap:tap, 0, 3)
kk <- kk / sum(kk)
n <- w

# 一维模糊的矩阵：Toeplitz，第 (i,j) 元是核在 |i-j| 处的值
Tmat <- outer(1:n, 1:n, function(i, j) {
  d <- abs(i - j)
  ifelse(d <= tap, kk[d + tap + 1], 0)     # 核的中心在 kk[tap+1]
})
Tmat <- Tmat / rowSums(Tmat)               # 边界上核被截断，按行归一化补回来
say("核：", paste(round(kk, 4), collapse = ", "), "，和为 ", round(sum(kk), 6))
say("T 是 ", n, "×", n, "，每行只有一个「凸起」，所以每行和为 1 —— 亮度的总和不变。")

# 挑一行「变化最剧烈」的：穿过字母和斜线，模糊的后果最明显
gd <- rowSums(abs(A[, -1] - A[, -w]))
ri <- which.max(gd)
x_true <- A[ri, ]
say("取第 ", ri, " 行（左邻右舍的差最大），共 ", w, " 个像素。")
b <- as.vector(Tmat %*% x_true)        # 模糊后的同一行
say("这一行：原图总亮度 ", round(sum(x_true), 1), "，模糊后 ", round(sum(b), 1), "（几乎相等）")
say("最暗处：原图 ", round(min(x_true), 1), "，模糊后 ", round(min(b), 1), " —— 黑线被周围的浅色摊薄了。")

# ===========================================================================
# 二、直接解：病态会把噪声放大
# ===========================================================================
say("\n== 4.2 直接求解：条件数说话 ==")

say("条件数 kappa(T) = ", format(kappa(Tmat), digits = 6),
    "（2-范数下最大与最小奇异值之比；eigen_rcond 给的是倒数）")
x_hat <- as.vector(solve(eigenMat(Tmat), eigenVec(b)))
say("无噪声时直接解：与原行的 RMSE = ", round(il_rmse(x_hat, x_true), 8), "（几乎是精确的）")

# 加一点点噪声：模拟相机噪声，标准差 1 个灰阶
set.seed(1)
bn <- b + rnorm(n, 0, 1)
x_noisy <- as.vector(solve(eigenMat(Tmat), eigenVec(bn)))
say("只加了 σ = 1 的噪声（相当于 ±1 个灰阶），解出来的图：")
say("  与原行的 RMSE = ", round(il_rmse(x_noisy, x_true), 2),
    "，最大值 ", round(max(abs(x_noisy)), 1), "（像素上限才 255）")
say("  1 个灰阶的噪声，解出来偏差 ", round(il_rmse(x_noisy, x_true), 0),
    " 个灰阶 —— 放大了四个数量级，这就是病态。")
say("  病态的根源：T 的最小奇异值 ", format(min(svd(Tmat, nu = 0, nv = 0)$d), scientific = TRUE),
    " 太小，b 里一点点扰动，除以它就被放大成巨浪。")

# ===========================================================================
# 三、正则化：给方程加一条「别太离谱」的约束
# ===========================================================================
say("\n== 4.3 正则化：(T'T + λI) x = T'b ==")

ridge <- function(lam) {
  M <- crossprod(Tmat) + lam * diag(n)
  as.vector(solve(eigenMat(M), eigenVec(as.vector(crossprod(Tmat, bn)))))
}
lams <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100)
say("λ 从 0.0001 扫到 100，看 RMSE 怎么变：")
best <- c(Inf, NA)
for (lam in lams) {
  xr <- ridge(lam)
  e <- il_rmse(xr, x_true)
  say("  λ = ", formatC(lam, width = 6), "：RMSE = ", round(e, 3),
      "，最大值 ", round(max(abs(xr)), 1))
  if (e < best[1]) best <- c(e, lam)
}
say("最好的 λ = ", best[2], "，RMSE = ", round(best[1], 3))

say("\nλ 在做什么：它把 (T'T)^{-1} 里那些特别小的特征值抬起来，")
say("  于是「噪声除以小特征值」这一步被挡住了；代价是解变平滑、细节变模糊。")
say("  这就是「偏差-方差权衡」在方程组上的样子。")

il_fig_sheet("img04-1d",
             list(matrix(x_true, 1), matrix(b, 1), matrix(pmax(pmin(x_noisy, 255), 0), 1),
                  matrix(ridge(1), 1), matrix(ridge(100), 1)),
             c("original row", "blurred (b = T x)", "naive solve (noise x100)",
               "ridge, lambda = 1", "ridge, lambda = 100"), ncol = 1, aspect = 0.18, panel = 900)

# ===========================================================================
# 四、整幅图：240 个方程组，共用一次分解
# ===========================================================================
say("\n== 4.4 整幅图：240 个方程组共用一份分解 ==")

# ROpenEigen 的分解对象可以复用：分解一次，解 240 个右端项。
# Bmat 的每一列是图像的一行（320 个数），共 240 列。
Bmat <- t(A)
blur2d <- t(Tmat %*% Bmat)                     # 正问题：b = T x
say("整幅图模糊后与原图的 RMSE = ", round(il_rmse(blur2d, A), 2),
    "（这是「模糊」本身，不是求解误差）")

set.seed(2)
Bnoise <- t(blur2d) + matrix(rnorm(length(Bmat), 0, 1), nrow = n)
lu <- eigen_lu(Tmat)                            # 分解一次
naive2d <- t(as.matrix(solve(lu, eigenMat(Bnoise))))
say("加 σ = 1 的噪声后直接解：RMSE = ", round(il_rmse(naive2d, A), 2),
    "，解的取值范围 ", paste(round(range(naive2d), 1), collapse = " .. "),
    "（像素只该在 0..255）")

ridge2d <- function(lam) {
  M <- eigen_lu(crossprod(Tmat) + lam * diag(n))
  t(as.matrix(solve(M, eigenMat(as.matrix(crossprod(Tmat, Bnoise))))))
}
# λ 的量级要跟 T'T 的对角线（约 1）比，不能照搬一维那张表
best2d <- c(Inf, NA)
for (lam in c(0.001, 0.003, 0.01, 0.03, 0.1, 1)) {
  e <- il_rmse(ridge2d(lam), A)
  if (e < best2d[1]) best2d <- c(e, lam)
}
reg2d <- ridge2d(best2d[2])
say("二维的 λ 扫一遍（0.001 … 1），最好的 λ = ", best2d[2],
    "：RMSE = ", round(best2d[1], 2),
    "，取值范围 ", paste(round(range(reg2d), 1), collapse = " .. "))

# 和 OpenCV 的高斯模糊对一下：同一个核，同一个结果
cvblur <- as_array(cv_GaussianBlur(il_cv(A), c(15, 15), 3), numeric = TRUE)
say("\n和 cv_GaussianBlur(15×15, σ=3) 对比：两者 RMSE = ",
    round(il_rmse(blur2d, cvblur), 3),
    "（矩阵乘出来的模糊和 OpenCV 的模糊是同一件事，差别只在边界处理）")

il_fig_sheet("img04-2d",
             list(A, blur2d, il_clamp(naive2d), il_clamp(reg2d)),
             c("original", "b = T x (blurred)", "naive solve: noise explodes",
               sprintf("ridge, lambda = %g", best2d[2])), ncol = 4)

# ===========================================================================
# 五、最小二乘：320 个方程，2 个未知数
# ===========================================================================
say("\n== 4.5 超定方程组：拿像素点拟合一条直线 ==")

# 测试图里有一条黑色斜线（见 il_chart）。每列找最黑的像素，得到一组点。
# 只保留「确实有黑色像素」的列：纯背景的列里每一点都一样黑，which.min()
# 会随便返回第一行 —— 那不是数据，用错数据比用错方法更致命。
cols_all <- seq_len(w)
has_dark <- vapply(cols_all, function(j) min(A[, j]) < 60, logical(1))
cols <- cols_all[has_dark & cols_all <= 150]
say("有 ", length(cols), " 列能取到线上的点（其余的列背景是纯色，取不出点）。")
y <- vapply(cols, function(j) which.min(A[, j]), integer(1))
# 拟合 y = a + b·x，写成超定方程 X c = y
X <- cbind(1, cols)
cf <- as.vector(eigen_pinv(eigenMat(X)) %*% eigenVec(as.numeric(y)))
say("拟合结果：行号 = ", round(cf[1], 4), " + ", round(cf[2], 6), " × 列号")
resid <- as.numeric(y) - as.vector(X %*% cf)
say("残差平方和 = ", round(sum(resid^2), 3), "，最大偏差 ", round(max(abs(resid)), 3), " 像素")

# 正规方程与 SVD 两条路，应当一致
cf_ne <- as.vector(solve(eigenMat(crossprod(X)), eigenVec(as.vector(crossprod(X, y)))))
say("正规方程给的系数与 SVD 差：", format(max(abs(cf_ne - cf)), scientific = TRUE))
say("（", length(cols), " 个方程、2 个未知数，方程远多于未知数 —— 这就是「超定」。）")
say("（病态时正规方程会失准，理由见实验七 7.2。）")

# ===========================================================================
# 六、实验习题
# ===========================================================================
say("\n== 4.6 习题 ==")
say("1. 把核宽度从 15 改成 31，kappa(T) 变成多少？为什么模糊越厉害越难还原？")
say("2. 4.2 节里噪声标准差是 1。改成 3、5，RMSE 怎么变？是线性增长吗？")
say("3. 用 eigen_sparse() 把 T 存成稀疏矩阵，再用 eigen_sparse_solve() 解，")
say("   结果和稠密解差多少？为什么这里稀疏解法特别合适？")
say("4. 4.5 节用的是「每列最黑的那个像素」。如果这条线是斜的很陡（接近竖直），")
say("   这种做法会出什么问题？该怎么改？")

say("\nimg04 结束：模糊是 Ax=b，病态是放大镜，正则化是刹车。")
