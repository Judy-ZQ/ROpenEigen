# ===========================================================================
# 实验五 · 特征值：把图像的「主方向」找出来
#
# 对应教材：第五章 §5.1–§5.3（特征值、相似、实对称矩阵的对角化）
# 建议学时：3
#
# 实验目的
#   1. 看到特征向量就是「被 A 作用后方向不变」的那些方向；
#   2. 用 PCA 把一张图的主轴找出来，画在图上；
#   3. 看到实对称矩阵的正交对角化：转过去、拉开、转回来；
#   4. 看到幂迭代：反复作用同一个矩阵，任何向量都会被拉向主特征方向。
#
# 运行：Rscript img05-eigen.R   或在 RStudio 里 source()
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

img <- il_read()
A <- il_gray(img)
h <- nrow(A); w <- ncol(A)

rot <- function(deg) {
  t <- deg * pi / 180
  matrix(c(cos(t), -sin(t), sin(t), cos(t)), 2, byrow = TRUE)
}
deg <- function(v) atan2(v[2], v[1]) * 180 / pi

# ===========================================================================
# 一、把图的「形状」变成一个 2×2 矩阵
# ===========================================================================
say("\n== 5.1 图像的二阶矩：一张 2×2 的协方差矩阵 ==")

# 把每个像素看成平面上的一点，亮度当权重，算加权协方差。
# 亮的地方「重」，所以协方差矩阵描述的是「亮块摊开的形状」。
xs <- rep(seq_len(w), each = h)
ys <- rep(seq_len(h), times = w)
wt <- as.vector(A)
wt <- wt / sum(wt)
mx <- sum(wt * xs); my <- sum(wt * ys)
dx <- xs - mx; dy <- ys - my
C <- matrix(c(sum(wt * dx * dx), sum(wt * dx * dy),
              sum(wt * dx * dy), sum(wt * dy * dy)), 2, byrow = TRUE)
say("亮度的重心（加权平均位置）= (", round(mx, 2), ", ", round(my, 2), ")")
say("协方差矩阵 C = ")
print(round(C, 3))

lam <- as.vector(eigen_values(eigenMat(C)))
V <- as.matrix(eigen_vectors(eigenMat(C)))
say("特征值：", paste(round(lam, 3), collapse = ", "), "（eigen_values 给的是升序）")
say("特征向量（列）：")
print(round(V, 4))
say("主轴方向 = ", round(deg(V[, 2]), 2), " 度（对应较大的特征值 ", round(lam[2], 3), "）")
say("trace(C) = ", round(sum(diag(C)), 4), "，两个特征值之和 = ", round(sum(lam), 4), "（相似不变量）")
say("det(C) = ", format(det(C), digits = 6), "，两个特征值之积 = ",
    format(prod(lam), digits = 6))

# 把两条轴画在图上：长度取 3 个标准差
overlay <- img
for (k in 1:2) {
  s <- 3 * sqrt(lam[k])
  v <- V[, k]
  overlay <- il_draw(overlay, mx - s * v[1], my - s * v[2],
                     mx + s * v[1], my + s * v[2], val = c(255, 0, 0), width = 2)
}
say("\n红线 = 两个特征方向，长度 = 3 倍标准差（sqrt(特征值)）。")
say("它们一定互相垂直：<v1, v2> = ", format(sum(V[, 1] * V[, 2]), scientific = TRUE),
    "（实对称矩阵的特征向量正交）")

# 协方差是二次型：把整幅点云旋转 R，新的协方差就是 R C R'。
# 它的特征值一个都不变，特征向量全部跟着转 R —— 这正好预演了 5.3 的对角化。
Rr <- rot(30)
C2 <- Rr %*% C %*% t(Rr)
lam2 <- as.vector(eigen_values(eigenMat(C2)))
V2 <- as.matrix(eigen_vectors(eigenMat(C2)))
say("把点云整体转 30 度：特征值 ", paste(round(lam2, 3), collapse = ", "),
    " 一个都没变；主轴从 ", round(deg(V[, 2]), 2), " 度变成 ", round(deg(V2[, 2]), 2), " 度。")
expect_near(lam2, lam, 1e-10, "旋转不该改变特征值")

il_fig_sheet("img05-axes", list(img, overlay),
             c("original", "eigen-directions of the brightness distribution"), ncol = 2)

# ===========================================================================
# 二、A v = λ v：方向不变，只伸缩
# ===========================================================================
say("\n== 5.2 沿特征方向拉伸，形状不变 ==")

S <- V %*% diag(c(1, 2.2)) %*% t(V)          # 沿主轴放大 2.2 倍，另一个方向不变
say("S = V·diag(1, 2.2)·V'，det S = ", round(eigen_det(eigenMat(S)), 6),
    "，特征值 = ", paste(round(as.vector(eigen_values(eigenMat(S))), 4), collapse = ", "))

# S 作用在图上：沿主轴拉长。因为 S 的特征方向就是 V，图像不会被「拧」。
il_fig_sheet("img05-stretch",
             list(A, il_warp_inv(A, S, out = c(240L, 400L))),
             c("original", "stretched along the main axis only (eigen-direction)"), ncol = 2)

# ===========================================================================
# 三、正交对角化：转过去、拉开、转回来
# ===========================================================================
say("\n== 5.3 Q' A Q = Λ：三步做完同一个变换 ==")

# 取一个对称矩阵 T，用它的特征分解
Ssym <- matrix(c(1.6, 0.6,
                 0.6, 1.0), 2, byrow = TRUE)
Q <- as.matrix(eigen_vectors(eigenMat(Ssym)))
L <- diag(as.vector(eigen_values(eigenMat(Ssym))))
say("T = [[1.6, 0.6], [0.6, 1.0]]，特征值 ", paste(round(diag(L), 4), collapse = ", "))
say("Q 正交：Q'Q - I 的最大偏差 = ", format(max(abs(crossprod(Q) - diag(2))), scientific = TRUE))
say("对角化还原：||Q Λ Q' - T|| = ", format(max(abs(Q %*% L %*% t(Q) - Ssym)), scientific = TRUE))

# 先在「点」上验证：这才是代数意义上的相等
P <- matrix(runif(200, -100, 100), 2)
say("200 个点分别用 T 和 Q Λ Q' 变换，最大坐标差 = ",
    format(max(abs(Ssym %*% P - Q %*% L %*% t(Q) %*% P)), scientific = TRUE))

# 再在「像素」上验证：三次重采样一定会有误差，但那不是矩阵的错
direct <- il_warp_inv(A, Ssym, out = c(300L, 380L))
stepwise <- il_warp_inv(il_warp_inv(il_warp_inv(A, t(Q), out = c(300L, 380L)), L, out = c(300L, 380L)),
                        Q, out = c(300L, 380L))
say("直接按 T 变换 与 「先转 Q'、沿轴拉伸、再转回 Q」：RMSE = ",
    round(il_rmse(direct, stepwise), 3))
say("  这点差别全部来自「三次最近邻重采样」，与 1e-16 的代数误差无关。")

il_fig_sheet("img05-diag",
             list(A, il_warp_inv(A, t(Q), out = c(300L, 380L)),
                  il_warp_inv(il_warp_inv(A, t(Q), out = c(300L, 380L)), L, out = c(300L, 380L)),
                  direct),
             c("original", "step 1: rotate into the eigenbasis (Q')",
               "step 2: scale along the axes (Lambda)",
               "T applied directly (after step 3: rotate back)"), ncol = 4)

# ===========================================================================
# 四、幂迭代：任何方向都会被拉向主特征方向
# ===========================================================================
say("\n== 5.4 幂迭代：为什么 A^k 会「失忆」 ==")

S2 <- Ssym
lam2 <- as.vector(eigen_values(eigenMat(S2)))
say("T 的两个特征值相差 ", round(lam2[2] / lam2[1], 4), " 倍。")
say("每作用一次，两个方向的比值就乘上这个数；作用 k 次就是 ", round(lam2[2] / lam2[1], 4), "^k。")

v <- c(1, 0)                                  # 随便挑一个初始方向
say("\n从 v0 = (1, 0) 开始迭代 v <- T v / ||T v||：")
for (k in 1:6) {
  v <- as.vector(S2 %*% v); v <- v / sqrt(sum(v^2))
  say("  k = ", k, "：方向 ", formatC(round(deg(v), 4), width = 8, format = "f"),
      " 度，与主特征方向的夹角 ",
      formatC(round(abs(deg(v) - deg(Q[, 2])), 6), width = 9, format = "f"), " 度")
}
say("角度每轮缩小 ", round(lam2[1] / lam2[2], 4), " 倍 —— 收敛速度由「次特征值/主特征值」决定。")

# ===========================================================================
# 五、反复作用：图像被拉成条纹
# ===========================================================================
say("\n== 5.5 把 A^k 打在图上：所有东西都朝主方向倒 ==")

Q2 <- rot(35)
S3 <- Q2 %*% diag(c(1.03, 0.97)) %*% t(Q2)     # 缓慢的各向异性拉伸
say("S = Q·diag(1.03, 0.97)·Q'，Q 转 35 度：每作用一次，35 度方向长 3%，垂直方向缩 3%。")
for (k in c(1, 8, 24, 60)) {
  Ak <- Q2 %*% diag(c(1.03^k, 0.97^k)) %*% t(Q2)
  say("  k = ", formatC(k, width = 2), "：两个方向的长度比 = ",
      round(1.03^k / 0.97^k, 2), " 倍")
}

il_fig_sheet("img05-iterate",
             list(A,
                  il_warp_inv(A, Q2 %*% diag(c(1.03^8, 0.97^8)) %*% t(Q2), out = c(200L, 420L)),
                  il_warp_inv(A, Q2 %*% diag(c(1.03^24, 0.97^24)) %*% t(Q2), out = c(200L, 420L)),
                  il_warp_inv(A, Q2 %*% diag(c(1.03^60, 0.97^60)) %*% t(Q2), out = c(200L, 420L))),
             c("original", "k = 8", "k = 24", "k = 60: almost a stripe"), ncol = 4)

say("这就是 det 和特征值的分工：det 管「面积变几倍」，特征值管「各方向分别变几倍」。")
say("k 足够大时，最大的那个特征值独占一切，图像塌成一维 —— 和实验一的 det = 0 是同一件事的极限形式。")

# ===========================================================================
# 六、实验习题
# ===========================================================================
say("\n== 5.6 习题 ==")
say("1. 把测试图换成一张竖条纹图（例如 (row %% 20 < 10) * 255），主轴方向是多少度？")
say("   为什么？再用 eigen_values() 看两个特征值差多少。")
say("2. 5.4 节里初始方向换成 (0, 1)，迭代 6 次后角度序列一样吗？")
say("3. 5.5 节如果把 1.03 改成 1.003，k 要多大才看得出条纹？")
say("   这说明了「幂迭代收敛慢」在数值上的原因。")
say("4. 为什么 eigen_eigen() 只能用于对称矩阵？用它算旋转矩阵会得到什么？")
say("   （提示：回看实验七的 5.1 节，或者直接试 [[0,-1],[1,0]]。）")

say("\nimg05 结束：特征方向是「不转向的方向」，对角化是「换一套坐标做同一件事」。")
