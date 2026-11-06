# ===========================================================================
# 实验二 · 矩阵运算：复合、逆、初等变换、秩
#
# 对应教材：第二章 §2.1–§2.6（运算、分块、秩、初等变换、逆矩阵）
# 建议学时：3
#
# 实验目的
#   1. 看到 AB 就是「先做 B 再做 A」，以及为什么 AB ≠ BA；
#   2. 看到逆矩阵真的把图变回去，误差只剩插值；
#   3. 把「左乘初等矩阵 = 对矩阵做初等行变换」直接打在像素上；
#   4. 看到秩：低秩图像长什么样，为什么秩小就可以压缩。
#
# 运行：Rscript img02-operations.R   或在 RStudio 里 source()
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
gray <- il_gray(img)
h <- nrow(gray); w <- ncol(gray)

rot <- function(deg) {
  t <- deg * pi / 180
  matrix(c(cos(t), -sin(t), sin(t), cos(t)), 2, byrow = TRUE)
}
say("测试图：", h, " × ", w, "，灰度矩阵的秩 = ", eigen_rank(eigenMat(gray)))
say("（测试图有大片纯色区域，行与行之间重复，所以秩小于 ", h,
    "；换成真实照片会更接近满秩。）")

# ===========================================================================
# 一、AB = 先 B 后 A
# ===========================================================================
say("\n== 2.1 矩阵乘就是变换的复合 ==")

R30 <- rot(30)
S2  <- matrix(c(2, 0, 0, 1), 2, byrow = TRUE)      # 横向拉长一倍
say("R30 = 旋转 30 度，S2 = 横向拉伸 2 倍")
say("det(R30) = ", round(eigen_det(eigenMat(R30)), 6), "，det(S2) = ", eigen_det(eigenMat(S2)),
    "，det(R30 %*% S2) = ", round(eigen_det(eigenMat(R30 %*% S2)), 6),
    " = 两者的乘积（det 是乘性的）")

# 顺序一：先拉伸，再旋转  ->  R30 %*% S2 %*% x
# 顺序二：先旋转，再拉伸  ->  S2 %*% R30 %*% x
img_RS <- il_warp_inv(gray, R30 %*% S2, out = c(330L, 560L))
img_SR <- il_warp_inv(gray, S2 %*% R30, out = c(330L, 560L))
say("先拉伸再旋转 与 先旋转再拉伸，两张图的 RMSE = ", round(il_rmse(img_RS, img_SR), 3),
    "（若可交换应当为 0）")

il_fig_sheet("img02-order",
             list(gray, img_RS, img_SR, 255 * (abs(img_RS - img_SR) > 40)),
             c("original", "rotate(R30) o stretch(S2)", "stretch(S2) o rotate(R30)",
               "where they differ (|diff| > 40)"), ncol = 4)

say("结论：矩阵乘法的顺序就是操作的顺序，而两个顺序一般不等价，")
say("      所以 AB != BA —— 图上是同一个变换套在另一个变换里面，还是套在外面。")

# ===========================================================================
# 二、逆矩阵把图变回去
# ===========================================================================
say("\n== 2.2 A^{-1} 真的能把图变回去吗 ==")

# 取一个「缩小 + 旋转」的 A：缩小保证整幅图都还留在画布里，
# 这样变回去时不会有内容被裁掉，误差就只剩下重采样。
A <- 0.85 * rot(20)
say("A = 缩小 0.85 倍 + 旋转 20 度，det A = ", round(eigen_det(eigenMat(A)), 6))

Ai <- as.matrix(solve(eigenMat(A)))
say("ROpenEigen 算出的 A^{-1}·A 与单位矩阵的最大偏差：",
    format(max(abs(Ai %*% A - diag(2))), scientific = TRUE))

moved <- il_warp_inv(gray, A, out = dim(gray))
canvas <- il_warp_inv(moved, Ai, out = dim(gray))
say("图先正变换再逆变换，整图 RMSE = ", round(il_rmse(canvas, gray), 3),
    "，PSNR = ", round(il_psnr(canvas, gray), 2), " dB —— 看起来很差？")
say("  但这个数主要是被黑框拉偏的：缩小 0.85 倍后，画布边缘有 ",
    round(100 * mean(moved < 1), 1), "% 的像素是补的黑色，")
say("  它们当然变不回来。只看中间的区域：RMSE = ",
    round(il_rmse(canvas[41:(h - 40), 41:(w - 40)], gray[41:(h - 40), 41:(w - 40)]), 3),
    "，PSNR = ",
    round(il_psnr(canvas[41:(h - 40), 41:(w - 40)], gray[41:(h - 40), 41:(w - 40)]), 2), " dB")
say("  而真正「变了」的像素只占 ", round(100 * mean(abs(canvas - gray) > 20), 2),
    "%，全在硬边缘上（取整造成的）。")
say("记住三层的区别：")
say("  A^{-1}A = I 精确到 1e-16 —— 线性代数意义上的可逆；")
say("  边缘像素差几十 —— 离散采样 + 最近邻取整的代价；")
say("  黑框那块差 255 —— 压根没有信息，跟矩阵无关。")

il_fig_sheet("img02-inverse",
             list(gray, moved, canvas, 255 * (abs(canvas - gray) > 20)),
             c("original", "A applied", "A^{-1} applied: back home",
               "residual (resampling)"), ncol = 4)

# ===========================================================================
# 三、初等矩阵 = 初等行变换（这次是字面意义上的「行」）
# ===========================================================================
say("\n== 2.3 左乘初等矩阵，就是给图像矩阵做初等行变换 ==")

# 行交换：单位矩阵的行倒过来，左乘就是把图像的行倒过来
Pswap <- diag(h)[h:1, ]
vflip <- Pswap %*% gray
say("行交换矩阵 P 左乘：图像上下翻转，与 gray[h:1, ] 的最大差 = ",
    max(abs(vflip - gray[h:1, ])))
expect(max(abs(vflip - gray[h:1, ])) < 1e-12, "行交换矩阵没把行倒过来")

# 行倍加：E[i, j] = k 表示「第 i 行 += k × 第 j 行」。
# 取 j = i - 12，就是「每一行叠上它上方 12 行的那一行」，于是出现一个
# 上移 12 像素的重影 —— 这就是教科书里那句话的字面执行。
Eadd <- diag(h)
Eadd[cbind(13:h, 1:(h - 12))] <- 0.6
ghost <- il_clamp(Eadd %*% gray)
say("行倍加 E[i, i-12] = 0.6：每个像素叠上了它上方 12 行的 60%，")
say("  于是图像出现向上偏移 12 像素的重影。全程只用了矩阵乘法。")

# 行倍乘：只放大第 100 行
Esc <- diag(h); Esc[100, 100] <- 3
bright <- il_clamp(Esc %*% gray)
say("行倍乘 E[100, 100] = 3：只有第 100 行被乘了 3 倍，")
say("  平均亮度从 ", round(mean(gray[100, ]), 2), " 变成 ", round(mean(bright[100, ]), 2),
    "，其余行不变：", max(abs(bright[-100, ] - gray[-100, ])) == 0)

il_fig_sheet("img02-elementary",
             list(gray, vflip, ghost, bright),
             c("original", "row swap = vertical flip", "row addition = ghosting",
               "row scaling = one bright line"), ncol = 4)

# ===========================================================================
# 四、秩：低秩图像长什么样
# ===========================================================================
say("\n== 2.4 秩决定了这张图「有多少种独立的东西」==")

# 秩 1 的图像：每一行都是同一个向量的倍数，所以整幅图只有一列信息
u <- seq(0, 1, length.out = w)
v <- sin(seq(0, 2.4 * pi, length.out = h))
rank1 <- 255 * outer(v, u) / max(abs(outer(v, u)))
say("outer(v, u) 造出来的图，秩 = ", eigen_rank(eigenMat(rank1)),
    "：每行都是 u 的倍数，只有一种「花样」。")
say("它看起来是一条条横纹 —— 因为第 i 行 = v[i]·u，只有明暗不同。")

# 截断 SVD：秩 k 的最佳近似（Eckart-Young），实验七会细讲
S <- eigen_svd(eigenMat(gray))
sv <- as.vector(S$singular_values())
U <- as.matrix(S$u()); V <- as.matrix(S$v())
approx <- function(k) U[, 1:k, drop = FALSE] %*% (sv[1:k] * t(V[, 1:k, drop = FALSE]))
say("\n灰度的奇异值：前 5 个是 ", paste(round(sv[1:5], 1), collapse = ", "))
say("用前 k 个奇异值重构（k = 1, 5, 20, 60）：")
for (k in c(1, 5, 20, 60)) {
  ak <- approx(k)
  say("  k = ", formatC(k, width = 2), "：RMSE = ", formatC(il_rmse(ak, gray), 8, format = "f"),
      "，秩 = ", eigen_rank(eigenMat(ak)))
}

il_fig_sheet("img02-rank",
             list(rank1, approx(1), approx(5), approx(20), approx(60), gray),
             c("rank 1 (outer product)", "grey, rank 1 approx", "rank 5",
               "rank 20", "rank 60", "original (rank 240)"), ncol = 3)

# ===========================================================================
# 五、手写映射 vs OpenCV 的双线性插值
# ===========================================================================
say("\n== 2.5 为什么库函数要「插值」 ==")

A2 <- rot(18)
hand <- il_warp_inv(gray, A2, out = dim(gray))       # 最近邻

# 同一件事交给 OpenCV：2×3 的仿射矩阵，平移项按「绕图像中心」补齐
cx <- (w + 1) / 2; cy <- (h + 1) / 2
M23 <- rbind(c(A2[1, ], cx - sum(A2[1, ] * c(cx, cy))),
             c(A2[2, ], cy - sum(A2[2, ] * c(cx, cy))))
cvres <- as_array(cv_warpAffine(il_cv(gray), Mat(M23), c(w, h),
                                INTER_LINEAR, BORDER_CONSTANT, 0), numeric = TRUE)

say("同一旋转：手写最近邻 vs cv_warpAffine（双线性）")
say("  两者的 RMSE = ", round(il_rmse(hand, cvres), 3))
say("  差异集中在边缘（|diff| > 30 的像素占 ",
    round(100 * mean(abs(hand - cvres) > 30), 2), "%）：")
say("  最近邻只会「搬」像素，双线性会「算」出中间值，所以斜边更平滑。")

il_fig_sheet("img02-interp",
             list(hand, cvres, 255 * (abs(hand - cvres) > 30)),
             c("nearest neighbour (hand-written)", "cv_warpAffine (bilinear)",
               "where they differ"), ncol = 3)

# ===========================================================================
# 六、实验习题
# ===========================================================================
say("\n== 2.6 习题 ==")
say("1. 找一个 2×2 的 A 使得 AB = BA 对所有 B 成立。用图像验证你的答案。")
say("2. 行倍加 E[i, j] = k 里的 k 取 0.4 和取 -0.4，图像为什么一个是拖影、")
say("   一个像浮雕？把 k 从 -1 扫到 1，找出图像「消失」的那个 k。")
say("3. 秩 1 的图像只有一列信息，那 rank(outer(v, u)) 一定等于 1 吗？")
say("   如果 v 或 u 里有零元呢？用 eigen_rank() 试。")
say("4. 2.5 节里把 INTER_LINEAR 换成 INTER_NEAREST，RMSE 应该变成多少？")
say("   换成 INTER_CUBIC 呢？自己跑一遍再解释。")

say("\nimg02 结束：乘法是复合，逆是回退，初等矩阵是行操作，秩是信息的维数。")
