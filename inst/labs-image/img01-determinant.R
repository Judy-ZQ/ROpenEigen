# ===========================================================================
# 实验一 · 行列式：面积、翻转、塌陷
#
# 对应教材：第一章 §1.1–§1.5，重点是 §1.5 行列式的应用
# 建议学时：2
#
# 实验目的
#   1. 看到「行列式 = 面积（体积）的缩放倍数」这句话长什么样；
#   2. 看到 det 的符号为什么叫「定向」—— 负号就是照镜子；
#   3. 看到 det = 0 时图像被压成一条线，以及为什么这时矩阵不可逆；
#   4. 把 3×3 行列式用起来：RGB 三个通道就是三维向量，3×3 矩阵就是调色。
#
# 运行：Rscript img01-determinant.R   或在 RStudio 里 source()
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

# 多边形面积（鞋带公式）：教材 §1.5 用的就是它。
poly_area <- function(P) {
  n <- ncol(P)
  abs(sum(P[1, ] * P[2, c(2:n, 1)] - P[1, c(2:n, 1)] * P[2, ])) / 2
}

img <- il_read()
gray <- il_gray(img)
say("测试图：", paste(dim(gray), collapse = " × "), "（高 × 宽）")
say("图会存到：", il_outdir())

# ===========================================================================
# 一、det 就是面积缩放倍数
# ===========================================================================
say("\n== 1.1 |det A| = 面积缩放倍数 ==")

unit <- rbind(c(0, 1, 1, 0), c(0, 0, 1, 1))     # 单位正方形的四个顶点
A <- matrix(c(2, 0.5,
              0, 1), 2, byrow = TRUE)

say("A = [[2, 0.5], [0, 1]]，det A = ", eigen_det(eigenMat(A)), "（手算 2×1 - 0.5×0 = 2）")
P <- A %*% unit
say("单位正方形 -> 平行四边形，顶点：")
print(round(t(P), 3))
say("鞋带公式算面积 = ", poly_area(P), "，|det A| = ", abs(eigen_det(eigenMat(A))))
expect_near(poly_area(P), abs(eigen_det(eigenMat(A))), 1e-12,
            "面积不等于 |det|")

# 用像素再量一次：一张白方块，变换后数亮像素。
sq <- matrix(0, 240, 240)
sq[81:160, 81:160] <- 255                     # 80×80 = 6400 个白像素
warped <- il_warp_inv(sq, A, out = c(300, 360))
n_white <- sum(warped > 128)
say("\n白方块 80×80 = ", 80 * 80, " 个像素，变换后亮像素 ", n_white,
    "，比值 = ", round(n_white / 6400, 4))
say("理论值 |det A| = 2，实测比值 ", round(n_white / 6400, 4),
    "（差的是边界上一个像素的取整）")
expect(abs(n_white / 6400 - 2) < 0.02, "像素测出来的面积比不对")

il_fig_sheet("img01-area",
             list(sq, warped, gray, il_warp_inv(gray, A, out = c(240L, 660L))),
             c("white square 80x80", "after A: area x2 (pixel count)",
               "the test image", "A applied to the image"), ncol = 4)

# ===========================================================================
# 二、符号 = 定向：负号就是照镜子
# ===========================================================================
say("\n== 1.2 det < 0 就是镜像 ==")

Mflip <- matrix(c(1, 0,
                  0, -1), 2, byrow = TRUE)
say("M = [[1, 0], [0, -1]]：上下翻转，det M = ", eigen_det(eigenMat(Mflip)))
say("面积不变（|det| = 1），但定向反了。图里的字母 A 会变成倒的：")

il_fig_sheet("img01-flip",
             list(gray, il_warp_inv(gray, Mflip, out = dim(gray)),
                  il_warp_inv(gray, matrix(c(-1, 0, 0, 1), 2, byrow = TRUE), out = dim(gray))),
             c("original", "flip vertically: det = -1", "flip horizontally: det = -1"), ncol = 3)

# 先上下翻、再左右翻，两个 det = -1 相乘得 +1：图像转 180 度，与直接转 180 度一致。
Mrot <- matrix(c(-1, 0,
                 0, -1), 2, byrow = TRUE)
twice <- il_warp_inv(il_warp_inv(gray, Mflip, out = dim(gray)),
                     matrix(c(-1, 0, 0, 1), 2, byrow = TRUE), out = dim(gray))
say("两次镜像（det 相乘 = +1）应当等于旋转 180 度，实测差：",
    round(il_rmse(twice, il_warp_inv(gray, Mrot, out = dim(gray))), 6))

# ===========================================================================
# 三、det = 0：图像被压成一条线
# ===========================================================================
say("\n== 1.3 det = 0：塌陷 ==")

S <- matrix(c(1, 1,
              1, 1), 2, byrow = TRUE)
say("S = [[1, 1], [1, 1]]，det S = ", eigen_det(eigenMat(S)),
    "，R(S) = ", eigen_rank(eigenMat(S)))
say("S 把平面上的每个点都送到直线 y = x 上 —— 二维塌成一维，")
say("所以不存在 S^{-1}：信息丢了，回不去。")

col <- il_warp_fwd(gray, S, out = dim(gray))
say("正映射：亮像素占 ", round(100 * mean(col > 0), 2), "%（大部分格子空着，因为像素被搬到一条线上）")
say("ROpenEigen 的 solve() 对奇异矩阵不报错，它返回：")
inv_try <- as.matrix(solve(eigenMat(S)))
say("  第一列 = (", paste(round(inv_try[, 1], 3), collapse = ", "), ") —— 全是 NaN/Inf，这就是「算不出来」的样子。")
say("（base R 的 solve() 会直接报错；两者行为不同，实验二的 lab02 里讨论过。）")

il_fig_sheet("img01-collapse",
             list(gray, col, il_warp_fwd(gray, matrix(c(1, 2, 2, 4), 2, byrow = TRUE), out = dim(gray))),
             c("original", "S = [[1,1],[1,1]]: collapses to a line",
               "T = [[1,2],[2,4]]: the same, rank 1"), ncol = 3)

# ===========================================================================
# 四、3×3 行列式：颜色就是三维向量
# ===========================================================================
say("\n== 1.4 3×3：RGB 是三维向量，矩阵就是调色 ==")

# 3×3 矩阵作用在「每个像素的颜色向量」上。这是最便宜的一次 3 维线性变换演示：
# 你看到的颜色变化，就是三个通道被重新线性组合的结果。
apply_color <- function(img, M) {
  h <- dim(img)[1]; w <- dim(img)[2]
  v <- matrix(as.vector(img), ncol = 3)          # 每行一个 RGB 向量
  out <- pmin(pmax(v %*% t(M), 0), 255)          # M v，再截断到 0..255
  array(as.vector(out), c(h, w, 3))
}

swapRB <- matrix(c(0, 0, 1,
                   0, 1, 0,
                   1, 0, 0), 3, byrow = TRUE)
say("交换 R 与 B：det = ", eigen_det(eigenMat(swapRB)),
    "（一个换位，负号 = 手性反了，但颜色体积不变）")

# 绕「灰度轴」(1,1,1) 旋转 120 度：cos = -1/2，sin = sqrt(3)/2
ct <- cos(2 * pi / 3); st <- sin(2 * pi / 3)
J <- matrix(1, 3, 3)
K <- matrix(c(0, -1, 1,
              1, 0, -1,
              -1, 1, 0), 3, byrow = TRUE) / sqrt(3)
hue <- ct * diag(3) + (1 - ct) * J / 3 + st * K
say("绕灰度轴转 120 度：det = ", round(eigen_det(eigenMat(hue)), 10),
    "（旋转，所以 det = 1，体积守恒）")
say("灰度轴不动：(1,1,1) 被映到 ", paste(round(hue %*% c(1, 1, 1), 10), collapse = ", "),
    " —— 灰色的还是灰色，只是红绿蓝换了位置。")

img_swap <- apply_color(img, swapRB)
img_hue <- apply_color(img, hue)
say("\n橙色那个像素 (224, 123, 57)：")
say("  交换 R/B -> (", paste(round(img_swap[183, 148, ]), collapse = ", "),
    ")，通道和 ", round(sum(img_swap[183, 148, ]), 1), "（每行和为 1，所以守恒）")
say("  转 120 度 -> (", paste(round(img_hue[183, 148, ]), collapse = ", "), ")")

# 秩为 1 的 3×3：三个通道都变成同一个灰度值 —— 彩色图变成灰度图
flat <- matrix(1 / 3, 3, 3)
say("\n[[1/3,1/3,1/3] 三行相同] 的秩 = ", eigen_rank(eigenMat(flat)),
    "，det = ", round(eigen_det(eigenMat(flat)), 12), "：它把颜色投影到灰度轴上。")
img_flat <- apply_color(img, flat)
say("  用它作用一次，图就变成灰的了（三个通道相等）。")

il_fig_sheet("img01-color",
             list(img, img_swap, img_hue, img_flat),
             c("original", "swap R and B: det = -1", "hue rotation: det = 1",
               "rank-1 projection: grayscale"), ncol = 4)

# ===========================================================================
# 五、共线判定：面积为零
# ===========================================================================
say("\n== 1.5 三个点共线 ⇔ 三角形面积为 0 ==")

tri <- function(p1, p2, p3) {
  M <- cbind(c(p1[1], p2[1], p3[1]), c(p1[2], p2[2], p3[2]), c(1, 1, 1))
  abs(eigen_det(eigenMat(M)))
}
say("(0,0), (2,1), (4,2) 共线：行列式 = ", tri(c(0, 0), c(2, 1), c(4, 2)), "，面积 = ", tri(c(0,0), c(2,1), c(4,2)) / 2)
say("(0,0), (2,1), (1,3) 不共线：行列式 = ", tri(c(0, 0), c(2, 1), c(1, 3)), "，面积 = ", tri(c(0,0), c(2,1), c(1,3)) / 2)
say("行列式在这里就是「两倍面积」：2·2.5 = 5。")

# ===========================================================================
# 六、实验习题
# ===========================================================================
say("\n== 1.6 习题 ==")
say("1. 取 A = [[3, 1], [1, 2]]，先算 det，再用像素数法量面积，两者差多少？")
say("   把白方块换成 160×160，误差会变大还是变小？为什么？")
say("2. det = -1 的矩阵一定是对称的镜像吗？举一个 det = -1 但不是镜像的例子，")
say("   用图像验证它确实把定向反了（提示：错切 [[1,1],[0,-1]]）。")
say("3. 找一个 3×3 的颜色矩阵，使 det = 0。图像会变成什么样？为什么？")
say("4. 用 il_warp_fwd() 做 rank 2 的 3×3 变换不现实，但 3×3 作用在颜色上是")
say("   rank 3 之外的例子：把 R、G、B 都换成 (R+G+B)/3，矩阵的秩是多少？")

say("\nimg01 结束：行列式是面积，符号是定向，零是塌陷。")
