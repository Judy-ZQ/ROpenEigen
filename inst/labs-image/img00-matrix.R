# ===========================================================================
# 实验零 · 图像就是矩阵
#
# 对应教材：全书通用（做后面所有图像实验之前先做这个）
# 建议学时：1
#
# 实验目的
#   1. 看清一张图在 R 里到底是几个矩阵，维度、取值范围、存的是什么；
#   2. 把「灰度化」认成一次线性组合 —— 三个通道的加权和；
#   3. 会做逐元素（逐像素）运算，并知道为什么必须截断；
#   4. 会在 R 数组和 cvMat 之间来回转，会显示、会存盘。
#
# 运行：Rscript img00-matrix.R   或在 RStudio 里 source()
#       （图默认存到当前目录的 imgfigs/，可用 LINALG_FIGS 改）
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

say("ROpenCV ", as.character(packageVersion("ROpenCV")),
    " | ROpenEigen ", as.character(packageVersion("ROpenEigen")))
say("图会存到：", il_outdir())

# ===========================================================================
# 一、一张图是三个矩阵
# ===========================================================================
say("\n== 0.1 一张彩色图 = 三个矩阵 ==")

img <- il_read()                       # 没有 LINALG_IMG 时用 il_chart() 画一张
say("尺寸 (高, 宽, 通道)：", paste(dim(img), collapse = " × "))
say("取值：", paste(range(img), collapse = " .. "), "（8 位无符号，0 = 黑，255 = 白）")

# 通道分离：手写一遍，再用库函数一遍，两者必须一样
R <- img[, , 1]; G <- img[, , 2]; B <- img[, , 3]
sp <- cv_split(il_cv(img))          # il_cv(): R 数组 -> 8 位 cvMat（见 0.4）
say("手写分离的 R 通道均值：", round(mean(R), 2))
say("cv_split() 第 1 个通道均值：", round(mean(as_array(sp[[1]], numeric = TRUE)), 2))
expect_near(mean(R), mean(as_array(sp[[1]], numeric = TRUE)), 1e-12,
            "手写分离与 cv_split 不一致")

# 挑一块「有变化」的地方打印 —— 全是常数的块看不出矩阵的味道。
say("\n找一块对比最剧烈的地方，打印 R 通道 6×8 的像素值（就是教材里那种矩阵）：")
# 在粗网格上找标准差最大的一块 6×8 窗口（数字变得最热闹的那种）
# 判据是「不同取值的个数」：边缘窗口只有两种值，虽然方差大但不好看；
# 灰度过渡带上一块窗口能有十几种值，那才像教材里的例子。
best <- c(-1, -1, 1L, 1L)
for (ii in seq(1L, nrow(R) - 6L, by = 6L)) {
  for (jj in seq(1L, ncol(R) - 8L, by = 6L)) {
    w <- as.vector(R[ii:(ii + 5), jj:(jj + 7)])
    nu <- length(unique(round(w)))
    if (nu > best[1] || (nu == best[1] && stats::sd(w) > best[2]))
      best <- c(nu, stats::sd(w), ii, jj)
  }
}
i0 <- best[3]; j0 <- best[4]
say("  第 ", i0, "–", i0 + 5, " 行、第 ", j0, "–", j0 + 7, " 列：")
print(round(R[i0:(i0 + 5), j0:(j0 + 7)]))

say("\n每一个通道单独拿出来都是一张灰度图：")
il_fig_sheet("img00-channels", list(img, R, G, B),
             c("original (RGB)", "R channel", "G channel", "B channel"), ncol = 4)

# ===========================================================================
# 二、灰度 = 三个矩阵的线性组合
# ===========================================================================
say("\n== 0.2 灰度化是一次线性组合 ==")

w <- c(0.299, 0.587, 0.114)            # ITU-R BT.601
gray_hand <- w[1] * R + w[2] * G + w[3] * B
gray_cv <- as_array(cv_cvtColor(il_cv(img), COLOR_RGB2GRAY), numeric = TRUE)

say("权重 w = (", paste(w, collapse = ", "), ")，且 sum(w) = ", sum(w))
say("手算灰度与 cv_cvtColor 的最大差：", max(abs(gray_hand - gray_cv)), "（差 1 以内是四舍五入）")
expect(max(abs(gray_hand - gray_cv)) <= 1.001, "两种灰度化差别超过 1")

# 这不是巧合：灰度化就是「把三个通道当坐标，做内积」。
# 换句话说，灰度图 = 彩色图在方向 w 上的投影。
say("把 (R,G,B) 看成三维向量，灰度就是它在 w 方向上的投影：")
spread <- pmax(R, G, B) - pmin(R, G, B)      # 找一个三通道差别最大的像素
kk <- which(spread == max(spread), arr.ind = TRUE)[1, ]
i0 <- kk[1]; j0 <- kk[2]
v <- c(R[i0, j0], G[i0, j0], B[i0, j0])
say("  像素 (", i0, ",", j0, ") 的颜色向量 = (", paste(v, collapse = ", "), ")")
say("  与 w 的内积 = ", round(sum(v * w), 3),
    "，灰度图上的值 = ", round(gray_hand[i0, j0], 3))
expect_near(sum(v * w), gray_hand[i0, j0], 1e-12, "灰度不是内积")

il_fig_sheet("img00-grey",
             list(gray_hand, 255 * (abs(gray_hand - gray_cv) > 0.5)),
             c("grey = 0.299R + 0.587G + 0.114B", "white = pixels where the two differ"), ncol = 2)

# ===========================================================================
# 三、逐像素运算：和矩阵的加法、数乘是同一件事
# ===========================================================================
say("\n== 0.3 变亮 / 提高对比度 / 反色 / gamma ==")

clamp <- function(m) pmin(pmax(m, 0), 255)
brighter <- clamp(gray_hand + 40)          # 加一个标量矩阵 40·J
contrast <- clamp(1.6 * (gray_hand - 128) + 128)   # 先减均值再放大再移回
invert <- 255 - gray_hand                  # 相当于 255·J - G
gamma <- 255 * (gray_hand / 255)^0.45      # 非线性，逐元素

say("变亮 +40：均值 ", round(mean(gray_hand), 2), " -> ", round(mean(brighter), 2),
    "（被截断的像素 ", round(100 * mean(gray_hand + 40 > 255), 1), "%）")
say("反色：两次反色回到原图，最大差 = ", format(max(abs(255 - invert - gray_hand)), scientific = TRUE))
expect(max(abs(255 - invert - gray_hand)) < 1e-12, "反色不是一一对应")
say("gamma 0.45：均值 ", round(mean(gamma), 2), "（提亮暗部，但不线性）")

il_fig_sheet("img00-pointwise", list(gray_hand, brighter, contrast, invert, gamma),
             c("original", "+40 (clipped)", "contrast x1.6", "255 - x", "gamma 0.45"), ncol = 3)

# ===========================================================================
# 四、R 数组 <-> cvMat
# ===========================================================================
say("\n== 0.4 两种对象怎么来回转，以及 depth 这件事 ==")

say("R 数组：", paste(dim(gray_hand), collapse = " × "), "，类型 ", typeof(gray_hand))

md <- Mat(gray_hand)
say("Mat(双精度矩阵)  -> depth = ", mat_info(md)[5], "（CV_64F）")
mr <- Mat({ r <- as.raw(round(gray_hand)); dim(r) <- dim(gray_hand); r })
say("Mat(raw 数组)     -> depth = ", mat_info(mr)[5], "（CV_8U，也就是图片的存法）")
say("il_cv(矩阵)       -> depth = ", mat_info(il_cv(gray_hand))[5], "（同一件事，包在工具里）")

say("这不是小细节：cvtColor / warpAffine / GaussianBlur 都只吃 8U、16U、32F，")
say("喂它们一个 CV_64F 会直接报 'Unsupported depth of input image'。")
say("所以本套实验的规矩是：算的时候用 R 的 double，交给 OpenCV 之前走 il_cv()。")

back <- as_array(il_cv(gray_hand), numeric = TRUE)
say("8 位往返（double -> 8U -> 整数）最大差：", max(abs(back - round(gray_hand))),
    "（取整带来的，不是 bug）")

ok <- cv_imwrite(file.path(il_outdir(), "img00-grey.png"), il_cv(gray_hand))
say("cv_imwrite() = ", ok, "，文件存在：",
    file.exists(file.path(il_outdir(), "img00-grey.png")))
say("读回来用 cv_imread()，它给的是 BGR —— 单通道图无所谓，彩色图记得翻通道。")

# ===========================================================================
# 五、实验习题
# ===========================================================================
say("\n== 0.5 习题 ==")
say("1. 用 w = (1/3, 1/3, 1/3) 做一次灰度化，和 0.299/0.587/0.114 的结果差多少？")
say("   为什么人眼对绿色最敏感这件事，会体现在权重上？")
say("2. 把 contrast 的 1.6 改成 3，调亮和调暗的两端会怎样？用直方图解释。")
say("3. 写一个函数，把图像左右翻转（不许用 cv_flip），再用它验证「翻转两次 = 恒等」。")
say("4. Mat() 一个 double 矩阵得到 CV_64F，直接 plot() 会报什么错？为什么？")
say("   （提示：as.raster 按 depth 决定归一化的分母。）")

say("\nimg00 结束：图像是矩阵，通道是矩阵的层，灰度是线性组合。")
