# ===========================================================================
# 实验七 · 奇异值分解：一张图能压缩到多小
#
# 对应教材：第六章（二次型的标准形）+ 第五章（特征值），综合应用
# 建议学时：2
#
# 实验目的
#   1. 用 SVD 把一张图拆成「一层层从粗到细」的图；
#   2. 看到秩 k 近似：k 就是「保留几层」，误差由丢掉的奇异值决定；
#   3. 算清压缩比：为什么 k 小的时候真能省空间；
#   4. 彩色图 = 三个灰度图，各做各的。
#
# 运行：Rscript img07-svd.R   或在 RStudio 里 source()
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
say("灰度图 ", h, " × ", w, "，一共 ", h * w, " 个像素，也就这么多数据。")

# ===========================================================================
# 一、拆开它
# ===========================================================================
say("\n== 7.1 SVD：A = U Σ V' ==")

S <- eigen_svd(eigenMat(A))
sv <- as.vector(S$singular_values())
U <- as.matrix(S$u())
V <- as.matrix(S$v())
say("奇异值个数 = ", length(sv), "（= min(", h, ", ", w, ")）")
say("前 10 个：", paste(round(sv[1:10], 1), collapse = ", "))

# 这张测试图有大片纯色背景，很多行完全一样，于是矩阵不满秩：
# 尾部的奇异值是真 0，直接拿去算「最大 / 最小」只会得到 Inf。
tol <- 1e-9 * sv[1]
nz <- sv[sv > tol]
say("\n非零奇异值只有 ", length(nz), " 个（后面 ", length(sv) - length(nz),
    " 个精确为 0），所以 rank(A) = ", length(nz), " 而不是 ", length(sv), "。")
say("原因：大片纯色让许多行向量线性相关 —— 秩是「这张图到底有多少独立信息」。")
expect(eigen_rank(eigenMat(A)) == length(nz), "秩和零奇异值的个数对不上")
say("非零奇异值里最大 / 最小 = ", format(sv[1] / nz[length(nz)], digits = 4),
    "，这才是有效的条件数（拿 0 作分母只会得到 Inf）。")
say("能量占比：前 1 个奇异值占 ", round(100 * sv[1]^2 / sum(sv^2), 2),
    "%，前 10 个占 ", round(100 * sum(sv[1:10]^2) / sum(sv^2), 2), "%。")
say("第一层就吃掉九成能量：图的「大结构」比「细节」重要得多。")

# 单独看第 1 层：一个外积，秩一定是 1
layer1 <- sv[1] * U[, 1, drop = FALSE] %*% t(V[, 1, drop = FALSE])
say("\n第 1 层 u1·σ1·v1' 的秩 = ", eigen_rank(eigenMat(layer1)),
    "，它就是「整幅图最粗的那一层」。")
say("第 1 层的像素范围 ", paste(round(range(layer1), 1), collapse = " .. "),
    "：底色浅，所以最粗的一层整片都是正的，看着还像张图。")
for (i in 2:3) {
  Li <- sv[i] * U[, i, drop = FALSE] %*% t(V[, i, drop = FALSE])
  say("第 ", i, " 层的像素范围 ", paste(round(range(Li), 1), collapse = " .. "),
      "：已经带负号了。")
}
say("层本身不是图，只是矩阵；正负互相抵消，加起来才回到 0..255。")

# ===========================================================================
# 二、秩 k 近似
# ===========================================================================
say("\n== 7.2 秩 k 近似：只留下前 k 层 ==")

approx <- function(k) {
  U[, 1:k, drop = FALSE] %*% (sv[1:k] * t(V[, 1:k, drop = FALSE]))
}

ks <- c(1, 2, 5, 10, 20, 50, 100, length(nz), 240)
tab <- data.frame(k = ks, rmse = NA_real_, psnr = NA_real_, storage = NA_integer_)
for (i in seq_along(ks)) {
  k <- ks[i]
  ak <- il_clamp(approx(k))
  tab$rmse[i] <- il_rmse(ak, A)
  tab$psnr[i] <- il_psnr(ak, A)
  tab$storage[i] <- k * (h + w + 1)          # σ 一行 + u 的 k 列 + v 的 k 列
}
tab$ratio <- round(h * w / tab$storage, 2)
print(tab, row.names = FALSE)

say("\n记法：秩 k 近似要存 k 个数（σ）+ k×", h, "（u）+ k×", w, "（v），")
say("     一共有 ", h, "×", w, " = ", h * w, " 个数。k 小的时候确实省。")
say("表里 k 从 ", length(nz), " 加到 240 一个像素都没好 —— 零奇异值加了也白加，")
say("rank(A) 已经是「能存下的全部信息」了。")

# Eckart-Young：误差的平方 = 被丢掉的那些奇异值的平方和
k <- 20
ak <- approx(k)
ey <- sqrt(sum(sv[(k + 1):length(sv)]^2))
say("\nEckart-Young 定理：k = ", k, " 时的 RMSE = sqrt(Σ_{i>k} σ_i^2 / (m·n))")
say("  实测（不截断）", round(il_rmse(ak, A), 6),
    "，理论值 ", round(ey / sqrt(h * w), 6), " —— 完全相等")
expect_near(il_rmse(ak, A), ey / sqrt(h * w), 1e-12, "Eckart-Young 不成立")
say("  截断到 0..255 之后变成 ", round(il_rmse(il_clamp(ak), A), 6),
    "（差值来自那些本来会被截掉的像素）")

il_fig_sheet("img07-rank",
             list(A, il_clamp(approx(1)), il_clamp(approx(5)), il_clamp(approx(20)),
                  il_clamp(approx(50)), il_clamp(approx(120))),
             c("original", "k = 1", "k = 5", "k = 20", "k = 50",
               sprintf("k = 120 (%.1f dB)", il_psnr(il_clamp(approx(120)), A))),
             ncol = 3)

# ===========================================================================
# 三、压缩比与视觉质量
# ===========================================================================
say("\n== 7.3 压缩比和「看起来还行」的平衡点 ==")

say("把「可接受」定成 PSNR ≥ 30 dB（经验上肉眼不太看得出），找最小的 k：")
ok <- tab$psnr >= 30
if (any(ok)) {
  kk <- min(tab$k[ok])
  say("  k = ", kk, " 时 PSNR = ", round(tab$psnr[tab$k == kk], 2),
      " dB，压缩比 ", tab$ratio[tab$k == kk], " : 1")
} else {
  say("  这份数据里所有 k 都到不了 30 dB")
}
say("k 再大，PSNR 涨得很慢，但存储线性增长 —— 这就是「边际收益递减」。")

# ===========================================================================
# 四、彩色图：三个通道各做一次
# ===========================================================================
say("\n== 7.4 彩色图：每个通道各拆一次 ==")

chan <- lapply(1:3, function(ch) {
  m <- img[, , ch]
  s <- eigen_svd(eigenMat(m))
  v <- as.vector(s$singular_values())
  uu <- as.matrix(s$u()); vv <- as.matrix(s$v())
  k <- 30
  list(k = k, sv = v, approx = uu[, 1:k, drop = FALSE] %*% (v[1:k] * t(vv[, 1:k, drop = FALSE])))
})
col_approx <- il_rgb(il_clamp(chan[[1]]$approx), il_clamp(chan[[2]]$approx),
                     il_clamp(chan[[3]]$approx))
say("三个通道的奇异值谱几乎一样（图像是灰调的）：")
for (ch in 1:3) say("  通道 ", ch, " 前 3 个奇异值：",
                    paste(round(chan[[ch]]$sv[1:3], 1), collapse = ", "))
say("k = 30 的彩色近似：RMSE = ", round(il_rmse(col_approx, img), 3),
    "，PSNR = ", round(il_psnr(col_approx, img), 2), " dB")
say("存储：3 × 30 × (", h, " + ", w, " + 1) = ", 3 * 30 * (h + w + 1),
    " 个数，原始 3 × ", h * w, " 个数，压缩比 ",
    round(3 * h * w / (3 * 30 * (h + w + 1)), 2), " : 1")

il_fig_sheet("img07-colour", list(img, col_approx),
             c("original", "k = 30 per channel"), ncol = 2)

# ===========================================================================
# 五、和「按行列砍一半」比一比
# ===========================================================================
say("\n== 7.5 对比：直接缩小一半 vs SVD 压缩 ==")

small <- as_array(cv_resize(il_cv(A), NULL, 0.5, 0.5), numeric = TRUE)
back <- as_array(cv_resize(il_cv(small), NULL, 2, 2), numeric = TRUE)
say("把图缩到一半再放大回来（相当于只存 25% 的像素）：")
say("  RMSE = ", round(il_rmse(back, A), 2), "，PSNR = ", round(il_psnr(back, A), 2), " dB")
say("SVD 的 k = 20（存 ", round(100 * 20 * (h + w + 1) / (h * w), 1),
    "% 的数）：RMSE = ", round(il_rmse(il_clamp(approx(20)), A), 2),
    "，PSNR = ", round(il_psnr(il_clamp(approx(20)), A), 2), " dB")
say("两者各有胜负：缩放会整体糊掉，SVD 保留大结构、丢掉细节（线条会断）。")

# ===========================================================================
# 六、实验习题
# ===========================================================================
say("\n== 7.6 习题 ==")
say("1. 把测试图换成纯棋盘（il_checker()），奇异值谱会变成什么样？为什么？")
say("2. k 从 1 加到 240，RMSE 单调下降吗？PSNR 呢？画出来。")
say("3. 用 SVD 重构时，为什么第 1 层看起来像「一团糊」？它保留了哪些信息？")
say("   （提示：u1 和 v1 各自只有 ", h, " 和 ", w, " 个数。）")
say("4. k 设成 rank(A) = ", length(nz), " 时重构误差只剩 ", format(il_rmse(approx(length(nz)), A), digits = 3),
    "，为什么不是精确的 0？（提示：浮点数只有 53 位有效二进制位。）")

say("\nimg07 结束：SVD 把图分成「从粗到细的层」，秩就是「留几层」。")
