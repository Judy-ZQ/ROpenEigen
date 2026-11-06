# ===========================================================================
# 实验六 · 二次型：椭圆、双曲线和「正定」到底在说什么
#
# 对应教材：第六章 §6.1–§6.3（二次型、标准形、正定）
# 建议学时：2
#
# 实验目的
#   1. 把 f(x) = x'Ax 画出来：正定是椭圆，不定是双曲线，半定是两条平行线；
#   2. 看到椭圆的两个半轴长就是 1/sqrt(特征值)，方向就是特征向量；
#   3. 用二次型「造」一张图（高斯斑点），A 一变形状就变；
#   4. 把惯性指数和正定判别连起来。
#
# 运行：Rscript img06-quadratic.R   或在 RStudio 里 source()
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

h <- 240; w <- 320
# 网格：以画布中心为原点，坐标范围 [-1.6, 1.6] × [-1.2, 1.2]
gx <- seq(-1.6, 1.6, length.out = w)
gy <- seq(-1.2, 1.2, length.out = h)
GX <- matrix(rep(gx, each = h), h, w)
GY <- matrix(rep(gy, times = w), h, w)

# 二次型 q(x) = x'Ax（A 对称）
qform <- function(A) A[1, 1] * GX^2 + 2 * A[1, 2] * GX * GY + A[2, 2] * GY^2

# 画等值线 q = 1：极坐标下 r = 1/sqrt(q(cosθ, sinθ))，一步到位。
# 某些方向上 q <= 0（不定型），那就没有交点 —— 于是画出两支双曲线。
level_curve <- function(img, A, val = c(255, 0, 0)) {
  th <- seq(0, 2 * pi, length.out = 720)
  a <- A[1, 1]; b <- A[1, 2]; c <- A[2, 2]
  qq <- a * cos(th)^2 + 2 * b * cos(th) * sin(th) + c * sin(th)^2
  ok <- qq > 1e-12
  r <- rep(NA_real_, length(th)); r[ok] <- 1 / sqrt(qq[ok])
  px <- (r * cos(th)) / (gx[2] - gx[1]) + (w + 1) / 2
  py <- -(r * sin(th)) / (gy[2] - gy[1]) + (h + 1) / 2
  for (i in seq_len(length(th) - 1)) {
    if (!is.na(px[i]) && !is.na(px[i + 1]))
      img <- il_draw(img, px[i], py[i], px[i + 1], py[i + 1], val = val, width = 2)
  }
  img
}

APD <- matrix(c(1.4, 0.5, 0.5, 1.0), 2, byrow = TRUE)     # 正定
AID <- matrix(c(1.0, 0.6, 0.6, -1.0), 2, byrow = TRUE)    # 不定
ASD <- matrix(c(1.0, 1.0, 1.0, 1.0), 2, byrow = TRUE)     # 半定（秩 1）

say("三个 2×2 对称矩阵，看它们的二次型")
for (nm in c("APD", "AID", "ASD")) {
  A <- get(nm)
  l <- as.vector(eigen_values(eigenMat(A)))
  say("  ", nm, " = [[", A[1, 1], ", ", A[1, 2], "], [", A[2, 1], ", ", A[2, 2], "]]",
      "：特征值 (", paste(round(l, 3), collapse = ", "), ")",
      "，det = ", round(eigen_det(eigenMat(A)), 3),
      "，惯性指数 = (", sum(l > 1e-9), ", ", sum(l < -1e-9), ")")
}

# 背景画成 f 的灰度（暗 = 值大），再叠上 f = 1 的红线
# 背景做成三通道，这样红线才画得上去（二维矩阵只能存灰度）
bg <- function(A) {
  v <- 255 - il_disp(log1p(abs(qform(A))))
  il_rgb(v, v, v)
}
p_ell <- level_curve(bg(APD), APD)
p_hyp <- level_curve(bg(AID), AID)
p_par <- level_curve(bg(ASD), ASD)

il_fig_sheet("img06-levels", list(p_ell, p_hyp, p_par),
             c("positive definite: ellipse", "indefinite: hyperbola",
               "semi-definite (rank 1): two parallel lines"), ncol = 3)

# ===========================================================================
# 一、半轴长 = 1/sqrt(特征值)
# ===========================================================================
say("\n== 6.1 椭圆的半轴：1 / sqrt(特征值) ==")

lam <- as.vector(eigen_values(eigenMat(APD)))
V <- as.matrix(eigen_vectors(eigenMat(APD)))
say("APD 的特征值 ", paste(round(lam, 4), collapse = ", "))
say("对应的半轴长 1/sqrt(λ) = ", paste(round(1 / sqrt(lam), 4), collapse = ", "))
say("面积 = π·a·b = π/sqrt(det) = ", round(pi / sqrt(eigen_det(eigenMat(APD))), 4))

# 数值验证：沿特征方向走到 q = 1 的距离
for (k in 1:2) {
  v <- V[, k]
  tt <- 1 / sqrt(as.numeric(t(v) %*% APD %*% v))
  say("  沿 v", k, " 方向：t = ", round(tt, 4), "，代入 q = ",
      round(as.numeric(tt^2 * t(v) %*% APD %*% v), 10))
}

# ===========================================================================
# 二、用二次型「造」一张图
# ===========================================================================
say("\n== 6.2 用 exp(-x'Ax) 造斑点：A 的形状就是斑点的形状 ==")

blob <- function(A, s = 1) {
  z <- exp(-qform(A) / s)
  255 * il_disp(z) / 255
}
# 四张：正定（细长）、正定转 90 度、正定（圆）、不定（马鞍）
A1 <- matrix(c(3, 0.8, 0.8, 0.7), 2, byrow = TRUE)
A2 <- matrix(c(0.7, 0.8, 0.8, 3), 2, byrow = TRUE)
A3 <- matrix(c(1.5, 0, 0, 1.5), 2, byrow = TRUE)
imgs <- list(blob(A1), blob(A2), blob(A3), il_disp(abs(1 - 2 * exp(-qform(AID)))))
say("A1 = [[3, 0.8], [0.8, 0.7]]：斑点又细又斜，主轴 = ",
    round(atan2(as.matrix(eigen_vectors(eigenMat(A1)))[2, 2],
                as.matrix(eigen_vectors(eigenMat(A1)))[1, 2]) * 180 / pi, 2), " 度")
say("  注意方向是反的：半轴长 = 1/sqrt(λ)，所以特征值大的方向斑点反而窄。")
say("A2 把 3 和 0.7 对调，斑点就转成另一个方向 —— 主对角元控制「摊得多开」。")
say("A3 = 1.5·I 是圆：两个特征值相等，方向没区别。")

il_fig_sheet("img06-blobs", imgs,
             c("A1: elongated", "A2: mirrored", "A3: circular (equal eigenvalues)",
               "indefinite: saddle"), ncol = 4)

# ===========================================================================
# 三、正定性判别：四种方法互相印证
# ===========================================================================
say("\n== 6.3 正定的四种判别法 ==")

minors <- function(A) {
  c(A[1, 1], A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1])
}
classify <- function(A) {
  l <- as.vector(eigen_values(eigenMat(A)))
  if (all(l > 1e-9)) "正定" else if (all(l < -1e-9)) "负定"
  else if (any(abs(l) <= 1e-9) && all(l >= -1e-9)) "半正定"
  else "不定"
}
for (nm in c("APD", "AID", "ASD", "A3")) {
  A <- get(nm)
  dm <- minors(A)
  say("  ", nm, "：顺序主子式 (", paste(round(dm, 3), collapse = ", "), ")",
      "；特征值 (", paste(round(as.vector(eigen_values(eigenMat(A))), 3), collapse = ", "),
      ")；Cholesky ", if (isTRUE(tryCatch({eigen_chol(eigenMat(A)); TRUE},
                                          error = function(e) FALSE))) "成功" else "失败",
      "；判定 ", classify(A))
}
say("\n四种说法是同一件事：顺序主子式全正 ⇔ 特征值全正 ⇔ Cholesky 分解成功 ⇔ 惯性指数 (2, 0)。")

# ===========================================================================
# 四、参数扫描：从椭圆到双曲线
# ===========================================================================
say("\n== 6.4 让一个参数走过去，看椭圆变成双曲线 ==")

# A(a) = [[a, 1], [1, a]]：det = a^2 - 1，a 从 1.6 降到 0.6
for (a in c(1.6, 1.1, 1.0, 0.8)) {
  Aa <- matrix(c(a, 1, 1, a), 2, byrow = TRUE)
  l <- as.vector(eigen_values(eigenMat(Aa)))
  say("  a = ", a, "：det = ", round(eigen_det(eigenMat(Aa)), 4),
      "，特征值 (", paste(round(l, 3), collapse = ", "), ") -> ", classify(Aa))
}
say("临界点 a = 1：det = 0，椭圆退化成两条平行线 —— 这正是「半正定但不正定」。")

# ===========================================================================
# 五、实验习题
# ===========================================================================
say("\n== 6.5 习题 ==")
say("1. 取 A = [[2, 1], [1, 2]]，画出 x'Ax = 1 的椭圆，量出两个半轴，")
say("   和 1/sqrt(特征值) 对一下。")
say("2. 把 6.2 的 blob() 换成 exp(-(x'Ax)^2)，斑点形状会变吗？为什么？")
say("3. 设计一个 A，使 x'Ax = 1 是一条真正的双曲线（两支），并算出它的渐近线")
say("   方向。（提示：q(方向) = 0 的方向就是渐近线。）")
say("4. 二次型 x'Ax 在约束 ||x|| = 1 下的最大值就是最大特征值（Rayleigh 商）。")
say("   用随机方向采样验证这一点，再说明它和主轴方向的关系。")

say("\nimg06 结束：二次型把矩阵变成了形状，特征值就是这把尺子。")
