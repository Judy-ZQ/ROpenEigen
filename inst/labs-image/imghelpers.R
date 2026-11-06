# ===========================================================================
# imghelpers.R —— 图像实验的公共工具
#
# 这套实验只做一件事：把「矩阵」画出来给你看。
# 所以工具也分三类：
#   1. 造图 / 读图：把一张图变成 R 里的矩阵（或三维数组）；
#   2. 画图：把矩阵画回屏幕，或拼成一张对照图存成 PNG；
#   3. 断言：和 inst/labs 一样的 expect_near()，数字必须自洽。
#
# 约定
#   * 彩色图统一用 RGB 顺序的 (h, w, 3) 数组，取值 0..255 的 double。
#     ROpenCV 的 imread() 给的是 BGR，读进来立刻翻转一次，之后不再管。
#   * 灰度图是 (h, w) 矩阵，取值 0..255 的 double。
#   * 图里的标题一律用 ASCII：Windows 的默认设备画中文会掉字。
# ===========================================================================

## ---------------------------------------------------------------- 输出与断言

# 图存哪里：默认当前目录下的 imgfigs/，可用环境变量 LINALG_FIGS 改。
il_outdir <- function() {
  d <- Sys.getenv("LINALG_FIGS", unset = "")
  if (!nzchar(d)) d <- file.path(getwd(), "imgfigs")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

say <- function(...) cat(..., "\n", sep = "")

# 打开一张待保存的 PNG。用完必须 il_done() 关掉。
il_fig <- function(name, width = 1000, height = 680, res = 110) {
  p <- file.path(il_outdir(), paste0(name, ".png"))
  grDevices::png(p, width = width, height = height, res = res)
  p
}

il_done <- function(p) {
  grDevices::dev.off()
  say("  图已保存：", p)
  invisible(p)
}

relerr <- function(a, b) {
  s <- max(abs(b))
  if (s == 0) s <- 1
  max(abs(a - b)) / s
}

expect <- function(ok, what) {
  if (!isTRUE(ok)) stop("断言失败：", what, call. = FALSE)
  invisible(TRUE)
}

expect_near <- function(a, b, tol = 1e-8, what = "") {
  e <- relerr(a, b)
  if (!is.finite(e) || e > tol)
    stop(sprintf("断言失败：%s（相对误差 %.3e > %.1e）", what, e, tol), call. = FALSE)
  invisible(TRUE)
}

## ---------------------------------------------------------------- 造图与读图

# 棋盘：m 个格子，每格 tile 像素。明暗两级的对比很硬，适合量面积、看形变。
il_checker <- function(m = 8L, tile = 24L, lo = 40, hi = 210) {
  n <- m * tile
  r <- ((row(matrix(0, n, n)) - 1L) %/% tile) + ((col(matrix(0, n, n)) - 1L) %/% tile)
  ifelse(r %% 2L == 0L, hi, lo)
}

# 用 R 的绘图设备画一张「测试图」，再用 ROpenCV 读回来。
# 这一步本身就是一条链：图形 -> PNG -> 解码 -> 矩阵。
il_chart <- function(w = 320L, h = 240L) {
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = w, height = h, bg = "white")
  op <- graphics::par(no.readonly = TRUE)
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::plot.window(c(0, w), c(0, h), xaxs = "i", yaxs = "i")
  # 背景用中间调：左边一块略暗的带，右边留白。中间调的好处是后面做
  # 「加 40」「乘 1.6」这类逐像素运算时，不会一动就全局截断。
  graphics::rect(0, 0, w, h, col = "grey75", border = NA)
  graphics::rect(0, 0, w / 2, h, col = "grey55", border = NA)
  # 三个几何体
  graphics::rect(w * 0.10, h * 0.55, w * 0.32, h * 0.85, col = "grey30", border = NA)
  th <- seq(0, 2 * pi, length.out = 200)
  graphics::polygon(w * 0.72 + w * 0.16 * cos(th), h * 0.30 + h * 0.20 * sin(th),
                    col = "grey38", border = NA)
  graphics::polygon(c(w * 0.15, w * 0.40, w * 0.27), c(h * 0.40, h * 0.40, h * 0.15),
                    col = "grey90", border = NA)
  # 一条斜线 + 一个字：高频细节，SVD 一截断就糊掉
  graphics::segments(w * 0.05, h * 0.92, w * 0.95, h * 0.55, lwd = 2, col = "black")
  graphics::text(w * 0.55, h * 0.75, "A", cex = 7, font = 2, col = "black")
  # 两块彩色：给「三个通道是三个不同的矩阵」留素材
  graphics::rect(w * 0.46, h * 0.05, w * 0.64, h * 0.24, col = "#E07B39", border = NA)
  graphics::polygon(w * 0.85 + w * 0.07 * cos(th), h * 0.14 + h * 0.09 * sin(th),
                    col = "#3C74A8", border = NA)
  graphics::par(op)
  grDevices::dev.off()
  on.exit(unlink(f))
  suppressPackageStartupMessages(requireNamespace("ROpenCV", quietly = TRUE))
  a <- il_read(f)
  a
}

# 读一张图，缩到长边不超过 max_side，转成 RGB 的 (h, w, 3) 数组。
# 默认的图是 il_chart()，可以用环境变量 LINALG_IMG 换成自己的照片：
#   LINALG_IMG=D:/photo.jpg Rscript img02-operations.R
il_read <- function(path = NULL, max_side = 360L) {
  if (is.null(path)) path <- Sys.getenv("LINALG_IMG", unset = "")
  if (!nzchar(path)) return(il_chart())
  m <- ROpenCV::cv_imread(path, ROpenCV::IMREAD_COLOR_BGR)
  if (any(dim(m) == 0L)) stop("读不到这张图：", path, call. = FALSE)
  d <- max(dim(m)[1:2])
  if (d > max_side) {
    s <- max_side / d
    m <- ROpenCV::cv_resize(m, NULL, s, s)
  }
  a <- ROpenCV::as_array(m, numeric = TRUE)
  a[, , 3:1, drop = FALSE]                    # BGR -> RGB，只翻这一次
}

# 灰度：教材里的「线性组合」，三个通道的加权和。权重是 ITU-R BT.601。
il_gray <- function(img) {
  if (is.matrix(img)) return(img)
  0.299 * img[, , 1] + 0.587 * img[, , 2] + 0.114 * img[, , 3]
}

# 由三张同尺寸的单通道矩阵拼成彩色图。
il_rgb <- function(r, g, b) {
  a <- array(0, c(dim(r), 3L))
  a[, , 1] <- r; a[, , 2] <- g; a[, , 3] <- b
  a
}

# R 的数值数组 -> 8 位 cvMat。
# 直接 Mat(双精度数组) 得到的是 CV_64F（depth = 6），cvtColor / warpAffine
# 这类函数会拒收；转成 raw 再 Mat() 才是 CV_8U（depth = 0）。
il_cv <- function(x) {
  a <- round(pmin(pmax(x, 0), 255))
  r <- as.raw(a)
  dim(r) <- dim(x)
  Mat(r)
}

## ---------------------------------------------------------------- 画图

# 把一个矩阵 / 数组 / cvMat 画到当前设备上。
il_show <- function(x, title = NULL, cex = 1.05) {
  if (inherits(x, "cvMat")) {
    graphics::plot(x)
  } else {
    a <- if (is.matrix(x)) pmin(pmax(x, 0), 255) / 255
         else pmin(pmax(x, 0), 255) / 255
    graphics::plot(grDevices::as.raster(a))
  }
  if (!is.null(title))
    graphics::title(main = title, cex.main = cex, font.main = 1)
  invisible(x)
}

# 拼图：把几张图排成网格存进当前设备（配合 il_fig()/il_done() 用）。
il_sheet <- function(imgs, titles = NULL, ncol = NULL) {
  n <- length(imgs)
  if (is.null(ncol)) ncol <- min(n, 3L)
  nrow <- ceiling(n / ncol)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::layout(matrix(seq_len(nrow * ncol), nrow, ncol, byrow = TRUE))
  for (i in seq_len(n)) {
    graphics::par(mar = c(0.5, 0.5, 1.8, 0.5))
    il_show(imgs[[i]], if (is.null(titles)) NULL else titles[i])
  }
  invisible(NULL)
}

# 拼图并直接存盘：按面板数算好画布尺寸，省得每处都手写宽高。
il_fig_sheet <- function(name, imgs, titles = NULL, ncol = NULL,
                         aspect = 0.78, panel = 320) {
  n <- length(imgs)
  if (is.null(ncol)) ncol <- min(n, 3L)
  nrow <- ceiling(n / ncol)
  p <- il_fig(name, width = panel * ncol,
              height = round(panel * aspect + 42) * nrow)
  il_sheet(imgs, titles, ncol)
  il_done(p)
}

## ---------------------------------------------------------------- 矩阵作用在图上

# 正映射（splat）：把源图每个像素按 x' = A x 搬到新位置。
# 奇异矩阵也能用 —— 这时所有点被压到一条线上，中间会留下空洞，
# 而这正是「行列式为 0」的样子。
il_warp_fwd <- function(img, A, out = dim(img)[1:2], offset = 0) {
  h <- dim(img)[1]; w <- dim(img)[2]
  g <- matrix(0, h, w)
  # 像素坐标先中心化，免得旋转时整幅图跑出画布
  xs <- seq_len(w) - (w + 1) / 2
  ys <- seq_len(h) - (h + 1) / 2
  # 点阵按 (h, w) 摆好：as.vector() 是列优先，行号（y）变化最快，
  # 正好和 as.vector(img) 的顺序对齐。
  P <- rbind(as.vector(outer(rep(1, h), xs)),
             as.vector(outer(ys, rep(1, w))))
  Q <- A %*% P
  qx <- round(Q[1, ] + (out[2] + 1) / 2 + offset)
  qy <- round(Q[2, ] + (out[1] + 1) / 2 + offset)
  ok <- qx >= 1 & qx <= out[2] & qy >= 1 & qy <= out[1]
  v <- if (is.matrix(img)) as.vector(img) else rowMeans(matrix(as.vector(img), ncol = 3))
  g[cbind(qy[ok], qx[ok])] <- v[ok]
  g
}

# 反映射（gather）：对输出图每个像素，问「你来自源图的哪里」，即用 A^{-1}。
# 只有可逆矩阵才做得到；采样用最近邻，所以放大时会有锯齿 —— 和
# cv_warpAffine 的双线性插值一对比，就能看出「插值」在补什么。
il_warp_inv <- function(img, A, out = dim(img)[1:2]) {
  Ai <- solve(A)
  h <- dim(img)[1]; w <- dim(img)[2]
  g <- if (is.matrix(img)) matrix(0, out[1], out[2])
       else array(0, c(out[1], out[2], 3L))
  xs <- seq_len(out[2]) - (out[2] + 1) / 2
  ys <- seq_len(out[1]) - (out[1] + 1) / 2
  P <- rbind(as.vector(outer(rep(1, out[1]), xs)),
             as.vector(outer(ys, rep(1, out[2]))))
  S <- Ai %*% P
  sx <- round(S[1, ] + (w + 1) / 2)
  sy <- round(S[2, ] + (h + 1) / 2)
  ok <- sx >= 1 & sx <= w & sy >= 1 & sy <= h
  idx <- rep(NA_integer_, length(sx))
  idx[ok] <- (sx[ok] - 1L) * h + sy[ok]          # R 是列优先：行号变化最快
  if (is.matrix(img)) {
    v <- as.vector(img)[idx]
    g[] <- ifelse(is.na(v), 0, v)
  } else {
    flat <- matrix(as.vector(img), ncol = 3L)[idx, , drop = FALSE]
    flat[is.na(flat)] <- 0
    g[, , 1] <- matrix(flat[, 1], out[1], out[2])
    g[, , 2] <- matrix(flat[, 2], out[1], out[2])
    g[, , 3] <- matrix(flat[, 3], out[1], out[2])
  }
  g
}


# 在图上画一条线段（直接改像素，省得纠结图形设备的坐标系）。
# 用作「把特征向量画在图上」这类叠加。
il_draw <- function(img, x0, y0, x1, y1, val = 0, width = 1) {
  n <- max(abs(c(x1 - x0, y1 - y0))) + 1
  for (d in seq_len(width) - 1) {
    xs <- round(seq(x0, x1, length.out = n)) + d
    ys <- round(seq(y0, y1, length.out = n))
    ok <- xs >= 1 & xs <= ncol(img) & ys >= 1 & ys <= nrow(img)
    if (is.matrix(img)) {
      img[cbind(ys[ok], xs[ok])] <- val
    } else {
      for (ch in 1:3) img[cbind(ys[ok], xs[ok], ch)] <- val[ch]
    }
  }
  img
}

## ---------------------------------------------------------------- 比较两张图

# 均方根误差：和教材里向量距离的定义完全一致，只是维数高得多。
il_rmse <- function(a, b) sqrt(mean((as.vector(a) - as.vector(b))^2))

# 两张图当作两个高维向量时的夹角余弦。
il_cos <- function(a, b) {
  x <- as.vector(a - mean(a)); y <- as.vector(b - mean(b))
  sum(x * y) / sqrt(sum(x^2) * sum(y^2))
}

il_psnr <- function(a, b) {
  e <- il_rmse(a, b)
  if (e == 0) return(Inf)
  20 * log10(255 / e)
}

# 截断到像素能表示的范围。逐像素运算之后必须做，否则会溢出/变负。
il_clamp <- function(m, lo = 0, hi = 255) pmin(pmax(m, lo), hi)

# 把任意实数矩阵映到 0..255，用于显示（不改变后续计算用的原数据）。
il_disp <- function(m) {
  r <- range(m, finite = TRUE)
  if (diff(r) == 0) return(matrix(0, nrow(m), ncol(m)))
  255 * (m - r[1]) / diff(r)
}
