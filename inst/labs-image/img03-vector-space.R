# ===========================================================================
# 实验三 · 图像是向量：内积、正交、秩、子空间
#
# 对应教材：第三章 §3.1–§3.4（内积、相关无关、最大无关组与秩、正交化）
# 建议学时：3
#
# 实验目的
#   1. 接受「一张图就是一个 76800 维向量」这件事，并把内积、范数、夹角算出来；
#   2. 看到线性组合就是「叠图」，正交就是「两张图互不相似」；
#   3. 看到投影：把一张图投影到两张图张成的平面上，剩下的是什么；
#   4. 把 Gram-Schmidt 打在图上：把「相关的部分」剥掉，只留正交残差。
#
# 运行：Rscript img03-vector-space.R   或在 RStudio 里 source()
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

# 造两张「有关系」的图：B 是 A 平移过的版本，C 是 A 模糊过的版本。
B <- A[, c(31:w, rep(w, 30))]                     # 整体左移 30 像素
C <- as_array(cv_GaussianBlur(il_cv(A), c(15, 15), 5), numeric = TRUE)

say("A、B、C 都是 ", h, " × ", w, " 的灰度矩阵，也就是 ", h * w, " 维向量。")
say("B = A 左移 30 像素，C = A 高斯模糊（σ = 5）。")
say("图会存到：", il_outdir())

# 把图当向量用的两个基本函数
dotv <- function(x, y) sum(as.vector(x) * as.vector(y))
normv <- function(x) sqrt(dotv(x, x))

# ===========================================================================
# 一、内积、范数、夹角
# ===========================================================================
say("\n== 3.1 内积与夹角 ==")

say("||A|| = ", round(normv(A), 2), "，||B|| = ", round(normv(B), 2))
say("<A, B> = ", format(dotv(A, B), big.mark = ","), "，<A, C> = ", format(dotv(A, C), big.mark = ","))
say("两个 76800 维向量照样有夹角：")
ca <- function(x, y) dotv(x, y) / (normv(x) * normv(y))
say("  cos(A, B) = ", round(ca(A, B), 6), " -> ", round(acos(ca(A, B)) * 180 / pi, 2), " 度")
say("  cos(A, C) = ", round(ca(A, C), 6), " -> ", round(acos(ca(A, C)) * 180 / pi, 2), " 度")
say("余弦越大越相似：模糊只改变一点（6.9 度），平移 30 像素后角度到 18.7 度。")
say("真正「不像」的两张图，余弦会接近 0 —— 下一节构造给你看。")

# ===========================================================================
# 二、线性组合 = 叠图
# ===========================================================================
say("\n== 3.2 线性组合就是「淡入淡出」 ==")

mix <- function(t) t * A + (1 - t) * B
say("0.5·A + 0.5·B 的某个像素：0.5×", round(A[60, 100], 1), " + 0.5×", round(B[60, 100], 1),
    " = ", round(mix(0.5)[60, 100], 1), "（逐像素加权，没有别的魔法）")
say("α 从 0 扫到 1，就是 B 渐渐变成 A。这一族图正好是 span{A, B} 里的一条线段。")

il_fig_sheet("img03-combo",
             list(B, mix(0.25), mix(0.5), mix(0.75), A),
             c("B (shifted)", "0.25A + 0.75B", "0.5A + 0.5B", "0.75A + 0.25B", "A"), ncol = 5)

# ===========================================================================
# 三、正交 = 不相似
# ===========================================================================
say("\n== 3.3 正交是什么意思 ==")

# (A - mean(A)) 与常向量 1 正交：这就是「把均值减掉」的几何含义
Am <- A - mean(A)
one <- matrix(1, h, w)
say("<A - mean(A), 1> = ", format(dotv(Am, one), scientific = TRUE),
    "（浮点误差级别，理论上是 0）")
say("所以「去均值」= 把向量投影到「与常向量正交」的子空间里。")

# 「反相」的图（255 - A）是不是与 A 正交？不是 —— 因为它俩共享同一个大均值。
neg <- 255 - A
say("cos(A, 255 - A) = ", round(ca(A, neg), 6), "：看着完全相反，余弦却是正的。")
say("原因：两张图都有一个很大的常数底（均值 ", round(mean(A), 1),
    " 与 ", round(mean(neg), 1), "），")
say("这个共同的分量把余弦拉高了。去掉均值再看：cos = ", round(ca(Am, neg - mean(neg)), 6),
    " —— 这才是「完全负相关」。")
say("所以「正交」讲的是内积为 0，和「看起来像不像」不是一回事。")

# ===========================================================================
# 四、Gram-Schmidt：把相关的部分剥掉
# ===========================================================================
say("\n== 3.4 Gram-Schmidt：B 里有多少是 A，剩下的是什么 ==")

c1 <- dotv(B, A) / dotv(A, A)
R <- B - c1 * A                                   # 残差
say("B 在 A 方向上的坐标 c = <B,A>/<A,A> = ", round(c1, 6))
say("残差 R = B - c·A，检验 <R, A> = ", format(dotv(R, A), scientific = TRUE), "（应为 0）")
say("||R|| / ||B|| = ", round(normv(R) / normv(B), 4),
    "：B 的长度有 ", round(100 * normv(R) / normv(B), 1), "% 是 A 这个方向给不了的。")

il_fig_sheet("img03-gs",
             list(A, B, c1 * A, il_disp(A + R)),
             c("A", "B", "the part of B along A", "residual R (B minus its A-part)"), ncol = 4)

# ===========================================================================
# 五、投影与子空间
# ===========================================================================
say("\n== 3.5 投影：C 在 span{A, B} 上能还原多少 ==")

# 先正交化 A、B，再用正交集做投影（和教材里「把基正交化再算坐标」一个道理）
q1 <- A / normv(A)
v2 <- B - dotv(B, q1) * q1
q2 <- v2 / normv(v2)
say("q1、q2 是 span{A, B} 的一组标准正交基：<q1,q2> = ",
    format(dotv(q1, q2), scientific = TRUE), "，||q1|| = ", round(normv(q1), 10),
    "，||q2|| = ", round(normv(q2), 10))

proj <- dotv(C, q1) * q1 + dotv(C, q2) * q2
res <- C - proj
say("C 的投影长度占 C 的 ", round(100 * normv(proj) / normv(C), 2), "%（能量按平方算就是 ",
    round(100 * dotv(proj, proj) / dotv(C, C), 2), "%）")
say("残差与 q1、q2 都正交：<res,q1> = ", format(dotv(res, q1), scientific = TRUE),
    "，<res,q2> = ", format(dotv(res, q2), scientific = TRUE))
say("投影误差 ||res|| / ||C|| = ", round(normv(res) / normv(C), 4),
    "：两张图张成的平面解释了模糊图的绝大部分，剩下这一点点就是「平面外的分量」。")

il_fig_sheet("img03-proj",
             list(C, il_disp(proj), il_disp(abs(res) * 3)),
             c("C (blurred)", "projection onto span{A,B}", "residual (x3 for visibility)"), ncol = 3)

# ===========================================================================
# 六、秩：到底需要几张图才够
# ===========================================================================
say("\n== 3.6 秩 = 最大无关组的个数 ==")

# 三张图的「行向量」组成的矩阵：这里的秩说明这三张图互相独立到什么程度
M3 <- rbind(as.vector(A), as.vector(B), as.vector(C))
say("把 A、B、C 摞成 3 × 76800 的矩阵，秩 = ", eigen_rank(eigenMat(M3)),
    "（= 3，说明三张图线性无关）")

# 但 C 是 A 模糊来的，它几乎落在 A 的附近：奇异值能看出「几乎相关」
sv <- as.vector(eigen_svd(eigenMat(M3))$singular_values())
say("它的三个奇异值：", paste(format(sv, digits = 6), collapse = ", "))
say("第三个比第一个小 ", format(sv[1] / sv[3], digits = 4),
    " 倍 —— 「线性无关」和「几乎线性相关」在数值上是两件事。")

# ===========================================================================
# 七、实验习题
# ===========================================================================
say("\n== 3.7 习题 ==")
say("1. 用 il_warp_inv() 造一张把 A 旋转 5 度的图，算 cos(A, A_rot)。")
say("   旋转角度变大时这个余弦怎么变？画出来。")
say("2. 找出一个 α，使 ||A - α·B|| 最小。它和 <A,B>/<B,B> 是什么关系？")
say("3. 把 3.6 节的第三张图换成 2A - B，秩变成多少？为什么？")
say("4. 3.3 节说「去均值 = 投影到与 1 正交的子空间」。那「除以 255 归一化」")
say("   在几何上是什么操作？（提示：数乘一个向量，不改变方向。）")

say("\nimg03 结束：图是向量，内积是相似度，正交是无关，投影是「最好的近似」。")
