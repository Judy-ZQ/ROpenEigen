# ===========================================================================
# 实验三 · 向量组的线性相关性及线性空间
#
# 对应教材：第三章
#   第一节 n 维向量              第五节 向量空间
#   第二节 向量组的线性相关与线性无关  第六节 线性空间与线性变换
#   第三节 向量组的最大线性无关组与秩  第七节 线性空间和线性变换的应用
#   第四节 正交向量组
# 建议学时：3 学时
#
# 实验目的
#   1. 把「线性相关」翻译成「齐次方程组有非零解」，用秩一步判定；
#   2. 会用行最简形找出最大无关组，并把其余向量用它表示出来；
#   3. 实现 Gram-Schmidt 正交化，并与 QR 分解对照；
#   4. 会用基与坐标、过渡矩阵做坐标变换；
#   5. 在线性变换层面验证 B = P^{-1}AP（同一变换在不同基下的矩阵）。
#
# 运行：Rscript lab03-vector-space.R   或在 RStudio 里 source()
# ===========================================================================

# 定位 helpers.R：脚本自己所在目录 -> 当前目录（及往上的包源码树）-> 已安装的包。
# 三种入口都能跑：Rscript labXX.R、在 RStudio 里 source()、包内回归测试。
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
  # (2) 当前目录，以及从当前目录往上 4 层找包源码树（<root>/inst/labs）
  d <- getwd()
  for (i in 0:4) {
    cand <- c(cand, d, file.path(d, "inst", "labs"), file.path(d, "labs"))
    d <- dirname(d)
  }
  # (3) 已安装的包（安装副本可能比源码树旧，所以放在最后）
  cand <- c(cand, system.file("labs", package = "ROpenEigen"))
  hit <- cand[file.exists(file.path(cand, "helpers.R"))]
  if (!length(hit))
    stop("找不到 helpers.R。三种跑法都支持：在 inst/labs 下 Rscript、",
         "在 RStudio 里 source()、或先安装 ROpenEigen。", call. = FALSE)
  hit[1]
})
source(file.path(.lab_dir, "helpers.R"))
suppressPackageStartupMessages(library(ROpenEigen))

# ===========================================================================
# 一、n 维向量与内积（教材 §3.1）
# ===========================================================================
say("== 3.1 向量与内积 ==")
a <- eigenVec(c(1, 2, 2)); b <- eigenVec(c(2, 1, 0))
say("a = (1,2,2)，b = (2,1,0)")
say("2a - 3b = ", paste(as.vector(2 * a - 3 * b), collapse = ", "))
say("内积 (a,b) = ", eigen_dot(a, b), "（= 1*2 + 2*1 + 2*0）")
say("长度 |a| = ", sqrt(sum(as.vector(a)^2)), "，|b| = ", sqrt(sum(as.vector(b)^2)))
say("eigen_norm(a, \"2\") = ", eigen_norm(a, "2"), "（与上一致）")

# 夹角余弦 cos θ = (a,b) / (|a||b|)
cos_t <- eigen_dot(a, b) / (eigen_norm(a, "2") * eigen_norm(b, "2"))
say("cos(a,b) = ", cos_t, "，夹角 = ", round(acos(cos_t) * 180 / pi, 4), " 度")

# 柯西-施瓦茨不等式 |(a,b)| <= |a||b|，随机检验
say("\n柯西-施瓦茨不等式随机检验（20 组）：")
set.seed(3)
ok <- TRUE
for (t in 1:20) {
  u <- rnorm(5); v <- rnorm(5)
  ok <- ok && abs(sum(u * v)) <= sqrt(sum(u^2)) * sqrt(sum(v^2)) + 1e-12
}
say("  全部满足：", ok)

# 正交：内积为 0
say("(1,0,0) 与 (0,1,0) 正交：", eigen_dot(eigenVec(c(1, 0, 0)), eigenVec(c(0, 1, 0))) == 0)

# ===========================================================================
# 二、线性相关与线性无关（教材 §3.2）
# ===========================================================================
say("\n== 3.2 线性相关性 ==")
# 判据：α1..αs 线性相关 ⇔ 存在不全为零的 k_i 使 Σ k_i α_i = 0
#       ⇔ 齐次方程组 (α1..αs) x = 0 有非零解 ⇔ R(α1..αs) < s
is_dependent <- function(...) {
  A <- do.call(cbind, lapply(list(...), as.vector))
  eigen_rank(eigenMat(A)) < ncol(A)
}

a1 <- c(1, 1, 1); a2 <- c(1, 2, 3); a3 <- c(3, 4, 5)
say("α1=(1,1,1) α2=(1,2,3) α3=(3,4,5)")
say("  R(α1,α2,α3) = ", eigen_rank(eigenMat(cbind(a1, a2, a3))),
    " < 3，所以线性相关")
say("  is_dependent() 判定：", is_dependent(a1, a2, a3))
# 找出相关的那组系数：解 (α1..α3)x = 0
say("  相关关系（零空间）：")
print(round(nullspace(cbind(a1, a2, a3)), 10))
say("  验证 k1α1+k2α2+k3α3 = ",
    paste(round(cbind(a1, a2, a3) %*% nullspace(cbind(a1, a2, a3)), 10), collapse = ", "))

# 换成线性无关的一组
b1 <- c(1, 0, 0); b2 <- c(1, 1, 0); b3 <- c(1, 1, 1)
say("β1=(1,0,0) β2=(1,1,0) β3=(1,1,1)：R = ",
    eigen_rank(eigenMat(cbind(b1, b2, b3))), "，线性无关")

# 三条常用结论
say("\n常用结论：")
say("  (1) 含零向量的向量组必相关：", is_dependent(c(1, 1), c(0, 0), c(2, 3)))
say("  (2) 两个向量相关 ⇔ 成比例：",
    is_dependent(c(1, 2, 3), c(2, 4, 6)),
    " / ", is_dependent(c(1, 2, 3), c(2, 4, 7)))
say("  (3) 部分相关则整体相关：若 α1,α2 相关，加入 α3 后仍相关：",
    is_dependent(c(1, 2), c(2, 4), c(7, 1)))

# 向量个数 > 维数时必相关（教材 §3.2 定理）
say("  4 个 3 维向量必相关：", is_dependent(c(1, 0, 0), c(0, 1, 0),
                                              c(0, 0, 1), c(1, 1, 1)))

# ===========================================================================
# 三、最大无关组与秩（教材 §3.3）
# ===========================================================================
say("\n== 3.3 最大无关组 ==")
# 取 4 个 3 维向量，秩为 2，最大无关组含 2 个向量
V <- cbind(c(1, 2, 3), c(2, 4, 6), c(1, 1, 1), c(3, 5, 7))
colnames(V) <- paste0("α", 1:4)
say("向量组（按列）:")
print(V)
R <- rref(V)
say("行最简形：")
print(R$rref)
say("主元列 = ", paste(R$pivots, collapse = ", "),
    " → 最大无关组是 {", paste0("α", R$pivots, collapse = ", "), "}")
say("向量组的秩 = ", R$rank, "，与 eigen_rank 一致：",
    R$rank == eigen_rank(eigenMat(V)))

# 行最简形的主元列还直接给出「其余向量用最大无关组表示」的系数：
# 第 j 列（非主元列）= Σ_i R[i, j] * 第 i 个主元列
say("\n其余向量用最大无关组表示：")
for (j in setdiff(seq_len(ncol(V)), R$pivots)) {
  coef <- R$rref[seq_along(R$pivots), j]
  say(sprintf("  α%d = %s", j,
              paste(sprintf("%.0f·α%d", coef, R$pivots), collapse = " + ")))
  expect_near(V[, j], V[, R$pivots] %*% coef, 1e-12, "线性表示系数不对")
}

# 增删向量对秩的影响：加一个能被表示的向量，秩不变
V2 <- cbind(V, V[, 1] + V[, 3])
say("加一个 α1+α3 之后，秩仍为 ", eigen_rank(eigenMat(V2)))

# ===========================================================================
# 四、正交向量组与 Gram-Schmidt（教材 §3.4）
# ===========================================================================
say("\n== 3.4 正交化 ==")
G <- cbind(c(1, 1, 1), c(1, 0, 1), c(0, 1, 1))
say("原始向量组（按列）：")
print(G)
gs <- gram_schmidt(G)
say("Schmidt 正交化并单位化后：")
print(round(gs$Q, 6))
say("两两正交（Q'Q 应为 E）：")
print(round(t(gs$Q) %*% gs$Q, 12))
expect_near(t(gs$Q) %*% gs$Q, diag(3), 1e-12, "正交化结果不是规范正交组")
expect_near(gs$Q %*% gs$R, G, 1e-12, "A = QR 不成立")

# 与库函数 QR 分解对照：Q 的列相同，但可能差一个符号（R 的对角元可正可负）
qr_res <- eigen_qr(eigenMat(G))
Qlib <- as.matrix(qr_res$matrixQ())
say("库 QR 的 Q：")
print(round(Qlib, 6))
say("两套 Q 的列向量绝对值相同：",
    relerr(abs(gs$Q), abs(Qlib)) < 1e-12,
    "（差的正负号无关紧要）")

# 正交矩阵的性质：Q'Q = QQ' = E，|Q| = ±1，且保持长度
say("Q'Q = E：", relerr(t(gs$Q) %*% gs$Q, diag(3)) < 1e-12)
say("|det(Q)| = ", abs(det(eigenMat(gs$Q))), "（应为 1）")
x <- c(3, -1, 2)
say("正交变换保持长度：|x| = ", sqrt(sum(x^2)),
    "，|Q'x| = ", sqrt(sum((t(gs$Q) %*% x)^2)))
y <- c(1, 2, -3)
say("正交变换保持内积：<x,y> = ", sum(x * y),
    "，<Q'x, Q'y> = ", sum((t(gs$Q) %*% x) * (t(gs$Q) %*% y)))

# ===========================================================================
# 五、向量空间：基、维数、坐标、过渡矩阵（教材 §3.5）
# ===========================================================================
say("\n== 3.5 向量空间 ==")
# (1) 齐次方程组的解空间就是 R^n 的子空间，维数 = n - R(A)
A <- matrix(c(1, 2, 3, 4, 2, 4, 6, 8, 1, 1, 1, 1), 3, 4, byrow = TRUE)
N <- nullspace(A)
say("A 的秩 = ", eigen_rank(eigenMat(A)), "，未知数 n = 4")
say("解空间维数 = n - R(A) = ", 4 - eigen_rank(eigenMat(A)),
    "，实际算出基础解系含 ", ncol(N), " 个向量")
say("基础解系：")
print(round(N, 6))
say("A %*% N 应当是零矩阵：", max(abs(A %*% N)) < 1e-12)

# (2) 基与坐标：同一向量在两组织基下的坐标不同，靠过渡矩阵联系
#     基 I：e1=(1,0,0), e2=(0,1,0), e3=(0,0,1)
#     基 II：f1=(1,1,0), f2=(0,1,1), f3=(1,0,1)
F <- cbind(c(1, 1, 0), c(0, 1, 1), c(1, 0, 1))
x <- c(2, 3, 1)
coord_F <- as.vector(solve(eigenMat(F), eigenVec(x)))     # x = F · coord
say("x = (2,3,1)")
say("在标准基下的坐标：", paste(x, collapse = ", "))
say("在基 {f1,f2,f3} 下的坐标：", paste(round(coord_F, 10), collapse = ", "))
expect_near(F %*% coord_F, x, 1e-12, "坐标与基不匹配")

# (3) 过渡矩阵：从基 I 到基 II 的过渡矩阵 C 满足 (f1,f2,f3) = (e1,e2,e3) C
#     这里 C 就是 F 本身（因为标准基就是单位矩阵）。坐标变换公式 x_I = C x_II。
say("过渡矩阵 C =")
print(F)
say("坐标变换 x_I = C x_II 成立：", relerr(F %*% coord_F, x) < 1e-12)

# ===========================================================================
# 六、线性变换与矩阵表示（教材 §3.6）
# ===========================================================================
say("\n== 3.6 线性变换的矩阵 ==")
# 教材 §3.7 的例子：多项式空间 P3 的基 α1=1, α2=x, α3=x^2, α4=x^3，
# 求导算子 D 在这组基下的矩阵。
#   D(α1)=0, D(α2)=α1, D(α3)=2α2, D(α4)=3α3
Dmat <- matrix(c(0, 1, 0, 0,
                 0, 0, 2, 0,
                 0, 0, 0, 3,
                 0, 0, 0, 0), 4, 4, byrow = TRUE)
say("D 在基 (1, x, x^2, x^3) 下的矩阵：")
print(Dmat)
# 验证：把系数向量 (c0,c1,c2,c3) 乘 Dmat，应当得到导数的系数 (c1, 2c2, 3c3, 0)
cfs <- c(2, -1, 4, 5)                       # 2 - x + 4x^2 + 5x^3
say("f(x) = 2 - x + 4x^2 + 5x^3 的系数：", paste(cfs, collapse = ", "))
say("矩阵法求导：", paste(as.vector(Dmat %*% cfs), collapse = ", "))
say("手算求导：系数应为  -1, 8, 15, 0")
expect_near(Dmat %*% cfs, c(-1, 8, 15, 0), 1e-12, "求导矩阵不对")

# 换基：改用 (α2, α3, α1, α4) 这组基，求 D 在新基下的矩阵。
# 新基用旧基表示：f1=α2, f2=α3, f3=α1, f4=α4，过渡矩阵 P 的列就是它们。
# 注意 P 的每一「列」是新基向量在旧基下的坐标：
#   f1 = α2 → (0,1,0,0)'，f2 = α3 → (0,0,1,0)'，f3 = α1 → (1,0,0,0)'，f4 = α4
P <- cbind(c(0, 1, 0, 0), c(0, 0, 1, 0), c(1, 0, 0, 0), c(0, 0, 0, 1))
Dnew <- solve(eigenMat(P)) %*% eigenMat(Dmat) %*% eigenMat(P)
say("换基后的矩阵 B = P^{-1} A P：")
print(round(as.matrix(Dnew), 10))
# 教材给出的结果：
Dtext <- matrix(c(0, 2, 0, 0,
                  0, 0, 0, 3,
                  1, 0, 0, 0,
                  0, 0, 0, 0), 4, 4, byrow = TRUE)
say("与教材给出的结果一致：", relerr(as.matrix(Dnew), Dtext) < 1e-12)

# ===========================================================================
# 七、应用：多项式空间的内积与正交化（教材 §3.7）
# ===========================================================================
say("\n== 3.7 多项式空间的内积 ==")
# 在 P2 = span{1, x, x^2} 上定义内积 (f,g) = ∫_0^1 f(x)g(x) dx。
# 基 1, x, x^2 的 Gram 矩阵 G_ij = ∫_0^1 x^(i-1) x^(j-1) dx = 1/(i+j-1)
Gram <- outer(1:3, 1:3, function(i, j) 1 / (i + j - 1))
say("基 {1, x, x^2} 的 Gram 矩阵（正好是 3 阶 Hilbert 矩阵）：")
print(round(Gram, 6))
# 这组基并不正交：(1, x) = 1/2 ≠ 0
say("(1, x) = ", Gram[1, 2], " ≠ 0，所以这组基不正交")

# 这里的「向量」是多项式在基 {1, x, x^2} 下的坐标，内积由 Gram 矩阵给出，
# 所以不能直接用上一个 gram_schmidt()（那个默认内积是坐标点积）。
Qp <- gram_schmidt_metric(diag(3), Gram)
say("正交化后各正交多项式在基 {1, x, x^2} 下的坐标（按列）：")
print(round(Qp, 6))
say("验证 Q'GQ = E：")
print(round(t(Qp) %*% Gram %*% Qp, 12))
# 换回多项式：q_k(x) = Σ_j Q[j,k] x^(j-1)
# 把坐标写成多项式
poly_str <- function(cf) {
  terms <- sprintf("%+.4f", cf)
  nm <- c("", "·x", "·x^2")
  paste0(paste0(terms, nm), collapse = " ")
}
say("于是正交多项式为：")
for (k in 1:3) say(sprintf("  q%d(x) = %s", k, poly_str(Qp[, k])))
say("验证两两正交：")
print(round(t(Qp) %*% Gram %*% Qp, 12))

# 把每个 q_k 除以它自己的最高次项系数，得到首一（monic）形式，
# 会发现它们正是移位 Legendre 多项式的首一形式：
#   1,  x - 1/2,  x^2 - x + 1/6
# 「首一」= 除以最高次非零项系数（q1、q2 的最高次不是二次）
monic <- function(cf) {
  k <- max(which(abs(cf) > 1e-12))
  cf / cf[k]
}
say("化为首一形式后：")
for (k in 1:3) say(sprintf("  q%d(x) ~ %s", k, poly_str(monic(Qp[, k]))))
shifted_monic <- list(c(1, 0, 0), c(-0.5, 1, 0), c(1 / 6, -1, 1))
say("与 1, x-1/2, x^2-x+1/6（乘以 6 就是 6x^2-6x+1）一致：",
    all(vapply(1:3, function(k) {
      relerr(matrix(monic(Qp[, k]), 3, 1), matrix(shifted_monic[[k]], 3, 1)) < 1e-10
    }, logical(1))))

# ===========================================================================
# 八、实验习题
# ===========================================================================
# 1. 判断 α1=(1,1,2,3)、α2=(1,-1,1,1)、α3=(1,3,3,5)、α4=(1,-5,3,-1) 的相关性，
#    求一个最大无关组，并把其余向量用它表示。（提示：rref(cbind(...))）
#
# 2. 设 α1=(1,1,0)、α2=(0,1,1)、α3=(1,0,1)。验证它们线性无关，并把
#    β=(2,3,1) 用它们表示。再换一组基算同一个 β 的坐标，用过渡矩阵核对。
#
# 3. 教材 §3.4 的 Schmidt 正交化公式里，每一步都要减去前面的投影。
#    如果用修改过的 Gram-Schmidt（MGS，每步重新计算投影）会怎样？
#    构造一个「列几乎平行」的矩阵，比较两种做法的 Q'Q 偏离单位阵多少。
#
# 4. 验证：n 维空间中任意 n+1 个向量必线性相关。用随机矩阵试 n = 3..8。
#
# 5. 教材 §3.7 用多项式空间讲了线性变换。请把「求导 + 乘 2」这个变换
#    (D+2I) 在基 (1,x,x^2,x^3) 下的矩阵写出来，并验证 (D+2I)f 的结果。
#
# 思考题
#   (1) 「向量组线性相关」与「齐次方程组有非零解」是同一件事的两种说法。
#       那么「向量组能表示 β」对应方程组的什么性质？
#   (2) 为什么最大无关组不唯一，但秩唯一？（提示：从行最简形看主元列。）
#   (3) 正交化把一组基变成规范正交基，代价是什么？（提示：条件数与舍入误差放大。）

say("\nlab03 结束。")
