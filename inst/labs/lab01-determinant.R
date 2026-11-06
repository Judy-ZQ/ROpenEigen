# ===========================================================================
# 实验一 · 行列式
#
# 对应教材：第一章 行列式
#   第一节 预备知识（排列、逆序数）
#   第二节 n 阶行列式的定义
#   第三节 行列式的性质
#   第四节 行列式的展开（余子式、代数余子式、伴随矩阵）
#   第五节 行列式的应用（向量积、几何判定）
# 建议学时：2 学时
#
# 实验目的
#   1. 把「逆序数定义」变成能跑的代码，看清行列式为什么是 n! 项之和；
#   2. 逐条验证教材第三节的七条性质，理解每条性质在计算上的用途；
#   3. 用代数余子式实现伴随矩阵，验证 A·adj(A) = |A|·E；
#   4. 用行列式解决几何判定问题（共线、共面、三平面交于一点、面积体积）；
#   5. 认识数值行列式的两个陷阱：阶乘级代价、溢出与下溢。
#
# 运行：Rscript lab01-determinant.R   或在 RStudio 里 source()
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
# 一、预备知识：排列与逆序数（教材 §1.1）
# ===========================================================================
say("== 1.1 逆序数 ==")
# 逆序数 τ(p)：对每个 p_i，数一数排在它前面又比它大的数有几个，然后相加。
say("τ(3421)   = ", inversions(c(3, 4, 2, 1)))
say("τ(517924) = ", inversions(c(5, 1, 7, 9, 2, 4)))
# 逆序数为偶数叫偶排列，奇数叫奇排列。行列式定义里 (-1)^τ 就是从这里来的。
say("奇偶性：τ(1234) = ", inversions(1:4), " → 偶排列")

# 3 阶行列式的 6 项，把每一项的排列、逆序数、符号、乘积因子都列出来
say("\n3 阶行列式按定义展开的 6 项：")
P <- permutations(3)
for (k in seq_len(nrow(P))) {
  p <- P[k, ]
  say(sprintf("  a1%d a2%d a3%d   τ=%d  sign=%+d", p[1], p[2], p[3],
              inversions(p), ifelse(inversions(p) %% 2 == 0, 1, -1)))
}

# ===========================================================================
# 二、定义：n! 项之和（教材 §1.2）
# ===========================================================================
say("\n== 1.2 按定义计算 ==")
A3 <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 10), 3, 3, byrow = TRUE)
say("A = ")
print(A3)
say("定义法   ：", det_by_definition(A3))
say("eigen_det：", det(eigenMat(A3)))
expect_near(det_by_definition(A3), det(eigenMat(A3)), 1e-12, "定义法与 eigen_det 不符")

# 定义法的代价是 n!。下面这张表说明为什么它只能用于教学：
say("\nn 与项数 n! 的关系：")
for (n in 2:10) say(sprintf("  n=%2d  %10.0f 项", n, factorial(n)))
say("n=15 已经是 1.3e12 项，所以真实计算必须用消元法（O(n^3)）。")

# 上三角行列式 = 主对角线之积（教材 §1.2 例）
U <- matrix(c(2, 3, -1, 0, 5, 4, 0, 0, 7), 3, 3, byrow = TRUE)
say("\n上三角矩阵：")
print(U)
say("det(U) = ", det(eigenMat(U)), "，对角线之积 = ", prod(diag(U)))

# ===========================================================================
# 三、性质逐条验证（教材 §1.3）
# ===========================================================================
say("\n== 1.3 行列式的性质 ==")
set.seed(11)
n <- 4
M <- matrix(sample(-4:4, n * n, replace = TRUE), n, n, byrow = TRUE)
d <- det(eigenMat(M))
say("取一个 4 阶整数矩阵 M，det(M) = ", d)
say("性质1  |M'| = |M| ：", relerr(det(eigenMat(t(M))), d) < 1e-12)

M2 <- M; M2[c(2, 3), ] <- M2[c(3, 2), ]        # 交换两行
say("性质2  交换两行变号：", relerr(det(eigenMat(M2)), -d) < 1e-12)

M3 <- M; M3[2, ] <- 5 * M3[2, ]                # 某行乘 5
say("性质3  一行乘 k，行列式乘 k：", relerr(det(eigenMat(M3)), 5 * d) < 1e-12)

say("性质3' |kM| = k^n|M|（n=", n, "）：",
    relerr(det(eigenMat(3 * M)), 3^n * d) < 1e-12)

M4 <- M; M4[2, ] <- M4[2, ] + 7 * M4[1, ]      # 一行加上另一行的倍数
say("性质4  行倍加不变：", relerr(det(eigenMat(M4)), d) < 1e-12)

# 有两行相同（或成比例）→ 行列式为 0
M5 <- M; M5[3, ] <- M5[1, ]
say("性质5  两行相同则 |M| = 0：", abs(det(eigenMat(M5))) < 1e-12)

N <- matrix(sample(-3:3, n * n, replace = TRUE), n, n, byrow = TRUE)
say("性质6  |MN| = |M||N|：",
    relerr(det(eigenMat(M %*% N)), det(eigenMat(M)) * det(eigenMat(N))) < 1e-10)

# 分块三角：| [[A, C], [0, D]] | = |A||D|（教材 §1.3 分块矩阵一节）
P1 <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
P2 <- matrix(c(5, 6, 7, 8), 2, 2, byrow = TRUE)
C <- matrix(c(9, 1, 2, 3), 2, 2, byrow = TRUE)
blk <- rbind(cbind(P1, C), cbind(matrix(0, 2, 2), P2))
say("性质7  分块三角 |[[A,C],[0,D]]| = |A||D|：",
    relerr(det(eigenMat(blk)), det(eigenMat(P1)) * det(eigenMat(P2))) < 1e-12)

# 可逆时 |A^{-1}| = 1/|A|
say("性质8  |M^{-1}| = 1/|M|：",
    relerr(det(eigenMat(solve(M))), 1 / d) < 1e-10)

# ===========================================================================
# 四、展开：余子式、代数余子式、伴随矩阵（教材 §1.4）
# ===========================================================================
say("\n== 1.4 按行展开与伴随矩阵 ==")
B <- matrix(c(1, 2, 3, 2, 5, 3, 1, 0, 8), 3, 3, byrow = TRUE)
say("B = ")
print(B)
say("按第 1 行展开：Σ a1j·A1j = ",
    sum(B[1, ] * cofactor_matrix(B)[1, ]))
say("按第 2 列展开：Σ ai2·Ai2 = ",
    sum(B[, 2] * cofactor_matrix(B)[, 2]))
say("det(B) = ", det(eigenMat(B)))

# 关键：用别的行的代数余子式去乘本行元素，结果必为 0（教材 §1.4 推论）。
# 这一条正是 A·adj(A) = |A|·E 的来源：对角线上是 |A|，其余位置全是 0。
say("错位展开 Σ a1j·A2j = ", sum(B[1, ] * cofactor_matrix(B)[2, ]), "（应为 0）")

adjB <- adjugate(B)
say("adj(B) = ")
print(adjB)
say("B %*% adj(B) 应当等于 |B|·E：")
print(round(as.matrix(eigenMat(B) %*% eigenMat(adjB)), 12))
expect_near(B %*% adjB, det(eigenMat(B)) * diag(3), 1e-12, "伴随矩阵恒等式不成立")

# 于是 B^{-1} = adj(B)/|B|，这就是教材里的伴随矩阵求逆公式
say("adj(B)/|B| 与 solve(B) 的相对误差：",
    format(relerr(adjB / det(eigenMat(B)), as.matrix(solve(eigenMat(B)))), scientific = TRUE))

# 克拉默法则放在实验二（教材把克拉默法则排在第二章第七节），这里先用伴随矩阵
# 把方程组解出来看一眼，正式讨论见 lab02。

# ===========================================================================
# 五、应用：向量积与几何判定（教材 §1.5）
# ===========================================================================
say("\n== 1.5 行列式的应用 ==")
# (1) 叉积的坐标式就是三个二阶行列式
a <- c(1, 2, 3); b <- c(-1, 0, 2)
cross_by_det <- c(det(matrix(c(a[2], a[3], b[2], b[3]), 2)),
                  -det(matrix(c(a[1], a[3], b[1], b[3]), 2)),
                  det(matrix(c(a[1], a[2], b[1], b[2]), 2)))
say("a × b（按行列式算） = ", paste(cross_by_det, collapse = ", "))
say("a × b（eigen_cross） = ", paste(round(as.vector(eigen_cross(eigenVec(a), eigenVec(b))), 12),
                                     collapse = ", "))
expect_near(cross_by_det, as.vector(eigen_cross(eigenVec(a), eigenVec(b))), 1e-12,
            "叉积的行列式表示不符")

# (2) 混合积 (a×b)·c = |a b c|，几何意义是平行六面体体积
cc <- c(2, 1, 1)
vol_det <- det(eigenMat(rbind(a, b, cc)))
vol_dot <- sum(cross_by_det * cc)
say("混合积（行列式） = ", vol_det, "，（点积） = ", vol_dot)
say("平行六面体体积 = |混合积| = ", abs(vol_det))
expect_near(vol_det, vol_dot, 1e-12, "混合积两种算法不一致")

# (3) 共面判定：三点共面（含原点）⇔ 行列式为 0
say("三向量共面？", abs(vol_det) < 1e-12, "（本例不共面）")
say("把 c 换成 a+b 之后共面？",
    abs(det(eigenMat(rbind(a, b, a + b)))) < 1e-12)

# (4) 三点共线：| x1 y1 1; x2 y2 1; x3 y3 1 | = 0，这就是教材 §1.5 例 3
#     里「过两点的直线方程」的行列式形式。
line_det <- function(p1, p2, p3) {
  det(eigenMat(rbind(c(p1, 1), c(p2, 1), c(p3, 1))))
}
P1 <- c(1, 2); P2 <- c(4, 6)
# 在直线上取点时，理论上行列式恰好是 0；浮点算出 1e-15 量级，
# 这就是实验零里说的「等于要用相对误差判断」。
say("过 (1,2)、(4,6) 的直线上取一点 (7,10)，行列式 = ",
    format(line_det(P1, P2, c(7, 10)), scientific = TRUE), "（数值零）")
say("不在直线上的 (7,11)，行列式 = ", line_det(P1, P2, c(7, 11)))

# (5) 三角形面积 = 1/2 |det[[x1,y1,1],[x2,y2,1],[x3,y3,1]]|
tri <- rbind(c(0, 0, 1), c(4, 0, 1), c(0, 3, 1))
say("三角形面积（行列式） = ", abs(det(eigenMat(tri))) / 2, "（应为 6）")

# (6) 两向量张成的平行四边形面积 = |det[[a1,a2],[b1,b2]]|
par_area <- abs(det(eigenMat(rbind(c(3, 1), c(1, 2)))))
say("平行四边形面积 = ", par_area)

# (7) 三平面交于一点的条件（教材 §1.5 例 2）：系数行列式 ≠ 0
#     x + y + z = 6, 2x - y + z = 3, x + 2y - z = 2
Apl <- matrix(c(1, 1, 1, 2, -1, 1, 1, 2, -1), 3, 3, byrow = TRUE)
bpl <- c(6, 3, 2)
say("三平面系数行列式 = ", det(eigenMat(Apl)), "（≠0 → 交于唯一一点）")
say("交点 = ", paste(round(as.vector(solve(eigenMat(Apl), eigenVec(bpl))), 10), collapse = ", "))

# 把第三个平面换成与前两个共线（不相交于一点）的情形
Apl2 <- Apl; Apl2[3, ] <- Apl2[1, ] + Apl2[2, ]
say("退化后系数行列式 = ", det(eigenMat(Apl2)), "（=0 → 不交于唯一一点）")

# ===========================================================================
# 六、数值提醒：溢出、下溢与该不该用行列式
# ===========================================================================
say("\n== 1.6 数值陷阱 ==")
# 对角元取 1e20，20 阶行列式就是 1e400，超过了 double 的上限 1.8e308
D <- diag(1e20, 20)
say("|1e20·E20| 直接算 = ", det(eigenMat(D)), "（溢出成 Inf）")
say("log|1e20·E20| = ", eigen_logdet(eigenMat(D)), "（对数版仍是有限数）")
say("20·log(1e20) = ", 20 * log(1e20))

# 下溢：Hilbert 矩阵阶数一高，det 掉到浮点下限以下就变 0
for (k in c(4, 8, 12, 16)) {
  Hk <- outer(1:k, 1:k, function(i, j) 1 / (i + j - 1))
  say(sprintf("  %2d 阶 Hilbert：det = %-12.3e logdet = %8.3f rcond = %.2e",
              k, det(eigenMat(Hk)), eigen_logdet(eigenMat(Hk)),
              eigen_rcond(eigenMat(Hk))))
}
say("结论：判断「行列式是否为 0」只在手算题里可靠；程序里要看秩与 rcond。")

# ===========================================================================
# 七、实验习题
# ===========================================================================
# 1.（教材习题一 第 1 题风格）求 τ(3421)、τ(517924)，并说明
#    「n 阶行列式中带负号的项有多少个」与 n 的关系。
#
# 2. 教材习题一第 3 题：5 阶行列式中含因子 a12·a23·a45 的项有几项？分别写出。
#    提示：固定了 3 个位置后，剩下两行取剩下两列，共 2! 种，逐一看符号。
#    验证代码：
#      part <- c(1, 2, NA, 4, NA)          # 第 i 行取第 part[i] 列
#      rest <- setdiff(1:5, c(2, 3, 5))
#      for (q in list(rest, rev(rest))) {
#        p <- part; p[3] <- q[1]; p[5] <- q[2]
#        cat(p, "τ =", inversions(p), "\n")
#      }
#
# 3. 教材习题一第 5 题：设
#        f(x) = | 1 1  1   1  |
#               | 1 2 -2   x  |
#               | 1 4  4  x^2 |
#               | 1 8 -8  x^3 |
#    求 f(x) = 0 的根。
#    这是节点为 1, 2, -2, x 的 Vandermonde 行列式，f(x) = ∏_{i<j}(t_j - t_i)。
#    用数值方法验证：
#      Vand <- function(tt) { n <- length(tt)
#        det(eigenMat(outer(0:(n-1), tt, function(p, t) t^p))) }
#      f <- function(x) Vand(c(1, 2, -2, x))
#      grid <- seq(-3, 3, by = 0.001)
#      v <- vapply(grid, f, numeric(1))
#      roots <- grid[which(diff(sign(v)) != 0)]      # 变号点
#      print(unique(round(roots, 2)))
#    再把 f(x) 与 12*(x-1)*(x-2)*(x+2) 对比，说明为什么恰好是这三个根。
#
# 4. Vandermonde 行列式 V(t1..tn) = ∏_{i<j}(t_j - t_i)。对 n = 4 用随机点验证，
#    并解释：为什么插值节点靠得很近时，解方程组会变得不可靠？
#    （提示：V 很小，但 b 里的舍入误差没有变小，所以解被放大。）
#
# 5. 用行列式判断四个点是否共面：给 4 个三维点 P1..P4，
#    det[P2-P1; P3-P1; P4-P1] = 0 ⇔ 共面。随机构造共面与不共面各一例验证。
#
# 思考题
#   (1) 教材用逆序数定义行列式。若改用「对角线法则」定义，能推广到 4 阶吗？
#       用代码算一个 4 阶矩阵的「对角线法则结果」，与真值比较。
#   (2) 为什么 ROpenEigen 提供 eigen_logdet()，而 det() 不算一个「更好的」接口？
#   (3) 三平面交于一点用 |A| ≠ 0 判定。若 |A| = 1e-14 但确实不为 0，
#       计算机给出的交点还可信吗？请设计一个实验说明。

say("\nlab01 结束。")
