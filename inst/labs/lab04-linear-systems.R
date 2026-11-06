# ===========================================================================
# 实验四 · 线性方程组
#
# 对应教材：第四章 线性方程组
#   第一节 线性方程组的同解定理
#   第二节 齐次线性方程组
#   第三节 非齐次线性方程组
#   第四节 线性方程组的应用
# 建议学时：3 学时
#
# 实验目的
#   1. 用增广矩阵的行最简形一次读出「有解 / 无解 / 无穷多解」；
#   2. 求齐次方程组的基础解系，验证解空间维数 = n - R(A)；
#   3. 写出非齐次方程组的通解 x = 特解 + 齐次通解，并验证它确实是全部解；
#   4. 处理不相容方程组：最小二乘解与最小范数解；
#   5. 完成两个应用算例：食堂就餐比例的马尔可夫模型、化学方程式配平。
#
# 运行：Rscript lab04-linear-systems.R   或在 RStudio 里 source()
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
# 一、同解定理与有解判定（教材 §4.1）
# ===========================================================================
say("== 4.1 同解定理 ==")
# 对增广矩阵做初等行变换，解集不变；所以把 [A|b] 化成行最简形就能读出一切。
A <- matrix(c(1, 2, -1, 2, 1, 1, -1, 1, 2), 3, 3, byrow = TRUE)
b <- c(2, 5, 3)
Ab <- cbind(A, b)
say("A ="); print(A)
say("b = ", paste(b, collapse = ", "))
R <- rref(Ab)
say("[A|b] 的行最简形：")
print(R$rref)

rA <- eigen_rank(eigenMat(A))
rAb <- eigen_rank(eigenMat(Ab))
say("R(A) = ", rA, "，R(A|b) = ", rAb, "，未知数 n = ", ncol(A))
say("判定：", if (rA < rAb) "无解" else if (rA == ncol(A)) "唯一解" else "无穷多解")
say("实际解：", paste(round(as.vector(solve(eigenMat(A), eigenVec(b))), 10), collapse = ", "))

# 三种情形各给一个例子，用同一套判据机器判一遍
cases <- list(
  list(name = "唯一解", A = matrix(c(2, 1, 1, 3), 2, 2, byrow = TRUE), b = c(3, 5)),
  list(name = "无穷多解", A = matrix(c(1, 2, 2, 4), 2, 2, byrow = TRUE), b = c(1, 2)),
  list(name = "无解", A = matrix(c(1, 2, 2, 4), 2, 2, byrow = TRUE), b = c(1, 3))
)
say("\n三种情形：")
for (cs in cases) {
  rr <- eigen_rank(eigenMat(cs$A)); ra <- eigen_rank(eigenMat(cbind(cs$A, cs$b)))
  verdict <- if (rr < ra) "无解" else if (rr == ncol(cs$A)) "唯一解" else "无穷多解"
  say(sprintf("  %-8s R(A)=%d R(A|b)=%d → %s", cs$name, rr, ra, verdict))
}

# ===========================================================================
# 二、齐次线性方程组（教材 §4.2）
# ===========================================================================
say("\n== 4.2 齐次方程组 ==")
# 第 3 行 = 第 1 行 + 第 2 行，所以秩为 2，自由变量有 2 个
Ah <- matrix(c(1, 2, 2, -1,
               2, 4, 8, 0,
               3, 6, 10, -1), 3, 4, byrow = TRUE)
say("A ="); print(Ah)
rh <- eigen_rank(eigenMat(Ah))
N <- nullspace(Ah)
say("R(A) = ", rh, "，n = ", ncol(Ah), "，解空间维数 = n - R(A) = ", ncol(Ah) - rh)
say("基础解系（按列）：")
print(round(N, 10))
say("A %*% N 应为零矩阵：", max(abs(Ah %*% N)) < 1e-12)
say("基础解系自身线性无关：", eigen_rank(eigenMat(N)) == ncol(N))

# 通解：x = c1 ξ1 + c2 ξ2，取几组 c 验证
say("\n通解 x = c1ξ1 + c2ξ2，随机取系数验证：")
set.seed(5)
for (t in 1:3) {
  cc <- round(runif(ncol(N), -2, 2), 1)
  x <- as.vector(N %*% cc)
  say(sprintf("  c = (%s) → x = (%s)，|Ax| = %.2e",
              paste(sprintf("%.1f", cc), collapse = ", "),
              paste(sprintf("%.1f", x), collapse = ", "),
              sqrt(sum((Ah %*% x)^2))))
}

# 若 R(A) = n，则只有零解
Az <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 10), 3, 3, byrow = TRUE)
say("满秩时基础解系为空：", ncol(nullspace(Az)) == 0,
    "（R(A) = ", eigen_rank(eigenMat(Az)), " = n = 3）")

# ===========================================================================
# 三、非齐次线性方程组（教材 §4.3）
# ===========================================================================
say("\n== 4.3 非齐次方程组 ==")
An <- matrix(c(1, 2, 2, -1,
               2, 4, 8, 0,
               3, 6, 10, -1), 3, 4, byrow = TRUE)
bn <- c(1, 2, 3)
say("A ="); print(An)
say("b = ", paste(bn, collapse = ", "))
say("R(A) = ", eigen_rank(eigenMat(An)),
    "，R(A|b) = ", eigen_rank(eigenMat(cbind(An, bn))), " → 有解")

xp <- particular_solution(An, bn)
say("一个特解 x* = ", paste(round(xp, 10), collapse = ", "))
say("A x* = ", paste(round(as.vector(An %*% xp), 10), collapse = ", "), "（应等于 b）")

Nh <- nullspace(An)
say("导出组（齐次）的基础解系含 ", ncol(Nh), " 个向量")
say("通解 x = x* + N c，随机取 c 验证：")
for (t in 1:3) {
  cc <- round(runif(ncol(Nh), -3, 3), 1)
  x <- xp + as.vector(Nh %*% cc)
  say(sprintf("  c = (%s) → max|Ax - b| = %.2e",
              paste(sprintf("%.1f", cc), collapse = ", "),
              max(abs(An %*% x - bn))))
}

# 用库函数解一遍对账（这里方程组有唯一性吗？4 个未知数 3 个方程 → 无穷多解）
say("库函数给出的一个解（最小范数解）：")
say("  ", paste(round(as.vector(eigen_pinv(eigenMat(An)) %*% eigenVec(bn)), 10), collapse = ", "))
say("  它的模 = ", sqrt(sum(as.vector(eigen_pinv(eigenMat(An)) %*% eigenVec(bn))^2)),
    "，小于等于其他特解的模（最小范数解的定义）")

# 无解的情形：rank 不相等
say("\n无解例子：b 改成 c(1, 2, 4)")
say("  R(A) = ", eigen_rank(eigenMat(An)),
    "，R(A|b) = ", eigen_rank(eigenMat(cbind(An, c(1, 2, 4)))), " → 无解")
# 注意：方程组无解时，ROpenEigen 的 solve() 不会报错——对长方形（超定）矩阵
# 它交给 Eigen 做 QR 最小二乘，返回的是「残差最小」的那个 x，而不是解。
# base R 的 solve() 遇到非方阵/奇异矩阵会直接拒绝。
b_no <- c(1, 2, 4)
r_ls <- solve(eigenMat(An), eigenVec(b_no))
say("  ROpenEigen 返回的 x = ", paste(round(as.vector(r_ls), 6), collapse = ", "))
say("  代回后 A x - b = ", paste(round(as.vector(An %*% as.vector(r_ls) - b_no), 6), collapse = ", "))
say("  它不是解（残差不为 0），而是最小二乘解 —— 这正是下一节的内容。")

# ===========================================================================
# 四、不相容方程组：最小二乘解（拓展）
# ===========================================================================
say("\n== 4.4 最小二乘 ==")
# 超定方程组：4 个方程 2 个未知数，一般无解，找一个使 ||Ax-b|| 最小的 x
A1 <- cbind(1, c(0, 1, 2, 3))
b1 <- c(1, 3, 5, 6)                     # 大致在直线 y = 1 + 1.8x 附近，末点有偏差
say("A ="); print(A1)
say("b = ", paste(b1, collapse = ", "))

# 方法一：正规方程 A'A x = A'b
x_ne <- as.vector(solve(eigenMat(t(A1) %*% A1), eigenVec(t(A1) %*% b1)))
# 方法二：QR 分解（数值上更稳）
x_qr <- as.vector(solve(eigenMat(A1), eigenVec(b1)))      # ROpenEigen 对超定做 QR 求解
# 方法三：SVD 的伪逆
x_svd <- as.vector(eigen_pinv(eigenMat(A1)) %*% eigenVec(b1))

say("正规方程解：", paste(round(x_ne, 8), collapse = ", "))
say("QR 解      ：", paste(round(x_qr, 8), collapse = ", "))
say("SVD 伪逆解 ：", paste(round(x_svd, 8), collapse = ", "))

# 三种解应当一致；再验证「残差平方和最小」这个定义
resid <- function(x) sum((A1 %*% x - b1)^2)
say("残差平方和：正规方程 ", round(resid(x_ne), 10),
    "，QR ", round(resid(x_qr), 10), "，SVD ", round(resid(x_svd), 10))
set.seed(9)
best <- min(replicate(2000, resid(x_ne + rnorm(2, 0, 0.5))))
say("在解附近随机扰动 2000 次，最小残差 = ", round(best, 6),
    "，不小于最优值：", best >= resid(x_ne) - 1e-12)

# 正规方程的条件数是原矩阵的平方，这正是它不如 QR/SVD 的原因
say("cond(A'A) = ", 1 / eigen_rcond(eigenMat(t(A1) %*% A1)),
    "，cond(A) = ", 1 / eigen_rcond(eigenMat(A1)),
    "（平方关系）")

# ===========================================================================
# 五、应用一：食堂就餐比例的马尔可夫模型（教材 §4.4 例 1）
# ===========================================================================
say("\n== 4.5 应用：食堂就餐比例 ==")
# 原在 A 就餐的学生，下月 1/4 留 A、3/4 去 B；原在 B 的，1/3 来 A、2/3 留 B。
# 于是 x_n = P^n x_0，其中 P = [[1/4, 1/3], [3/4, 2/3]]
P <- matrix(c(1/4, 1/3,
              3/4, 2/3), 2, 2, byrow = TRUE)
say("转移矩阵 P ="); print(round(P, 6))
say("P 的列和都是 1（概率守恒）：", paste(colSums(P), collapse = ", "))

# 矩阵幂：x_n = P^n x_0。写成一个循环，别手写一长串乘法。
matpow <- function(M, k) {
  R <- diag(nrow(M)); Mk <- M
  while (k > 0) {
    if (k %% 2 == 1) R <- R %*% Mk
    Mk <- Mk %*% Mk; k <- k %/% 2
  }
  R
}
x0 <- c(0.6, 0.4)
say("初始比例 x0 = (0.6, 0.4)")
for (n in c(1, 2, 3, 5, 10, 30)) {
  xn <- as.vector(matpow(P, n) %*% x0)
  say(sprintf("  第 %2d 个月后：A = %.6f，B = %.6f", n, xn[1], xn[2]))
}

# 稳态：P 的不动点，解 (P - E)x = 0 且 x1 + x2 = 1
# 把「和为 1」作为附加方程，替换掉一个线性相关的方程
Sys <- rbind(P - diag(2), c(1, 1))
rhs <- c(0, 0, 1)
xs <- rref(Sys)$rref
say("稳态方程 (P-E)x = 0，x1+x2 = 1 的解：")
xsol <- particular_solution(Sys, rhs)
say("  x* = (", paste(round(xsol, 10), collapse = ", "), ")")
say("  即 A 占 ", round(xsol[1], 10), " = 4/13，B 占 ", round(xsol[2], 10), " = 9/13")
say("  验证 P x* = x*：", relerr(P %*% xsol, xsol) < 1e-12)

# ===========================================================================
# 六、应用二：化学方程式配平
# ===========================================================================
say("\n== 4.6 应用：配平化学方程式 ==")
# 丙烷燃烧：x1 C3H8 + x2 O2 → x3 CO2 + x4 H2O
# 按元素守恒列方程（生成物移到左边取负）：
#   C: 3x1 - x3 = 0
#   H: 8x1 - 2x4 = 0
#   O: 2x2 - 2x3 - x4 = 0
Achem <- matrix(c(3, 0, -1, 0,
                  8, 0, 0, -2,
                  0, 2, -2, -1), 3, 4, byrow = TRUE)
say("元素守恒方程组 A x = 0，A =")
print(Achem)
Nc <- nullspace(Achem)
say("基础解系：", paste(round(Nc[, 1], 10), collapse = ", "))
# 化成最小正整数解
v <- Nc[, 1]
v <- v / min(abs(v[v != 0]))
mult <- 1
while (any(abs(v * mult - round(v * mult)) > 1e-9)) mult <- mult + 1
v <- round(v * mult)
say("最小正整数解：x = ", paste(v, collapse = ", "))
say("即 ", v[1], " C3H8 + ", v[2], " O2 → ", v[3], " CO2 + ", v[4], " H2O")
say("逐元素核对：C ", 3 * v[1], " = ", v[3],
    "；H ", 8 * v[1], " = ", 2 * v[4],
    "；O ", 2 * v[2], " = ", 2 * v[3] + v[4])

# ===========================================================================
# 七、实验习题
# ===========================================================================
# 1. 解方程组并写出通解：
#       x1 + x2 - 3x3 - x4 = 1
#       x1 - x2 - x3 + 3x4 = 3
#       3x1 + x2 - 5x3 - 3x4 = 5
#    要求：分别用 rref(增广) 和 nullspace() + particular_solution() 两条路径，
#    并说明解空间的维数为什么是 n - R(A)。
#
# 2. 讨论 λ 取何值时下列方程组有唯一解、无解、无穷多解：
#       λx1 + x2 + x3 = 1
#       x1 + λx2 + x3 = 1
#       x1 + x2 + λx3 = 1
#    提示：先算 |A| = (λ-1)^2(λ+2)，再分别取 λ = 1、λ = -2、其他值用代码验证。
#
# 3. 最小二乘：给 5 个点，拟合二次多项式 y = a + bx + cx^2。
#    用正规方程与 SVD 两种方法求系数，比较残差与条件数。
#    再把 x 的取值改成 0, 1, 2, ..., 4 与 0, 0.01, 0.02, ..., 0.04，
#    观察条件数的变化，解释为什么「自变量尺度」会影响结果。
#
# 4. 教材 §4.4 的食堂模型：若学校新建第三食堂 C，转移矩阵改为 3 阶，
#    自己设定一组合理的转移概率（列和为 1），求稳态分布，
#    并说明为什么稳态与初始分布无关。
#
# 5. 配平：KMnO4 + HCl → KCl + MnCl2 + Cl2 + H2O。
#    （提示：6 个未知数，按 K、Mn、O、H、Cl 列 5 个方程，基础解系 1 维。）
#
# 思考题
#   (1) 为什么「R(A) = R(A|b)」就保证有解？从列空间的角度说明。
#   (2) 通解写成 x* + Nc，其中的 x* 不唯一。不同的 x* 之间差什么？
#   (3) 最小二乘解在最优点处残差 r = Ax - b 与 A 的列空间垂直。用代码验证这一点。

say("\nlab04 结束。")
