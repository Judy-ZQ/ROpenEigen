# ===========================================================================
# 实验二 · 矩阵
#
# 对应教材：第二章 矩阵
#   第一节 矩阵的概念        第五节 初等变换与初等方阵
#   第二节 矩阵的基本运算    第六节 逆矩阵
#   第三节 分块矩阵          第七节 克拉默法则
#   第四节 矩阵的秩          第八节 矩阵的应用（投入产出模型）
# 建议学时：3 学时
#
# 实验目的
#   1. 熟悉矩阵运算，特别是乘法不满足交换律、消去律这两件事；
#   2. 会用分块矩阵做乘法，理解它和整体乘法结果相同；
#   3. 把「秩 = 最高阶非零子式的阶数」算出来，验证秩的常用性质；
#   4. 把初等行变换实现成代码，用它求秩、求逆、解方程组；
#   5. 用克拉默法则与投入产出模型各做一个完整算例。
#
# 运行：Rscript lab02-matrix.R   或在 RStudio 里 source()
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
# 一、概念与特殊矩阵（教材 §2.1）
# ===========================================================================
say("== 2.1 矩阵的概念 ==")
A <- eigenMat(matrix(c(1, 2, 3, 4, 5, 6), 2, 3, byrow = TRUE))
say("A 是 ", nrow(A), " x ", ncol(A), " 矩阵")
say("零矩阵 eigen_zeros(2,3)：")
print(eigen_zeros(2, 3))
say("单位矩阵 eigen_identity(3)：")
print(eigen_identity(3))
say("对角矩阵 diag(c(2,5,7))：")
print(eigenMat(diag(c(2, 5, 7))))

# 方阵的迹：主对角线之和（教材 §2.2 定义）
say("tr(diag(2,5,7)) = ", eigen_trace(eigenMat(diag(c(2, 5, 7)))))

# ===========================================================================
# 二、基本运算（教材 §2.2）
# ===========================================================================
say("\n== 2.2 基本运算 ==")
A <- eigenMat(matrix(c(1, 2, 3, 4, 5, 6), 2, 3, byrow = TRUE))
B <- eigenMat(matrix(c(6, 5, 4, 3, 2, 1), 2, 3, byrow = TRUE))
say("A + B =")
print(A + B)
say("2A - 3B =")
print(2 * A - 3 * B)

# 乘法要求内维相等：2x3 乘 3x2 得 2x2
C <- eigenMat(matrix(c(1, 0, 1, 0, 1, 1), 3, 2, byrow = TRUE))
say("A(2x3) %*% C(3x2) =")
print(A %*% C)

# (1) 乘法不可交换。先取两个同阶方阵
M <- eigenMat(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
N <- eigenMat(matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE))
say("MN ="); print(M %*% N)
say("NM ="); print(N %*% M)
say("MN == NM ？", relerr(as.matrix(M %*% N), as.matrix(N %*% M)) < 1e-12)

# (2) 没有消去律：AB = 0 推不出 A = 0 或 B = 0
Z1 <- eigenMat(matrix(c(1, 1, 1, 1), 2, 2, byrow = TRUE))
Z2 <- eigenMat(matrix(c(1, -1, -1, 1), 2, 2, byrow = TRUE))
say("A = [[1,1],[1,1]]，B = [[1,-1],[-1,1]]，AB =")
print(Z1 %*% Z2)
say("A、B 都不是零矩阵，但 AB = 0")

# (3) 转置的反序律 (AB)' = B'A'
say("(AB)' 与 B'A' 相等：",
    relerr(as.matrix(t(M %*% N)), as.matrix(t(N) %*% t(M))) < 1e-12)

# (4) 方阵的幂与矩阵多项式：f(A) = A^2 - 3A + 2E
A2 <- eigenMat(matrix(c(1, 2, 0, 1), 2, 2, byrow = TRUE))
say("A = [[1,2],[0,1]]")
say("A^2 ="); print(A2 %*% A2)
say("A^3 ="); print(A2 %*% A2 %*% A2)
fA <- A2 %*% A2 - 3 * A2 + 2 * eigen_identity(2)
say("f(A) = A^2 - 3A + 2E ="); print(fA)

# (5) 任意方阵都能写成一个对称阵与一个反对称阵之和
S <- 0.5 * (M + t(M)); K <- 0.5 * (M - t(M))
say("A = (A+A')/2 + (A-A')/2 还原误差：",
    format(relerr(as.matrix(S + K), as.matrix(M)), scientific = TRUE))
say("(A+A')/2 对称：", relerr(as.matrix(S), as.matrix(t(S))) < 1e-12,
    "；(A-A')/2 的迹为 0：", abs(eigen_trace(K)) < 1e-12)

# ===========================================================================
# 三、分块矩阵（教材 §2.3）
# ===========================================================================
say("\n== 2.3 分块矩阵 ==")
# 把 4x4 的 A 按 2x2 分成四块，用分块乘法算 A^2，再和整体相乘对比
set.seed(7)
A4 <- matrix(sample(-3:3, 16, replace = TRUE), 4, 4, byrow = TRUE)
A11 <- A4[1:2, 1:2]; A12 <- A4[1:2, 3:4]
A21 <- A4[3:4, 1:2]; A22 <- A4[3:4, 3:4]
blk_product <- rbind(
  cbind(A11 %*% A11 + A12 %*% A21, A11 %*% A12 + A12 %*% A22),
  cbind(A21 %*% A11 + A22 %*% A21, A21 %*% A12 + A22 %*% A22)
)
say("A4 =")
print(A4)
say("分块乘法的结果与整体 A4 %*% A4 的相对误差：",
    format(relerr(blk_product, A4 %*% A4), scientific = TRUE))

# 分块对角矩阵的幂：对角块各自求幂
D1 <- diag(c(2, 3)); D2 <- matrix(c(1, 1, 0, 1), 2, 2, byrow = TRUE)
BD <- rbind(cbind(D1, matrix(0, 2, 2)), cbind(matrix(0, 2, 2), D2))
BD2 <- BD %*% BD
say("分块对角阵的平方，右下块等于 D2^2：",
    relerr(BD2[3:4, 3:4], D2 %*% D2) < 1e-12)
say("左上块等于 D1^2：", relerr(BD2[1:2, 1:2], D1 %*% D1) < 1e-12)

# ===========================================================================
# 四、矩阵的秩（教材 §2.4）
# ===========================================================================
say("\n== 2.4 矩阵的秩 ==")
# 定义：最高阶非零子式的阶数。先按定义硬算，再和 eigen_rank 对账。
max_minor_order <- function(A, tol = 1e-9) {
  A <- as.matrix(A); m <- nrow(A); n <- ncol(A)
  for (k in min(m, n):1) {
    rs <- combn(m, k); cs <- combn(n, k)
    for (i in seq_len(ncol(rs))) for (j in seq_len(ncol(cs))) {
      if (abs(det(A[rs[, i], cs[, j], drop = FALSE])) > tol) return(k)
    }
  }
  0L
}
Ar <- matrix(c(1, 2, 2, 4, -2, 6, 1, 2, -1, 3, 3, 6), 3, 4, byrow = TRUE)
say("Ar =")
print(Ar)
say("按定义（枚举子式）求秩 = ", max_minor_order(Ar))
say("eigen_rank 求秩      = ", eigen_rank(eigenMat(Ar)))

# 秩的常用性质，逐条用随机矩阵验证
say("\n秩的性质（随机 20 组，全部满足才算通过）：")
set.seed(2024)
ok1 <- ok2 <- ok3 <- TRUE
for (t in 1:20) {
  m <- sample(2:5, 1); n <- sample(2:5, 1); k <- sample(2:4, 1)
  X <- matrix(sample(-2:2, m * n, TRUE), m, n)
  Y <- matrix(sample(-2:2, n * k, TRUE), n, k)
  ok1 <- ok1 && eigen_rank(eigenMat(X)) <= min(m, n)
  ok2 <- ok2 && eigen_rank(eigenMat(X %*% Y)) <= min(eigen_rank(eigenMat(X)),
                                                      eigen_rank(eigenMat(Y)))
  ok3 <- ok3 && eigen_rank(eigenMat(t(X))) == eigen_rank(eigenMat(X))
}
say("  R(A) <= min(m, n)：", ok1)
say("  R(AB) <= min(R(A), R(B))：", ok2)
say("  R(A') = R(A)：", ok3)

# 教材 §4.2 例 9 的结论：若 AB = 0，则 R(A) + R(B) <= n
say("\nAB = 0 时 R(A) + R(B) <= n 的检验：")
for (t in 1:5) {
  n <- 4
  Bm <- matrix(sample(-2:2, n * n, TRUE), n, n)
  # 把 B 的列换成 A 的零空间的向量，保证 AB = 0
  Nsp <- nullspace(Bm)
  if (ncol(Nsp) == 0) next
  Am <- t(Nsp)                                   # A 的行取零空间向量
  say(sprintf("  R(A)=%d R(B)=%d 和=%d <= n=%d ：%s",
              eigen_rank(eigenMat(Am)), eigen_rank(eigenMat(Bm)),
              eigen_rank(eigenMat(Am)) + eigen_rank(eigenMat(Bm)), n,
              eigen_rank(eigenMat(Am)) + eigen_rank(eigenMat(Bm)) <= n))
}

# ===========================================================================
# 五、初等变换与初等方阵（教材 §2.5）
# ===========================================================================
say("\n== 2.5 初等变换 ==")
# 三种初等行变换都可以写成「左乘一个初等方阵」：
#   E1 交换两行、E2 某行乘非零数、E3 某行加上另一行的 k 倍
E1 <- matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 1), 3, 3, byrow = TRUE)
E2 <- diag(c(1, 4, 1))
E3 <- diag(3); E3[3, 1] <- 5                      # r3 <- r3 + 5 r1
T <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 10), 3, 3, byrow = TRUE)
say("T ="); print(T)
say("E1 %*% T（交换第 1、2 行）："); print(E1 %*% T)
say("E2 %*% T（第 2 行乘 4）："); print(E2 %*% T)
say("E3 %*% T（第 3 行加第 1 行的 5 倍）："); print(E3 %*% T)
say("三个初等方阵都非奇异：", det(eigenMat(E1)) != 0, det(eigenMat(E2)) != 0,
    det(eigenMat(E3)) != 0)

# 用初等行变换把 T 化成行阶梯形 / 行最简形，同时读出秩
R <- rref(T)
say("T 的行最简形：")
print(R$rref)
say("主元所在列：", paste(R$pivots, collapse = ", "), "；秩 = ", R$rank)
say("与 eigen_rank 一致：", R$rank == eigen_rank(eigenMat(T)))

# ===========================================================================
# 六、逆矩阵（教材 §2.6）
# ===========================================================================
say("\n== 2.6 逆矩阵 ==")
Ainv <- matrix(c(1, 2, 3, 2, 5, 3, 1, 0, 8), 3, 3, byrow = TRUE)
say("A ="); print(Ainv)
say("|A| = ", det(eigenMat(Ainv)), "（非零 → 可逆）")

# 方法一：伴随矩阵法 A^{-1} = adj(A)/|A|
inv_adj <- adjugate(Ainv) / det(eigenMat(Ainv))
# 方法二：初等变换法 [A | E] -> [E | A^{-1}]
aug <- rref(cbind(Ainv, diag(3)))
inv_gauss <- aug$rref[, 4:6]
# 方法三：库函数
inv_lib <- as.matrix(solve(eigenMat(Ainv)))
say("伴随矩阵法与库函数的相对误差：", format(relerr(inv_adj, inv_lib), scientific = TRUE))
say("初等变换法与库函数的相对误差：", format(relerr(inv_gauss, inv_lib), scientific = TRUE))
say("A %*% A^{-1} 与 E 的偏差：", format(relerr(Ainv %*% inv_lib, diag(3)), scientific = TRUE))

# 逆的性质
Bm <- matrix(c(2, 0, 1, 1, 1, 0, 0, 1, 3), 3, 3, byrow = TRUE)
say("(AB)^{-1} = B^{-1}A^{-1}：",
    relerr(solve(eigenMat(Ainv %*% Bm)), as.matrix(solve(eigenMat(Bm)) %*% solve(eigenMat(Ainv)))) < 1e-10)
say("(A')^{-1} = (A^{-1})'：",
    relerr(solve(eigenMat(t(Ainv))), t(solve(eigenMat(Ainv)))) < 1e-10)

# 不可逆的情形——这里有一处必须记牢的差异：
#   base R 的 solve() 遇到奇异矩阵会报 "system is computationally singular"；
#   ROpenEigen 的 solve() 把活交给 Eigen 的 LU，而 LU 不做奇异检测，
#   于是它照样算下去，返回一张 NaN/Inf 的矩阵。不报错，但结果全是垃圾。
Sing <- matrix(c(1, 2, 3, 2, 4, 6, 1, 1, 1), 3, 3, byrow = TRUE)
say("奇异矩阵 |A| = ", det(eigenMat(Sing)), "，rcond = ",
    format(eigen_rcond(eigenMat(Sing)), scientific = TRUE))
say("直接 solve() 得到的：")
print(as.matrix(solve(eigenMat(Sing))))
say("base R 的 solve() 则会直接报错：")
r0 <- try(solve(Sing), silent = TRUE)
say("  -> ", if (inherits(r0, "try-error")) "报错（正确行为）" else "没有报错")
# 所以求逆前先看 rcond；helpers.R 里的 safe_inverse() 就是这道闸门
say("safe_inverse() 的反应：")
r <- try(safe_inverse(Sing), silent = TRUE)
say("  -> ", if (inherits(r, "try-error")) conditionMessage(attr(r, "condition")) else "算出了结果")

# ===========================================================================
# 七、克拉默法则（教材 §2.7）
# ===========================================================================
say("\n== 2.7 克拉默法则 ==")
# 若 |A| ≠ 0，则 Ax = b 的唯一解是 x_j = |A_j| / |A|，A_j 是把第 j 列换成 b
cramer <- function(A, b) {
  D <- det(eigenMat(A))
  if (abs(D) < 1e-12) stop("系数行列式为 0，克拉默法则失效")
  n <- ncol(A)
  vapply(seq_len(n), function(j) {
    Aj <- A; Aj[, j] <- b
    det(eigenMat(Aj)) / D
  }, numeric(1))
}
Ac <- matrix(c(2, 1, -1, 1, -1, 1, 1, 1, 2), 3, 3, byrow = TRUE)
bc <- c(4, 1, 5)
say("A ="); print(Ac)
say("b = ", paste(bc, collapse = ", "), "；|A| = ", det(eigenMat(Ac)))
say("克拉默法则：", paste(round(cramer(Ac, bc), 10), collapse = ", "))
say("库函数求解：", paste(round(as.vector(solve(eigenMat(Ac), eigenVec(bc))), 10), collapse = ", "))
say("回代验证 Ax - b 的最大偏差：",
    format(max(abs(Ac %*% cramer(Ac, bc) - bc)), scientific = TRUE))

# 行列式为 0 时法则失效：方程组可能无解，也可能有无穷多解
As <- matrix(c(1, 2, 3, 2, 4, 6, 1, 1, 1), 3, 3, byrow = TRUE)
say("奇异情形 |A| = ", det(eigenMat(As)))
say("  b = c(1,2,3) 时 rank(A)=", eigen_rank(eigenMat(As)),
    "，rank(A|b)=", eigen_rank(eigenMat(cbind(As, c(1, 2, 3)))),
    " → 两者相等，无穷多解")
say("  b = c(1,3,1) 时 rank(A)=", eigen_rank(eigenMat(As)),
    "，rank(A|b)=", eigen_rank(eigenMat(cbind(As, c(1, 3, 1)))),
    " → 不相等，无解")

# ===========================================================================
# 八、应用：投入产出模型（教材 §2.8 定义 2.8.1）
# ===========================================================================
say("\n== 2.8 投入产出模型 ==")
# 三个部门：农业、工业、服务业。a_ij = 第 j 部门生产 1 元产品，直接消耗第 i 部门的产品量。
Aio <- matrix(c(0.20, 0.30, 0.10,
                0.30, 0.20, 0.20,
                0.10, 0.20, 0.30), 3, 3, byrow = TRUE)
Y <- c(100, 80, 60)                    # 最终产品（消费 + 积累 + 出口）
say("直接消耗系数矩阵 A =")
print(Aio)
say("各列之和（应小于 1，模型才有意义）：", paste(colSums(Aio), collapse = ", "))

# 平衡方程 X = A X + Y  →  X = (E - A)^{-1} Y
EminusA <- diag(3) - Aio
say("|E - A| = ", det(eigenMat(EminusA)))
X <- as.vector(solve(eigenMat(EminusA), eigenVec(Y)))
say("总产出 X = ", paste(round(X, 4), collapse = ", "))
say("平衡方程残差 max|X - (AX + Y)| = ",
    format(max(abs(X - (Aio %*% X + Y))), scientific = TRUE))

# 完全消耗系数矩阵 B = (E - A)^{-1} - E：b_ij 是 j 部门多生产 1 元最终产品，
# 对 i 部门的完全消耗（含间接消耗）。
Bfull <- as.matrix(solve(eigenMat(EminusA))) - diag(3)
say("完全消耗系数矩阵 B =")
print(round(Bfull, 4))
say("B 的每一列都大于相等位置的 A 列（间接消耗必然增加总消耗）：",
    all(Bfull >= Aio - 1e-12))

# 若最终产品需求变化 ΔY，总产出的变化 ΔX = (E - A)^{-1} ΔY
dY <- c(10, 0, 0)
dX <- as.vector(solve(eigenMat(EminusA), eigenVec(dY)))
say("农业最终产品增加 10 时，三个部门总产出分别增加：",
    paste(round(dX, 4), collapse = ", "))
say("注意 dX 的农业分量 ", round(dX[1], 4),
    " 大于 10：农业为了多交 10 元的最终产品，还得给本部门和其他部门提供中间产品，",
    "这就是乘数效应，它不是简单相加。")

# ===========================================================================
# 九、实验习题
# ===========================================================================
# 1. 举出两个 2 阶方阵 A、B，使 (A+B)^2 ≠ A^2 + 2AB + B^2，并说明展开式在
#    什么条件下才等于右边。用代码验证你的答案。
#
# 2. 教材 §2.4 用「最高阶非零子式」定义秩。对
#        A = [[1, 2, 2, 4, -2, 6], [1, 2, -1, 3, 3, 6]]
#    求出所有 2 阶子式，数一数非零的有几个，再解释为什么 R(A) = 2。
#
# 3. 用初等行变换求下面矩阵的逆，并与 solve() 对比：
#        A = [[ 1, 2, 3], [2, 5, 3], [1, 0, 8]]
#    再回答：为什么 [A|E] 化成 [E|A^{-1}] 这一步要求 A 可逆？
#
# 4. 克拉默法则的代价是 n+1 个 n 阶行列式，即 O(n·n!)。写一个计时实验，
#    比较 n = 6, 8, 10 时克拉默法则与 solve() 的耗时（可用 system.time）。
#
# 5. 投入产出模型：若某部门的技术改造使 a_11 从 0.20 降到 0.15，
#    在最终产品 Y 不变的前提下，总产出 X 会怎么变？用上面的算例算一遍。
#
# 思考题
#   (1) 教材指出 AB = 0 不能推出 A = 0 或 B = 0。从线性变换的角度看，
#       这句话在说什么？（提示：B 的列都在 A 的零空间里。）
#   (2) 为什么程序里判断矩阵可逆不看 |A| 是否等于 0，而看 rcond 或秩？
#   (3) 分块乘法 A11B11 + A12B21 的项数与整体乘法一样多，那分块的意义在哪？
#       （提示：稀疏、并行、以及把高阶矩阵化为低阶问题。）

say("\nlab02 结束。")
