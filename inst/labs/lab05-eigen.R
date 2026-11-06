# ===========================================================================
# 实验五 · 矩阵的相似对角化
#
# 对应教材：第五章
#   第一节 方阵的特征值与特征向量
#   第二节 相似矩阵
#   第三节 实对称矩阵的相似对角化
#   第四节 特征值理论的应用（差分方程、动态模型）
# 建议学时：3 学时
#
# 实验目的
#   1. 用 Av = λv 直接验证特征值与特征向量，并用特征多项式核对；
#   2. 分清「代数重数」与「几何重数」，据此判断能不能对角化；
#   3. 验证相似不变量（行列式、迹、特征多项式、秩）；
#   4. 完成实对称矩阵的正交对角化 Q'ΛQ，并说明特征向量为什么不唯一；
#   5. 用对角化求差分方程 x_{k+1} = A x_k 的通项与长期行为。
#
# 运行：Rscript lab05-eigen.R   或在 RStudio 里 source()
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
# 一、特征值与特征向量（教材 §5.1）
# ===========================================================================
say("== 5.1 特征值与特征向量 ==")
A <- matrix(c(3, -1, -1, 3), 2, 2, byrow = TRUE)
say("A ="); print(A)

# 定义：Av = λv。先手算再让程序算，两边对账。
say("特征多项式 |λE - A| = (λ-3)^2 - 1 = λ^2 - 6λ + 8 = (λ-2)(λ-4)")
say("所以特征值 λ1 = 2，λ2 = 4")
ev <- eigen_eigen(eigenMat(A))
lam <- as.vector(ev$values()); V <- as.matrix(ev$vectors())
say("eigen_eigen 给出：λ = ", paste(round(lam, 10), collapse = ", "),
    "（升序，这是 SelfAdjointEigenSolver 的约定）")
say("对应的特征向量（按列）：")
print(round(V, 10))

# 逐个验证 Av = λv：这是唯一不会骗人的判据
say("验证 A v_i - λ_i v_i 的最大模：")
for (i in seq_along(lam)) {
  say(sprintf("  i=%d：%.2e", i, max(abs(A %*% V[, i] - lam[i] * V[, i]))))
}

# 迹与行列式的两条恒等式（教材 §5.1 性质）
say("tr(A) = ", sum(diag(A)), " = Σλ = ", sum(lam))
say("|A| = ", det(eigenMat(A)), " = Πλ = ", prod(lam))

# 一般（非对称）矩阵的复特征值：旋转 90°
N <- matrix(c(0, -1, 1, 0), 2, 2, byrow = TRUE)
say("\n旋转 90° 的矩阵 [[0,-1],[1,0]]：")
say("  eigen_eigen（只读下三角，要求对称）：", paste(round(as.vector(eigen_eigen(eigenMat(N))$values()), 6), collapse = ", "),
    " ← 这是错的结果")
g <- eigen_geigen(eigenMat(N))
say("  eigen_geigen（一般矩阵，返回复数）：", paste(round(as.vector(g$values()), 6), collapse = ", "))
say("  特征值是 ±i：", relerr(Re(as.vector(g$values())), c(0, 0)) < 1e-12,
    "，虚部 ", paste(round(Im(as.vector(g$values())), 10), collapse = ", "))

# ===========================================================================
# 二、特征多项式（教材 §5.1）
# ===========================================================================
say("\n== 5.2 特征多项式 ==")
# |λE - A| = λ^n + c1 λ^(n-1) + ... + cn。Faddeev-LeVerrier 递推算系数。
B <- matrix(c(1, 2, 3, 2, 5, 3, 1, 0, 8), 3, 3, byrow = TRUE)
say("B ="); print(B)
cf <- charpoly_faddeev(B)
say("特征多项式系数（降幂）：", paste(round(cf, 10), collapse = ", "))
say("即 φ(λ) = λ^3 + (", round(cf[2], 8), ")λ^2 + (", round(cf[3], 8),
    ")λ + (", round(cf[4], 8), ")")
say("用 det(λE - B) 在 λ = 5 处核对：",
    round(det(eigenMat(5 * diag(3) - B)), 8), " vs ",
    round(sum(cf * 5^(3:0)), 8))

# 特征多项式的根就是特征值；用 R 的 polyroot 求根，与 eigen_geigen 对比
rts <- polyroot(rev(cf))
say("polyroot 求得的根：", paste(round(rts, 8), collapse = ", "))
say("eigen_geigen 给出的：", paste(round(as.vector(eigen_geigen(eigenMat(B))$values()), 8),
                                   collapse = ", "))

# ===========================================================================
# 三、能不能对角化：代数重数 vs 几何重数（教材 §5.1、§5.2）
# ===========================================================================
say("\n== 5.3 可对角化的判定 ==")
# Jordan 块 J = [[3,1],[0,3]]：λ = 3 是二重根，但特征向量只有 1 个
J <- matrix(c(3, 1, 0, 3), 2, 2, byrow = TRUE)
say("J = [[3,1],[0,3]]")
lamJ <- as.vector(eigen_geigen(eigenMat(J))$values())
say("  特征值：", paste(round(lamJ, 10), collapse = ", "), " → 代数重数 2")
say("  rank(J - 3E) = ", eigen_rank(eigenMat(J - 3 * diag(2))),
    "，特征子空间维数 = n - rank = ", 2 - eigen_rank(eigenMat(J - 3 * diag(2))),
    " → 几何重数 1")
say("  几何重数 < 代数重数 → 不可对角化")
say("  证据：nullspace(J - 3E) 只有 ", ncol(nullspace(J - 3 * diag(2))), " 个向量")

# 可对角化的例子：重根但特征子空间够了
K <- matrix(c(3, 0, 0, 3), 2, 2, byrow = TRUE)
say("\nK = 3E：λ = 3 的代数重数 2，rank(K-3E) = 0 → 几何重数 2 → 可对角化")

# 一般判定：n 阶方阵可对角化 ⇔ 每个特征值的几何重数 = 代数重数
say("\n一般判定函数：")
is_diagonalizable <- function(M, tol = 1e-8) {
  lam <- as.vector(eigen_geigen(eigenMat(M))$values())
  # 把相近的特征值聚成一组，比较每组的多重数与零空间维数
  for (l in unique(round(lam, 6))) {
    alg <- sum(abs(lam - l) < tol)
    geo <- ncol(nullspace(as.matrix(M) - l * diag(nrow(M))))
    if (geo < alg) return(FALSE)
  }
  TRUE
}
say("  J 可对角化？", is_diagonalizable(J))
say("  K 可对角化？", is_diagonalizable(K))
say("  A = [[3,-1],[-1,3]] 可对角化？", is_diagonalizable(A))

# ===========================================================================
# 四、相似矩阵与不变量（教材 §5.2）
# ===========================================================================
say("\n== 5.4 相似矩阵 ==")
# 取一个可逆 P，构造 B = P^{-1} A P，验证不变量
P <- matrix(c(1, 2, 1, 1), 2, 2, byrow = TRUE)
Bs <- as.matrix(solve(eigenMat(P)) %*% eigenMat(A) %*% eigenMat(P))
say("P = [[1,2],[1,1]]，B = P^{-1} A P =")
print(round(Bs, 10))
say("A 与 B 不相等：", relerr(Bs, A) > 1e-6)
say("相似不变量：")
say("  |A| = ", round(det(eigenMat(A)), 10), "，|B| = ", round(det(eigenMat(Bs)), 10))
say("  tr(A) = ", round(eigen_trace(eigenMat(A)), 10), "，tr(B) = ", round(eigen_trace(eigenMat(Bs)), 10))
say("  R(A) = ", eigen_rank(eigenMat(A)), "，R(B) = ", eigen_rank(eigenMat(Bs)))
say("  特征值相同：",
    relerr(sort(as.vector(eigen_eigen(eigenMat(A))$values())),
           sort(as.vector(eigen_eigen(eigenMat(Bs))$values()))) < 1e-12)

# 相似但不等：A 与对角阵 Λ 相似 ⇔ A 可对角化
Lam <- diag(lam)
say("A 与 diag(2,4) 相似：A = V Λ V^{-1} 还原误差 ",
    format(relerr(V %*% Lam %*% solve(eigenMat(V)), A), scientific = TRUE))

# ===========================================================================
# 五、实对称矩阵的正交对角化（教材 §5.3）
# ===========================================================================
say("\n== 5.5 实对称矩阵的正交对角化 ==")
# 教材 §5.3 的例子：A = [[4,2,2],[2,4,2],[2,2,4]]，特征值 8, 2, 2
S <- matrix(c(4, 2, 2,
              2, 4, 2,
              2, 2, 4), 3, 3, byrow = TRUE)
say("S ="); print(S)
es <- eigen_eigen(eigenMat(S))
lamS <- as.vector(es$values()); Q <- as.matrix(es$vectors())
say("特征值：", paste(round(lamS, 10), collapse = ", "), "（即 2, 2, 8）")
say("特征向量（按列，已是规范正交）：")
print(round(Q, 10))
say("Q'Q = E：", relerr(t(Q) %*% Q, diag(3)) < 1e-12)
say("Q'SQ = Λ：", relerr(t(Q) %*% S %*% Q, diag(lamS)) < 1e-12)
say("S = QΛQ' 还原误差：", format(relerr(Q %*% diag(lamS) %*% t(Q), S), scientific = TRUE))

# 教材手算时取 p1 = (1,1,1)/√3（λ=8），p2、p3 是 λ=2 的特征子空间里的一组正交基。
# 程序给出的 Q 是同一子空间里的另一组正交基，两者不一定逐列相同。
p1 <- c(1, 1, 1) / sqrt(3)
say("手工取的 p1 = ", paste(round(p1, 6), collapse = ", "))
say("它确实是 λ=8 的特征向量：", max(abs(S %*% p1 - 8 * p1)) < 1e-12)
say("与程序给出的 λ=8 的特征向量夹角余弦（绝对值）应接近 1：",
    round(abs(sum(p1 * Q[, 3])), 12))

# 重特征值的特征子空间：λ=2 对应二维子空间，任取一组正交基都可以
idx2 <- which(abs(lamS - 2) < 1e-8)
say("λ=2 的特征子空间维数 = ", length(idx2),
    "（代数重数也是 2，所以可正交对角化）")
say("该子空间内任意两个向量的内积为 0：",
    abs(sum(Q[, idx2[1]] * Q[, idx2[2]])) < 1e-12)

# 随机对称矩阵的通用检验
say("\n随机 5 阶对称矩阵的正交对角化检验（10 组）：")
set.seed(17)
allok <- TRUE
for (t in 1:10) {
  M <- matrix(rnorm(25), 5, 5); M <- (M + t(M)) / 2
  em <- eigen_eigen(eigenMat(M))
  Qm <- as.matrix(em$vectors()); lm <- as.vector(em$values())
  allok <- allok && relerr(t(Qm) %*% M %*% Qm, diag(lm)) < 1e-10 &&
    relerr(t(Qm) %*% Qm, diag(5)) < 1e-12
}
say("  全部通过：", allok)

# ===========================================================================
# 六、应用：差分方程与长期行为（教材 §5.4）
# ===========================================================================
say("\n== 5.6 应用：差分方程 x_{k+1} = A x_k ==")
# 若 A = P Λ P^{-1}，则 x_k = A^k x_0 = P Λ^k P^{-1} x_0 —— 幂运算变成数的幂。
M2 <- matrix(c(0.9, 0.2,
               0.1, 1.05), 2, 2, byrow = TRUE)
say("M ="); print(M2)
em2 <- eigen_geigen(eigenMat(M2))
lm2 <- Re(as.vector(em2$values()))
say("特征值：", paste(round(lm2, 10), collapse = ", "))
say("主特征值（模最大）= ", round(max(abs(lm2)), 10),
    " → 长期每年的增长率趋近这个数")

x0 <- c(100, 50)
say("x0 = (100, 50)，用 A^k x0 逐年迭代：")
xk <- x0
for (k in 1:8) {
  xk <- as.vector(M2 %*% xk)
  if (k <= 5 || k == 8) say(sprintf("  k=%d：(%8.3f, %8.3f)  比值 = %.6f",
                                    k, xk[1], xk[2], xk[1] / xk[2]))
}

# 与谱分解的结果对照：x_k = V Λ^k V^{-1} x_0，幂运算变成数的幂
# （特征向量可能带复数，所以整条链路都按复数算，最后取实部）
ev2 <- eigen_geigen(eigenMat(M2))
Vc <- as.matrix(ev2$vectors())
lam_c <- as.vector(ev2$values())
c2 <- as.vector(solve(Vc) %*% x0)              # x0 在特征向量基下的坐标
say("按特征向量分解的系数：", paste(round(Re(c2), 6), collapse = ", "))
xk_spec <- Re(Vc %*% (lam_c^8 * c2))
say("第 8 步 VΛ^8V^{-1}x0 = ", paste(round(xk_spec, 6), collapse = ", "))
say("与逐步迭代之差：", format(max(abs(xk_spec - xk)), scientific = TRUE))

# Fibonacci：[[F_{n+1}],[F_n]] = [[1,1],[1,0]]^n [[1],[0]]
say("\nFibonacci 用矩阵幂：")
Fm <- matrix(c(1, 1, 1, 0), 2, 2, byrow = TRUE)
fib <- function(n) as.vector(matpow_local(Fm, n) %*% c(1, 0))[2]
matpow_local <- function(M, k) {           # 快速幂
  R <- diag(nrow(M)); Mk <- M
  while (k > 0) { if (k %% 2 == 1) R <- R %*% Mk; Mk <- Mk %*% Mk; k <- k %/% 2 }
  R
}
say("  前 10 项：", paste(vapply(0:9, fib, numeric(1)), collapse = ", "))
# 通项由特征值给出：λ = (1±√5)/2，F_n = (λ1^n - λ2^n)/√5
phi <- (1 + sqrt(5)) / 2; psi <- (1 - sqrt(5)) / 2
say("  Binet 公式 F_10 = ", round((phi^10 - psi^10) / sqrt(5), 10),
    "，矩阵法 = ", fib(10))
say("  主特征值 φ = ", round(phi, 10), "，所以 F 的增长率趋于 φ")

# ===========================================================================
# 七、数值提醒：特征多项式为什么不能用来求特征值
# ===========================================================================
say("\n== 5.7 为什么不要展开特征多项式 ==")
# Wilkinson 多项式：根为 1, 2, ..., 10，系数稍微一动，根就大动
p <- poly_from_roots(1:10)                 # 升幂系数
comp <- function(pp) {                     # 由系数造友矩阵
  n <- length(pp) - 1
  A <- matrix(0, n, n)
  A[2:n, 1:(n - 1)] <- diag(n - 1)
  A[, n] <- -pp[1:n] / pp[n + 1]
  A
}
ev0 <- sort(Re(as.vector(eigen_geigen(eigenMat(comp(p)))$values())))
say("友矩阵的特征值：", paste(round(ev0, 6), collapse = ", "))
p2 <- p; p2[2] <- p2[2] * (1 + 1e-8)       # 只把 λ^9 的系数扰动 1e-8（相对）
ev1 <- sort(Re(as.vector(eigen_geigen(eigenMat(comp(p2)))$values())))
say("把其中一个系数扰动 1e-8 后：", paste(round(ev1, 6), collapse = ", "))
say("根的最大位移 = ", format(max(abs(ev1 - ev0)), scientific = TRUE),
    "（放大了约 ", round(max(abs(ev1 - ev0)) / 1e-8), " 倍）")
say("结论：先展开多项式再求根，会把误差放大几个数量级；")
say("      正确做法是直接对矩阵做 QR 迭代（eigen_geigen 内部就是 Schur 分解）。")

# ===========================================================================
# 八、实验习题
# ===========================================================================
# 1. 求 A = [[1,-1],[2,4]] 的特征值与特征向量，验证 Av = λv 与 tr、det 两条恒等式。
#
# 2. 讨论 a 取何值时 A = [[1,a],[0,2]] 可对角化。用 is_diagonalizable() 在
#    a = 0, 1, -3 上验证你的结论。
#
# 3. 证明：若 A 可逆，则 A 的特征值是 1/λ_i，特征向量不变。用代码验证，
#    并说明为什么「A 可逆」等价于「没有零特征值」。
#
# 4. 教材 §5.3 的例子中，λ=2 是二重根。请自己在该特征子空间里另取一组
#    正交单位基，写出对应的正交矩阵 Q，验证 Q'SQ = Λ 仍然成立，
#    并解释「正交矩阵不唯一」不影响对角化的结论。
#
# 5. 差分方程：设某城市人口按 x_{k+1} = [[0.9,0.1],[0.1,0.9]] x_k 在城区与
#    郊区之间迁移。求极限分布，并算出多少年后城区人口降到初始偏差的一半。
#    （提示：另一个特征值是 0.8，偏差按 0.8^k 衰减，解 0.8^k = 0.5。）
#
# 思考题
#   (1) 为什么实对称矩阵一定能正交对角化，而一般矩阵不一定？
#   (2) 相似不变量里为什么没有「元素之和」？举一个反例。
#   (3) 特征值对扰动的敏感程度由什么决定？（提示：左右特征向量的夹角。）

say("\nlab05 结束。")
