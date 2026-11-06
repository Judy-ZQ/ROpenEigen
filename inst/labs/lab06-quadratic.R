# ===========================================================================
# 实验六 · 二次型
#
# 对应教材：第六章
#   第一节 二次型及其标准形
#   第二节 化二次型为标准形（正交变换法、配方法）
#   第三节 正定二次型（顺序主子式、惯性定理）
#   第四节 二次型的应用（二次曲面、极值）
# 建议学时：3 学时
#
# 实验目的
#   1. 会在二次型与对称矩阵之间来回转换，并说明为什么必须对称化；
#   2. 用正交变换把二次型化成标准形，验证变换前后函数值不变；
#   3. 用配方法做同一个问题，比较两种标准形的异同（惯性定理）；
#   4. 掌握正定、负定、半定、不定的四种判别法，互相印证；
#   5. 把二次型用到二次曲面分类与多元函数极值上。
#
# 运行：Rscript lab06-quadratic.R   或在 RStudio 里 source()
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
# 一、二次型与对称矩阵（教材 §6.1）
# ===========================================================================
say("== 6.1 二次型的矩阵表示 ==")
# 教材 §6.3 例 1 的二次型：
#   f = 5x1^2 + x2^2 + 5x3^2 + 4x1x2 - 8x1x3 - 4x2x3
# 交叉项系数要平分给 a_ij 和 a_ji。
A <- quad_matrix(3, list(c(1, 1, 5), c(2, 2, 1), c(3, 3, 5),
                         c(1, 2, 4), c(1, 3, -8), c(2, 3, -4)))
say("二次型的矩阵 A =")
print(A)
say("A 对称：", relerr(A, t(A)) == 0)

# 代入验证：f(x) 用矩阵算和用代数式算必须一致
x <- c(1, -2, 3)
f_matrix <- quad_eval(A, x)
f_algebra <- 5 * x[1]^2 + x[2]^2 + 5 * x[3]^2 +
  4 * x[1] * x[2] - 8 * x[1] * x[3] - 4 * x[2] * x[3]
say("在 x = (1, -2, 3) 处：矩阵法 = ", f_matrix, "，代数式 = ", f_algebra)
expect_near(f_matrix, f_algebra, 1e-12, "二次型的矩阵表示不对")

# 二次型的矩阵只取对称的那一半：任何一个 A 都能换成 (A+A')/2
An <- matrix(c(5, 4, -8, 0, 1, -4, 0, 0, 5), 3, 3, byrow = TRUE)  # 上三角写法
say("上三角写法的矩阵乘 x 两次：", quad_eval(An, x))
say("换成 (A+A')/2 之后：", quad_eval((An + t(An)) / 2, x),
    "（与对称矩阵一致）")

# ===========================================================================
# 二、正交变换法化标准形（教材 §6.2）
# ===========================================================================
say("\n== 6.2 正交变换法 ==")
e <- eigen_eigen(eigenMat(A))
lam <- as.vector(e$values()); Q <- as.matrix(e$vectors())
say("特征值：", paste(round(lam, 10), collapse = ", "))
say("正交矩阵 Q（按列）：")
print(round(Q, 8))
say("Q'AQ = Λ：", relerr(t(Q) %*% A %*% Q, diag(lam)) < 1e-12)
say("标准形：f = ", paste(sprintf("%.6f y%d^2", lam, 1:3), collapse = " + "))

# 关键验证：作正交变换 x = Qy 后，函数值不变
say("\n变换前后函数值必须相同（随机取 5 组 y）：")
set.seed(4)
for (t in 1:5) {
  y <- round(rnorm(3), 2)
  x <- as.vector(Q %*% y)
  say(sprintf("  y = (%s) → f(x) = %.8f，Σλy^2 = %.8f",
              paste(sprintf("%.2f", y), collapse = ", "),
              quad_eval(A, x), sum(lam * y^2)))
}
v <- c(3, -1, 2)
say("正交变换只旋转、不改变长度：|Qv| = ", sqrt(sum((Q %*% v)^2)),
    "，|v| = ", sqrt(sum(v^2)))

# ===========================================================================
# 三、配方法（教材 §6.2）
# ===========================================================================
say("\n== 6.3 配方法 ==")
# 例：f = x1^2 + 2x1x2 + 2x2^2 + 2x2x3 + x3^2
# 手算：
#   f = (x1 + x2)^2 + (x2^2 + 2x2x3 + x3^2) = (x1+x2)^2 + (x2+x3)^2
# 令 y1 = x1+x2, y2 = x2+x3, y3 = x3
B <- quad_matrix(3, list(c(1, 1, 1), c(2, 2, 2), c(3, 3, 1),
                         c(1, 2, 2), c(2, 3, 2)))
say("B ="); print(B)
say("配方：f = (x1+x2)^2 + (x2+x3)^2")

# 由 y 反解 x，得到变换矩阵 C（x = C y）
#   y1 = x1+x2, y2 = x2+x3, y3 = x3
#   → x3 = y3, x2 = y2-y3, x1 = y1-y2+y3
C <- matrix(c(1, -1, 1,
              0, 1, -1,
              0, 0, 1), 3, 3, byrow = TRUE)
say("变换矩阵 C（x = C y）：")
print(C)
say("C'BC =")
print(round(t(C) %*% B %*% C, 10))
say("标准形：f = 1·y1^2 + 1·y2^2 + 0·y3^2（y3 缺失，因为 B 不满秩）")
expect_near(t(C) %*% B %*% C, diag(c(1, 1, 0)), 1e-12, "配方法的合同变换不对")

# 随机取 y 验证函数值不变（注意：这里的变换不是正交变换，不保长度）
say("\n配方法同样保持函数值（随机 3 组 y）：")
for (t in 1:3) {
  y <- round(rnorm(3), 2)
  x <- as.vector(C %*% y)
  say(sprintf("  y = (%s) → f(Cy) = %.8f，Σd y^2 = %.8f",
              paste(sprintf("%.2f", y), collapse = ", "),
              quad_eval(B, x), sum(c(1, 1, 0) * y^2)))
}

# ===========================================================================
# 四、惯性定理：标准形不唯一，惯性指数唯一
# ===========================================================================
say("\n== 6.4 惯性定理 ==")
# 对同一个二次型 B（上一节那个）分别用两种方法化标准形：
#   正交变换法 → 系数是特征值
#   配方法     → 系数是 1, 1, 0
lamB <- as.vector(eigen_eigen(eigenMat(B))$values())
say("对同一个 B 求两次标准形：")
say("  正交变换法：", paste(sprintf("%.10f", lamB), collapse = ", "))
say("  配方法    ：1, 1, 0")
say("  系数完全不同，但正、负、零的个数必须相同：")
say("    正交变换：正 ", sum(lamB > 1e-10), "，负 ", sum(lamB < -1e-10),
    "，零 ", sum(abs(lamB) <= 1e-10))
say("    配方法  ：正 2，负 0，零 1")
say("  两边的惯性指数都是 (2, 0, 1)，这正是惯性定理。")
say("  配方法给出 1、1、0，正交变换给出 0、1、3 —— 两者都是标准形，系数不同。")
say("  所以「标准形」不唯一；只有「规范形」（系数只取 ±1 和 0）才由惯性指数唯一确定。")

# 规范形：把标准形的非零系数都除以自己的绝对值
canon <- ifelse(abs(lamB) > 1e-10, sign(lamB), 0)
say("  正交变换法的规范形系数：", paste(canon, collapse = ", "),
    "；配方法的规范形系数：1, 1, 0")

# ===========================================================================
# 五、正定二次型的四种判别法（教材 §6.3）
# ===========================================================================
say("\n== 6.5 正定性判别 ==")

# 判别法一：顺序主子式全大于 0（教材 §6.3 定理 6.3.2）
say("对 A 算顺序主子式：", paste(round(leading_minors(A), 10), collapse = ", "))
say("  全为正 → A 正定")

# 判别法二：特征值全大于 0
say("对 A 算特征值：", paste(round(lam, 10), collapse = ", "))
say("  全为正 → A 正定")

# 判别法三：Cholesky 分解能成功（存在可逆下三角 L，A = L L'）
r <- try(eigen_chol(eigenMat(A)), silent = TRUE)
say("Cholesky 分解：", if (inherits(r, "try-error")) "失败 → 不正定" else "成功 → 正定")
if (!inherits(r, "try-error")) {
  L <- as.matrix(r)
  say("  L L' 还原 A 的误差：", format(relerr(L %*% t(L), A), scientific = TRUE))
}

# 判别法四：惯性指数 p = n（教材 §6.3 定理 6.3.3）
say("惯性指数：" , paste(names(inertia(A)), inertia(A), sep = "=", collapse = "，"))

# 把四种判别法打包成一个函数，对若干矩阵一起跑
classify_quad <- function(M, tol = 1e-10) {
  M <- (M + t(M)) / 2
  n <- nrow(M); inm <- inertia(M, tol)
  ch <- !inherits(try(eigen_chol(eigenMat(M)), silent = TRUE), "try-error")
  if (inm["positive"] == n) "正定"
  else if (inm["negative"] == n) "负定"
  else if (inm["zero"] > 0 && (inm["positive"] == 0 || inm["negative"] == 0)) "半定"
  else "不定"
}
say("\n四种情形的对照：")
tests <- list(
  list(name = "正定  [[2,1],[1,2]]", M = matrix(c(2, 1, 1, 2), 2, 2, byrow = TRUE)),
  list(name = "负定  [[-2,-1],[-1,-2]]", M = matrix(c(-2, -1, -1, -2), 2, 2, byrow = TRUE)),
  list(name = "不定  [[1,2],[2,1]]", M = matrix(c(1, 2, 2, 1), 2, 2, byrow = TRUE)),
  list(name = "半正定 [[1,1],[1,1]]", M = matrix(c(1, 1, 1, 1), 2, 2, byrow = TRUE))
)
for (tt in tests) {
  say(sprintf("  %-26s → %s（顺序主子式：%s）", tt$name, classify_quad(tt$M),
              paste(round(leading_minors(tt$M), 6), collapse = ", ")))
}

# 教材 §6.3 例 2：讨论 a 取何值时 f = a x^2 + a y^2 + z^2 + 2xy + 2yz 正定
say("\n教材 §6.3 例 2：f = a x^2 + a y^2 + z^2 + 2xy + 2yz")
quad_a <- function(a) matrix(c(a, 1, 0,
                               1, a, 1,
                               0, 1, 1), 3, 3, byrow = TRUE)
say("  顺序主子式：a11 = a，|A2| = a^2 - 1，|A| = a^2 - a - 1")
say("  （注意 |A| 是一次项为 0 的三阶行列式，手算得 a(a-1) - 1 = a^2 - a - 1）")
say("  三个都要 > 0：a > 0，a > 1，a^2 - a - 1 > 0")
root <- uniroot(function(a) a^2 - a - 1, c(1, 3))$root
say("  a^2 - a - 1 = 0 的正根 = (1+√5)/2 = ", round(root, 10),
    "（黄金比 (1+√5)/2 = ", round((1 + sqrt(5)) / 2, 10), "）")
say("  结论：a > (1+√5)/2 时正定；a = 1 时 |A2| = 0，需另作判断")
for (a in c(1.0, 1.4, 1.618, 1.7, 2.0)) {
  say(sprintf("    a = %.4f：主子式 = (%7.4f, %7.4f, %8.4f) → %s",
              a, a, a^2 - 1, a^2 - a - 1, classify_quad(quad_a(a))))
}

# ===========================================================================
# 六、应用一：二次曲面分类（教材 §6.4）
# ===========================================================================
say("\n== 6.6 应用：二次曲面 ==")
# 二次曲面 F(x,y,z) = [x y z] A [x y z]' = 1 的类型，由 A 的特征值符号决定：
#   三个正 → 椭球面；两正一负 → 单叶双曲面；一正两负 → 双叶双曲面
surface_type <- function(M) {
  nm <- inertia((M + t(M)) / 2)
  p <- nm["positive"]; q <- nm["negative"]
  if (p == 3) "椭球面" else if (p == 2 && q == 1) "单叶双曲面"
  else if (p == 1 && q == 2) "双叶双曲面" else "退化/其他"
}
surfaces <- list(
  list(name = "x^2 + 2y^2 + 3z^2 = 1", M = diag(c(1, 2, 3))),
  list(name = "x^2 + y^2 - z^2 = 1", M = diag(c(1, 1, -1))),
  list(name = "x^2 - y^2 - z^2 = 1", M = diag(c(1, -1, -1)))
)
for (sf in surfaces) {
  say(sprintf("  %-24s 特征值 (%s) → %s", sf$name,
              paste(round(as.vector(eigen(sf$M, symmetric = TRUE,
                                          only.values = TRUE)$values), 4),
                    collapse = ", "),
              surface_type(sf$M)))
}

# 带交叉项的曲面：先正交变换消去交叉项，再判断类型
Mx <- matrix(c(2, 1, 0,
               1, 2, 0,
               0, 0, 1), 3, 3, byrow = TRUE)
say("\n带交叉项的曲面 [x y z] M [x y z]' = 1，M =")
print(Mx)
say("  特征值：", paste(round(as.vector(eigen(Mx, symmetric = TRUE,
                                              only.values = TRUE)$values), 6), collapse = ", "))
say("  → 标准形 ", paste(sprintf("%.4f u%d^2", sort(as.vector(eigen(Mx, symmetric = TRUE,
                                                                     only.values = TRUE)$values)), 1:3),
                          collapse = " + "), " = 1")
say("  → ", surface_type(Mx))
say("  （旋转不改变曲面的类型，只改变坐标轴的方向）")

# ===========================================================================
# 七、应用二：多元函数的极值
# ===========================================================================
say("\n== 6.7 应用：Hessian 与极值 ==")
# f(x,y) = x^2 + 3xy + y^2 在原点：Hessian 就是二次型的矩阵
H <- matrix(c(2, 3, 3, 2), 2, 2, byrow = TRUE)
say("f = x^2 + 3xy + y^2 的 Hessian H =")
print(H)
say("  H 的特征值：", paste(round(as.vector(eigen(H, symmetric = TRUE,
                                                 only.values = TRUE)$values), 6), collapse = ", "))
say("  一正一负 → 鞍点（不定）")

# 比较一个极小值点
H2 <- matrix(c(2, 1, 1, 4), 2, 2, byrow = TRUE)
say("\ng = x^2 + xy + 2y^2 的 Hessian H =")
print(H2)
say("  特征值全正 → 严格极小值点")
say("  用 Cholesky 判断：", if (!inherits(try(eigen_chol(eigenMat(H2)), silent = TRUE),
                                          "try-error")) "正定（极小）" else "非正定")

# 数值印证：在原点附近随机撒点，看函数值是否都比 f(0) 大
f2 <- function(v) quad_eval(H2, v)
set.seed(21)
pts <- matrix(rnorm(2000), 1000, 2) * 0.3
vals <- apply(pts, 1, f2)
say("  在原点附近撒 1000 个点，f 的最小值 = ", format(min(vals), scientific = TRUE),
    "（≥ 0，说明原点确为极小）")
f1 <- function(v) quad_eval(H, v)
vals1 <- apply(pts, 1, f1)
say("  鞍点情形：最小值 = ", round(min(vals1), 6), "，最大值 = ", round(max(vals1), 6),
    "（有正有负）")

# ===========================================================================
# 八、实验习题
# ===========================================================================
# 1. 写出二次型 f = x1^2 + 2x2^2 + 3x3^2 - 4x1x2 - 6x2x3 的矩阵，
#    用正交变换化为标准形，写出所用的正交矩阵，并判断是否正定。
#
# 2. 用配方法把上题化成标准形，比较两种标准形的系数，验证惯性定理。
#
# 3. 教材 §6.3 定理 6.3.3：负定的充要条件是负惯性指数等于 n；
#    半正定的充要条件是正惯性指数等于秩。用 classify_quad() 各造一个例子验证。
#
# 4. 讨论 a, b 满足什么条件时，f = x1^2 + x2^2 + 5x3^2 + 2a x1x2 - 2x1x3 + 4x2x3
#    正定。（提示：写出矩阵，算三个顺序主子式，用 uniroot 找临界值。）
#
# 5. 二次曲面 x^2 + y^2 + z^2 + 2xy = 1 是什么曲面？
#    先正交变换，再判断，最后用随机点验证（在曲面上取点，检查是否满足方程）。
#
# 思考题
#   (1) 正交变换法与配方法都能化标准形，为什么工程上更常用正交变换？
#       （提示：正交变换的条件数是 1，配方法的变换矩阵可能很病态。）
#   (2) 二次型的标准形不唯一，但惯性指数唯一。这句话对「用正交变换得到的
#       标准形」还成立吗？请举例说明。
#   (3) 正定矩阵一定是满秩的，反过来对吗？半正定矩阵的秩与惯性指数什么关系？

say("\nlab06 结束。")
