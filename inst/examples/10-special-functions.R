# ===========================================================================
# 10 · 特殊函数（unsupported/SpecialFunctions）
#
# 对应 EigenExample：
#   examples/Cwise_erf.cpp
#   examples/Cwise_erfc.cpp
#   examples/Cwise_lgamma.cpp
#
# 运行：Rscript 10-special-functions.R
#
# EigenExample 里只有这三个片段用到 SpecialFunctions 模块，但它们代表的是
# 整个 unsupported/Eigen/SpecialFunctions：同一个头文件里的 igamma、zeta、
# betainc、bessel 系列在 ROpenEigen 里也一并封装了，本章一并给出，并在每条
# 上标出它在 Eigen 里的名字。
#
# 与 03 章那一组 Cwise 逐元素函数的区别：那些是 Dense 模块自带的，
# 这些来自 unsupported 模块，必须 #include <unsupported/Eigen/SpecialFunctions>。
# ===========================================================================

suppressPackageStartupMessages(library(ROpenEigen))

v <- c(-0.5, 2, 0, -7)

# ---- Cwise_erf.cpp ----------------------------------------------------------
# cpp:  Array4d v(-0.5, 2, 0, -7);  v.erf()   ->  -0.5205 0.995322 0 -1
print(eigen_erf(v))

# ---- Cwise_erfc.cpp ---------------------------------------------------------
# cpp:  v.erfc()  ->  1.5205 0.00467773 1 2
print(eigen_erfc(v))

# ---- Cwise_lgamma.cpp -------------------------------------------------------
# cpp:  Array4d v(0.5, 10, 0, -1);  v.lgamma()  ->  0.572365 12.8018 inf inf
# R 的 lgamma() 对 0 和负整数同样返回 Inf，两者一致。
print(eigen_gammaln(c(0.5, 10, 0, -1)))
print(lgamma(c(0.5, 10, 0, -1)))          # base R 对照

# ---- 其余同模块函数（cpp 未单列，同一个头文件）-----------------------------
# 每条都给出 R 的对照，方便确认定义一致。eigen_* 返回的始终是 eigenMat，
# 要数字就 as.vector()。
cat("digamma(1)      :", as.vector(eigen_digamma(1)),
    " base:", digamma(1), "\n")
cat("ndtri(0.975)    :", as.vector(eigen_ndtri(0.975)),
    " base qnorm:", qnorm(0.975), "\n")
cat("igamma(2, 3)    :", as.vector(eigen_igamma(2, 3)),
    " base pgamma:", pgamma(3, 2), "\n")
cat("igammac(2, 3)   :", as.vector(eigen_igammac(2, 3)),
    " base 上尾:", pgamma(3, 2, lower.tail = FALSE), "\n")
cat("betainc(2,3,0.4):", as.vector(eigen_betainc(2, 3, 0.4)),
    " base pbeta:", pbeta(0.4, 2, 3), "\n")
cat("polygamma(1, 2) :", as.vector(eigen_polygamma(1, 2)),
    " base psigamma:", psigamma(2, 1), "\n")

# zeta：Eigen 实现的是 Hurwitz 级数，不做解析延拓 —— x <= 1 直接给 NaN，
# 而 R 的 zeta() 在 0 处会返回 -1/2。这个差异被测试钉住了。
# base R 没有 zeta()，这里的数值只能与数学值对照：zeta(2) = pi^2/6 ≈ 1.644934。
cat("zeta(2)         :", as.vector(eigen_zeta(2)),
    " pi^2/6 =", pi^2 / 6, "\n")
cat("zeta(0)         :", as.vector(eigen_zeta(0)), "（不做解析延拓，故为 NaN）\n")

# Bessel 家族：type 参数就是 Eigen 的十二个函数名，末尾的 e 表示指数缩放版。
cat("besselJ(0, 1)   :", as.vector(eigen_bessel(1, "j0")),
    " base:", besselJ(1, 0), "\n")
cat("besselI(0, 1)   :", as.vector(eigen_bessel(1, "i0")),
    " base:", besselI(1, 0), "\n")
cat("besselK(1, 2)   :", as.vector(eigen_bessel(2, "k1")),
    " base:", besselK(2, 1), "\n")

# 逐元素语义：给矩阵就按元素算，形状不变。
m <- matrix(c(0.5, 1, 2, 3), 2, 2, byrow = TRUE)
print(eigen_erf(m))
print(eigen_gammaln(m))

cat("\n10-special-functions.R 结束\n")
