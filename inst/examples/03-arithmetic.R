# ===========================================================================
# 03 · 算术运算与逐元素函数
#
# 对应 EigenExample：
#   examples/tut_arithmetic_add_sub.cpp
#   examples/tut_arithmetic_scalar_mul_div.cpp
#   examples/tut_arithmetic_matrix_mul.cpp
#   examples/tut_arithmetic_dot_cross.cpp
#   examples/tut_arithmetic_redux_basic.cpp
#   examples/Tutorial_ArrayClass_{addition,mult,cwise_other,interop,interop_matrix,
#                                 accessors}.cpp
#   snippets/tut_arithmetic_transpose_{aliasing,conjugate,inplace}.cpp
#   snippets/MatrixBase_{array,array_const}.cpp
#   snippets/Cwise_*.cpp                                        （52 个）
#   snippets/MatrixBase_cwise{Abs,Abs2,Arg,Equal,Inverse,Max,Min,
#                             NotEqual,Product,Quotient,Sign,Sqrt}.cpp
#      —— 与 Cwise_* 是同一条 API 的成员函数写法，R 侧写法完全相同，
#         所以只在下面的 Cwise 小节里各给一次例子，不再重复。
#
# 运行：Rscript 03-arithmetic.R
#
# 本章的核心对照：Eigen 有 Matrix 与 Array 两套语义，Matrix 的 * 是矩阵乘，
# Array 的 * 是逐元素乘；R 只有一套 —— * 永远逐元素，矩阵乘写 %*%。
# 所以 Eigen 的 Array 片段在 R 里全都能用最朴素的写法表达。
# ===========================================================================

suppressPackageStartupMessages(library(ROpenEigen))

# ===========================================================================
# 一、基本算术
# ===========================================================================

# ---- tut_arithmetic_add_sub.cpp ---------------------------------------------
# cpp 里 a(0,0)=1, a(1,0)=2, a(0,1)=3, a(1,1)=4，按行写就是 1 3 / 2 4。
a <- eigenMat(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
b <- eigenMat(matrix(c(2, 3, 1, 4), 2, 2, byrow = TRUE))

print(a + b)
print(a - b)
print(a + b)                                   # cpp 的 a += b 之后 a 就是 a+b
print(-eigenVec(c(1, 2, 3)) + eigenVec(c(1, 0, 0)) - eigenVec(c(1, 2, 3)))

# ---- tut_arithmetic_scalar_mul_div.cpp --------------------------------------
# cpp:  a * 2.5 / ...  R 的标量乘除作用于 eigenMat 的每个系数。
a <- eigenMat(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
print(a * 2.5)
print(a / 2)
print(eigenVec(c(1, 2, 3)) * 0.1)
print(eigenVec(c(1, 2, 3)) * 2)

# ---- tut_arithmetic_matrix_mul.cpp ------------------------------------------
# 矩阵乘在 R 里是 %*%，并且 eigenMat 与 numeric 向量可以混着乘。
mat <- eigenMat(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
u <- c(-1, 1)
v <- c(2, 0)

print(mat %*% mat)
print(mat %*% u)                    # 矩阵 × 列向量
print(t(u) %*% mat)                 # 行向量 × 矩阵
print(as.numeric(t(u) %*% v))       # 内积
print(eigen_outer(eigenMat(u), eigenMat(v)))   # 外积 u v^T

# ---- tut_arithmetic_dot_cross.cpp -------------------------------------------
# 点积用 eigen_dot()，叉积用 eigen_cross()。注意 eigen_dot() **不取共轭**，
# 与 Eigen 的 MatrixBase::dot() 相反（见 README「Dot product」一节）。
p <- eigenVec(c(1, 2, 3)); q <- eigenVec(c(0, 1, 2))
print(eigen_dot(p, q))
print(eigen_cross(p, q))

# ---- tut_arithmetic_redux_basic.cpp -----------------------------------------
# 归约：sum/prod/mean/min/max/trace。sum()、mean() 等直接派发；min/max 在
# base R 里语义不同（有 na.rm、返回标量），所以包提供的是 eigen_which_min() 与
# eigen_sum() 一族，外加 eigen_trace()。
mat <- eigenMat(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
cat("sum:  ", eigen_sum(mat), "\n")
cat("prod: ", eigen_prod(mat), "\n")
cat("mean: ", eigen_mean(mat), "\n")
cat("min:  ", min(as.vector(mat)), "\n")
cat("max:  ", max(as.vector(mat)), "\n")
cat("trace:", eigen_trace(mat), "\n")

# ---- tut_arithmetic_transpose_aliasing.cpp ----------------------------------
# cpp 里 a = a.transpose() 会踩别名，触发 Eigen 断言（本目录跑这个 cpp 时
# 直接 abort）。R 里 a <- t(a) 永远安全：t() 生成新对象，赋值是重新绑定。
a <- eigenMat(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
a <- t(a)
print(a)

# ---- tut_arithmetic_transpose_conjugate.cpp ---------------------------------
# 转置 t()、共轭 Conj()、共轭转置 eigen_adjoint()。
a <- eigenMat(matrix(c(1 + 2i, 3 - 1i, 0 + 1i, 2 + 0i), 2, 2, byrow = TRUE))
print(a)
print(t(a))
print(Conj(a))
print(eigen_adjoint(a))

# ---- tut_arithmetic_transpose_inplace.cpp -----------------------------------
# cpp 的 transposeInPlace() 是原地转置。R 里 t() 已经生成新对象，
# 不需要（也没有）原地版本。
a <- eigenMat(matrix(1:6, 2, 3, byrow = TRUE))
print(t(a))

# ---- MatrixBase_transpose.cpp / MatrixBase_adjoint.cpp ----------------------
# cpp 的 transpose() 与 adjoint()（共轭转置）。R 里分别是 t() 与 eigen_adjoint()，
# 实数矩阵上两者相同，复数矩阵上差一个共轭。
m <- eigenMat(matrix(1:6, 2, 3, byrow = TRUE))
print(t(m))
z <- eigenVec(c(1 + 2i, 3 - 1i))
print(t(z))            # 只转置
print(eigen_adjoint(z)) # 共轭转置

# ===========================================================================
# 二、Matrix 与 Array 的互转
# ===========================================================================

# ---- Tutorial_ArrayClass_addition.cpp ---------------------------------------
# cpp 用 Array 才能 a + b 逐元素；R 的 + 本来就是逐元素。
a <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
b <- matrix(c(5, 6, 7, 8), 2, 2, byrow = TRUE)
print(eigenMat(a) + eigenMat(b))

# ---- Tutorial_ArrayClass_mult.cpp -------------------------------------------
# cpp:  ArrayXXf a(2,2), b(2,2);  a * b  是逐元素乘。
print(eigenMat(a) * eigenMat(b))

# ---- Tutorial_ArrayClass_cwise_other.cpp ------------------------------------
# 逐元素的其他运算：sqrt、比较、条件。
c <- eigenMat(matrix(c(1, 4, 9, 16), 2, 2, byrow = TRUE))
print(sqrt(c))
print(c > 4)

# ---- Tutorial_ArrayClass_interop.cpp ----------------------------------------
# cpp 演示 Matrix 与 Array 之间的 .array() / .matrix() 互转。R 里两者没有
# 类型之分：同一个 eigenMat，* 是逐元素，%*% 是矩阵乘。
m <- eigenMat(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
print(m * m)      # Array 语义
print(m %*% m)    # Matrix 语义

# ---- Tutorial_ArrayClass_interop_matrix.cpp ---------------------------------
# cpp:  m.array() * n.array() 之后再 .matrix()。R 里一步到位。
n <- eigenMat(matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE))
print(m * n)

# ---- Tutorial_ArrayClass_accessors.cpp --------------------------------------
# cpp 用 operator() 逐元素读写 Array。R 里是 [i]（线性下标）或 [i, j]。
arr <- eigenMat(matrix(0, 2, 3))
for (k in 1:6) arr[k] <- k * 10
print(arr)

# ---- MatrixBase_array.cpp / MatrixBase_array_const.cpp ----------------------
# cpp 的 .array() 得到一个 Array 视图，.array().abs() 之类才有逐元素语义。
# R 没有这层包装：eigenMat 上的 abs() 已经逐元素。
print(abs(eigenMat(matrix(c(-1, 2, -3, 4), 2, 2, byrow = TRUE))))

# ===========================================================================
# 三、Cwise_* 逐元素函数
#
# Eigen 的 cwiseXxx() 全部是逐元素，R 里对应 base R 的同名向量化函数，
# 直接作用在 eigenMat 上即可（Ops/Math 两组都派发）。下面按 cpp 文件逐条列。
# 每条的 "cpp:" 行是原始调用，紧跟的 R 行是等价写法。
# ===========================================================================

v <- eigenVec(c(1, -2, -3))
w <- eigenVec(c(4, 2, 3))

# ---- Cwise_abs.cpp ----------------------------------------------------------
# cpp:  v.abs()            R: abs(v)
print(abs(v))

# ---- Cwise_abs2.cpp ---------------------------------------------------------
# cpp:  v.abs2()           R: abs(v)^2  （等价于 v*v，但复数下是 |z|^2）
print(abs(v)^2)

# ---- Cwise_square.cpp -------------------------------------------------------
# cpp:  v.square()         R: v^2
print(eigenVec(c(2, 3, 4))^2)

# ---- Cwise_cube.cpp ---------------------------------------------------------
# cpp:  v.cube()           R: v^3
print(eigenVec(c(2, 3, 4))^3)

# ---- Cwise_sqrt.cpp ---------------------------------------------------------
# cpp:  v.sqrt()           R: sqrt(v)
print(sqrt(eigenVec(c(1, 2, 4))))

# ---- Cwise_cbrt.cpp ---------------------------------------------------------
# cpp:  v.cbrt()           R: sign(v) * abs(v)^(1/3)  （负数也要能开）
print(sign(eigenVec(c(1, 2, 4))) * abs(eigenVec(c(1, 2, 4)))^(1 / 3))

# ---- Cwise_inverse.cpp ------------------------------------------------------
# cpp:  v.inverse()        R: 1/v
print(1 / eigenVec(c(2, 3, 4)))

# ---- Cwise_sign.cpp ---------------------------------------------------------
# cpp:  v.sign()           R: sign(v)
print(sign(eigenVec(c(-3, 5, 0))))

# ---- Cwise_exp.cpp / Cwise_log.cpp / Cwise_log10.cpp ------------------------
# cpp:  v.exp(), v.log(), log10(v)
print(exp(eigenVec(c(1, 2, 3))))
print(log(eigenVec(c(1, 2, 3))))
print(log10(eigenVec(c(-1, 0, 1, 2))))     # 负数的 log10 是 NaN，0 是 -Inf

# ---- Cwise_pow.cpp / Cwise_scalar_power_array.cpp / Cwise_array_power_array.cpp
# cpp:  v.pow(0.333333), pow(10, e), x^e
print(eigenVec(c(8, 27, 64))^0.333333)
print(10^eigenVec(c(2, -3, 1 / 3)))
print(eigenVec(c(8, 25, 3))^eigenVec(c(1 / 3, 0.5, 2)))

# ---- Cwise_sin.cpp / Cwise_cos.cpp / Cwise_tan.cpp --------------------------
# cpp:  v.sin(), v.cos(), v.tan()      参数是 π, π/2, π/3
v <- eigenVec(c(pi, pi / 2, pi / 3))
print(sin(v)); print(cos(v)); print(tan(v))

# ---- Cwise_asin.cpp / Cwise_acos.cpp / Cwise_atan.cpp -----------------------
print(asin(eigenVec(c(0, sqrt(2) / 2, 1))))
print(acos(eigenVec(c(0, sqrt(2) / 2, 1))))
print(atan(eigen_linspace(0, 1, 5)))

# ---- Cwise_sinh.cpp / Cwise_cosh.cpp / Cwise_tanh.cpp -----------------------
v <- eigen_linspace(0, 1, 5)
print(sinh(v)); print(cosh(v)); print(tanh(v))

# ---- Cwise_array_atan2_array.cpp --------------------------------------------
# cpp:  y.atan2(x)，即逐元素 atan2(y, x)。
# atan2 在 R 里是 primitive，不走 S3 派发，所以不能直接喂 eigenMat；
# 转成普通向量再算，这也是唯一一处必须"落地"再算的逐元素运算。
y <- c(8, -25, 3); x <- c(1 / 3, 0.5, -2)
print(atan2(y, x))

# ---- Cwise_ceil.cpp / Cwise_floor.cpp / Cwise_rint.cpp / Cwise_round.cpp / Cwise_trunc.cpp
v <- eigen_linspace(-2, 2, 7)
print(ceiling(v))
print(floor(v))
print(trunc(v))
# R 的 round() 是 IEEE 的 round-half-to-even，与 Eigen 的 rint() 一致；
# Eigen 的 round() 是 half-away-from-zero，R 里没有直接对应（要自己写）：
print(round(v))
cat("Eigen round(0.5)=1, rint(0.5)=0；R round(0.5)=", round(0.5), "\n")

# ---- Cwise_plus.cpp / Cwise_minus.cpp ---------------------------------------
print(eigenVec(c(1, 2, 3)) + 5)
print(eigenVec(c(1, 2, 3)) - 5)

# ---- Cwise_plus_equal.cpp / Cwise_minus_equal.cpp / Cwise_times_equal.cpp / Cwise_slash_equal.cpp
# cpp 的 v += 5 之类在 R 里就是重新赋值（R 没有原地复合赋值）。
v <- eigenVec(c(1, 2, 3)); print(v + 5)
v <- eigenVec(c(1, 2, 3)); print(v - 5)
print(eigenVec(c(1, 2, 3)) * eigenVec(c(2, 3, 0)))
print(eigenVec(c(3, 2, 4)) / eigenVec(c(5, 4, 2)))

# ---- Cwise_product.cpp / Cwise_quotient.cpp ---------------------------------
# cpp 的 a * b / a / b 在 Array 上就是逐元素，R 里同样。
print(eigenVec(c(1, 2, 3)) * eigenVec(c(4, 5, 6)))
print(eigenVec(c(2, 3, 4)) / eigenVec(c(4, 2, 3)))

# ---- Cwise_min.cpp / Cwise_max.cpp ------------------------------------------
# cpp:  v.min(w), v.max(w)  -> 逐元素
print(pmin(eigenVec(c(2, 3, 4)), eigenVec(c(4, 2, 3))))
print(pmax(eigenVec(c(2, 3, 4)), eigenVec(c(4, 2, 3))))

# ---- Cwise_equal_equal.cpp / Cwise_not_equal.cpp ---------------------------
print(eigenVec(c(1, 2, 3)) == eigenVec(c(3, 2, 1)))
print(eigenVec(c(1, 2, 3)) != eigenVec(c(3, 2, 1)))

# ---- Cwise_less.cpp / Cwise_less_equal.cpp / Cwise_greater.cpp / Cwise_greater_equal.cpp
print(eigenVec(c(1, 2, 3)) <  eigenVec(c(3, 2, 1)))
print(eigenVec(c(1, 2, 3)) <= eigenVec(c(3, 2, 1)))
print(eigenVec(c(1, 2, 3)) >  eigenVec(c(3, 2, 1)))
print(eigenVec(c(1, 2, 3)) >= eigenVec(c(3, 2, 1)))

# ---- Cwise_boolean_and.cpp / Cwise_boolean_or.cpp / Cwise_boolean_not.cpp ---
v <- eigenVec(c(-1, 2, 1)); w <- eigenVec(c(-3, 2, 3))
print((v < w) & (v < 0))
print((v < w) | (v < 0))
v <- eigenVec(c(1, NaN, Inf)); print(!is.finite(v))

# ---- Cwise_isFinite.cpp / Cwise_isInf.cpp / Cwise_isNaN.cpp ----------------
v <- eigenVec(c(1, NaN, Inf))
print(is.finite(v))
print(is.infinite(v))
print(is.na(v))

# ---- Cwise_arg.cpp ----------------------------------------------------------
# cpp:  arg(v)             R: Arg(v)
print(Arg(eigenVec(c(0.0327601 - 0.726553i, 0.362112 - 0.15168i))))

cat("\n03-arithmetic.R 结束\n")
