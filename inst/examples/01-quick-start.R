# ===========================================================================
# 01 · 快速起步
#
# 对应 EigenExample/examples/ 与 snippets/ 里的入门示例：
#   examples/QuickStart_example.cpp
#   examples/QuickStart_example2_dynamic.cpp
#   examples/QuickStart_example2_fixed.cpp
#   examples/Tutorial_simple_example_dynamic_size.cpp
#   examples/Tutorial_simple_example_fixed_size.cpp
#   snippets/Tutorial_commainit_01.cpp
#   snippets/Tutorial_commainit_01b.cpp
#   snippets/Tutorial_commainit_02.cpp
#   snippets/Tutorial_AdvancedInitialization_{ThreeWays,CommaTemporary,Join,Block,
#                                            LinSpaced,Zero}.cpp
#
# 运行：Rscript 01-quick-start.R
#
# 两条转写规则，全书通用：
#   1. Eigen 的逗号初始化 `m << 1, 2, 3, 4;` 是**按行**填的，而 R 的 matrix()
#      默认按列。所以凡是要和 cpp 输出逐字对上的地方，一律写 byrow = TRUE。
#   2. Eigen 的对象一旦构造就有固定形状，R 的 eigenMat 也一样（构造时定形状），
#      但 R 里没有 Matrix3d / MatrixXd 的静态类型之分 —— 尺寸是运行期的值。
#      cpp 里 dynamic 与 fixed 两个版本，在 R 里是同一段代码，区别只在于有没有
#      把尺寸写成常量。
#
# 涉及随机数的例子，R 的 eigen_random() 走 R 的 RNG，Eigen 的 MatrixXd::Random()
# 走 std::rand，两者序列不同。所以这类例子里数字和 .cpp 的输出不会一样，
# 这里用 set.seed() 固定 R 侧的结果，便于重复。
# ===========================================================================

suppressPackageStartupMessages(library(ROpenEigen))

# ---- QuickStart_example.cpp -------------------------------------------------
# cpp:
#   MatrixXd m(2, 2);
#   m(0, 0) = 3;  m(1, 0) = 2.5;
#   m(0, 1) = -1; m(1, 1) = m(1, 0) + m(0, 1);
#   cout << m;
# 输出： 3  -1
#       2.5 1.5
#
# R 的 [i, j] 是 1-based，对应 cpp 的 m(i-1, j-1)。eigenMat 只在构造时拷贝一次，
# 之后每个元素都在 C++ 侧，逐元素赋值会经过 R <- C++ 的边界，所以更常见的写法
# 是先在 R 里把矩阵拼好（下面这样），再交给 eigenMat()。

m <- matrix(0, 2, 2)
m[1, 1] <- 3
m[2, 1] <- 2.5
m[1, 2] <- -1
m[2, 2] <- m[2, 1] + m[1, 2]
em <- eigenMat(m)
print(em)

# 也可以在 eigenMat 上原地改，语义与 cpp 的 m(i,j) = ... 一致：
em2 <- eigenMat(matrix(0, 2, 2))
em2[1, 1] <- 3
em2[2, 1] <- 2.5
em2[1, 2] <- -1
em2[2, 2] <- as.matrix(em2)[2, 1] + as.matrix(em2)[1, 2]
print(em2)

# ---- QuickStart_example2_dynamic.cpp ----------------------------------------
# cpp:
#   MatrixXd m = MatrixXd::Random(3, 3);
#   m = (m + MatrixXd::Constant(3, 3, 1.2)) * 50;
#   VectorXd v(3); v << 1, 2, 3;
#   cout << m * v;
#
# eigen_random() 产生 [-1, 1] 的均匀分布，与 Eigen 的 Random() 同一个区间。

set.seed(2026)
m <- (eigen_random(3, 3) + eigen_constant(1.2, 3, 3)) * 50
print(m)

v <- eigenVec(c(1, 2, 3))
print(m %*% v)

# ---- QuickStart_example2_fixed.cpp ------------------------------------------
# cpp 用 Matrix3d / Vector3d，尺寸写死在类型里。R 只有一种矩阵类型，
# 这段代码与上面完全相同 —— 这正是 R 侧不需要区分 fixed/dynamic 的原因。
# 如果确实想断言尺寸，在 R 里用 stopifnot(dim(m) == c(3, 3)) 比类型系统更直接。

set.seed(2026)
m3 <- (eigen_random(3, 3) + eigen_constant(1.2, 3, 3)) * 50
stopifnot(identical(dim(m3), c(3L, 3L)))
print(m3 %*% eigenVec(c(1, 2, 3)))

# ---- Tutorial_simple_example_dynamic_size.cpp -------------------------------
# cpp 在循环里构造 size x (size+1) 的整数矩阵，逐元素写 m(i, j) = i + j * size。
# R 侧用 outer() 一步得到同样的内容，再按需截取。

for (size in 1:4) {
  m <- outer(0:(size - 1), 0:size, function(i, j) i + j * size)
  print(eigenMat(m))
  cat("\n")
}

# ---- Tutorial_simple_example_fixed_size.cpp ---------------------------------
# cpp 用 Matrix3i，循环里只改数值不改形状。R 里就是复用同一个 eigenMat。

m <- eigenMat(matrix(0L, 3, 4))
for (j in 0:3) {
  for (i in 0:2) {
    m[i + 1, j + 1] <- i + 3 * j
  }
}
print(m)

# ---- Tutorial_commainit_01.cpp ----------------------------------------------
# cpp:  Matrix3f m;  m << 1,2,3, 4,5,6, 7,8,9;
# 逗号初始化按行填，等价于 R 的 byrow = TRUE。

m <- eigenMat(matrix(1:9, 3, 3, byrow = TRUE))
print(m)

# ---- Tutorial_commainit_01b.cpp ---------------------------------------------
# cpp 分三次填：第 0 行、block(1,0,2,2)、col(2).tail(2)。
# R 里对应 m[1, ]、m[2:3, 1:2]、m[2:3, 3]。

m <- eigenMat(matrix(0, 3, 3))
m[1, ] <- c(1, 2, 3)                 # m.row(0)
m[2:3, 1:2] <- matrix(c(4, 5, 7, 8), 2, byrow = TRUE)   # m.block(1,0,2,2)
m[2:3, 3] <- c(6, 9)                 # m.col(2).tail(2)
print(m)

# ---- Tutorial_commainit_02.cpp ----------------------------------------------
# cpp 用 3x3 字面量、零块、单位块拼出 5x5。
# R 侧用 rbind/cbind 拼装 —— eigen_rbind()/eigen_cbind() 返回 eigenMat，
# 但这里只是构造，直接拼 base R 矩阵更直观。

top_left  <- matrix(1:9, 3, 3, byrow = TRUE)
top_right <- matrix(0, 3, 2)
bot_left  <- matrix(0, 2, 3)
bot_right <- diag(2)
print(eigenMat(rbind(cbind(top_left, top_right), cbind(bot_left, bot_right))))

# ---- Tutorial_AdvancedInitialization_ThreeWays.cpp --------------------------
# 同一个分块对角矩阵，cpp 用三种写法：整块赋值、块上 setXxx()、逗号拼接。
# R 里三种写法都成立，而且都是拷贝语义，不存在别名问题。

size <- 6
half <- size / 2

# 写法一：整块赋值（右值是 base R 矩阵）
mat1 <- matrix(0, size, size)
mat1[1:half, 1:half]           <- matrix(0, half, half)
mat1[1:half, (half + 1):size]  <- diag(half)
mat1[(half + 1):size, 1:half]  <- diag(half)
mat1[(half + 1):size, (half + 1):size] <- matrix(0, half, half)
print(eigenMat(mat1))

# 写法二：整块赋值（右值是 eigenMat，直接由 eigen_xxx() 造出来）
mat2 <- eigenMat(matrix(0, size, size))
mat2[1:half, 1:half]           <- eigen_zeros(half, half)
mat2[1:half, (half + 1):size]  <- eigen_identity(half)
mat2[(half + 1):size, 1:half]  <- eigen_identity(half)
mat2[(half + 1):size, (half + 1):size] <- eigen_zeros(half, half)
print(mat2)

# 写法三：一次拼装，用 eigen_cbind()/eigen_rbind() 把四块粘起来
mat3 <- eigen_rbind(
  eigen_cbind(eigen_zeros(half, half), eigen_identity(half)),
  eigen_cbind(eigen_identity(half), eigen_zeros(half, half))
)
print(mat3)

# 三种写法结果一致：
stopifnot(identical(as.matrix(mat2), mat1), identical(as.matrix(mat3), mat1))

# ---- Tutorial_AdvancedInitialization_CommaTemporary.cpp ---------------------
# cpp 里 (MatrixXf(2,2) << 0,1,1,0).finished() 是一个**临时**矩阵，马上用于乘法。
# R 里临时矩阵就是普通对象，没有 finished() 这一步。

set.seed(1)
mat <- eigen_random(2, 3)
print(mat)
flip <- eigenMat(matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE))
print(flip)
print(flip %*% mat)                  # 交换两行

# ---- Tutorial_AdvancedInitialization_Join.cpp -------------------------------
# cpp 用逗号把两个行向量接成一个：joined << vec1, vec2。
# R 侧就是 c() 或 eigen_cbind()。

vec1 <- c(1, 2, 3)
vec2 <- c(1, 4, 9, 16)
print(eigenVec(c(vec1, vec2)))

# ---- Tutorial_AdvancedInitialization_Block.cpp ------------------------------
# cpp:  matB << matA, matA/10, matA/10, matA;   2x2 块拼成 4x4

matA <- matrix(1:4, 2, 2, byrow = TRUE)
matB <- rbind(cbind(matA, matA / 10), cbind(matA / 10, matA))
print(eigenMat(matB))

# ---- Tutorial_AdvancedInitialization_LinSpaced.cpp --------------------------
# cpp 用 ArrayXXf 建 10x4 表：角度、弧度、正弦、余弦。
# R 侧 eigen_linspace(from, to, n) 对应 ArrayXf::LinSpaced。

deg <- eigen_linspace(0, 90, 10)
rad <- deg * (pi / 180)
tab <- cbind(as.matrix(deg), as.matrix(rad),
             sin(as.matrix(rad)), cos(as.matrix(rad)))
colnames(tab) <- c("Degrees", "Radians", "Sine", "Cosine")
print(round(tab, 6))

# ---- Tutorial_AdvancedInitialization_Zero.cpp -------------------------------
# cpp:  ArrayXXf a = ArrayXXf::Zero(2, 3);
# R:    eigen_zeros(2, 3)

print(eigen_zeros(2, 3))
print(eigen_ones(2, 3))
print(eigen_identity(3))

cat("\n01-quick-start.R 结束\n")
