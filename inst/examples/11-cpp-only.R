# ===========================================================================
# 11 · 只属于 C++ 的写法：别名、表达式模板、模板元编程
#
# 对应 EigenExample：
#   snippets/TopicAliasing_{mult1,mult2,mult3,mult4,mult5,block,block_correct,
#                          cwise}.cpp
#   snippets/function_taking_ref.cpp / function_taking_eigenbase.cpp
#   snippets/TemplateKeyword_simple.cpp / TemplateKeyword_flexible.cpp
#   examples/CustomizingEigen_Inheritance.cpp
#   examples/nullary_indexing.cpp
#   examples/make_circulant.cpp（含 .entry/.evaluator/.expression/.main/.preamble/
#                                .traits 六个片段）/ make_circulant2.cpp
#   snippets/IOFormat.cpp
#   special_examples/random_cpp11.cpp
#
# 运行：Rscript 11-cpp-only.R
#
# 这一章的 .cpp 在 R 侧**没有对等物**，原因各不相同，逐条写在下面。
# 与其硬造一个"看起来像"的 R 版本，不如说清楚为什么不需要 —— 这些片段的
# 教训大多是"Eigen 的表达式模板带来的坑"，而 R 根本不用那套机制。
# ===========================================================================

suppressPackageStartupMessages(library(ROpenEigen))

# ---- TopicAliasing_mult1..5 / _block / _block_correct / _cwise --------------
# cpp 的主题是**别名**：m = m * m 会让 Eigen 在读取右值时把左值改坏，
# 除非用 .eval() 强制先算成临时量，或者用 .noalias() 声明不别名。
# 这在 R 里不会发生：每个运算符都返回新对象，赋值是重新绑定。
m <- eigenMat(matrix(c(2, 0, 0, 2), 2, 2, byrow = TRUE))
m <- m %*% m                       # TopicAliasing_mult3 的 noalias() 那种写法
print(m)                           # 4 4 / 4 4？不，是 2x2 对角 4

a <- eigenMat(matrix(1:9, 3, 3, byrow = TRUE))
a[1:2, 1:2] <- a[2:3, 2:3]         # TopicAliasing_block：右值先求值，安全
print(a)

# ---- function_taking_ref.cpp -----------------------------------------------
# cpp 用 Eigen::Ref<const MatrixXf> 接收"任意矩阵表达式"，避免模板膨胀。
# R 没有按引用传参这回事：函数参数拿到的是对象，eigenMat 内部的指针由
# 引用计数管理，传进函数不会复制数据。
inv_cond <- function(a) {
  s <- as.vector(as.matrix(eigen_svd(eigenMat(a))$singular_values()))
  s[length(s)] / s[1]
}
set.seed(9)
mm <- matrix(runif(16, -1, 1), 4, 4)
cat("inv_cond(m)          :", inv_cond(mm), "\n")
cat("inv_cond(m[1:3,1:3]) :", inv_cond(mm[1:3, 1:3]), "\n")
cat("inv_cond(m + I)      :", inv_cond(mm + diag(4)), "\n")

# ---- function_taking_eigenbase.cpp -----------------------------------------
# cpp 用 EigenBase<Derived> 接收"任何能报出尺寸的东西"，包括 asDiagonal() 这种
# 伪表达式。R 里没有伪表达式，尺寸永远是确定的：dim() 直接给。
print_size <- function(x) cat("size:", length(x), " (", nrow(x), ",", ncol(x), ")\n")
print_size(eigenVec(c(1, 2, 3)))
print_size(eigenMat(diag(c(1, 2, 3))))

# ---- TemplateKeyword_simple.cpp / TemplateKeyword_flexible.cpp -------------
# cpp 的 temple 关键字（dst.template triangularView<Upper>()）是 C++ 语法要求，
# 与线性代数无关。R 里"把上三角部分从 src 拷到 dst"就是一次带条件的赋值。
src <- matrix(1, 4, 4)
dst <- matrix(0, 4, 4)
up <- upper.tri(dst, diag = TRUE)
dst[up] <- src[up]
print(eigenMat(dst))

# ---- CustomizingEigen_Inheritance.cpp --------------------------------------
# cpp 从 Eigen::Matrix 派生一个类，加上自己的方法。R 的对应做法是给 eigenMat
# 写 S3 方法（或用 eigenMat 包一层），不需要继承：S3 派发本来就是开放的。
# 这里演示"用 eigenMat 承载的类"最简单的一种：
my_matrix <- function(x) structure(eigenMat(x), class = c("myMatrix", "eigenMat"))
print.myMatrix <- function(x, ...) { cat("myMatrix:", "\n"); print(as.matrix(x)) }
print(my_matrix(matrix(1:4, 2)))

# ---- nullary_indexing.cpp / random_cpp11.cpp -------------------------------
# cpp 用 NullaryExpr(rows, cols, functor) 按需生成每个元素（不先造矩阵再填）。
# R 里最接近的是 outer()/expand.grid + apply，或者直接向量化计算再 matrix()。
# random_cpp11.cpp 用 std::poisson_distribution 生成泊松随机数，R 里就是 rpois()。
set.seed(1)
v <- rpois(10, 4.1)
print(v)

# nullary_indexing 想表达的是"元素值由下标决定"：R 里 outer(行, 列, f) 正是这个。
f <- function(i, j) i * 10 + j
print(eigenMat(outer(0:2, 0:3, f)))

# ---- make_circulant.cpp / make_circulant2.cpp ------------------------------
# cpp 定义一个 circulant_functor，让 Eigen 惰性生成循环矩阵（每个元素取
# 向量第 (row - col) mod n 个）。R 里写成一个普通函数更直接：
make_circulant <- function(v) {
  n <- length(v)
  idx <- outer(seq_len(n), seq_len(n), function(i, j) ((i - j) %% n) + 1L)
  matrix(v[idx], n, n)
}
print(eigenMat(make_circulant(c(1, 2, 3))))

# ---- IOFormat.cpp ----------------------------------------------------------
# cpp 用 IOFormat 控制矩阵打印的分隔符、对齐、精度。R 里对应的是 print()
# 的参数（digits）与 format()/write.table() 的输出格式，没有"矩阵对象的
# 格式属性"这一层。
m <- matrix(c(1.111111, 2, 3.33333, 4, 5, 6, 7, 8.888888, 9), 3, 3, byrow = TRUE)
print(round(m, 3))                                  # CleanFmt 那种 3 位小数
cat(paste(sprintf("%.6f", as.vector(t(m))), collapse = ", "), "\n")  # OctaveFmt 风格

# ===========================================================================
# 尾声：demos/mix_eigen_and_c —— 这层 C 包装在 R 里是自动的
# ---- demos/mix_eigen_and_c/binary_library.cpp -----------------------------
# cpp 用一个 extern "C" 的不透明句柄把 MatrixXd 包成 C 可用的结构体：
#   MatrixXd_new(m,n) / MatrixXd_set_coeff / MatrixXd_multiply / MatrixXd_delete
# 内部靠 reinterpret_cast 在 C_MatrixXd* 与 MatrixXd* 之间来回转。
# ROpenEigen 的 src/ 里有一模一样的一层：eigenMat 就是那个不透明句柄，
# 只是"句柄"换成了 R 的 EXTPTRSXP，"delete"换成了外部指针的终结器（GC 触发），
# 所以用户在 R 里看不到 new/delete，也不会忘记 delete。
m <- eigenMat(matrix(c(1, 2, 3, 4), 2, 2))
cat("句柄类型:", typeof(m), " class:", paste(class(m), collapse = "/"), "\n")
cat("它是不可直接读写的裸指针，取内容必须 as.matrix()：\n")
print(as.matrix(m))
print(m %*% m)                     # MatrixXd_multiply(matrix1, matrix1, matrix2)

# ---- demo_Map_MatrixXd 的原地改写 ----------------------------------------
# cpp 的 Map_MatrixXd_add(map, map, map) 会**改写调用方的 double 数组**
# （数组元素当场翻倍）。R 的 Map 层不存在，且 R 是写时复制：
# `x2 <- x * 2` 一定得到新对象，x 本身不变。
x <- 0:4
x2 <- x * 2
cat("原数组未被改写:", paste(x, collapse = " "), "\n")
cat("新对象:", paste(x2, collapse = " "), "\n")

cat("\n11-cpp-only.R 结束\n")
