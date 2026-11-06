# ===========================================================================
# 02 · 矩阵与向量：构造、尺寸、存储顺序
#
# 对应 EigenExample：
#   snippets/Matrix_initializer_list_23_cxx11.cpp
#   snippets/Matrix_initializer_list_vector_cxx11.cpp
#   snippets/Matrix_variadic_ctor_cxx11.cpp
#   snippets/Array_initializer_list_23_cxx11.cpp
#   snippets/Array_initializer_list_vector_cxx11.cpp
#   snippets/Array_variadic_ctor_cxx11.cpp
#   snippets/Matrix_resize_int.cpp
#   snippets/Matrix_resize_int_int.cpp
#   snippets/Matrix_resize_int_NoChange.cpp
#   snippets/Matrix_resize_NoChange_int.cpp
#   examples/tut_matrix_resize.cpp
#   examples/tut_matrix_resize_fixed_size.cpp
#   examples/tut_matrix_assignment_resizing.cpp
#   snippets/Matrix_set{Constant,Identity,Ones,Random,Zero}_int*.cpp
#   snippets/MatrixBase_{set,setZero,setOnes,setIdentity,setRandom}.cpp
#   snippets/MatrixBase_{ones,zero,identity,random}.cpp
#   snippets/MatrixBase_{ones,zero,identity,random}_int*.cpp
#   snippets/Map_{simple,general_stride,inner_stride,outer_stride,placement_new}.cpp
#   snippets/Matrix_Map_stride.cpp
#   snippets/Tutorial_Map_{rowmajor,using}.cpp
#   snippets/TopicStorageOrders_example.cpp
#   snippets/MatrixBase_reshaped_{auto,fixed,int_int,to_vector}.cpp
#   snippets/Tutorial_reshaped_vs_resize_{1,2}.cpp
#   snippets/Tutorial_ReshapeMat2Mat.cpp / Tutorial_ReshapeMat2Vec.cpp
#   snippets/MatrixBase_cast.cpp
#   snippets/Tutorial_range_for_loop_{1d,2d}_cxx11.cpp
#   snippets/Tutorial_std_sort.cpp / Tutorial_std_sort_rows_cxx11.cpp
#
# 运行：Rscript 02-matrix-class.R
# ===========================================================================

suppressPackageStartupMessages(library(ROpenEigen))

# ---- Matrix_initializer_list_23_cxx11.cpp -----------------------------------
# cpp:  MatrixXd m{{1, 2, 3}, {4, 5, 6}};      // 2x3，按行写
m <- eigenMat(matrix(c(1, 2, 3, 4, 5, 6), 2, 3, byrow = TRUE))
print(m)

# ---- Matrix_initializer_list_vector_cxx11.cpp -------------------------------
# cpp:  VectorXi v{{1, 2}};                    // 列向量
print(eigenVec(c(1L, 2L)))

# ---- Matrix_variadic_ctor_cxx11.cpp -----------------------------------------
# cpp:  Matrix<int, 1, 6> a(1,2,3,4,5,6);  Matrix<int, 3, 1> b{1,2,3};
# 参数个数固定形状的构造在 R 里就是给 c() 一个向量，形状由 nrow/ncol 决定。
print(eigenMat(matrix(c(1, 2, 3, 4, 5, 6), 1, 6)))   # 行向量
print(eigenMat(matrix(c(1, 2, 3), 3, 1)))            # 列向量

# ---- Array_*_cxx11.cpp ------------------------------------------------------
# Array 与 Matrix 的区别只在默认运算（系数运算 vs 线性代数）。R 里所有矩阵
# 默认就是逐元素的 +、*、/，所以 Array 的三种写法与上面 Matrix 的完全一样。
print(eigenMat(matrix(c(1, 2, 3, 3, 4, 5), 2, 3, byrow = TRUE)))
print(eigenVec(c(1L, 2L, 3L, 4L, 5L)))

# ---- Matrix_resize_int.cpp --------------------------------------------------
# cpp:  VectorXd v(10);  v.resize(3);       -> 3x1
#       RowVector3d w;   w.resize(3);      -> 仍然是 1x3（定长类型，resize 无效）
# Eigen 的 resize 可以改元素个数，R 的 dim<- 不行：它只重排、不增删，元素数
# 必须对得上。要改元素个数，R 里就是取子集（或重新构造）。
# 这个 cpp 例子只打印尺寸，所以两种写法的结果一样。

v <- eigenVec(seq_len(10))
v <- eigenVec(v[seq_len(3)])          # 尺寸变小 = 取前三个
print(dim(v))

w <- eigenMat(matrix(0, 1, 3))
dim(w) <- c(3L, 1L)          # 元素数不变，形状可以随便改；cpp 里的 RowVector3d 不会
print(dim(w))

# ---- Matrix_resize_int_int.cpp ----------------------------------------------
# cpp 把 2x3 改成 3x2，注释说这是 "conservative resizing"，系数不变，只是重新
# 解释形状；改成 2x2 则**不是** conservative，cpp 里内容变成未初始化。
# R 的 dim<- 走的是前一种语义；后一种 R 直接拒绝，不会给你一块未初始化的内存。

m <- eigenMat(matrix(1:6, 2, 3, byrow = TRUE))
print(m)
dim(m) <- c(3L, 2L)          # 保守重排：列主序的 1 4 2 5 3 6 摊成 3x2
print(m)

# 长度不匹配时 R 会报错，这一点比 cpp 安全：
print(tryCatch({ dim(m) <- c(2L, 2L); "改掉了" },
               error = function(e) paste("拒绝:", conditionMessage(e))))

# ---- Matrix_resize_int_NoChange.cpp / Matrix_resize_NoChange_int.cpp --------
# cpp 的 NoChange 常量表示"这一维别动"。R 的矩阵没有这种半指定尺寸的接口：
# 元素数一变就得新建，于是"另一维保持不变"变成直接写出那个值。
print(dim(eigen_zeros(5, 4)))   # resize(5, NoChange)，原矩阵 3x4
print(dim(eigen_zeros(3, 5)))   # resize(NoChange, 5)

# ---- tut_matrix_resize.cpp --------------------------------------------------
# cpp 打印 rows/cols/size。R 里分别是 dim()[1]、dim()[2]、length()。
m <- eigen_zeros(4, 3)          # cpp 里是 2x5 的 m.resize(4,3)
cat("size:", dim(m)[1], "x", dim(m)[2], " coefficients:", length(m), "\n")
v <- eigenVec(numeric(5))       # cpp 里是 2 维的 v.resize(5)
cat("v size:", length(v), " as matrix:", dim(v)[1], "x", dim(v)[2], "\n")

# ---- tut_matrix_resize_fixed_size.cpp ---------------------------------------
# cpp 里 Matrix4d 的 resize(4,4) 是空操作。R 里没有定长类型，dim<- 到同样
# 尺寸也只是重排（内容保持不变），效果相同。
m <- eigenMat(matrix(1:16, 4, 4))
dim(m) <- c(4L, 4L)
cat("m is", dim(m)[1], "x", dim(m)[2], "\n")

# ---- tut_matrix_assignment_resizing.cpp -------------------------------------
# cpp 里 a = b 会让 a 跟着 b 变尺寸（Eigen 的赋值是深拷贝 + 自动 resize）。
# R 的赋值是**绑定**：a <- b 之后 a 就是 b，连拷贝都没有，直到有人修改。
a <- eigenMat(matrix(0, 2, 2))
cat("a is", dim(a)[1], "x", dim(a)[2], "\n")
b <- eigenMat(matrix(0, 3, 3))
a <- b
cat("a is now", dim(a)[1], "x", dim(a)[2], "\n")

# ---- Matrix_setXxx_int*.cpp -------------------------------------------------
# cpp 的 MatrixXd::Zero(r,c)、::Ones、::Identity、::Constant、::Random 都是
# 静态构造。R 侧是同名的一族函数，参数顺序统一为 (值, 行, 列)。

print(eigen_zeros(2, 3))
print(eigen_ones(1, 4))
print(eigen_constant(2.5, 2, 2))
print(eigen_identity(3))
set.seed(7); print(round(eigen_random(2, 2), 4))

# ---- MatrixBase_set.cpp -----------------------------------------------------
# cpp 依次用逗号初始化、块初始化、再整块重写。最后一步在 R 里就是逐块赋值。

m1 <- matrix(1:9, 3, 3, byrow = TRUE)
m2 <- diag(3)
m2[1:2, 1:2] <- matrix(c(10, 11, 12, 13), 2, byrow = TRUE)
print(eigenMat(m2))

v1 <- c(14, 15)
m2[1, ] <- c(v1[1], v1[2], 16)
m2[2:3, 1] <- v1
m2[2:3, 2:3] <- m1[2:3, 2:3]
print(eigenMat(m2))

# ---- MatrixBase_setZero / setOnes / setIdentity / setRandom -----------------
# cpp 对某一行/整块调用 setXxx()。R 里是给下标赋相应的常量矩阵。
m <- eigenMat(matrix(1:16, 4, 4, byrow = TRUE))
m[2, ] <- 0
print(m)

m <- eigenMat(matrix(1:16, 4, 4, byrow = TRUE))
m[2, ] <- 1
print(m)

m <- eigenMat(matrix(0, 4, 4))
m[1:3, 1:3] <- eigen_identity(3)
print(m)

set.seed(11)
m <- eigen_random(4, 4)
print(round(m, 4))

# ---- MatrixBase_{ones,zero,identity,random}_int*.cpp -----------------------
# 这些片段是"不带参数的静态构造"，等价于 R 里给 nrow/ncol 传标量。
print(eigen_ones(3, 3))
print(eigen_zeros(2, 2))
print(eigen_identity(2))
print(length(eigen_random(5, 1)))     # VectorXd::Random(5)

# ---- Map_simple.cpp ---------------------------------------------------------
# cpp 的 Map 是"把一段裸内存当成矩阵看"，零拷贝。R 里没有裸内存这一层：
# eigenMat(x) 一定拷贝。能对应的是"同一串数据换一种形状/顺序解读"，
# 用 matrix(x, nrow, ncol, byrow=) 表达。
x <- 0:8
print(eigenMat(matrix(x, 3, 3)))                     # 列优先：0 3 6 / 1 4 7 / 2 5 8
print(eigenMat(matrix(x, 3, 3, byrow = TRUE)))       # 行优先：0 1 2 / 3 4 5 / 6 7 8

# ---- Map_inner_stride.cpp / Map_outer_stride.cpp / Map_general_stride.cpp ---
# 这三种 stride 在 R 里就是带步长的下标：seq(1, by = 2, length.out = n)。
x <- 0:11
print(x[seq(1, by = 2, length.out = 6)])             # InnerStride<2>
# OuterStride<>(4)：行内连续、行间隔 4，等价于取下标 0,1,2 / 4,5,6 / 8,9,10
print(matrix(x[c(0:2, 4:6, 8:10) + 1], 3, 3))

# ---- Matrix_Map_stride.cpp --------------------------------------------------
# cpp 从 A(1,1) 开始按 Stride<8,2> 取 2x2：元素是 A(1,1) A(1,3) / A(3,1) A(3,3)，
# 下标从 0 数。R 里直接写这两个下标更清楚。

A <- matrix(1:16, 4, 4, byrow = TRUE)
print(eigenMat(A[c(2, 4), c(2, 4)]))

# ---- Map_placement_new.cpp --------------------------------------------------
# cpp 用 placement new 让同一个 Map 对象改指向另一段内存。R 里就是重新赋值，
# 没有"就地改指向"这个概念（也不需要）。
v <- c(1, 2, 3, 4)
print(eigenMat(matrix(v, 1, 4)))
v <- c(5, 6, 7, 8, 9)
print(eigenMat(matrix(v, 1, 5)))

# ---- Tutorial_Map_rowmajor.cpp ----------------------------------------------
# 同一段 8 个整数，按列优先和按行优先读出来是两回事。
x <- 0:7
cat("Column-major:\n"); print(matrix(x, 2, 4))
cat("Row-major:\n");    print(matrix(x, 2, 4, byrow = TRUE))

# ---- Tutorial_Map_using.cpp -------------------------------------------------
# cpp 的要点是 Map 与 Map<const> 共享内存：m2map(3) = 7 会改到 m2。
# R 没有零拷贝视图，eigenMat(m2) 是独立副本，改副本不影响原件 —— 这一点
# 与 cpp 相反，写代码时值得记一笔。

set.seed(3)
m2 <- eigen_random(1, 5)
map <- eigenMat(as.matrix(m2))       # 拷贝
map[4] <- 7
cat("原对象第 4 个系数:", round(as.matrix(m2)[4], 4), "\n")
cat("副本第 4 个系数:", as.matrix(map)[4], "\n")

# ---- tut_matrix_coefficient_accessors.cpp -----------------------------------
# cpp 逐个元素写入 m(i,j) 与 v(i)。R 的 [i, j] 是 1-based，与 cpp 的 0-based 相差 1。
m <- matrix(0, 2, 2)
m[1, 1] <- 3; m[2, 1] <- 2.5; m[1, 2] <- -1
m[2, 2] <- m[2, 1] + m[1, 2]
print(eigenMat(m))
v <- numeric(2); v[1] <- 4; v[2] <- v[1] - 1
print(eigenVec(v))

# ---- DenseBase_LinSpaced.cpp / DenseBase_LinSpacedInt.cpp / ----------------
#      DenseBase_setLinSpaced.cpp
# cpp:  VectorXi::LinSpaced(4, 7, 10) / VectorXd::LinSpaced(5, 0, 1)
# R:    eigen_linspace(from, to, n)，整数端点时步长仍是整数
print(eigen_linspace(7, 10, 4))
print(eigen_linspace(0, 1, 5))

# _setLinSpaced 是"原地重填"，R 里就是重新赋值一个新向量：
v <- eigen_linspace(7, 10, 4)
print(v)

# ---- DenseBase_LinSpaced_seq_deprecated.cpp --------------------------------
# cpp 里那个 LinSpaced(size, low, high) 的旧重载已经废弃，语义与上面一致。
# R 侧只有一个 eigen_linspace()，不存在废弃版本。

# ---- TopicStorageOrders_example.cpp -----------------------------------------
# cpp 打印同一矩阵的列优先存储序列与行优先存储序列。R 的 as.vector() 给的是
# 列优先序列；要行优先，先转置再取。

A <- matrix(c(8, 2, 2, 9, 9, 1, 4, 4, 3, 5, 4, 5), 3, 4, byrow = TRUE)
print(eigenMat(A))
cat("列优先存储:", as.vector(A), "\n")
cat("行优先存储:", as.vector(t(A)), "\n")

# ---- MatrixBase_reshaped_int_int.cpp ----------------------------------------
# cpp:  m.reshaped(2, 8) 返回**视图**，不改 m。R 的 dim<- 是原地改，要视图
# 就先拷贝一份再改形状（eigenMat 的拷贝是 C++ 侧的深拷贝）。

set.seed(5)
m <- eigenMat(matrix(1:16, 4, 4, byrow = TRUE))
print(m)
view <- eigenMat(as.matrix(m))
dim(view) <- c(2L, 8L)
print(view)
cat("原矩阵尺寸未变:", dim(m)[1], "x", dim(m)[2], "\n")

# ---- MatrixBase_reshaped_to_vector.cpp --------------------------------------
# cpp 比较 reshaped() 与 reshaped<RowMajor>() 两种展平顺序。
m <- matrix(1:16, 4, 4, byrow = TRUE)
cat("按列展平:", as.vector(m), "\n")
cat("按行展平:", as.vector(t(m)), "\n")

# ---- Tutorial_reshaped_vs_resize_1.cpp / _2.cpp -----------------------------
# cpp 用同一个矩阵演示：reshaped() 不改原对象，resize() 改。R 里对应
# eigen_reshape()（返回新对象）与 dim<-（原地改）。
m <- eigenMat(matrix(1:16, 4, 4, byrow = TRUE))
print(eigen_reshape(m, 2, 8))        # 新对象
print(dim(m))                        # 原对象还是 4x4

# ---- Tutorial_ReshapeMat2Mat.cpp --------------------------------------------
# 2x6 的列优先数据重排成 6x2，就是 R 的 dim<-（数据不动，形状变）。
M1 <- matrix(1:12, 2, 6, byrow = TRUE)
M2 <- M1
dim(M2) <- c(6L, 2L)
print(eigenMat(M2))

# ---- Tutorial_ReshapeMat2Vec.cpp --------------------------------------------
# cpp 里列优先矩阵按列展平得到 1 4 7 2 5 8 3 6 9；行优先矩阵按行展平得到
# 1 2 3 4 5 6 7 8 9。R 里 as.vector() 就是列优先展平，行优先转置即可。
M1 <- matrix(1:9, 3, 3, byrow = TRUE)     # 内容是 1..9，行优先存放
cat("列优先展平:", as.vector(M1), "\n")
cat("行优先展平:", as.vector(t(M1)), "\n")

# ---- MatrixBase_cast.cpp ----------------------------------------------------
# cpp 的 cast<double>() 是把 float 矩阵的元素类型换掉。R 里只有 double 与
# complex 两种存储，类型转换在构造之前用 as.numeric()/as.integer() 完成。
md <- diag(2) * 0.45
mf <- diag(2) * 1
print(eigenMat(md + mf))                  # 逐元素相加，相当于 md + mf.cast<double>()
print(eigenMat(matrix(as.integer(c(1.7, 2.2)), 2, 1)))   # 截断成整数

# ---- Tutorial_range_for_loop_1d_cxx11.cpp / --------------------------------
#      Tutorial_range_for_loop_2d_cxx11.cpp /
#      MatrixBase_colwise_iterator_cxx11.cpp
# cpp 用 range-for 遍历系数（1d 走 vector，2d 走 reshaped()，都是列优先顺序）。
# R 里就是 for 循环或 apply 系列，顺序同样是列优先。
m <- matrix(1:6, 2, 3, byrow = TRUE)
for (x in as.vector(m)) cat(x, " ")      # 1d：列优先
cat("\n")
for (j in seq_len(ncol(m))) for (i in seq_len(nrow(m))) cat(m[i, j], " ")
cat("\n")

# ---- Tutorial_std_sort.cpp / Tutorial_std_sort_rows_cxx11.cpp ---------------
# cpp 对整块排序、对每一行分别排序。R 里是 sort() 与 apply(..., 1, sort)。
A <- matrix(c(3, 1, 2, 6, 5, 4), 2, 3, byrow = TRUE)
print(eigenMat(matrix(sort(as.vector(A)), 2, 3)))       # 整体排序后重新成矩阵
print(eigenMat(t(apply(A, 1, sort))))                   # 每行各自排序

cat("\n02-matrix-class.R 结束\n")
