# ===========================================================================
# 04 · 块、切片与视图
#
# 对应 EigenExample：
#   examples/Tutorial_BlockOperations_{print_block,block_assignment,colrow,corner,
#                                     vector}.cpp
#   snippets/MatrixBase_block_int_int.cpp / _int_int_int_int.cpp / fixedBlock_int_int
#   snippets/MatrixBase_{col,row}.cpp
#   snippets/MatrixBase_{topRows,bottomRows,leftCols,rightCols}_int.cpp
#   snippets/MatrixBase_{topLeftCorner,bottomLeftCorner,topRightCorner,
#                        bottomRightCorner}_int_int.cpp
#   snippets/MatrixBase_template_int_*.cpp                      （16 个）
#   snippets/MatrixBase_segment_int_int.cpp / start_int.cpp / end_int.cpp
#   snippets/MatrixBase_{diagonal,diagonal_int,diagonal_template_int,asDiagonal}.cpp
#   snippets/MatrixBase_{replicate,replicate_int_int,reverse,select}.cpp
#   snippets/MatrixBase_{eval,noalias,applyOnTheLeft,applyOnTheRight}.cpp
#   snippets/MatrixBase_triangularView.cpp / selfadjointView.cpp
#   snippets/MatrixBase_homogeneous.cpp / MatrixBase_hnormalized.cpp
#   snippets/VectorwiseOp_homogeneous.cpp / DirectionWise_hnormalized.cpp
#   snippets/class_{Block,FixedBlock,VectorBlock,FixedVectorBlock,Reshaped,
#                   FixedReshaped,CwiseUnaryOp,CwiseUnaryOp_ptrfun,CwiseBinaryOp}.cpp
#   snippets/Slicing_{arrayexpr,custom_padding_cxx11,rawarray_cxx11,stdvector_cxx11}.cpp
#   snippets/Tutorial_SlicingCol.cpp / Tutorial_SlicingVec.cpp
#
# 运行：Rscript 04-blocks.R
#
# 两条对照原则：
#   1. Eigen 的 block()/corner()/segment() 全都能写成 R 的下标表达式，
#      而且 R 的下标**直接就是拷贝**，没有"视图别名"这一层。
#   2. Eigen 的模板参数形式（block<2,2>、bottomRows<2>）与运行期参数形式
#      （block(2,2)、bottomRows(2)）在 R 里是同一个东西：R 没有编译期尺寸。
#      所以 template_int_* 那一批片段统一折进运行期形式的例子里。
# ===========================================================================

suppressPackageStartupMessages(library(ROpenEigen))

m4 <- matrix(1:16, 4, 4, byrow = TRUE)     # Eigen 的 m << 1..16 按行填

# ---- Tutorial_BlockOperations_print_block.cpp -------------------------------
# cpp:  m.block<2,2>(1,1)   /  m.block(0, 0, i, i)  /  m.block(0, 0, j1, j2)
# R:    m[2:3, 2:3]         /  m[1:i, 1:i]          /  m[1:j1, 1:j2]
# 注意 Eigen 的 block(起始行, 起始列, 行数, 列数) 是 0-based，R 是 1-based。
print(eigenMat(m4[2:3, 2:3]))

for (i in 1:3) {
  cat("Block of size", i, "x", i, "\n")
  print(eigenMat(m4[1:i, 1:i]))
}

j1 <- 2; j2 <- 3
print(eigenMat(m4[1:j1, 1:j2]))

# ---- MatrixBase_block_int_int.cpp -------------------------------------------
# cpp:  m.block<2,2>(1,1) = 0 —— 对块赋值
# R:    一次下标赋值（赋值也是拷贝语义，原矩阵真的被改）
m <- eigenMat(m4)
m[2:3, 2:3] <- 0
print(m)

# ---- MatrixBase_block_int_int_int_int.cpp -----------------------------------
# cpp 的 block(i, j, p, q) 四个参数版，与上面等价，只是尺寸是运行期变量。
i <- 1; j <- 1; p <- 2; q <- 2
print(eigenMat(m4[(i + 1):(i + p), (j + 1):(j + q)]))

# ---- MatrixBase_fixedBlock_int_int.cpp --------------------------------------
# cpp 的 fixedBlock<2,2>(i,j) 是"模板尺寸 + 运行期起点"，R 里没有区别：
print(eigenMat(m4[2:3, 1:2]))

# ---- Tutorial_BlockOperations_block_assignment.cpp --------------------------
# cpp 把 2x2 块写进 4x4 的中间，再把右下 2x3 块拷到左上。
# R 里同一件事就是两次下标赋值；而且右值先被求值，不存在 cpp 注释里担心的别名。
a <- matrix(0.6, 4, 4)
a[2:3, 2:3] <- matrix(c(1, 2, 3, 4), 2, byrow = TRUE)
print(eigenMat(a))
a[1:2, 1:3] <- a[3:4, 2:4]
print(eigenMat(a))

# ---- Tutorial_BlockOperations_colrow.cpp ------------------------------------
# cpp:  m.row(1) 与 m.col(2) += 3 * m.col(0)
m <- matrix(1:9, 3, 3, byrow = TRUE)
print(eigenMat(m[2, , drop = FALSE]))          # 第 2 行（0-based 的 row(1)）
m[, 3] <- m[, 3] + 3 * m[, 1]
print(eigenMat(m))

# ---- MatrixBase_col.cpp / MatrixBase_row.cpp --------------------------------
m <- diag(3)
m[, 2] <- c(4, 5, 6)
print(eigenMat(m))
print(eigenMat(m[1, ]))

# ---- Tutorial_BlockOperations_corner.cpp ------------------------------------
# cpp:  leftCols(2), bottomRows<2>(), topLeftCorner(2,2), topRightCorner<2,3>()
print(eigenMat(m4[, 1:2]))
print(eigenMat(m4[3:4, ]))
print(eigenMat(m4[1:2, 1:2]))
print(eigenMat(m4[1:2, 3:4]))

# ---- MatrixBase_topRows_int / bottomRows_int / leftCols_int / rightCols_int --
# 这四个是最常用的角/边，R 里就是取前 n 行 / 后 n 行 / 前 n 列 / 后 n 列。
n <- 2
print(eigenMat(m4[1:n, ]))                       # topRows(2)
print(eigenMat(m4[(nrow(m4) - n + 1):nrow(m4), ]))  # bottomRows(2)
print(eigenMat(m4[, 1:n]))                       # leftCols(2)
print(eigenMat(m4[, (ncol(m4) - n + 1):ncol(m4)]))  # rightCols(2)

# ---- MatrixBase_*Corner_int_int --------------------------------------------
# 四个角：topLeft / bottomLeft / topRight / bottomRight，参数是行数、列数。
p <- 2; q <- 2
print(eigenMat(m4[1:p, 1:q]))
print(eigenMat(m4[(nrow(m4) - p + 1):nrow(m4), 1:q]))
print(eigenMat(m4[1:p, (ncol(m4) - q + 1):ncol(m4)]))
print(eigenMat(m4[(nrow(m4) - p + 1):nrow(m4), (ncol(m4) - q + 1):ncol(m4)]))

# ---- MatrixBase_template_int_{topRows,bottomRows,leftCols,rightCols,
#                              topLeftCorner,bottomLeftCorner,topRightCorner,
#                              bottomRightCorner,segment,end,start,int_block,...}
# 这一批 16 个片段是同一件事的模板参数写法。R 没有编译期尺寸，统一用上面
# 运行期写法即可，不再重复。

# ---- Tutorial_BlockOperations_vector.cpp ------------------------------------
# cpp 对向量取 head / tail / segment。
v <- c(1, 2, 3, 4, 5, 6)
print(v[1:2])                       # head(2)
print(v[(length(v) - 1):length(v)]) # tail(2)
print(v[3:4])                       # segment(2, 2) —— 从下标 2（0-based）起 2 个

# ---- MatrixBase_segment_int_int.cpp / MatrixBase_start_int.cpp / ------------
#      MatrixBase_end_int.cpp
# segment(i, n) 是从第 i 个（0-based）起取 n 个；start()/end() 是首尾元素的
# 单元素块。R 里就是 v[i+1] 与 v[length(v)]。
v <- eigenVec(c(0, 1, 2, 3, 4, 5))
print(v[3:5])          # segment(2, 3)
print(v[1])            # start()
print(v[length(v)])    # end()

# ---- DenseBase_middleCols_int.cpp / DenseBase_middleRows_int.cpp -----------
# cpp 的 middleCols(i, n) / middleRows(i, n) 是"从第 i 条起取 n 条"。
# R 里就是 m[, (i+1):(i+n)] 与 m[(i+1):(i+n), ]。
A <- matrix(1:25, 5, 5, byrow = TRUE)          # cpp 里是 5x5 随机矩阵
print(eigenMat(A[, 2:4]))                       # A.middleCols(1, 3)
print(eigenMat(A[3:4, ]))                       # A.middleRows(2, 2)

# ---- DenseBase_template_int_middleCols.cpp / -------------------------------
#      DenseBase_template_int_middleRows.cpp
# 模板参数版本（middleCols<3>(i)）在 R 里没有区别，同上一段。

# ---- MatrixBase_diagonal.cpp / diagonal_int.cpp / diagonal_template_int.cpp -
# cpp 的 diagonal() 取主对角线，diagonal(k) 取第 k 条（正为上，负为下）。
# R 侧对应物是 eigen_diag()，它**只做提取**（返回 k x 1），
# 与 base R 的 diag() 在矩阵输入下一致；建对角矩阵是另一个方向。
m <- matrix(1:16, 4, 4, byrow = TRUE)
print(eigen_diag(eigenMat(m)))                  # diagonal()
print(diag(m))                                  # base R 对照
print(m[cbind(1:3, 2:4)])                       # diagonal(1)，上对角线
print(m[cbind(2:4, 1:3)])                       # diagonal(-1)，下对角线

# ---- MatrixBase_asDiagonal.cpp ----------------------------------------------
# cpp:  v.asDiagonal() 把向量变成对角矩阵。R 里是 base R 的 diag(v)
# —— 注意方向相反：eigen_diag() 是"矩阵 → 对角线"，diag() 是双向的。
print(eigenMat(diag(c(1, 2, 3))))

# ---- MatrixBase_replicate.cpp / replicate_int_int.cpp -----------------------
# cpp:  v.replicate(1, 3) 把行向量重复 3 次；m.replicate(2, 2) 各方向 2 次。
# R:    eigen_rep(x, nrow, ncol)
v <- eigenVec(c(1, 2, 3))
print(eigen_rep(v, 1, 3))
print(eigen_rep(eigenMat(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)), 2, 2))

# ---- MatrixBase_reverse.cpp -------------------------------------------------
# cpp 的 reverse() 按**存储顺序**（列优先）反转整个矩阵，返回的是一个视图，
# 所以 m.reverse()(1,0) = 4 会真的改到 m。R 里 rev() 只作用于向量，
# 矩阵要先 as.vector 再按原尺寸放回；改不到原对象（拷贝语义）。
m <- matrix(1:12, 3, 4, byrow = TRUE)
print(eigenMat(matrix(rev(as.vector(m)), 3, 4)))

# ---- MatrixBase_select.cpp --------------------------------------------------
# cpp:  m = (m.array() >= 5).select(-m, m)
# R:    ifelse(条件, -m, m)
m <- matrix(1:9, 3, 3, byrow = TRUE)
print(eigenMat(ifelse(m >= 5, -m, m)))

# ---- MatrixBase_eval.cpp / noalias.cpp --------------------------------------
# cpp 的 .eval() 强制求值、.noalias() 声明"右值不别名、可以直接写"。R 里
# 每次运算都生成新对象，赋值也立即完成，这两个修饰都没有对应物。
a <- eigenMat(matrix(1:4, 2, 2, byrow = TRUE))
print(a %*% a)          # 不需要 eval()
print(a)                # 也不需要 noalias()：左值右值绝不会是同一块内存

# ---- MatrixBase_applyOnTheLeft.cpp / applyOnTheRight.cpp --------------------
# cpp:  m.applyOnTheLeft(A) 等价于 m = A * m；applyOnTheRight(B) 等价于 m = m * B。
# R 里直接写矩阵乘，再赋值。
m <- eigenMat(matrix(1:6, 2, 3, byrow = TRUE))
A <- eigenMat(matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE))
print(A %*% m)
print(m %*% eigen_identity(3))

# ---- MatrixBase_triangularView.cpp ------------------------------------------
# cpp 的 triangularView<Upper>() 把另一半当成 0。R 里用 upper.tri()/lower.tri()
# 显式清零，或者直接用 backsolve()/forwardsolve() 解三角系统。
m <- matrix(c(2, 1, 1, 3), 2, 2, byrow = TRUE)
upper <- m; upper[lower.tri(upper)] <- 0
print(eigenMat(upper))
b <- c(5, 7)
print(backsolve(upper, b))              # 上三角求解，与 cpp 的 .solve() 对应

# ---- MatrixBase_selfadjointView.cpp -----------------------------------------
# cpp 的 selfadjointView<Lower>() 只读下三角、当作对称矩阵。R 里没有这层
# 惰性视图，标准做法是手工补全对称部分再算。
lower <- matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE)   # 只下三角有意义
sym <- lower + t(lower) - diag(diag(lower))
print(eigenMat(sym))

# ---- MatrixBase_homogeneous.cpp / hnormalized.cpp ---------------------------
# cpp 的 homogeneous() 在末尾补 1，hnormalized() 是前 n-1 个分量除以最后一个。
v <- c(1, 2, 3)
print(c(v, 1))                          # homogeneous()
h <- c(2, 4, 6, 2)
print(h[1:3] / h[4])                    # hnormalized()

# ---- VectorwiseOp_homogeneous.cpp / DirectionWise_hnormalized.cpp -----------
# 对每一列/每一行做齐次化，R 里是 rbind 一行 1、或者按行归一化。
m <- matrix(1:6, 2, 3, byrow = TRUE)
print(eigenMat(rbind(m, 1)))            # 每列补 1（colwise().homogeneous()）
print(eigenMat(m / rep(colSums(m), each = nrow(m))))   # 每列按列和归一化

# ---- MatrixBase_colwise.cpp / rowwise.cpp / colwise_iterator_cxx11.cpp ------
# cpp 的 colwise()/rowwise() 返回"按列/按行"的视图，后面接归约或广播（见下一章）。
# 单纯取列/行在 R 里是 m[, j] / m[i, ]，遍历就是 apply(m, 2, f) / apply(m, 1, f)。
m <- matrix(1:6, 2, 3, byrow = TRUE)
for (j in seq_len(ncol(m))) cat("col", j, ":", m[, j], "\n")
for (i in seq_len(nrow(m))) cat("row", i, ":", m[i, ], "\n")

# ---- class_Block.cpp / class_FixedBlock.cpp / class_VectorBlock.cpp / -------
# ---- class_FixedVectorBlock.cpp / class_Reshaped.cpp / class_FixedReshaped.cpp
# ---- class_CwiseUnaryOp.cpp / class_CwiseUnaryOp_ptrfun.cpp / class_CwiseBinaryOp.cpp
# 这些片段演示的是"块/运算都是**表达式**，可以嵌套、可以延迟求值"。
# R 没有表达式模板：m[1:2, 1:2] 立刻就是一块新内存，abs(m) 也立刻算完。
# 所以这一组在 R 里没有专门写法，上面的下标表达式已经覆盖了它们的语义。
print(eigenMat(abs(m4[1:2, 1:2] + m4[3:4, 3:4])))

# ---- Slicing_arrayexpr.cpp / Slicing_rawarray_cxx11.cpp / -------------------
# ---- Slicing_stdvector_cxx11.cpp / Slicing_custom_padding_cxx11.cpp ---------
# cpp 这几段的核心是"对一段连续内存按任意起点/步长切片"。R 里就是带步长的
# seq()：x[seq(from, by = step, length.out = n)]。
x <- 0:19
print(x[seq(1, by = 3, length.out = 4)])        # 每 3 个取 1 个

# ---- Tutorial_SlicingCol.cpp ------------------------------------------------
# cpp 把 3x8 矩阵每隔 3 列取一列（列优先存储下等价于隔 3 个元素取 1 个）。
M1 <- matrix(1:24, 3, 8, byrow = TRUE)
print(eigenMat(M1[, seq(1, by = 3, length.out = 3)]))

# ---- Tutorial_SlicingVec.cpp ------------------------------------------------
v <- c(0, 1, 2, 3, 4, 5, 6, 7)
print(v[seq(1, by = 2, length.out = 4)])

cat("\n04-blocks.R 结束\n")
