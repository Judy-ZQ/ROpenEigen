# ROpenEigen 用例教程

这套脚本把 Eigen 官方 `EigenExample/` 目录里的 **377 个 `.cpp` / `.h`**
（examples 70、snippets 297、special_examples 3、demos 7 组）逐个改写成
ROpenEigen 的 R 写法，按主题分成 11 章。每个 cpp 都能在某章里找到对应的
代码段，章首注释里写着它到底对上了哪些文件。

目的不是"把 C++ 翻译成 R"这么简单，而是把**哪些概念在 R 里有对等物、
哪些没有、没有的该怎么绕**讲清楚。有对等物的部分给出可直接运行的 R 代码；
没有的（表达式模板、模板元编程、编译期尺寸……）集中收在第 11 章逐条说明。

## 怎么跑

```bash
# 装好包之后（R CMD INSTALL ROpenEigen_0.2.0.tar.gz）
Rscript run-all.R              # 顺序跑完 11 章
Rscript run-all.R 03 07        # 只跑第 3、7 章
Rscript run-all.R --quiet      # 只要每章一行的摘要
Rscript 04-blocks.R            # 单章直接跑

# 包已安装但不想找路径时
Rscript "$(Rscript -e 'cat(system.file("examples","run-all.R",package="ROpenEigen"))')"
```

唯一依赖是 ROpenEigen 本身。第 07 章用到 `Matrix::Schur` 做对照，代码里带
`requireNamespace()` 守卫，没装 Matrix 时该小段自动跳过，不影响其他部分。

## 章节一览

| 章 | 主题 | 覆盖文件 | 代表片段 |
|----|------|---------:|----------|
| 01 | 快速开始：构造、逗号初始化、固定/动态尺寸 | 14 | `QuickStart_example`, `Tutorial_commainit_02` |
| 02 | 矩阵与向量：尺寸、存储顺序、Map、reshape、cast | 66 | `Matrix_resize_int_int`, `TopicStorageOrders_example` |
| 03 | 算术运算与逐元素函数（Cwise 全家桶） | 85 | `tut_arithmetic_matrix_mul`, `Cwise_*` 52 个 |
| 04 | 块、切片与视图 | 75 | `Tutorial_BlockOperations_block_assignment` |
| 05 | 归约、访问者与广播 | 37 | `Tutorial_ReductionsVisitorsBroadcasting_colwise` |
| 06 | 线性方程组求解 | 23 | `TutorialLinAlgExSolveColPivHouseholderQR` |
| 07 | 矩阵分解：LU / QR / Cholesky / SVD / Schur / QZ | 50 | `LLT_example`, `JacobiSVD_basic` |
| 08 | 稀疏矩阵与迭代解法 | 5 | `Tutorial_sparse_example`, `BiCGSTAB_step_by_step` |
| 09 | 几何：四元数、轴角、变换、umeyama | 6 | `AngleAxis_mimic_euler`, `demos/opengl/trackball.cpp` |
| 10 | 特殊函数（unsupported/SpecialFunctions） | 3 | `Cwise_erf`, `Cwise_lgamma` |
| 11 | 只有 C++ 才有的东西（逐条说明为什么） | 19 | `TopicAliasing_mult1`, `function_taking_ref` |

合计 383 处引用、去重后 377 个文件（6 个文件跨越两章，如 `Cwise_erf` 既在
03 章的 Cwise 表里，也在 10 章单独展开）。

## 动手前先记住六条

这六条能解释教程里 90% 的"为什么这么写"。

1. **下标**：Eigen 0-based，R 1-based。`v(i)` 对应 `v[i + 1]`。
2. **填充顺序**：Eigen 的逗号初始化**按行**填，R 的 `matrix()` 默认**按列**填，
   所以一律写 `byrow = TRUE`。
3. **两套语义 vs 一套**：Eigen 的 `Matrix * Matrix` 是矩阵乘、`Array * Array`
   是逐元素；R 只有一套 —— `*` 永远逐元素，矩阵乘写 `%*%`。
   于是 Eigen 的 Array 片段在 R 里全都能用最朴素的写法表达。
4. **值语义**：Eigen 有表达式模板、别名问题、`eval()`/`noalias()`；R 每次
   运算都产出新对象，`as.matrix()` 是深拷贝。所以一整类"别名陷阱"在 R 里
   根本不存在。
5. **编译期 vs 运行期**：`block<2,2>(i,j)` 与 `block(i,j,2,2)` 在 R 里是同一个
   东西，`template_int_*` 那一批片段统一折进运行期写法的例子里。
6. **视图 vs 拷贝**：Eigen 的 `block()`、`col()`、`diagonal()` 是视图（改它
   就改了原矩阵）；R 的下标是拷贝，想改原矩阵要写 `m[i, j] <- value`。

## API 速查

左边是 Eigen，右边是 ROpenEigen（或 base R）。这张表按章分组，详细用法看对应章节。

**构造与生成**

| Eigen | R |
|-------|---|
| `MatrixXd m(2, 3);` | `eigenMat(matrix(0, 2, 3))` |
| `MatrixXd::Zero(2, 3)` | `eigen_zeros(2, 3)` |
| `MatrixXd::Ones(2, 3)` | `eigen_ones(2, 3)` |
| `MatrixXd::Constant(2, 3, 7)` | `eigen_constant(7, 2, 3)` |
| `MatrixXd::Identity(3, 3)` | `eigen_identity(3)` |
| `MatrixXd::Random(2, 3)` | `eigen_random(2, 3)` |
| `VectorXd::LinSpaced(5, 0, 1)` | `eigen_linspace(0, 1, 5)` |
| `m << 1, 2, 3, 4;`（按行） | `eigenMat(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))` |
| `v << 1, 2, 3;` | `eigenVec(c(1, 2, 3))` |

**形状与存储**

| Eigen | R |
|-------|---|
| `m.rows()` / `m.cols()` | `nrow(m)` / `ncol(m)` |
| `m.size()` | `length(m)` |
| `m.resize(r, c)` | `dim(m) <- c(r, c)`（不能改元素总数） |
| `m.setZero()` / `setOnes()` / `setIdentity()` | 重新赋值 `eigen_zeros()` 等 |
| `m.transpose()` | `t(m)` |
| `m.adjoint()` | `eigen_adjoint(m)` / `Conj(t(m))` |
| `m.reshaped(r, c)` | `eigen_reshape(m, r, c)` |
| `v.replicate(2, 3)` | `eigen_rep(v, 2, 3)` |
| `MatrixXd::Map(p, r, c)` | 无（R 没有裸内存视图，见 02 章） |

**块与切片**

| Eigen | R |
|-------|---|
| `m.block(i, j, p, q)` | `m[(i+1):(i+p), (j+1):(j+q)]` / `eigen_block(m, ...)` |
| `m.col(j)` / `m.row(i)` | `m[, j + 1]` / `m[i + 1, ]` |
| `v.segment(i, n)` | `v[(i+1):(i+n)]` |
| `m.topRows(n)` / `bottomRows(n)` | `m[1:n, ]` / `m[(nrow(m)-n+1):nrow(m), ]` |
| `m.leftCols(n)` / `rightCols(n)` | `m[, 1:n]` / `m[, (ncol(m)-n+1):ncol(m)]` |
| `m.topLeftCorner(p, q)` | `m[1:p, 1:q]` |
| `m.diagonal()` | `eigen_diag(m)`（只提取，不构造） |
| `v.asDiagonal()` | `eigenMat(diag(v))` |
| `m.reverse()` | `m[nrow(m):1, ncol(m):1]` / `rev(v)` |
| `m.triangularView<Upper>().solve(b)` | `backsolve(m, b)`（下三角用 `forwardsolve`） |
| `m/selfadjointView<Lower>()` | 先把下三角补成对称：`m + t(m) - diag(diag(m))` |

**归约与内积**

| Eigen | R |
|-------|---|
| `m.sum()` / `prod()` / `mean()` | `sum(m)` / `prod(m)` / `mean(m)` |
| `m.minCoeff()` / `maxCoeff()` | `min(m)` / `max(m)` |
| `m.norm()`（Frobenius） | `eigen_norm(m, "F")`（默认 `"O"` 是 1-范数） |
| `m.operatorNorm()`（谱范数） | `eigen_norm(m, "2")` |
| `m.squaredNorm()` | `sum(m^2)` |
| `m.trace()` | `eigen_trace(m)` |
| `m.all()` / `m.any()` | `all(m)` / `any(m)` |
| `m.colwise().sum()` | `colSums(m)` / `apply(m, 2, sum)` |
| `m.rowwise().maxCoeff()` | `apply(m, 1, max)` |
| `v.dot(w)` | `eigen_dot(v, w)` / `sum(v * w)` |
| `v.cross(w)` | `eigen_cross(v, w)`（2 维返回标量 z） |
| `m.cwiseProduct(n)` | `m * n` |

**求解与分解**

| Eigen | R |
|-------|---|
| `A.inverse()` | `eigen_inverse(A)` / `solve(A)` |
| `A.determinant()` | `eigen_det(A)` / `det(A)` |
| `A.llt().solve(b)` | `eigen_llt(A)` / `eigen_chol(A)` / `eigen_solve(A, b)` |
| `A.ldlt().solve(b)` | `eigen_ldlt(A)` |
| `A.partialPivLu().solve(b)` | `eigen_lu(A)` |
| `A.fullPivLu()` | `eigen_full_lu(A)` |
| `A.householderQr()` / `colPivHouseholderQr()` | `eigen_qr(A)` / `eigen_colpiv_qr(A)` |
| `A.jacobiSvd()` | `eigen_svd(A)` |
| `A.selfadjointView<Lower>().eigenvalues()` | `eigen_eigen(A)` |
| `A.eigenvalues()` | `eigen_values(A)` |
| `A.eigenvectors()` | `eigen_vectors(A)` |
| `A.jacobiSvd().singularValues()` | `eigen_singular_values(A)` |
| `A.rank()` / `A.rcond()` | `eigen_rank(A)` / `eigen_rcond(A)` |

**稀疏**

| Eigen | R |
|-------|---|
| `SparseMatrix<double> A(m, n); A.insert(i,j) = x;` | `eigen_sparse(i, j, x, nrow = m, ncol = n)` |
| `A.nonZeros()` | `eigen_sparse_nnz(A)` |
| `A.coeff(i, j)` | `eigen_sparse_get(A, i, j)` |
| `SparseMatrix` + `SimplicialLLT` / `BiCGSTAB` / `GMRES` … | `eigen_sparse_solve(A, b, method = "llt")` 等 11 种 |

**几何与特殊函数**

| Eigen | R |
|-------|---|
| `Quaterniond(w, x, y, z)` | `eigen_quat(w, x, y, z)` |
| `AngleAxisd(angle, axis)` | `eigen_angle_axis(angle, axis)` |
| `q.slerp(t, r)` | `eigen_quat_slerp(t, q, r)` |
| `Transform<double,3,Affine>` | `eigen_transform()` + `$rotate()` / `$translate()` |
| `umeyama(src, dst)` | `eigen_umeyama(src, dst)` |
| `x.erf()` / `x.lgamma()` | `eigen_erf(x)` / `eigen_gammaln(x)` |
| `igamma(a, x)` / `zeta(x)` / `bessel_*` | `eigen_igamma` / `eigen_zeta` / `eigen_bessel` |
| `A.exp()` / `log()` / `sqrt()` / `pow(p)`（矩阵函数） | `eigen_expm` / `eigen_logm` / `eigen_sqrtm` / `eigen_powm` |

## 没有对等物的部分

集中在第 11 章，一条一条说明，不是"跳过"：

- **表达式模板与别名**（`TopicAliasing_*`）：R 每次运算都产出新对象，
  `a = a * a` 里的右值先算完才赋值，别名问题天然不存在。
- **`eval()` / `noalias()`**：没有惰性表达式，也就没有强制求值和跳过临时量。
- **`Ref<>` / `EigenBase<>` / 模板元编程**（`function_taking_*`、`TemplateKeyword_*`）：
  泛型在 R 里靠 S3 派发，写不出"接受任意 Eigen 表达式"的签名。
- **CRTP 自定义表达式**（`CustomizingEigen_Inheritance`、`nullary_indexing`、
  `make_circulant`）：R 里用普通函数返回矩阵即可，代价是失去惰性求值。
- **编译期尺寸**（`MatrixBase_template_int_*`、`Matrix_setXxx_int*`）：
  固定尺寸在 R 里没有性能意义，全部折进运行期写法。
- **`Map` 与 placement new**：R 没有裸内存视图，`eigenMat()` 一定拷贝。
- **Qt / OpenGL 的 GUI 层**（`demos/mandelbrot`、`demos/opengl`）：只取其中的
  数学内核改写（05 章末尾的 Mandelbrot 迭代、09 章末尾的 icosphere 细分）。
- **`demos/mix_eigen_and_c` 的 C 包装**：这一层在 ROpenEigen 内部已经存在
  （`eigenMat` 就是那个不透明句柄，`delete` 交给 GC 终结器），用户不用手写。

## 数字从哪来

教程里出现的每个数值（特征值、误差、迭代步数、棱长）都来自**实际运行**，
不是从 Eigen 文档抄的期望值。写这批脚本时先做了一份 C++ 侧基准：
把 297 个 snippet 逐个包成 `main()` 编译运行、71 个 example 直接编译，
输出留在 `scratch/eigendemo/`、`scratch/snip/`，再拿 R 的结果去对。
对不上的地方一律以实跑为准（也因此发现并修掉了两个 R 侧真实 bug：
向量右值铺块、`eigenMat` 右值赋值）。

想自己核对某一章：直接 `Rscript 03-arithmetic.R`，输出就是基准。

## 目录

```
inst/examples/
├── README.md              ← 本文件
├── run-all.R              ← 顺序执行入口
├── 01-quick-start.R       ← 快速开始
├── 02-matrix-class.R      ← 矩阵与向量
├── 03-arithmetic.R        ← 算术与逐元素
├── 04-blocks.R            ← 块与切片
├── 05-reductions.R        ← 归约与广播
├── 06-linear-solve.R      ← 线性求解
├── 07-decompositions.R    ← 矩阵分解
├── 08-sparse.R            ← 稀疏
├── 09-geometry.R          ← 几何
├── 10-special-functions.R ← 特殊函数
└── 11-cpp-only.R          ← 只有 C++ 才有的东西
```
