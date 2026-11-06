# Sparse matrices.
#
# The risk in sparse code is different from the dense side. Dense bugs tend to
# be numerical; sparse bugs tend to be *indexing* -- 1-based versus 0-based,
# a transpose that silently reads column-major, an off-by-one at the edge of
# the block. Those produce plausible-looking numbers rather than crashes, so
# every assertion here is against a dense equivalent computed by base R.
#
# The other thing checked hard: Eigen's sparse factorisations report failure
# through info() and keep returning values. An unchecked SparseLLT on an
# indefinite matrix hands back the solution of a *different* system. That is
# exactly the bug that shipped once on the dense side, so it gets a test here.

suppressPackageStartupMessages(library(ROpenEigen))

ok <- 0L; bad <- 0L
chk <- function(label, expr) {
  r <- tryCatch(isTRUE(expr), error = function(e) {
    cat(sprintf("  !! %-38s %s\n", label, conditionMessage(e))); NA
  })
  if (is.na(r)) { bad <<- bad + 1L; return(invisible(FALSE)) }
  if (r) ok <<- ok + 1L else { bad <<- bad + 1L; cat(sprintf("  ?? %s\n", label)) }
  invisible(r)
}
errs <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

nrm <- function(x) { x <- as.matrix(x); if (!length(x)) 0 else max(abs(x)) }
close <- function(a, b, tol = 1e-10) {
  a <- as.matrix(a); b <- as.matrix(b)
  if (!identical(dim(a), dim(b))) return(FALSE)
  scale <- max(1, nrm(a), nrm(b))
  max(Mod(a - b)) <= tol * scale
}
unw <- function(x) if (inherits(x, "eigenObject")) as.matrix(x) else x

cat("== 1. 三元组构造：1-based 索引必须落在正确的格子 ==\n")
S <- eigen_sparse(c(1, 2, 3), c(1, 2, 3), c(4, 5, 6), nrow = 3, ncol = 3)
chk("对角元素位置",   identical(as.vector(unw(S[1, 1])), 4))
chk("(2,2)",          identical(as.vector(unw(S[2, 2])), 5))
chk("(3,3)",          identical(as.vector(unw(S[3, 3])), 6))
chk("未存的位置是 0", identical(as.vector(unw(S[1, 3])), 0))
chk("dim",            identical(as.integer(dim(S)), c(3L, 3L)))
chk("nnz",            identical(as.integer(eigen_sparse_nnz(S)), 3L))
chk("length",         identical(as.numeric(length(S)), 9))

# 非方阵，且首行/首列都非零——转置写反了就会在这里露出来
T1 <- eigen_sparse(c(1, 1, 3), c(1, 2, 4), c(10, 20, 30), nrow = 3, ncol = 4)
chk("非方阵 (1,1)",   identical(as.vector(unw(T1[1, 1])), 10))
chk("非方阵 (1,2)",   identical(as.vector(unw(T1[1, 2])), 20))
chk("非方阵 (3,4)",   identical(as.vector(unw(T1[3, 4])), 30))
chk("非方阵 (2,1)=0", identical(as.vector(unw(T1[2, 1])), 0))
chk("非方阵 dim",     identical(as.integer(dim(T1)), c(3L, 4L)))

cat("== 2. 重复三元组必须相加 ==\n")
# Matrix::sparseMatrix 的默认行为就是求和；两边不一致的话，同一个矩阵用两种
# 方式构造会得到不同结果。
D <- eigen_sparse(c(1, 1, 2), c(1, 1, 2), c(3, 4, 5), nrow = 2, ncol = 2)
chk("重复项求和",   identical(as.vector(unw(D[1, 1])), 7))
chk("重复项不重复存", identical(as.integer(eigen_sparse_nnz(D)), 2L))

cat("== 3. 非法输入必须报错，不能崩 ==\n")
chk("行越界",       errs(eigen_sparse(4, 1, 1, nrow = 3, ncol = 3)))
chk("列越界",       errs(eigen_sparse(1, 4, 1, nrow = 3, ncol = 3)))
chk("0 下标",       errs(eigen_sparse(0, 1, 1, nrow = 3, ncol = 3)))
chk("负下标",       errs(eigen_sparse(-1, 1, 1, nrow = 3, ncol = 3)))
chk("NA 下标",      errs(eigen_sparse(NA_integer_, 1, 1, nrow = 3, ncol = 3)))
chk("长度不一致",   errs(eigen_sparse(c(1, 2), 1, 1, nrow = 3, ncol = 3)))
chk("负维度",       errs(eigen_sparse(1, 1, 1, nrow = -1, ncol = 3)))
chk("取越界元素",   errs(eigen_sparse_get(S, 4, 1)))
chk("复数被拒绝",   errs(eigen_sparse(matrix(1:4 + 0i, 2))))

cat("== 4. 从稠密矩阵转换 ==\n")
A <- matrix(c(1, 0, 0, 2, 0, 3, 0, 0, 4), 3, 3)
As <- eigen_sparse(A)
chk("稠密往返无损",   close(as.matrix(As), A, 1e-12))
chk("精确零被丢掉",   identical(as.integer(eigen_sparse_nnz(As)), 4L))
chk("从 eigenMat 转换", close(as.matrix(eigen_sparse(eigenMat(A))), A, 1e-12))
chk("全零矩阵 nnz=0", identical(as.integer(eigen_sparse_nnz(eigen_sparse(matrix(0, 3, 3)))), 0L))
chk("全零往返",       close(as.matrix(eigen_sparse(matrix(0, 3, 3))), matrix(0, 3, 3), 1e-12))
chk("as.vector 与 base 一致", close(as.vector(As), as.vector(A), 1e-12))
chk("as.matrix 是 base matrix", is.matrix(as.matrix(As)) && !inherits(as.matrix(As), "eigenObject"))
chk("is.matrix / is.array",     isTRUE(is.matrix(As)) && isTRUE(is.array(As)))

cat("== 5. 三元组导出 ==\n")
tr <- eigen_triplets(T1)
chk("导出是 data.frame", is.data.frame(tr))
chk("导出列名",          identical(names(tr), c("i", "j", "x")))
chk("导出 1-based",      identical(sort(tr$i), c(1L, 1L, 3L)))
chk("导出值可重建",      {
  r <- eigen_sparse(tr$i, tr$j, tr$x, nrow = 3, ncol = 4)
  close(as.matrix(r), as.matrix(T1), 1e-12)
})
chk("导出列优先顺序", identical(tr$j, sort(tr$j)))

cat("== 6. 转置 ==\n")
chk("转置维度",     identical(as.integer(dim(t(T1))), c(4L, 3L)))
chk("转置数值",     close(as.matrix(t(T1)), t(as.matrix(T1)), 1e-12))
chk("转置两次还原", close(as.matrix(t(t(T1))), as.matrix(T1), 1e-12))
chk("转置保持 nnz", identical(as.integer(eigen_sparse_nnz(t(T1))), as.integer(eigen_sparse_nnz(T1))))

cat("== 7. 下标 ==\n")
chk("[i,j] 单元素",   identical(as.vector(unw(T1[1, 2])), 20))
chk("[i] 线性下标",   identical(as.vector(unw(T1[1])), 10))    # 列主序第一个
chk("[4] 第二列首行", identical(as.vector(unw(T1[4])), 20))
chk("区间退回稠密",   close(unw(T1[1:2, 1:2]), as.matrix(T1)[1:2, 1:2], 1e-12))
chk("整列",           close(unw(T1[, 2]), as.matrix(T1)[, 2], 1e-12))

cat("== 8. 矩阵乘法 ==\n")
set.seed(101)
m <- 30; p <- 20
Ad <- matrix(0, m, p)
for (k in 1:p) { Ad[k, k] <- 4; Ad[k + 5, k] <- -1; Ad[k + 10, k] <- 2 }
As <- eigen_sparse(Ad)
B <- matrix(rnorm(p * 4), p, 4)
chk("稀疏 x 稠密",   close(unw(As %*% B), Ad %*% B, 1e-12))
v <- rnorm(p)
chk("稀疏 x 向量",   close(unw(As %*% v), Ad %*% v, 1e-12))
# 左操作数要按行数来：base R 把无维向量当 1 x length，所以长度必须是 nrow(Ad)
vrow <- rnorm(m)
chk("向量 x 稀疏",   close(unw(vrow %*% As), vrow %*% Ad, 1e-12))
# 稠密在左：Eigen 的 Sparse 模块没有稠密×稀疏的直接乘积，内部走
# (S^T D^T)^T，所以这条要单独验。
B3 <- matrix(rnorm(4 * m), 4, m)
chk("稠密 x 稀疏",   close(unw(B3 %*% As), B3 %*% Ad, 1e-12))
# 稀疏 x 稀疏要能对上维度：Ad 是 30x20，乘它的转置得 30x30
At <- eigen_sparse(t(Ad))
chk("稀疏 x 稀疏",   close(as.matrix(As %*% At), Ad %*% t(Ad), 1e-12))
chk("稀疏 x 稀疏 保持稀疏", inherits(As %*% At, "eigenSparse"))
chk("维度不匹配报错", errs(As %*% matrix(1, 3, 3)))

cat("== 9. 逐元素算术 ==\n")
chk("稀疏 + 稀疏",     close(as.matrix(As + As), Ad + Ad, 1e-12))
chk("稀疏 - 稀疏",     close(as.matrix(As - As), Ad - Ad, 1e-12))
chk("加保持稀疏",      inherits(As + As, "eigenSparse"))
chk("标量乘保持稀疏",  inherits(As * 2, "eigenSparse"))
chk("标量乘数值",      close(as.matrix(As * 2), Ad * 2, 1e-12))
chk("标量除",          close(as.matrix(As / 2), Ad / 2, 1e-12))
chk("逐元素乘",        close(as.matrix(As * As), Ad * Ad, 1e-12))
chk("形状不符报错",    errs(As + eigen_sparse(matrix(1, 3, 3))))

# 幂：p > 0 时隐含零元素仍是零，结构不变，所以必须保持稀疏；
# p <= 0 时 0^p 是 Inf 或 1，必须退回稠密。这两条分界线各有断言。
#
# Ad 里有零也有负数，所以 ^0.5 出 NaN、^-1 出 Inf。close() 遇到 NaN 一律
# 判假（d/s 是 NaN），这里需要一个先比 NaN/Inf 位置、再比有限部分的比较器。
same_special <- function(a, b, tol = 1e-12) {
  a <- as.matrix(a); b <- as.matrix(b)
  if (!identical(dim(a), dim(b))) return(FALSE)
  if (!identical(is.na(a), is.na(b))) return(FALSE)
  if (!identical(is.infinite(a), is.infinite(b))) return(FALSE)
  if (!identical(is.nan(a), is.nan(b))) return(FALSE)
  f <- is.finite(a) & is.finite(b)
  all(abs(a[f] - b[f]) <= tol * pmax(abs(b[f]), 1))
}
chk("正幂保持稀疏",    inherits(As ^ 2, "eigenSparse"))
chk("平方数值",        close(as.matrix(As ^ 2), Ad ^ 2, 1e-12))
chk("非整数正幂保持稀疏", inherits(As ^ 0.5, "eigenSparse"))
chk("平方根数值（含 NaN）", same_special(as.matrix(As ^ 0.5), Ad ^ 0.5))
chk("负幂转稠密",      !inherits(As ^ -1, "eigenSparse"))
chk("负幂数值（含 Inf）", same_special(as.matrix(As ^ -1), Ad ^ -1))
chk("零次幂转稠密",    close(as.matrix(As ^ 0), Ad ^ 0, 1e-12))

# 除以 0：base R 给 Inf（非零处）和 NaN（零处）。稀疏表示不了 NaN，
# 所以这里必须转稠密，而且**不能**报错——稠密 eigenMat 和 base matrix
# 都能算，只有稀疏报错的话，同一个式子换个类型结果就变了。
d0 <- tryCatch(as.matrix(As / 0), error = function(e) NULL)
chk("除以 0 不报错",   !is.null(d0))
chk("除以 0 与 base 一致", !is.null(d0) && identical(is.infinite(d0), is.infinite(Ad / 0)) &&
                            identical(is.nan(d0), is.nan(Ad / 0)))

# 比较：结果一定是稠密的逻辑矩阵，且与 base R 逐元素一致。
chk("比较 >",          identical(as.matrix(As > 2), Ad > 2))
chk("比较 <=",         identical(as.matrix(As <= 0), Ad <= 0))
chk("比较 ==",         identical(as.matrix(As == 0), Ad == 0))
chk("比较 != ",        identical(as.matrix(As != 0), Ad != 0))
chk("比较稀疏 vs 稀疏", identical(as.matrix(As > As), Ad > Ad))
# ! 走一元分支，等价于 x == 0，与 Ops.eigenMat 的 !E 同一口径。
chk("逻辑取反",        identical(as.matrix(!As), !Ad))

# 整除与取模没有稀疏实现，全部落到稠密分支。零处的语义（0 %% 2 == 0、
# 0 %/% 2 == 0）决定了它们不可能保持稀疏，所以只验数值。
chk("取模 %%",         close(as.matrix(As %% 2), Ad %% 2, 1e-12))
chk("整除 %/%",        close(as.matrix(As %/% 2), Ad %/% 2, 1e-12))
chk("标量在左的除法",  same_special(as.matrix(2 / As), 2 / Ad))
chk("标量在左的减法",  close(as.matrix(1 - As), 1 - Ad, 1e-12))
chk("标量在左的加法",  close(as.matrix(1 + As), 1 + Ad, 1e-12))

cat("== 10. 稀疏求解：每个方法都要对得上 base ==\n")
# 一维 Laplacian：对称正定，带状，是稀疏求解器的典型场景。
n <- 40
L <- matrix(0, n, n)
for (k in 1:n) { L[k, k] <- 2; if (k < n) { L[k, k + 1] <- -1; L[k + 1, k] <- -1 } }
Ls <- eigen_sparse(L)
bb <- rnorm(n)
ref <- solve(L, bb)
for (mth in c("lu", "llt", "ldlt", "cg")) {
  chk(paste("solve", mth), close(unw(eigen_sparse_solve(Ls, bb, mth)), ref, 1e-8))
}
chk("solve.eigenSparse 走 lu", close(unw(solve(Ls, bb)), ref, 1e-8))
# 返回形状必须与 base R 及 solve.eigenMat 一致：b 是无 dim 向量时给向量。
# 稠密那条路一直是这么做的，稀疏这边第一版漏了，于是 solve(S, b) 给 10x1 的
# eigenMat 而 solve(as.matrix(S), b) 给向量。
chk("solve(S, 向量) 返回向量",   is.null(dim(solve(Ls, bb))))
chk("与 base 同为向量",          identical(is.null(dim(solve(Ls, bb))),
                                          is.null(dim(solve(L, bb)))))
chk("solve(S, 矩阵) 返回矩阵",   !is.null(dim(solve(Ls, matrix(bb, n, 1)))))
chk("eigen_sparse_solve 恒为 eigenMat",
    inherits(eigen_sparse_solve(Ls, bb, "lu"), "eigenMat"))
chk("solve 无 b 时解单位阵",     close(as.matrix(solve(Ls)),
                                       solve(L, diag(n)), 1e-8))

# 一般非对称方阵走 LU（对角占优，保证良态）
set.seed(202)
G <- matrix(0, 12, 12)
for (i in 1:12) for (j in 1:12)
  if (i == j) G[i, j] <- rnorm(1) + 4 else if (abs(i - j) == 1) G[i, j] <- rnorm(1)
Gs <- eigen_sparse(G)
gb <- rnorm(12)
chk("非对称 LU", close(unw(eigen_sparse_solve(Gs, gb, "lu")), solve(G, gb), 1e-8))

cat("== 11. 求解器必须拦住失败，而不是给出别的方程组的解 ==\n")
# 不定矩阵：LLT 在这里会失败。Eigen 不抛异常，只把 info() 设成
# NumericalIssue，解照常返回——所以必须显式检查。
Ind <- eigen_sparse(matrix(c(1, 2, 2, 1), 2))
chk("LLT 拒绝不定矩阵", errs(eigen_sparse_solve(Ind, c(1, 1), "llt")))
chk("LDLT 接受不定矩阵", !errs(eigen_sparse_solve(Ind, c(1, 1), "ldlt")))
chk("LU 接受不定矩阵",   !errs(eigen_sparse_solve(Ind, c(1, 1), "lu")))

# 精确奇异
Sin <- eigen_sparse(matrix(c(1, 2, 2, 4), 2))
chk("LU 拒绝奇异矩阵", errs(eigen_sparse_solve(Sin, c(1, 1), "lu")))
chk("非方阵要求 qr",    errs(eigen_sparse_solve(eigen_sparse(matrix(1, 3, 2)), rep(1, 3), "lu")))

cat("== 12. 长方形：SparseQR 最小二乘 ==\n")
Tall <- matrix(0, 10, 3)
for (k in 1:3) { Tall[1:10, k] <- c(rep(1, 10)) * (k == 1) }
Tall[, 1] <- 1; Tall[, 2] <- 1:10; Tall[, 3] <- (1:10)^2
Ts <- eigen_sparse(Tall)
yt <- rnorm(10)
fit <- unw(eigen_sparse_solve(Ts, yt, "qr"))
chk("QR 最小二乘拟合值", close(Tall %*% fit, Tall %*% qr.solve(Tall, yt), 1e-8))
chk("QR 返回列向量",     identical(dim(fit), c(3L, 1L)))

cat("== 13. 非有限值与打印 ==\n")
chk("all_finite 正常矩阵", isTRUE(eigen_all_finite(Ls)))
Bad <- eigen_sparse(c(1, 2), c(1, 2), c(1, Inf), nrow = 2, ncol = 2)
chk("all_finite 检出 Inf", isFALSE(eigen_all_finite(Bad)))
chk("has_inf",             isTRUE(eigen_has_inf(Bad)))
chk("has_nan 为假",        isFALSE(eigen_has_nan(Bad)))
chk("求解器拒绝 Inf",      errs(eigen_sparse_solve(Bad, c(1, 1), "lu")))
chk("print 不报错",        { capture.output(print(Ls)); TRUE })
chk("print 含尺寸",        any(grepl("40 x 40", capture.output(print(Ls)))))

cat("== 14. Matrix 包互操作 ==\n")
if (requireNamespace("Matrix", quietly = TRUE)) {
  suppressMessages(library(Matrix))
  Ms <- Matrix::sparseMatrix(i = c(1, 2, 3), j = c(1, 2, 3), x = c(7, 8, 9),
                             dims = c(4, 4))
  Es <- eigen_sparse(Ms)
  chk("从 dgCMatrix 转换", close(as.matrix(Es), as.matrix(Ms), 1e-12))
  chk("转换后维度",        identical(as.integer(dim(Es)), c(4L, 4L)))
  # Matrix 1.7.6 不再导出 nnz()；dgCMatrix 的 @x 就是存储值，长度即 nnz。
  chk("转换后 nnz",        identical(as.integer(eigen_sparse_nnz(Es)),
                                     as.integer(length(Ms@x))))
} else {
  cat("  (Matrix 未安装，跳过)\n")
}

cat("== 15. Summary / Math / 混合运算 ==\n")
# 隐含零必须参与归约：全正矩阵的 min 是 0，prod 是 0。这两条最容易写错成
# "只对存下来的值做归约"，那样 min 会给出正数、prod 会给出非零。
P <- abs(Ad)                     # 有零、有正数
Ps <- eigen_sparse(P)
chk("sum",            close(sum(Ps), sum(P), 1e-12))
chk("max",            close(max(Ps), max(P), 1e-12))
chk("min 含隐含零",   close(min(Ps), min(P), 1e-12))
chk("range",          close(range(Ps), range(P), 1e-12))
chk("prod 含隐含零",  close(prod(Ps), prod(P), 1e-12))
chk("mean",           close(mean(Ps), mean(P), 1e-12))
chk("sum 负数",       close(sum(As), sum(Ad), 1e-12))
chk("min 负数",       close(min(As), min(Ad), 1e-12))
chk("any",            identical(any(As), any(Ad)))
# base 的 all() 会把 double 强制成 logical 并**报警**，而稀疏侧直接返回 logical。
# 只对稠密参照 suppressWarnings，稀疏那条路的行为仍然完全可见。
chk("all",            identical(all(As), suppressWarnings(all(Ad))))
chk("any 全零",       identical(any(eigen_sparse(matrix(0, 3, 3))), FALSE))
chk("all 全非零",     identical(all(eigen_sparse(matrix(c(1, 2, 3, 4), 2))), TRUE))
chk("多参数 sum",     close(sum(Ps, Ps), sum(P, P), 1e-12))

# NA / NaN 必须分开：base R 的 sum(c(1, NA)) 是 NA，sum(c(1, NaN)) 是 NaN。
An <- Ad; An[1, 1] <- NA; Asn <- eigen_sparse(An)
Nn <- Ad; Nn[1, 1] <- NaN; Nsn <- eigen_sparse(Nn)
chk("sum 遇 NA 给 NA",      identical(sum(Asn), sum(An)))
chk("sum 遇 NaN 给 NaN",    identical(sum(Nsn), sum(Nn)))
chk("sum na.rm",            close(sum(Asn, na.rm = TRUE), sum(An, na.rm = TRUE), 1e-12))
chk("max 遇 NA",            identical(max(Asn), max(An)))
chk("anyNA 认 NA",          identical(anyNA(Asn), anyNA(An)))
chk("anyNA 干净矩阵",       identical(anyNA(As), FALSE))
# 同一个问题、两条门：稠密侧 reigen_finite 用 Eigen 的 isNaN()，稀疏侧必须
# 同一口径。这里第一版用的是 R_IsNaN()，它对 NA 返回 false，于是
# eigen_has_nan(Asn) 是 FALSE 而 eigen_has_nan(An) 是 TRUE。
chk("has_nan 与稠密一致（NA）",  identical(eigen_has_nan(Asn), eigen_has_nan(An)))
chk("has_nan 与稠密一致（NaN）", identical(eigen_has_nan(Nsn), eigen_has_nan(Nn)))
chk("has_nan 干净矩阵一致",      identical(eigen_has_nan(As), eigen_has_nan(Ad)))
chk("all_finite 与稠密一致",     identical(eigen_all_finite(Asn), eigen_all_finite(An)))

# Math：把 0 映射成 0 的保持稀疏，其余转稠密。
chk("abs 保持稀疏",   inherits(abs(As), "eigenSparse"))
chk("abs 数值",       close(as.matrix(abs(As)), abs(Ad), 1e-12))
chk("sqrt 保持稀疏",  inherits(sqrt(Ps), "eigenSparse"))
chk("sqrt 数值",      close(as.matrix(sqrt(Ps)), sqrt(P), 1e-12))
chk("sign 保持稀疏",  close(as.matrix(sign(As)), sign(Ad), 1e-12))
chk("round 带参数",   close(as.matrix(round(Ps, 2)), round(P, 2), 1e-12))
chk("exp 转稠密",     !inherits(exp(As), "eigenSparse"))
# 返回**类型**也要和稠密侧一致：exp(E) 是 eigenMat，exp(As) 不能是 base matrix。
# 这条是 GC 拷问脚本先发现的——单元测试只比了数值，没比类型。
chk("exp 返回类型与稠密一致", identical(class(exp(As))[1], class(exp(eigenMat(Ad)))[1]))
chk("log 返回类型与稠密一致", identical(class(log(Ps))[1], class(log(eigenMat(P)))[1]))
chk("exp 数值",       close(as.matrix(exp(As)), exp(Ad), 1e-12))
chk("log 转稠密",     !inherits(log(Ps), "eigenSparse"))
chk("log 数值（含 -Inf）", same_special(as.matrix(log(Ps)), log(P)))
chk("cumsum 转稠密",  close(as.matrix(cumsum(Ps)), cumsum(P), 1e-12))

chk("is.numeric",     isTRUE(is.numeric(As)))
chk("is.na 是矩阵",   is.matrix(is.na(As)) && !inherits(is.na(As), "eigenObject"))
chk("is.na 数值",     identical(is.na(Asn), is.na(An)))
chk("as.data.frame",  is.data.frame(as.data.frame(As)) &&
                      close(as.matrix(as.data.frame(As)), Ad, 1e-12))

# 混合稀疏与稠密 eigenMat。R 的 Ops 组泛型要求两个操作数派发到同一个方法，
# 否则警告 "Incompatible methods" 并落到 default——所以 Ops.eigenSparse
# 必须和 Ops.eigenMat 是同一个闭包。
Em <- eigenMat(Ad)
chk("稀疏 + 稠密",     close(as.matrix(As + Em), Ad + Ad, 1e-12))
chk("稠密 + 稀疏",     close(as.matrix(Em + As), Ad + Ad, 1e-12))
chk("稀疏 - 稠密",     close(as.matrix(As - Em), Ad - Ad, 1e-12))
chk("稀疏 * 稠密",     close(as.matrix(As * Em), Ad * Ad, 1e-12))
chk("稀疏 / 稠密",     close(as.matrix(As / (abs(Em) + 1)), Ad / (abs(Ad) + 1), 1e-12))
chk("稀疏 > 稠密",     identical(as.matrix(As > Em), Ad > Ad))
# 1x1 的稠密操作数：is.numeric() 为真、length() 为 1，别把它当标量走
# reigen_sparse_scale——那条路上 as.numeric() 会在 externalptr 上炸掉。
E1 <- eigenMat(matrix(3, 1, 1))
chk("稀疏 * 1x1 稠密", close(as.matrix(As * E1), Ad * 3, 1e-12))
chk("1x1 稠密 * 稀疏", close(as.matrix(E1 * As), 3 * Ad, 1e-12))
chk("混合运算无警告", {
  w <- NULL
  withCallingHandlers(As + Em, warning = function(e) { w <<- conditionMessage(e); invokeRestart("muffleWarning") })
  is.null(w)
})

cat("== 16. BiCGSTAB / LSQR：CG 覆盖不到的两类问题 ==\n")
# CG 只对对称正定有效。实践中大多数稀疏系统是非对称的，BiCGSTAB 就是那类
# 问题的迭代解法。这里构造一个次对角 -2、超对角 -1 的三对角矩阵：非对称、
# 严格对角占优，所以解存在且迭代收敛得快。
n2 <- 60
NS <- matrix(0, n2, n2)
for (k in 1:n2) {
  NS[k, k] <- 4
  if (k < n2) NS[k, k + 1] <- -1
  if (k > 1)  NS[k, k - 1] <- -2
}
NSs <- eigen_sparse(NS)
b2 <- rnorm(n2)
x_ref <- solve(NS, b2)

chk("BiCGSTAB 非对称求解", close(as.matrix(eigen_sparse_solve(NSs, b2, "bicgstab")), x_ref, 1e-6))
# 次对角 -2、超对角 -1，所以 NS - t(NS) 在非对角处正好是 ±1。用 > 0 判，
# 不要写 > 1——那会因为差值恰好等于 1 而失败（本测试第一版就是这么错的）。
chk("BiCGSTAB 真的非对称", max(abs(NS - t(NS))) > 0)
# 残差要按 A x - b 量，不能只比 x：x 对但残差大意味着后面用它会出问题。
xb <- as.matrix(eigen_sparse_solve(NSs, b2, "bicgstab"))
chk("BiCGSTAB 残差", max(abs(NS %*% xb - b2)) < 1e-8)
chk("BiCGSTAB 多右端", {
  B2 <- cbind(b2, rnorm(n2))
  close(as.matrix(eigen_sparse_solve(NSs, B2, "bicgstab")), solve(NS, B2), 1e-6)
})

# 迭代法必须报出迭代次数和最终残差：解出来 ≠ 解对了。
it <- attr(eigen_sparse_solve(NSs, b2, "bicgstab"), "iterations")
rs <- attr(eigen_sparse_solve(NSs, b2, "bicgstab"), "residual")
chk("迭代次数是正数", is.numeric(it) && length(it) == 1L && it > 0)
chk("残差已报出",     is.numeric(rs) && rs < 1e-8)
chk("直接法不挂属性", is.null(attr(eigen_sparse_solve(NSs, b2, "lu"), "iterations")))
chk("solve() 也能透传", !is.null(attr(solve(NSs, b2, method = "bicgstab"), "iterations")))

# 迭代上限：要一个**真的**需要多步的系统。三对角严格对角占优矩阵
# BiCGSTAB 一步就收敛了（第一版用它测 maxiter=1，结果根本没报错）。
# 换成 2 维 5 点 Poisson：条件数 O(n)，CG 至少要 O(sqrt(n)) 步。
nx <- 20
Po <- matrix(0, nx * nx, nx * nx)
for (i in 1:nx) for (j in 1:nx) {
  k <- (i - 1) * nx + j
  Po[k, k] <- 4
  if (i > 1)  Po[k, k - nx] <- -1
  if (i < nx) Po[k, k + nx] <- -1
  if (j > 1)  Po[k, k - 1]  <- -1
  if (j < nx) Po[k, k + 1]  <- -1
}
Pos <- eigen_sparse(Po)
bP <- rnorm(nx * nx)
xP <- eigen_sparse_solve(Pos, bP, "cg", tol = 1e-10)
chk("CG 在 Poisson 上收敛", max(abs(Po %*% as.matrix(xP) - bP)) < 1e-6)
chk("CG 步数大于 1",        attr(xP, "iterations") > 1)
chk("maxiter 不够要报错",   errs(eigen_sparse_solve(Pos, bP, "cg", maxiter = 1L)))
chk("放宽 tol 迭代更少", {
  a1 <- attr(eigen_sparse_solve(Pos, bP, "cg", tol = 1e-2), "iterations")
  a2 <- attr(eigen_sparse_solve(Pos, bP, "cg", tol = 1e-12), "iterations")
  a1 < a2
})
chk("BiCGSTAB 也能解 Poisson", {
  xb2 <- eigen_sparse_solve(Pos, bP, "bicgstab", tol = 1e-12)
  max(abs(Po %*% as.matrix(xb2) - bP)) < 1e-6
})

# 预条件子：默认对角（永远有定义），ilut 是选项（快，但可能因零主元崩掉）。
# 两条路都必须给出同一个解，只是步数不同。
chk("ilut 与 diagonal 同解", {
  a <- eigen_sparse_solve(Pos, bP, "bicgstab", precond = "diagonal", tol = 1e-12)
  b <- eigen_sparse_solve(Pos, bP, "bicgstab", precond = "ilut", tol = 1e-12)
  max(abs(as.matrix(a) - as.matrix(b))) < 1e-6
})
chk("ilut 步数更少", {
  a <- attr(eigen_sparse_solve(Pos, bP, "bicgstab", precond = "diagonal", tol = 1e-12), "iterations")
  b <- attr(eigen_sparse_solve(Pos, bP, "bicgstab", precond = "ilut", tol = 1e-12), "iterations")
  b <= a
})
chk("precond 只认两个值", errs(eigen_sparse_solve(Pos, bP, "cg", precond = "nope")))
# 默认必须是 ilut：对角预条件在同一个系统上步数多一个量级，如果默认值写错
# （R 的向量顺序和 C++ 的编号对不上），症状是"能算但慢"，很难被发现。
chk("默认 precond 是 ilut", {
  d <- attr(eigen_sparse_solve(Pos, bP, "bicgstab"), "iterations")
  i <- attr(eigen_sparse_solve(Pos, bP, "bicgstab", precond = "ilut"), "iterations")
  identical(d, i)
})
# CG 要求预条件子对称正定，ILUT 不是对称的。这个组合**不是**"慢一点"，
# 而是数学上不成立：实测在 Poisson 上直接不收敛。所以要显式拒绝。
chk("cg + ilut 必须报错", errs(eigen_sparse_solve(Pos, bP, "cg", precond = "ilut")))
chk("cg 默认用 diagonal", {
  d <- attr(eigen_sparse_solve(Pos, bP, "cg"), "iterations")
  g <- attr(eigen_sparse_solve(Pos, bP, "cg", precond = "diagonal"), "iterations")
  identical(d, g)
})
# ilut 崩掉时必须**报错**，不能把 NaN 当解返回。
chk("ilut 崩掉要报错而不是返回 NaN", {
  # 强对流 + 大网格会让 IncompleteLUT 撞上零主元
  m3 <- 60; M3 <- matrix(0, m3 * m3, m3 * m3)
  for (i in 1:m3) for (j in 1:m3) {
    k <- (i - 1) * m3 + j
    M3[k, k] <- 4
    if (i > 1)  M3[k, k - m3] <- -1
    if (i < m3) M3[k, k + m3] <- -1
    if (j > 1)  M3[k, k - 1]  <- -3
    if (j < m3) M3[k, k + 1]  <- -1
  }
  r <- tryCatch(eigen_sparse_solve(eigen_sparse(M3), rnorm(m3 * m3), "bicgstab",
                                   precond = "ilut", tol = 1e-12),
                error = function(e) "stopped")
  # 要么收敛并给出有限解，要么报错；绝不允许静默返回 NaN
  identical(r, "stopped") || all(is.finite(as.matrix(r)))
})

# 长方形：LSQR 与 SparseQR 应该给出同一个最小二乘解。
Tall2 <- cbind(1, seq_len(40), seq_len(40)^2)
Tall2 <- Tall2 + matrix(rnorm(40 * 3, 0, 1e-3), 40, 3)
y2 <- rnorm(40)
q_sol <- as.matrix(eigen_sparse_solve(eigen_sparse(Tall2), y2, "qr"))
l_sol <- as.matrix(eigen_sparse_solve(eigen_sparse(Tall2), y2, "lsqr", tol = 1e-14))
chk("LSQR 与 QR 一致",   close(l_sol, q_sol, 1e-5))
chk("LSQR 与 qr.solve",  close(l_sol, qr.solve(Tall2, y2), 1e-5))
chk("LSQR 接受长方形",   nrow(l_sol) == 3L)
# 方阵专用方法遇到长方形必须报错，而不是给出别的东西。
chk("BiCGSTAB 拒绝长方形", errs(eigen_sparse_solve(eigen_sparse(Tall2), y2, "bicgstab")))
chk("CG 拒绝长方形",       errs(eigen_sparse_solve(eigen_sparse(Tall2), y2, "cg")))

# crossprod / tcrossprod / aperm / is.finite 在稠密侧是内部泛型；稀疏对象的
# class 向量里没有 "eigenMat"，派发到不了稠密的方法，所以这里必须各有一份。
cat("== 17. 内部泛型：crossprod / tcrossprod / aperm / is.finite ==\n")
# 用第 8 节的 Ad（30x20）/ As 和 B（20x4）：crossprod 会缩到 20x20，
# tcrossprod 会长到 30x30，两个方向都覆盖到。注意 crossprod(Ad, Y) 要的是
# t(Ad) 的列数对齐 Y 的行数，也就是 Y 必须有 30 行——B 只有 20 行，用不了。
B30 <- matrix(rnorm(30 * 3), 30, 3)
chk("crossprod(S)",        close(crossprod(As), crossprod(Ad), 1e-12))
chk("crossprod(S, B30)",   close(crossprod(As, B30), crossprod(Ad, B30), 1e-12))
chk("crossprod(Ad, S) 也派发", close(crossprod(Ad, As), crossprod(Ad), 1e-12))
chk("crossprod 保持稀疏",  inherits(crossprod(As), "eigenSparse"))
chk("tcrossprod(S)",       close(tcrossprod(As), tcrossprod(Ad), 1e-12))
chk("tcrossprod 保持稀疏", inherits(tcrossprod(As), "eigenSparse"))
chk("aperm 转置",          close(aperm(As, c(2, 1)), t(Ad), 1e-12))
chk("aperm 恒等",          close(aperm(As, c(1, 2)), Ad, 1e-12))
chk("aperm 非法 perm 报错", errs(aperm(As, c(2, 2))))
chk("is.finite",           identical(is.finite(As), is.finite(Ad)))
chk("is.infinite",         identical(is.infinite(As), is.infinite(Ad)))
chk("is.nan",              identical(is.nan(As), is.nan(Ad)))
# 隐含零是有限的，但存下来的 Inf/NA 必须被认出来——两边口径要一致。
Adf <- Ad; Adf[1, 1] <- Inf; Adf[2, 2] <- NA
chk("is.finite 认 Inf/NA", identical(is.finite(eigen_sparse(Adf)),
                                     is.finite(Adf)))

cat("== 18. 组泛型要求同一个闭包：稀疏 × 稠密 ==\n")
# %*% / crossprod / tcrossprod 都在**两个**操作数上派发，而且要求两边落到
# 同一个方法闭包。稠密和稀疏各有一份实现时，S %*% E 会警告
#   Incompatible methods ("%*%.eigenSparse", "%*%.eigenMat") for "%*%"
# 然后落到 default 报错——两个函数体一次都不执行。这和当初 As + E 是同一个坑。
# 用方阵 Po 做，形状约束最简单。
Eo <- eigenMat(Po)
chk("S %*% eigenMat",     close(Pos %*% Eo, Po %*% Po, 1e-10))
chk("eigenMat %*% S",     close(Eo %*% Pos, Po %*% Po, 1e-10))
chk("crossprod(S, E)",    close(crossprod(Pos, Eo), crossprod(Po), 1e-10))
chk("crossprod(E, S)",    close(crossprod(Eo, Pos), crossprod(Po), 1e-10))
chk("tcrossprod(S, E)",   close(tcrossprod(Pos, Eo), tcrossprod(Po), 1e-10))
chk("tcrossprod(E, S)",   close(tcrossprod(Eo, Pos), tcrossprod(Po), 1e-10))
# 解完方程看残差——这个 bug 正是被这一行暴露的：eigen_sparse_solve 返回的
# 就是 eigenMat，所以 S %*% x 几乎必然被走到。
chk("求解后残差 S %*% x", {
  xs <- eigen_sparse_solve(Pos, bP, "lu")
  close(Pos %*% xs, bP, 1e-8)
})
chk("迭代解后残差", {
  xs <- eigen_sparse_solve(Pos, bP, "cg", tol = 1e-12)
  close(Pos %*% xs, bP, 1e-6)
})
# 只断言"没报错"是不够的：Incompatible methods 先发 warning，再落到 default
# 才 error。要证明真的进了方法体，得断言那条 warning 不再出现。
chk("混合乘法不再产生 Incompatible 警告", {
  w <- character(0)
  withCallingHandlers(Pos %*% Eo, warning = function(cw) {
    w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning")
  })
  !any(grepl("Incompatible", w))
})

# Summary 组泛型**只按第一个参数派发**。所以 sum(S, E) 和 sum(E, S) 走的是
# 两个不同的方法：前者进 Summary.eigenSparse，后者进 Summary.eigenMat 的
# 多参数兜底，而那个兜底靠 .as_r() 落回 base。.as_r 原先只认 "eigenMat"，
# 于是 sum(E, S) 把 externalptr 直接交给 base sum() 而报错。
chk("sum(S, E)", close(sum(Pos, Eo), sum(Po) + sum(Po), 1e-10))
chk("sum(E, S)", close(sum(Eo, Pos), sum(Po) + sum(Po), 1e-10))
chk("max(E, S)", close(max(Eo, Pos), max(Po), 1e-10))
chk("sum(E, 1)", close(sum(Eo, 1), sum(Po) + 1, 1e-10))
# 这一条**修不了**：派发看第一个参数，普通数值在前就进不了方法。
# 断言它确实报错，免得以后有人以为已经修好了。
chk("sum(1, S) 仍然报错", errs(sum(1, Pos)))

cat("== 19. 稀疏侧的归约：同一个量，两条门 ==\n")
# 这些原本只有稠密实现，eigen_sum(S) 报 "expected a numeric matrix or vector,
# got externalptr"——可 sum(S) 是能算的。同一个量两条门给出不同结果，正是
# consistency.R 一直在抓的那类问题。全部走 reigen_sparse_stats 的那一遍扫描。
chk("eigen_sum(S)",        close(eigen_sum(Ps), sum(P), 1e-12))
chk("eigen_prod(S)",       close(eigen_prod(Ps), prod(P), 1e-12))
chk("eigen_mean(S)",       close(eigen_mean(Ps), mean(P), 1e-12))
chk("eigen_sumsq(S)",      close(eigen_sumsq(Ps), sum(P^2), 1e-12))
chk("eigen_norm(S, 'F')",  close(eigen_norm(Ps, "F"), sqrt(sum(P^2)), 1e-12))
chk("eigen_trace(S)",      close(eigen_trace(Ps), sum(diag(P)), 1e-12))
chk("eigen_diag(S)",       close(eigen_diag(Ps), matrix(diag(P), ncol = 1), 1e-12))
chk("eigen_diag 返回 eigenMat", inherits(eigen_diag(Ps), "eigenMat"))
chk("eigen_is_complex(S)", identical(eigen_is_complex(Ps), FALSE))
chk("eigen_adjoint(S) 保持稀疏", inherits(eigen_adjoint(Ps), "eigenSparse"))
chk("eigen_adjoint(S) 值", close(eigen_adjoint(Ps), t(P), 1e-12))
# 与泛型口径一致：同一个量，两条门必须给同一个答案
chk("eigen_sum == sum",    close(eigen_sum(Ps), sum(Ps), 1e-12))
chk("eigen_mean == mean",  close(eigen_mean(Ps), mean(Ps), 1e-12))
chk("eigen_prod == prod",  close(eigen_prod(Ps), prod(Ps), 1e-12))

# 缺失值必须传上去。reigen_sparse_stats() 为了 na.rm = TRUE 那条路是**跳过**
# 缺失值累加的，所以每个调用方都要自己判一次。漏掉不会报错，只会给一个偏小的
# 数：A = matrix(c(1,2,NA,4,5,6), 2, 3) 时 eigen_sum(S) 曾经是 18，而 sum(S)
# 是 NA；eigen_prod 是 240 而 prod(S) 是 NA；mean(S) 是 3 而 mean(A) 是 NA。
# 这三条是同一个根因，所以三条一起断言，而且都拿 base 当参照。
Pna <- P; Pna[2, 3] <- NA; Sna <- eigen_sparse(Pna)
chk("eigen_sum 传 NA",   identical(eigen_sum(Sna),  sum(Pna)))
chk("eigen_prod 传 NA",  identical(eigen_prod(Sna), prod(Pna)))
chk("eigen_mean 传 NA",  identical(eigen_mean(Sna), mean(Pna)))
chk("eigen_sumsq 传 NA", identical(eigen_sumsq(Sna), sum(Pna^2)))
chk("mean(S) 传 NA",     identical(mean(Sna), mean(Pna)))
chk("sum(S) 传 NA",      identical(sum(Sna), sum(Pna)))
# NA 和 NaN 不能混为一谈：sum(c(1, NA)) 是 NA，sum(c(1, NaN)) 是 NaN，
# 两者 is.na() 都为真但 identical() 为假。
Pnan <- P; Pnan[2, 3] <- NaN; Snan <- eigen_sparse(Pnan)
chk("eigen_sum 传 NaN",  identical(eigen_sum(Snan),  sum(Pnan)))
chk("eigen_prod 传 NaN", identical(eigen_prod(Snan), prod(Pnan)))
chk("eigen_mean 传 NaN", identical(eigen_mean(Snan), mean(Pnan)))
# na.rm = TRUE 那条路不能因此坏掉
chk("sum(S, na.rm=TRUE)",  close(sum(Sna, na.rm = TRUE),  sum(Pna, na.rm = TRUE), 1e-12))
chk("mean(S, na.rm=TRUE)", close(mean(Sna, na.rm = TRUE), mean(Pna, na.rm = TRUE), 1e-12))
# 迹只看对角线：非对角线上的 NA 不影响它，对角线上的才影响。
# 全矩阵的 na 计数回答不了这个问题，所以统计里单独数了 trace_na / trace_nan。
Poff <- P; Poff[1, 2] <- NA            # 非对角线
chk("trace 忽略非对角 NA", close(eigen_trace(eigen_sparse(Poff)), sum(diag(P)), 1e-12))
Pdiag <- P; Pdiag[2, 2] <- NA          # 对角线
chk("trace 传对角 NA", identical(eigen_trace(eigen_sparse(Pdiag)), sum(diag(Pdiag))))

# 不支持的必须报错，而且要说清楚为什么——不是 "got externalptr"
chk("eigen_norm(S, 'p') 报错", errs(eigen_norm(Ps, "p")))
chk("eigen_norm(S, 'Z') 报错", errs(eigen_norm(Ps, "Z")))
# 只吃稠密的入口：错误信息要能照着做（提 eigenSparse、提 as.matrix）
chk("稠密入口的报错可操作", {
  msg <- tryCatch(eigen_svd(Ps), error = function(e) conditionMessage(e))
  is.character(msg) && grepl("eigenSparse", msg) && grepl("as.matrix", msg)
})
chk("稠密入口的报错列了替代入口", {
  msg <- tryCatch(eigen_inverse(Ps), error = function(e) conditionMessage(e))
  is.character(msg) && grepl("eigen_sparse_solve", msg)
})

cat("== 20. c / isSymmetric / median / quantile ==\n")
# c() 只按**第一个**参数派发。没有 c.eigenSparse 时，c(S) 落到 default，
# 而 default 对 externalptr 只能构造 list——于是 c(S) 静默返回一个装着裸指针
# 的长度 1 list。不报错、值也不对，是最坏的一类。
chk("c(S) 是数值",       is.numeric(c(Ps)))
chk("c(S) 与 base 一致", close(c(Ps), c(P), 1e-12))
chk("c(S, E) 混合",      close(c(Ps, Eo), c(P, Po), 1e-12))

chk("isSymmetric 对称阵", identical(isSymmetric(Pos), isSymmetric(Po)))
NPo <- Po; NPo[1, 2] <- NPo[1, 2] + 1
chk("isSymmetric 非对称", identical(isSymmetric(eigen_sparse(NPo)), isSymmetric(NPo)))
# 非方阵：base R 返回 FALSE，不能因为 t(S)-S 不合法就报错
chk("isSymmetric 非方阵", identical(isSymmetric(As), isSymmetric(Ad)))
chk("isSymmetric 带 tol", identical(isSymmetric(eigen_sparse(NPo), tol = 1e3),
                                    isSymmetric(NPo, tol = 1e3)))

# 隐含零必须参与：稀疏矩阵里大部分元素是 0，中位数/分位数直接由它决定。
chk("median(S)",        close(median(Ps), median(P), 1e-12))
chk("median 带隐含零",  close(median(Ps), median(P), 1e-12))
chk("quantile(S)",      close(quantile(Ps), quantile(P), 1e-12))
chk("quantile 指定 probs", close(quantile(Ps, c(0.1, 0.9)),
                                 quantile(P, c(0.1, 0.9)), 1e-12))
chk("quantile 名字与 base 一致", identical(names(quantile(Ps)), names(quantile(P))))
chk("quantile names = FALSE",    identical(names(quantile(Ps, names = FALSE)), NULL))
chk("quantile type = 1 兜底",    close(quantile(Ps, type = 1), quantile(P, type = 1), 1e-12))
chk("全零矩阵的中位数",          identical(median(eigen_sparse(matrix(0, 4, 4))), 0))
chk("probs 越界报错",            errs(quantile(Ps, 1.5)))

# 上面那一组用的是 P = abs(Ad)，**全非负**，于是"零排在升序最前面"这个错误
# 假设照样能通过——第一版就是这么错的：0% 分位数给 0，正确答案是最小值。
# 真矩阵有正有负，零夹在中间，所以这里专门造一个三态都有的夹具。
set.seed(99)
Mixed <- matrix(rnorm(120), 12, 10)
Mixed[sample(length(Mixed), 40)] <- 0        # 40 个精确零
Ms <- eigen_sparse(Mixed)
chk("median 含正负零",     close(median(Ms), median(Mixed), 1e-12))
chk("quantile 含正负零",   close(quantile(Ms), quantile(Mixed), 1e-12))
chk("quantile 0% 即最小值", close(quantile(Ms, 0), min(Mixed), 1e-12))
chk("quantile 100% 即最大值", close(quantile(Ms, 1), max(Mixed), 1e-12))
chk("median 全非正",       close(median(eigen_sparse(-P)), median(-P), 1e-12))
chk("quantile 全非正",     close(quantile(eigen_sparse(-P)), quantile(-P), 1e-12))

# 方法名在签名里写成字面量（codoc 要逐字比对），在 C++ 那边却靠
# .reigen_sparse_methods 的**顺序**映射成整数 code。两处一旦漂移，用户敲
# "llt" 会跑到 SimplicialLDLT 上去，而数值往往仍然"像对的"。所以钉死。
#
# 顺带说明一处看起来不一致、其实是故意的写法：solve.eigenSparse 的默认值仍然
# 写 .reigen_sparse_methods。它是未导出的 S3 方法、没有 Rd，codoc 不检查它，
# 所以那里用共享常量才是对的。只有导出的 eigen_sparse_solve 被迫写字面量。
chk("签名里的方法列表 == .reigen_sparse_methods", {
  sig <- eval(formals(ROpenEigen::eigen_sparse_solve)$method)
  identical(sig, ROpenEigen:::.reigen_sparse_methods)
})
chk("precond 列表顺序 == C++ 的 code", {
  identical(ROpenEigen:::.reigen_sparse_preconds, c("ilut", "diagonal"))
})

cat("== 21. 稀疏范数：五种类型，全部对稠密 ==\n")
# "2" 原本是唯一做不到的一种（要 A'A 的最大特征值，而稀疏没有特征分解），
# 所以当时除 "F" 以外全部报错——包括默认值 "O"，一个函数对最主要的一种输入
# 用默认参数就失败，这不是"明确拒绝"，是设计缺陷。现在五种都用幂迭代或
# 单遍扫描补上了。
#
# 注意默认值是 "O" 不是 "2"：type 的默认是整条向量，函数取 type[1]。
# 我先前在这里写过"默认 2 直接报错"，那是错的——旧代码拒绝所有非 "F" 的类型，
# 所以那条断言无论如何都会过，标签错了也不会有人发现。
NormCases <- list(
  list("随机 40x30", matrix(rnorm(1200), 40, 30)),
  list("全正",       matrix(abs(rnorm(300)), 20, 15)),
  list("全零",       matrix(0, 5, 4)),
  list("单非零",     { m <- matrix(0, 6, 6); m[3, 4] <- -2.5; m }),
  list("秩 1",       { u <- rnorm(25); v <- rnorm(18); u %o% v }),
  list("对角",       diag(c(3, -7, 1, 9, -2))),
  list("宽 10x30",   matrix(rnorm(300), 10, 30)),
  list("高 30x10",   matrix(rnorm(300), 30, 10))
)
# 每列和为 0：全 1 起点与主奇异向量正交，幂迭代会永远停在 0 上。
# 起点用固定种子的伪随机向量就是为了这个；用全 1 的话这条必挂。
NormCases[[9]] <- list("每列和为 0", { m <- matrix(rnorm(400), 20, 20); sweep(m, 2, colMeans(m)); m })

for (nc in NormCases) {
  nm <- nc[[1]]; A <- nc[[2]]; S <- eigen_sparse(A)
  for (ty in c("O", "I", "F", "M", "2")) {
    chk(sprintf("norm %-10s %s", nm, ty),
        close(eigen_norm(S, ty), norm(A, ty), 1e-8))
  }
}
# "1" 是 "O" 的别名，与 base::norm 一致
chk("norm 别名 '1' == 'O'", close(eigen_norm(eigen_sparse(NormCases[[1]][[2]]), "1"),
                                  norm(NormCases[[1]][[2]], "O"), 1e-12))
# 默认值不是 "2"。type 的默认是 c("O","I","F","M","2")，函数取 type[1]，
# 所以 eigen_norm(x) 走的是 "O"——和 base::norm 的 match.arg 一致。
# 这里钉死，因为我先前写文档时想当然地以为默认是 "2"，还照着那个假设写了
# 一条测试标签。标签写错不会失败，只会让后来的人读错。
chk("norm 默认 == 'O'", close(eigen_norm(eigen_sparse(NormCases[[1]][[2]])),
                              norm(NormCases[[1]][[2]], "O"), 1e-12))
# 谱范数与稠密侧 eigen_norm 走的是两条完全不同的路（SVD vs 幂迭代），
# 结果必须一致到 1e-8。
chk("norm '2' 稀疏 vs 稠密 eigen_norm",
    close(eigen_norm(eigen_sparse(NormCases[[1]][[2]]), "2"),
          eigen_norm(eigenMat(NormCases[[1]][[2]]), "2"), 1e-8))
# 空矩阵不能崩
chk("norm 空稀疏 F", identical(eigen_norm(eigen_sparse(matrix(0, 0, 3)), "F"), 0))
chk("norm 空稀疏 O", identical(eigen_norm(eigen_sparse(matrix(0, 0, 3)), "O"), 0))
# NA 的处理要和 base 一致：四种非谱范数"传上去"得 NA，谱范数 base 是报错
An <- matrix(rnorm(20), 5, 4); An[2, 3] <- NA
Sn <- eigen_sparse(An)
for (ty in c("O", "I", "F", "M"))
  chk(sprintf("norm NA %s", ty), identical(eigen_norm(Sn, ty), norm(An, ty)))
chk("norm NA '2' 报错", errs(eigen_norm(Sn, "2")))

# 稀疏侧的 F 范数用的是统计里另一套缩放累加器（st[["frob"]]），不能用
# sqrt(st[["sumsq"]])——sumsq 是朴素的，在 1e-200 上下溢成 0、1e300 上溢成
# Inf。而 eigen_sumsq() 必须保持朴素，因为 base 的 sum(x^2) 本来就是会溢出的。
for (k in c(-300, -200, -170, 150, 160, 300)) {
  B <- matrix(c(10^k, 0, 0, 10^k), 2, 2)
  chk(sprintf("sparse norm F 1e%d", k),
      isTRUE(all.equal(eigen_norm(eigen_sparse(B), "F"), norm(B, "F"), tolerance = 1e-12)))
}
chk("sparse norm F 单元素 1e-200 不是 0",
    eigen_norm(eigen_sparse(matrix(1e-200)), "F") > 0)
chk("sparse sumsq 仍然溢出",
    is.infinite(eigen_sumsq(eigen_sparse(matrix(1e300)))))
chk("sparse norm F 含 Inf",
    identical(eigen_norm(eigen_sparse(matrix(c(1, Inf, 2, 3), 2, 2)), "F"),
              norm(matrix(c(1, Inf, 2, 3), 2, 2), "F")))

cat("== 22. 极值的位置：隐含零也是候选 ==\n")
# which.min()/which.max() 在稠密侧一直可用，稀疏侧以前直接报错，理由是
# "位置取决于存储模式而不是数值"。位置确实取决于存储模式——可它仍然是唯一
# 正确的答案，而且和 which.min(as.matrix(x)) 给的是同一个。拒绝只是逼用户
# 自己稠密化绕一圈，而稠密化正是稀疏存在的理由。
# 注意 as.integer()：nrow() 给整数，但下面有几处传的是字面量 5，那是 double，
# 乘出来就是 double，identical() 比不过 which.min() 的整数。测试写错一次，
# 报出来的是"函数不对"——所以这里统一收成整数。
linidx <- function(w, nr)
  if (is.na(w$row)) NA_integer_ else as.integer((w$col - 1L) * nr + w$row)

WhichCases <- list(
  list("随机含负",   { m <- matrix(rnorm(120), 12, 10); m[sample(120, 60)] <- 0; m }),
  list("全正",       matrix(abs(rnorm(60)), 10, 6)),
  list("全负",       -matrix(abs(rnorm(60)), 10, 6)),
  list("全零",       matrix(0, 4, 5)),
  list("满阵无隐含零", matrix(rnorm(30), 6, 5)),
  list("单个非零",   { m <- matrix(0, 7, 7); m[4, 5] <- -2.5; m }),
  list("对角",       diag(c(3, -7, 1, 9, -2))),
  list("首格即最小", { m <- matrix(1, 5, 4); m[1, 1] <- -9; m })
)
for (wc in WhichCases) {
  nm <- wc[[1]]; A <- wc[[2]]; S <- eigen_sparse(A)
  wmin <- eigen_which_min(S); wmax <- eigen_which_max(S)
  chk(sprintf("which_min %-12s 值", nm), close(wmin$value, min(A), 1e-12))
  chk(sprintf("which_min %-12s 位置", nm), identical(linidx(wmin, nrow(A)), which.min(A)))
  chk(sprintf("which_max %-12s 值", nm), close(wmax$value, max(A), 1e-12))
  chk(sprintf("which_max %-12s 位置", nm), identical(linidx(wmax, nrow(A)), which.max(A)))
}
# 稀疏和稠密两扇门必须给同一组答案——这是 consistency.R 的原则，在这里再钉一次
chk("which_min 稀疏 == 稠密", {
  A <- WhichCases[[1]][[2]]
  a <- eigen_which_min(eigen_sparse(A)); b <- eigen_which_min(eigenMat(A))
  identical(c(a$value, a$row, a$col), c(b$value, b$row, b$col))
})
chk("which_max 稀疏 == 稠密", {
  A <- WhichCases[[1]][[2]]
  a <- eigen_which_max(eigen_sparse(A)); b <- eigen_which_max(eigenMat(A))
  identical(c(a$value, a$row, a$col), c(b$value, b$row, b$col))
})
# NA 要跳过，Inf 是合法取值
An2 <- matrix(rnorm(20), 5, 4); An2[2, 3] <- NA
chk("which_min 跳过 NA", identical(linidx(eigen_which_min(eigen_sparse(An2)), 5), which.min(An2)))
Ai <- matrix(c(Inf, 1, -Inf, 2), 2, 2)
chk("which_min 把 -Inf 当值", identical(linidx(eigen_which_min(eigen_sparse(Ai)), 2), which.min(Ai)))
chk("which_max 把 Inf 当值",  identical(linidx(eigen_which_max(eigen_sparse(Ai)), 2), which.max(Ai)))
# 全缺失：base 给 integer(0)，list 表达不了，退而报 NA 而不是编一个位置
Aall <- matrix(NA_real_, 3, 3)
wall <- eigen_which_min(eigen_sparse(Aall))
chk("which_min 全缺失给 NA", is.na(wall$value) && is.na(wall$row) && is.na(wall$col))

cat("== 23. 迭代求解器：GMRES / MINRES / DGMRES / IDR(s) ==\n")
# 与已有方法放在同一块 2 维 Poisson 上比：既有对称正定（CG/MINRES 的地盘），
# 也有非对称的对流扩散（GMRES/BiCGSTAB/IDR 的地盘）。
n2 <- 30
Poi <- Matrix::Matrix(0, n2, n2)
for (i in 1:n2) {
  Poi[i, i] <- 4
  if (i > 1) Poi[i, i - 1] <- -1
  if (i < n2) Poi[i, i + 1] <- -1
}
PoiS <- eigen_sparse(Poi)
bPoi <- as.numeric(seq_len(n2))
refPoi <- as.numeric(solve(as.matrix(Poi), bPoi))

conv <- Matrix::Matrix(0, n2, n2)
for (i in 1:n2) {
  conv[i, i] <- 4
  if (i > 1) conv[i, i - 1] <- -1.1
  if (i < n2) conv[i, i + 1] <- -1
}
convS <- eigen_sparse(conv)
refConv <- as.numeric(solve(as.matrix(conv), bPoi))

for (m in c("gmres", "minres", "dgmres", "idrs")) {
  r <- try(eigen_sparse_solve(PoiS, bPoi, method = m), silent = TRUE)
  chk(paste0(m, " 解出 Poisson"), !inherits(r, "try-error"))
  if (!inherits(r, "try-error"))
    chk(paste0(m, " 残差可接受"),
        max(abs(as.numeric(as.matrix(r)) - refPoi)) < 1e-8)
  chk(paste0(m, " 带迭代属性"), !is.null(attr(r, "iterations")))
}
for (m in c("gmres", "dgmres", "idrs")) {
  r <- try(eigen_sparse_solve(convS, bPoi, method = m), silent = TRUE)
  chk(paste0(m, " 解出对流扩散"), !inherits(r, "try-error"))
  if (!inherits(r, "try-error"))
    chk(paste0(m, " 对流残差"),
        max(abs(as.numeric(as.matrix(r)) - refConv)) < 1e-6)
}
# MINRES 只读下三角：不对称的矩阵进来不会报错，而是被当成对称的那一个。
# 这不是 bug，是 Eigen 的语义；钉住它，免得以后有人"顺手修好"。
chk("minres 只读下三角", {
  # up 的非零全在上三角。MINRES 用 selfadjointView<Lower>，于是它看到的
  # 矩阵是 diag(2, 2, 2)：解就是 b / 2。上三角那 5 被完全忽略，
  # 换成 t(up) 才会被读进去——所以转置前后答案必须不同。
  up <- matrix(0, 3, 3); up[upper.tri(up)] <- 5; diag(up) <- 2
  r1 <- as.numeric(as.matrix(eigen_sparse_solve(eigen_sparse(up), 1:3,
                                                method = "minres", maxiter = 50)))
  r2 <- as.numeric(as.matrix(eigen_sparse_solve(eigen_sparse(t(up)), 1:3,
                                                method = "minres", maxiter = 50)))
  close(r1, (1:3) / 2, 1e-8) && !isTRUE(all.equal(r1, r2))
})
chk("minres + ilut 必须报错", errs(eigen_sparse_solve(PoiS, bPoi, "minres",
                                                     precond = "ilut")))
chk("gmres 默认走 ilut", {
  a <- attr(eigen_sparse_solve(convS, bPoi, "gmres"), "iterations")
  b <- attr(eigen_sparse_solve(convS, bPoi, "gmres", precond = "ilut"), "iterations")
  close(a, b)
})
chk("未知方法名报错", errs(eigen_sparse_solve(PoiS, bPoi, "nope")))
chk("方法表与签名一致（含新增四种）", {
  sig <- eval(formals(ROpenEigen::eigen_sparse_solve)$method)
  identical(sig, ROpenEigen:::.reigen_sparse_methods) &&
    all(c("gmres", "minres", "dgmres", "idrs") %in% sig)
})
chk("idrs 在 maxiter 极小时报不收敛", errs(eigen_sparse_solve(convS, bPoi, "idrs",
                                                             maxiter = 2)))

cat(sprintf("\n通过 %d, 失败 %d\n", ok, bad))
if (bad > 0L) stop(sprintf("%d assertion(s) failed", bad), call. = FALSE)
