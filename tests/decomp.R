# Decomposition contracts.
#
# fuzz.R checks the matrix-level behaviour (arithmetic, indexing, svd/pinv/det)
# but only touches two of the nine decompositions. The other seven were covered
# by hand-picked examples in core.R, which proves they work on the example and
# says nothing about the contract.
#
# Each decomposition is defined by an identity. If the identity holds on random
# input of assorted shapes and conditioning, the factorisation is right; if it
# does not, no amount of spot-checking `$solve()` will reveal it. Reconstruction
# is the strongest single check because it exercises every factor at once.
#
# Sign and ordering are deliberately not asserted for eigenvectors: V and -V are
# both valid, and Eigen sorts ascending while base R sorts descending.

suppressPackageStartupMessages(library(ROpenEigen))

ok <- 0L; bad <- 0L
chk <- function(label, expr) {
  r <- tryCatch(isTRUE(expr), error = function(e) {
    cat(sprintf("  !! %-32s %s\n", label, conditionMessage(e))); NA
  })
  if (is.na(r)) { bad <<- bad + 1L; return(invisible(FALSE)) }
  if (r) ok <<- ok + 1L else { bad <<- bad + 1L; cat(sprintf("  ?? %s\n", label)) }
  invisible(r)
}

# 判据按 ||A|| 缩放：重建是后向稳定的，误差量级 ~ eps * ||A||，用绝对容差
# 在量级大的矩阵上会误报。
nrm <- function(x) {
  x <- as.matrix(x)
  if (length(x) == 0L) return(0)
  if (is.complex(x)) max(Mod(x)) else max(abs(x))
}
close <- function(a, b, tol = 1e-10) {
  a <- as.matrix(a); b <- as.matrix(b)
  if (!identical(dim(a), dim(b))) return(FALSE)
  scale <- max(1, nrm(a), nrm(b))
  max(Mod(a - b)) <= tol * scale
}
orthonormal <- function(q, tol = 1e-10) {
  q <- as.matrix(q)
  g <- Conj(t(q)) %*% q
  close(g, diag(ncol(q)), tol)
}
unw <- function(x) if (inherits(x, "eigenObject")) as.matrix(x) else x

cat("== 1. 实矩阵：每个分解的重建恒等式 ==\n")
for (n in c(1L, 2L, 3L, 5L, 8L)) {
  set.seed(1000 + n)
  A <- matrix(rnorm(n * n), n, n)
  S <- crossprod(A) + n * diag(n)          # 对称正定
  tag <- sprintf("n=%d", n)
  Es <- eigenMat(S); E <- eigenMat(A)

  chk(paste("llt 重建", tag),       close(eigen_llt(Es)$reconstructed(), S))
  chk(paste("llt L%*%t(L)", tag),   close(unw(eigen_llt(Es)$matrixL()) %*% t(unw(eigen_llt(Es)$matrixL())), S))
  chk(paste("ldlt 重建", tag),      close(eigen_ldlt(Es)$reconstructed(), S))
  chk(paste("lu 重建", tag),        close(eigen_lu(E)$reconstructed(), A))
  chk(paste("fulllu 重建", tag),    close(eigen_full_lu(E)$reconstructed(), A))
  chk(paste("qr 重建", tag),        close(eigen_qr(E)$reconstructed(), A))
  chk(paste("qr Q 正交", tag),      orthonormal(eigen_qr(E)$matrixQ()))
  chk(paste("colpiv 重建", tag),    close(eigen_colpiv_qr(E)$reconstructed(), A))
  chk(paste("colpiv Q 正交", tag),  orthonormal(eigen_colpiv_qr(E)$matrixQ()))
  chk(paste("svd 重建", tag),       close(eigen_svd(E)$reconstructed(), A))
  chk(paste("svd U 正交", tag),     orthonormal(eigen_svd(E)$u()))
  chk(paste("svd V 正交", tag),     orthonormal(eigen_svd(E)$v()))

  # LU：P %*% A == L %*% U
  lu <- eigen_lu(E)
  chk(paste("lu P*A=L*U", tag),
      close(unw(lu$permutationP()) %*% A, unw(lu$matrixL()) %*% unw(lu$matrixU())))
}

cat("== 2. 特征分解 ==\n")
for (n in c(1L, 2L, 3L, 6L)) {
  set.seed(2000 + n)
  A <- matrix(rnorm(n * n), n, n); S <- crossprod(A) + n * diag(n)
  Es <- eigenMat(S)
  d <- eigen_eigen(Es)
  V <- unw(d$vectors()); lam <- as.vector(unw(d$values()))
  chk(sprintf("eigen A=V L V' n=%d", n), close(V %*% diag(lam, n) %*% t(V), S))
  chk(sprintf("eigen V 正交 n=%d", n),  orthonormal(V))
  # Eigen 升序、base R 降序，比较排序后的值
  chk(sprintf("eigen 值对照 base n=%d", n),
      isTRUE(all.equal(sort(lam), sort(eigen(S, only.values = TRUE)$values), tolerance = 1e-8)))

  # 非对称：A V = V diag(lambda)
  E <- eigenMat(A); g <- eigen_geigen(E)
  Vc <- as.matrix(g$vectors()); lc <- as.vector(as.matrix(g$values()))
  chk(sprintf("geigen A V=V L n=%d", n),
      close(as.matrix(A) %*% Vc, Vc %*% diag(lc, n), 1e-8))
}

cat("== 3. 复数 ==\n")
for (n in c(2L, 4L)) {
  set.seed(3000 + n)
  R <- matrix(complex(real = rnorm(n * n), imaginary = rnorm(n * n)), n, n)
  H <- R %*% Conj(t(R)) + n * diag(n)        # Hermitian 正定
  Eh <- eigenMat(H); Ec <- eigenMat(R)

  chk(sprintf("复数 llt 重建 n=%d", n),   close(eigen_llt(Eh)$reconstructed(), H))
  chk(sprintf("复数 llt L L^H n=%d", n),
      close(unw(eigen_llt(Eh)$matrixL()) %*% Conj(t(unw(eigen_llt(Eh)$matrixL()))), H))
  chk(sprintf("复数 ldlt 重建 n=%d", n),  close(eigen_ldlt(Eh)$reconstructed(), H))
  chk(sprintf("复数 lu 重建 n=%d", n),    close(eigen_lu(Ec)$reconstructed(), R))
  chk(sprintf("复数 qr 重建 n=%d", n),    close(eigen_qr(Ec)$reconstructed(), R))
  chk(sprintf("复数 svd 重建 n=%d", n),   close(eigen_svd(Ec)$reconstructed(), R))
  chk(sprintf("复数 svd U 正交 n=%d", n), orthonormal(eigen_svd(Ec)$u()))
  chk(sprintf("复数 eigen 重建 n=%d", n),
      { d <- eigen_eigen(Eh); V <- as.matrix(d$vectors()); l <- as.vector(as.matrix(d$values()))
        close(V %*% diag(l, n) %*% Conj(t(V)), H, 1e-8) })
}

cat("== 4. 长方形 ==\n")
for (dims in list(c(5L, 3L), c(3L, 5L), c(6L, 6L))) {
  nr <- dims[1]; nc <- dims[2]
  set.seed(4000 + nr * 10 + nc)
  A <- matrix(rnorm(nr * nc), nr, nc); E <- eigenMat(A)
  tag <- sprintf("%dx%d", nr, nc)
  chk(paste("svd 重建", tag),     close(eigen_svd(E)$reconstructed(), A))
  chk(paste("qr 重建", tag),      close(eigen_qr(E)$reconstructed(), A))
  chk(paste("colpiv 重建", tag),  close(eigen_colpiv_qr(E)$reconstructed(), A))
  chk(paste("fulllu 重建", tag),  close(eigen_full_lu(E)$reconstructed(), A))
  # 最小二乘：拟合值应与 base qr.solve 一致（超定系统残差非零，比的是拟合值）
  b <- matrix(rnorm(nr), nr, 1)
  sol <- as.matrix(eigen_svd(E)$solve(b))
  ref <- try(qr.solve(A, b), silent = TRUE)
  if (!inherits(ref, "try-error"))
    chk(paste("svd 最小二乘拟合值", tag), close(A %*% sol, A %*% ref, 1e-8))
}

cat("== 5. 秩亏与奇异 ==\n")
set.seed(5000)
X <- matrix(rnorm(20), 5, 4); A <- X %*% matrix(rnorm(16), 4, 4)   # rank <= 4
E <- eigenMat(A)
chk("秩亏 colpiv rank", identical(as.integer(eigen_colpiv_qr(E)$rank()), 4L))
chk("秩亏 svd rank",    identical(as.integer(eigen_svd(E)$rank()), 4L))
chk("秩亏 fulllu rank", identical(as.integer(eigen_full_lu(E)$rank()), 4L))
chk("秩亏 svd 重建",    close(eigen_svd(E)$reconstructed(), A))
chk("秩亏 colpiv 重建", close(eigen_colpiv_qr(E)$reconstructed(), A))
Sing <- eigenMat(matrix(c(1, 2, 2, 4), 2))
chk("奇异 det 为 0",    abs(eigen_det(Sing)) < 1e-12)
chk("奇异 rcond 为 0",  eigen_rcond(Sing) < 1e-12)

cat("== 6. 分解对象与 base R 的数值一致性 ==\n")
set.seed(6000)
A <- matrix(rnorm(36), 6, 6); S <- crossprod(A) + 6 * diag(6)
Es <- eigenMat(S); E <- eigenMat(A)
chk("llt logdet",    close(eigen_llt(Es)$logdet(), as.vector(determinant(S)$modulus), 1e-10))
chk("ldlt logdet",   close(eigen_ldlt(Es)$logdet(), as.vector(determinant(S)$modulus), 1e-10))
chk("svd 奇异值",    close(as.vector(unw(eigen_svd(E)$singular_values())), svd(A)$d, 1e-10))
chk("qr log|det|",   close(eigen_qr(E)$log_abs_det(), sum(log(abs(svd(A)$d))), 1e-10))
chk("llt rcond",     close(eigen_llt(Es)$rcond(), eigen_rcond(Es), 1e-8))
chk("ldlt rcond",    close(eigen_ldlt(Es)$rcond(), eigen_rcond(Es), 1e-8))
chk("fulllu det",    close(as.vector(unw(eigen_full_lu(E)$det())), det(A), 1e-8))
chk("fulllu rcond 与 eigen_rcond 同口径",
    close(eigen_full_lu(E)$rcond(), eigen_rcond(E), 1e-8))

# LDLT 的存在理由就是不定矩阵：它不要求正定，重构出来必须还是原矩阵，
# rcond 也必须和 eigen_rcond() 同一个数（这条曾经不成立——Eigen 的
# LDLT::rcond() 是另一套估计）。
Ind <- matrix(c(1, 2, 2, 1), 2)
Ei <- eigenMat(Ind)
chk("ldlt 不定矩阵重建", close(eigen_ldlt(Ei)$reconstructed(), Ind, 1e-12))
chk("ldlt 不定矩阵 rcond", close(eigen_ldlt(Ei)$rcond(), eigen_rcond(Ei), 1e-8))

cat("== 7. SVD 双后端（Jacobi / BDC）==\n")
# 后端可以换，分解结果不能换。这条契约是整个 method 参数成立的前提：
# 如果三个后端给出的奇异值/秩/伪逆不一致，那 method 就不该暴露给用户。
set.seed(7000)
Big <- matrix(rnorm(40 * 40), 40, 40)
Eb <- eigenMat(Big)
svj <- eigen_svd(Eb, "jacobi")
svb <- eigen_svd(Eb, "bdc")
chk("auto 大矩阵选 bdc",      identical(eigen_svd(Eb)$backend(), "bdc"))
chk("auto 小矩阵选 jacobi",   identical(eigen_svd(eigenMat(matrix(1:4, 2)))$backend(), "jacobi"))
chk("显式 jacobi",            identical(svj$backend(), "jacobi"))
chk("显式 bdc",               identical(svb$backend(), "bdc"))
chk("两后端奇异值一致",
    close(as.vector(unw(svj$singular_values())), as.vector(unw(svb$singular_values())), 1e-12))
chk("两后端秩一致",           identical(as.integer(svj$rank()), as.integer(svb$rank())))
chk("两后端 threshold 一致",  close(svj$threshold(), svb$threshold(), 0))
chk("bdc 重建",               close(svb$reconstructed(), Big))
chk("bdc 伪逆",               close(svb$pseudo_inverse(), svj$pseudo_inverse(), 1e-10))
chk("bdc U 正交",             close(crossprod(as.matrix(svb$u())), diag(40), 1e-10))
chk("bdc V 正交",             close(crossprod(as.matrix(svb$v())), diag(40), 1e-10))

# thin U/V：bdc 也必须只给 min(n,p) 列，否则 n >> p 时会多出几 GB 零空间。
Tall <- matrix(rnorm(30 * 5), 30, 5)
Et <- eigenMat(Tall)
chk("bdc 高瘦 U 尺寸",  identical(dim(eigen_svd(Et, "bdc")$u()), c(30L, 5L)))
chk("bdc 高瘦 V 尺寸",  identical(dim(eigen_svd(Et, "bdc")$v()), c(5L, 5L)))
chk("bdc 高瘦重建",     close(eigen_svd(Et, "bdc")$reconstructed(), Tall))

# 宽矩阵（rows < cols）是另一条形状分支：BDCSVD 内部按 diagSize 方阵计算，
# 转置由它自己处理，这里必须实证而不是假设。
Wide <- matrix(rnorm(5 * 30), 5, 30)
Ew <- eigenMat(Wide)
chk("bdc 矮胖 U 尺寸",  identical(dim(eigen_svd(Ew, "bdc")$u()), c(5L, 5L)))
chk("bdc 矮胖 V 尺寸",  identical(dim(eigen_svd(Ew, "bdc")$v()), c(30L, 5L)))
chk("bdc 矮胖重建",     close(eigen_svd(Ew, "bdc")$reconstructed(), Wide))
chk("bdc 矮胖与 jacobi 奇异值一致",
    close(as.vector(unw(eigen_svd(Ew, "bdc")$singular_values())),
          as.vector(unw(eigen_svd(Ew, "jacobi")$singular_values())), 1e-12))

# 秩亏：BDCSVD 用 SVDBase 的同一个 rank()，截断口径必须和 Jacobi 一致。
Xr <- matrix(rnorm(24), 6, 4) %*% matrix(rnorm(16), 4, 4)
Er <- eigenMat(Xr)
br <- rnorm(6)   # 同一个右端项，否则比的不是同一个解
chk("bdc 秩亏 rank",     identical(as.integer(eigen_svd(Er, "bdc")$rank()), 4L))
chk("bdc 秩亏伪逆有限",  all(is.finite(as.matrix(eigen_svd(Er, "bdc")$pseudo_inverse()))))
chk("bdc 秩亏最小二乘",
    close(Xr %*% as.matrix(eigen_svd(Er, "bdc")$solve(br)),
          Xr %*% as.matrix(eigen_svd(Er, "jacobi")$solve(br)), 1e-6))

# 复数走的是另一份实例化，必须单独验。
Z <- matrix(rnorm(400) + 1i * rnorm(400), 20, 20)
Ez <- eigenMat(Z)
chk("bdc 复数重建",     close(eigen_svd(Ez, "bdc")$reconstructed(), Z))
chk("bdc 复数与 jacobi 一致",
    close(as.vector(unw(eigen_svd(Ez, "bdc")$singular_values())),
          as.vector(unw(eigen_svd(Ez, "jacobi")$singular_values())), 1e-12))

chk("非法 method 报错",
    inherits(try(eigen_svd(Eb, "nope"), silent = TRUE), "try-error"))

cat(sprintf("\n通过 %d, 失败 %d\n", ok, bad))
if (bad > 0L) stop(sprintf("%d assertion(s) failed", bad), call. = FALSE)
