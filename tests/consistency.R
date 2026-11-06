# Cross-entry-point consistency.
#
# Every other test file compares ROpenEigen against base R, or a decomposition
# against its own defining identity. Neither catches the failure mode that has
# actually bitten this package three times:
#
#   the same quantity, computed two ways, giving two different answers.
#
# Examples that shipped before being caught:
#   * `det(E)` returned exp(det(A)) while `eigen_det(E)` returned det(A)
#   * `llt$rcond()` and `eigen_rcond(E)` returned different numbers
#   * `eigen_llt()` accepted an indefinite matrix and factorised a *different*
#     matrix, so `$solve()` was confidently wrong
#
# A per-function test cannot see any of those, because each door is individually
# "working". This file opens all the doors and compares them to each other.
#
# Rule: if two functions claim to compute the same mathematical object, they must
# agree. Anything else is a bug in at least one of them.

suppressPackageStartupMessages(library(ROpenEigen))

ok <- 0L; bad <- 0L
chk <- function(label, expr) {
  r <- tryCatch(isTRUE(expr), error = function(e) {
    cat(sprintf("  !! %-34s %s\n", label, conditionMessage(e))); NA
  })
  if (is.na(r)) { bad <<- bad + 1L; return(invisible(FALSE)) }
  if (r) ok <<- ok + 1L else { bad <<- bad + 1L; cat(sprintf("  ?? %s\n", label)) }
  invisible(r)
}

nrm <- function(x) {
  x <- as.matrix(x)
  if (length(x) == 0L) return(0)
  if (is.complex(x)) max(Mod(x)) else max(abs(x))
}
close <- function(a, b, tol = 1e-9) {
  a <- as.matrix(a); b <- as.matrix(b)
  if (!identical(dim(a), dim(b))) return(FALSE)
  scale <- max(1, nrm(a), nrm(b))
  max(Mod(a - b)) <= tol * scale
}
unw <- function(x) if (inherits(x, "eigenObject")) as.matrix(x) else x
# 标量口径：det 这类量在良态矩阵上可能相差几个量级，用相对误差比绝对误差稳。
rel <- function(a, b, tol = 1e-8) {
  a <- as.numeric(a); b <- as.numeric(b)
  scale <- max(1, abs(a), abs(b))
  abs(a - b) <= tol * scale
}

# 四种性态：良态、病态、秩亏、奇异。跨入口一致性在这些矩阵上都该成立；
# 只在"恰好卡在数值秩边界"的矩阵上才可能分歧，所以那种不进来。
set.seed(2024)
n <- 6
G  <- matrix(rnorm(n * n), n, n)
SPD <- crossprod(G) + n * diag(n)              # 对称正定、良态
ILL <- {                                        # 病态但非奇异
  q <- qr.Q(qr(matrix(rnorm(n * n), n, n)))
  q %*% diag(10^seq(0, -12, length.out = n)) %*% t(q)
}
RD  <- matrix(rnorm(n * 3), n, 3) %*% matrix(rnorm(3 * n), 3, n)   # 秩 <= 3
SIN <- matrix(c(1, 2, 2, 4), 2)                # 精确奇异

mats <- list(
  list(name = "SPD",   A = SPD, sq = TRUE),
  list(name = "ill",   A = ILL, sq = TRUE),
  list(name = "rdef",  A = RD,  sq = TRUE),
  list(name = "sing",  A = SIN, sq = TRUE)
)

cat("== 1. det：每个入口必须给同一个数 ==\n")
for (m in mats) {
  A <- m$A; E <- eigenMat(A); tag <- m$name
  ref <- det(A)                                  # LAPACK
  chk(paste("det() == LAPACK", tag),      rel(det(E), ref))
  chk(paste("eigen_det == LAPACK", tag),  rel(eigen_det(E), ref))
  chk(paste("full_lu$det == LAPACK", tag), rel(unw(eigen_full_lu(E)$det()), ref))
  chk(paste("lu$det == LAPACK", tag),     rel(unw(eigen_lu(E)$det()), ref))
  # exp(logdet) 必须还原 |det|；奇异矩阵 logdet = -Inf，跳过这一条
  if (abs(ref) > 1e-8)
    chk(paste("exp(eigen_logdet) == |det|", tag),
        rel(exp(eigen_logdet(E)), abs(ref), 1e-7))
}

cat("== 2. log|det|：每个入口必须给同一个数 ==\n")
# log|det| 是病态量。cond = 1e12 时各算法的结果差 5e-8（绝对值），而 log 值本身
# 是 -83，用 1e-9 的**相对**判据会把这点噪声判成失败。所以这里用「绝对 + 相对」
# 混合容差：既要相对得上，也允许 1e-5 的绝对抖动。
ok_logdet <- function(a, b) abs(a - b) <= 1e-5 * max(1, abs(a), abs(b))
for (m in mats[1:2]) {
  A <- m$A; E <- eigenMat(A); tag <- m$name
  ref <- as.vector(determinant(A)$modulus)
  chk(paste("eigen_logdet", tag),        ok_logdet(eigen_logdet(E), ref))
  chk(paste("qr log_abs_det", tag),      ok_logdet(eigen_qr(E)$log_abs_det(), ref))
  chk(paste("colpiv log_abs_det", tag),  ok_logdet(eigen_colpiv_qr(E)$log_abs_det(), ref))
  chk(paste("determinant()$modulus", tag),
      ok_logdet(as.vector(determinant(E)$modulus), ref))
}
# 秩亏矩阵的 log|det| 真值是 -Inf，任何有限值都只是数值噪声（实测 LAPACK -108.4、
# 我们 -108.0、QR -103.1，三个数彼此差到 5，但都远小于 0）。这种矩阵上要求
# 「相等」没有数学意义，只能断言它确实塌到很负的位置。
Erd <- eigenMat(RD)
chk("秩亏 log|det| 塌到 -30 以下",     eigen_logdet(Erd) < -30)
chk("秩亏 qr log_abs_det 也塌了",     eigen_qr(Erd)$log_abs_det() < -30)
chk("秩亏 det 接近 0",                abs(eigen_det(Erd)) < 1e-6)
# 正定矩阵上 LLT / LDLT 的 logdet 也必须是同一个数
Es <- eigenMat(SPD)
refs <- as.vector(determinant(SPD)$modulus)
chk("llt$logdet == LAPACK",  rel(eigen_llt(Es)$logdet(), refs, 1e-9))
chk("ldlt$logdet == LAPACK", rel(eigen_ldlt(Es)$logdet(), refs, 1e-9))
chk("llt$logdet == eigen_logdet",  rel(eigen_llt(Es)$logdet(), eigen_logdet(Es), 1e-9))
chk("ldlt$logdet == eigen_logdet", rel(eigen_ldlt(Es)$logdet(), eigen_logdet(Es), 1e-9))
chk("llt$logdet == ldlt$logdet",
    rel(eigen_llt(Es)$logdet(), eigen_ldlt(Es)$logdet(), 1e-9))

cat("== 3. rank：每个入口必须给同一个数 ==\n")
# 注意这里**不含** ill。它的 sigma_min/sigma_max = 1e-12，正好落在两个默认阈值
# 中间：R 的 qr() 用 tol = 1e-7 判秩，Eigen 用 diagSize*eps ≈ 1.3e-15。于是
# base 说 rank 4，我们说 6——两边都对，只是容差不同。这种矩阵不能用来验一致性，
# 否则测的是「谁家的默认阈值」，不是「我们内部是否自洽」。
for (m in mats[c(1, 3, 4)]) {
  A <- m$A; E <- eigenMat(A); tag <- m$name
  ref <- as.integer(qr(A)$rank)                  # LAPACK
  chk(paste("eigen_rank", tag),     identical(as.integer(eigen_rank(E)), ref))
  chk(paste("full_lu$rank", tag),   identical(as.integer(eigen_full_lu(E)$rank()), ref))
  chk(paste("colpiv$rank", tag),    identical(as.integer(eigen_colpiv_qr(E)$rank()), ref))
  chk(paste("svd$rank", tag),       identical(as.integer(eigen_svd(E)$rank()), ref))
}
# 把上面这个容差差异钉成契约：eigen_rank() 用 Eigen 的紧阈值，比 base 的 qr()
# 宽松阈值更"敢于"认为满秩。这不该悄悄变。
chk("ill: eigen_rank 用紧阈值判满秩", identical(as.integer(eigen_rank(eigenMat(ILL))), 6L))
chk("ill: base qr() 用 1e-7 判亏秩",  as.integer(qr(ILL)$rank) < 6L)
# 但显式给 tol 时两边应该能对上
chk("显式 tol 后与 base 一致", identical(as.integer(eigen_rank(eigenMat(ILL), tol = 1e-7)), 4L))

cat("== 4. rcond：每个入口必须给同一个数 ==\n")
# 注意 base 的 rcond(A) **不是**同一个量：它是 1-norm 条件数的倒数
# （1/(||A||_1 * ||A^-1||_1)），我们的是 2-norm（sigma_min/sigma_max）。
# 实测 SPD 上 base 给 0.1693、我们给 0.2476，两个都对，只是口径不同。
# 真正该对上的是 1/kappa(A, exact = TRUE)，那也是 2-norm。
for (m in mats[1:2]) {
  A <- m$A; E <- eigenMat(A); tag <- m$name
  ref <- eigen_rcond(E)
  chk(paste("rcond == 1/kappa(exact)", tag), rel(ref, 1 / kappa(A, exact = TRUE), 1e-8))
  chk(paste("full_lu$rcond", tag), rel(eigen_full_lu(E)$rcond(), ref, 1e-8))
  svv <- as.vector(unw(eigen_svd(E)$singular_values()))
  chk(paste("rcond == sigma ratio", tag), rel(ref, svv[length(svv)] / svv[1], 1e-8))
}
# 秩亏矩阵的 rcond 是 1/eps 量级，相对比较没有意义（两个数都是噪声），
# 只断言它确实被判成数值奇异。
chk("秩亏 rcond < 1e-10", eigen_rcond(eigenMat(RD)) < 1e-10)
chk("奇异 rcond < 1e-14", eigen_rcond(eigenMat(SIN)) < 1e-14)
Es <- eigenMat(SPD)
chk("llt$rcond == eigen_rcond",  rel(eigen_llt(Es)$rcond(), eigen_rcond(Es), 1e-8))
chk("ldlt$rcond == eigen_rcond", rel(eigen_ldlt(Es)$rcond(), eigen_rcond(Es), 1e-8))

cat("== 5. solve：每个入口必须给同一个解 ==\n")
b <- rnorm(n)
for (nm in c("lu", "full_lu", "qr", "svd")) {
  obj <- switch(nm, lu = eigen_lu, full_lu = eigen_full_lu, qr = eigen_qr, svd = eigen_svd)
  chk(paste(nm, "$solve == base solve"), close(unw(obj(eigenMat(SPD))$solve(b)), solve(SPD, b)))
}
chk("llt$solve == base",  close(unw(eigen_llt(Es)$solve(b)), solve(SPD, b)))
chk("ldlt$solve == base", close(unw(eigen_ldlt(Es)$solve(b)), solve(SPD, b)))
chk("eigen_solve == base", close(unw(eigen_solve(eigenMat(SPD), b)), solve(SPD, b)))
chk("solve(eigenMat) == base", close(unw(solve(eigenMat(SPD), b)), solve(SPD, b)))
# 全部与 llt 的解对齐：close() 是传递的，都贴着同一个参照就等于两两一致。
sol <- lapply(list(llt = eigen_llt(Es), ldlt = eigen_ldlt(Es),
                   lu = eigen_lu(eigenMat(SPD)), qr = eigen_qr(eigenMat(SPD)),
                   svd = eigen_svd(eigenMat(SPD))), function(o) unw(o$solve(b)))
for (k in 2:length(sol)) chk(paste("solve 与 llt 一致:", names(sol)[k]), close(sol[[k]], sol[[1]]))

cat("== 6. inverse：每个入口必须给同一个逆 ==\n")
refinv <- solve(SPD)
chk("eigen_inverse == base",       close(unw(eigen_inverse(eigenMat(SPD))), refinv))
chk("lu$inverse == base",          close(unw(eigen_lu(eigenMat(SPD))$inverse()), refinv))
chk("full_lu$inverse == base",     close(unw(eigen_full_lu(eigenMat(SPD))$inverse()), refinv))
chk("svd$pseudo_inverse == base",  close(unw(eigen_svd(eigenMat(SPD))$pseudo_inverse()), refinv))
chk("eigen_pinv == eigen_inverse", close(unw(eigen_pinv(eigenMat(SPD))), unw(eigen_inverse(eigenMat(SPD)))))
chk("pinv == svd$pseudo_inverse",
    close(unw(eigen_pinv(eigenMat(SPD))), unw(eigen_svd(eigenMat(SPD))$pseudo_inverse())))

cat("== 7. 一次性与对象式入口一致 ==\n")
E <- eigenMat(G); Eg <- eigenMat(SPD)
chk("eigen_singular_values == svd$singular_values",
    close(unw(eigen_singular_values(E)), unw(eigen_svd(E)$singular_values())))
chk("eigen_values == eigen$values",
    close(unw(eigen_values(Eg)), unw(eigen_eigen(Eg)$values())))
chk("eigen_vectors == eigen$vectors",
    close(unw(eigen_vectors(Eg)), unw(eigen_eigen(Eg)$vectors())))
chk("eigen_chol == llt$matrixL",
    close(unw(eigen_chol(Eg)), unw(eigen_llt(Eg)$matrixL())))
chk("eigen_det == determinant()",
    rel(eigen_det(E), as.vector(det(E))))
chk("eigen_trace == sum(diag)", rel(eigen_trace(E), sum(diag(G))))
chk("eigen_trace == sum(eigen_values)", rel(eigen_trace(Eg), sum(unw(eigen_values(Eg)))))

cat("== 8. 范数：每个口径都必须对得上 ==\n")
for (ty in c("1", "2", "F", "I")) {
  chk(paste("eigen_norm", ty, "== base"), rel(eigen_norm(E, ty), norm(G, ty), 1e-9))
}
sv <- as.vector(unw(eigen_svd(E)$singular_values()))
chk("norm-2 == max sigma", rel(eigen_norm(E, "2"), sv[1], 1e-9))
chk("norm-F == sqrt(sum sigma^2)", rel(eigen_norm(E, "F"), sqrt(sum(sv^2)), 1e-9))
chk("max sigma == sqrt(max eigenvalue of A^T A)",
    rel(sv[1], sqrt(max(as.vector(unw(eigen_values(eigenMat(crossprod(G))))))), 1e-7))

cat("== 9. 伪逆满足 Moore-Penrose 四条件 ==\n")
# 秩亏矩阵最能考验截断口径：四条恒等式同时成立，说明零空间处理正确。
# 长方形必须单独测：只有那里 A^+ 真的不等于 A^-1，四条恒等式才是有效判据
# （方阵上它们退化成 A P = I，测不出截断是否合理）。
set.seed(9001)
RC <- list(
  list(nm = "SPD",   A = SPD),
  list(nm = "rdef",  A = RD),
  list(nm = "sing",  A = SIN),
  list(nm = "tall",  A = matrix(rnorm(15), 5, 3)),   # 5x3，超定
  list(nm = "wide",  A = matrix(rnorm(15), 3, 5)),   # 3x5，欠定
  list(nm = "rtall", A = matrix(rnorm(24), 6, 4) %*% matrix(rnorm(16), 4, 4))  # 秩亏长方
)
for (m in RC) {
  A <- m$A; P <- unw(eigen_pinv(eigenMat(A))); tag <- m$nm
  chk(paste("A P A == A", tag),    close(A %*% P %*% A, A, 1e-8))
  chk(paste("P A P == P", tag),    close(P %*% A %*% P, P, 1e-8))
  chk(paste("(A P)' == A P", tag), close(t(A %*% P), A %*% P, 1e-8))
  chk(paste("(P A)' == P A", tag), close(t(P %*% A), P %*% A, 1e-8))
  # 形状契约：P 必须是 ncol x nrow
  chk(paste("pinv 形状", tag), identical(dim(P), rev(dim(A))))
}

cat("== 10. 转置不改变任何不变量 ==\n")
for (m in mats[1:3]) {
  A <- m$A; E <- eigenMat(A); Et <- t(E); tag <- m$name
  chk(paste("det(t(E)) == det(E)", tag), rel(eigen_det(Et), eigen_det(E), 1e-8))
  chk(paste("trace(t(E)) == trace(E)", tag), rel(eigen_trace(Et), eigen_trace(E)))
  chk(paste("sigma(t(E)) == sigma(E)", tag),
      close(unw(eigen_svd(Et)$singular_values()), unw(eigen_svd(E)$singular_values()), 1e-9))
  chk(paste("norm-2(t(E)) == norm-2(E)", tag),
      rel(eigen_norm(Et, "2"), eigen_norm(E, "2"), 1e-9))
  chk(paste("rank(t(E)) == rank(E)", tag), identical(as.integer(eigen_rank(Et)), as.integer(eigen_rank(E))))
}

cat("== 11. 正定矩阵：特征值与奇异值必须对得上 ==\n")
# 对正定矩阵 sigma_i = lambda_i（不是 lambda_i^2！平方关系是 A'A 的谱，
# 别把 lambda 和 sigma^2 混了——这正是我第一版测试写错的地方）。
Es <- eigenMat(SPD)
ev <- sort(as.vector(unw(eigen_values(Es))))
sg <- sort(as.vector(unw(eigen_svd(Es)$singular_values())))
chk("eigen_values == sort(sigma)",  close(ev, sg, 1e-8))
chk("eigen_values 全正",            all(ev > 0))
chk("eigen_geigen 给同一组特征值",
    close(ev, sort(Re(as.vector(unw(eigen_geigen(Es)$values())))), 1e-8))
chk("eigen_gvalues == eigen_values",
    close(unw(eigen_gvalues(Es)), unw(eigen_values(Es)), 1e-8))
# LDLT 的 D 是对角主元，不是特征值（实测 9.08 vs 6.96）；sum(D) 也不是 trace
# （A = P'LDL'P，trace(A) = trace(D L'L)）。真正成立的是 prod(D) = det(A)，
# 因为 det(L) = 1。所以只能拿行列式比，不能拿迹比。
ld <- as.vector(unw(eigen_ldlt(Es)$vectorD()))
chk("prod(ldlt$vectorD) == det",
    rel(prod(ld), eigen_det(Es), 1e-8))
chk("sum(log|ldlt$vectorD|) == log|det|",
    rel(sum(log(abs(ld))), eigen_logdet(Es), 1e-8))
chk("sum(eigen_values) == trace",   rel(sum(ev), eigen_trace(Es), 1e-9))
chk("trace == sum(diag(A))",        rel(eigen_trace(Es), sum(diag(SPD)), 1e-12))

cat("== 12. 一次性函数 vs base R 的同义函数 ==\n")
# 每个 eigen_* 一次性函数都有一个 base R 对应物。它们是同一把门的两种开法，
# 结果必须一致；不一致就说明至少有一边的语义写歪了。
A <- matrix(c(1, 3, 5, 2, 4, 6), 2, 3)     # 与 probe 用同一布局，便于核对
Ea <- eigenMat(A)
Z <- matrix(c(1 + 2i, 2, 3, 4 - 2i), 2)
Ez <- eigenMat(Z)

chk("eigen_sum == sum",        rel(eigen_sum(Ea), sum(A), 1e-12))
chk("eigen_sum 复数",          close(eigen_sum(Ez), sum(Z), 1e-12))
chk("eigen_prod == prod",      rel(eigen_prod(Ea), prod(A), 1e-12))
chk("eigen_mean == mean",      rel(eigen_mean(Ea), mean(A), 1e-12))
chk("eigen_sumsq == sum(abs^2)", rel(eigen_sumsq(Ez), sum(abs(Z)^2), 1e-12))
chk("eigen_sumsq == sum(A^2) 实数", rel(eigen_sumsq(Ea), sum(A^2), 1e-12))

# dot 的约定：sum(a*b)，**不取共轭**。实数下与 base 完全一致；
# 复数下刻意不跟 Eigen 的 MatrixBase::dot()（那个会共轭第一个参数）。
u <- c(1, 2, 3); v <- c(4, 5, 6)
chk("eigen_dot == sum(a*b) 实数", rel(eigen_dot(u, v), sum(u * v), 1e-12))
a <- c(1 + 2i, 3 - 1i); b <- c(2 - 1i, 1 + 1i)
chk("eigen_dot == sum(a*b) 复数", close(eigen_dot(a, b), sum(a * b), 1e-12))
# 把这个"刻意的不一致"钉住：复数下 dot(x,x) 不是平方范数，而是 sum(x*x)。
chk("dot(z,z) == sum(z*z)",       close(eigen_dot(Z, Z), sum(Z * Z), 1e-12))
chk("dot(z,z) != sumsq(z) 复数",  !isTRUE(all.equal(eigen_dot(Z, Z), eigen_sumsq(Ez))))
# 文档推荐用 eigen_adjoint(a) %*% b 取 Hermitian 内积，那就得真的等于 sum(Conj(a)*b)。
chk("adjoint(a) %*% b == Hermitian 内积",
    close(unw(eigen_adjoint(a) %*% b), sum(Conj(a) * b), 1e-12))

e1 <- c(1, 0, 0); e2 <- c(0, 1, 0)
chk("eigen_cross(e1,e2) == e3", close(unw(eigen_cross(e1, e2)), c(0, 0, 1), 1e-12))
chk("cross 与两输入正交", {
  w <- as.vector(unw(eigen_cross(u, v)))
  abs(sum(w * u)) < 1e-10 && abs(sum(w * v)) < 1e-10
})
chk("eigen_outer == %o%", close(unw(eigen_outer(u, v)), u %o% v, 1e-12))

# which_min 返回 list(value, row, col)，不是线性下标——转换后才能比。
wm <- unw(eigen_which_min(Ea))
chk("which_min value", rel(wm$value, min(A), 1e-12))
chk("which_min 位置",  identical(c(wm$row, wm$col), as.integer(which(A == min(A), arr.ind = TRUE)[1, ])))
wM <- unw(eigen_which_max(Ea))
chk("which_max value", rel(wM$value, max(A), 1e-12))
chk("which_max 位置",  identical(c(wM$row, wM$col), as.integer(which(A == max(A), arr.ind = TRUE)[1, ])))

chk("eigen_approx 相等为真", isTRUE(eigen_approx(Ea, eigenMat(A))))
chk("eigen_approx 不同为假", isFALSE(eigen_approx(Ea, eigenMat(A + 1))))

chk("eigen_reshape == matrix()", close(unw(eigen_reshape(Ea, 3, 2)), matrix(A, 3, 2), 1e-12))
chk("eigen_block 行列块",        close(unw(eigen_block(Ea, 1:2, 2:3)), A[1:2, 2:3], 1e-12))
chk("eigen_block 只有行",        close(unw(eigen_block(Ea, 2)), A[2, , drop = FALSE], 1e-12))
chk("eigen_diag == diag",        close(unw(eigen_diag(Ea)), as.matrix(diag(A)), 1e-12))

chk("eigen_identity == diag",    close(unw(eigen_identity(4)), diag(4), 1e-12))
chk("eigen_zeros",               close(unw(eigen_zeros(3, 2)), matrix(0, 3, 2), 1e-12))
chk("eigen_ones",                close(unw(eigen_ones(3, 2)), matrix(1, 3, 2), 1e-12))
chk("eigen_constant",            close(unw(eigen_constant(7, 2, 3)), matrix(7, 2, 3), 1e-12))

chk("eigen_rbind == rbind",      close(unw(eigen_rbind(Ea, Ea)), rbind(A, A), 1e-12))
chk("eigen_cbind == cbind",      close(unw(eigen_cbind(Ea, Ea)), cbind(A, A), 1e-12))
chk("eigen_bind along=row",      close(unw(eigen_bind(Ea, Ea, along = "row")), rbind(A, A), 1e-12))
chk("eigen_bind along=col",      close(unw(eigen_bind(Ea, Ea, along = "col")), cbind(A, A), 1e-12))
chk("eigen_rep 横向平铺",        close(unw(eigen_rep(Ea, 1, 2)), cbind(A, A), 1e-12))
chk("eigen_rep 纵向平铺",        close(unw(eigen_rep(Ea, 2, 1)), rbind(A, A), 1e-12))

chk("eigen_adjoint == Conj(t())", close(unw(eigen_adjoint(Ez)), Conj(t(Z)), 1e-12))
chk("eigen_linspace == seq",     close(unw(eigen_linspace(0, 1, 5)),
                                       as.matrix(seq(0, 1, length.out = 5)), 1e-12))

N <- eigenMat(matrix(c(1, NA, Inf, 4), 2))
chk("all_finite 有 NA/Inf 为假", isFALSE(eigen_all_finite(N)))
chk("has_nan 为真",              isTRUE(eigen_has_nan(N)))
chk("has_inf 为真",              isTRUE(eigen_has_inf(N)))
chk("all_finite 干净矩阵为真",   isTRUE(eigen_all_finite(Ea)))
chk("is_complex 实数",           isFALSE(eigen_is_complex(Ea)))
chk("is_complex 复数",           isTRUE(eigen_is_complex(Ez)))

cat("== 13. 泛型入口与一次性函数必须一致 ==\n")
# sum(E) 走 Summary 组，eigen_sum(E) 直接进 C++。两条路算的是同一个和。
chk("sum(E) == eigen_sum(E)",   rel(sum(Ea), eigen_sum(Ea), 1e-12))
chk("prod(E) == eigen_prod(E)", rel(prod(Ea), eigen_prod(Ea), 1e-12))
chk("mean(E) == eigen_mean(E)", rel(mean(Ea), eigen_mean(Ea), 1e-12))
chk("t(E) == eigen_adjoint(E) 实数", close(unw(t(Ea)), unw(eigen_adjoint(Ea)), 1e-12))
chk("diag(E) == eigen_diag(E)", close(unw(diag(Ea)), unw(eigen_diag(Ea)), 1e-12))
chk("rbind(E,E) == eigen_rbind", close(unw(rbind(Ea, Ea)), unw(eigen_rbind(Ea, Ea)), 1e-12))
chk("cbind(E,E) == eigen_cbind", close(unw(cbind(Ea, Ea)), unw(eigen_cbind(Ea, Ea)), 1e-12))
chk("dim(E) == dim(A)",         identical(as.integer(dim(Ea)), as.integer(dim(A))))
chk("nrow/ncol",                identical(c(nrow(Ea), ncol(Ea)), c(2L, 3L)))
chk("length(E) == length(A)",   identical(as.integer(length(Ea)), as.integer(length(A))))
chk("as.vector(E) == as.vector(A)", close(as.vector(Ea), as.vector(A), 1e-12))
chk("as.matrix(E) 往返无损",     identical(as.matrix(eigenMat(A)), A))
chk("sum(E) 带 na.rm 走 base",   rel(sum(N, na.rm = TRUE), 5, 1e-12))

cat(sprintf("\n通过 %d, 失败 %d\n", ok, bad))
if (bad > 0L) stop(sprintf("%d assertion(s) failed", bad), call. = FALSE)
