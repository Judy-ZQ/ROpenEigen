# Base R compatibility surface.
#
# The det() bug was not "wrong answer when called" -- it was "silently wrong
# answer from a function nobody knew was dispatching". This file pins down the
# two halves of that contract:
#
#   same()  functions that must return exactly what base R returns
#   refuses()  functions we do not support, which must error rather than
#           quietly return something else
#
# The second half looks odd as a test, but it is the important one: if a future
# change makes one of these start returning a value, that is a behaviour change
# and this file should be updated deliberately, not by accident.
#
# Deliberately not covered: c(1, E), unlist(E). c() only dispatches on its first
# argument, and unlist() is not a generic, so neither can be hooked.

suppressPackageStartupMessages(library(ROpenEigen))

ok <- 0L; bad <- 0L
chk <- function(label, expr) {
  r <- tryCatch(isTRUE(expr), error = function(e) {
    cat(sprintf("  !! %-30s %s\n", label, conditionMessage(e))); NA
  })
  if (is.na(r)) { bad <<- bad + 1L; return(invisible(FALSE)) }
  if (r) ok <<- ok + 1L else { bad <<- bad + 1L; cat(sprintf("  ?? %s\n", label)) }
  invisible(r)
}

set.seed(20261103)
m  <- matrix(rnorm(24), 4, 6)
ms <- matrix(rnorm(16), 4, 4); ms <- ms + t(ms) + 4 * diag(4)
v  <- matrix(rnorm(4), 4, 1)
E  <- eigenMat(m); Es <- eigenMat(ms); Ev <- eigenMat(v)

# 结果若是 eigenMat，先落回 base 再比较；返回 base 对象的原样比较。
unwrap <- function(x) if (inherits(x, "eigenObject")) as.matrix(x) else x
same <- function(a, b) {
  a <- unwrap(a); b <- unwrap(b)
  if (is.list(a) || is.list(b)) return(identical(a, b))
  isTRUE(all.equal(a, b, tolerance = 1e-8, check.attributes = FALSE))
}

cat("== A. 与 base R 结果一致 ==\n")

# 不是泛型、但因为内部调用 as.matrix() 而"恰好能工作"的一批
chk("qr()$rank",        same(qr(E)$rank, qr(m)$rank))
chk("qr()$q",           same(qr(E)$q, qr(m)$q))
chk("svd()$d",          same(svd(E)$d, svd(m)$d))
chk("eigen()$values",   same(eigen(Es)$values, eigen(ms)$values))
chk("kappa()",          same(kappa(E), kappa(m)))
chk("rcond()",          same(rcond(E), rcond(m)))
chk("det()",            same(det(Es), det(ms)))
chk("max.col()",        same(max.col(E), max.col(m)))
chk("apply()",          same(apply(E, 1, sum), apply(m, 1, sum)))
chk("scale()",          same(scale(E), scale(m)))
chk("sweep()",          same(sweep(E, 2, 1), sweep(m, 2, 1)))
chk("prop.table()",     same(prop.table(E), prop.table(m)))
chk("head()",           same(head(E), head(m)))
chk("tail()",           same(tail(E), tail(m)))
chk("upper.tri()",      same(upper.tri(E), upper.tri(m)))
chk("lower.tri()",      same(lower.tri(E), lower.tri(m)))
chk("kronecker()",      same(kronecker(Es, Es), kronecker(ms, ms)))
chk("outer()",          same(outer(Ev, Ev), outer(v, v)))

# 本轮新挂上的 S3 方法
chk("is.array()",       is.array(E))
chk("is.matrix()",      is.matrix(E))
chk("diag()",           same(diag(eigenMat(matrix(c(1, 2, 3, 4), 2))), c(1, 4)))
chk("diag() 长方形",    same(diag(E), diag(m)))
chk("isSymmetric()",    same(isSymmetric(Es), isSymmetric(ms)))
chk("median()",         same(median(E), median(m)))
chk("quantile()",       same(quantile(E, c(.25, .75)), quantile(m, c(.25, .75))))
chk("sort()",           same(sort(E), sort(m)))
chk("rev()",            same(rev(E), rev(m)))
chk("unique() 按行",    same(unique(eigenMat(rbind(m[1, ], m[1, ], m[2, ]))),
                             unique(rbind(m[1, ], m[1, ], m[2, ]))))
chk("duplicated()",     same(duplicated(eigenMat(rbind(m[1, ], m[1, ]))),
                             duplicated(rbind(m[1, ], m[1, ]))))
chk("as.data.frame()",  same(as.data.frame(E), as.data.frame(m)))
chk("diff()",           same(diff(E), diff(m)))
chk("rbind()",          same(rbind(E, E), rbind(m, m)))
chk("cbind()",          same(cbind(E, E), cbind(m, m)))
chk("c() 展平",         same(c(E, 1), c(m, 1)))
chk("! 等于 == 0",      same(!eigenMat(matrix(c(0, 1, 2, 0), 2)), !matrix(c(0, 1, 2, 0), 2)))
chk("range()",          same(range(E), range(m)))

cat("== B. 不支持，必须报错（不能静默给别的结果）==\n")
refuses <- function(label, expr) chk(label, inherits(try(expr, silent = TRUE), "try-error"))

# 这几个**确实**挂不上方法：body 里没有 UseMethod，也不在 Ops/Summary/Math
# 这些组泛型里，dispatch 根本不存在。只能用 eigen_* 的等价物。
refuses("rowSums()",    rowSums(E))
refuses("colSums()",    colSums(E))
refuses("rowMeans()",   rowMeans(E))
refuses("colMeans()",   colMeans(E))
refuses("norm()",       norm(E, "O"))
refuses("sd()",         sd(E))
refuses("var()",        var(E))
refuses("which.min()",  which.min(E))
refuses("which.max()",  which.max(E))

# 这一段是修一个**判断错误**留下的。crossprod / tcrossprod / aperm /
# is.finite 原先就列在上面那份名单里，理由是"body 里没有 UseMethod"。
# 那个理由不成立：它们是**内部泛型**，dispatch 发生在 .Internal 内部，
# 读 body 看不出来。实测（scratch/probe_dispatch2.R）它们都会派发。
# 教训：判断"能不能挂方法"只能靠实测派发，不能靠读源码。
cat("== B2. 内部泛型：crossprod / tcrossprod / aperm / is.finite ==\n")
chk("crossprod(E)",         same(crossprod(E), crossprod(m)))
chk("crossprod(E, m)",      same(crossprod(E, m), crossprod(m)))
chk("crossprod(E, v)",      same(crossprod(E, v), crossprod(m, v)))
chk("crossprod(m, E) 也派发", same(crossprod(m, E), crossprod(m)))
chk("crossprod(v)",         same(crossprod(v), crossprod(v)))
chk("crossprod 保持 eigenMat", inherits(crossprod(E), "eigenMat"))
chk("tcrossprod(E)",        same(tcrossprod(E), tcrossprod(m)))
chk("tcrossprod(E, m)",     same(tcrossprod(E, m), tcrossprod(m)))
chk("tcrossprod(m, E) 也派发", same(tcrossprod(m, E), tcrossprod(m)))
chk("tcrossprod(v)",        same(tcrossprod(v), tcrossprod(v)))
# aperm：二维只有两种合法排列
chk("aperm 转置",           same(aperm(E, c(2, 1)), t(m)))
chk("aperm 恒等",           same(aperm(E, c(1, 2)), m))
chk("aperm 默认 perm 即转置", same(aperm(E), t(m)))
chk("aperm 非法 perm 报错",  inherits(try(aperm(E, c(1, 1)), silent = TRUE), "try-error"))
# is.finite 系列返回**逻辑矩阵**（eigenMat 只存 double，装不下逻辑值）
mfin <- m; mfin[1, 1] <- NA; mfin[2, 2] <- Inf; mfin[3, 3] <- NaN
Efin <- eigenMat(mfin)
chk("is.finite",   identical(is.finite(Efin), is.finite(mfin)))
chk("is.infinite", identical(is.infinite(Efin), is.infinite(mfin)))
chk("is.nan",      identical(is.nan(Efin), is.nan(mfin)))
chk("is.finite 返回逻辑矩阵", is.logical(is.finite(Efin)))

# 替代路径存在，确保用户不会走投无路
cat("== C. 替代路径可用 ==\n")
chk("crossprod 用 %*%",  same(t(as.matrix(E)) %*% as.matrix(E), crossprod(m)))
chk("rowSums 用 eigen_sum", same(eigen_sum(E), sum(m)))
# eigen_which_min() 返回 list(value, row, col)（1 基行列），which.min() 返回线性
# 下标；把后者换算成行列再比。
chk("which.min 替代", {
  w <- eigen_which_min(E); k <- which.min(m)
  w$value == min(m) && w$row == (k - 1L) %% nrow(m) + 1L && w$col == (k - 1L) %/% nrow(m) + 1L
})
chk("norm 替代",         same(eigen_norm(E, "O"), norm(m, "O")))
chk("is.finite 替代",    isTRUE(eigen_all_finite(E)))

cat("== D. which.min/which.max 与 base 逐位置一致 ==\n")
# 这里曾经有一个只在**列主序第一格**发作的 bug。Eigen 的 minCoeff(&r,&c) 拿
# 第一个元素当初始最优，而它的比较器 numext::mini(b,a) 等价于 (b<a)?b:a——
# 初始 a 是 NaN 时，后面每个 b<NaN 都是 false，NaN 就一直留着。结果：
# 3x3 全 1、把其中一格换成 NA，只有 NA 在 (1,1) 那一种情况返回 NA@(1,1)，
# 其余 8 种都正确跳过 NA。而 which.min(c(NA,1)) 是 2。
#
# 这种 bug 的特点是"只在特定位置上错"，所以测试必须遍历所有位置，不能只挑
# 一个。下面是 3x3 全 1、NA 依次放在每一格。
wlin <- function(w, nr)
  if (is.na(w$row)) NA_integer_ else as.integer((w$col - 1L) * nr + w$row)
for (k in 1:9) {
  A <- matrix(1, 3, 3)
  A[(k - 1) %% 3 + 1, (k - 1) %/% 3 + 1] <- NA
  E2 <- eigenMat(A)
  chk(sprintf("NA 在 (%d,%d) 时 which.min 与 base 一致", (k - 1) %% 3 + 1, (k - 1) %/% 3 + 1),
      identical(wlin(eigen_which_min(E2), 3), which.min(A)))
  chk(sprintf("NA 在 (%d,%d) 时 which.max 与 base 一致", (k - 1) %% 3 + 1, (k - 1) %/% 3 + 1),
      identical(wlin(eigen_which_max(E2), 3), which.max(A)))
}
# NaN 走同一条路径（ISNAN 同时认 NA 和 NaN）
chk("NaN 在首格", identical(wlin(eigen_which_min(eigenMat(matrix(c(NaN, 1, 2, 3), 2, 2))), 2),
                           which.min(matrix(c(NaN, 1, 2, 3), 2, 2))))
# Inf 不是缺失值：which.min(c(Inf,1)) 是 2，which.min(c(-Inf,1)) 是 1
Ainf <- matrix(c(Inf, 1, -Inf, 2), 2, 2)
chk("Inf 是值不是缺失", identical(wlin(eigen_which_min(eigenMat(Ainf)), 2), which.min(Ainf)))
chk("-Inf 是值不是缺失", identical(wlin(eigen_which_max(eigenMat(Ainf)), 2), which.max(Ainf)))
# 全部缺失：base 给 integer(0)，list 表达不了，给 NA 而不是编一个位置
Aall <- matrix(NA_real_, 2, 2)
wall <- eigen_which_min(eigenMat(Aall))
chk("全缺失给 NA", is.na(wall$value) && is.na(wall$row) && is.na(wall$col))
# 并列时取列主序第一个——和 which.min 的约定一致
chk("并列取第一个", {
  A <- matrix(c(0, 0, 0, 5, 5, 5), 3, 2)
  identical(wlin(eigen_which_min(eigenMat(A)), 3), which.min(A))
})

cat(sprintf("\n通过 %d, 失败 %d\n", ok, bad))
if (bad > 0L) stop(sprintf("%d assertion(s) failed", bad), call. = FALSE)
