suppressPackageStartupMessages(library(ROpenEigen))

ok <- 0L; bad <- 0L; fails <- character(0)
chk <- function(label, expr) {
  r <- tryCatch(isTRUE(expr), error = function(e) {
    bad <<- bad + 1L; fails <<- c(fails, paste0(label, " [error] ", conditionMessage(e)))
    return(NA)
  })
  if (is.na(r)) return(invisible(FALSE))
  if (r) ok <<- ok + 1L else { bad <<- bad + 1L; fails <<- c(fails, label) }
  invisible(r)
}
eq  <- function(a, b, tol = 1e-10) isTRUE(all.equal(a, b, tolerance = tol))
err <- function(e) inherits(try(e, silent = TRUE), "try-error")

cat("== A. 空与单元素 ==\n")
z0 <- eigenMat(matrix(numeric(0), 0, 0))
chk("0x0 dim",        eq(dim(z0), c(0L, 0L)))
chk("0x0 as.matrix",  eq(as.matrix(z0), matrix(numeric(0), 0, 0)))
chk("0x0 length",     eq(length(z0), 0L))
chk("0x0 sum",        eq(sum(z0), 0))
chk("0x0 print",      is.character(capture.output(print(z0))))
z03 <- eigenMat(matrix(numeric(0), 0, 3))
chk("0x3 dim",        eq(dim(z03), c(0L, 3L)))
chk("0x3 cbind 3x1",  eq(dim(eigen_cbind(z03, eigen_zeros(0, 1))), c(0L, 4L)))
one <- eigenMat(5)
chk("标量 dim",       eq(dim(one), c(1L, 1L)))
chk("标量 值",         eq(as.matrix(one)[1, 1], 5))
chk("标量 det",       eq(eigen_det(one), 5))
chk("标量 solve",     eq(as.vector(as.matrix(solve(one, 3))), 0.6))
chk("标量 t",         eq(dim(t(one)), c(1L, 1L)))

cat("== B. 转置 vs 共轭转置（复数）==\n")
C <- eigenMat(matrix(c(1+2i, 3-1i, 0+1i, 2+0i), 2, 2))
chk("t() 不共轭",      eq(as.matrix(t(C)), t(as.matrix(C))))
chk("t() 与 base 一致", eq(as.matrix(t(C)), t(matrix(c(1+2i,3-1i,0+1i,2+0i),2,2))))
chk("adjoint 共轭",    eq(as.matrix(eigen_adjoint(C)), Conj(t(as.matrix(C)))))
chk("adjoint != t",    !eq(as.matrix(eigen_adjoint(C)), as.matrix(t(C))))
chk("Conj",           eq(as.matrix(Conj(C)), Conj(as.matrix(C))))
chk("Mod",            eq(as.matrix(Mod(C)), Mod(as.matrix(C))))
chk("Arg",            eq(as.matrix(Arg(C)), Arg(as.matrix(C))))
chk("Re/Im",          eq(as.matrix(Re(C)), Re(as.matrix(C))) && eq(as.matrix(Im(C)), Im(as.matrix(C))))
chk("复数 det",       eq(eigen_det(C), { m <- as.matrix(C); m[1,1]*m[2,2] - m[1,2]*m[2,1] }))
# base R 在 modulus 上挂了 logarithm 属性，而且 exp() 会把它带下去，
# 所以比较时用 as.vector() 剥掉，否则 all.equal 报的是属性差异不是数值差异。
chk("复数 determinant", eq(as.vector(determinant(C)$modulus), log(Mod(eigen_det(C)))))
# base R 对复矩阵的 determinant() 直接报错，sign 的约定由本包定义：
# 单位复数 det/|det|，使 det = sign * exp(modulus) 仍然成立。
chk("复数 determinant 重构",
    eq(as.vector(determinant(C)$sign * exp(determinant(C)$modulus)), eigen_det(C)))
chk("复数 trace",     eq(eigen_trace(C), sum(diag(as.matrix(C)))))
chk("复数 t(t(C)) 往返", eq(as.matrix(t(t(C))), as.matrix(C)))

cat("== C. reshape / rep / block ==\n")
M <- matrix(1:12, 3, 4)
E <- eigenMat(M)
chk("reshape 3x4->4x3", eq(as.matrix(eigen_reshape(E, 4, 3)), matrix(1:12, 4, 3)))
chk("reshape 3x4->2x6", eq(as.matrix(eigen_reshape(E, 2, 6)), matrix(1:12, 2, 6)))
chk("reshape 元素数不符报错", err(eigen_reshape(E, 5, 5)))
chk("rep 2x2",        eq(as.matrix(eigen_rep(E, 2, 2)), M[rep(1:3, 2), rep(1:4, 2), drop = FALSE]))
chk("block 全取",      eq(as.matrix(eigen_block(E)), M))
chk("block 2:3, 2:4",  eq(as.matrix(eigen_block(E, 2:3, 2:4)), M[2:3, 2:4, drop = FALSE]))
chk("block 单行",      eq(as.matrix(eigen_block(E, 2, 1:4)), M[2, , drop = FALSE]))
chk("block 单列",      eq(dim(eigen_block(E, 1:3, 3)), c(3L, 1L)))

cat("== D. 下标边界 ==\n")
chk("[负下标]",   eq(as.matrix(E[-1, ]), M[-1, , drop = FALSE]))
chk("[逻辑下标]", eq(as.matrix(E[c(TRUE,FALSE,TRUE), ]), M[c(TRUE,FALSE,TRUE), , drop = FALSE]))
chk("[越界报错]", err(E[10, 1]))
chk("[0 下标返回空]", eq(length(E[0, 1]), length(M[0, 1])))
chk("[<- 标量]",  { X <- E; X[1, 1] <- 99; eq(as.matrix(X)[1,1], 99) })
chk("[<- 不改原件]", eq(as.matrix(E)[1,1], 1))
chk("[<- 整块]",  { X <- E; X[1:2, 1:2] <- matrix(c(0,0,0,0), 2); eq(sum(as.matrix(X)[1:2,1:2]), 0) })
# 赋值右侧是向量时要按目标块的形状铺开（曾报 "value is 3x1, block is 1x3"）
chk("[<- 行赋向量]", { X <- eigenMat(matrix(0, 3, 3)); X[1, ] <- c(1, 2, 3)
                       eq(as.matrix(X)[1, ], c(1, 2, 3)) })
chk("[<- 列赋向量]", { X <- eigenMat(matrix(0, 3, 3)); X[, 2] <- c(4, 5, 6)
                       eq(as.matrix(X)[, 2], c(4, 5, 6)) })
chk("[<- 标量铺满块]", { X <- eigenMat(matrix(0, 3, 3)); X[2:3, 2:3] <- 7
                         eq(sum(as.matrix(X)[2:3, 2:3]), 28) })
# 右侧是 eigenMat 时要先取内容再赋值（曾报 "number of items to replace ..."）
chk("[<- 右值是 eigenMat]", { X <- eigenMat(matrix(0, 4, 4)); X[1:2, 1:2] <- eigen_identity(2)
                               eq(as.matrix(X)[1:2, 1:2], diag(2)) })
chk("[<- 右值是 eigenMat 不改右值]", { Y <- eigen_identity(3)
                                       eq(as.matrix(Y), diag(3)) })

cat("== E. 形状不匹配必须报错 ==\n")
A2 <- eigenMat(matrix(1:4, 2, 2)); A3 <- eigenMat(matrix(1:6, 3, 2))
chk("2x2 + 3x2 报错", err(A2 + A3))
chk("2x2 %*% 2x2 ok", eq(as.matrix(A2 %*% A2), matrix(1:4,2,2) %*% matrix(1:4,2,2)))
chk("3x2 %*% 2x2 ok", eq(as.matrix(A3 %*% A2), matrix(1:6,3,2) %*% matrix(1:4,2,2)))
chk("2x2 %*% 3x2 报错", err(A2 %*% A3))
chk("不同长度 dot 报错", err(eigen_dot(eigenVec(1:3), eigenVec(1:2))))
chk("2 维 cross 返回标量 z", eq(as.vector(eigen_cross(eigenVec(c(1, 0)), eigenVec(c(0, 1)))), 1))
chk("2 维 cross 反对称",     eq(as.vector(eigen_cross(eigenVec(c(0, 1)), eigenVec(c(1, 0)))), -1))
chk("长度不等 cross 报错",   err(eigen_cross(eigenVec(1:3), eigenVec(1:2))))
chk("长度 4 cross 报错",     err(eigen_cross(eigenVec(1:4), eigenVec(1:4))))

cat("== F. 奇异 / 秩亏 ==\n")
Sg <- eigenMat(matrix(c(1, 2, 2, 4), 2))
chk("奇异 det ≈ 0",    abs(eigen_det(Sg)) < 1e-10)
chk("奇异 rank = 1",   eq(eigen_rank(Sg), 1))
chk("奇异 rcond ≈ 0",  eigen_rcond(Sg) < 1e-10)
chk("奇异 inverse 不崩", is.finite(eigen_rcond(Sg)))
chk("秩亏 pinv 与 ginv 一致", eq(as.matrix(eigen_pinv(Sg)), MASS::ginv(matrix(c(1,2,2,4),2)), tol = 1e-8))
chk("秩亏 svd$solve 有限", all(is.finite(as.matrix(eigen_svd(Sg)$solve(c(1, 1))))))

cat("== G. 非对称特征值 ==\n")
NS <- matrix(c(0, -1, 1, 0), 2)
ge <- eigen_geigen(eigenMat(NS))
chk("geigen 复数",     eigen_is_complex(ge$values()))
chk("geigen 值 ±i",    eq(sort(Im(as.vector(as.matrix(ge$values())))), c(-1, 1)))
chk("geigen 重建",     max(Mod(as.matrix(ge$vectors() %*% diag(as.vector(as.matrix(ge$values()))) %*% solve(as.matrix(ge$vectors()))) - NS)) < 1e-8)
chk("geigen 与 base 一致",
    eq(sort(round(Im(eigen(NS)$values), 8)), sort(round(Im(as.vector(as.matrix(ge$values()))), 8))))

cat("== H. 生成器边界 ==\n")
chk("linspace n=1",   eq(as.vector(as.matrix(eigen_linspace(0, 10, 1))), 0))
chk("linspace n=2",   eq(as.vector(as.matrix(eigen_linspace(0, 10, 2))), c(0, 10)))
chk("identity 1",     eq(as.matrix(eigen_identity(1)), matrix(1, 1, 1)))
chk("zeros 0x0",      eq(dim(eigen_zeros(0, 0)), c(0L, 0L)))
chk("constant 复数？", eq(as.matrix(eigen_constant(2, 2, 2)), matrix(2, 2, 2)))
set.seed(42); r1 <- as.matrix(eigen_random(3, 3))
set.seed(42); r2 <- as.matrix(eigen_random(3, 3))
chk("random 可复现",  eq(r1, r2))
chk("random 范围",    all(r1 >= 0 & r1 <= 1))

cat("== I. Summary / Math 的 na.rm ==\n")
N <- eigenMat(matrix(c(1, NA, 3, 4), 2, 2))
chk("sum 有 NA",      is.na(sum(N)))
chk("sum na.rm",      eq(sum(N, na.rm = TRUE), 8))
chk("mean na.rm",     eq(mean(N, na.rm = TRUE), 8/3))
chk("anyNA",          anyNA(N))
chk("is.na 矩阵",     eq(as.matrix(is.na(N)), is.na(matrix(c(1,NA,3,4),2,2))))
chk("Math 有 NA",     is.na(as.matrix(log(N))[2, 1]))
chk("log(x, base)",   eq(as.matrix(log(eigenMat(matrix(c(8, 16), 1, 2)), 2)), matrix(c(3, 4), 1, 2)))

cat("== J. 分解对象的基本契约 ==\n")
SPD <- eigenMat(crossprod(matrix(c(1,2,3,4,5,6,7,8,10), 3)))
for (nm in c("llt","ldlt","lu","full_lu","qr","colpiv_qr","svd","eigen")) {
  f <- get(paste0("eigen_", nm))
  obj <- f(SPD)
  chk(paste0(nm, ": names() 非空"), length(names(obj)) > 0)
  chk(paste0(nm, ": print 不崩"),   is.character(capture.output(print(obj))))
  chk(paste0(nm, ": $ 未知方法报错"), err(obj$no_such_method))
}
chk("llt 非法输入报错",   err(eigen_llt(eigenMat(matrix(1:6, 2, 3)))))
# 非正定：Eigen 只在 info() 里记一个 NumericalIssue，solve()/reconstructed()
# 照样返回值——返回的是另一个矩阵的分解（实测重构出 [[1,2],[2,5]]）。必须报错。
chk("llt 非正定报错",     err(eigen_llt(eigenMat(matrix(c(1, 2, 2, 1), 2)))))
chk("chol 非正定报错",    err(eigen_chol(eigenMat(matrix(c(1, 2, 2, 1), 2)))))
chk("ldlt 接受不定矩阵",  !err(eigen_ldlt(eigenMat(matrix(c(1, 2, 2, 1), 2)))))
# 正定的情形不受影响
chk("llt 正定可用",       !err(eigen_llt(eigenMat(matrix(c(4, 1, 1, 3), 2)))))
chk("eigen 非法输入报错", err(eigen_eigen(eigenMat(matrix(1:6, 2, 3)))))

cat("== K. 极端量级：Frobenius 范数不能溢出或下溢 ==\n")
# sqrt(sum(x^2)) 有两个真实的坏结局，而且都不报错：
#   元素 1e-200 -> 0     （平方下溢，真值 1.414e-200）
#   元素 1e300 -> Inf    （平方上溢，真值 1.414e+300）
# 返回 0 尤其糟：它看起来完全正常。base 的 norm(x, "F") 走 LAPACK 的 dlange，
# 本来就是缩放的，所以它在这两个量级上都对。eigen_norm 现在也缩放。
for (k in c(-300, -250, -200, -170, -160, 150, 155, 160, 200, 300)) {
  B <- matrix(c(10^k, 0, 0, 10^k), 2, 2)
  chk(sprintf("norm F 1e%d", k),
      isTRUE(all.equal(eigen_norm(eigenMat(B), "F"), norm(B, "F"), tolerance = 1e-12)))
}
# 单个元素的下溢：1e-200 的 F 范数就是 1e-200，不是 0
chk("norm F 单元素 1e-200 不是 0",
    isTRUE(all.equal(eigen_norm(eigenMat(matrix(1e-200)), "F"), 1e-200, tolerance = 1e-12)))
chk("norm F 单元素 1e-200 非零", eigen_norm(eigenMat(matrix(1e-200)), "F") > 0)
# 缩放不能改变正常量级的结果
chk("norm F 正常量级不变",
    isTRUE(all.equal(eigen_norm(eigenMat(matrix(c(3, 4, 0, 0), 2, 2)), "F"), 5,
                     tolerance = 1e-14)))
# Inf / NaN 的语义要和 base 一致：base 对 Inf 给 Inf，对 NaN 给 NaN
Ainf <- matrix(c(1, Inf, 2, 3), 2, 2)
chk("norm F 含 Inf", identical(eigen_norm(eigenMat(Ainf), "F"), norm(Ainf, "F")))
chk("norm F 含 Inf 是 Inf", is.infinite(eigen_norm(eigenMat(Ainf), "F")))
Anan <- matrix(c(1, NaN, 2, 3), 2, 2)
chk("norm F 含 NaN", identical(eigen_norm(eigenMat(Anan), "F"), norm(Anan, "F")))
# sumsq 保持朴素：base 的 sum(x^2) 在 1e300 上**就是** Inf，不该被"修好"
chk("sumsq 仍然溢出", is.infinite(eigen_sumsq(eigenMat(matrix(1e300)))))

cat("\n")
cat(sprintf("通过 %d, 失败 %d\n", ok, bad))
if (bad > 0L) { cat("失败项：\n"); cat(paste0("  - ", fails, collapse = "\n"), "\n") }
if (bad > 0L) stop(sprintf("%d assertion(s) failed", bad), call. = FALSE)
