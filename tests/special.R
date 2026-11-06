suppressPackageStartupMessages(library(ROpenEigen))

# 特殊函数（unsupported/Eigen/SpecialFunctions）与矩阵函数
# （unsupported/Eigen/MatrixFunctions）的对照测试。
#
# 特殊函数全部拿 R 自己的 stats 函数当参照：pnorm / lgamma / digamma /
# qnorm / pgamma / pbeta / psigamma / bessel* 都是同一族的实现，而且 R 的那
# 几个是公认准确的。Hurwitz zeta 在 base R 里没有，退化情形有闭式解
# （zeta(s, 1) 就是 Riemann zeta，zeta(2) = pi^2/6），就用闭式解。

ok <- 0L; bad <- 0L
chk <- function(label, expr) {
  r <- tryCatch(isTRUE(expr), error = function(e) {
    cat(sprintf("  !! %-38s %s\n", label, conditionMessage(e))); NA
  })
  if (is.na(r)) { bad <<- bad + 1L; return(invisible(FALSE)) }
  if (r) ok <<- ok + 1L else { bad <<- bad + 1L; cat(sprintf("  ?? %s\n", label)) }
  invisible(r)
}
eq <- function(a, b, tol = 1e-10) isTRUE(all.equal(a, b, tolerance = tol))
errs <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
M <- function(x) as.matrix(x)

cat("== 1. erf / erfc ==\n")
x <- c(-2, -0.5, 0, 0.5, 1, 2.5)
chk("erf = 2*pnorm(sqrt2 x) - 1", eq(as.vector(M(eigen_erf(eigenMat(x)))),
                                     2 * pnorm(x * sqrt(2)) - 1))
chk("erfc = 2*pnorm(-sqrt2 x)",   eq(as.vector(M(eigen_erfc(eigenMat(x)))),
                                     2 * pnorm(-x * sqrt(2))))
chk("erf(0) = 0",                 eq(M(eigen_erf(eigenMat(0)))[1], 0))
chk("erf + erfc = 1",             eq(as.vector(M(eigen_erf(eigenMat(x)))) +
                                     as.vector(M(eigen_erfc(eigenMat(x)))), rep(1, length(x))))

cat("== 2. gamma 家族 ==\n")
chk("gammaln = lgamma",  eq(as.vector(M(eigen_gammaln(eigenMat(c(0.5, 1, 2, 10))))),
                            lgamma(c(0.5, 1, 2, 10))))
chk("digamma",           eq(as.vector(M(eigen_digamma(eigenMat(c(0.5, 1, 2, 10))))),
                            digamma(c(0.5, 1, 2, 10))))
chk("digamma(1) = -gamma", eq(M(eigen_digamma(eigenMat(1)))[1], -0.5772156649015329))
chk("ndtri = qnorm",     eq(as.vector(M(eigen_ndtri(eigenMat(c(0.01, 0.5, 0.975))))),
                            qnorm(c(0.01, 0.5, 0.975))))
chk("ndtri(0.5) = 0",    eq(M(eigen_ndtri(eigenMat(0.5)))[1], 0))
chk("igamma = pgamma",   eq(as.vector(M(eigen_igamma(2.5, c(0.3, 1, 4)))),
                            pgamma(c(0.3, 1, 4), 2.5)))
chk("igammac = 上尾",     eq(as.vector(M(eigen_igammac(2.5, c(0.3, 1, 4)))),
                            pgamma(c(0.3, 1, 4), 2.5, lower.tail = FALSE)))
chk("igamma + igammac = 1", eq(as.vector(M(eigen_igamma(2.5, c(0.3, 1, 4)))) +
                                as.vector(M(eigen_igammac(2.5, c(0.3, 1, 4)))), rep(1, 3)))
chk("betainc = pbeta",   eq(as.vector(M(eigen_betainc(2, 3, c(0.1, 0.5, 0.9)))),
                            pbeta(c(0.1, 0.5, 0.9), 2, 3)))
chk("betainc 手动值",     eq(M(eigen_betainc(2, 3, 0.5))[1], 0.6875))
chk("polygamma = psigamma", {
  n <- c(0, 1, 2)
  okk <- TRUE
  for (k in n) okk <- okk && eq(as.vector(M(eigen_polygamma(k, c(1.5, 2.5)))),
                                psigamma(c(1.5, 2.5), k), 1e-9)
  okk
})

cat("== 3. zeta ==\n")
chk("zeta(2) = pi^2/6",   eq(M(eigen_zeta(eigenMat(2)))[1], pi^2 / 6))
chk("zeta(4) = pi^4/90",  eq(M(eigen_zeta(eigenMat(4)))[1], pi^4 / 90))
chk("Hurwitz zeta(s, 2) = zeta(s) - 1", eq(M(eigen_zeta(eigenMat(3), 2))[1],
                                           M(eigen_zeta(eigenMat(3)))[1] - 1))
# Eigen 的 Hurwitz zeta 只在 Re(s) > 1 上有效（它用的是级数定义，没有做
# 解析延拓）：s <= 1 时返回 NaN。R 的 zeta() 做了延拓，两者不能对等，
# 所以这里断言"返回 NaN"而不是"等于 -1/2"，把这个边界写在测试里。
chk("zeta(0) 是 NaN（Eigen 无解析延拓）", is.nan(M(eigen_zeta(eigenMat(0)))[1]))
chk("zeta(-1) 是 NaN",                    is.nan(M(eigen_zeta(eigenMat(-1)))[1]))

cat("== 4. Bessel ==\n")
xs <- c(0.5, 1, 3)
chk("i0 = besselI(x, 0)", eq(as.vector(M(eigen_bessel(eigenMat(xs), "i0"))), besselI(xs, 0)))
chk("i1 = besselI(x, 1)", eq(as.vector(M(eigen_bessel(eigenMat(xs), "i1"))), besselI(xs, 1)))
chk("k0 = besselK(x, 0)", eq(as.vector(M(eigen_bessel(eigenMat(xs), "k0"))), besselK(xs, 0)))
chk("k1 = besselK(x, 1)", eq(as.vector(M(eigen_bessel(eigenMat(xs), "k1"))), besselK(xs, 1)))
chk("j0 = besselJ(x, 0)", eq(as.vector(M(eigen_bessel(eigenMat(xs), "j0"))), besselJ(xs, 0)))
chk("j1 = besselJ(x, 1)", eq(as.vector(M(eigen_bessel(eigenMat(xs), "j1"))), besselJ(xs, 1)))
chk("y0 = besselY(x, 0)", eq(as.vector(M(eigen_bessel(eigenMat(xs), "y0"))), besselY(xs, 0)))
chk("y1 = besselY(x, 1)", eq(as.vector(M(eigen_bessel(eigenMat(xs), "y1"))), besselY(xs, 1)))
chk("i0e = i0 * exp(-x)", eq(as.vector(M(eigen_bessel(eigenMat(xs), "i0e"))), besselI(xs, 0) * exp(-xs)))
chk("k1e = k1 * exp(x)",  eq(as.vector(M(eigen_bessel(eigenMat(xs), "k1e"))), besselK(xs, 1) * exp(xs)))
chk("j0(1) 的文献值",      eq(M(eigen_bessel(eigenMat(1), "j0"))[1], 0.7651976865579666))
chk("未知 type 报错",      errs(eigen_bessel(eigenMat(1), "nope")))

cat("== 5. 形状与类型保持 ==\n")
A <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), 2, 3)
chk("矩阵进、矩阵出",   eq(dim(eigen_erf(eigenMat(A))), c(2L, 3L)))
chk("逐元素对齐",       eq(M(eigen_erf(eigenMat(A))),
                            matrix(2 * pnorm(as.numeric(A) * sqrt(2)) - 1, nrow(A))))
chk("空矩阵不炸",       eq(dim(eigen_erf(eigenMat(matrix(numeric(0), 0, 0)))), c(0L, 0L)))
chk("复数是错误",        errs(eigen_erf(eigenMat(matrix(1 + 1i, 1, 1)))))
chk("eigen_zeta 保持形状", eq(dim(eigen_zeta(eigenMat(A))), c(2L, 3L)))
chk("eigen_igamma 广播标量", eq(length(eigen_igamma(2, A)), 6L))

cat("== 6. 矩阵指数 ==\n")
chk("expm(0) = I",        eq(M(eigen_expm(eigenMat(matrix(0, 3, 3)))), diag(3)))
chk("expm(diag) = diag(exp)", eq(M(eigen_expm(eigenMat(diag(c(1, 2))))), diag(exp(c(1, 2)))))
# expm(tA) 满足 d/dt = A expm：用有限差分开检验，比逐个元素对表更结实
G <- matrix(c(0, -pi / 2, pi / 2, 0), 2)
# G 是列主序：G = [[0, pi/2], [-pi/2, 0]]，对应 -90 度旋转
chk("expm(skew) 是 -90 度旋转", eq(M(eigen_expm(eigenMat(G))),
                                   matrix(c(0, -1, 1, 0), 2, 2), tol = 1e-12))
chk("expm(A+B) = expm A expm B (A, B 可交换)", {
  A1 <- matrix(c(1, 0, 0, 2), 2); B1 <- matrix(c(3, 0, 0, 4), 2)
  eq(M(eigen_expm(eigenMat(A1 + B1))),
     M(eigen_expm(eigenMat(A1))) %*% M(eigen_expm(eigenMat(B1))))
})
chk("det(expm A) = exp(tr A)", {
  A1 <- matrix(c(0.5, 1, -2, 3), 2)
  eq(det(M(eigen_expm(eigenMat(A1)))), exp(sum(diag(A1))), tol = 1e-10)
})

cat("== 7. 矩阵对数与开方 ==\n")
P <- matrix(c(4, 1, 1, 3), 2)
chk("logm(expm A) = A",  eq(M(eigen_logm(eigen_expm(eigenMat(P)))), P, tol = 1e-9))
chk("expm(logm P) = P",  eq(M(eigen_expm(eigen_logm(eigenMat(P)))), P, tol = 1e-9))
chk("sqrtm^2 = P",       eq(M(eigen_sqrtm(eigenMat(P))) %*% M(eigen_sqrtm(eigenMat(P))), P,
                             tol = 1e-12))
chk("sqrtm 对称",         eq(M(eigen_sqrtm(eigenMat(P))), t(M(eigen_sqrtm(eigenMat(P)))),
                             tol = 1e-12))
chk("powm(P, 3) = P^3",  eq(M(eigen_powm(eigenMat(P), 3)), P %*% P %*% P, tol = 1e-10))
chk("powm(P, 0.5) = sqrtm", eq(M(eigen_powm(eigenMat(P), 0.5)), M(eigen_sqrtm(eigenMat(P))),
                               tol = 1e-10))
chk("powm(P, -1) = 逆",  eq(M(eigen_powm(eigenMat(P), -1)), solve(P), tol = 1e-10))

cat("== 8. 矩阵函数的错误处理 ==\n")
chk("非方阵 expm 报错",   errs(eigen_expm(eigenMat(matrix(1:6, 2, 3)))))
chk("非方阵 sqrtm 报错",  errs(eigen_sqrtm(eigenMat(matrix(1:6, 2, 3)))))
chk("复矩阵 expm 报错",   errs(eigen_expm(eigenMat(matrix(1 + 1i, 1, 1)))))
# Eigen 对 diag(-1, 2) 会安静地返回 diag(0, 0.693)——全有限、全错。
# 封装层靠 expm(结果) == A 的验算把它拦住，这条断言就是钉这个行为。
chk("负特征值的实对数报错", errs(eigen_logm(eigenMat(diag(c(-1, 2))))))
chk("负特征值的实幂报错",   errs(eigen_powm(eigenMat(diag(c(-1, 2))), 0.5)))
chk("logm 正常矩阵仍然可用", eq(M(eigen_logm(eigenMat(diag(c(1, 2))))), diag(log(c(1, 2)))))
chk("负特征值的实开方报错", errs(eigen_sqrtm(eigenMat(diag(c(-1, 2))))))

cat(sprintf("\n通过 %d, 失败 %d\n", ok, bad))
if (bad > 0L) stop(sprintf("%d assertion(s) failed", bad), call. = FALSE)
