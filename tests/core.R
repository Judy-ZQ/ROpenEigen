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
eq <- function(a, b, tol = 1e-10) isTRUE(all.equal(a, b, tolerance = tol))

cat("== 1. 构造与往返 ==\n")
m <- matrix(c(4, 1, 1, 3), 2)
A <- eigenMat(m)
chk("dim",            eq(dim(A), c(2L, 2L)))
chk("as.matrix 往返",  eq(as.matrix(A), m))
chk("length",         eq(length(A), 4L))
chk("is_eigen",       is_eigen(A))
chk("is.matrix",      is.matrix(A))
chk("is.numeric",     is.numeric(A))
chk("eigenMat(矩阵) 幂等", identical(eigenMat(A), A))
chk("向量 -> 列向量",  eq(as.matrix(eigenMat(1:3)), matrix(1:3, 3, 1)))
chk("eigenVec",       eq(dim(eigenVec(1:3)), c(3L, 1L)))
chk("eigenVec 类",     inherits(eigenVec(1:3), "eigenMat"))

cat("== 2. 下标 ==\n")
B <- eigenMat(matrix(1:12, 3, 4))
chk("B[2,3]",   eq(B[2, 3], 8))
chk("B[1,]",    eq(B[1, ], c(1, 4, 7, 10)))
chk("B[,2]",    eq(B[, 2], c(4, 5, 6)))
chk("B[1:2,3:4]", eq(B[1:2, 3:4], matrix(1:12, 3, 4)[1:2, 3:4]))
chk("B[-1,]",   eq(B[-1, ], matrix(1:12, 3, 4)[-1, ]))
B2 <- B; B2[2, 2] <- 99
chk("B[2,2] <- 99", eq(B2[2, 2], 99))
chk("原对象未变",      eq(B[2, 2], 5))
B3 <- B; B3[1, ] <- 0
chk("整行赋值",        eq(B3[1, ], c(0, 0, 0, 0)))

# 线性下标：x[k] 取的是按列展平后的第 k 个元素，不是第 k 行。这一组之前完全
# 没有覆盖，于是 E[1] 静默返回了第一行、E[5] 直接报越界、diag(E) 跟着挂。
M34 <- matrix(1:12, 3, 4); E34 <- eigenMat(M34)
chk("E[1] 是元素",       eq(E34[1], M34[1]))
chk("E[5] 是元素",       eq(E34[5], M34[5]))
chk("E[c(1,5)]",         eq(E34[c(1, 5)], M34[c(1, 5)]))
chk("E[1:12] 展平",      eq(E34[1:12], M34[1:12]))
chk("E[13] 越界为 NA",   eq(E34[13], M34[13]))
chk("E[0] 为空",         eq(E34[0], M34[0]))
chk("E[cbind(1:2,1:2)]", eq(E34[cbind(1:2, 1:2)], M34[cbind(1:2, 1:2)]))
chk("diag(E34)",         eq(diag(E34), diag(M34)))
chk("diag 方阵",         eq(diag(eigenMat(matrix(c(1,2,3,4), 2))), c(1, 4)))
E34b <- eigenMat(M34); E34b[5] <- 99
chk("E[5] <- 99 写元素",  eq(as.matrix(E34b), { z <- M34; z[5] <- 99; z }))
E34c <- eigenMat(M34); E34c[3] <- 0
chk("E[3] <- 0 写元素",   eq(as.matrix(E34c), { z <- M34; z[3] <- 0; z }))
chk("diag<- 赋值",        eq(as.matrix({ z <- eigenMat(matrix(0, 2, 2)); diag(z) <- c(7, 8); z }),
                             { z <- matrix(0, 2, 2); diag(z) <- c(7, 8); z }))
chk("x[] <- 0",           eq(as.matrix({ z <- eigenMat(M34); z[] <- 0; z }), M34 * 0))

cat("== 3. 算术 ==\n")
chk("A + 1",     eq(as.matrix(A + 1), m + 1))
chk("1 + A",     eq(as.matrix(1 + A), m + 1))
chk("A + m",     eq(as.matrix(A + m), m + m))
chk("A * A",     eq(as.matrix(A * A), m * m))
chk("A / 2",     eq(as.matrix(A / 2), m / 2))
chk("A ^ 2",     eq(as.matrix(A^2), m^2))
chk("-A",        eq(as.matrix(-A), -m))
chk("A == m",    eq(A == m, m == m))
chk("A > 2",     eq(A > 2, m > 2))
chk("%*%",       eq(as.matrix(A %*% A), m %*% m))
chk("非方阵 %*%", eq(as.matrix(B %*% t(B)), matrix(1:12, 3, 4) %*% t(matrix(1:12, 3, 4))))

cat("== 4. Math / Summary ==\n")
chk("abs",    eq(as.matrix(abs(eigenMat(-2:2))), matrix(abs(-2:2), 5, 1)))
chk("sqrt",   eq(as.matrix(sqrt(eigenMat(c(1, 4, 9)))), matrix(c(1, 2, 3), 3, 1)))
chk("exp",    eq(as.matrix(exp(eigenMat(1))), matrix(exp(1), 1, 1)))
chk("log(8,2)", eq(as.vector(log(eigenMat(8), base = 2)), 3))
chk("sum",    eq(sum(B), sum(1:12)))
chk("prod",   eq(prod(eigenMat(c(1, 2, 3))), 6))
chk("min/max", eq(c(min(B), max(B)), c(1, 12)))
chk("mean",   eq(mean(B), mean(1:12)))
chk("any/all", eq(c(any(eigenMat(c(0, 1))), all(eigenMat(c(1, 2)))), c(TRUE, TRUE)))

cat("== 5. 线性代数 ==\n")
chk("t",          eq(as.matrix(t(B)), t(matrix(1:12, 3, 4))))
chk("solve",      eq(solve(A, c(1, 2)), solve(m, c(1, 2))))
chk("solve 无 b",  eq(as.matrix(solve(A)), solve(m)))
# eigen_solve 对两种入参都要能用：eigenMat 走一次性求解，分解对象走 $solve。
# 之前这里零覆盖，而 eigenMat 也带 "eigenObject" 类，于是 eigen_solve(A, b)
# 掉进 $.eigenObject 报 "eigenMat has no method 'solve'"。
chk("eigen_solve 矩阵",   eq(as.vector(as.matrix(eigen_solve(A, c(1, 2)))), solve(m, c(1, 2))))
chk("eigen_solve 分解",   eq(as.vector(as.matrix(eigen_solve(eigen_llt(A), c(1, 2)))),
                             solve(m, c(1, 2))))
chk("eigen_det",  eq(eigen_det(A), det(m)))
# determinant() 的 modulus 是 log|det|，不是 |det|。这里必须对照 base R 的
# determinant()，而不是 det()：写成 modulus == det(m) 会让「返回 |det|」这个
# bug 一路绿灯（历史上正是如此），而且 det(E) 会变成 exp(det(A))。
chk("determinant modulus", eq(determinant(A)$modulus, determinant(m)$modulus))
chk("determinant 类", identical(class(determinant(A)), "det"))
chk("determinant sign", eq(determinant(A)$sign, 1))
chk("det() 走 determinant", eq(det(A), det(m)))
# 负行列式：sign 必须是 -1，否则 det() 和 determinant() 一起错。
Dneg <- eigenMat(matrix(c(0, 1, 1, 0), 2))
chk("determinant 负号", eq(determinant(Dneg)$sign, -1))
chk("det() 负行列式",  eq(det(Dneg), det(matrix(c(0, 1, 1, 0), 2))))
# 奇异：log|det| = -Inf，det() 应为 0 而不是 NaN。
Dsing <- eigenMat(matrix(c(1, 2, 2, 4), 2))
chk("determinant 奇异 modulus", is.infinite(determinant(Dsing)$modulus))
chk("det() 奇异为 0", eq(det(Dsing), 0))
# 大矩阵：det 溢出，logdet 不能。base R 的 determinant() 同样不溢出。
set.seed(42)
Big <- eigenMat(matrix(rnorm(400), 20, 20) + diag(20))
chk("logdet 不溢出", is.finite(eigen_logdet(Big)))
chk("logdet 对照 base R",
    eq(eigen_logdet(Big), as.vector(determinant(as.matrix(Big))$modulus), 1e-8))
chk("eigen_norm F", eq(eigen_norm(A, "F"), norm(m, "F")))
chk("eigen_norm 2", eq(eigen_norm(A, "2"), norm(m, "2")))
chk("eigen_trace", eq(eigen_trace(A), sum(diag(m))))
chk("eigen_inverse", eq(as.matrix(eigen_inverse(A)), solve(m)))
chk("eigen_rank",  eq(eigen_rank(B), 2))
chk("eigen_dot",   eq(eigen_dot(eigenVec(1:3), eigenVec(1:3)), 14))
chk("eigen_dot 标量", identical(length(eigen_dot(eigenVec(1:3), eigenVec(1:3))), 1L) &&
                       is.null(dim(eigen_dot(eigenVec(1:3), eigenVec(1:3)))))
chk("eigen_dot 不复共轭", eq(eigen_dot(eigenVec(c(1i, 1)), eigenVec(c(1, 1))), 1 + 1i))
chk("eigen_dot 长度不符报错", inherits(try(eigen_dot(eigenVec(1:3), eigenVec(1:2)), silent = TRUE), "try-error"))
chk("eigen_cross", eq(as.vector(as.matrix(eigen_cross(eigenVec(c(1,0,0)), eigenVec(c(0,1,0))))), c(0,0,1)))

cat("== 6. 生成器 ==\n")
chk("identity",  eq(as.matrix(eigen_identity(3)), diag(3)))
chk("zeros",     eq(as.matrix(eigen_zeros(2, 3)), matrix(0, 2, 3)))
chk("ones",      eq(as.matrix(eigen_ones(2, 3)), matrix(1, 2, 3)))
chk("linspace",  eq(as.vector(as.matrix(eigen_linspace(0, 1, 5))), seq(0, 1, length.out = 5)))
chk("diag 取",    eq(as.vector(as.matrix(eigen_diag(A))), diag(m)))
set.seed(42); r1 <- as.matrix(eigen_random(3, 3))
set.seed(42); r2 <- as.matrix(eigen_random(3, 3))
chk("random 可复现", eq(r1, r2))
chk("random 范围",  all(r1 >= 0 & r1 <= 1))
chk("rbind",     eq(as.matrix(eigen_rbind(A, A)), rbind(m, m)))
chk("cbind",     eq(as.matrix(eigen_cbind(A, A)), cbind(m, m)))
chk("reshape",   eq(as.matrix(eigen_reshape(B, 4, 3)), matrix(1:12, 4, 3)))
chk("rep",       eq(as.matrix(eigen_rep(A, 2, 2)), kronecker(matrix(1, 2, 2), m)))

cat("== 7. 复数 ==\n")
z <- matrix(c(1+1i, 2-1i, 0+2i, 3), 2)
Z <- eigenMat(z)
chk("复数构造",  isTRUE(eigen_is_complex(Z)))
chk("复数往返",  eq(as.matrix(Z), z))
chk("Re",       eq(as.matrix(Re(Z)), Re(z)))
chk("Im",       eq(as.matrix(Im(Z)), Im(z)))
chk("Mod",      eq(as.matrix(Mod(Z)), Mod(z)))
chk("Conj",     eq(as.matrix(Conj(Z)), Conj(z)))
chk("复数 %*%",  eq(as.matrix(Z %*% Z), z %*% z))

cat("== 8. 分解 ==\n")
S <- matrix(c(4, 2, 2, 3), 2)          # 对称正定
s <- eigenMat(S)
llt <- eigen_llt(s)
chk("llt is_ok",   isTRUE(llt$is_ok()))
chk("llt$solve",   eq(as.vector(as.matrix(llt$solve(c(1, 2)))), solve(S, c(1, 2))))
chk("llt 无参",     eq(as.matrix(llt$solve()), solve(S)))
chk("llt$matrixL", eq(as.matrix(llt$matrixL()) %*% t(as.matrix(llt$matrixL())), S))
chk("llt$logdet",  eq(llt$logdet(), log(det(S))))
chk("solve(llt,b)", eq(solve(llt, c(1, 2)), solve(S, c(1, 2))))
chk("names()",      all(c("solve", "matrixL") %in% names(llt)))

ldlt <- eigen_ldlt(s)
chk("ldlt$solve",  eq(as.vector(as.matrix(ldlt$solve(c(1, 2)))), solve(S, c(1, 2))))
chk("ldlt 重建",    eq(as.matrix(ldlt$matrixL()) %*% diag(as.vector(as.matrix(ldlt$vectorD()))) %*% t(as.matrix(ldlt$matrixL())), S))

lu <- eigen_lu(s)
chk("lu$solve",    eq(as.vector(as.matrix(lu$solve(c(1, 2)))), solve(S, c(1, 2))))
chk("lu$det",      eq(as.vector(lu$det()), det(S)))
chk("lu$inverse",  eq(as.matrix(lu$inverse()), solve(S)))

flu <- eigen_full_lu(B)
chk("fulllu$rank", eq(flu$rank(), 2))
# FullPivLU 的 det/inverse 只对方阵有定义（Eigen 内部 assert，NDEBUG 下越界），
# 长方形上必须报错而不是崩。rcond 走 SVD，长方形可用。
chk("fulllu$det 拒绝长方形",     inherits(try(flu$det(),     silent = TRUE), "try-error"))
chk("fulllu$inverse 拒绝长方形", inherits(try(flu$inverse(), silent = TRUE), "try-error"))
chk("fulllu$rcond 长方形可用",   eq(flu$rcond(), eigen_rcond(B), tol = 1e-10))
chk("fulllu$rcond 方阵",         eigen_full_lu(s)$rcond() > 0)
chk("eigen_rcond 长方形 = smin/smax",
    eq(eigen_rcond(B), { dd <- svd(matrix(1:12, 3, 4), nu = 0, nv = 0)$d; dd[3] / dd[1] }, tol = 1e-8))
chk("eigen_rcond 方阵 = 1/kappa", eq(eigen_rcond(s), 1 / kappa(S, exact = TRUE), tol = 1e-8))
chk("eigen_rcond 奇异 = 0",      eigen_rcond(eigenMat(matrix(c(1, 2, 2, 4), 2))) < 1e-12)

qr <- eigen_qr(B)
chk("qr$matrixQ 正交", eq(as.matrix(t(qr$matrixQ()) %*% qr$matrixQ()), diag(3)))
chk("qr 重建",      eq(as.matrix(qr$matrixQ() %*% qr$matrixR()), matrix(1:12, 3, 4), tol = 1e-8))

cp <- eigen_colpiv_qr(B)
chk("colpiv rank",  eq(cp$rank(), 2))
chk("colpiv 重建",   eq(as.matrix(cp$matrixQ() %*% cp$matrixR() %*% t(cp$permutationP())), matrix(1:12, 3, 4), tol = 1e-8))

sv <- eigen_svd(B)
chk("svd 奇异值",    eq(as.vector(as.matrix(sv$singular_values())), svd(matrix(1:12,3,4), nu=0, nv=0)$d, tol = 1e-8))
chk("svd 重建",      eq(as.matrix(sv$reconstructed()), matrix(1:12, 3, 4), tol = 1e-8))
B3 <- matrix(c(1, 2, 3, 4, 5, 6), 3, 2)
chk("svd$solve 超定", eq(as.vector(as.matrix(eigen_svd(eigenMat(B3))$solve(c(1,2,3)))), qr.solve(B3, c(1,2,3)), tol = 1e-8))
chk("eigen_pinv",    eq(as.matrix(eigen_pinv(B)), MASS::ginv(matrix(1:12,3,4)), tol = 1e-8))
chk("pinv 满足 Penrose", max(abs(matrix(1:12,3,4) %*% as.matrix(eigen_pinv(B)) %*% matrix(1:12,3,4) - matrix(1:12,3,4))) < 1e-8)

ev <- eigen_eigen(s)
chk("eigen 值",      eq(sort(as.vector(as.matrix(ev$values()))), sort(eigen(S, symmetric=TRUE, only.values=TRUE)$values), tol = 1e-8))
chk("eigen 重建",    eq(as.matrix(ev$vectors() %*% diag(as.vector(as.matrix(ev$values()))) %*% t(ev$vectors())), S, tol = 1e-8))
chk("eigen$solve",   eq(as.vector(as.matrix(ev$solve(c(1, 2)))), solve(S, c(1, 2)), tol = 1e-8))

G <- matrix(c(0, -1, 1, 0), 2)          # 旋转矩阵，特征值 ±i
ge <- eigen_geigen(eigenMat(G))
chk("geigen 是复数",  isTRUE(eigen_is_complex(ge$values())))
chk("geigen 值",      eq(sort(Im(as.vector(as.matrix(ge$values())))), sort(c(-1, 1)), tol = 1e-8))

cat("== 9. print ==\n")
chk("print eigenMat",  is.character(capture.output(print(A))))
chk("print 分解",       is.character(capture.output(print(llt))))

cat(sprintf("\n通过 %d, 失败 %d\n", ok, bad))
if (bad > 0L) stop(sprintf("%d assertion(s) failed", bad), call. = FALSE)
