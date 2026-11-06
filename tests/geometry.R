suppressPackageStartupMessages(library(ROpenEigen))

# 几何模块（Eigen/Geometry）的对照测试。
#
# 这里的每一个断言都对着 Eigen 自己的语义写，而不是对着本包的实现写：
# 四元数是 (w, x, y, z)、角用弧度、q1 * q2 先转 q2、变换是 4x4 齐次矩阵。
# 只要有一条对不上，就说明封装层把 Eigen 的意思改掉了——那是比崩溃更坏的错，
# 因为它看起来能跑。

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

cat("== 1. 四元数：构造与系数 ==\n")
q <- eigen_quat_from_angle_axis(pi / 2, c(0, 0, 1))
chk("类名",            identical(class(q), c("eigenQuat", "eigenObject")))
chk("绕 z 轴 90 度的系数", eq(unname(eigen_quat_coeffs(q)), c(cos(pi / 4), 0, 0, sin(pi / 4))))
chk("单位四元数",       eq(unname(eigen_quat_coeffs(eigen_quat())), c(1, 0, 0, 0)))
chk("构造时归一化",     eq(unname(eigen_quat_coeffs(eigen_quat(2, 0, 0, 0))), c(1, 0, 0, 0)))
chk("norm",            eq(q$norm(), 1))
chk("squaredNorm",     eq(q$squaredNorm(), 1))
chk("零四元数是错误",   errs(eigen_quat(0, 0, 0, 0)))
chk("非数值输入是错误", errs(eigen_quat("a")))

cat("== 2. 四元数：旋转与矩阵 ==\n")
chk("转 x 轴得到 y 轴",  eq(as.numeric(eigen_quat_rotate(q, c(1, 0, 0))), c(0, 1, 0)))
chk("转 z 轴不变",       eq(as.numeric(eigen_quat_rotate(q, c(0, 0, 1))), c(0, 0, 1)))
chk("模长保持",          eq(sum(eigen_quat_rotate(q, c(1, 2, 3))^2), 14))
chk("旋转矩阵",          eq(as.matrix(q), matrix(c(0, 1, 0, -1, 0, 0, 0, 0, 1), 3, 3)))
chk("矩阵 -> 四元数往返", eq(abs(eigen_quat_coeffs(eigen_quat_from_matrix(as.matrix(q)))),
                             abs(eigen_quat_coeffs(q))))
chk("from_matrix 接受任意 3x3", is.finite(sum(as.matrix(
  eigen_quat_from_matrix(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 10), 3, 3))))))

cat("== 3. 四元数：复合、求逆、内积 ==\n")
r <- eigen_quat_from_angle_axis(pi / 2, c(1, 0, 0))
chk("逆 * 自身 = 1",   eq(unname(eigen_quat_coeffs(q$inverse() %*% q)), c(1, 0, 0, 0)))
chk("共轭 = 逆",       eq(unname(eigen_quat_coeffs(q$conjugate())), unname(eigen_quat_coeffs(q$inverse()))))
chk("dot(q, q) = 1",   eq(q$dot(q), 1))
chk("自己到自己的角距为 0", eq(q$angularDistance(q), 0))
chk("dot(q, -1) = -cos(45)", eq(q$dot(eigen_quat(-1, 0, 0, 0)), -sqrt(0.5)))
# q1 %*% q2 的旋转矩阵必须等于 R1 %*% R2：这是"顺序"唯一的硬证据
comp <- q %*% r
chk("复合矩阵 = 矩阵复合", eq(as.matrix(comp), as.matrix(q) %*% as.matrix(r)))
chk("复合两次 90 度 = 180 度", {
  two <- q %*% q
  eq(as.matrix(two), diag(c(-1, -1, 1)))
})
chk("与标量相乘是错误", errs(q %*% 2))

cat("== 4. 四元数：slerp 与轴角互转 ==\n")
mid <- eigen_quat_slerp(eigen_quat(), q, 0.5)
chk("slerp 到一半 = 45 度", eq(unname(eigen_quat_coeffs(mid)),
                               c(cos(pi / 8), 0, 0, sin(pi / 8))))
chk("slerp t=0 是起点",     eq(unname(eigen_quat_coeffs(eigen_quat_slerp(eigen_quat(), q, 0))),
                               c(1, 0, 0, 0)))
chk("slerp t=1 是终点",     eq(abs(eigen_quat_coeffs(eigen_quat_slerp(eigen_quat(), q, 1))),
                               abs(eigen_quat_coeffs(q))))
chk("slerp t 越界是错误",    errs(eigen_quat_slerp(eigen_quat(), q, 1.5)))
aa <- q$toAngleAxis()
chk("toAngleAxis 角度",      eq(aa$angle(), pi / 2))
chk("toAngleAxis 轴",        eq(abs(aa$axis()), c(0, 0, 1)))
chk("轴角转四元数一致",       eq(abs(eigen_quat_coeffs(eigen_angle_axis(pi / 2, c(0, 0, 1))$toQuaternion())),
                                abs(eigen_quat_coeffs(q))))
chk("两向量之间",            eq(as.numeric(eigen_quat_rotate(
  eigen_quat_from_two_vectors(c(1, 0, 0), c(0, 1, 0)), c(1, 0, 0))), c(0, 1, 0)))

cat("== 5. 轴角 ==\n")
a2 <- eigen_angle_axis(pi / 3, c(0, 0, 2))   # 轴会被归一化
chk("类名",        identical(class(a2), c("eigenAngleAxis", "eigenObject")))
chk("角度保留",     eq(a2$angle(), pi / 3))
chk("轴归一化",     eq(a2$axis(), c(0, 0, 1)))
chk("矩阵是旋转矩阵", eq(as.matrix(a2) %*% t(as.matrix(a2)), diag(3)))
chk("行列式为 1",    eq(det(as.matrix(a2)), 1))
chk("零轴是错误",    errs(eigen_angle_axis(1, c(0, 0, 0))))

cat("== 6. 仿射变换 ==\n")
I <- eigen_transform()
chk("类名",           identical(class(I), c("eigenTransform", "eigenObject")))
chk("单位变换",        eq(as.matrix(I), diag(4)))
t <- eigen_translation(c(1, 2, 3))
homog <- function(Rm = diag(3), tr = c(0, 0, 0)) {
  m <- diag(4); m[1:3, 1:3] <- Rm; m[1:3, 4] <- tr; m
}
chk("平移矩阵",        eq(as.matrix(t), homog(tr = c(1, 2, 3))))
chk("变换点",          eq(eigen_transform_points(t, rbind(c(0, 0, 0), c(1, 1, 1))),
                          rbind(c(1, 2, 3), c(2, 3, 4))))
chk("求逆",            eq(as.matrix(t$inverse()), homog(tr = c(-1, -2, -3))))
chk("值语义：t 未被修改", eq(as.matrix(t)[1:3, 4], c(1, 2, 3)))
chk("链式结果不回写",    eq(as.matrix(t$translate(c(1, 1, 1)))[1:3, 4], c(2, 3, 4)))
chk("缩放矩阵",         eq(as.matrix(eigen_scaling(c(1, 2, 3)))[1:3, 1:3], diag(c(1, 2, 3))))
chk("从四元数构造",      eq(as.matrix(eigen_transform(q))[1:3, 1:3], as.matrix(q)))
chk("从轴角构造",        eq(as.matrix(eigen_transform(a2))[1:3, 1:3], as.matrix(a2)))
chk("从 4x4 构造",       eq(as.matrix(eigen_transform(as.matrix(t))), as.matrix(t)))
chk("eigen_transform 幂等", identical(eigen_transform(t), t))
chk("rotate 后的 translate 走局部坐标系", {
  tt <- eigen_transform()$rotate(q)$translate(c(1, 0, 0))
  eq(as.matrix(tt)[1:3, 4], c(0, 1, 0))      # q 把 x 轴转到了 y 轴
})
chk("先 translate 后 rotate 则平移不变", {
  tt <- eigen_transform()$translate(c(1, 0, 0))$rotate(q)
  eq(as.matrix(tt)[1:3, 4], c(1, 0, 0))
})
chk("compose = 矩阵乘",  eq(as.matrix(t$compose(t)), as.matrix(t) %*% as.matrix(t)))
chk("t1 %*% t2",         eq(as.matrix(t %*% t), as.matrix(t) %*% as.matrix(t)))
chk("翻译点三维校验",     errs(eigen_transform_points(t, matrix(1:4, 2, 2))))

cat("== 7. Umeyama 对齐 ==\n")
src <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
Rd <- as.matrix(eigen_quat_from_angle_axis(pi / 2, c(0, 0, 1)))
dst <- t(Rd %*% t(src)) + matrix(c(1, 2, 3), 4, 3, byrow = TRUE)
fit <- eigen_umeyama(src, dst)
chk("还原旋转",   eq(as.matrix(fit)[1:3, 1:3], Rd, tol = 1e-8))
chk("还原平移",   eq(as.matrix(fit)[1:3, 4], c(1, 2, 3), tol = 1e-8))
chk("点对点吻合", eq(eigen_transform_points(fit, src), dst, tol = 1e-8))
# 带尺度：dst = 2 * R * src
dst2 <- 2 * t(Rd %*% t(src))
chk("含尺度时缩放为 2", eq(as.matrix(eigen_umeyama(src, dst2))[1, 1]^2 +
                            as.matrix(eigen_umeyama(src, dst2))[2, 1]^2 +
                            as.matrix(eigen_umeyama(src, dst2))[3, 1]^2, 4, tol = 1e-8))
chk("不允许尺度时缩放为 1", eq(sum(as.matrix(eigen_umeyama(src, dst2, scaling = FALSE))[1:3, 1]^2),
                               1, tol = 1e-8))
chk("点太少是错误",  errs(eigen_umeyama(src[1:2, ], dst[1:2, ])))
chk("维度不一致是错误", errs(eigen_umeyama(src, dst[, 1:2])))

cat("== 8. Rotation2D ==\n")
r2 <- eigen_rotation2d(pi / 2)
chk("类名",     identical(class(r2), c("eigenRotation2D", "eigenObject")))
chk("角度",     eq(r2$angle(), pi / 2))
chk("矩阵",     eq(as.matrix(r2), matrix(c(0, 1, -1, 0), 2, 2)))
chk("求逆矩阵", eq(as.matrix(r2$inverse()), t(as.matrix(r2))))

cat("== 9. $ 方法表与打印 ==\n")
chk("四元数方法名", identical(names(q), c("coeffs", "matrix", "inverse", "conjugate",
                                          "normalize", "norm", "squaredNorm", "dot",
                                          "angularDistance", "slerp", "rotate",
                                          "toAngleAxis", "toTransform")))
chk("变换方法名",   identical(names(I), c("matrix", "inverse", "translation",
                                          "translate", "rotate", "scale", "points",
                                          "compose")))
chk("未知方法是错误", errs(q$nope()))
chk("print 四元数",  is.character(capture.output(print(q))))
chk("print 轴角",    is.character(capture.output(print(a2))))
chk("print 变换",    is.character(capture.output(print(t))))
chk("print Rotation2D", is.character(capture.output(print(r2))))
chk("as.matrix 一致", identical(as.matrix(q), q$matrix()))

cat(sprintf("\n通过 %d, 失败 %d\n", ok, bad))
if (bad > 0L) stop(sprintf("%d assertion(s) failed", bad), call. = FALSE)
