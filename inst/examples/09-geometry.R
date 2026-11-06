# ===========================================================================
# 09 · 几何：旋转、四元数、仿射变换
#
# 对应 EigenExample：
#   snippets/AngleAxis_mimic_euler.cpp
#   demos/opengl/trackball.cpp          （轨迹球：把 2-D 拖动映射成旋转）
#   demos/opengl/quaternion_demo.cpp    （四元数插值、欧拉角↔四元数）
#   demos/opengl/camera.h               （相机：绕目标旋转 / 本地旋转 / 平移）
#
# 运行：Rscript 09-geometry.R
#
# 这三个 demo 都是 OpenGL 程序，Windows/CI 上跑不起来（要 GLUT/Qt），但里面的
# 数学是纯粹的 Eigen 几何运算，与绘图无关。下面把数学部分逐段搬过来：拖动
# 轨迹球、四元数 slerp、相机绕点旋转，结果都是可验证的数字。
# ===========================================================================

suppressPackageStartupMessages(library(ROpenEigen))

# ---- AngleAxis_mimic_euler.cpp ----------------------------------------------
# cpp:  三个 AngleAxis 依次相乘得到旋转矩阵，再检查 isUnitary()
m <- eigen_angle_axis_matrix(eigen_angle_axis(0.25 * pi, c(1, 0, 0))) %*%
     eigen_angle_axis_matrix(eigen_angle_axis(0.5 * pi, c(0, 1, 0))) %*%
     eigen_angle_axis_matrix(eigen_angle_axis(0.33 * pi, c(0, 0, 1)))
print(round(as.matrix(m), 6))
cat("is unitary:", all(abs(t(as.matrix(m)) %*% as.matrix(m) - diag(3)) < 1e-6), "\n")

# 用四元数复合写同一件事，结果相同（q1 %*% q2 与 R1 %*% R2 同序）：
q <- eigen_quat_from_angle_axis(0.25 * pi, c(1, 0, 0)) %*%
     eigen_quat_from_angle_axis(0.5 * pi, c(0, 1, 0)) %*%
     eigen_quat_from_angle_axis(0.33 * pi, c(0, 0, 1))
cat("四元数与矩阵的差:",
    max(abs(as.matrix(eigen_quat_matrix(q)) - as.matrix(m))), "\n")

# ---- trackball.cpp：mapToSphere ---------------------------------------------
# cpp 把屏幕坐标 (x, y) 映射到单位球面：x、y 各自按视口归一化后取正弦，
# z 用 sqrt(1 - x² - y²)，超出球外则为 0。这段与 OpenGL 无关，直接照搬。
map_to_sphere <- function(px, py, w, h) {
  x <- (px - 0.5 * w) / w
  y <- (0.5 * h - py) / h
  sx <- sin(pi * x * 0.5)
  sy <- sin(pi * y * 0.5)
  s2 <- sx^2 + sy^2
  c(x = sx, y = sy, z = if (s2 < 1) sqrt(1 - s2) else 0)
}

w <- 800; h <- 600
p1 <- map_to_sphere(400, 300, w, h)     # 视口中心 → 球顶
p2 <- map_to_sphere(600, 300, w, h)     # 向右拖
cat("球心处:", round(p1, 6), " 拖动后:", round(p2, 6), "\n")

# ---- trackball.cpp：由两点求旋转 --------------------------------------------
# cpp:
#   axis  = last.cross(new).normalized();
#   angle = 2 * acos(last.dot(new));
#   相机绕目标转 q，或是本地转 -q
track <- function(last, new) {
  axis <- eigen_cross(eigenVec(last), eigenVec(new))
  n <- sqrt(sum(as.vector(axis)^2))
  if (n < 1e-12) return(NULL)                    # 两点几乎重合，不转
  axis <- as.vector(axis) / n
  cos_angle <- sum(last * new)
  if (abs(cos_angle) >= 1) return(NULL)
  angle <- 2 * acos(cos_angle)
  list(axis = axis, angle = angle)
}

r <- track(p1, p2)
cat("旋转轴:", round(r$axis, 6), " 角度:", round(r$angle, 6), "\n")

# 验证：把 p1 绕这个轴转这个角度，应当落到 p2 上
q <- eigen_quat_from_angle_axis(r$angle, r$axis)
cat("转到 p2 的误差:", max(abs(eigen_quat_rotate(q, p1) - p2)), "\n")

# ---- quaternion_demo.cpp：四元数 slerp --------------------------------------
# cpp:  Quaternionf lerp(float t, const Quaternionf& a, const Quaternionf& b)
#         { return a.slerp(t, b); }
# R:    eigen_quat_slerp(a, b, t)，t 是 [0, 1] 的插值参数。
qa <- eigen_quat_from_angle_axis(0, c(0, 0, 1))            # 单位旋转
qb <- eigen_quat_from_angle_axis(pi / 2, c(0, 0, 1))       # 绕 z 转 90°
for (t in c(0, 0.25, 0.5, 0.75, 1)) {
  qt <- eigen_quat_slerp(qa, qb, t)
  aa <- eigen_quat_to_angle_axis(qt)
  cat(sprintf("t=%.2f  角度=%.6f  轴=%s\n", t, aa$angle(),
              paste(round(aa$axis(), 6), collapse = " ")))
}

# slerp 走的是大圆，匀速；四元数点乘得到的夹角就是两姿态之间的角距离：
cat("角距离:", eigen_quat_angular_distance(qa, qb), "（应为 π/2）\n")

# ---- quaternion_demo.cpp：欧拉角 → 四元数 ----------------------------------
# cpp 把 (绕 Y 的 dx, 绕 X 的 -dy) 两个 AngleAxis 乘起来。R 里同样写法：
euler_to_quat <- function(dx, dy) {
  eigen_quat_from_angle_axis(dx * pi, c(0, 1, 0)) %*%
    eigen_quat_from_angle_axis(-dy * pi, c(1, 0, 0))
}
q <- euler_to_quat(0.1, 0.2)
print(round(as.matrix(eigen_quat_matrix(q)), 6))

# ---- camera.h：绕目标旋转 / 本地旋转 ---------------------------------------
# cpp 的 rotateAroundTarget(q)：把相机位置相对目标点做旋转；
# localRotate(q)：在相机自身坐标系里转。两者的差别只在参考系。
target <- c(0, 0, 0)
eye    <- c(0, 0, 5)

rot <- eigen_quat_from_angle_axis(pi / 6, c(0, 1, 0))     # 绕 Y 转 30°
around <- eigen_quat_rotate(rot, eye - target) + target    # 绕目标转
cat("绕目标转到:", round(around, 6), "\n")
cat("到目标的距离不变:", sqrt(sum((around - target)^2)), "\n")

# ---- 变换对象：t$rotate(q)$translate(v) 的语义 -----------------------------
# ROpenEigen 的 eigenTransform 是引用对象，$ 方法返回**新**对象（值语义）。
# 注意 t$rotate(q)$translate(v) 的平移发生在**旋转后的坐标系**里，
# 与 Eigen::Transform::translate() 一致：存下来的平移量是 R·v 而不是 v。
t <- eigen_transform()
t2 <- t$rotate(eigen_quat_from_angle_axis(pi / 2, c(0, 0, 1)))$translate(c(1, 0, 0))
print(round(as.matrix(t2), 6))
cat("平移量被旋转成了 R·v =", round(t2$translation(), 6), "\n")

# ---- umeyama：两组点的最小二乘刚体对齐 -------------------------------------
# cpp 里没有对应片段，但它是 Geometry 模块提供的关键算法，这里补一个可验证的例子。
src <- matrix(c(0, 0, 0,
                1, 0, 0,
                0, 1, 0,
                0, 0, 1), 4, 3, byrow = TRUE)
Rtrue <- as.matrix(eigen_rotation2d(pi / 3))
R3 <- diag(3); R3[1:2, 1:2] <- Rtrue
dst <- sweep(src %*% t(R3), 2, c(1, 2, 3), "+")      # 旋转 + 平移

u <- eigen_umeyama(src, dst)
cat("恢复的旋转误差:", max(abs(as.matrix(u)[1:3, 1:3] - R3)),
    " 平移:", round(as.matrix(u)[1:3, 4], 6), "\n")

# ===========================================================================
# 尾声：两个 OpenGL demo 里真正的 Eigen 部分
# ---- demos/opengl/icosphere.cpp / icosphere.h -----------------------------
# cpp 用 `Map<Vector3f>(vdata[i])` 把 12 组 GLfloat[3] 直接当成向量用（零拷贝），
# 再递归 _subdivide() 把 20 个三角面细分到球面网格。
# R 没有裸内存视图（见 02 章 Map 一节），但"顶点表 + 归一化"这一步是纯几何，
# 可以直接算，顺带验证正二十面体的棱长与顶点数。
phi <- (1 + sqrt(5)) / 2                       # 黄金比例，c 里写死成宏 X/Z
ico <- matrix(c(-1,  phi, 0,   1,  phi, 0,
                -1, -phi, 0,   1, -phi, 0,
                 0, -1,  phi,  0,  1,  phi,
                 0, -1, -phi,  0,  1, -phi,
                 phi, 0, -1,   phi, 0,  1,
                -phi, 0, -1,  -phi, 0,  1), 12, 3, byrow = TRUE)
ico <- ico / sqrt(rowSums(ico^2))              # 逐行归一化 → 落在单位球上
cat("顶点数:", nrow(ico), "\n")
cat("全部落在单位球上（最大偏差）:",
    format(max(abs(sqrt(rowSums(ico^2)) - 1)), scientific = TRUE), "\n")
cat("棱长（相邻两顶点距离）:", round(sqrt(sum((ico[1, ] - ico[6, ])^2)), 6), "\n")
cat("levels 层细分后的三角面数 20*4^levels:", sapply(0:3, function(k) 20 * 4^k), "\n")

# 一层细分：每条边取中点再归一化，面数 x4。这正是 _subdivide() 的数学。
subdivide <- function(v, faces) {
  mids <- t(apply(faces, 1, function(f) {
    m <- (v[f[1], ] + v[f[2], ]) / 2
    m / sqrt(sum(m^2))
  }))
  rbind(v, mids)
}
faces <- matrix(c(1, 6, 2), 1, 3)              # 取一个面做示范
v2 <- subdivide(ico, faces)
cat("细分后新增顶点也在单位球上:", max(abs(sqrt(rowSums(v2^2)) - 1)) < 1e-15, "\n")
cat("新顶点到两端点等距:",
    round(sqrt(sum((v2[13, ] - ico[1, ])^2)), 6) ==
    round(sqrt(sum((v2[13, ] - ico[6, ])^2)), 6), "\n")

# ---- demos/opengl/gpuhelper.cpp / gpuhelper.h -----------------------------
# cpp 里 Eigen 只出现在两处：`typedef Vector4f Color` 和
# `pushMatrix(const Matrix4f&, GLenum)`（转成 float[16] 交给 glLoadMatrixf），
# 其余全是 glPushMatrix / glOrtho 这类 GL 状态机调用。R 侧对应的只有
# "造一个 4x4 变换矩阵"这一半：eigen_transform() 的 as.matrix() 就是 4x4。
M4 <- as.matrix(eigen_transform())
cat("4x4 变换矩阵的对角线:", diag(M4), "\n")
cat("identity 分支:", identical(dim(as.matrix(eigen_identity(4))), c(4L, 4L)), "\n")

cat("\n09-geometry.R 结束\n")
