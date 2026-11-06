# 几何对象：四元数、轴角、仿射变换、二维旋转
#
# 与 eigenMat 一样，这些对象是 C++ 侧的引用对象（EXTPTRSXP），R 只拿指针。
# 它们没有 dim 属性——四元数不是 4x1 矩阵，变换也不是 4x4 矩阵的别名，
# 硬给它们挂 dim 会让 is.array() 之类的东西产生误解。要数据就 as.matrix()。
#
# 方法用 $ 取用（q$slerp(r, 0.5)、t$inverse()），与分解对象同一套人机接口；
# 同时每个方法都有一个同名的 eigen_ 前缀函数，方便在管道里用。

# ---------------------------------------------------------------- 四元数

#' Quaternions, axis-angle rotations and affine transforms
#'
#' Thin R wrappers around Eigen's `Geometry` module.  Rotations are kept on the
#' C++ side as reference objects, so slerp, composition and inverse are the
#' Eigen implementations rather than hand-rolled R arithmetic.
#'
#' [`eigen_quat()`] takes the coefficients in the conventional `(w, x, y, z)`
#' order and normalises them; a zero quaternion is an error, since it does not
#' represent a rotation.  `eigen_quat_from_matrix()` accepts any 3x3 matrix and
#' projects it onto the nearest rotation, as Eigen does.
#'
#' The transform constructors all return a **new** object: `t$translate(v)`
#' does not modify `t`, because R users expect value semantics and a shared
#' mutable transform would be a trap.
#'
#' @param angle Rotation angle in radians.
#' @param axis A length-3 vector; it is normalised, and a zero axis is an error.
#' @param m For `eigen_quat_from_matrix()`, a 3x3 rotation matrix.  For
#'   `eigen_transform()`, a 4x4 homogeneous matrix.
#' @param a,b For `eigen_quat_slerp()`, the two end points; for
#'   `eigen_quat_from_two_vectors()`, two non-zero 3-vectors.
#' @param t Interpolation parameter in `[0, 1]`.
#' @param x A geometry object (or, for `eigen_transform_points()`, the transform).
#' @param v A length-3 translation / scaling vector or a point to rotate.
#' @param pts An `n x 3` matrix of points, one per row.
#' @param src,dst `n x 3` point sets with `n >= 3`; row `i` of `src` is matched
#'   with row `i` of `dst`.
#' @param scaling Whether `eigen_umeyama()` may include a scale factor.
#' @return A `eigenQuat`, `eigenAngleAxis`, `eigenTransform` or `eigenRotation2D`
#'   object.  Use `as.matrix()` to get the rotation/homogeneous matrix, or
#'   `names()` to list the available `$` methods.
#' @name eigen-geometry
NULL

#' Quaternion from its four coefficients
#'
#' Builds the quaternion `w + x i + y j + z k`.  The coefficients are
#' normalised, so any non-zero quadruple is accepted; a zero quaternion is an
#' error, since it does not represent a rotation.
#'
#' This lives in its own help page only because `w`, `x`, `y` and `z` cannot
#' share a page with the other constructors, whose first argument is also called
#' `x` -- `R CMD check` rejects the duplicated `\argument` entry.
#'
#' @param w,x,y,z Quaternion coefficients, `(w, x, y, z)`.
#' @return An `eigenQuat` object.
#' @seealso [`eigen-geometry`] for the rest of the `Geometry` wrappers.
#' @examples
#' q <- eigen_quat(cos(pi / 4), 0, 0, sin(pi / 4))
#' eigen_quat_coeffs(q)
#' @export
eigen_quat <- function(w = 1, x = 0, y = 0, z = 0) {
  reigen_quat_new(c(w, x, y, z))
}

#' @rdname eigen-geometry
#' @export
eigen_quat_from_matrix <- function(m) reigen_quat_from_matrix(.as_eigen(m))

#' @rdname eigen-geometry
#' @export
eigen_quat_from_angle_axis <- function(angle, axis) {
  reigen_quat_from_angle_axis(angle, axis)
}

#' @rdname eigen-geometry
#' @export
eigen_quat_from_two_vectors <- function(a, b) reigen_quat_from_two_vectors(a, b)

#' @rdname eigen-geometry
#' @export
eigen_quat_coeffs <- function(x) reigen_quat_coeffs(x)

#' @rdname eigen-geometry
#' @export
eigen_quat_matrix <- function(x) reigen_quat_matrix(x)

#' @rdname eigen-geometry
#' @export
eigen_quat_slerp <- function(a, b, t) reigen_quat_slerp(a, b, t)

#' @rdname eigen-geometry
#' @export
eigen_quat_rotate <- function(x, v) reigen_quat_rotate(x, v)

#' @rdname eigen-geometry
#' @export
eigen_quat_compose <- function(a, b) reigen_quat_compose(a, b)

#' @rdname eigen-geometry
#' @export
eigen_quat_inverse <- function(x) reigen_quat_inverse(x)

#' @rdname eigen-geometry
#' @export
eigen_quat_conjugate <- function(x) reigen_quat_conjugate(x)

#' @rdname eigen-geometry
#' @export
eigen_quat_normalize <- function(x) reigen_quat_normalized(x)

#' @rdname eigen-geometry
#' @export
eigen_quat_dot <- function(a, b) reigen_quat_dot(a, b)

#' @rdname eigen-geometry
#' @export
eigen_quat_angular_distance <- function(a, b) reigen_quat_angular_distance(a, b)

#' @rdname eigen-geometry
#' @export
eigen_quat_to_angle_axis <- function(x) reigen_quat_to_angle_axis(x)

#' @export
print.eigenQuat <- function(x, ...) {
  cf <- reigen_quat_coeffs(x)
  cat(sprintf("<eigenQuat> w=%g x=%g y=%g z=%g  (norm %g)\n",
              cf[1], cf[2], cf[3], cf[4], reigen_quat_norm(x)))
  invisible(x)
}

#' @export
as.matrix.eigenQuat <- function(x, ...) reigen_quat_matrix(x)

# ---------------------------------------------------------------- 轴角

#' @rdname eigen-geometry
#' @export
eigen_angle_axis <- function(angle, axis) reigen_angle_axis_new(angle, axis)

#' @rdname eigen-geometry
#' @export
eigen_angle_axis_matrix <- function(x) reigen_angle_axis_matrix(x)

#' @export
print.eigenAngleAxis <- function(x, ...) {
  cat(sprintf("<eigenAngleAxis> angle=%g rad, axis=(%s)\n",
              reigen_angle_axis_angle(x),
              paste(format(reigen_angle_axis_axis(x), digits = 4), collapse = ", ")))
  invisible(x)
}

#' @export
as.matrix.eigenAngleAxis <- function(x, ...) reigen_angle_axis_matrix(x)

# ---------------------------------------------------------------- 仿射变换

#' @rdname eigen-geometry
#' @export
eigen_transform <- function(m = NULL) {
  if (is.null(m)) return(reigen_transform_identity())
  if (inherits(m, "eigenTransform")) return(m)
  if (inherits(m, "eigenQuat")) return(reigen_transform_from_quat(m))
  if (inherits(m, "eigenAngleAxis")) return(reigen_transform_from_angle_axis(m))
  reigen_transform_from_matrix(.as_eigen(m))
}

#' @rdname eigen-geometry
#' @export
eigen_translation <- function(v) reigen_transform_translation(v)

#' @rdname eigen-geometry
#' @export
eigen_scaling <- function(v) reigen_transform_scaling(v)

#' @rdname eigen-geometry
#' @export
eigen_transform_matrix <- function(x) reigen_transform_matrix(x)

#' @rdname eigen-geometry
#' @export
eigen_transform_inverse <- function(x) reigen_transform_inverse(x)

#' @rdname eigen-geometry
#' @export
eigen_transform_points <- function(x, pts) reigen_transform_points(x, pts)

#' @rdname eigen-geometry
#' @export
eigen_umeyama <- function(src, dst, scaling = TRUE) {
  reigen_umeyama(src, dst, isTRUE(scaling))
}

#' @export
print.eigenTransform <- function(x, ...) {
  m <- reigen_transform_matrix(x)
  cat("<eigenTransform> 4 x 4 affine\n")
  print(m, digits = 6)
  invisible(x)
}

#' @export
as.matrix.eigenTransform <- function(x, ...) reigen_transform_matrix(x)

# q1 %*% q2：与矩阵乘法的顺序一致，先作用 q2 再作用 q1。
# 不提供 *.eigenQuat：四元数还有标量乘这一层歧义，%*% 是唯一不会误读的写法。
#' @export
`%*%.eigenQuat` <- function(x, y) {
  if (!inherits(y, "eigenQuat")) {
    stop("%*% on an eigenQuat needs another eigenQuat on the right",
         call. = FALSE)
  }
  reigen_quat_compose(x, y)
}

# t1 %*% t2：与矩阵乘法的顺序一致，先作用 t2 再作用 t1。
#' @export
`%*%.eigenTransform` <- function(x, y) {
  if (!inherits(y, "eigenTransform")) {
    stop("%*% on an eigenTransform needs another eigenTransform on the right",
         call. = FALSE)
  }
  reigen_transform_compose(x, y)
}

# ---------------------------------------------------------------- Rotation2D

#' @rdname eigen-geometry
#' @export
eigen_rotation2d <- function(angle) reigen_rotation2d_new(angle)

#' @export
print.eigenRotation2D <- function(x, ...) {
  cat(sprintf("<eigenRotation2D> angle=%g rad\n", reigen_rotation2d_angle(x)))
  invisible(x)
}

#' @export
as.matrix.eigenRotation2D <- function(x, ...) reigen_rotation2d_matrix(x)

# ---------------------------------------------------------------- $ 方法表

.reigen_reg("eigenQuat", m = list(
  coeffs          = function(self) reigen_quat_coeffs(self),
  matrix          = function(self) reigen_quat_matrix(self),
  inverse         = function(self) reigen_quat_inverse(self),
  conjugate       = function(self) reigen_quat_conjugate(self),
  normalize       = function(self) reigen_quat_normalized(self),
  norm            = function(self) reigen_quat_norm(self),
  squaredNorm     = function(self) reigen_quat_squared_norm(self),
  dot             = function(self, other) reigen_quat_dot(self, other),
  angularDistance = function(self, other) reigen_quat_angular_distance(self, other),
  slerp           = function(self, other, t) reigen_quat_slerp(self, other, t),
  rotate          = function(self, v) reigen_quat_rotate(self, v),
  toAngleAxis     = function(self) reigen_quat_to_angle_axis(self),
  toTransform     = function(self) reigen_transform_from_quat(self)
))

.reigen_reg("eigenAngleAxis", m = list(
  angle       = function(self) reigen_angle_axis_angle(self),
  axis        = function(self) reigen_angle_axis_axis(self),
  matrix      = function(self) reigen_angle_axis_matrix(self),
  toQuaternion = function(self) reigen_angle_axis_to_quat(self),
  toTransform = function(self) reigen_transform_from_angle_axis(self)
))

.reigen_reg("eigenTransform", m = list(
  matrix      = function(self) reigen_transform_matrix(self),
  inverse     = function(self) reigen_transform_inverse(self),
  translation = function(self) reigen_transform_translation_of(self),
  translate   = function(self, v) reigen_transform_translate(self, v),
  rotate      = function(self, q) {
    if (inherits(q, "eigenAngleAxis")) reigen_transform_rotate_angle_axis(self, q)
    else reigen_transform_rotate_quat(self, q)
  },
  scale       = function(self, v) reigen_transform_scale(self, v),
  points      = function(self, pts) reigen_transform_points(self, pts),
  compose     = function(self, other) reigen_transform_compose(self, other)
))

.reigen_reg("eigenRotation2D", m = list(
  angle   = function(self) reigen_rotation2d_angle(self),
  matrix  = function(self) reigen_rotation2d_matrix(self),
  inverse = function(self) reigen_rotation2d_inverse(self)
))
