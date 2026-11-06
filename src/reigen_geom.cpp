// reigen_geom.cpp —— Eigen Geometry 模块的 R 侧承载
//
// 与稠密矩阵一样，这些对象是 EXTPTRSXP，数据留在 C++ 侧，R 只拿指针。
// 差别在于它们没有「形状」：四元数是 4 个数、变换是 4x4 矩阵，但都不是
// R 意义上的矩阵，所以不挂 dim 属性，as.matrix() 才落回 R。
//
// 为什么不用数值向量表示四元数：Eigen 的 Quaterniond 有自己的运算
// （slerp、angularDistance、FromTwoVectors），这些都不是逐元素的；
// 一旦拆成 4 个 double，就只能靠调用方自己拼回去，等于没封装。
//
// 维度固定为 double / 3D。Eigen 的 Geometry 也支持 float 和 2D，但 R 的
// 数值只有 double，而 2D 情形用 AngleAxis / Rotation2D 已经够用；
// 为每种标量再展开一遍会让这一层的符号翻四倍，收益却接近零。

#include "reigen.h"
#include <Eigen/Geometry>

using namespace reigen;

namespace {

struct QuatObj {
  Eigen::Quaterniond q;
};

struct AAObj {
  Eigen::AngleAxisd aa;
};

struct TransObj {
  Eigen::Transform<double, 3, Eigen::Affine> t;
};

struct Rot2Obj {
  Eigen::Rotation2Dd r;
};

// R 侧的 3 维向量：接受长度 3 的向量，或 3x1 / 1x3 矩阵，或 3x3 里的一列。
// 列主序下 data()[0..2] 正好是前三个元素，所以三种形状都能直接取。
Eigen::Vector3d vec3(SEXP x, const char* what) {
  Mat m = from_sexp(x, what);
  if (m.cplx) Rcpp::stop("%s must be real, not complex", what);
  if (m.r.size() != 3)
    Rcpp::stop("%s must have exactly 3 elements, got %lld", what,
               (long long)m.r.size());
  Eigen::Vector3d v;
  for (int i = 0; i < 3; ++i) v[i] = m.r.data()[i];
  return v;
}

// n x 3 的点集，每行一个点。Eigen 的 umeyama() 要的是每列一个点，
// 所以调用方负责转置——这里只按 R 的直觉收 n x 3。
Eigen::MatrixXd points_nx3(SEXP x, const char* what) {
  Mat m = from_sexp(x, what);
  if (m.cplx) Rcpp::stop("%s must be real, not complex", what);
  if (m.r.cols() != 3)
    Rcpp::stop("%s must have 3 columns (one point per row), got %lld columns",
               what, (long long)m.r.cols());
  return m.r;
}

inline SEXP wrap_quat(const Eigen::Quaterniond& q) {
  QuatObj* p = new QuatObj();
  p->q = q;
  return wrap_obj(p, "eigenQuat");
}

inline SEXP wrap_aa(const Eigen::AngleAxisd& aa) {
  AAObj* p = new AAObj();
  p->aa = aa;
  return wrap_obj(p, "eigenAngleAxis");
}

inline SEXP wrap_trans(const Eigen::Transform<double, 3, Eigen::Affine>& t) {
  TransObj* p = new TransObj();
  p->t = t;
  return wrap_obj(p, "eigenTransform");
}

inline SEXP wrap_rot2(const Eigen::Rotation2Dd& r) {
  Rot2Obj* p = new Rot2Obj();
  p->r = r;
  return wrap_obj(p, "eigenRotation2D");
}

inline SEXP mat3_to_sexp(const Eigen::Matrix3d& m) {
  return to_sexp(make(Eigen::MatrixXd(m)));
}

}  // namespace

// ============================================================== 四元数

// [[Rcpp::export]]
SEXP reigen_quat_identity() { return wrap_quat(Eigen::Quaterniond::Identity()); }

// [[Rcpp::export]]
SEXP reigen_quat_new(SEXP wxyz) {
  Mat m = from_sexp(wxyz, "q");
  if (m.cplx || m.r.size() != 4)
    Rcpp::stop("q must be a real vector of length 4: c(w, x, y, z)");
  const double* p = m.r.data();
  Eigen::Quaterniond q(p[0], p[1], p[2], p[3]);
  if (q.norm() == 0.0)
    Rcpp::stop("q must not be the zero quaternion (its norm is 0)");
  return wrap_quat(q.normalized());
}

// [[Rcpp::export]]
SEXP reigen_quat_from_matrix(SEXP m) {
  Mat M = from_sexp(m, "m");
  if (M.cplx || M.r.rows() != 3 || M.r.cols() != 3)
    Rcpp::stop("m must be a real 3x3 rotation matrix");
  // 传进来的是旋转矩阵，但 R 用户可能给一个带缩放/反射的矩阵；Eigen 的
  // 转换假定正交，这里只做形状检查，剩下的交给它。
  // 必须是定长的 Matrix3d：Quaternion 从矩阵构造走的是
  // internal::quaternionbase_assign_impl，它只为 3x3 的定长矩阵特化，
  // 传 MatrixXd 进来会在模板深处报 "incomplete type"。
  Eigen::Matrix3d R = M.r;
  Eigen::Quaterniond q(R);
  q.normalize();
  return wrap_quat(q);
}

// [[Rcpp::export]]
SEXP reigen_quat_from_angle_axis(SEXP angle, SEXP axis) {
  const double a = Rcpp::as<double>(angle);
  Eigen::Vector3d ax = vec3(axis, "axis");
  if (ax.norm() == 0.0) Rcpp::stop("axis must not be the zero vector");
  return wrap_quat(Eigen::Quaterniond(Eigen::AngleAxisd(a, ax.normalized())));
}

// [[Rcpp::export]]
SEXP reigen_quat_coeffs(SEXP o) {
  QuatObj* p = obj_ptr<QuatObj>(o, "eigenQuat", "q");
  Rcpp::NumericVector out(4);
  out[0] = p->q.w();
  out[1] = p->q.x();
  out[2] = p->q.y();
  out[3] = p->q.z();
  out.names() = Rcpp::CharacterVector::create("w", "x", "y", "z");
  return out;
}

// [[Rcpp::export]]
SEXP reigen_quat_matrix(SEXP o) {
  QuatObj* p = obj_ptr<QuatObj>(o, "eigenQuat", "q");
  return mat3_to_sexp(p->q.toRotationMatrix());
}

// [[Rcpp::export]]
SEXP reigen_quat_conjugate(SEXP o) {
  QuatObj* p = obj_ptr<QuatObj>(o, "eigenQuat", "q");
  return wrap_quat(p->q.conjugate());
}

// [[Rcpp::export]]
SEXP reigen_quat_inverse(SEXP o) {
  QuatObj* p = obj_ptr<QuatObj>(o, "eigenQuat", "q");
  return wrap_quat(p->q.inverse());
}

// [[Rcpp::export]]
SEXP reigen_quat_normalized(SEXP o) {
  QuatObj* p = obj_ptr<QuatObj>(o, "eigenQuat", "q");
  if (p->q.norm() == 0.0) Rcpp::stop("cannot normalise a zero quaternion");
  return wrap_quat(p->q.normalized());
}

// [[Rcpp::export]]
SEXP reigen_quat_norm(SEXP o) {
  QuatObj* p = obj_ptr<QuatObj>(o, "eigenQuat", "q");
  return Rcpp::wrap(p->q.norm());
}

// [[Rcpp::export]]
SEXP reigen_quat_squared_norm(SEXP o) {
  QuatObj* p = obj_ptr<QuatObj>(o, "eigenQuat", "q");
  return Rcpp::wrap(p->q.squaredNorm());
}

// [[Rcpp::export]]
SEXP reigen_quat_dot(SEXP a, SEXP b) {
  QuatObj* pa = obj_ptr<QuatObj>(a, "eigenQuat", "a");
  QuatObj* pb = obj_ptr<QuatObj>(b, "eigenQuat", "b");
  return Rcpp::wrap(pa->q.dot(pb->q));
}

// 两个朝向之间的夹角（弧度）。Eigen 的 angularDistance() 返回 [0, pi]，
// 与「旋转角」不同：q 与 -q 表示同一个朝向，夹角为 0。
// [[Rcpp::export]]
SEXP reigen_quat_angular_distance(SEXP a, SEXP b) {
  QuatObj* pa = obj_ptr<QuatObj>(a, "eigenQuat", "a");
  QuatObj* pb = obj_ptr<QuatObj>(b, "eigenQuat", "b");
  return Rcpp::wrap(pa->q.angularDistance(pb->q));
}

// [[Rcpp::export]]
SEXP reigen_quat_slerp(SEXP a, SEXP b, SEXP t) {
  QuatObj* pa = obj_ptr<QuatObj>(a, "eigenQuat", "a");
  QuatObj* pb = obj_ptr<QuatObj>(b, "eigenQuat", "b");
  const double u = Rcpp::as<double>(t);
  if (u < 0.0 || u > 1.0) Rcpp::stop("t must lie in [0, 1], got %g", u);
  return wrap_quat(pa->q.slerp(u, pb->q));
}

// 四元数乘法 q1 * q2：先转 q2 再转 q1，与矩阵乘法同序。
// R 侧接成 %*%，这样 q1 %*% q2 的读法和 t1 %*% t2 一致。
// [[Rcpp::export]]
SEXP reigen_quat_compose(SEXP a, SEXP b) {
  QuatObj* pa = obj_ptr<QuatObj>(a, "eigenQuat", "a");
  QuatObj* pb = obj_ptr<QuatObj>(b, "eigenQuat", "b");
  return wrap_quat(pa->q * pb->q);
}

// [[Rcpp::export]]
SEXP reigen_quat_to_angle_axis(SEXP o) {  QuatObj* p = obj_ptr<QuatObj>(o, "eigenQuat", "q");
  return wrap_aa(Eigen::AngleAxisd(p->q));
}

// 用四元数旋转一个 3 维向量。Eigen 里是 q * v，返回旋转后的向量。
// [[Rcpp::export]]
SEXP reigen_quat_rotate(SEXP o, SEXP v) {
  QuatObj* p = obj_ptr<QuatObj>(o, "eigenQuat", "q");
  Eigen::Vector3d x = vec3(v, "v");
  Eigen::Vector3d y = p->q * x;
  Rcpp::NumericVector out(3);
  for (int i = 0; i < 3; ++i) out[i] = y[i];
  return out;
}

// 两个向量之间的最短旋转（FromTwoVectors）。示例里常用它把某个方向
// 对齐到另一个方向。
// [[Rcpp::export]]
SEXP reigen_quat_from_two_vectors(SEXP a, SEXP b) {
  Eigen::Vector3d x = vec3(a, "a");
  Eigen::Vector3d y = vec3(b, "b");
  if (x.norm() == 0.0 || y.norm() == 0.0)
    Rcpp::stop("both vectors must be non-zero");
  return wrap_quat(Eigen::Quaterniond::FromTwoVectors(x, y));
}

// ============================================================== 轴角

// [[Rcpp::export]]
SEXP reigen_angle_axis_new(SEXP angle, SEXP axis) {
  const double a = Rcpp::as<double>(angle);
  Eigen::Vector3d ax = vec3(axis, "axis");
  if (ax.norm() == 0.0) Rcpp::stop("axis must not be the zero vector");
  return wrap_aa(Eigen::AngleAxisd(a, ax.normalized()));
}

// [[Rcpp::export]]
SEXP reigen_angle_axis_angle(SEXP o) {
  AAObj* p = obj_ptr<AAObj>(o, "eigenAngleAxis", "aa");
  return Rcpp::wrap(p->aa.angle());
}

// [[Rcpp::export]]
SEXP reigen_angle_axis_axis(SEXP o) {
  AAObj* p = obj_ptr<AAObj>(o, "eigenAngleAxis", "aa");
  Rcpp::NumericVector out(3);
  for (int i = 0; i < 3; ++i) out[i] = p->aa.axis()[i];
  return out;
}

// [[Rcpp::export]]
SEXP reigen_angle_axis_matrix(SEXP o) {
  AAObj* p = obj_ptr<AAObj>(o, "eigenAngleAxis", "aa");
  return mat3_to_sexp(Eigen::Matrix3d(p->aa.toRotationMatrix()));
}

// [[Rcpp::export]]
SEXP reigen_angle_axis_to_quat(SEXP o) {
  AAObj* p = obj_ptr<AAObj>(o, "eigenAngleAxis", "aa");
  return wrap_quat(Eigen::Quaterniond(p->aa));
}

// ============================================================== 仿射变换

// [[Rcpp::export]]
SEXP reigen_transform_identity() {
  return wrap_trans(Eigen::Transform<double, 3, Eigen::Affine>::Identity());
}

// [[Rcpp::export]]
SEXP reigen_transform_from_matrix(SEXP m) {
  Mat M = from_sexp(m, "m");
  if (M.cplx || M.r.rows() != 4 || M.r.cols() != 4)
    Rcpp::stop("m must be a real 4x4 homogeneous matrix");
  Eigen::Transform<double, 3, Eigen::Affine> t;
  t.matrix() = M.r;
  return wrap_trans(t);
}

// [[Rcpp::export]]
SEXP reigen_transform_translation(SEXP v) {
  Eigen::Vector3d x = vec3(v, "v");
  Eigen::Transform<double, 3, Eigen::Affine> t =
      Eigen::Transform<double, 3, Eigen::Affine>::Identity();
  t.translate(x);
  return wrap_trans(t);
}

// [[Rcpp::export]]
SEXP reigen_transform_from_quat(SEXP q) {
  QuatObj* p = obj_ptr<QuatObj>(q, "eigenQuat", "q");
  Eigen::Transform<double, 3, Eigen::Affine> t =
      Eigen::Transform<double, 3, Eigen::Affine>::Identity();
  t.rotate(p->q);
  return wrap_trans(t);
}

// [[Rcpp::export]]
SEXP reigen_transform_from_angle_axis(SEXP aa) {
  AAObj* p = obj_ptr<AAObj>(aa, "eigenAngleAxis", "aa");
  Eigen::Transform<double, 3, Eigen::Affine> t =
      Eigen::Transform<double, 3, Eigen::Affine>::Identity();
  t.rotate(p->aa);
  return wrap_trans(t);
}

// [[Rcpp::export]]
SEXP reigen_transform_scaling(SEXP s) {
  Eigen::Vector3d x = vec3(s, "s");
  Eigen::Transform<double, 3, Eigen::Affine> t =
      Eigen::Transform<double, 3, Eigen::Affine>::Identity();
  t.scale(x);
  return wrap_trans(t);
}

// [[Rcpp::export]]
SEXP reigen_transform_matrix(SEXP o) {
  TransObj* p = obj_ptr<TransObj>(o, "eigenTransform", "t");
  Eigen::MatrixXd m = p->t.matrix();
  return to_sexp(make(m));
}

// [[Rcpp::export]]
SEXP reigen_transform_inverse(SEXP o) {
  TransObj* p = obj_ptr<TransObj>(o, "eigenTransform", "t");
  return wrap_trans(p->t.inverse());
}

// 变换复合：返回 a * b，即先作用 b 再作用 a（与矩阵乘法同序）。
// [[Rcpp::export]]
SEXP reigen_transform_compose(SEXP a, SEXP b) {
  TransObj* pa = obj_ptr<TransObj>(a, "eigenTransform", "a");
  TransObj* pb = obj_ptr<TransObj>(b, "eigenTransform", "b");
  return wrap_trans(pa->t * pb->t);
}

// n x 3 的点集整体变换，返回 n x 3。
// [[Rcpp::export]]
SEXP reigen_transform_points(SEXP o, SEXP pts) {
  TransObj* p = obj_ptr<TransObj>(o, "eigenTransform", "t");
  Eigen::MatrixXd P = points_nx3(pts, "pts");
  Eigen::MatrixXd out = P;
  for (Eigen::Index i = 0; i < P.rows(); ++i) {
    Eigen::Vector3d v(P(i, 0), P(i, 1), P(i, 2));
    Eigen::Vector3d w = p->t * v;
    for (int j = 0; j < 3; ++j) out(i, j) = w[j];
  }
  return to_sexp(make(out));
}

// 以下三个返回新对象而不是就地修改：R 的值语义要求 A 被改过之后，
// 原来的 A 不能跟着变。就地修改会让 `B <- t; t$translate(v)` 变成
// 两个名字指向同一个变换，与 R 用户的预期相反。
// [[Rcpp::export]]
SEXP reigen_transform_translate(SEXP o, SEXP v) {
  TransObj* p = obj_ptr<TransObj>(o, "eigenTransform", "t");
  Eigen::Transform<double, 3, Eigen::Affine> t = p->t;
  t.translate(vec3(v, "v"));
  return wrap_trans(t);
}

// [[Rcpp::export]]
SEXP reigen_transform_rotate_quat(SEXP o, SEXP q) {
  TransObj* p = obj_ptr<TransObj>(o, "eigenTransform", "t");
  QuatObj* pq = obj_ptr<QuatObj>(q, "eigenQuat", "q");
  Eigen::Transform<double, 3, Eigen::Affine> t = p->t;
  t.rotate(pq->q);
  return wrap_trans(t);
}

// [[Rcpp::export]]
SEXP reigen_transform_rotate_angle_axis(SEXP o, SEXP aa) {
  TransObj* p = obj_ptr<TransObj>(o, "eigenTransform", "t");
  AAObj* pa = obj_ptr<AAObj>(aa, "eigenAngleAxis", "aa");
  Eigen::Transform<double, 3, Eigen::Affine> t = p->t;
  t.rotate(pa->aa);
  return wrap_trans(t);
}

// [[Rcpp::export]]
SEXP reigen_transform_scale(SEXP o, SEXP s) {
  TransObj* p = obj_ptr<TransObj>(o, "eigenTransform", "t");
  Eigen::Transform<double, 3, Eigen::Affine> t = p->t;
  t.scale(vec3(s, "s"));
  return wrap_trans(t);
}

// [[Rcpp::export]]
SEXP reigen_transform_translation_of(SEXP o) {
  TransObj* p = obj_ptr<TransObj>(o, "eigenTransform", "t");
  Eigen::Vector3d v = p->t.translation();
  Rcpp::NumericVector out(3);
  for (int i = 0; i < 3; ++i) out[i] = v[i];
  return out;
}

// ============================================================== Rotation2D

// [[Rcpp::export]]
SEXP reigen_rotation2d_new(SEXP angle) {
  return wrap_rot2(Eigen::Rotation2Dd(Rcpp::as<double>(angle)));
}

// [[Rcpp::export]]
SEXP reigen_rotation2d_angle(SEXP o) {
  Rot2Obj* p = obj_ptr<Rot2Obj>(o, "eigenRotation2D", "r");
  return Rcpp::wrap(p->r.angle());
}

// [[Rcpp::export]]
SEXP reigen_rotation2d_matrix(SEXP o) {
  Rot2Obj* p = obj_ptr<Rot2Obj>(o, "eigenRotation2D", "r");
  return to_sexp(make(Eigen::MatrixXd(p->r.toRotationMatrix())));
}

// [[Rcpp::export]]
SEXP reigen_rotation2d_inverse(SEXP o) {
  Rot2Obj* p = obj_ptr<Rot2Obj>(o, "eigenRotation2D", "r");
  return wrap_rot2(p->r.inverse());
}

// ============================================================== umeyama

// 最小二乘意义下的相似变换估计：找 s、R、t 使 dst ≈ s·R·src + t。
// 输入按 R 的直觉收 n x 3（每行一个点），内部转成 Eigen 要的 3 x n。
// 这是 Eigen 的 umeyama()，与 OpenCV 的 estimateAffine3D 不同：它不剔除外点。
// [[Rcpp::export]]
SEXP reigen_umeyama(SEXP src, SEXP dst, SEXP with_scaling) {
  Eigen::MatrixXd S = points_nx3(src, "src");
  Eigen::MatrixXd D = points_nx3(dst, "dst");
  if (S.rows() != D.rows())
    Rcpp::stop("src and dst must have the same number of points (%lld vs %lld)",
               (long long)S.rows(), (long long)D.rows());
  if (S.rows() < 3)
    Rcpp::stop("umeyama needs at least 3 point pairs, got %lld",
               (long long)S.rows());
  const bool ws = Rcpp::as<bool>(with_scaling);
  Eigen::Matrix4d m = Eigen::umeyama(S.transpose(), D.transpose(), ws);
  Eigen::Transform<double, 3, Eigen::Affine> t;
  t.matrix() = m;
  return wrap_trans(t);
}
