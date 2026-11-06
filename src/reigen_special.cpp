// reigen_special.cpp —— 特殊函数与矩阵函数
//
// 这两块在 Eigen 里都是 unsupported/ 下的模块，但它们的实现质量与 Dense
// 同级（SpecialFunctions 内部还做了矢量化），只是不承诺 API 稳定。
// 本包随附了整棵 unsupported/，所以这里直接把它们接出来。
//
// 特殊函数一律是逐元素的（Eigen 的 Array API），输入是矩阵也只是把它
// 摊平成 n 个数；输出保持原来的形状。这和 R 的 sin()/gamma() 一致。
//
// 复数不支持：Eigen 的这些函数只对实数标量有实现（betainc 之类甚至只
// 支持 double）。强行对 complex 编译会在模板深处报一串看不懂的错误，
// 不如在这里给一句能看懂的话。

#include "reigen.h"

#include <unsupported/Eigen/SpecialFunctions>
#include <unsupported/Eigen/MatrixFunctions>

using namespace reigen;

namespace {

// 形状保持不变地把一元函数作用到每个元素上。
// 直接复用 m.r 的列主序存储：MatrixXd 与 ArrayXd 的内存布局相同，
// Map 过来就是零拷贝。
template <class F>
Mat map_real(const Mat& m, F&& f) {
  if (m.cplx)
    Rcpp::stop("this function is defined for real matrices only; "
               "take Re() or Im() first");
  const Eigen::Index nr = m.r.rows(), nc = m.r.cols();
  const Eigen::Index n = m.r.size();
  Eigen::Map<const Eigen::ArrayXd> a(m.r.data(), n);
  Eigen::ArrayXd out(n);
  for (Eigen::Index i = 0; i < n; ++i) out[i] = f(a[i]);
  Mat r;
  r.r.resize(nr, nc);
  if (n > 0) std::memcpy(r.r.data(), out.data(), (size_t)n * sizeof(double));
  return r;
}

inline double num1(SEXP x, const char* what) {
  Rcpp::NumericVector v(x);
  if (v.size() != 1) Rcpp::stop("%s must be a single number", what);
  return v[0];
}

}  // namespace

// ============================================================ 一元特殊函数

// code: 0 erf, 1 erfc, 2 lgamma, 3 digamma, 4 ndtri
//       10..13 i0 i0e i1 i1e, 14..17 k0 k0e k1 k1e,
//       18..21 j0 j1 y0 y1
// [[Rcpp::export]]
SEXP reigen_special(SEXP x, int code) {
  Mat X = from_sexp(x, "x");
  switch (code) {
    case 0:  return wrap_mat(map_real(X, [](double v) { return std::erf(v); }));
    case 1:  return wrap_mat(map_real(X, [](double v) { return std::erfc(v); }));
    case 2:  return wrap_mat(map_real(X, [](double v) { return std::lgamma(v); }));
    case 3:  return wrap_mat(map_real(X, [](double v) { return Eigen::numext::digamma(v); }));
    case 4:  return wrap_mat(map_real(X, [](double v) { return Eigen::numext::ndtri(v); }));
    case 10: return wrap_mat(map_real(X, [](double v) { return Eigen::numext::bessel_i0(v); }));
    case 11: return wrap_mat(map_real(X, [](double v) { return Eigen::numext::bessel_i0e(v); }));
    case 12: return wrap_mat(map_real(X, [](double v) { return Eigen::numext::bessel_i1(v); }));
    case 13: return wrap_mat(map_real(X, [](double v) { return Eigen::numext::bessel_i1e(v); }));
    case 14: return wrap_mat(map_real(X, [](double v) { return Eigen::numext::bessel_k0(v); }));
    case 15: return wrap_mat(map_real(X, [](double v) { return Eigen::numext::bessel_k0e(v); }));
    case 16: return wrap_mat(map_real(X, [](double v) { return Eigen::numext::bessel_k1(v); }));
    case 17: return wrap_mat(map_real(X, [](double v) { return Eigen::numext::bessel_k1e(v); }));
    case 18: return wrap_mat(map_real(X, [](double v) { return Eigen::numext::bessel_j0(v); }));
    case 19: return wrap_mat(map_real(X, [](double v) { return Eigen::numext::bessel_j1(v); }));
    case 20: return wrap_mat(map_real(X, [](double v) { return Eigen::numext::bessel_y0(v); }));
    case 21: return wrap_mat(map_real(X, [](double v) { return Eigen::numext::bessel_y1(v); }));
    default: Rcpp::stop("unknown special-function code %d", code);
  }
}

// ============================================================ 二元特殊函数

// code: 0 igamma(a, x), 1 igammac(a, x), 2 zeta(x, q), 3 polygamma(n, x)
// [[Rcpp::export]]
SEXP reigen_special2(SEXP a, SEXP b, int code) {
  Mat A = from_sexp(a, "a");
  Mat B = from_sexp(b, "b");
  if (A.cplx || B.cplx) Rcpp::stop("special functions are defined for real input only");

  // 标量广播，与 Ops 组泛型的口径一致：一边是 1x1 时重复到另一边的形状。
  const bool as = (A.r.size() == 1), bs = (B.r.size() == 1);
  Eigen::Index nr, nc;
  if (as && bs)      { nr = 1;        nc = 1; }
  else if (as)       { nr = B.r.rows(); nc = B.r.cols(); }
  else if (bs)       { nr = A.r.rows(); nc = A.r.cols(); }
  else {
    if (A.r.rows() != B.r.rows() || A.r.cols() != B.r.cols())
      Rcpp::stop("non-conformable operands: %lldx%lld vs %lldx%lld",
                 (long long)A.r.rows(), (long long)A.r.cols(),
                 (long long)B.r.rows(), (long long)B.r.cols());
    nr = A.r.rows(); nc = A.r.cols();
  }

  const Eigen::Index n = nr * nc;
  Eigen::ArrayXd x(n), y(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    x[i] = A.r.data()[as ? 0 : i];
    y[i] = B.r.data()[bs ? 0 : i];
  }
  Eigen::ArrayXd out(n);
  switch (code) {
    case 0: out = Eigen::igamma(x, y);  break;
    case 1: out = Eigen::igammac(x, y); break;
    case 2: out = Eigen::zeta(x, y);    break;
    case 3: out = Eigen::polygamma(x, y); break;
    default: Rcpp::stop("unknown special-function code %d", code);
  }
  Mat r;
  r.r.resize(nr, nc);
  if (n > 0) std::memcpy(r.r.data(), out.data(), (size_t)n * sizeof(double));
  return wrap_mat(r);
}

// betainc 是三元函数，单独开口：Eigen 的签名是 (a, b, x)，
// 容易与 R 的 pbeta(x, a, b) 弄反，所以参数名写死成 a/b/x。
// [[Rcpp::export]]
SEXP reigen_betainc(SEXP a, SEXP b, SEXP x) {
  const double A = num1(a, "a");
  const double B = num1(b, "b");
  Mat X = from_sexp(x, "x");
  return wrap_mat(map_real(X, [&](double v) { return Eigen::numext::betainc(A, B, v); }));
}

// ============================================================ 矩阵函数

// code: 0 exp, 1 log, 2 sqrt, 3 pow(p)
// 这些不是逐元素的：exp(A) 是矩阵指数 sum A^k/k!，不是每个元素取 exp。
// 名字里带 m 是刻意的提醒。
//
// 这里不信任 Eigen 的返回值，而是**回头验算一遍**——原因是实测出来的：
// diag(-1, 2) 的实对数应该是复矩阵，Eigen 却安静地返回 diag(0, 0.693)，
// 全是有限数，看不出任何异常。只查 allFinite() 拦不住这种答案。
// 所以 log 验 expm(out) == A，sqrt 验 out^2 == A，非整数幂验
// expm(p * logm(A)) 且 logm(A) 自身先过一遍验算。整数幂直接连乘对照，
// 那是最便宜也最硬的检查。

static double rel_diff(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y) {
  // 相对 Frobenius 误差，用 max(1, ||Y||) 归一：A 接近 0 时不把绝对误差
  // 放大成相对误差。
  const double denom = std::max(1.0, Y.norm());
  return (X - Y).norm() / denom;
}

static Eigen::MatrixXd int_power(const Eigen::MatrixXd& A, long k) {
  Eigen::MatrixXd r = Eigen::MatrixXd::Identity(A.rows(), A.cols());
  for (long i = 0; i < k; ++i) r = r * A;
  return r;
}

// [[Rcpp::export]]
SEXP reigen_matrix_function(SEXP x, int code, double p) {
  Mat X = from_sexp(x, "x");
  if (X.cplx)
    Rcpp::stop("matrix functions are implemented for real matrices only");
  if (X.r.rows() != X.r.cols())
    Rcpp::stop("matrix functions need a square matrix, got %lldx%lld",
               (long long)X.r.rows(), (long long)X.r.cols());
  const Eigen::MatrixXd A = X.r;
  Eigen::MatrixXd out;

  const char* no_real =
      " -- the matrix probably has a negative real eigenvalue (or a complex "
      "eigenvalue pair), for which no real log/sqrt/power exists; use a "
      "different function, or split the matrix into real Schur blocks.";

  switch (code) {
    case 0:
      out = A.exp();          // 指数对任何实矩阵都有定义，不需要验算
      break;
    case 1: {
      out = A.log();
      if (!out.allFinite() || rel_diff(out.exp(), A) > 1e-7)
        Rcpp::stop("logm did not verify: expm(logm(A)) differs from A%s", no_real);
      break;
    }
    case 2: {
      out = A.sqrt();
      if (!out.allFinite() || rel_diff(out * out, A) > 1e-7)
        Rcpp::stop("sqrtm did not verify: sqrtm(A)^2 differs from A%s", no_real);
      break;
    }
    case 3: {
      out = A.pow(p);
      if (!out.allFinite())
        Rcpp::stop("powm produced non-finite entries%s", no_real);
      const double rounded = std::floor(p + 0.5);
      if (std::fabs(p - rounded) < 1e-12 && std::fabs(p) <= 64.0) {
        const long k = (long)std::fabs(rounded);
        if (rounded >= 0.0) {
          if (rel_diff(out, int_power(A, k)) > 1e-8)
            Rcpp::stop("powm(A, %g) did not verify against repeated "
                       "multiplication%s", p, no_real);
        } else {
          const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(A.rows(), A.cols());
          if (rel_diff(out * int_power(A, k), I) > 1e-8)
            Rcpp::stop("powm(A, %g) did not verify against the inverse%s", p, no_real);
        }
      } else {
        // 非整数幂：Eigen 内部就是 expm(p * logm(A))。既然 log 本身可能
        // 悄悄给出错答案，就先自己算一遍并验算，再比对结果。
        const Eigen::MatrixXd L = A.log();
        if (!L.allFinite() || rel_diff(L.exp(), A) > 1e-7)
          Rcpp::stop("powm(A, %g) needs a real logarithm, but A has none%s",
                     p, no_real);
        if (rel_diff(out, (p * L).exp()) > 1e-6)
          Rcpp::stop("powm(A, %g) did not verify against expm(p * logm(A))%s",
                     p, no_real);
      }
      break;
    }
    default:
      Rcpp::stop("unknown matrix-function code %d", code);
  }
  if (!out.allFinite())
    Rcpp::stop("matrix function produced non-finite entries%s", no_real);
  return wrap_mat(make(out));
}
