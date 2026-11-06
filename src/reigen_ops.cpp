// reigen_ops.cpp —— 逐元素运算、归约、范数、求解
//
// 分派靠 op 整数码而不是「一个函数一个符号」：Eigen 的算子都是模板，
// 每个都要实例化两份（double / complex），符号数量会爆炸。一个 switch
// 的开销相对于一次 n^2 的遍历可以忽略。
//
// 码表见 R 侧 R/ops.R 的注释，两边必须一致。

#include "reigen.h"

using namespace reigen;

// ============================================================ 逐元素二元

// op: 0 +, 1 -, 2 *, 3 /, 4 ^, 5 pmin, 6 pmax, 7 atan2, 8 hypot
template <class S>
static Mat cwise2_t(const Mat& a, const Mat& b, int op) {
  typedef Eigen::Matrix<S, -1, -1> M;
  typedef Eigen::Array<S, -1, -1>  A_;
  M X = to_eigen<S>(a);
  M Y = to_eigen<S>(b);

  // 标量广播要两边都判。R 里 1 + A 和 A + 1 都是合法的，只判 b
  // 会把 1 + A 判成「1x1 与 2x2 不齐」。
  const bool as = (X.size() == 1);
  const bool bs = (Y.size() == 1);
  Eigen::Index nr, nc;
  if (as && bs)      { nr = 1;         nc = 1; }
  else if (as)       { nr = Y.rows();  nc = Y.cols(); }
  else if (bs)       { nr = X.rows();  nc = X.cols(); }
  else {
    if (X.rows() != Y.rows() || X.cols() != Y.cols())
      Rcpp::stop("non-conformable elementwise operands: %lldx%lld vs %lldx%lld",
                 (long long)X.rows(), (long long)X.cols(),
                 (long long)Y.rows(), (long long)Y.cols());
    nr = X.rows(); nc = X.cols();
  }
  if (as && (X.rows() != nr || X.cols() != nc)) X = M::Constant(nr, nc, X.data()[0]);
  if (bs && (Y.rows() != nr || Y.cols() != nc)) Y = M::Constant(nr, nc, Y.data()[0]);

  A_ xa = X.array(), ya = Y.array();
  A_ ra;
  switch (op) {
    case 0: ra = xa + ya; break;
    case 1: ra = xa - ya; break;
    case 2: ra = xa * ya; break;
    case 3: ra = xa / ya; break;
    case 4: ra = xa.pow(ya); break;
    // min/max 需要全序，std::complex 没有。用 if constexpr 把这两支
    // 从复数实例化里彻底剔除——否则即使运行期走不到，模板展开时
    // 仍会去找 operator<，编译直接失败。
    case 5:
      if constexpr (std::is_same<S, double>::value) ra = xa.min(ya);
      else Rcpp::stop("pmin is not defined for complex numbers");
      break;
    case 6:
      if constexpr (std::is_same<S, double>::value) ra = xa.max(ya);
      else Rcpp::stop("pmax is not defined for complex numbers");
      break;
    default: Rcpp::stop("operator %d is not defined for this type", op);
  }
  return make(M(ra));
}

// atan2 / hypot 走循环：Eigen 的二元自定义函数要用 binaryExpr 拼，
// 代码比一次线性遍历还长，而这两个函数也不在热路径上。
static Mat atan2_loop(const Mat& a, const Mat& b) {
  Eigen::MatrixXd X = as_real(a, "atan2"), Y;
  if (nelem(b) == 1) Y = Eigen::MatrixXd::Constant(X.rows(), X.cols(), as_real(b, "atan2").data()[0]);
  else {
    Y = as_real(b, "atan2");
    if (Y.rows() != X.rows() || Y.cols() != X.cols()) Rcpp::stop("non-conformable atan2 operands");
  }
  Eigen::MatrixXd R(X.rows(), X.cols());
  for (Eigen::Index i = 0; i < X.size(); ++i)
    R.data()[i] = std::atan2(X.data()[i], Y.data()[i]);
  return make(R);
}

static Mat hypot_loop(const Mat& a, const Mat& b) {
  Eigen::MatrixXd X = as_real(a, "hypot"), Y;
  if (nelem(b) == 1) Y = Eigen::MatrixXd::Constant(X.rows(), X.cols(), as_real(b, "hypot").data()[0]);
  else {
    Y = as_real(b, "hypot");
    if (Y.rows() != X.rows() || Y.cols() != X.cols()) Rcpp::stop("non-conformable hypot operands");
  }
  Eigen::MatrixXd R(X.rows(), X.cols());
  for (Eigen::Index i = 0; i < X.size(); ++i)
    R.data()[i] = std::hypot(X.data()[i], Y.data()[i]);
  return make(R);
}

// [[Rcpp::export]]
SEXP reigen_cwise2(SEXP a, SEXP b, int op) {
  Mat A = from_sexp(a, "a");
  Mat B = from_sexp(b, "b");
  if (op == 7) return wrap_mat(atan2_loop(A, B));
  if (op == 8) return wrap_mat(hypot_loop(A, B));
  const bool cz = A.cplx || B.cplx;
  if (cz && (op == 5 || op == 6))
    Rcpp::stop("pmin/pmax are not defined for complex numbers");
  return wrap_mat(pick_cplx(cz, [&](auto tag) {
    typedef typename scalar_of<decltype(tag)>::type S;
    return cwise2_t<S>(A, B, op);
  }));
}

// ============================================================ 逐元素一元

// op: 0 abs, 1 sqrt, 2 exp, 3 log, 4 log10, 5 sin, 6 cos, 7 tan, 8 asin,
//     9 acos, 10 atan, 11 sinh, 12 cosh, 13 tanh, 14 round, 15 floor,
//     16 ceil, 17 sign, 18 square, 19 cube, 20 inverse, 21 rsqrt,
//     22 log1p, 23 expm1, 24 conj, 25 Re, 26 Im, 27 Arg, 28 abs2
//
// 复数输入下 0/26/27/28 返回实矩阵，25 返回实部，其余保持复数。

template <class S>
static Mat cwise1_real_t(const Eigen::Matrix<S, -1, -1>& X, int op) {
  typedef Eigen::Matrix<S, -1, -1> M;
  M R = X;
  switch (op) {
    case 0:  R = X.cwiseAbs(); break;
    case 1:  R = X.cwiseSqrt(); break;
    case 2:  R = X.array().exp().matrix(); break;
    case 3:  R = X.array().log().matrix(); break;
    case 4:  R = X.array().log10().matrix(); break;
    case 5:  R = X.array().sin().matrix(); break;
    case 6:  R = X.array().cos().matrix(); break;
    case 7:  R = X.array().tan().matrix(); break;
    case 8:  R = X.array().asin().matrix(); break;
    case 9:  R = X.array().acos().matrix(); break;
    case 10: R = X.array().atan().matrix(); break;
    case 11: R = X.array().sinh().matrix(); break;
    case 12: R = X.array().cosh().matrix(); break;
    case 13: R = X.array().tanh().matrix(); break;
    case 14: R = X.array().round().matrix(); break;
    case 15: R = X.array().floor().matrix(); break;
    case 16: R = X.array().ceil().matrix(); break;
    case 17: R = X.array().sign().matrix(); break;
    case 18: R = X.array().square().matrix(); break;
    case 19: R = X.array().cube().matrix(); break;
    case 20: R = X.array().inverse().matrix(); break;
    case 21: R = X.array().rsqrt().matrix(); break;
    case 22: R = X.array().log1p().matrix(); break;
    case 23: R = X.array().expm1().matrix(); break;
    case 24: break;  // conj: 实矩阵不变
    case 25: break;  // Re:   实矩阵不变
    case 28: R = X.array().square().matrix(); break;
    default: Rcpp::stop("operation %d is not defined for real matrices", op);
  }
  return make(R);
}

static Mat cwise1_cplx(const Mat& m, int op) {
  Eigen::MatrixXcd X = as_cplx(m);
  switch (op) {
    case 0:  return make(Eigen::MatrixXd(X.cwiseAbs()));
    case 28: return make(Eigen::MatrixXd(X.cwiseAbs2()));
    case 27: return make(Eigen::MatrixXd(X.array().arg().matrix()));
    case 25: return make(Eigen::MatrixXd(X.real()));
    case 26: return make(Eigen::MatrixXd(X.imag()));
    case 24: return make(Eigen::MatrixXcd(X.conjugate()));
    case 1:  return make(Eigen::MatrixXcd(X.cwiseSqrt()));
    case 2:  return make(Eigen::MatrixXcd(X.array().exp().matrix()));
    case 3:  return make(Eigen::MatrixXcd(X.array().log().matrix()));
    case 5:  return make(Eigen::MatrixXcd(X.array().sin().matrix()));
    case 6:  return make(Eigen::MatrixXcd(X.array().cos().matrix()));
    case 7:  return make(Eigen::MatrixXcd(X.array().tan().matrix()));
    case 10: return make(Eigen::MatrixXcd(X.array().atan().matrix()));
    case 11: return make(Eigen::MatrixXcd(X.array().sinh().matrix()));
    case 12: return make(Eigen::MatrixXcd(X.array().cosh().matrix()));
    case 13: return make(Eigen::MatrixXcd(X.array().tanh().matrix()));
    case 18: return make(Eigen::MatrixXcd(X.array().square().matrix()));
    case 19: return make(Eigen::MatrixXcd(X.array().cube().matrix()));
    case 20: return make(Eigen::MatrixXcd(X.array().inverse().matrix()));
    default:
      Rcpp::stop("operation %d is not defined for complex matrices", op);
  }
}

// [[Rcpp::export]]
SEXP reigen_cwise1(SEXP x, int op) {
  Mat X = from_sexp(x, "x");
  if (X.cplx) return wrap_mat(cwise1_cplx(X, op));
  return wrap_mat(cwise1_real_t<double>(X.r, op));
}

// ============================================================ 矩阵乘法

// [[Rcpp::export]]
SEXP reigen_matmul(SEXP a, SEXP b) {
  Mat A = from_sexp(a, "a");
  Mat B = from_sexp(b, "b");
  // 注意这里的 %%*%%：Rcpp::stop 走 tinyformat，会把 % 当格式符解析。
  // 直接写 %*% 的话 tinyformat 看到 %*（动态宽度）就抛
  // "Too many conversion specifiers in format string"，
  // 真正的诊断信息反而被这条内部错误盖掉——用户看到的是 tinyformat 的抱怨，
  // 而不是"non-conformable arguments"。要输出字面 % 必须写 %%。
  if (cols(A) != rows(B))
    Rcpp::stop("non-conformable arguments: %lldx%lld %%*%% %lldx%lld",
               (long long)rows(A), (long long)cols(A),
               (long long)rows(B), (long long)cols(B));
  return wrap_mat(pick_cplx(A.cplx || B.cplx, [&](auto tag) {
    typedef typename scalar_of<decltype(tag)>::type S;
    typedef Eigen::Matrix<S, -1, -1> M;
    M X = to_eigen<S>(A), Y = to_eigen<S>(B);
    return make(M(X * Y));
  }));
}

// ============================================================ 归约

// op: 0 sum, 1 prod, 2 mean, 3 min, 4 max, 5 sum of squares,
//     6 all, 7 any, 8 count of positives
// [[Rcpp::export]]
SEXP reigen_reduce(SEXP x, int op) {
  Mat X = from_sexp(x, "x");
  const Eigen::Index n = nelem(X);

  if (X.cplx) {
    Eigen::MatrixXcd A = X.c;
    switch (op) {
      case 0: return wrap_cplx(A.sum());
      case 1: return wrap_cplx(A.prod());
      case 2: return wrap_cplx(n ? A.sum() / (double)n : cplx_t(NA_REAL, 0));
      case 5: return Rcpp::wrap(A.cwiseAbs2().sum());
      case 6: return Rcpp::wrap(A.cwiseAbs2().sum() != 0.0);
      case 7: return Rcpp::wrap(A.cwiseAbs2().sum() != 0.0);
      default:
        if (op == 3 || op == 4)
          Rcpp::stop("min() and max() are not defined for complex matrices; "
                     "compare Mod(x) or Re(x) instead");
        Rcpp::stop("operation %d is not defined for complex matrices", op);
    }
  }

  const Eigen::MatrixXd& A = X.r;
  if (n == 0 && (op == 3 || op == 4))
    Rcpp::stop("min/max of an empty matrix");
  switch (op) {
    case 0: return Rcpp::wrap(A.sum());
    case 1: return Rcpp::wrap(A.prod());
    case 2: return Rcpp::wrap(n ? A.sum() / (double)n : NA_REAL);
    case 3: return Rcpp::wrap(A.minCoeff());
    case 4: return Rcpp::wrap(A.maxCoeff());
    case 5: return Rcpp::wrap(A.cwiseAbs2().sum());
    case 6: return Rcpp::wrap(A.size() == 0 ? true : (A.array() != 0.0).all());
    case 7: return Rcpp::wrap(A.size() > 0 ? (A.array() != 0.0).any() : false);
    case 8: return Rcpp::wrap((double)(A.array() > 0.0).count());
    default: Rcpp::stop("unknown reduction code %d", op);
  }
}

// 返回 list(value, row, col)，行列 1-based。
//
// 不用 Eigen 的 minCoeff(&r,&c)，因为它和 base 的 which.min() 在 NA 上不一致。
// minCoeff 把第一个元素当作初始"最优"，而它的比较器 numext::mini(b,a) 等价于
// (b < a) ? b : a —— 初始 a 是 NaN 时，后面每个 b < NaN 都是 false，NaN 就被
// 一直留着。实测：3x3 全 1、把其中一个位置换成 NA，只有 NA 在列主序第一位
// 那一种情况返回 NA@(1,1)，其余 8 种都正常跳过。也就是说这是个"只在第一个
// 元素上发作"的 bug，而 which.min(c(NA,1)) 是 2。
//
// 所以自己扫一遍，列主序、跳过 NA/NaN、严格小于才替换，第一个出现的位置就是
// 答案。Inf 不算缺失：which.min(c(Inf,1)) 是 2，which.min(c(-Inf,1)) 是 1，
// 所以用 ISNAN() 而不是 !R_FINITE()。
//
// 全部缺失时 base 的 which.min() 返回 integer(0)，list 表达不了；这里返回
// value/row/col 全 NA，至少不会像以前那样谎报一个位置。
// [[Rcpp::export]]
SEXP reigen_which(SEXP x, int op) {
  Mat X = from_sexp(x, "x");
  if (X.cplx) Rcpp::stop("min/max index is not defined for complex matrices");
  const Eigen::MatrixXd& A = X.r;
  if (A.size() == 0) Rcpp::stop("min/max of an empty matrix");
  const bool want_max = (op == 1);
  bool found = false;
  double best = 0.0;
  Eigen::Index br = 0, bc = 0;
  for (Eigen::Index c = 0; c < A.cols(); ++c)
    for (Eigen::Index r = 0; r < A.rows(); ++r) {
      const double v = A(r, c);
      if (ISNAN(v)) continue;
      if (!found || (want_max ? v > best : v < best)) {
        best = v; br = r; bc = c; found = true;
      }
    }
  if (!found)
    return Rcpp::List::create(Rcpp::Named("value") = NA_REAL,
                              Rcpp::Named("row")   = NA_INTEGER,
                              Rcpp::Named("col")   = NA_INTEGER);
  return Rcpp::List::create(Rcpp::Named("value") = best,
                            Rcpp::Named("row")   = (int)br + 1,
                            Rcpp::Named("col")   = (int)bc + 1);
}

// ============================================================ 范数

// Frobenius 范数，先按最大元素缩放再求平方和。
//
// 朴素的 sqrt(sum(x^2)) 有两个真实的坏结局，而且都不报错：
//   元素 1e-200 -> 返回 0      （平方下溢；真值是 1.414e-200）
//   元素 1e300  -> 返回 Inf    （平方上溢；真值是 1.414e+300）
// 返回 0 尤其糟——它看起来完全正常，而 0 和 1.4e-200 在后续计算里天差地别。
//
// base 的 norm(x, "F") 走 LAPACK 的 dlange/zlange，本来就是缩放的，所以那两个
// 量级它都是对的。这里用同样的做法：记下 max|x|，累加 (x/max)^2，最后乘回
// max。O(n)，不分配临时矩阵。
//
// Inf / NaN 要单独处理。base 对 Inf 给 Inf、对 NaN 给 NaN；而缩放算法在
// max = Inf 时会算出 Inf/Inf = NaN，所以不能让它进循环。NaN 优先于 Inf，
// 与 base 一致。
template <typename Derived>
double frobenius_scaled(const Eigen::MatrixBase<Derived>& A) {
  const typename Derived::Scalar* d = A.derived().data();
  const Eigen::Index n = A.size();
  double mx = 0.0, ssq = 1.0;   // LAPACK 的初值：全零时 mx*sqrt(ssq) = 0
  bool inf = false, nan = false;
  for (Eigen::Index i = 0; i < n; ++i) {
    const double av = std::abs(d[i]);
    if (ISNAN(av)) { nan = true; continue; }
    if (av == R_PosInf) { inf = true; continue; }
    if (av == 0.0) continue;
    if (av > mx) { const double r = mx / av; ssq = 1.0 + ssq * r * r; mx = av; }
    else         { const double r = av / mx; ssq += r * r; }
  }
  if (nan) return R_NaN;
  if (inf) return R_PosInf;
  return mx * std::sqrt(ssq);
}

// type: 0 "F", 1 "O"/"1"（最大列绝对和）, 2 "I"（最大行绝对和）,
//       3 "M"（最大绝对值）, 4 "2"（谱范数）, 5 一般 p 范数（展平）
// [[Rcpp::export]]
SEXP reigen_norm(SEXP x, int type, double p) {
  Mat X = from_sexp(x, "x");
  if (X.cplx) {
    Eigen::MatrixXcd A = X.c;
    switch (type) {
      case 0: return Rcpp::wrap(frobenius_scaled(A));
      case 1: { double best = 0; for (Eigen::Index j = 0; j < A.cols(); ++j) best = std::max(best, A.col(j).cwiseAbs().sum()); return Rcpp::wrap(best); }
      case 2: { double best = 0; for (Eigen::Index i = 0; i < A.rows(); ++i) best = std::max(best, A.row(i).cwiseAbs().sum()); return Rcpp::wrap(best); }
      case 3: return Rcpp::wrap(A.size() ? A.cwiseAbs().maxCoeff() : 0.0);
      case 4: {
        if (A.size() == 0) return Rcpp::wrap(0.0);
        return Rcpp::wrap(Eigen::JacobiSVD<Eigen::MatrixXcd>(A).singularValues()(0));
      }
      case 5: {
        if (p <= 0) Rcpp::stop("p must be positive");
        double s = 0;
        for (Eigen::Index i = 0; i < A.size(); ++i) s += std::pow(std::abs(A.data()[i]), p);
        return Rcpp::wrap(std::pow(s, 1.0 / p));
      }
      default: Rcpp::stop("unknown norm type %d", type);
    }
  }
  const Eigen::MatrixXd& A = X.r;
  switch (type) {
    case 0: return Rcpp::wrap(frobenius_scaled(A));
    case 1: { double best = 0; for (Eigen::Index j = 0; j < A.cols(); ++j) best = std::max(best, A.col(j).cwiseAbs().sum()); return Rcpp::wrap(best); }
    case 2: { double best = 0; for (Eigen::Index i = 0; i < A.rows(); ++i) best = std::max(best, A.row(i).cwiseAbs().sum()); return Rcpp::wrap(best); }
    case 3: return Rcpp::wrap(A.size() ? A.cwiseAbs().maxCoeff() : 0.0);
    case 4: {
      if (A.size() == 0) return Rcpp::wrap(0.0);
      return Rcpp::wrap(Eigen::JacobiSVD<Eigen::MatrixXd>(A).singularValues()(0));
    }
    case 5: {
      if (p <= 0) Rcpp::stop("p must be positive");
      double s = 0;
      for (Eigen::Index i = 0; i < A.size(); ++i) s += std::pow(std::abs(A.data()[i]), p);
      return Rcpp::wrap(std::pow(s, 1.0 / p));
    }
    default: Rcpp::stop("unknown norm type %d", type);
  }
}

// ============================================================ 内积 / 叉积

// 与 R 的 sum(a * b) 一致：不取共轭。Eigen 的 .dot() 对复数会共轭第一个
// 参数，而且 .dot() 是向量方法，对矩阵会触发静态断言，所以这里手写。
// [[Rcpp::export]]
SEXP reigen_dot(SEXP a, SEXP b) {
  Mat A = from_sexp(a, "a");
  Mat B = from_sexp(b, "b");
  if (nelem(A) != nelem(B))
    Rcpp::stop("dot: lengths differ (%lld vs %lld)", (long long)nelem(A), (long long)nelem(B));
  Mat out;
  if (A.cplx || B.cplx) {
    Eigen::MatrixXcd X = as_cplx(A), Y = as_cplx(B);
    out.cplx = true;
    out.c.resize(1, 1);
    out.c(0, 0) = (X.array() * Y.array()).sum();
  } else {
    out.r.resize(1, 1);
    out.r(0, 0) = (A.r.array() * B.r.array()).sum();
  }
  return wrap_mat(std::move(out));
}

// [[Rcpp::export]]
SEXP reigen_cross(SEXP a, SEXP b) {
  Mat A = from_sexp(a, "a");
  Mat B = from_sexp(b, "b");
  if (A.cplx || B.cplx) Rcpp::stop("cross product is defined for real vectors only");
  // 二维叉积不是一个向量，而是一个标量（z 分量）：Eigen 5 的
  // Vector2d::cross 就是 double。这里照它的语义返回长度 1 的数值，
  // 而不是硬塞成 3 维向量的第三个分量——后者会让 length() 说谎。
  if (nelem(A) == 2 && nelem(B) == 2)
    return Rcpp::wrap(A.r.data()[0] * B.r.data()[1] - A.r.data()[1] * B.r.data()[0]);
  if (nelem(A) != 3 || nelem(B) != 3)
    Rcpp::stop("cross product requires two length-3 vectors, or two length-2 vectors");
  Eigen::Vector3d x, y;
  for (int i = 0; i < 3; ++i) { x(i) = A.r.data()[i]; y(i) = B.r.data()[i]; }
  Eigen::Vector3d z = x.cross(y);
  Eigen::MatrixXd out(3, 1);
  out.col(0) = z;
  return wrap_mat(make(out));
}

// ============================================================ 线性代数

// [[Rcpp::export]]
SEXP reigen_det(SEXP x) {
  Mat X = from_sexp(x, "x");
  if (rows(X) != cols(X)) Rcpp::stop("determinant requires a square matrix");
  if (X.cplx) return wrap_cplx(X.c.determinant());
  return Rcpp::wrap(X.r.determinant());
}

// log|det|，大矩阵上比 det() 安全：det 很容易溢出成 Inf/0。
// 注意不能写成 log(abs(determinant()))——那一步已经溢出了，再取对数救不回来。
// [[Rcpp::export]]
SEXP reigen_logdet(SEXP x) {
  Mat X = from_sexp(x, "x");
  if (rows(X) != cols(X)) Rcpp::stop("logdet requires a square matrix");
  double logabs = 0.0;
  cplx_t phase;
  if (X.cplx) det_parts(X.c, logabs, phase);
  else        det_parts(X.r, logabs, phase);
  return Rcpp::wrap(logabs);
}

// determinant() 需要 log|det| 和符号两样东西，一次算完。
// 返回 c(log|det|, sign_re, sign_im)：实数矩阵 sign_im 恒为 0，
// 复数矩阵的 sign 是单位复数（det = sign * exp(modulus) 仍然成立）。
// [[Rcpp::export]]
SEXP reigen_det_parts(SEXP x) {
  Mat X = from_sexp(x, "x");
  if (rows(X) != cols(X)) Rcpp::stop("determinant requires a square matrix");
  double logabs = 0.0;
  cplx_t phase(1.0, 0.0);
  if (X.cplx) det_parts(X.c, logabs, phase);
  else        det_parts(X.r, logabs, phase);
  return Rcpp::NumericVector::create(logabs, phase.real(), phase.imag());
}

// 秩和条件数是标量，不包成 eigenMat——否则 R 侧拿到的是一个外部指针，
// 用户还得再 as.vector(as.matrix(...)) 才能用，没有道理。
template <class S>
static double rank_t(const Mat& X, double tol) {
  typedef Eigen::Matrix<S, -1, -1> M;
  M A = to_eigen<S>(X);
  Eigen::FullPivLU<M> lu(A);
  if (tol > 0) lu.setThreshold(tol);
  return (double)lu.rank();
}

// 倒数条件数：sigma_min / sigma_max。方阵、长方形都成立，奇异矩阵给 0。
// 实现见 reigen.h 的 rcond_svd()（那里记了为什么不用 FullPivLU::rcond()）。
template <class S>
static double rcond_t(const Mat& X) {
  return rcond_svd(to_eigen<S>(X));
}

// [[Rcpp::export]]
double reigen_rank(SEXP x, double tol) {
  Mat X = from_sexp(x, "x");
  return X.cplx ? rank_t<cplx_t>(X, tol) : rank_t<double>(X, tol);
}

// [[Rcpp::export]]
double reigen_rcond(SEXP x) {
  Mat X = from_sexp(x, "x");
  return X.cplx ? rcond_t<cplx_t>(X) : rcond_t<double>(X);
}

// [[Rcpp::export]]
SEXP reigen_inverse(SEXP x) {
  Mat X = from_sexp(x, "x");
  if (rows(X) != cols(X)) Rcpp::stop("inverse requires a square matrix");
  return wrap_mat(pick_cplx(X.cplx, [&](auto tag) {
    typedef typename scalar_of<decltype(tag)>::type S;
    typedef Eigen::Matrix<S, -1, -1> M;
    M A = to_eigen<S>(X);
    return make(M(Eigen::FullPivLU<M>(A).inverse()));
  }));
}

// 方阵走 LU；非方阵走 ColPivHouseholderQR 的最小二乘解。
// 这条分叉让 reigen_solve(A, b) 在 A 是 n x p (n != p) 时也能给出
// 与 R 的 qr.solve() 一致的答案，而不是直接报「非方阵」。
// [[Rcpp::export]]
SEXP reigen_solve(SEXP a, SEXP b) {
  Mat A = from_sexp(a, "a");
  Mat B = from_sexp(b, "b");
  if (rows(A) != rows(B))
    Rcpp::stop("solve: %lld rows in a but %lld in b", (long long)rows(A), (long long)rows(B));
  return wrap_mat(pick_cplx(A.cplx || B.cplx, [&](auto tag) {
    typedef typename scalar_of<decltype(tag)>::type S;
    typedef Eigen::Matrix<S, -1, -1> M;
    M X = to_eigen<S>(A), Y = to_eigen<S>(B);
    if (X.rows() == X.cols())
      return make(M(X.partialPivLu().solve(Y)));
    return make(M(X.colPivHouseholderQr().solve(Y)));
  }));
}

// ============================================================ 近似 / 有限性

// [[Rcpp::export]]
bool reigen_approx(SEXP a, SEXP b, double tol) {
  Mat A = from_sexp(a, "a");
  Mat B = from_sexp(b, "b");
  if (rows(A) != rows(B) || cols(A) != cols(B)) return false;
  if (A.cplx || B.cplx) {
    Eigen::MatrixXcd X = as_cplx(A), Y = as_cplx(B);
    return (X - Y).norm() <= tol * std::min(X.norm(), Y.norm());
  }
  return A.r.isApprox(B.r, tol);
}

// which: 0 allFinite, 1 hasNaN, 2 hasInf
// [[Rcpp::export]]
bool reigen_finite(SEXP x, int which) {
  Mat X = from_sexp(x, "x");
  if (X.cplx) {
    const Eigen::MatrixXcd& A = X.c;
    switch (which) {
      case 0: return A.array().isFinite().all();
      case 1: return A.array().isNaN().any();
      default: {
        for (Eigen::Index i = 0; i < A.size(); ++i) {
          if (std::isinf(A.data()[i].real()) || std::isinf(A.data()[i].imag())) return true;
        }
        return false;
      }
    }
  }
  const Eigen::MatrixXd& A = X.r;
  switch (which) {
    case 0: return A.array().isFinite().all();
    case 1: return A.array().isNaN().any();
    default: return A.array().isInf().any();
  }
}
