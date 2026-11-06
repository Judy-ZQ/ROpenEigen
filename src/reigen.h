// reigen.h —— ROpenEigen 转换层（手写，不生成）
//
// 与 ROpenCV 的 ropecv.h 同构，区别只在承载类型：
//
//   1. 矩阵在 R 侧是 EXTPTRSXP，class = c("eigenMat", "eigenObject")，
//      指针指向下面的 Mat。释放靠 XPtr 的 finalizer（delete），不需要注册
//      R 的 finalizer，也不需要 .Call 里手工管理生命周期。
//   2. 实矩阵用 Eigen::MatrixXd，复矩阵用 Eigen::MatrixXcd。一个 Mat 只装一种，
//      用 cplx 标记。int / logical 输入一律升为 double——Eigen 里 MatrixXi 和
//      MatrixXd 是两套不互通的类型，若保留整型，每个运算都要写两遍，
//      而 R 用户对「1:4 变成 double」的容忍度远高于对「一半函数不支持整型」的容忍度。
//      这一点与 RcppEigen 的 as<MatrixXd> 一致。
//   3. R 矩阵与 Eigen 矩阵都是列主序，所以 R -> Mat 是一次 memcpy，
//      不是逐元素循环。复数因为 Rcomplex 与 std::complex 的布局只是「事实上」
//      相同（标准不保证），走循环，不走 memcpy。
//   4. R 的 GC 看不到 Eigen 的堆分配。若不自记账，R 会一直以为内存充裕，
//      finalizer 迟迟不跑。所以 wrap_mat 累加分配字节数，超过 128MB 主动
//      R_gc() 一次——这就是 ROpenCV 里那套压力阀的移植。

#pragma once

#define R_NO_REMAP

#include <Rcpp.h>
#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace reigen {

typedef std::complex<double> cplx_t;

// ============================================================== GC 压力阀

inline size_t& gc_bytes() {
  static size_t v = 0;
  return v;
}

inline void gc_tick(size_t bytes) {
  gc_bytes() += bytes;
  if (gc_bytes() >= (size_t)128 * 1024 * 1024) {
    gc_bytes() = 0;
    R_gc();
  }
}

// ============================================================== 矩阵载体

struct Mat {
  bool cplx;
  Eigen::MatrixXd  r;
  Eigen::MatrixXcd c;
  Mat() : cplx(false) {}
};

inline Eigen::Index rows(const Mat& m)  { return m.cplx ? m.c.rows() : m.r.rows(); }
inline Eigen::Index cols(const Mat& m)  { return m.cplx ? m.c.cols() : m.r.cols(); }
inline Eigen::Index nelem(const Mat& m) { return rows(m) * cols(m); }

inline Mat make(const Eigen::MatrixXd& x)  { Mat m; m.r = x; return m; }
inline Mat make(const Eigen::MatrixXcd& x) { Mat m; m.cplx = true; m.c = x; return m; }

inline Eigen::MatrixXd as_real(const Mat& m, const char* what) {
  if (m.cplx)
    Rcpp::stop("%s: complex matrix is not supported here; take Re() or Im() first", what);
  return m.r;
}

inline Eigen::MatrixXcd as_cplx(const Mat& m) {
  return m.cplx ? m.c : m.r.template cast<cplx_t>();
}

// 倒数条件数 = sigma_min / sigma_max（2-范数意义），方阵和长方形都成立，
// 奇异矩阵给 0。
//
// 为什么不用 FullPivLU::rcond()：它内部走 rcond_estimate_helper()，那里有
// eigen_assert(rows == cols)，而 R 编译包默认 -DNDEBUG 把断言编掉，长方形上
// 是越界读（随机段错误）。而且对恰好奇异的矩阵它给出的估计不可靠——实测
// 一个 2x2 奇异矩阵返回 1，而正确值是 0。SVD 慢一点，但定义清楚、不会崩。
template <class M>
inline double rcond_svd(const M& A) {
  if (A.rows() == 0 || A.cols() == 0) return 1.0;
  Eigen::JacobiSVD<M> svd(A);
  const Eigen::Index k = svd.singularValues().rows();
  if (k == 0) return 1.0;
  const double smax = (double)svd.singularValues()(0);
  const double smin = (double)svd.singularValues()(k - 1);
  return smax > 0.0 ? smin / smax : 0.0;
}

// det = sign(P) * prod(U_ii)。这里不真的把乘积算出来，而是逐项取对数再相加，
// 顺便把符号从「U 对角线的相位 × 置换的奇偶」凑出来。
//
// 为什么不能直接 determinant()：500x500 的随机矩阵行列式在 1e±300 量级，
// 中间乘积立刻溢出成 Inf 或下溢成 0，于是 log|det| 变成 Inf 或 -Inf，
// 而正确答案是个正常的有限数。base R 的 determinant() 同样是 LU 累加对数，
// 所以这里必须照做，否则同一份数据在 base R 和本包里给出两个答案。
//
// 置换的符号：把 P 的 indices 看成置换，符号 = (-1)^(n - 圈数)。
// Eigen 把 m_det_p 藏成私有成员，只能这样自己算。
template <class M>
inline void det_parts(const M& A, double& logabs, cplx_t& phase) {
  typedef typename M::Scalar S;
  logabs = 0.0;
  phase = cplx_t(1.0, 0.0);
  if (A.rows() == 0) return;

  Eigen::PartialPivLU<M> lu(A);
  const M& LU = lu.matrixLU();
  for (Eigen::Index i = 0; i < LU.rows(); ++i) {
    const S d = LU(i, i);
    const double a = (double)std::abs(d);
    if (a == 0.0) {
      // 奇异：log|det| = -Inf，符号无意义。base R 在这种情况下仍给 ±1
      // （它取的是置换的符号），这里同样返回 1，免得 det() 变成 0 * 0。
      logabs = -std::numeric_limits<double>::infinity();
      phase = cplx_t(1.0, 0.0);
      return;
    }
    logabs += std::log(a);
    phase *= cplx_t(d) / a;     // 复数取单位相位，实数取 ±1
  }

  const auto& idx = lu.permutationP().indices();
  const int n = (int)idx.size();
  std::vector<char> seen((size_t)n, 0);
  int cycles = 0;
  for (int i = 0; i < n; ++i) {
    if (seen[(size_t)i]) continue;
    ++cycles;
    int j = i;
    while (!seen[(size_t)j]) { seen[(size_t)j] = 1; j = (int)idx[j]; }
  }
  if (((n - cycles) & 1) != 0) phase = -phase;
}

// 按 cplx 标记把匿名函数分派到 double / complex 两个实现上。
// 用法： pick_cplx(m.cplx, [&](auto tag) { using S = typename scalar_of<decltype(tag)>::type; ... })
template <class F>
inline Mat pick_cplx(bool c, F&& f) {
  return c ? f(std::true_type{}) : f(std::false_type{});
}

template <class Tag>          struct scalar_of   { typedef double  type; };
template <>                   struct scalar_of<std::true_type>  { typedef cplx_t type; };

template <class S>            struct mat_of      { typedef Eigen::Matrix<S, -1, -1> type; };

// Mat -> Eigen::Matrix<S>
template <class S> inline Eigen::Matrix<S, -1, -1> to_eigen(const Mat& m);
template <> inline Eigen::MatrixXd   to_eigen<double>(const Mat& m) { return as_real(m, "to_eigen"); }
template <> inline Eigen::MatrixXcd to_eigen<cplx_t>(const Mat& m) { return as_cplx(m); }

// ============================================================== 对象包装
//
// 这几个函数拿到的 x 往往是刚 wrap 出来、R 侧还没有任何引用的对象。
// 函数体里只要发生一次 R 对象分配（Rcpp::CharacterVector::create、
// IntegerVector::create、Rf_install……），GC 就可能跑起来；那一刻 x 没有
// 任何保护，会被回收——XPtr 的 finalizer 连 C++ 对象一起 delete——随后
// Rf_setAttrib 就是往已释放的内存里写。这种堆损坏不会当场报错，而是在
// 后面某次随机分配时炸成段错误，极难定位（本包就为此吃过一次）。
// 所以：只要参数是"裸 SEXP 且后面还要分配"，一律先 PROTECT。

inline SEXP set_class(SEXP x, Rcpp::CharacterVector cl) {
  PROTECT(x);
  Rf_setAttrib(x, R_ClassSymbol, cl);
  UNPROTECT(1);
  return x;
}

inline SEXP set_class(SEXP x, const char* cls) {
  PROTECT(x);
  SEXP r = set_class(x, Rcpp::CharacterVector::create(cls, "eigenObject"));
  UNPROTECT(1);
  return r;
}

// eigenVec 要排在 eigenMat 前面，S3 派发才会先找 eigenVec 的方法、
// 找不到再顺着 class 向量落到 eigenMat。少了中间这档，eigenVec 就只剩
// eigenObject 兜底，dim()/print() 全部失灵。
inline SEXP set_class_vec(SEXP x) {
  PROTECT(x);
  SEXP r = set_class(x, Rcpp::CharacterVector::create("eigenVec", "eigenMat", "eigenObject"));
  UNPROTECT(1);
  return r;
}

template <class T>
inline SEXP wrap_obj(T* p, const char* cls, size_t bytes = 0) {
  Rcpp::XPtr<T> xp(p);
  SEXP x = xp;
  PROTECT(x);
  set_class(x, cls);
  if (bytes) gc_tick(bytes);   // 可能调 R_gc()，x 必须已被保护
  UNPROTECT(1);
  return x;
}

// Rcpp 没有 std::complex<double> 的 wrap，这里手工造长度 1 的 complex 向量。
inline SEXP wrap_cplx(cplx_t v) {
  Rcpp::ComplexVector out(1);
  Rcomplex* p = reinterpret_cast<Rcomplex*>(out.begin());
  p[0].r = v.real();
  p[0].i = v.imag();
  return out;
}

// 把 1x1 的 Mat 转成 R 标量。
inline SEXP scalar_of_mat(const Mat& m) {
  if (m.cplx) return wrap_cplx(m.c.data()[0]);
  return Rcpp::wrap(m.r.data()[0]);
}

inline SEXP wrap_mat(Mat m) {
  size_t bytes = (size_t)nelem(m) * (m.cplx ? 16u : 8u);
  return wrap_obj(new Mat(std::move(m)), "eigenMat", bytes);
}

// 分解对象自己没有形状可推，构造时把原矩阵的尺寸挂在自定义属性
// "eigenDim" 上。这样 R 侧 dim()、以及 $solve() 在缺省右端项时
// 造单位阵，都不用再回 C++ 问一次。
// 用自定义名字而不是 "dim"：XPtr 一旦有了 dim 属性，is.array() 会
// 当成数组，很多泛型会按数组对待它，那是假象。
//
// 符号缓存成函数内 static：Rf_install 首次调用会分配，之后直接返回。
// 不放在全局，是因为全局静态初始化时 R 还没起来。
inline SEXP dim_symbol() {
  static SEXP s = Rf_install("eigenDim");
  return s;
}

inline SEXP set_dim_attr(SEXP x, int nr, int nc) {
  PROTECT(x);
  Rf_setAttrib(x, dim_symbol(), Rcpp::IntegerVector::create(nr, nc));
  UNPROTECT(1);
  return x;
}

inline SEXP get_dim_attr(SEXP x) {
  PROTECT(x);
  SEXP r = Rf_getAttrib(x, dim_symbol());
  UNPROTECT(1);
  return r;
}

inline bool is_eigen(SEXP x) {
  return TYPEOF(x) == EXTPTRSXP && Rf_inherits(x, "eigenMat");
}

inline Mat* mat_ptr(SEXP x, const char* what = "x") {
  if (!is_eigen(x))
    Rcpp::stop("%s must be an eigenMat object, got %s", what, Rf_type2char(TYPEOF(x)));
  Mat* p = static_cast<Mat*>(R_ExternalPtrAddr(x));
  if (!p) Rcpp::stop("%s: the eigenMat pointer is null (object already released)", what);
  return p;
}

template <class T>
inline T* obj_ptr(SEXP x, const char* cls, const char* what) {
  if (TYPEOF(x) != EXTPTRSXP || !Rf_inherits(x, cls))
    Rcpp::stop("%s must be a %s object, got %s", what, cls, Rf_type2char(TYPEOF(x)));
  T* p = static_cast<T*>(R_ExternalPtrAddr(x));
  if (!p) Rcpp::stop("%s: the %s pointer is null (object already released)", what, cls);
  return p;
}

// ============================================================== R -> Mat

// 接受：eigenMat、R 矩阵、无 dim 的数值/逻辑/复数向量、单个数字。
// 无 dim 的一律当列向量（n x 1）——与 Eigen 的 VectorXd 语义对齐。
inline Mat from_sexp(SEXP x, const char* what = "x") {
  if (is_eigen(x)) return *mat_ptr(x, what);

  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  Eigen::Index nr, nc;
  if (dim != R_NilValue) {
    if (Rf_length(dim) != 2)
      Rcpp::stop("%s: only 2-d arrays are accepted, dim has length %d",
                 what, (int)Rf_length(dim));
    nr = (Eigen::Index)INTEGER(dim)[0];
    nc = (Eigen::Index)INTEGER(dim)[1];
  } else {
    nr = (Eigen::Index)Rf_xlength(x);
    nc = 1;
  }

  const Eigen::Index n = nr * nc;
  Mat m;

  if (Rf_isComplex(x)) {
    m.cplx = true;
    m.c.resize(nr, nc);
    const Rcomplex* p = COMPLEX(x);
    for (Eigen::Index i = 0; i < n; ++i)
      m.c.data()[i] = cplx_t(p[i].r, p[i].i);
    return m;
  }

  if (!(Rf_isReal(x) || Rf_isInteger(x) || Rf_isLogical(x)))
    Rcpp::stop("%s: expected a numeric matrix or vector, got %s",
               what, Rf_type2char(TYPEOF(x)));

  m.cplx = false;
  m.r.resize(nr, nc);
  if (Rf_isReal(x)) {
    if (n > 0) std::memcpy(m.r.data(), REAL(x), (size_t)n * sizeof(double));
  } else if (Rf_isInteger(x)) {
    const int* p = INTEGER(x);
    for (Eigen::Index i = 0; i < n; ++i) m.r.data()[i] = (double)p[i];
  } else {
    const int* p = LOGICAL(x);
    for (Eigen::Index i = 0; i < n; ++i)
      m.r.data()[i] = (p[i] == NA_LOGICAL) ? NA_REAL : (double)p[i];
  }
  return m;
}

// ============================================================== Mat -> R

// 分解专用：拒绝含 NaN/Inf 的输入。
// base R 的 svd()/qr()/solve()/eigen() 遇到缺失值都是直接报错；如果这里放过去，
// Eigen 会安静地产出 NaN，用户拿到一坨垃圾却没有任何提示。
// 只给分解用，不给逐元素运算/归约用——那些场景 NA 应该自然传播。
inline Mat from_sexp_finite(SEXP a, const char* what) {
  Mat m = from_sexp(a, what);
  const bool ok = m.cplx ? m.c.allFinite() : m.r.allFinite();
  if (!ok)
    Rcpp::stop("%s: matrix contains NaN or Inf; base R's decompositions reject "
               "such input too -- drop or impute the missing values first", what);
  return m;
}

inline SEXP to_sexp(const Mat& m) {
  const Eigen::Index nr = rows(m), nc = cols(m);
  if (m.cplx) {
    Rcpp::ComplexMatrix out((int)nr, (int)nc);
    Rcomplex* p = reinterpret_cast<Rcomplex*>(out.begin());
    const Eigen::Index n = nr * nc;
    for (Eigen::Index i = 0; i < n; ++i) {
      p[i].r = m.c.data()[i].real();
      p[i].i = m.c.data()[i].imag();
    }
    return out;
  }
  Rcpp::NumericMatrix out((int)nr, (int)nc);
  const Eigen::Index n = nr * nc;
  if (n > 0) std::memcpy(out.begin(), m.r.data(), (size_t)n * sizeof(double));
  return out;
}

// 1x1 取出成 R 标量（复数返回 complex 长度 1）。
inline SEXP to_scalar(SEXP x) {
  const Mat& m = *mat_ptr(x);
  if (nelem(m) != 1)
    Rcpp::stop("expected a 1x1 eigenMat, got %lldx%lld",
               (long long)rows(m), (long long)cols(m));
  return scalar_of_mat(m);
}

// 索引换算：R 是 1-based，C++ 是 0-based；0 表示「不指定」，由调用方展开。
inline Eigen::Index idx1(int i, Eigen::Index n, const char* what) {
  if (i < 1 || (Eigen::Index)i > n)
    Rcpp::stop("%s index %d out of range [1, %lld]", what, i, (long long)n);
  return (Eigen::Index)(i - 1);
}

}  // namespace reigen
