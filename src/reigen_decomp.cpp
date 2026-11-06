// reigen_decomp.cpp —— 矩阵分解
//
// 设计：每个分解在 R 侧是一个 EXTPTRSXP，class = c("eigenLLT", "eigenObject")
// 之类。分解算一次，结果留在 C++ 侧；后续 $solve() / $matrixL() 反复用，
// 不再重算——这是 Eigen 相对「每次调用都重新分解」的全部意义。
//
// 每个分解都写成 template<S> 的实现，double 和 complex<double> 各存一份成员，
// 用 cplx 标记走哪边。模板只写一遍，实例化两遍；这比给复数单独开一套类
// 省一半代码，也比「复数一律转成实部再算」诚实。
//
// 与 ROpenCV 的对应关系：cv 那边的分解是 cv::solve 之类的一次性函数，
// 没有可复用的对象。Eigen 的分解对象本身有状态，所以这里多了一层
// XPtr 包住 solver，这正是「结合 OpenCV 经验」时需要加出来的部分。

#include "reigen.h"

using namespace reigen;

// 取 Impl 里的矩阵类型 / 标量类型，省掉一长串 decay_t<decltype(...)>。
#define REIGEN_M(d) typename std::decay_t<decltype(d)>::M
#define REIGEN_S(d) typename std::decay_t<decltype(d)>::Scalar

// 两个存储（实 / 复）各一份的分解对象。
#define REIGEN_DUAL_OBJ(ObjT, ImplT) \
  struct ObjT { bool cplx = false; ImplT<double> r; ImplT<cplx_t> c; }

// 按 cplx 把 lambda 分派到 r 或 c 上。
template <class Obj, class F>
static SEXP with_obj(SEXP x, const char* cls, F&& f) {
  Obj* p = obj_ptr<Obj>(x, cls, cls);
  return p->cplx ? f(p->c) : f(p->r);
}

// 方阵校验：大部分分解要求方阵，非方阵只有 QR / SVD 有意义。
static void need_square(const Mat& A, const char* what) {
  if (rows(A) != cols(A))
    Rcpp::stop("%s requires a square matrix, got %lldx%lld", what,
               (long long)rows(A), (long long)cols(A));
}

// ============================================================== Cholesky LLT
//
// A = L L^H，要求 A 对称正定。只读下三角，所以 A 的上三角填什么值都不影响
// 结果——这一点和 Eigen 的约定一致，也意味着「传进来的矩阵不对称」不会报错，
// 只会静默地当成对称矩阵的下三角。文档里必须写清楚。

template <class S>
struct LLTImpl {
  typedef S Scalar;
  typedef Eigen::Matrix<S, -1, -1> M;
  Eigen::LLT<M> d;

  void compute(const M& a) { d.compute(a); }
  Mat solve(const Mat& b) const { return make(M(d.solve(to_eigen<S>(b)))); }
  Mat matrixL() const { return make(M(d.matrixL())); }
  Mat matrixU() const { return make(M(d.matrixU())); }
  Mat reconstructed() const { return make(M(d.reconstructedMatrix())); }
  // 和 eigen_rcond() 同一口径（sigma_min / sigma_max）。Eigen 的 LLT::rcond()
  // 用的是另一套估计（对 diag(L) 做比较，量纲上其实是 1/cond(L)^2 那一类近似），
  // 同一个矩阵走 $rcond() 和 eigen_rcond() 会给出两个明显不同的数。
  // 代价是多一次 SVD；rcond 本来就是诊断用的，不值得为它牺牲一致性。
  double rcond() const { return rcond_svd(M(d.reconstructedMatrix())); }
  bool ok() const { return d.info() == Eigen::Success; }

  // LLT 没有 logdet()：det(A) = prod(diag(L))^2，取对数累加即可。
  // 用 abs 是因为复数下 diag(L) 是复数，而正定矩阵的行列式必为正实数。
  double logdet() const {
    M L = d.matrixL();
    double s = 0;
    for (Eigen::Index i = 0; i < L.rows(); ++i) s += std::log(std::abs(L(i, i)));
    return 2.0 * s;
  }
};

REIGEN_DUAL_OBJ(LLTObj, LLTImpl);

// [[Rcpp::export]]
SEXP reigen_llt_new(SEXP a) {
  Mat A = from_sexp_finite(a, "a");
  need_square(A, "LLT");
  LLTObj* p = new LLTObj();
  p->cplx = A.cplx;
  if (A.cplx) p->c.compute(A.c); else p->r.compute(A.r);

  // 非正定矩阵上 Eigen 的 LLT 不抛异常，只在 info() 里留一个 NumericalIssue，
  // 而 solve()/matrixL()/reconstructed() 照样有返回值——返回的是**另一个矩阵**
  // 的分解。实测 2x2 不定矩阵 [[1,2],[2,1]] 重构出来是 [[1,2],[2,5]]。
  // 这种"能算、但算的不是你的矩阵"必须当场拦下；base R 的 chol() 也是直接报错。
  // 注意这里 p 还没交给 XPtr，stop 之前必须自己 delete，否则泄漏。
  const bool ok = p->cplx ? p->c.ok() : p->r.ok();
  if (!ok) {
    delete p;
    Rcpp::stop("matrix is not positive definite; Cholesky failed. "
               "For a symmetric matrix that is merely indefinite, use "
               "eigen_ldlt() instead -- LDL^T does not require definiteness.");
  }

  return set_dim_attr(wrap_obj(p, "eigenLLT"), (int)rows(A), (int)cols(A));
}

// [[Rcpp::export]]
SEXP reigen_llt_solve(SEXP o, SEXP b) {
  Mat B = from_sexp(b, "b");
  return with_obj<LLTObj>(o, "eigenLLT", [&](const auto& d) {
    return wrap_mat(d.solve(B));
  });
}

// [[Rcpp::export]]
SEXP reigen_llt_matrixL(SEXP o) {
  return with_obj<LLTObj>(o, "eigenLLT", [&](const auto& d) { return wrap_mat(d.matrixL()); });
}

// [[Rcpp::export]]
SEXP reigen_llt_matrixU(SEXP o) {
  return with_obj<LLTObj>(o, "eigenLLT", [&](const auto& d) { return wrap_mat(d.matrixU()); });
}

// [[Rcpp::export]]
SEXP reigen_llt_reconstructed(SEXP o) {
  return with_obj<LLTObj>(o, "eigenLLT", [&](const auto& d) { return wrap_mat(d.reconstructed()); });
}

// [[Rcpp::export]]
double reigen_llt_logdet(SEXP o) {
  LLTObj* p = obj_ptr<LLTObj>(o, "eigenLLT", "eigenLLT");
  return p->cplx ? p->c.logdet() : p->r.logdet();
}

// [[Rcpp::export]]
double reigen_llt_rcond(SEXP o) {
  LLTObj* p = obj_ptr<LLTObj>(o, "eigenLLT", "eigenLLT");
  return p->cplx ? p->c.rcond() : p->r.rcond();
}

// [[Rcpp::export]]
bool reigen_llt_ok(SEXP o) {
  LLTObj* p = obj_ptr<LLTObj>(o, "eigenLLT", "eigenLLT");
  return p->cplx ? p->c.ok() : p->r.ok();
}

// ============================================================== LDLT
//
// A = P^T L D L^T P。比 LLT 宽松：不要求正定，只要对称。D 的对角线带符号，
// 所以 logdet 要按符号累加（这里返回 log|det|，与 reigen_logdet 一致）。

template <class S>
struct LDLTImpl {
  typedef S Scalar;
  typedef Eigen::Matrix<S, -1, -1> M;
  Eigen::LDLT<M> d;

  void compute(const M& a) { d.compute(a); }
  Mat solve(const Mat& b) const { return make(M(d.solve(to_eigen<S>(b)))); }
  Mat matrixL() const { return make(M(d.matrixL())); }
  Mat matrixU() const { return make(M(d.matrixU())); }
  Mat vectorD() const { return make(M(d.vectorD())); }
  Mat reconstructed() const { return make(M(d.reconstructedMatrix())); }
  // 同 LLTImpl::rcond：统一走 sigma_min / sigma_max，和 eigen_rcond() 对齐。
  double rcond() const { return rcond_svd(M(d.reconstructedMatrix())); }
  bool ok() const { return d.info() == Eigen::Success; }

  double logdet() const {
    M D = d.vectorD();
    double s = 0;
    for (Eigen::Index i = 0; i < D.rows(); ++i) s += std::log(std::abs(D(i, 0)));
    return s;
  }
};

REIGEN_DUAL_OBJ(LDLTObj, LDLTImpl);

// [[Rcpp::export]]
SEXP reigen_ldlt_new(SEXP a) {
  Mat A = from_sexp_finite(a, "a");
  need_square(A, "LDLT");
  LDLTObj* p = new LDLTObj();
  p->cplx = A.cplx;
  if (A.cplx) p->c.compute(A.c); else p->r.compute(A.r);
  return set_dim_attr(wrap_obj(p, "eigenLDLT"), (int)rows(A), (int)cols(A));
}

// [[Rcpp::export]]
SEXP reigen_ldlt_solve(SEXP o, SEXP b) {
  Mat B = from_sexp(b, "b");
  return with_obj<LDLTObj>(o, "eigenLDLT", [&](const auto& d) { return wrap_mat(d.solve(B)); });
}

// [[Rcpp::export]]
SEXP reigen_ldlt_matrixL(SEXP o) {
  return with_obj<LDLTObj>(o, "eigenLDLT", [&](const auto& d) { return wrap_mat(d.matrixL()); });
}

// [[Rcpp::export]]
SEXP reigen_ldlt_matrixU(SEXP o) {
  return with_obj<LDLTObj>(o, "eigenLDLT", [&](const auto& d) { return wrap_mat(d.matrixU()); });
}

// [[Rcpp::export]]
SEXP reigen_ldlt_vectorD(SEXP o) {
  return with_obj<LDLTObj>(o, "eigenLDLT", [&](const auto& d) { return wrap_mat(d.vectorD()); });
}

// [[Rcpp::export]]
SEXP reigen_ldlt_reconstructed(SEXP o) {
  return with_obj<LDLTObj>(o, "eigenLDLT", [&](const auto& d) { return wrap_mat(d.reconstructed()); });
}

// [[Rcpp::export]]
double reigen_ldlt_logdet(SEXP o) {
  LDLTObj* p = obj_ptr<LDLTObj>(o, "eigenLDLT", "eigenLDLT");
  return p->cplx ? p->c.logdet() : p->r.logdet();
}

// [[Rcpp::export]]
double reigen_ldlt_rcond(SEXP o) {
  LDLTObj* p = obj_ptr<LDLTObj>(o, "eigenLDLT", "eigenLDLT");
  return p->cplx ? p->c.rcond() : p->r.rcond();
}

// [[Rcpp::export]]
bool reigen_ldlt_ok(SEXP o) {
  LDLTObj* p = obj_ptr<LDLTObj>(o, "eigenLDLT", "eigenLDLT");
  return p->cplx ? p->c.ok() : p->r.ok();
}

// ============================================================== PartialPivLU

template <class S>
struct LUImpl {
  typedef S Scalar;
  typedef Eigen::Matrix<S, -1, -1> M;
  Eigen::PartialPivLU<M> d;

  void compute(const M& a) { d.compute(a); }
  Mat solve(const Mat& b) const { return make(M(d.solve(to_eigen<S>(b)))); }
  Mat inverse() const { return make(M(d.inverse())); }
  Mat reconstructed() const { return make(M(d.reconstructedMatrix())); }

  // PartialPivLU / FullPivLU 只暴露 matrixLU()，不单独给 L 和 U。
  // L 是单位下三角（对角线隐含 1，赋值时不会被写入），所以先铺单位阵；
  // U 是上三角，先铺零阵。
  Mat matrixL() const {
    M L = M::Identity(d.matrixLU().rows(), d.matrixLU().rows());
    L = d.matrixLU().template triangularView<Eigen::UnitLower>();
    return make(L);
  }

  Mat matrixU() const {
    M U = M::Zero(d.matrixLU().rows(), d.matrixLU().cols());
    U = d.matrixLU().template triangularView<Eigen::Upper>();
    return make(U);
  }

  Mat permutationP() const {
    M P = M::Identity(d.matrixLU().rows(), d.matrixLU().rows());
    P = d.permutationP();
    return make(P);
  }

  SEXP det() const {
    if constexpr (std::is_same<S, double>::value) return Rcpp::wrap(d.determinant());
    else return wrap_cplx(d.determinant());
  }
};

REIGEN_DUAL_OBJ(LUObj, LUImpl);

// [[Rcpp::export]]
SEXP reigen_lu_new(SEXP a) {
  Mat A = from_sexp_finite(a, "a");
  need_square(A, "LU");
  LUObj* p = new LUObj();
  p->cplx = A.cplx;
  if (A.cplx) p->c.compute(A.c); else p->r.compute(A.r);
  return set_dim_attr(wrap_obj(p, "eigenLU"), (int)rows(A), (int)cols(A));
}

// [[Rcpp::export]]
SEXP reigen_lu_solve(SEXP o, SEXP b) {
  Mat B = from_sexp(b, "b");
  return with_obj<LUObj>(o, "eigenLU", [&](const auto& d) { return wrap_mat(d.solve(B)); });
}

// [[Rcpp::export]]
SEXP reigen_lu_matrixL(SEXP o) {
  return with_obj<LUObj>(o, "eigenLU", [&](const auto& d) { return wrap_mat(d.matrixL()); });
}

// [[Rcpp::export]]
SEXP reigen_lu_matrixU(SEXP o) {
  return with_obj<LUObj>(o, "eigenLU", [&](const auto& d) { return wrap_mat(d.matrixU()); });
}

// [[Rcpp::export]]
SEXP reigen_lu_permutationP(SEXP o) {
  return with_obj<LUObj>(o, "eigenLU", [&](const auto& d) { return wrap_mat(d.permutationP()); });
}

// [[Rcpp::export]]
SEXP reigen_lu_inverse(SEXP o) {
  return with_obj<LUObj>(o, "eigenLU", [&](const auto& d) { return wrap_mat(d.inverse()); });
}

// [[Rcpp::export]]
SEXP reigen_lu_reconstructed(SEXP o) {
  return with_obj<LUObj>(o, "eigenLU", [&](const auto& d) { return wrap_mat(d.reconstructed()); });
}

// [[Rcpp::export]]
SEXP reigen_lu_det(SEXP o) {
  return with_obj<LUObj>(o, "eigenLU", [&](const auto& d) { return d.det(); });
}

// ============================================================== FullPivLU
//
// 比 PartialPivLU 慢，但能给出数值秩、条件数、是否可逆——判断奇异与否
// 应该用它，而不是看 PartialPivLU 的解是否出现 NaN。

template <class S>
struct FullLUImpl {
  typedef S Scalar;
  typedef Eigen::Matrix<S, -1, -1> M;
  Eigen::FullPivLU<M> d;

  void compute(const M& a) { d.compute(a); }
  Mat solve(const Mat& b) const { return make(M(d.solve(to_eigen<S>(b)))); }
  Mat inverse() const { return make(M(d.inverse())); }
  Mat reconstructed() const { return make(M(d.reconstructedMatrix())); }

  // PartialPivLU / FullPivLU 只暴露 matrixLU()，不单独给 L 和 U。
  // L 是单位下三角（对角线隐含 1，赋值时不会被写入），所以先铺单位阵；
  // U 是上三角，先铺零阵。
  Mat matrixL() const {
    M L = M::Identity(d.matrixLU().rows(), d.matrixLU().rows());
    L = d.matrixLU().template triangularView<Eigen::UnitLower>();
    return make(L);
  }

  Mat matrixU() const {
    M U = M::Zero(d.matrixLU().rows(), d.matrixLU().cols());
    U = d.matrixLU().template triangularView<Eigen::Upper>();
    return make(U);
  }
  double rank() const { return (double)d.rank(); }
  // 和 eigen_rcond() 保持同一口径（sigma_min / sigma_max），否则同一个矩阵
  // 走 $rcond() 和 eigen_rcond() 会给出两个数。代价是要重构一次 A 再做 SVD。
  double rcond() const { return rcond_svd(M(d.reconstructedMatrix())); }
  bool invertible() const { return d.isInvertible(); }
  Eigen::Index rows() const { return d.matrixLU().rows(); }
  Eigen::Index cols() const { return d.matrixLU().cols(); }

  Mat permutationP() const {
    M P = M::Identity(d.matrixLU().rows(), d.matrixLU().rows());
    P = d.permutationP();
    return make(P);
  }
  Mat permutationQ() const {
    M Q = M::Identity(d.matrixLU().cols(), d.matrixLU().cols());
    Q = d.permutationQ();
    return make(Q);
  }

  SEXP det() const {
    if constexpr (std::is_same<S, double>::value) return Rcpp::wrap(d.determinant());
    else return wrap_cplx(d.determinant());
  }
};

REIGEN_DUAL_OBJ(FullLUObj, FullLUImpl);

// FullPivLU 是少数能吃非方阵的分解（rank/solve/matrixL/matrixU/permutation
// 都合法），但 determinant()、inverse()、rcond() 在 Eigen 内部只对方阵有定义：
// 它们都带 eigen_assert(rows == cols)，而 R 编译包时默认 -DNDEBUG 会把
// eigen_assert 编掉，非方阵上直接读到越界内存 —— 症状是随机位置的段错误。
// 所以这三个方法必须自己先拦。
template <class Impl>
static void need_square_fulllu(const Impl& d, const char* what) {
  if (d.rows() != d.cols())
    Rcpp::stop("%s() requires a square matrix, got %lldx%lld; "
               "use eigen_rcond(x), eigen_svd(x)$singular_values() or eigen_rank(x) instead",
               what, (long long)d.rows(), (long long)d.cols());
}

// [[Rcpp::export]]
SEXP reigen_fulllu_new(SEXP a) {
  Mat A = from_sexp_finite(a, "a");
  // FullPivLU 是少数支持非方阵的分解：矩阵可以有 m != n，
  // rank() / rcond() / solve() 都按最小二乘或秩揭示的意义工作。
  FullLUObj* p = new FullLUObj();
  p->cplx = A.cplx;
  if (A.cplx) p->c.compute(A.c); else p->r.compute(A.r);
  return set_dim_attr(wrap_obj(p, "eigenFullPivLU"), (int)rows(A), (int)cols(A));
}

// [[Rcpp::export]]
SEXP reigen_fulllu_solve(SEXP o, SEXP b) {
  Mat B = from_sexp(b, "b");
  return with_obj<FullLUObj>(o, "eigenFullPivLU", [&](const auto& d) { return wrap_mat(d.solve(B)); });
}

// [[Rcpp::export]]
SEXP reigen_fulllu_matrixL(SEXP o) {
  return with_obj<FullLUObj>(o, "eigenFullPivLU", [&](const auto& d) { return wrap_mat(d.matrixL()); });
}

// [[Rcpp::export]]
SEXP reigen_fulllu_matrixU(SEXP o) {
  return with_obj<FullLUObj>(o, "eigenFullPivLU", [&](const auto& d) { return wrap_mat(d.matrixU()); });
}

// [[Rcpp::export]]
SEXP reigen_fulllu_permutationP(SEXP o) {
  return with_obj<FullLUObj>(o, "eigenFullPivLU", [&](const auto& d) { return wrap_mat(d.permutationP()); });
}

// [[Rcpp::export]]
SEXP reigen_fulllu_permutationQ(SEXP o) {
  return with_obj<FullLUObj>(o, "eigenFullPivLU", [&](const auto& d) { return wrap_mat(d.permutationQ()); });
}

// [[Rcpp::export]]
SEXP reigen_fulllu_inverse(SEXP o) {
  FullLUObj* p = obj_ptr<FullLUObj>(o, "eigenFullPivLU", "eigenFullPivLU");
  if (p->cplx) { need_square_fulllu(p->c, "inverse"); return wrap_mat(p->c.inverse()); }
  need_square_fulllu(p->r, "inverse");
  return wrap_mat(p->r.inverse());
}

// [[Rcpp::export]]
SEXP reigen_fulllu_reconstructed(SEXP o) {
  return with_obj<FullLUObj>(o, "eigenFullPivLU", [&](const auto& d) { return wrap_mat(d.reconstructed()); });
}

// [[Rcpp::export]]
SEXP reigen_fulllu_det(SEXP o) {
  FullLUObj* p = obj_ptr<FullLUObj>(o, "eigenFullPivLU", "eigenFullPivLU");
  if (p->cplx) { need_square_fulllu(p->c, "det"); return p->c.det(); }
  need_square_fulllu(p->r, "det");
  return p->r.det();
}

// [[Rcpp::export]]
double reigen_fulllu_rank(SEXP o) {
  FullLUObj* p = obj_ptr<FullLUObj>(o, "eigenFullPivLU", "eigenFullPivLU");
  return p->cplx ? p->c.rank() : p->r.rank();
}

// [[Rcpp::export]]
double reigen_fulllu_rcond(SEXP o) {
  FullLUObj* p = obj_ptr<FullLUObj>(o, "eigenFullPivLU", "eigenFullPivLU");
  // 走 SVD（见 FullLUImpl::rcond），长方形也成立，不再需要方阵校验。
  return p->cplx ? p->c.rcond() : p->r.rcond();
}

// [[Rcpp::export]]
bool reigen_fulllu_invertible(SEXP o) {
  FullLUObj* p = obj_ptr<FullLUObj>(o, "eigenFullPivLU", "eigenFullPivLU");
  return p->cplx ? p->c.invertible() : p->r.invertible();
}

// ============================================================== HouseholderQR

template <class S>
struct QRImpl {
  typedef S Scalar;
  typedef Eigen::Matrix<S, -1, -1> M;
  Eigen::HouseholderQR<M> d;

  void compute(const M& a) { d.compute(a); }
  Mat solve(const Mat& b) const { return make(M(d.solve(to_eigen<S>(b)))); }

  Mat matrixQ() const {
    M I = M::Identity(d.rows(), d.rows());
    return make(M(d.householderQ() * I));
  }

  // matrixQR() 与 A 同形；只取上三角，下三角显式清零，
  // 免得把 Householder 向量当成 R 的一部分交出去。
  Mat matrixR() const {
    M R = M::Zero(d.rows(), d.cols());
    R = d.matrixQR().template triangularView<Eigen::Upper>();
    return make(R);
  }

  Mat reconstructed() const {
    M Q = M::Identity(d.rows(), d.rows());
    Q = d.householderQ() * Q;
    M R = M::Zero(d.rows(), d.cols());
    R = d.matrixQR().template triangularView<Eigen::Upper>();
    return make(M(Q * R));
  }

  double log_abs_det() const { return (double)d.logAbsDeterminant(); }
};

REIGEN_DUAL_OBJ(QRObj, QRImpl);

// [[Rcpp::export]]
SEXP reigen_qr_new(SEXP a) {
  Mat A = from_sexp_finite(a, "a");
  QRObj* p = new QRObj();
  p->cplx = A.cplx;
  if (A.cplx) p->c.compute(A.c); else p->r.compute(A.r);
  return set_dim_attr(wrap_obj(p, "eigenQR"), (int)rows(A), (int)cols(A));
}

// [[Rcpp::export]]
SEXP reigen_qr_solve(SEXP o, SEXP b) {
  Mat B = from_sexp(b, "b");
  return with_obj<QRObj>(o, "eigenQR", [&](const auto& d) { return wrap_mat(d.solve(B)); });
}

// [[Rcpp::export]]
SEXP reigen_qr_matrixQ(SEXP o) {
  return with_obj<QRObj>(o, "eigenQR", [&](const auto& d) { return wrap_mat(d.matrixQ()); });
}

// [[Rcpp::export]]
SEXP reigen_qr_matrixR(SEXP o) {
  return with_obj<QRObj>(o, "eigenQR", [&](const auto& d) { return wrap_mat(d.matrixR()); });
}

// [[Rcpp::export]]
SEXP reigen_qr_reconstructed(SEXP o) {
  return with_obj<QRObj>(o, "eigenQR", [&](const auto& d) { return wrap_mat(d.reconstructed()); });
}

// [[Rcpp::export]]
double reigen_qr_log_abs_det(SEXP o) {
  QRObj* p = obj_ptr<QRObj>(o, "eigenQR", "eigenQR");
  return p->cplx ? p->c.log_abs_det() : p->r.log_abs_det();
}

// ============================================================== ColPivHouseholderQR
//
// 带列主元，能处理秩亏的最小二乘。Eigen 文档明确推荐：不知道矩阵是否满秩时
// 用 ColPivHouseholderQR 而不是 HouseholderQR。

template <class S>
struct ColPivQRImpl {
  typedef S Scalar;
  typedef Eigen::Matrix<S, -1, -1> M;
  Eigen::ColPivHouseholderQR<M> d;

  void compute(const M& a) { d.compute(a); }
  Mat solve(const Mat& b) const { return make(M(d.solve(to_eigen<S>(b)))); }
  double rank() const { return (double)d.rank(); }
  double max_pivot() const { return (double)std::abs(d.maxPivot()); }
  double kernel_dim() const { return (double)d.dimensionOfKernel(); }
  bool invertible() const { return d.isInvertible(); }
  double log_abs_det() const { return (double)d.logAbsDeterminant(); }

  Mat matrixQ() const {
    M I = M::Identity(d.rows(), d.rows());
    return make(M(d.householderQ() * I));
  }

  Mat matrixR() const {
    M R = M::Zero(d.rows(), d.cols());
    R = d.matrixQR().template triangularView<Eigen::Upper>();
    return make(R);
  }

  Mat permutationP() const {
    M P = M::Identity(d.cols(), d.cols());
    P = d.colsPermutation();
    return make(P);
  }

  Mat reconstructed() const {
    M Q = M::Identity(d.rows(), d.rows());
    Q = d.householderQ() * Q;
    M R = M::Zero(d.rows(), d.cols());
    R = d.matrixQR().template triangularView<Eigen::Upper>();
    M P = M::Identity(d.cols(), d.cols());
    P = d.colsPermutation();
    return make(M(Q * R * P.transpose()));
  }
};

REIGEN_DUAL_OBJ(ColPivQRObj, ColPivQRImpl);

// [[Rcpp::export]]
SEXP reigen_colpivqr_new(SEXP a) {
  Mat A = from_sexp_finite(a, "a");
  ColPivQRObj* p = new ColPivQRObj();
  p->cplx = A.cplx;
  if (A.cplx) p->c.compute(A.c); else p->r.compute(A.r);
  return set_dim_attr(wrap_obj(p, "eigenColPivQR"), (int)rows(A), (int)cols(A));
}

// [[Rcpp::export]]
SEXP reigen_colpivqr_solve(SEXP o, SEXP b) {
  Mat B = from_sexp(b, "b");
  return with_obj<ColPivQRObj>(o, "eigenColPivQR", [&](const auto& d) { return wrap_mat(d.solve(B)); });
}

// [[Rcpp::export]]
SEXP reigen_colpivqr_matrixQ(SEXP o) {
  return with_obj<ColPivQRObj>(o, "eigenColPivQR", [&](const auto& d) { return wrap_mat(d.matrixQ()); });
}

// [[Rcpp::export]]
SEXP reigen_colpivqr_matrixR(SEXP o) {
  return with_obj<ColPivQRObj>(o, "eigenColPivQR", [&](const auto& d) { return wrap_mat(d.matrixR()); });
}

// [[Rcpp::export]]
SEXP reigen_colpivqr_permutationP(SEXP o) {
  return with_obj<ColPivQRObj>(o, "eigenColPivQR", [&](const auto& d) { return wrap_mat(d.permutationP()); });
}

// [[Rcpp::export]]
SEXP reigen_colpivqr_reconstructed(SEXP o) {
  return with_obj<ColPivQRObj>(o, "eigenColPivQR", [&](const auto& d) { return wrap_mat(d.reconstructed()); });
}

// [[Rcpp::export]]
double reigen_colpivqr_rank(SEXP o) {
  ColPivQRObj* p = obj_ptr<ColPivQRObj>(o, "eigenColPivQR", "eigenColPivQR");
  return p->cplx ? p->c.rank() : p->r.rank();
}

// ColPivHouseholderQR 没有 rcond()。要条件数请用 eigen_full_lu()，
// 那里走 FullPivLU::rcond()；这里能给的奇异程度指标是 max_pivot 和 rank。

// [[Rcpp::export]]
double reigen_colpivqr_max_pivot(SEXP o) {
  ColPivQRObj* p = obj_ptr<ColPivQRObj>(o, "eigenColPivQR", "eigenColPivQR");
  return p->cplx ? p->c.max_pivot() : p->r.max_pivot();
}

// [[Rcpp::export]]
double reigen_colpivqr_kernel_dim(SEXP o) {
  ColPivQRObj* p = obj_ptr<ColPivQRObj>(o, "eigenColPivQR", "eigenColPivQR");
  return p->cplx ? p->c.kernel_dim() : p->r.kernel_dim();
}

// [[Rcpp::export]]
bool reigen_colpivqr_invertible(SEXP o) {
  ColPivQRObj* p = obj_ptr<ColPivQRObj>(o, "eigenColPivQR", "eigenColPivQR");
  return p->cplx ? p->c.invertible() : p->r.invertible();
}

// [[Rcpp::export]]
double reigen_colpivqr_log_abs_det(SEXP o) {
  ColPivQRObj* p = obj_ptr<ColPivQRObj>(o, "eigenColPivQR", "eigenColPivQR");
  return p->cplx ? p->c.log_abs_det() : p->r.log_abs_det();
}

// ============================================================== SVD
//
// 用 ThinU | ThinV：U 和 V 只算到 min(n,p) 列。全 U 在 n >> p 时
// 会多出几 GB 的零空间，几乎从来不是用户想要的。
//
// 两套后端：JacobiSVD（精度好、O(n^3) 但常数极大）和 BDCSVD（分治，n>=200 时
// 快一个数量级）。两者都继承 SVDBase，threshold() / rank() / nonzeroSingularValues()
// 是同一个实现（SVDBase.h），所以换后端不会改变秩判定和伪逆的截断口径——
// 这是敢让它们互换的前提。实测 n=400 时奇异值相对误差 1e-13，threshold 完全相同。
//
// 这里用虚函数在两者之间二选一：每次方法调用多一次间接跳转，相对 O(n^3) 的
// 分解成本可以忽略，而它让 SVDObj 不必按 4 种（实/复 x 2 后端）展开。

template <class S>
struct SVDCore {
  typedef Eigen::Matrix<S, -1, -1> M;
  virtual ~SVDCore() {}
  virtual Mat solve(const Mat& b) const = 0;
  virtual Mat singularValues() const = 0;
  virtual Mat matrixU() const = 0;
  virtual Mat matrixV() const = 0;
  virtual Mat pseudoInverse() const = 0;
  virtual Mat reconstructed() const = 0;
  virtual double rank() const = 0;
  virtual double nonzero() const = 0;
  virtual double threshold() const = 0;
};

// 两个后端暴露的接口完全一致（都来自 SVDBase），所以主体只有一份实现。
template <class S, class Solver>
struct SVDCoreImpl : SVDCore<S> {
  typedef typename SVDCore<S>::M M;
  Solver d;

  // Eigen 5 把 ComputeThinU/ComputeThinV 从 compute() 的运行时参数搬到了
  // 类模板参数上：旧的 d.compute(a, options) 仍能用，但每次实例化都会
  // 冒 -Wdeprecated-declarations。这里按新写法把选项写进类型，compute 只收矩阵。
  void compute(const M& a) { d.compute(a); }

  Mat solve(const Mat& b) const override { return make(M(d.solve(to_eigen<S>(b)))); }
  Mat singularValues() const override { return make(M(M(d.singularValues()))); }
  Mat matrixU() const override { return make(M(d.matrixU())); }
  Mat matrixV() const override { return make(M(d.matrixV())); }
  double rank() const override { return (double)d.rank(); }
  double nonzero() const override { return (double)d.nonzeroSingularValues(); }
  double threshold() const override { return (double)d.threshold(); }

  // 两个后端都没有 pseudoInverse()。按 V * diag(1/s) * U^H 自己拼。
  //
  // 关键：threshold() 返回的是**相对**阈值（epsilon * diagSize），
  // 要和最大奇异值相乘才是绝对阈值——rank() 内部就是这么判的。
  // 直接拿 threshold() 当绝对值比，会把 1e-15 量级的零奇异值当成有效值
  // 求倒数，伪逆直接炸成 1e14。
  Mat pseudoInverse() const override {
    const Eigen::Index r = d.singularValues().rows();
    if (r == 0) return make(M(d.matrixV() * M::Zero(0, 0) * d.matrixU().adjoint()));
    const double abs_tol = (double)d.threshold() * (double)d.singularValues()(0);
    M Sinv = M::Zero(r, r);
    for (Eigen::Index i = 0; i < r; ++i) {
      const double s = (double)d.singularValues()(i);
      Sinv(i, i) = (s > abs_tol) ? (1.0 / s) : 0.0;
    }
    return make(M(d.matrixV() * Sinv * d.matrixU().adjoint()));
  }

  Mat reconstructed() const override {
    const Eigen::Index r = d.singularValues().rows();
    M Sig = M::Zero(r, r);
    for (Eigen::Index i = 0; i < r; ++i) Sig(i, i) = d.singularValues()(i);
    return make(M(d.matrixU() * Sig * d.matrixV().adjoint()));
  }
};

// "auto" 的门槛：低于这个尺寸 JacobiSVD 本来就快到无从优化，走久经考验的老路；
// 高于它 BDCSVD 的优势迅速拉开（n=100 约 4x，n=200 约 8x，n=400 约 16x）。
static const Eigen::Index k_bdc_min_dim = 16;

template <class S>
struct SVDImpl {
  typedef S Scalar;
  typedef Eigen::Matrix<S, -1, -1> M;

  std::unique_ptr<SVDCore<S>> core;
  bool bdc = false;

  // method: 0 = auto, 1 = jacobi, 2 = bdc
  void compute(const M& a, int method) {
    const Eigen::Index diag = std::min(a.rows(), a.cols());
    bdc = (method == 2) || (method == 0 && diag >= k_bdc_min_dim);
    // 用 unique_ptr 而不是裸 new：compute() 里抛 bad_alloc 时裸指针会漏，
    // 和 LLT 那边 stop 之前必须 delete 是同一个道理。
    if (bdc) {
      typedef Eigen::BDCSVD<M, Eigen::ComputeThinU | Eigen::ComputeThinV> Solver;
      std::unique_ptr<SVDCoreImpl<S, Solver>> p(new SVDCoreImpl<S, Solver>());
      p->compute(a);
      core = std::move(p);
    } else {
      typedef Eigen::JacobiSVD<M, Eigen::ComputeThinU | Eigen::ComputeThinV> Solver;
      std::unique_ptr<SVDCoreImpl<S, Solver>> p(new SVDCoreImpl<S, Solver>());
      p->compute(a);
      core = std::move(p);
    }
  }

  // 对象一旦构造完，core 就必然非空；这里只是让调用点少写一次判空。
  const SVDCore<S>& get() const {
    if (!core) Rcpp::stop("internal error: SVD object was not initialised");
    return *core;
  }

  Mat solve(const Mat& b) const { return get().solve(b); }
  Mat singularValues() const { return get().singularValues(); }
  Mat matrixU() const { return get().matrixU(); }
  Mat matrixV() const { return get().matrixV(); }
  Mat pseudoInverse() const { return get().pseudoInverse(); }
  Mat reconstructed() const { return get().reconstructed(); }
  double rank() const { return get().rank(); }
  double nonzero() const { return get().nonzero(); }
  double threshold() const { return get().threshold(); }
};

REIGEN_DUAL_OBJ(SVDObj, SVDImpl);

// [[Rcpp::export]]
SEXP reigen_svd_new(SEXP a, int method) {
  Mat A = from_sexp_finite(a, "a");
  if (rows(A) == 0 || cols(A) == 0) Rcpp::stop("SVD of an empty matrix");
  if (method < 0 || method > 2) Rcpp::stop("unknown SVD method code %d", method);
  SVDObj* p = new SVDObj();
  p->cplx = A.cplx;
  if (A.cplx) p->c.compute(A.c, method); else p->r.compute(A.r, method);
  return set_dim_attr(wrap_obj(p, "eigenSVD"), (int)rows(A), (int)cols(A));
}

// [[Rcpp::export]]
SEXP reigen_svd_backend(SEXP o) {
  SVDObj* p = obj_ptr<SVDObj>(o, "eigenSVD", "eigenSVD");
  return Rcpp::CharacterVector::create(p->cplx ? (p->c.bdc ? "bdc" : "jacobi")
                                               : (p->r.bdc ? "bdc" : "jacobi"));
}

// [[Rcpp::export]]
SEXP reigen_svd_solve(SEXP o, SEXP b) {
  Mat B = from_sexp(b, "b");
  return with_obj<SVDObj>(o, "eigenSVD", [&](const auto& d) { return wrap_mat(d.solve(B)); });
}

// [[Rcpp::export]]
SEXP reigen_svd_singular_values(SEXP o) {
  return with_obj<SVDObj>(o, "eigenSVD", [&](const auto& d) { return wrap_mat(d.singularValues()); });
}

// [[Rcpp::export]]
SEXP reigen_svd_matrixU(SEXP o) {
  return with_obj<SVDObj>(o, "eigenSVD", [&](const auto& d) { return wrap_mat(d.matrixU()); });
}

// [[Rcpp::export]]
SEXP reigen_svd_matrixV(SEXP o) {
  return with_obj<SVDObj>(o, "eigenSVD", [&](const auto& d) { return wrap_mat(d.matrixV()); });
}

// [[Rcpp::export]]
SEXP reigen_svd_pseudo_inverse(SEXP o) {
  return with_obj<SVDObj>(o, "eigenSVD", [&](const auto& d) { return wrap_mat(d.pseudoInverse()); });
}

// [[Rcpp::export]]
SEXP reigen_svd_reconstructed(SEXP o) {
  return with_obj<SVDObj>(o, "eigenSVD", [&](const auto& d) { return wrap_mat(d.reconstructed()); });
}

// [[Rcpp::export]]
double reigen_svd_rank(SEXP o) {
  SVDObj* p = obj_ptr<SVDObj>(o, "eigenSVD", "eigenSVD");
  return p->cplx ? p->c.rank() : p->r.rank();
}

// [[Rcpp::export]]
double reigen_svd_nonzero(SEXP o) {
  SVDObj* p = obj_ptr<SVDObj>(o, "eigenSVD", "eigenSVD");
  return p->cplx ? p->c.nonzero() : p->r.nonzero();
}

// [[Rcpp::export]]
double reigen_svd_threshold(SEXP o) {
  SVDObj* p = obj_ptr<SVDObj>(o, "eigenSVD", "eigenSVD");
  return p->cplx ? p->c.threshold() : p->r.threshold();
}

// ============================================================== SelfAdjointEigenSolver
//
// 对称（Hermitian）矩阵的特征分解：A = V D V^H，D 是实特征值。
// Eigen 只读下三角，上三角被忽略——传进来的矩阵不对称不会报错。
// 一般矩阵请用 ComplexEigenSolver。

template <class S>
struct SelfAdjEigImpl {
  typedef S Scalar;
  typedef Eigen::Matrix<S, -1, -1> M;
  typedef Eigen::SelfAdjointEigenSolver<M> Solver;
  Solver d;

  void compute(const M& a) { d.compute(a); }
  Mat eigenvalues() const { return make(Eigen::MatrixXd(d.eigenvalues())); }
  Mat eigenvectors() const { return make(M(d.eigenvectors())); }
  bool ok() const { return d.info() == Eigen::Success; }

  // 谱分解求解：x = V D^-1 V^H b。比走 LU 稳，但只在 A 对称时成立。
  Mat solve(const Mat& b) const {
    M B = to_eigen<S>(b);
    const Eigen::Index n = d.eigenvalues().rows();
    Eigen::MatrixXd inv(n, 1);
    for (Eigen::Index i = 0; i < n; ++i) {
      const double v = (double)d.eigenvalues()(i);
      inv(i, 0) = (v == 0.0) ? 0.0 : 1.0 / v;
    }
    M V = d.eigenvectors();
    return make(M(V * inv.asDiagonal() * V.adjoint() * B));
  }
};

REIGEN_DUAL_OBJ(SelfAdjEigObj, SelfAdjEigImpl);

// [[Rcpp::export]]
SEXP reigen_eigen_new(SEXP a) {
  Mat A = from_sexp_finite(a, "a");
  need_square(A, "SelfAdjointEigenSolver");
  SelfAdjEigObj* p = new SelfAdjEigObj();
  p->cplx = A.cplx;
  if (A.cplx) p->c.compute(A.c); else p->r.compute(A.r);
  return set_dim_attr(wrap_obj(p, "eigenSelfAdjointEigen"), (int)rows(A), (int)cols(A));
}

// [[Rcpp::export]]
SEXP reigen_eigen_values(SEXP o) {
  return with_obj<SelfAdjEigObj>(o, "eigenSelfAdjointEigen",
                                 [&](const auto& d) { return wrap_mat(d.eigenvalues()); });
}

// [[Rcpp::export]]
SEXP reigen_eigen_vectors(SEXP o) {
  return with_obj<SelfAdjEigObj>(o, "eigenSelfAdjointEigen",
                                 [&](const auto& d) { return wrap_mat(d.eigenvectors()); });
}

// [[Rcpp::export]]
SEXP reigen_eigen_solve(SEXP o, SEXP b) {
  Mat B = from_sexp(b, "b");
  return with_obj<SelfAdjEigObj>(o, "eigenSelfAdjointEigen",
                                 [&](const auto& d) { return wrap_mat(d.solve(B)); });
}

// [[Rcpp::export]]
bool reigen_eigen_ok(SEXP o) {
  SelfAdjEigObj* p = obj_ptr<SelfAdjEigObj>(o, "eigenSelfAdjointEigen", "eigenSelfAdjointEigen");
  return p->cplx ? p->c.ok() : p->r.ok();
}

// ============================================================== ComplexEigenSolver
//
// 一般（非对称）矩阵的特征分解。实矩阵也走这条：输入转成复数，
// 输出恒为复数。Eigen 的 EigenSolver<MatrixXd> 能给实矩阵更准的
// 实特征值，但接口要分两套，收益不值这个复杂度——特征值本来就是复数域的事。

struct CplxEigObj {
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> d;
};

// [[Rcpp::export]]
SEXP reigen_ceigen_new(SEXP a) {
  Mat A = from_sexp_finite(a, "a");
  need_square(A, "ComplexEigenSolver");
  CplxEigObj* p = new CplxEigObj();
  p->d.compute(as_cplx(A));
  return set_dim_attr(wrap_obj(p, "eigenComplexEigen"), (int)rows(A), (int)cols(A));
}

// [[Rcpp::export]]
SEXP reigen_ceigen_values(SEXP o) {
  CplxEigObj* p = obj_ptr<CplxEigObj>(o, "eigenComplexEigen", "eigenComplexEigen");
  return wrap_mat(make(Eigen::MatrixXcd(p->d.eigenvalues())));
}

// [[Rcpp::export]]
SEXP reigen_ceigen_vectors(SEXP o) {
  CplxEigObj* p = obj_ptr<CplxEigObj>(o, "eigenComplexEigen", "eigenComplexEigen");
  return wrap_mat(make(Eigen::MatrixXcd(p->d.eigenvectors())));
}

// [[Rcpp::export]]
bool reigen_ceigen_ok(SEXP o) {
  CplxEigObj* p = obj_ptr<CplxEigObj>(o, "eigenComplexEigen", "eigenComplexEigen");
  return p->d.info() == Eigen::Success;
}
