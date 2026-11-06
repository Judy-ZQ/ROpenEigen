// reigen_sparse.cpp —— 稀疏矩阵
//
// 与稠密侧的关系：eigenMat 包 Eigen::MatrixXd，eigenSparse 包
// Eigen::SparseMatrix<double>。两者都是 EXTPTRSXP，都带 "eigenObject"
// 作为公共尾类，所以 $ 方法和 print 可以共用一套。
//
// 只做实数。复数稀疏（SparseMatrix<complex<double>>）Eigen 本身支持，
// 但转换层、求解器和测试都要再来一遍，暂时不开口子；传复数进来会明确报错，
// 而不是静默取实部——后者会悄悄丢掉虚部，是最糟的一种"能跑"。
//
// 索引约定：R 侧一律 1-based（和 Matrix::sparseMatrix 一致），进 C++ 立刻
// 转成 Eigen 的 0-based。这是稀疏代码里最容易错的地方，所以转换只发生在
// 两个函数里（reigen_sparse_triplet / reigen_sparse_triplets），别处一律
// 用 Eigen 自己的坐标。

#include "reigen.h"

#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>
// GMRES / MINRES / DGMRES / IDR(s) 只在 unsupported 里，稠密头文件不含它们。
#include <unsupported/Eigen/IterativeSolvers>

#include <vector>

using namespace reigen;

namespace {

// 列主序，与 R 的 dgCMatrix 一致。行主序会让 as.matrix() 的转换慢一截，
// 而且 Matrix 包交换数据时还要转置。
typedef Eigen::SparseMatrix<double> SpMat;

struct SparseObj {
  SpMat m;
};

inline SparseObj* sp_ptr(SEXP x, const char* what) {
  return obj_ptr<SparseObj>(x, "eigenSparse", what);
}

// 包成 R 对象，并把尺寸挂到 "eigenDim" 属性上（与分解对象同一套约定）。
inline SEXP wrap_sparse(SpMat m) {
  SparseObj* p = new SparseObj();
  p->m = std::move(m);
  // 存了多少个非零，就占多少 (double + int)。这是给 gc_tick 的粗略估计。
  const size_t bytes = (size_t)p->m.nonZeros() * (sizeof(double) + sizeof(int)) + sizeof(SparseObj);
  SEXP r = wrap_obj(p, "eigenSparse", bytes);
  return set_dim_attr(r, (int)p->m.rows(), (int)p->m.cols());
}

// 稀疏矩阵的形状校验。Eigen 的 + / - / * 在 NDEBUG 下不做检查（eigen_assert
// 被编译掉了），维度不对会直接越界读，所以每一处都要自己判。
inline void need_same_shape(const SpMat& a, const SpMat& b, const char* op) {
  if (a.rows() != b.rows() || a.cols() != b.cols())
    Rcpp::stop("%s: dimensions differ (%lldx%lld vs %lldx%lld)", op,
               (long long)a.rows(), (long long)a.cols(),
               (long long)b.rows(), (long long)b.cols());
}

// 求解前拒绝非有限值。稠密侧是靠 from_sexp_finite() 在入口处拦的，稀疏矩阵
// 走的是自己的转换路径，绕过了那道关——所以这里必须单独查一遍。
// 不查的话 NaN/Inf 会一路传进分解，SparseLU 只会给出一个 NaN 解。
inline void need_finite(const SpMat& m, const char* what) {
  for (Eigen::Index c = 0; c < m.outerSize(); ++c)
    for (SpMat::InnerIterator it(m, c); it; ++it)
      if (!std::isfinite(it.value()))
        Rcpp::stop("%s: sparse matrix contains NaN or Inf; solvers reject such "
                   "input, as base R's solve() does -- drop or impute the "
                   "missing values first", what);
}

}  // namespace

// ============================================================ 构造

// 1-based 三元组 -> 稀疏矩阵。重复的 (i, j) 会被**相加**，与
// Matrix::sparseMatrix() 的默认行为一致（那边叫 "use.last" 的反面）。
// [[Rcpp::export]]
SEXP reigen_sparse_triplet(SEXP i_, SEXP j_, SEXP x_, SEXP nrow_, SEXP ncol_) {
  const int nr = Rcpp::as<int>(nrow_);
  const int nc = Rcpp::as<int>(ncol_);
  if (nr < 0 || nc < 0)
    Rcpp::stop("nrow and ncol must be non-negative, got %d and %d", nr, nc);

  Rcpp::IntegerVector ii(i_), jj(j_);
  Rcpp::NumericVector xx(x_);
  const int n = (int)xx.size();
  if ((int)ii.size() != n || (int)jj.size() != n)
    Rcpp::stop("i, j and x must have the same length, got %d, %d and %d",
               (int)ii.size(), (int)jj.size(), n);

  std::vector<Eigen::Triplet<double>> trip;
  trip.reserve((size_t)n);
  for (int k = 0; k < n; ++k) {
    if (ii[k] == NA_INTEGER || jj[k] == NA_INTEGER)
      Rcpp::stop("i and j must not contain NA (element %d)", k + 1);
    const int r = (int)ii[k] - 1;      // R 1-based -> Eigen 0-based
    const int c = (int)jj[k] - 1;
    if (r < 0 || r >= nr || c < 0 || c >= nc)
      Rcpp::stop("index out of bounds at element %d: (%d, %d) is not inside "
                 "[1, %d] x [1, %d]", k + 1, (int)ii[k], (int)jj[k], nr, nc);
    trip.push_back(Eigen::Triplet<double>(r, c, (double)xx[k]));
  }

  SpMat m(nr, nc);
  m.setFromTriplets(trip.begin(), trip.end());
  m.makeCompressed();
  return wrap_sparse(std::move(m));
}

// 稠密 -> 稀疏，丢掉**精确**的 0（与 Matrix(A, sparse = TRUE) 一致）。
// [[Rcpp::export]]
SEXP reigen_sparse_from_dense(SEXP a) {
  Mat A = from_sexp(a, "a");
  if (A.cplx)
    Rcpp::stop("complex sparse matrices are not supported yet; "
               "split into Re() and Im() and keep two real sparse matrices");

  std::vector<Eigen::Triplet<double>> trip;
  for (Eigen::Index c = 0; c < cols(A); ++c)
    for (Eigen::Index r = 0; r < rows(A); ++r) {
      const double v = A.r(r, c);
      // NaN != 0 为真，所以 NaN 会被保留下来，不会被当成 0 丢掉。
      if (v != 0.0) trip.push_back(Eigen::Triplet<double>((int)r, (int)c, v));
    }

  SpMat m((int)rows(A), (int)cols(A));
  m.setFromTriplets(trip.begin(), trip.end());
  m.makeCompressed();
  return wrap_sparse(std::move(m));
}

// ============================================================ 查询

// [[Rcpp::export]]
Rcpp::IntegerVector reigen_sparse_dim(SEXP o) {
  SparseObj* p = sp_ptr(o, "x");
  return Rcpp::IntegerVector::create((int)p->m.rows(), (int)p->m.cols());
}

// 存储的条目数。注意是**存储的**，显式零也算在内——这与 Matrix 包的 nnz()
// 一致，而不是"数值上非零的个数"。
// [[Rcpp::export]]
double reigen_sparse_nnz(SEXP o) {
  return (double)sp_ptr(o, "x")->m.nonZeros();
}

// [[Rcpp::export]]
SEXP reigen_sparse_as_matrix(SEXP o) {
  SparseObj* p = sp_ptr(o, "x");
  return wrap_mat(make(Eigen::MatrixXd(p->m)));
}

// 导出成 1-based 三元组，列优先顺序（与 Matrix 的 summary() 一致）。
// [[Rcpp::export]]
SEXP reigen_sparse_triplets(SEXP o) {
  SparseObj* p = sp_ptr(o, "x");
  const SpMat& m = p->m;
  const Eigen::Index n = m.nonZeros();
  Rcpp::IntegerVector ii((R_xlen_t)n), jj((R_xlen_t)n);
  Rcpp::NumericVector xx((R_xlen_t)n);
  Eigen::Index k = 0;
  for (Eigen::Index c = 0; c < m.outerSize(); ++c)
    for (SpMat::InnerIterator it(m, c); it; ++it) {
      ii[(R_xlen_t)k] = (int)it.row() + 1;    // Eigen 0-based -> R 1-based
      jj[(R_xlen_t)k] = (int)it.col() + 1;
      xx[(R_xlen_t)k] = it.value();
      ++k;
    }
  return Rcpp::List::create(Rcpp::Named("i") = ii,
                            Rcpp::Named("j") = jj,
                            Rcpp::Named("x") = xx);
}

// 单个元素，1-based。越界报错而不是返回 0——静默返回 0 会把"取错了"伪装成
// "那里正好是 0"。
// [[Rcpp::export]]
double reigen_sparse_get(SEXP o, int i, int j) {
  SparseObj* p = sp_ptr(o, "x");
  if (i < 1 || i > p->m.rows() || j < 1 || j > p->m.cols())
    Rcpp::stop("index out of bounds: (%d, %d) is not inside [1, %lld] x [1, %lld]",
               i, j, (long long)p->m.rows(), (long long)p->m.cols());
  return p->m.coeff(i - 1, j - 1);
}

// 非有限值检测。code: 0 全部有限, 1 有 NaN, 2 有 Inf —— 与稠密侧
// reigen_finite() 的编码一致，R 侧三个谓词共用。
// 稀疏矩阵只检查**存储的**条目：没存的位置按定义就是 0，必然有限。
//
// code 1 用 std::isnan 而不是 R_IsNaN：R_IsNaN() 对 NA 返回 false（R 的约定是
// R_IsNA 认 NA、R_IsNaN 只认"非 NA 的 NaN"），而稠密侧 reigen_finite 用的是
// Eigen 的 isNaN()，也就是 std::isnan，对 NA 一样返回 true。两边口径必须一致，
// 否则 eigen_has_nan(As) 和 eigen_has_nan(as.matrix(As)) 会对同一个矩阵给出
// 不同答案——这正是 consistency.R 存在的理由。
// [[Rcpp::export]]
bool reigen_sparse_finite(SEXP o, int code) {
  SparseObj* p = sp_ptr(o, "x");
  for (Eigen::Index c = 0; c < p->m.outerSize(); ++c)
    for (SpMat::InnerIterator it(p->m, c); it; ++it) {
      const double v = it.value();
      if (code == 0 && !std::isfinite(v)) return false;
      if (code == 1 && std::isnan(v)) return true;
      if (code == 2 && std::isinf(v)) return true;
    }
  return code == 0;
}

// ============================================================ 运算

// [[Rcpp::export]]
SEXP reigen_sparse_transpose(SEXP o) {
  SparseObj* p = sp_ptr(o, "x");
  SpMat out = p->m.transpose();   // 转置表达式，赋给 SpMat 才真正物化
  out.makeCompressed();
  return wrap_sparse(std::move(out));
}

// 稀疏 x 稠密 -> 稠密
// [[Rcpp::export]]
SEXP reigen_sparse_matmul(SEXP a, SEXP b) {
  SparseObj* pa = sp_ptr(a, "a");
  Mat B = from_sexp(b, "b");
  if (B.cplx)
    Rcpp::stop("complex operands are not supported in sparse products yet");
  if (pa->m.cols() != rows(B))
    Rcpp::stop("non-conformable arguments: %lldx%lld %%*%% %lldx%lld",
               (long long)pa->m.rows(), (long long)pa->m.cols(),
               (long long)rows(B), (long long)cols(B));
  return wrap_mat(make(Eigen::MatrixXd(pa->m * B.r)));
}

// 稀疏 x 稀疏 -> 稀疏
// [[Rcpp::export]]
SEXP reigen_sparse_matmul_sparse(SEXP a, SEXP b) {
  SparseObj* pa = sp_ptr(a, "a");
  SparseObj* pb = sp_ptr(b, "b");
  if (pa->m.cols() != pb->m.rows())
    Rcpp::stop("non-conformable arguments: %lldx%lld %%*%% %lldx%lld",
               (long long)pa->m.rows(), (long long)pa->m.cols(),
               (long long)pb->m.rows(), (long long)pb->m.cols());
  SpMat out = pa->m * pb->m;
  out.makeCompressed();
  return wrap_sparse(std::move(out));
}

// [[Rcpp::export]]
SEXP reigen_sparse_add(SEXP a, SEXP b) {
  SparseObj* pa = sp_ptr(a, "a");
  SparseObj* pb = sp_ptr(b, "b");
  need_same_shape(pa->m, pb->m, "+");
  SpMat out = pa->m + pb->m;
  out.makeCompressed();
  return wrap_sparse(std::move(out));
}

// [[Rcpp::export]]
SEXP reigen_sparse_sub(SEXP a, SEXP b) {
  SparseObj* pa = sp_ptr(a, "a");
  SparseObj* pb = sp_ptr(b, "b");
  need_same_shape(pa->m, pb->m, "-");
  SpMat out = pa->m - pb->m;
  out.makeCompressed();
  return wrap_sparse(std::move(out));
}

// [[Rcpp::export]]
SEXP reigen_sparse_scale(SEXP a, double s) {
  SparseObj* pa = sp_ptr(a, "a");
  SpMat out = pa->m * s;
  out.makeCompressed();
  return wrap_sparse(std::move(out));
}

// Summary 组要的东西一次算完。关键是**隐含零也要参与**：
//   prod 只要有一个隐含零就是 0；
//   max/min 要把 0 当候选值（全正矩阵的 min 是 0，不是最小的存值）。
// 显式存下来的零和隐含零在数值上没有区别，所以 nz 只统计真正非零的存值。
//
// NA 和 NaN 要分开数。base R 的 sum(c(1, NA)) 是 NA，sum(c(1, NaN)) 是 NaN，
// 两者 is.na() 都为真但 identical() 为假。NA_real_ 在 C++ 里也是一个 NaN，
// 只有 R_IsNA() 能把它和普通 NaN 区分开，所以不能只用 std::isnan()。
// [[Rcpp::export]]
Rcpp::NumericVector reigen_sparse_stats(SEXP a) {
  const SpMat& m = sp_ptr(a, "a")->m;
  const double total = (double)m.rows() * (double)m.cols();
  double sum = 0.0, prod = 1.0, vmin = 0.0, vmax = 0.0;
  double nnz = 0.0, nz = 0.0, nna = 0.0, nnan = 0.0;
  double sumsq = 0.0, trace = 0.0;
  double trace_na = 0.0, trace_nan = 0.0;
  // Frobenius 范数单独用一套缩放的累加器，理由见 frobenius_scaled()：
  // 朴素的平方和在 1e-200 上给 0、在 1e300 上给 Inf，两者都不报错。
  // 注意 sumsq 保持朴素——它对应的是 base 的 sum(x^2)，那个**就是**会溢出的，
  // 所以 eigen_sumsq() 应当和 base 一样给 Inf，只有范数需要缩放。
  double fmax = 0.0, fssq = 1.0;
  bool f_inf = false;
  bool have = false;
  for (Eigen::Index c = 0; c < m.outerSize(); ++c)
    for (SpMat::InnerIterator it(m, c); it; ++it) {
      const double v = it.value();
      ++nnz;
      // 迹只碰主对角线，所以对角线上的缺失值要单独数。全矩阵的 na/nan 计数
      // 回答不了"这条对角线上有没有缺失"——非对角线上的 NA 对迹没有影响。
      const bool on_diag = (it.row() == c);
      if (R_IsNA(v))  { ++nna;  if (on_diag) ++trace_na;  continue; }
      if (R_IsNaN(v)) { ++nnan; if (on_diag) ++trace_nan; continue; }
      sum += v;
      prod *= v;
      sumsq += v * v;
      const double av = std::abs(v);
      if (av == R_PosInf) f_inf = true;
      else if (av != 0.0) {
        if (av > fmax) { const double r = fmax / av; fssq = 1.0 + fssq * r * r; fmax = av; }
        else           { const double r = av / fmax; fssq += r * r; }
      }
      // 列主序：it.row() 是行号、c 是列号（本类型是 column-major）。
      if (on_diag) trace += v;
      if (v != 0.0) ++nz;
      if (!have) { vmin = vmax = v; have = true; }
      else { if (v < vmin) vmin = v; if (v > vmax) vmax = v; }
    }

  // 未被显式存储的位置就是零。存在隐含零时：
  //   prod 归零；min/max 必须把 0 拉进候选。
  // sumsq 和 trace 不受影响：隐含零对平方和没有贡献，也不在主对角线上
  // （主对角线上的零本来就"是"零，加不加一样）。
  if (total - nnz > 0.0) {
    prod = 0.0;
    if (!have) { vmin = vmax = 0.0; }
    else { if (0.0 < vmin) vmin = 0.0; if (0.0 > vmax) vmax = 0.0; }
  }

  return Rcpp::NumericVector::create(
    Rcpp::Named("nnz")   = nnz,
    Rcpp::Named("total") = total,
    Rcpp::Named("sum")   = sum,
    Rcpp::Named("prod")  = prod,
    Rcpp::Named("min")   = vmin,
    Rcpp::Named("max")   = vmax,
    Rcpp::Named("nz")    = nz,
    Rcpp::Named("na")    = nna,
    Rcpp::Named("nan")   = nnan,
    Rcpp::Named("sumsq") = sumsq,
    Rcpp::Named("frob")  = f_inf ? R_PosInf : fmax * std::sqrt(fssq),
    Rcpp::Named("trace") = trace,
    Rcpp::Named("trace_na")  = trace_na,
    Rcpp::Named("trace_nan") = trace_nan);
}

// 极值的位置。op: 0 = min，1 = max。返回 list(value, row, col)，1-based。
//
// 稀疏版的难点全在**隐含零**上：矩阵的大部分格子没存，那些格子的值是 0，
// 而它们完全有权成为最小值。所以候选集合是"存下来的有限值 ∪ {0}（如果存在
// 隐含零）"，而不是"存下来的有限值"。
//
// 找位置分两步。第一步扫一遍得到极值。第二步按列主序找第一个等于极值的位置：
// 逐列走，列内按行升序；每个存储项之前先看有没有"空档"——上一个存储行号与
// 当前行号之间有间隔，那第一个空档就是隐含零。空档只在极值等于 0 时才算答案，
// 因为非零的极值不可能出现在空档里。整列走完还要看列尾的空档。
//
// 复杂度 O(nnz + cols)，不稠密化。跟 which.min(as.matrix(x)) 的结果一致——
// 这正是这类函数该守的契约。
//
// 关于"这答案有用吗"：对一个大部分为零的矩阵，最小值几乎总是 0，位置几乎
// 总是列主序第一个未被存储的格子。信息量不大，但它是**对的**，而且和稠密
// 那一扇门给的是同一个答案。以前这里直接报错，理由是"位置取决于存储模式
// 而不是数值"——位置确实取决于存储模式，可它是唯一正确的答案，拒绝它只是
// 让用户自己去 as.matrix() 绕一圈。
// [[Rcpp::export]]
SEXP reigen_sparse_which(SEXP a, int op) {
  const SpMat& m = sp_ptr(a, "a")->m;
  const bool want_max = (op == 1);
  const double total = (double)m.rows() * (double)m.cols();
  const bool has_implicit = (total > (double)m.nonZeros());

  // 第一步：极值。ISNAN 而不是 !R_FINITE，Inf 是合法取值。
  bool found = false;
  double best = 0.0;
  for (Eigen::Index c = 0; c < m.outerSize(); ++c)
    for (SpMat::InnerIterator it(m, c); it; ++it) {
      const double v = it.value();
      if (ISNAN(v)) continue;
      if (!found || (want_max ? v > best : v < best)) { best = v; found = true; }
    }
  if (has_implicit) {
    // 0 是候选。隐含零是有限的，所以"全都缺失"在有隐含零时不可能成立。
    if (!found || (want_max ? 0.0 > best : 0.0 < best)) best = 0.0;
    found = true;
  }
  if (!found)  // 一个有限值都没有：全缺失且无隐含零（即整矩阵存满且全 NA）
    return Rcpp::List::create(Rcpp::Named("value") = NA_REAL,
                              Rcpp::Named("row")   = NA_INTEGER,
                              Rcpp::Named("col")   = NA_INTEGER);

  // 第二步：列主序第一个等于 best 的位置。
  for (Eigen::Index c = 0; c < m.cols(); ++c) {
    Eigen::Index prev = -1;               // 上一个存储行号，用于找空档
    for (SpMat::InnerIterator it(m, c); it; ++it) {
      const Eigen::Index r = it.row();
      if (has_implicit && best == 0.0 && prev + 1 < r)
        return Rcpp::List::create(Rcpp::Named("value") = best,
                                  Rcpp::Named("row")   = (int)(prev + 1) + 1,
                                  Rcpp::Named("col")   = (int)c + 1);
      const double v = it.value();
      if (!ISNAN(v) && v == best)
        return Rcpp::List::create(Rcpp::Named("value") = best,
                                  Rcpp::Named("row")   = (int)r + 1,
                                  Rcpp::Named("col")   = (int)c + 1);
      prev = r;
    }
    if (has_implicit && best == 0.0 && prev + 1 < m.rows())
      return Rcpp::List::create(Rcpp::Named("value") = best,
                                Rcpp::Named("row")   = (int)(prev + 1) + 1,
                                Rcpp::Named("col")   = (int)c + 1);
  }
  // 只有 has_implicit 为假时才可能走到这里（极值必定是某个存储值）。
  return Rcpp::List::create(Rcpp::Named("value") = NA_REAL,
                            Rcpp::Named("row")   = NA_INTEGER,
                            Rcpp::Named("col")   = NA_INTEGER);
}

// O 范数（最大列绝对值和）与 I 范数（最大行绝对值和），一次扫描同时算出来。
//
// 为什么不用 base 的 tabulate()：它没有 weights 参数（我记错了，第一次就是
// 这么写错的）。能用的是 tapply()/rowsum()，两者都要先按分组排序，是
// O(nnz log nnz)，而且 tapply 会为每个分组建一个 list 元素——nnz 上百万时
// 光是分配就够呛。这里是一遍遍历、两个累加器，O(nnz)，只多一个长度 rows 的
// 向量。
//
// NA 的处理跟 base::norm() 一致：只要有一个缺失值，两种范数都是 NA。
// 分开返回两个数是因为它们共用同一遍扫描，分开算等于白扫一次。
// [[Rcpp::export]]
Rcpp::NumericVector reigen_sparse_norm_oi(SEXP a) {
  const SpMat& m = sp_ptr(a, "a")->m;
  if (m.nonZeros() == 0) return Rcpp::NumericVector::create(0.0, 0.0);
  std::vector<double> rowabs((size_t)m.rows(), 0.0);
  double maxcol = 0.0;
  for (Eigen::Index c = 0; c < m.outerSize(); ++c) {
    double cs = 0.0;
    for (SpMat::InnerIterator it(m, c); it; ++it) {
      const double raw = it.value();
      if (R_IsNA(raw) || R_IsNaN(raw))
        return Rcpp::NumericVector::create(NA_REAL, NA_REAL);
      const double v = std::abs(raw);
      cs += v;
      rowabs[(size_t)it.row()] += v;
    }
    if (cs > maxcol) maxcol = cs;
  }
  double maxrow = 0.0;
  for (size_t i = 0; i < rowabs.size(); ++i)
    if (rowabs[i] > maxrow) maxrow = rowabs[i];
  return Rcpp::NumericVector::create(maxcol, maxrow);
}

// 谱范数（最大奇异值），幂迭代。
//
// 为什么不能像 F 范数那样一次扫描算出来：σ_max 是 AᵀA 的最大特征值的平方根，
// 而稀疏矩阵没有特征分解（Eigen 的 SelfAdjointEigenSolver 只吃稠密）。要么
// 稠密化——那正是稀疏矩阵想避免的；要么迭代。
//
// 迭代的是 AᵀA：z = Aᵀ(A v)，v ← z/‖z‖。收敛后 ‖z‖ 就是 λ_max，
// σ_max = √λ_max。这里用 ‖z‖ 而不是 vᵀ(AᵀA)v 当估计值：v 有误差 ε 时
// λ 的误差是 O(ε²)（Rayleigh 商的性质），所以按 ‖z‖ 判收敛，实际精度比
// 向量的收敛速度还好一倍。
//
// 起点用固定种子的伪随机向量，不用全 1。全 1 起点有个真实的坑：如果主奇异
// 向量与全 1 向量正交（例如矩阵每列和为 0），第一轮就投影成 0，迭代永远停在
// 0 上——一个 0 的答案，看起来完全合法。
//
// 返回值：收敛则返回 σ_max；maxiter 内没收敛返回 NA，由 R 侧报错说明。
// 不用 warning：R 侧的 norm() 契约是"要么给数，要么报错"，悄悄给一个
// 没收敛的数是最糟的第三种结果。
// [[Rcpp::export]]
double reigen_sparse_norm2(SEXP a, double tol, int maxiter) {
  const SpMat& m = sp_ptr(a, "a")->m;
  const Eigen::Index n = m.cols();
  if (n == 0 || m.rows() == 0) return 0.0;
  if (m.nonZeros() == 0) return 0.0;

  // 确定性起点：同一个矩阵永远得到同一个答案，测试才能断言到 1e-10。
  Eigen::VectorXd v(n);
  unsigned long long s = 88172645463325252ULL;
  for (Eigen::Index i = 0; i < n; ++i) {
    s ^= s << 13; s ^= s >> 7; s ^= s << 17;          // xorshift64
    v[i] = ((double)(s >> 11) / 9007199254740992.0) - 0.5;  // [-0.5, 0.5)
  }
  double nv = v.norm();
  if (nv == 0.0) { v.setOnes(); nv = std::sqrt((double)n); }
  v /= nv;

  double lambda = 0.0;
  for (int it = 0; it < maxiter; ++it) {
    const Eigen::VectorXd u = m * v;          // A v
    const Eigen::VectorXd z = m.transpose() * u;  // AᵀA v
    const double nz = z.norm();
    if (nz == 0.0) return 0.0;                // 整列为零
    v = z / nz;
    // 相对判据。nz 就是 λ_max 的当前估计，用相对误差才和矩阵的尺度无关。
    if (std::abs(nz - lambda) <= tol * std::max(1.0, std::abs(nz)))
      return std::sqrt(std::max(0.0, nz));
    lambda = nz;
  }
  return NA_REAL;
}

// 主对角线。Eigen 的稀疏矩阵没有现成的 diag()，coeff() 是 O(log nnz)，
// 取一条对角线正好合适。隐含零自然读成 0。
// [[Rcpp::export]]
Rcpp::NumericVector reigen_sparse_diag(SEXP a) {
  const SpMat& m = sp_ptr(a, "a")->m;
  const Eigen::Index k = std::min(m.rows(), m.cols());
  Rcpp::NumericVector out(k);
  for (Eigen::Index i = 0; i < k; ++i) out[i] = m.coeff(i, i);
  return out;
}


// 逐元素幂，只对**存下来的**元素做变换。前提是 p > 0：隐含的零元素
// 0^p 在 p>0 时仍然是 0，所以非零结构完全不变，结果还是稀疏的。
// p <= 0 时 0^p 是 Inf（p<0）或 1（p==0），整个矩阵几乎处处非零，
// 再硬撑稀疏只会又慢又占内存——那种情况由 R 侧退回稠密，不要在这里做。
// [[Rcpp::export]]
SEXP reigen_sparse_pow(SEXP a, double p) {
  if (!(p > 0.0))
    Rcpp::stop("reigen_sparse_pow: exponent must be positive, got %g", p);
  SparseObj* pa = sp_ptr(a, "a");
  SpMat out = pa->m;                       // 拷贝：绝不改动入参
  for (Eigen::Index c = 0; c < out.outerSize(); ++c)
    for (SpMat::InnerIterator it(out, c); it; ++it)
      it.valueRef() = std::pow(it.value(), p);
  out.makeCompressed();
  return wrap_sparse(std::move(out));
}

// 逐元素乘（不是矩阵乘）。Eigen 的 cwiseProduct 只在两侧都是 SparseMatrix
// 时返回稀疏结果；维度不匹配在 NDEBUG 下同样不报错，所以自己判。
// [[Rcpp::export]]
SEXP reigen_sparse_cwise(SEXP a, SEXP b) {
  SparseObj* pa = sp_ptr(a, "a");
  SparseObj* pb = sp_ptr(b, "b");
  need_same_shape(pa->m, pb->m, "*");
  SpMat out = pa->m.cwiseProduct(pb->m);
  out.makeCompressed();
  return wrap_sparse(std::move(out));
}

// ============================================================ 求解
//
// method: 0 SparseLU                       通用方阵，直接法
//         1 SimplicialLLT                  对称正定
//         2 SimplicialLDLT                 对称不定
//         3 SparseQR                       可长方形（最小二乘）
//         4 ConjugateGradient              对称正定，迭代
//         5 BiCGSTAB                       一般方阵（非对称），迭代
//         6 LeastSquaresConjugateGradient  长方形最小二乘，迭代
//         7 GMRES                          一般方阵，迭代，重启式
//         8 MINRES                         对称（可不定），迭代
//         9 DGMRES                         一般方阵，GMRES + 放缩
//        10 IDRS                           IDR(s)，一般方阵，迭代
//
// 全部都要查 info()。Eigen 的分解**不会抛异常**：失败只是把 info() 设成
// NumericalIssue，解仍然照常返回——返回的是另一个方程组的解。稠密侧
// eigen_llt() 就是栽在这上面（见 reigen_decomp.cpp 的注释）。
//
// 迭代法还多一层坑：没收敛时 info() 是 NoConvergence，但解照样返回，
// 而且**看起来很正常**。所以迭代法必须把 iterations()/error() 一起交出去，
// 让用户能看到它到底走了几步、残差多大。
//
// precond: 0 IncompleteLUT（默认）, 1 对角。只影响 4、5 以及 7~10。
//
// 默认给 ILUT 是量出来的，不是猜的。48x48 的 2 维对流扩散（次对角 -1.1）：
// 对角预条件要 1352 步，ILUT 要 99 步；对流再强一点（-1.2），对角预条件
// 干脆不收敛，ILUT 1925 步仍能到。ILUT 做的是**无主元**的不完全 LU，遇到
// 零主元会产出 NaN（不是崩溃），所以它并非永远可用——那种情况下会报
// 不收敛，让用户退回对角预条件或直接法。

template <class Solver>
static Mat run_solver(Solver& s, const SpMat& A, const Eigen::MatrixXd& B,
                      const char* name) {
  s.compute(A);
  if (s.info() != Eigen::Success)
    Rcpp::stop("%s failed to factorise the matrix (Eigen::info() = %d). "
               "The matrix may be singular, rank deficient, or not positive "
               "definite for this method; try method = \"lu\" or \"qr\".",
               name, (int)s.info());
  Eigen::MatrixXd X = s.solve(B);
  if (s.info() != Eigen::Success)
    Rcpp::stop("%s failed during solve (Eigen::info() = %d)", name, (int)s.info());
  return make(X);
}

// 直接法的返回：没有迭代次数可报，两个字段置 NA。
//
// wrap_mat() 返回的是裸 SEXP，构造 List 的过程中可能触发 GC，所以要先
// PROTECT 住——这是 wrap_obj() 里 gc_tick() 会调 R_gc() 之后新增的规矩。
static Rcpp::List plain_result(const Mat& m) {
  SEXP x = wrap_mat(m);
  PROTECT(x);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("x") = x,
                                      Rcpp::Named("iterations") = NA_REAL,
                                      Rcpp::Named("error") = NA_REAL);
  UNPROTECT(1);
  return out;
}

template <class Solver>
static Rcpp::List run_iterative(Solver& s, const SpMat& A,
                                const Eigen::MatrixXd& B,
                                double tol, int maxiter, const char* name) {
  s.setTolerance(tol);
  s.setMaxIterations(maxiter);
  s.compute(A);
  if (s.info() != Eigen::Success)
    Rcpp::stop("%s failed to initialise (Eigen::info() = %d)", name, (int)s.info());
  Eigen::MatrixXd X = s.solve(B);
  const int iters = (int)s.iterations();
  const double err = (double)s.error();
  if (s.info() != Eigen::Success)
    Rcpp::stop("%s did not converge in %d iterations (relative residual %g). "
               "Raise 'maxiter', loosen 'tol', or use a direct method.",
               name, iters, err);
  SEXP x = wrap_mat(make(X));
  PROTECT(x);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("x") = x,
                                      Rcpp::Named("iterations") = iters,
                                      Rcpp::Named("error") = err);
  UNPROTECT(1);
  return out;
}

// [[Rcpp::export]]
Rcpp::List reigen_sparse_solve(SEXP a, SEXP b, int method, double tol,
                               int maxiter, int precond) {
  SparseObj* pa = sp_ptr(a, "a");
  Mat B = from_sexp_finite(b, "b");
  if (B.cplx)
    Rcpp::stop("complex sparse solve is not supported yet");
  const SpMat& A = pa->m;

  if (method < 0 || method > 10)
    Rcpp::stop("unknown sparse solver code %d", method);

  need_finite(A, "a");

  // SparseQR 和 LeastSquaresConjugateGradient 接受长方形；其余要求方阵。
  const bool needs_square = (method != 3 && method != 6);
  if (needs_square && A.rows() != A.cols())
    Rcpp::stop("this solver requires a square matrix, got %lldx%lld; "
               "use method = \"qr\" or \"lsqr\" for least squares",
               (long long)A.rows(), (long long)A.cols());
  if (A.rows() != rows(B))
    Rcpp::stop("dimensions do not match: A is %lldx%lld but b has %lld rows",
               (long long)A.rows(), (long long)A.cols(), (long long)rows(B));

  if (!(tol > 0.0)) tol = 1e-12;
  if (maxiter <= 0) maxiter = (int)std::max<Eigen::Index>(1000, A.rows() * 10);

  const Eigen::MatrixXd& Bm = B.r;
  const bool use_ilut = (precond == 0);   // 0 = ilut（默认）, 1 = diagonal

  switch (method) {
    case 0: {
      Eigen::SparseLU<SpMat> s;
      return plain_result(run_solver(s, A, Bm, "SparseLU"));
    }
    case 1: {
      Eigen::SimplicialLLT<SpMat> s;
      return plain_result(run_solver(s, A, Bm, "SimplicialLLT"));
    }
    case 2: {
      Eigen::SimplicialLDLT<SpMat> s;
      return plain_result(run_solver(s, A, Bm, "SimplicialLDLT"));
    }
    case 3: {
      Eigen::SparseQR<SpMat, Eigen::COLAMDOrdering<int>> s;
      return plain_result(run_solver(s, A, Bm, "SparseQR"));
    }
    case 4: {
      // CG 只接受**对称正定**的预条件子。IncompleteLUT 不是对称的，配上去
      // 收敛性就没有任何保证——实测在 2 维 Poisson 上直接不收敛。所以这里
      // 只用对角预条件，precond 参数对 CG 无效（R 侧传 ilut 会先报错）。
      Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> s;
      return run_iterative(s, A, Bm, tol, maxiter, "ConjugateGradient");
    }
    case 5: {
      // BiCGSTAB 处理非对称系统——这是 CG 做不到的那一类，而实践中大多数
      // 稀疏系统恰恰是非对称的。
      if (use_ilut) {
        Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double>> s;
        return run_iterative(s, A, Bm, tol, maxiter, "BiCGSTAB");
      }
      Eigen::BiCGSTAB<SpMat> s;
      return run_iterative(s, A, Bm, tol, maxiter, "BiCGSTAB");
    }
    case 6: {
      Eigen::LeastSquaresConjugateGradient<SpMat> s;
      return run_iterative(s, A, Bm, tol, maxiter,
                           "LeastSquaresConjugateGradient");
    }
    case 7: {
      // GMRES 的残差单调下降，能处理 CG 和 BiCGSTAB 都不好办的病态非对称
      // 系统；代价是每步都要存一整个 Krylov 基，内存随迭代数线性增长。
      // 这里不设显式的重启长度，用 Eigen 的默认行为。
      if (use_ilut) {
        Eigen::GMRES<SpMat, Eigen::IncompleteLUT<double>> s;
        return run_iterative(s, A, Bm, tol, maxiter, "GMRES");
      }
      Eigen::GMRES<SpMat> s;
      return run_iterative(s, A, Bm, tol, maxiter, "GMRES");
    }
    case 8: {
      // MINRES 只要求对称（不要求正定），这是它相对 CG 的全部价值。
      // 预条件子同样必须对称，所以和 CG 一样只给对角。
      // 默认的 UpLo = Lower：只读下三角，上三角被忽略——传进来的矩阵若
      // 不对称，它不会报错，而是把你当成了下三角那个对称矩阵。
      Eigen::MINRES<SpMat, Eigen::Lower,
                    Eigen::DiagonalPreconditioner<double>> s;
      return run_iterative(s, A, Bm, tol, maxiter, "MINRES");
    }
    case 9: {
      // DGMRES 在 GMRES 之外多做一步对角放缩，对行尺度相差悬殊的系统
      // 明显更稳（例如同时含单位系数和 1e6 系数的方程组）。
      if (use_ilut) {
        Eigen::DGMRES<SpMat, Eigen::IncompleteLUT<double>> s;
        return run_iterative(s, A, Bm, tol, maxiter, "DGMRES");
      }
      Eigen::DGMRES<SpMat> s;
      return run_iterative(s, A, Bm, tol, maxiter, "DGMRES");
    }
    default: {
      // IDR(s)，s = 4（Eigen 默认值 1 在小问题上省内存，但在难收敛的
      // 不定系统上要的步数明显更多；4 是论文推荐的折中）。
      // IDR 没有给 Preconditioner_ 默认值，必须写全两个模板参数。
      Eigen::IDRS<SpMat, Eigen::DiagonalPreconditioner<double>> s;
      s.setS(4);
      return run_iterative(s, A, Bm, tol, maxiter, "IDRS");
    }
  }
}

// 求解器给出的残差，方便 R 侧判断 CG 是否可信。
// [[Rcpp::export]]
double reigen_sparse_cg_error(SEXP a, SEXP b) {
  SparseObj* pa = sp_ptr(a, "a");
  Mat B = from_sexp_finite(b, "b");
  if (B.cplx) Rcpp::stop("complex sparse solve is not supported yet");
  Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> s;
  s.setTolerance(1e-12);
  s.setMaxIterations((int)std::max<Eigen::Index>(1000, pa->m.rows() * 10));
  s.compute(pa->m);
  s.solve(B.r);
  return (double)s.error();
}
