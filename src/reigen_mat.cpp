// reigen_mat.cpp —— 构造、转换、下标、拼接、生成器
//
// 这一层只做「形状与存取」，不做数值。所有数值运算在 reigen_ops.cpp，
// 分解在 reigen_decomp.cpp。

#include "reigen.h"

using namespace reigen;

// ---------------------------------------------------------------- 构造

// [[Rcpp::export]]
SEXP reigen_mat_new(SEXP x) {
  return wrap_mat(from_sexp(x));
}

// 列向量：n x 1。class 多一档 eigenVec，S3 派发仍能找到 eigenMat 的方法。
// [[Rcpp::export]]
SEXP reigen_vec_new(SEXP x) {
  Mat m = from_sexp(x);
  Mat out;
  out.cplx = m.cplx;
  if (m.cplx) out.c = Eigen::Map<Eigen::MatrixXcd>(m.c.data(), (Eigen::Index)nelem(m), 1);
  else        out.r = Eigen::Map<Eigen::MatrixXd>(m.r.data(), (Eigen::Index)nelem(m), 1);
  Rcpp::XPtr<Mat> xp(new Mat(std::move(out)));
  set_class_vec((SEXP)xp);
  return xp;
}

// [[Rcpp::export]]
SEXP reigen_mat_to_r(SEXP x) {
  return to_sexp(*mat_ptr(x));
}

// [[Rcpp::export]]
SEXP reigen_mat_dim(SEXP x) {
  const Mat& m = *mat_ptr(x);
  return Rcpp::IntegerVector::create((int)rows(m), (int)cols(m));
}

// [[Rcpp::export]]
bool reigen_mat_is_complex(SEXP x) {
  return mat_ptr(x)->cplx;
}

// [[Rcpp::export]]
SEXP reigen_mat_reshape(SEXP x, int nr, int nc) {
  const Mat& m = *mat_ptr(x);
  const Eigen::Index n = nelem(m);
  if ((Eigen::Index)nr * nc != n)
    Rcpp::stop("cannot reshape %lld elements into %dx%d", (long long)n, nr, nc);
  Mat out;
  out.cplx = m.cplx;
  if (m.cplx)
    out.c = Eigen::Map<const Eigen::MatrixXcd>(m.c.data(), nr, nc);
  else
    out.r = Eigen::Map<const Eigen::MatrixXd>(m.r.data(), nr, nc);
  return wrap_mat(std::move(out));
}

// ---------------------------------------------------------------- 元素存取

// [[Rcpp::export]]
SEXP reigen_mat_get(SEXP x, int i, int j) {
  const Mat& m = *mat_ptr(x);
  Eigen::Index r = idx1(i, rows(m), "row");
  Eigen::Index c = idx1(j, cols(m), "column");
  if (m.cplx) return wrap_cplx(m.c(r, c));
  return Rcpp::wrap(m.r(r, c));
}

// 返回新对象，不原地改——R 的语义是值语义，调用方用 A[i,j] <- v 时
// R 的 [<- 本身就会把结果重新赋给 A。
// [[Rcpp::export]]
SEXP reigen_mat_set(SEXP x, int i, int j, SEXP v) {
  Mat out = *mat_ptr(x);
  Eigen::Index r = idx1(i, rows(out), "row");
  Eigen::Index c = idx1(j, cols(out), "column");
  Mat val = from_sexp(v, "value");
  if (val.cplx || out.cplx) {
    if (!out.cplx) { out.c = out.r.template cast<cplx_t>(); out.cplx = true; out.r.resize(0, 0); }
    out.c(r, c) = as_cplx(val).data()[0];
  } else {
    out.r(r, c) = val.r.data()[0];
  }
  return wrap_mat(std::move(out));
}

// [[Rcpp::export]]
SEXP reigen_mat_row(SEXP x, int i) {
  const Mat& m = *mat_ptr(x);
  Eigen::Index r = idx1(i, rows(m), "row");
  if (m.cplx) return wrap_mat(make(Eigen::MatrixXcd(m.c.row(r).transpose())));
  return wrap_mat(make(Eigen::MatrixXd(m.r.row(r).transpose())));
}

// [[Rcpp::export]]
SEXP reigen_mat_col(SEXP x, int j) {
  const Mat& m = *mat_ptr(x);
  Eigen::Index c = idx1(j, cols(m), "column");
  if (m.cplx) return wrap_mat(make(Eigen::MatrixXcd(m.c.col(c))));
  return wrap_mat(make(Eigen::MatrixXd(m.r.col(c))));
}

// r0..r1 / c0..c1 均为 1-based 闭区间；传 0 表示「到端点」。
// [[Rcpp::export]]
SEXP reigen_mat_block(SEXP x, int r0, int r1, int c0, int c1) {
  const Mat& m = *mat_ptr(x);
  const Eigen::Index nr = rows(m), nc = cols(m);
  Eigen::Index a = (r0 == 0) ? 0 : idx1(r0, nr, "row");
  Eigen::Index b = (r1 == 0) ? nr - 1 : idx1(r1, nr, "row");
  Eigen::Index c = (c0 == 0) ? 0 : idx1(c0, nc, "column");
  Eigen::Index d = (c1 == 0) ? nc - 1 : idx1(c1, nc, "column");
  if (b < a || d < c) Rcpp::stop("empty or reversed block range");
  if (m.cplx) return wrap_mat(make(Eigen::MatrixXcd(m.c.block(a, c, b - a + 1, d - c + 1))));
  return wrap_mat(make(Eigen::MatrixXd(m.r.block(a, c, b - a + 1, d - c + 1))));
}

// [[Rcpp::export]]
SEXP reigen_mat_block_set(SEXP x, int r0, int r1, int c0, int c1, SEXP v) {
  Mat out = *mat_ptr(x);
  const Eigen::Index nr = rows(out), nc = cols(out);
  Eigen::Index a = (r0 == 0) ? 0 : idx1(r0, nr, "row");
  Eigen::Index b = (r1 == 0) ? nr - 1 : idx1(r1, nr, "row");
  Eigen::Index c = (c0 == 0) ? 0 : idx1(c0, nc, "column");
  Eigen::Index d = (c1 == 0) ? nc - 1 : idx1(c1, nc, "column");
  Mat val = from_sexp(v, "value");
  const Eigen::Index bh = b - a + 1, bw = d - c + 1;
  if (nelem(val) == 1) {
    if (val.cplx || out.cplx) {
      if (!out.cplx) { out.c = out.r.template cast<cplx_t>(); out.cplx = true; out.r.resize(0, 0); }
      out.c.block(a, c, bh, bw).setConstant(as_cplx(val).data()[0]);
    } else {
      out.r.block(a, c, bh, bw).setConstant(val.r.data()[0]);
    }
    return wrap_mat(std::move(out));
  }
  if (rows(val) != bh || cols(val) != bw)
    Rcpp::stop("value is %lldx%lld, block is %lldx%lld",
               (long long)rows(val), (long long)cols(val), (long long)bh, (long long)bw);
  if (val.cplx || out.cplx) {
    if (!out.cplx) { out.c = out.r.template cast<cplx_t>(); out.cplx = true; out.r.resize(0, 0); }
    out.c.block(a, c, bh, bw) = as_cplx(val);
  } else {
    out.r.block(a, c, bh, bw) = val.r;
  }
  return wrap_mat(std::move(out));
}

// ---------------------------------------------------------------- 对角线

// [[Rcpp::export]]
SEXP reigen_mat_diag(SEXP x) {
  const Mat& m = *mat_ptr(x);
  if (m.cplx) return wrap_mat(make(Eigen::MatrixXcd(m.c.diagonal())));
  return wrap_mat(make(Eigen::MatrixXd(m.r.diagonal())));
}

// [[Rcpp::export]]
SEXP reigen_mat_diag_set(SEXP x, SEXP v) {
  Mat out = *mat_ptr(x);
  Mat val = from_sexp(v, "value");
  const Eigen::Index n = std::min(rows(out), cols(out));
  if (val.cplx || out.cplx) {
    if (!out.cplx) { out.c = out.r.template cast<cplx_t>(); out.cplx = true; out.r.resize(0, 0); }
    Eigen::MatrixXcd vv = as_cplx(val);
    if (vv.size() == 1) {
      const cplx_t s = vv.data()[0];
      for (Eigen::Index i = 0; i < n; ++i) out.c(i, i) = s;
    } else {
      if (vv.size() != n)
        Rcpp::stop("diagonal length %lld does not match %lld", (long long)vv.size(), (long long)n);
      for (Eigen::Index i = 0; i < n; ++i) out.c(i, i) = vv.data()[i];
    }
  } else {
    if (nelem(val) == 1) {
      const double s = val.r.data()[0];
      for (Eigen::Index i = 0; i < n; ++i) out.r(i, i) = s;
    } else {
      if (nelem(val) != n)
        Rcpp::stop("diagonal length %lld does not match %lld", (long long)nelem(val), (long long)n);
      for (Eigen::Index i = 0; i < n; ++i) out.r(i, i) = val.r.data()[i];
    }
  }
  return wrap_mat(std::move(out));
}

// [[Rcpp::export]]
SEXP reigen_mat_trace(SEXP x) {
  const Mat& m = *mat_ptr(x);
  if (m.cplx) return wrap_cplx(m.c.trace());
  return Rcpp::wrap(m.r.trace());
}

// ---------------------------------------------------------------- 拼接

// side: 1 = rbind（纵向），2 = cbind（横向）
// [[Rcpp::export]]
SEXP reigen_mat_bind(SEXP a, SEXP b, int side) {
  Mat A = from_sexp(a, "a");
  Mat B = from_sexp(b, "b");
  const bool cz = A.cplx || B.cplx;
  if (cz) {
    Eigen::MatrixXcd X = as_cplx(A), Y = as_cplx(B);
    Eigen::MatrixXcd out;
    if (side == 1) {
      if (X.cols() != Y.cols()) Rcpp::stop("rbind: %lld vs %lld columns", (long long)X.cols(), (long long)Y.cols());
      out.resize(X.rows() + Y.rows(), X.cols());
      out << X, Y;
    } else {
      if (X.rows() != Y.rows()) Rcpp::stop("cbind: %lld vs %lld rows", (long long)X.rows(), (long long)Y.rows());
      out.resize(X.rows(), X.cols() + Y.cols());
      out << X, Y;
    }
    return wrap_mat(make(out));
  }
  Eigen::MatrixXd X = A.r, Y = B.r;
  Eigen::MatrixXd out;
  if (side == 1) {
    if (X.cols() != Y.cols()) Rcpp::stop("rbind: %lld vs %lld columns", (long long)X.cols(), (long long)Y.cols());
    out.resize(X.rows() + Y.rows(), X.cols());
    out << X, Y;
  } else {
    if (X.rows() != Y.rows()) Rcpp::stop("cbind: %lld vs %lld rows", (long long)X.rows(), (long long)Y.rows());
    out.resize(X.rows(), X.cols() + Y.cols());
    out << X, Y;
  }
  return wrap_mat(make(out));
}

// [[Rcpp::export]]
SEXP reigen_mat_rep(SEXP x, int nr, int nc) {
  const Mat& m = *mat_ptr(x);
  if (nr < 1 || nc < 1) Rcpp::stop("nr and nc must be >= 1");
  if (m.cplx) return wrap_mat(make(Eigen::MatrixXcd(m.c.replicate(nr, nc))));
  return wrap_mat(make(Eigen::MatrixXd(m.r.replicate(nr, nc))));
}

// t()：只转置，不取共轭。
// [[Rcpp::export]]
SEXP reigen_mat_transpose(SEXP x) {
  const Mat& m = *mat_ptr(x);
  if (m.cplx) return wrap_mat(make(Eigen::MatrixXcd(m.c.transpose())));
  return wrap_mat(make(Eigen::MatrixXd(m.r.transpose())));
}

// 共轭转置（Eigen 的 adjoint）。实矩阵退化为普通转置。
// [[Rcpp::export]]
SEXP reigen_mat_adjoint(SEXP x) {
  const Mat& m = *mat_ptr(x);
  if (m.cplx) return wrap_mat(make(Eigen::MatrixXcd(m.c.adjoint())));
  return wrap_mat(make(Eigen::MatrixXd(m.r.transpose())));
}

// ---------------------------------------------------------------- 生成器

// [[Rcpp::export]]
SEXP reigen_mat_identity(int n) {
  if (n < 0) Rcpp::stop("n must be >= 0");
  return wrap_mat(make(Eigen::MatrixXd(Eigen::MatrixXd::Identity(n, n))));
}

// [[Rcpp::export]]
SEXP reigen_mat_zeros(int nr, int nc) {
  if (nr < 0 || nc < 0) Rcpp::stop("nr and nc must be >= 0");
  return wrap_mat(make(Eigen::MatrixXd(Eigen::MatrixXd::Zero(nr, nc))));
}

// [[Rcpp::export]]
SEXP reigen_mat_ones(int nr, int nc) {
  if (nr < 0 || nc < 0) Rcpp::stop("nr and nc must be >= 0");
  return wrap_mat(make(Eigen::MatrixXd(Eigen::MatrixXd::Ones(nr, nc))));
}

// [[Rcpp::export]]
SEXP reigen_mat_constant(int nr, int nc, double v) {
  if (nr < 0 || nc < 0) Rcpp::stop("nr and nc must be >= 0");
  return wrap_mat(make(Eigen::MatrixXd(Eigen::MatrixXd::Constant(nr, nc, v))));
}

// 用 R 的 RNG，保证 set.seed() 可复现——Eigen 自带的 Random 走 std::rand，
// 不受 set.seed 控制，在可复现的脚本里是个陷阱。
// [[Rcpp::export]]
SEXP reigen_mat_random(int nr, int nc) {
  if (nr < 0 || nc < 0) Rcpp::stop("nr and nc must be >= 0");
  const Eigen::Index n = (Eigen::Index)nr * nc;
  Eigen::MatrixXd out(nr, nc);
  Rcpp::NumericVector u = Rcpp::runif((int)n);
  for (Eigen::Index i = 0; i < n; ++i) out.data()[i] = u[i];
  return wrap_mat(make(out));
}

// 与 seq(lo, hi, length.out = n) 一致：端点是闭的。
// [[Rcpp::export]]
SEXP reigen_mat_linspace(double lo, double hi, int n) {
  if (n < 1) Rcpp::stop("n must be >= 1");
  Eigen::MatrixXd out(n, 1);
  if (n == 1) { out(0, 0) = lo; return wrap_mat(make(out)); }
  const double step = (hi - lo) / (double)(n - 1);
  for (int i = 0; i < n; ++i) out(i, 0) = lo + step * i;
  out(n - 1, 0) = hi;  // 抹掉累积误差
  return wrap_mat(make(out));
}
