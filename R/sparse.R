# 稀疏矩阵
#
# eigenSparse 包的是 Eigen::SparseMatrix<double>，列主序，与 Matrix 包的
# dgCMatrix 同构。索引在 R 侧一律 1-based，与 Matrix::sparseMatrix() 一致。
#
# 只做实数：复数稀疏会明确报错，而不是静默取实部。
#
# 为什么要它：稠密分解是 O(n^3)，一个 10^5 阶的带状矩阵根本存不下。稀疏求解
# 器（SparseLU / SimplicialLLT / CG）只对结构做工作，这才是 Eigen 相对 base R
# 的核心优势——base R 里没有稀疏线性代数。

# ---------------------------------------------------------------- 构造

#' Sparse matrices
#'
#' `eigen_sparse()` builds a real sparse matrix in compressed column form. It
#' has two calling forms:
#'
#' * `eigen_sparse(i, j, x)` from triplets, where `i` and `j` are **1-based**
#'   row and column indices, exactly as [Matrix::sparseMatrix()] takes them.
#'   Duplicate `(i, j)` pairs are **summed**.
#' * `eigen_sparse(x)` from a dense matrix, `eigenMat`, or (if the Matrix
#'   package is installed) a `sparseMatrix`. Only exact zeros are dropped.
#'
#' Complex sparse matrices are not supported: complex input stops with an error
#' rather than silently discarding the imaginary part.
#'
#' `as.numeric()` does **not** work on an `eigenSparse`: it is a primitive and
#' does not dispatch on S3 classes, so it sees only the external pointer. Use
#' `as.vector()` or `as.numeric(as.matrix(x))`. The same caveat applies to
#' `as.numeric()` on `eigenMat`.
#'
#' @section Arithmetic:
#' Operators keep a result sparse only when the result *is* sparse. That is a
#' property of the operation, not a preference:
#'
#' * `S + S`, `S - S`, `S * S` (elementwise), `S * scalar`, `S / scalar` and
#'   `S ^ p` for `p > 0` stay sparse — an implicit zero stays zero under all of
#'   them.
#' * `S ^ p` for `p <= 0` returns a dense `eigenMat`, because `0^p` is `Inf`
#'   (`p < 0`) or `1` (`p == 0`), so the result has no zeros left to exploit.
#' * `S / 0` follows base R and returns a dense result with `Inf` where the
#'   entry was non-zero and `NaN` where it was zero. It does **not** raise an
#'   error; base `matrix / 0` does not either.
#' * Comparisons (`==`, `<`, `>`, `<=`, `>=`, `!=`) and `!` return a plain
#'   logical matrix, as they do for `eigenMat`.
#' * Mixing a sparse and a dense operand returns dense.
#'
#' @section Reductions and math:
#' `sum()`, `prod()`, `max()`, `min()`, `range()`, `any()`, `all()`, `mean()`,
#' `is.na()`, `anyNA()` and `as.data.frame()` all work on an `eigenSparse` and
#' agree with base R. The reductions are computed from the stored entries in
#' `O(nnz)`; they do not densify. **Implicit zeros take part**: for an all
#' positive sparse matrix, `min()` is `0` and `prod()` is `0`, exactly as for
#' the dense matrix.
#'
#' `NA` and `NaN` stay distinct, as in base R: `sum()` of a matrix containing
#' `NA` is `NA`, one containing `NaN` is `NaN`. `na.rm = TRUE` drops both.
#'
#' The `Math` group keeps a result sparse when the operation maps zero to zero
#' (`abs`, `sqrt`, `sign`, `floor`, `round`, `sin`, `tan`, `log1p`, ...). It
#' returns a dense `eigenMat` for operations that do not (`exp`, `log`, `cos`,
#' `gamma`, `cumsum`, ...), since those either give a non-zero value at zero or
#' depend on the whole column.
#'
#' @param i,j Integer vectors of 1-based row and column indices. `i` may also
#'   be a dense matrix or an existing sparse matrix, in which case `j` and `x`
#'   must be missing.
#' @param x Numeric vector of values. Non-finite values are stored as given;
#'   the solvers reject them.
#' @param nrow,ncol Dimensions. Default to `max(i)` and `max(j)`; supply them
#'   explicitly to keep empty trailing rows or columns.
#' @return An `eigenSparse` object.
#' @name eigen_sparse
#' @examples
#' \dontrun{
#' S <- eigen_sparse(c(1, 2, 3), c(1, 2, 3), c(4, 5, 6), nrow = 3, ncol = 3)
#' dim(S); eigen_sparse_nnz(S)
#' as.matrix(S)
#' }
NULL

#' @rdname eigen_sparse
#' @export
eigen_sparse <- function(i, j = NULL, x = NULL, nrow = NULL, ncol = NULL) {
  # 只有一个参数：从稠密矩阵 / dgCMatrix / eigenMat 转换。
  if (is.null(j) && is.null(x)) {
    if (inherits(i, "sparseMatrix") || inherits(i, "Matrix")) return(.sparse_from_matrix(i))
    return(reigen_sparse_from_dense(.as_eigen(i)))
  }
  if (is.null(j) || is.null(x))
    stop("give either a single matrix-like object, or all three of i, j and x",
         call. = FALSE)

  # dgCMatrix 的三元组也走这条路（见 .sparse_from_matrix）。
  i <- as.integer(i)
  j <- as.integer(j)
  x <- as.numeric(x)

  # nrow/ncol 缺省时取下标的最大值。空输入没有下标可推，必须显式给。
  if (is.null(nrow)) nrow <- if (length(i)) max(i) else 0L
  if (is.null(ncol)) ncol <- if (length(j)) max(j) else 0L

  reigen_sparse_triplet(i, j, x, as.integer(nrow), as.integer(ncol))
}

# dgCMatrix -> 三元组 -> eigenSparse。summary() 给的是 1-based i/j/x，
# 正好是 reigen_sparse_triplet() 要的格式。
#
# as() 是 methods 包的 S4 泛型，不是 base 的 as.* 系列，所以必须写全
# methods::as —— 只写 as(m, "TsparseMatrix") 时 R CMD check 会报
# "no visible global function definition for 'as'"。
.sparse_from_matrix <- function(m) {
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("converting a Matrix object requires the Matrix package", call. = FALSE)
  s <- methods::as(m, "TsparseMatrix")
  sm <- Matrix::summary(s)
  if (!is.numeric(sm$x))
    stop("only numeric sparse matrices are supported, got ",
         class(sm$x)[1], call. = FALSE)
  reigen_sparse_triplet(as.integer(sm$i), as.integer(sm$j), as.numeric(sm$x),
                        as.integer(nrow(s)), as.integer(ncol(s)))
}

#' @export
dim.eigenSparse <- function(x) reigen_sparse_dim(x)

# XPtr 没有 "dim" 属性，base 的 length() 会返回 1。矩阵的 length 是元素总数，
# 这里要自己补上，否则 length(S) 与 length(as.matrix(S)) 不一致。
# 用 double 而不是 integer：50000 x 50000 的元素总数超过 .Machine$integer.max，
# 与 base R 对长向量的处理保持一致。
#' @export
length.eigenSparse <- function(x) prod(as.double(reigen_sparse_dim(x)))

#' @export
as.matrix.eigenSparse <- function(x, ...) as.matrix(reigen_sparse_as_matrix(x))

#' @export
as.vector.eigenSparse <- function(x, mode = "any") as.vector(as.matrix(x), mode = mode)

# 这里**不要**写 as.numeric.eigenSparse。as.numeric 是 primitive，不做 S3 派发，
# 写了也不会被调用（本文件第一版就写了，结果 as.numeric(S) 仍然报
# "cannot coerce type 'externalptr' to vector of type 'double'"）。
# 要用 as.vector(S) 或 as.numeric(as.matrix(S))。这一点在 ?eigen_sparse 里写明了。

#' @export
t.eigenSparse <- function(x) reigen_sparse_transpose(x)

# crossprod / tcrossprod / aperm / is.finite 这一组在稀疏侧也必须各写一份：
# 稀疏对象的 class 向量里没有 "eigenMat"，派发根本到不了稠密的方法（和
# Summary / Math 缺方法时一样的坑）。实现直接复用已经测过的 t() 和 %*%，
# 稀疏 × 稀疏仍然保持稀疏，稀疏 × 稠密自然退化成 eigenMat。
#' @export
crossprod.eigenSparse <- function(x, y = NULL, ...) {
  # 和 %*% 一样：真正的派发靠文件末尾的 `crossprod.eigenSparse <-
  # crossprod.eigenMat`。这里只是让 S3method 注册有个名字。
  .crossprod(x, y)
}

#' @export
tcrossprod.eigenSparse <- function(x, y = NULL, ...) {
  .tcrossprod(x, y)
}

#' @export
aperm.eigenSparse <- function(a, perm, ...) {
  d <- dim(a)
  if (missing(perm)) perm <- rev(seq_along(d))
  perm <- as.integer(perm)
  if (length(perm) != length(d) || !identical(sort(perm), seq_along(d)))
    stop("'perm' must be a permutation of seq_along(dim(a))", call. = FALSE)
  if (identical(perm, rev(seq_along(d)))) return(t(a))
  a
}

# 隐含零是有限的，所以 is.finite(S) 基本恒为 TRUE——但仍然要真的算一遍，
# 因为存下来的值里可能有 NaN/Inf。返回稠密逻辑矩阵，与 is.na.eigenSparse
# 同一口径。
#' @export
is.finite.eigenSparse <- function(x) is.finite(as.matrix(x))

#' @export
is.infinite.eigenSparse <- function(x) is.infinite(as.matrix(x))

#' @export
is.nan.eigenSparse <- function(x) is.nan(as.matrix(x))

# c() 只按**第一个**参数派发。没有 c.eigenSparse 时，c(S) 会落到 default，
# 而 default 对一个 externalptr 只能构造 list——于是 c(S) 静默返回一个
# 装着裸指针的长度 1 list，既不报错也不对。这是最坏的一类：值错、不响。
#' @export
c.eigenSparse <- function(..., recursive = FALSE) {
  args <- lapply(list(...), function(a) if (inherits(a, "eigenObject")) as.vector(a) else a)
  do.call(c, c(args, list(recursive = recursive)))
}

# isSymmetric 是真泛型，但稀疏对象没有 "eigenMat"，派发到不了稠密那份。
# 不稠密化：t(S) - S 仍然是稀疏的，直接读它的极值。
# 非方阵要先挡掉——base R 对非方阵返回 FALSE，而 t(S) - S 会直接报
# "non-conformable"，那是另一种答案。
#' @export
isSymmetric.eigenSparse <- function(object, tol = 100 * .Machine$double.eps, ...) {
  d <- dim(object)
  if (d[1] != d[2]) return(FALSE)
  st <- reigen_sparse_stats(t(object) - object)
  if (st[["na"]] + st[["nan"]] > 0) return(NA)
  max(abs(c(st[["min"]], st[["max"]]))) <= tol
}

# 稀疏矩阵的"全部元素"= 存下来的非零值 + (total - nnz) 个隐含零。分位数只
# 需要知道"排序后第 k 个是多少"，于是整件事是 O(nnz log nnz)。
#
# 关键：零**不在**排序结果的最前面。升序排列是
#     负数 ... 负数 | 零 ... 零 | 正数 ... 正数
# 所以必须把负值和正值分开，零夹在中间。第一版把零当成"排在最前的 zeros 个"，
# 结果 0% 分位数给了 0 而正确答案是最小值（负数）。这个错误在"全非负"的测试
# 矩阵上完全看不出来——那是测试夹具太客气，不是实现正确。
.sparse_quantile <- function(x, probs, na.rm, names_flag, type) {
  if (!identical(as.integer(type), 7L))
    return(quantile(as.matrix(x), probs = probs, na.rm = na.rm,
                    names = names_flag, type = type))
  st <- reigen_sparse_stats(x)
  probs <- as.numeric(probs)
  if (!length(probs) || any(!is.finite(probs)) || any(probs < 0 | probs > 1))
    stop("'probs' must be numbers in [0, 1]", call. = FALSE)
  if (st[["na"]] + st[["nan"]] > 0) {
    # base R 在 na.rm = FALSE 时直接报错；na.rm = TRUE 时退回稠密，这条路径
    # 少见，不值得为它再写一遍"跳过 NA 的第 k 个"。
    if (!isTRUE(na.rm))
      stop("missing values and NaN's not allowed if 'na.rm' is FALSE", call. = FALSE)
    return(quantile(as.matrix(x), probs = probs, na.rm = TRUE,
                    names = names_flag, type = 7))
  }
  v   <- reigen_sparse_triplets(x)$x
  neg <- sort(v[v < 0])
  pos <- sort(v[v > 0])
  nneg  <- length(neg)
  zeros <- st[["total"]] - nneg - length(pos)
  kth <- function(k) {
    # 必须对**任何** k 都给出有限值。p = 1 时 h = n、q = n、r = 0，第二个
    # 项仍然会去取第 n+1 个，而那一项乘的是 0——只要取到 NA，0 * NA 还是 NA，
    # 100% 分位数就变成 NA 了。所以越界一律夹住，空的 pos 一律给 0。
    out <- numeric(length(k))
    lo <- k <= nneg
    hi <- k > nneg + zeros
    if (any(lo)) out[lo] <- neg[pmin(nneg, k[lo])]
    if (any(hi)) {
      idx <- k[hi] - nneg - zeros
      out[hi] <- if (length(pos)) pos[pmin(length(pos), idx)] else 0
    }
    out
  }
  h <- (st[["total"]] - 1) * probs + 1
  q <- floor(h); r <- h - q
  out <- kth(q) + r * (kth(q + 1) - kth(q))
  # 名字和 base 完全一致：名字只取决于 probs，不取决于数据。
  if (isTRUE(names_flag)) names(out) <- names(quantile(0, probs = probs))
  out
}

#' @export
median.eigenSparse <- function(x, na.rm = FALSE, ...)
  .sparse_quantile(x, 0.5, na.rm, names_flag = FALSE, type = 7)

#' @export
quantile.eigenSparse <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE,
                                 names = TRUE, type = 7, ...)
  .sparse_quantile(x, probs, na.rm, names_flag = names, type = type)

# ---------------------------------------------------------------- 查询

#' @rdname eigen_sparse
#' @export
eigen_sparse_nnz <- function(x) reigen_sparse_nnz(x)

#' @rdname eigen_sparse
#' @export
eigen_triplets <- function(x) {
  out <- reigen_sparse_triplets(x)
  data.frame(i = out$i, j = out$j, x = out$x)
}

#' @rdname eigen_sparse
#' @export
eigen_sparse_get <- function(x, i, j) reigen_sparse_get(x, as.integer(i), as.integer(j))

#' @export
print.eigenSparse <- function(x, ...) {
  d <- reigen_sparse_dim(x)
  n <- reigen_sparse_nnz(x)
  total <- as.numeric(d[1]) * as.numeric(d[2])
  cat(sprintf("<eigenSparse> %d x %d, %s stored entr%s (%.1f%% dense)\n",
              d[1], d[2], format(n, big.mark = ","),
              if (n == 1) "y" else "ies",
              if (total > 0) 100 * n / total else 0))
  invisible(x)
}

#' @export
is.matrix.eigenSparse <- function(x) TRUE

#' @export
is.array.eigenSparse <- function(x) TRUE

# ---------------------------------------------------------------- 下标
#
# 单元素走 C++ 的 O(log n) 查找；区间退回 as.matrix() 再取。后者会物化整个
# 稠密矩阵，对大规模矩阵很浪费，但语义与 base R 完全一致——宁可慢，也不要
# 一个"看起来对、边界差一个"的稀疏切片。

#' @export
`[.eigenSparse` <- function(x, i, j, ..., drop = TRUE) {
  nidx <- .reigen_nindex(nargs(), missing(drop))

  if (nidx == 1L) {
    if (missing(i)) return(x)
    # 线性下标：先算 (行, 列) 再取。列主序，与 base R 一致。
    if (.is_pos_index(i)) {
      d <- dim(x); k <- as.integer(i)
      if (k <= d[1L] * d[2L])
        return(reigen_sparse_get(x, (k - 1L) %% d[1L] + 1L, (k - 1L) %/% d[1L] + 1L))
    }
    return(as.matrix(x)[i])
  }

  if (missing(i) && missing(j)) return(x)

  # 单个元素是最常见的用法，不走稠密化。
  if (!missing(i) && !missing(j) &&
      length(i) == 1L && length(j) == 1L &&
      .is_pos_index(i) && .is_pos_index(j))
    return(reigen_sparse_get(x, as.integer(i), as.integer(j)))

  m <- as.matrix(x)
  if (missing(j))      m[i, , drop = drop]
  else if (missing(i)) m[, j, drop = drop]
  else                 m[i, j, drop = drop]
}

# ---------------------------------------------------------------- 算术
#
# 稀疏 + 稀疏保持稀疏；一旦掺进稠密操作数就整体转稠密（结果本来就几乎处处
# 非零，硬撑稀疏只会更慢更占内存）。
#
# 实现放在 .ops_sparse() 里，而 Ops.eigenSparse 只是 Ops.eigenMat 的另一个
# 名字（见文件末尾）。这不是风格问题：R 的 Ops 组泛型要求两个操作数派发到
# **同一个**方法，否则警告 "Incompatible methods" 并落到 default——
# As + E 会直接报 "non-numeric argument to binary operator"。让两个类共用一个
# 闭包是唯一能让混合运算走到我们自己代码里的办法。

.ops_sparse <- function(e1, e2, op, unary) {

  # 一元：e2 缺失。这里必须是 missing(e2)，不要先算一个 "one" 标志再判空——
  # 那样很容易把 unary 和 binary 判反（本函数第一版就是这么错的：As + As
  # 走了一元分支直接返回 e1，As - As 变成了 -As）。
  if (unary) {
    if (op == "-") return(reigen_sparse_scale(e1, -1))
    if (op == "+") return(e1)
    # ! 对数值矩阵等价于 x == 0，与 Ops.eigenMat 的 !E 保持一致。
    if (op == "!") return(as.matrix(e1) == 0)
    stop(sprintf("unary '%s' is not defined for eigenSparse", op), call. = FALSE)
  }

  s1 <- inherits(e1, "eigenSparse")
  s2 <- inherits(e2, "eigenSparse")

  # 标量乘 / 除：保持稀疏。
  #
  # 这里必须排除 eigenObject：is.numeric(eigenMat) 是 TRUE，length() 也正常，
  # 所以 As * eigenMat(1x1) 会落进这个分支，然后 as.numeric() 在 externalptr
  # 上炸掉（as.numeric 是 primitive，不派发）。1x1 的稠密操作数交给下面的
  # 稠密分支去广播，语义与 base R 的 A * matrix(3, 1, 1) 一致。
  if (s1 && is.numeric(e2) && length(e2) == 1L && !inherits(e2, "eigenObject")) {
    if (op == "*") return(reigen_sparse_scale(e1, as.numeric(e2)))
    if (op == "/" && e2 != 0) return(reigen_sparse_scale(e1, 1 / as.numeric(e2)))
    # 正幂保持稀疏：隐含的零元素 0^p 在 p>0 时仍是 0，结构不变。
    # p <= 0 不能这么做（0^-1 是 Inf、0^0 是 1），落到下面的稠密分支。
    if (op == "^" && e2 > 0) return(reigen_sparse_pow(e1, as.numeric(e2)))
    # 除以 0 特意**不**在这里拦：base R 的 A/0 给出 Inf 和 NaN（零元素变
    # NaN），稀疏表示不了，必须转稠密。早期版本在这里 stop("division by
    # zero")，那是错的——稠密 eigenMat 和 base matrix 都能算，只有稀疏报错。
  }
  if (s2 && is.numeric(e1) && length(e1) == 1L && !inherits(e1, "eigenObject") &&
      op == "*")
    return(reigen_sparse_scale(e2, as.numeric(e1)))

  # 稀疏 op 稀疏
  if (s1 && s2) {
    if (op == "+") return(reigen_sparse_add(e1, e2))
    if (op == "-") return(reigen_sparse_sub(e1, e2))
    if (op == "*") return(reigen_sparse_cwise(e1, e2))
  }

  # 其余一律转稠密，交给 eigenMat 的 Ops 处理。
  #
  # 只转**稀疏的那一侧**。这里两边都过 .as_eigen() 是错的（本函数第一版就是
  # 这么写的）：标量 2 会变成 1x1 的 eigenMat，As > 2 落到 base R 的 > 时
  # 就成了 3x3 对 1x1，报 "non-conformable arrays"，而稠密侧 E > 2 正常。
  d1 <- if (s1) eigenMat(as.matrix(e1)) else e1
  d2 <- if (s2) eigenMat(as.matrix(e2)) else e2
  get(op, envir = baseenv())(d1, d2)
}

# Ops.eigenSparse 必须是**同一个闭包**，不能是另一份实现。R 的 Ops 组泛型在
# 两个操作数派发到不同方法时，会警告 "Incompatible methods" 并落到 default，
# 于是 As + E 报 "non-numeric argument to binary operator"，而两边的代码一次
# 都不会被调用。指向同一个函数之后，混合运算才能进到 .ops_sparse() 里。
# （NAMESPACE 里的 S3method(Ops, eigenSparse) 仍然要保留：派发需要这个名字。）
Ops.eigenSparse <- Ops.eigenMat

#' @export
`%*%.eigenSparse` <- function(x, y) {
  # %*% 会在**两个**操作数上派发，所以这个方法既可能是"稀疏在左"，
  # 也可能是"稀疏在右"。四种组合都得处理。
  .matmul(x, y)
}
# 上面这个定义只是为了 NAMESPACE 里的 S3method(%*%, eigenSparse) 有个名字可
# 注册——真正的派发要求**同一个闭包**。见文件末尾的 `%*%.eigenSparse <-
# %*%.eigenMat`：没有那一句，S %*% E 会警告 "Incompatible methods" 并落到
# default 报错，两个函数体一次都不会执行。
# 顺序必须和 C++ 里的 method 编号一致（reigen_sparse.cpp 的 switch）。
.reigen_sparse_methods <- c("lu", "llt", "ldlt", "qr", "cg", "bicgstab", "lsqr",
                            "gmres", "minres", "dgmres", "idrs")
# 预条件子：0 = IncompleteLUT，1 = 对角。顺序必须和 C++ 一致。
.reigen_sparse_preconds <- c("ilut", "diagonal")

# precond 的默认值取决于 method，所以不能写成形式参数默认值：
#   bicgstab -> ilut（实测步数少一个量级）
#   cg       -> diagonal（CG 要求预条件子对称正定，ILUT 不是）
# 显式传 cg + ilut 是**错的**，直接报错而不是让它悄悄不收敛。
.reigen_resolve_precond <- function(precond, method) {
  # 只接受对称预条件子的方法：CG 要求正定，MINRES 只要求对称，两者都
  # 配不上 IncompleteLUT（它不对称）。
  symmetric_only <- method %in% c("cg", "minres")
  if (is.null(precond)) {
    if (method %in% c("bicgstab", "gmres", "dgmres", "idrs")) return("ilut")
    return("diagonal")            # 直接法与对称方法：对角或无所谓
  }
  precond <- match.arg(precond, .reigen_sparse_preconds)
  if (symmetric_only && precond == "ilut")
    # R 不拼接相邻的字符串字面量（C 会），所以这里必须用 paste0。
    stop(paste0(if (method == "cg") "ConjugateGradient" else "MINRES",
                " requires a symmetric preconditioner, and IncompleteLUT is ",
                "not symmetric -- pairing them is not a valid combination ",
                "and does not converge. Use precond = \"diagonal\" (the ",
                "default for \"", method, "\"), or method = \"bicgstab\"."),
         call. = FALSE)
  precond
}

#' @export
solve.eigenSparse <- function(a, b, ...,
                              method = .reigen_sparse_methods,
                              tol = 1e-12, maxiter = NULL,
                              precond = NULL) {
  # 与 solve.eigenMat 保持同一套返回形状约定：
  #   solve(A, b) 在 b 是无 dim 向量时返回**向量**（base R 的行为），
  #   而 eigen_* 一次性函数一律返回 eigenMat。
  # 少了这一条，solve(S, b) 给 10x1 的 eigenMat、solve(as.matrix(S), b) 给向量，
  # 同一个泛型在稠密和稀疏两条路上返回不同形状。
  if (missing(b)) {
    d <- dim(a)
    if (d[1] != d[2]) stop("solve() without 'b' needs a square matrix", call. = FALSE)
    return(eigen_sparse_solve(a, eigen_identity(d[1]), method = method,
                              tol = tol, maxiter = maxiter, precond = precond))
  }
  drop <- is.null(dim(b)) && !inherits(b, "eigenMat")
  r <- eigen_sparse_solve(a, b, method = method, tol = tol, maxiter = maxiter,
                          precond = precond)
  if (!drop) return(r)
  # 降到向量时把迭代诊断一起带过去。as.vector() 会丢掉属性，而
  # solve(S, b, method = "bicgstab") 恰恰是最自然的调用方式——诊断信息
  # 不能只对 eigenMat 形式的返回值可见。
  v <- as.vector(r)
  attr(v, "iterations") <- attr(r, "iterations")
  attr(v, "residual") <- attr(r, "residual")
  v
}

#' Solve a sparse system
#'
#' Factorises once and solves. Unlike [base::solve()], which densifies, this
#' never forms the dense matrix -- that is the entire point.
#'
#' Every solver checks `Eigen::info()`. Eigen's sparse factorisations do not
#' throw on failure: they set `info()` to `NumericalIssue` and keep returning
#' values, so an unchecked call silently solves a *different* system. Every
#' method here stops instead. Direct methods are `"lu"`, `"llt"`, `"ldlt"` and
#' `"qr"`; the other seven are iterative.
#'
#' The iterative methods need a second check that the direct ones do not: when
#' they run out of iterations `info()` becomes `NoConvergence`, but a solution
#' is still returned and it looks perfectly ordinary. They therefore stop too,
#' and on success they attach the iteration count and final relative residual as
#' attributes, so you can see how much work the answer cost.
#'
#' @param a An `eigenSparse` matrix.
#' @param b Right-hand side; a vector or matrix.
#' @param method One of:
#'   * `"lu"` -- `SparseLU`, the general-purpose square solver.
#'   * `"llt"` -- `SimplicialLLT`, for symmetric positive definite matrices.
#'     It reads only the lower triangle, so a non-symmetric matrix is silently
#'     treated as its symmetric part; what it does reject is a non-positive
#'     pivot, which is how an indefinite matrix announces itself.
#'   * `"ldlt"` -- `SimplicialLDLT`, for symmetric indefinite matrices. Same
#'     lower-triangle-only reading as `"llt"`, and it accepts negative pivots.
#'   * `"qr"` -- `SparseQR`, least squares; accepts a rectangular matrix.
#'   * `"cg"` -- `ConjugateGradient`, iterative, symmetric positive definite.
#'   * `"bicgstab"` -- `BiCGSTAB`, iterative, for **non-symmetric** square
#'     systems. Use this when `"cg"` does not apply; it is the usual choice for
#'     convection-dominated problems.
#'   * `"lsqr"` -- `LeastSquaresConjugateGradient`, iterative, rectangular.
#'   * `"gmres"` -- `GMRES`, iterative, for general non-symmetric systems. Its
#'     residual decreases monotonically, which is what makes it worth the memory:
#'     one Krylov vector per iteration.
#'   * `"minres"` -- `MINRES`, iterative, for **symmetric** systems (positive
#'     definite or not). Only the lower triangle is read, so a non-symmetric
#'     matrix is silently treated as its symmetric part.
#'   * `"dgmres"` -- `DGMRES`, `GMRES` with an extra diagonal scaling; use it
#'     when the rows differ wildly in scale.
#'   * `"idrs"` -- `IDR(s)` with `s = 4`, iterative, general. Often converges
#'     where `"bicgstab"` does not, at a fixed memory cost.
#' @param tol Relative residual to stop at. Only used by the iterative methods.
#' @param maxiter Iteration cap. Only used by the iterative methods; the default
#'   is `max(1000, 10 * nrow(a))`.
#' @param precond Preconditioner for `"cg"` and `"bicgstab"`; ignored by the
#'   other methods. Defaults depend on the method:
#'   * `"bicgstab"` defaults to `"ilut"` (`IncompleteLUT`). It converges in far
#'     fewer iterations: on a 48x48 convection grid, 99 against 1352 for
#'     `"diagonal"`, and at slightly stronger convection `"diagonal"` does not
#'     converge at all. It performs an incomplete LU **without pivoting**, so it
#'     can hit a zero pivot and produce `NaN`; the solve then stops with a
#'     non-convergence error rather than returning the `NaN`.
#'   * `"gmres"`, `"dgmres"` and `"idrs"` also default to `"ilut"`.
#'   * `"cg"` and `"minres"` default to `"diagonal"` (`DiagonalPreconditioner`),
#'     and `precond = "ilut"` is an **error** for both: they require a symmetric
#'     preconditioner, and IncompleteLUT is not symmetric. The combination does
#'     not converge.
#' @return An `eigenMat` holding the solution. For the iterative methods it
#'   carries `iterations` and `residual` attributes.
#' @export
eigen_sparse_solve <- function(a, b,
                               method = c("lu", "llt", "ldlt", "qr", "cg",
                                          "bicgstab", "lsqr", "gmres",
                                          "minres", "dgmres", "idrs"),
                               tol = 1e-12, maxiter = NULL,
                               precond = NULL) {
  # 形式参数里写**字面量**，不写 .reigen_sparse_methods：R CMD check 的 codoc
  # 会逐字比对函数签名和 Rd 的 \usage，写成内部对象名就报 WARNING。
  # 代价是这份列表出现两次，所以 tests/sparse.R 里有一条断言两者相同。
  method <- match.arg(method)
  precond <- .reigen_resolve_precond(precond, method)
  code <- match(method, .reigen_sparse_methods) - 1L
  r <- reigen_sparse_solve(a, .as_eigen(b), code, as.numeric(tol),
                           if (is.null(maxiter)) 0L else as.integer(maxiter),
                           match(precond, .reigen_sparse_preconds) - 1L)
  out <- r$x
  # 迭代法把迭代次数和最终相对残差挂成属性；直接法两个都是 NA，不挂。
  # 迭代法"解出来了"和"解对了"是两件事，用户需要能看见它走了多少步。
  if (!is.na(r$iterations)) {
    attr(out, "iterations") <- r$iterations
    attr(out, "residual") <- r$error
  }
  out
}

# ---------------------------------------------------------------- 归约
#
# 稀疏对象必须自己实现 Summary / Math：类向量是 c("eigenSparse",
# "eigenObject")，其中没有 "eigenMat"，所以派发**不会**落到 eigenMat 的方法
# 上，而是落到 base 的 default——sum(As) 于是报
# "invalid 'type' (externalptr) of argument"。稠密侧这些都能用，稀疏侧不能用
# 是一种没有理由的不一致。

#' @export
Summary.eigenSparse <- function(..., na.rm = FALSE) {
  xs <- list(...)
  g <- .Generic

  # 多个参数：各自归约再合并。prod 的积、sum 的和、max 的最大值都能这样合并；
  # range 要拆成 min/max 两个分量再合。
  #
  # 注意这里走 .summary_sparse() 而不是递归调用 Summary.eigenSparse()：
  # .Generic 是组泛型派发时才绑进方法帧的，直接调用时它根本不存在，
  # 递归会立刻报 "object '.Generic' not found"。
  if (length(xs) > 1L) {
    parts <- lapply(xs, function(z) {
      if (inherits(z, "eigenSparse")) .summary_sparse(z, g, na.rm)
      else do.call(g, list(z, na.rm = na.rm))
    })
    if (g == "range") return(range(unlist(parts), na.rm = na.rm))
    return(do.call(g, c(parts, list(na.rm = na.rm))))
  }
  if (length(xs) == 0L)
    stop(sprintf("Summary '%s' needs at least one argument", g), call. = FALSE)

  .summary_sparse(xs[[1]], g, na.rm)
}

.summary_sparse <- function(x, g, na.rm) {
  st <- reigen_sparse_stats(x)
  n_na  <- st[["na"]]
  n_nan <- st[["nan"]]
  has_na <- n_na > 0 || n_nan > 0
  # base R 区分 NA 和 NaN：sum(c(1, NA)) 是 NA，sum(c(1, NaN)) 是 NaN。
  # 两者 is.na() 都为真，但 identical() 为假，所以不能一律返回 NA_real_。
  bad <- if (n_na > 0) NA_real_ else NaN

  if (has_na && !na.rm && g %in% c("sum", "prod", "max", "min"))
    return(bad)
  if (has_na && !na.rm && g == "range")
    return(c(bad, bad))

  switch(g,
    sum   = st[["sum"]],
    prod  = st[["prod"]],
    max   = st[["max"]],
    min   = st[["min"]],
    range = c(st[["min"]], st[["max"]]),
    # any/all 的 NA 语义：有决定性答案就返回它，否则 NA。
    #   any(c(NA, 1)) 是 TRUE，any(c(NA, 0)) 是 NA
    #   all(c(NA, 0)) 是 FALSE，all(c(NA, 1)) 是 NA
    any = if (st[["nz"]] > 0) TRUE else if (has_na) NA else FALSE,
    all = {
      live <- st[["total"]] - n_na - n_nan
      if (st[["nz"]] < live) FALSE else if (has_na) NA else TRUE
    },
    stop(sprintf("Summary '%s' is not defined for eigenSparse", g), call. = FALSE))
}

#' @export
mean.eigenSparse <- function(x, ...) {
  if (length(list(...)) > 0L) return(mean(as.matrix(x), ...))
  # 以前这里直接写 sum / length，缺失值被 stats 那遍跳过，于是 mean(S) 给出
  # 一个偏小的数，而 mean(as.matrix(S)) 是 NA。用同一个 .sparse_na() 判断，
  # 和 eigen_mean()、sum() 保持一致。
  st <- reigen_sparse_stats(x)
  bad <- .sparse_na(st)
  if (!is.null(bad)) return(bad)
  st[["sum"]] / length(x)
}

# 只有把 0 映射成 0 的运算能保持稀疏结构。其余（exp、log、cos、gamma、
# cumsum……）在零处给出非零值，或者依赖整列，只能转稠密。
# 注意这个集合**漏掉**某一项只会慢、不会错：没列进来的运算走稠密，结果一样。
.reigen_sparse_zero_preserving <- c(
  "abs", "sign", "sqrt", "floor", "ceiling", "trunc", "round", "signif",
  "sin", "tan", "sinpi", "tanpi", "asin", "atan", "sinh", "tanh",
  "asinh", "atanh", "log1p", "expm1"
)

.sparse_map <- function(x, f) {
  tt <- eigen_triplets(x)
  d <- dim(x)
  eigen_sparse(tt$i, tt$j, f(tt$x), nrow = d[1L], ncol = d[2L])
}

#' @export
Math.eigenSparse <- function(x, ...) {
  op <- .Generic
  if (length(list(...)) > 0L)
    return(do.call(op, c(list(as.matrix(x)), list(...))))
  if (!(op %in% .reigen_sparse_zero_preserving))
    # 稠密化之后交给 **eigenMat 的 Math 方法**，而不是直接在 base matrix 上算。
    # 后者会让 exp(As) 返回 base matrix 而 exp(E) 返回 eigenMat——同一个运算，
    # 换个输入类型就换了返回类型。带额外参数的调用（log(As, 10)）在稠密侧本来
    # 就落回 base matrix，所以上面那条保持原样。
    return(do.call(op, list(eigenMat(as.matrix(x)))))
  .sparse_map(x, function(v) do.call(op, list(v)))
}

#' @export
is.numeric.eigenSparse <- function(x) TRUE

# is.na 的结果是一个逻辑矩阵，稀疏表示不了（零处是 FALSE，但 TRUE 的位置
# 可能遍布全矩阵），与稠密侧一样返回 base matrix。
#' @export
is.na.eigenSparse <- function(x) is.na(as.matrix(x))

#' @export
anyNA.eigenSparse <- function(x, recursive = FALSE) reigen_sparse_finite(x, 1L)

#' @export
as.data.frame.eigenSparse <- function(x, row.names = NULL, optional = FALSE, ...)
  as.data.frame(as.matrix(x), row.names = row.names, optional = optional, ...)

# ---------------------------------------------------------------- 同一个闭包

# 组泛型要求两个操作数派发到**同一个**方法闭包，否则 R 警告
#   Incompatible methods ("a.eigenSparse", "a.eigenMat") for "a"
# 并落到 default 直接报错，两个函数体一次都不会执行。
#
# Ops 早就修过（见上面的 Ops.eigenSparse <- Ops.eigenMat）。%*%、crossprod、
# tcrossprod 是同一类问题，但漏了——发现它的是 README 的可执行块：
#   x <- eigen_sparse_solve(L, b, "cg");  L %*% x     # 看残差
# 这一步报 "requires numeric/complex matrix/vector arguments"。而"解完方程
# 看残差"恰恰是稀疏求解之后最自然的动作，eigen_sparse_solve 返回的正是
# eigenMat，所以这条路径几乎必然被走到。
#
# 三行赋值必须放在本文件**最后**：上面各自的函数定义只是给 NAMESPACE 里的
# S3method 一个可注册的名字，真正生效的是这里的重定向。
`%*%.eigenSparse`      <- `%*%.eigenMat`
crossprod.eigenSparse  <- crossprod.eigenMat
tcrossprod.eigenSparse <- tcrossprod.eigenMat
