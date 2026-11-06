# eigenMat：Eigen 稠密矩阵在 R 侧的表现
#
# 对象是 EXTPTRSXP，class = c("eigenMat", "eigenObject")。数据在 C++ 侧，
# R 只持有一个指针；as.matrix() 才落回 R 的矩阵。这与 ROpenCV 的 cvMat 一致。
#
# 这里把 base R 的矩阵惯用法全部接过来：dim / length / t / solve / det /
# norm / %*% / Ops / Math / Summary / Re / Im / [ / [<-。目标是不用背新 API：
# 凡是能写成 A %*% B 的，就不该逼用户写成 eigen_matmul(A, B)。

#' Construct an Eigen matrix
#'
#' @param x A numeric, integer, logical or complex matrix or vector, or an
#'   existing `eigenMat` (returned unchanged).
#' @param nrow,ncol Optional dimensions, interpreted as in [matrix()].
#' @return An `eigenMat` reference object.
#' @export
eigenMat <- function(x = 0, nrow = NULL, ncol = NULL) {
  if (inherits(x, "eigenMat") && is.null(nrow) && is.null(ncol)) return(x)
  if (is.data.frame(x)) x <- data.matrix(x)
  if (!is.null(nrow) || !is.null(ncol)) {
    if (is.null(ncol)) ncol <- ceiling(length(x) / nrow)
    if (is.null(nrow)) nrow <- ceiling(length(x) / ncol)
    x <- matrix(x, nrow = nrow, ncol = ncol)
  }
  reigen_mat_new(x)
}

#' Construct an Eigen column vector
#'
#' Same storage as [eigenMat()], but the class vector carries `eigenVec` first
#' so that `[` and `print` can treat it as a vector.
#'
#' @param x A numeric or complex vector.
#' @return An `eigenVec` object.
#' @export
eigenVec <- function(x = numeric(0)) {
  reigen_vec_new(x)
}

#' @export
dim.eigenMat <- function(x) reigen_mat_dim(x)

#' @export
`dim<-.eigenMat` <- function(x, value) {
  if (length(value) != 2L) stop("dim must be length 2", call. = FALSE)
  reigen_mat_reshape(x, as.integer(value[1]), as.integer(value[2]))
}

#' @export
length.eigenMat <- function(x) {
  d <- dim(x)
  as.integer(d[1] * d[2])
}

#' @export
as.matrix.eigenMat <- function(x, ...) reigen_mat_to_r(x)

#' @export
as.vector.eigenMat <- function(x, mode = "any") {
  as.vector(reigen_mat_to_r(x), mode = mode)
}

# c() 只对第一个参数做 S3 派发。没有这个方法时 c(E, 1) 返回一个 list
# （base R 对非向量对象就是这样），而 c(m, 1) 返回展平后的数值向量——
# 又一处"不报错但结果不一样"。注意 c(1, E) 派发不到这里，那是 c() 的规则。
#' @export
c.eigenMat <- function(..., recursive = FALSE) {
  args <- lapply(list(...), function(a) if (inherits(a, "eigenMat")) as.vector(a) else a)
  do.call(c, c(args, list(recursive = recursive)))
}

#' @export
is.matrix.eigenMat <- function(x) TRUE

# 矩阵本来就是数组。base R 的 is.array(matrix) 返回 TRUE；这里若返回 FALSE，
# 会和 is.matrix() 自相矛盾，而且 outer() 之类的函数会据此走错分支
# （outer 对数组和向量分别处理，见其源码里的 is.array 判断）。
#' @export
is.array.eigenMat <- function(x) TRUE

#' @export
is.numeric.eigenMat <- function(x) !reigen_mat_is_complex(x)

# ---------------------------------------------------------------- 转置

#' @export
t.eigenMat <- function(x) reigen_mat_transpose(x)

# ---------------------------------------------------------------- 打印

#' @export
print.eigenMat <- function(x, ..., max_rows = 10L, max_cols = 8L) {
  d <- dim(x)
  cat(sprintf("<eigenMat> %d x %d %s\n", d[1], d[2],
              if (reigen_mat_is_complex(x)) "complex" else "double"))
  if (d[1] == 0L || d[2] == 0L) return(invisible(x))

  m <- as.matrix(x)
  if (d[1] > max_rows) {
    h <- max_rows %/% 2L
    ri <- c(seq_len(h), seq.int(d[1] - (max_rows - h) + 1L, d[1]))
    m <- m[ri, , drop = FALSE]
  }
  if (d[2] > max_cols) {
    h <- max_cols %/% 2L
    ci <- c(seq_len(h), seq.int(d[2] - (max_cols - h) + 1L, d[2]))
    m <- m[, ci, drop = FALSE]
  }
  print(m, ...)
  if (d[1] > max_rows || d[2] > max_cols)
    cat(sprintf("  (%d x %d total; use as.matrix() for the whole matrix)\n", d[1], d[2]))
  invisible(x)
}

# ---------------------------------------------------------------- 下标

.eigen_row_out <- function(y) as.vector(reigen_mat_to_r(y))
.eigen_col_out <- function(y) as.vector(reigen_mat_to_r(y))

.is_pos_index <- function(i) {
  is.numeric(i) && length(i) == 1L && !is.na(i) && i >= 1 && i == as.integer(i)
}

# R 靠参数个数区分 x[i] 与 x[i, ]：在方法里两者都是 j missing，差别只在
# 传了几个下标参数。漏掉这个判断，x[i] 会被当成行下标——于是 E[1] 返回第一行
# （4 个数而不是 1 个），E[5] 直接报越界，diag(E) 也跟着挂：diag 内部取对角线
# 用的是纯线性下标 x[1 + 0:(m-1)*(nrow+1)]。
# nargs() 把 x 自己、已传的下标、以及（若显式给出）drop 都算进去。
.reigen_nindex <- function(nargs, drop_missing) nargs - 1L - !drop_missing

# base R 里 m[rows, cols] <- v 是把 v 按**列主序**铺进这个块，v 可以是向量。
# C++ 侧的 block_set 要的是一块形状确定的矩阵，所以向量必须先在 R 里铺好。
# 少了这一步，`m[1, ] <- 1:3` 会被当成 3x1 去填 1x3 的块，报
# "value is 3x1, block is 1x3" —— 而同样一行在 base R 矩阵上是合法的。
.reigen_block_value <- function(value, nr, nc) {
  if (is.matrix(value)) return(value)
  matrix(rep_len(value, nr * nc), nr, nc)
}

#' @export
`[.eigenMat` <- function(x, i, j, ..., drop = TRUE) {
  nidx <- .reigen_nindex(nargs(), missing(drop))

  # 线性下标 x[k]：按列主序展平后的第 k 个元素，与 base R 完全一致。
  # 矩阵下标 x[cbind(i, j)]、负下标、逻辑下标、越界（返回 NA）都走 as.matrix。
  if (nidx == 1L) {
    if (missing(i)) return(x)
    if (.is_pos_index(i)) {
      d <- dim(x); k <- as.integer(i)
      if (k <= d[1L] * d[2L])
        return(reigen_mat_get(x, (k - 1L) %% d[1L] + 1L, (k - 1L) %/% d[1L] + 1L))
    }
    return(as.matrix(x)[i])
  }

  if (missing(i) && missing(j)) return(x)          # x[] 或 x[,]

  # 快路径：单行 / 单列，不落回 R 矩阵，也不整表拷贝。
  if (missing(j) && .is_pos_index(i)) return(.eigen_row_out(reigen_mat_row(x, as.integer(i))))
  if (missing(i) && .is_pos_index(j)) return(.eigen_col_out(reigen_mat_col(x, as.integer(j))))

  # 其余（负下标、逻辑下标、范围、两维同取）交给 base R，语义完全一致。
  m <- as.matrix(x)
  if (missing(j))      m[i, , drop = drop]
  else if (missing(i)) m[, j, drop = drop]
  else                 m[i, j, drop = drop]
}

#' @export
`[<-.eigenMat` <- function(x, i, j, ..., value) {
  # 右值也可能是 eigenMat：m[1:2, 1:2] <- eigen_identity(2)。base R 的子赋值
  # 先看 length(value) 再决定怎么铺，而 externalptr 的 length 是 1，于是报
  # "number of items to replace is not a multiple of replacement length"。
  # 在 C++ 侧的指针外面套着 dim，看起来像矩阵，实际不是——先摊成真矩阵。
  if (inherits(value, "eigenMat")) value <- as.matrix(value)
  # nargs() 里已经含 x 和 value，剩下的才是下标个数：E[k] <- v 是 3 -> 1 个，
  # E[i, j] <- v 是 4 -> 2 个，E[] <- v 是 2 -> 0 个。
  nidx <- nargs() - 2L
  d <- dim(x)

  if (nidx <= 0L) {                              # x[] <- v
    m <- as.matrix(x); m[] <- value
    return(eigenMat(m))
  }

  # 线性下标赋值 E[k] <- v 写的是第 k 个元素，不是第 k 行。
  if (nidx == 1L) {
    if (missing(i)) { m <- as.matrix(x); m[] <- value; return(eigenMat(m)) }
    if (.is_pos_index(i)) {
      k <- as.integer(i)
      if (k <= d[1L] * d[2L])
        return(reigen_mat_set(x, (k - 1L) %% d[1L] + 1L, (k - 1L) %/% d[1L] + 1L, value))
    }
    m <- as.matrix(x); m[i] <- value
    return(eigenMat(m))
  }

  if (missing(j) && .is_pos_index(i))
    return(reigen_mat_block_set(x, as.integer(i), as.integer(i), 1L, d[2],
                                .reigen_block_value(value, 1L, d[2])))

  if (missing(i) && .is_pos_index(j))
    return(reigen_mat_block_set(x, 1L, d[1], as.integer(j), as.integer(j),
                                .reigen_block_value(value, d[1], 1L)))

  if (!missing(i) && !missing(j) && .is_pos_index(i) && .is_pos_index(j))
    return(reigen_mat_set(x, as.integer(i), as.integer(j), value))

  m <- as.matrix(x)
  if (missing(j))      m[i, ...] <- value
  else if (missing(i)) m[, j, ...] <- value
  else                 m[i, j, ...] <- value
  eigenMat(m)
}

# ---------------------------------------------------------- 交回 base R
#
# 这些函数在矩阵上的语义要么是"按行"、要么是"当成一个向量看"，在 C++ 里
# 重写一遍没有收益，还容易写错。统一做法：转成 base 矩阵（或向量）后原样
# 调用同名函数，参数照传——语义由 base R 保证。
#
# 返回的也是 base 对象，不是 eigenMat：这样 unique(E) 与 unique(m) 的
# identical() 为 TRUE。要留在 Eigen 侧的场景有 eigen_* 系列函数。

#' @export
sort.eigenMat <- function(x, decreasing = FALSE, na.last = NA, ...)
  sort(as.vector(x), decreasing = decreasing, na.last = na.last, ...)

#' @export
rev.eigenMat <- function(x) rev(as.vector(x))

#' @export
median.eigenMat <- function(x, na.rm = FALSE, ...)
  median(as.vector(x), na.rm = na.rm, ...)

#' @export
quantile.eigenMat <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE,
                              names = TRUE, type = 7, ...)
  quantile(as.vector(x), probs = probs, na.rm = na.rm, names = names,
           type = type, ...)

# which.min/which.max 不是 S3 泛型（base R 里没有 which.min.default，也没有
# UseMethod），挂方法不会被调用——实测仍然报 "cannot unclass an external
# pointer"。所以只能用 eigen_which_min()/eigen_which_max()。

# unique()/duplicated() 在矩阵上是**按行**的：base R 的 unique(matrix) 返回
# 去重后的行组成的矩阵，不是去重后的元素。所以这里必须转成矩阵交给 base，
# 自己用 as.vector() 去重会得到完全不同（而且静默错误）的结果。
#' @export
unique.eigenMat <- function(x, incomparables = FALSE, ...)
  unique(as.matrix(x), incomparables = incomparables, ...)

#' @export
duplicated.eigenMat <- function(x, incomparables = FALSE, ...)
  duplicated(as.matrix(x), incomparables = incomparables, ...)

#' @export
as.data.frame.eigenMat <- function(x, row.names = NULL, optional = FALSE, ...)
  as.data.frame(as.matrix(x), row.names = row.names, optional = optional, ...)

#' @export
isSymmetric.eigenMat <- function(object, tol = 100 * .Machine$double.eps, ...)
  isSymmetric(as.matrix(object), tol = tol, ...)

# diff() 是泛型，但 diff.default 内部对 x 做了 unclass()，EXTPTRSXP 上会报
# "cannot unclass an external pointer"。矩阵上 diff 是相邻行相减。
#' @export
diff.eigenMat <- function(x, lag = 1L, differences = 1L, ...)
  diff(as.matrix(x), lag = lag, differences = differences, ...)

# rbind() 在 EXTPTRSXP 上会报 "cannot coerce type 'externalptr' to vector of
# type 'list'"（cbind 反而能跑，因为两者的内部实现不同）。补上方法让两边一致，
# 结果仍然是 base 矩阵，与 cbind 的行为对齐。
#' @export
rbind.eigenMat <- function(..., deparse.level = 1)
  do.call(rbind, c(lapply(list(...), .as_r), list(deparse.level = deparse.level)))

#' @export
cbind.eigenMat <- function(..., deparse.level = 1)
  do.call(cbind, c(lapply(list(...), .as_r), list(deparse.level = deparse.level)))

# ---------------------------------------------------------------- 算术

# .Generic 到 C++ 算子码的映射，与 src/reigen_ops.cpp 的 switch 一致。
.reigen_cwise2_code <- c("+" = 0L, "-" = 1L, "*" = 2L, "/" = 3L, "^" = 4L)

#' @export
Ops.eigenMat <- function(e1, e2) {
  op <- .Generic

  # 稀疏操作数交给 .ops_sparse（R/sparse.R）。为什么这个分支要长在 eigenMat
  # 的方法里，而不是让 eigenSparse 自己有个独立方法：R 的 Ops 组泛型会在两个
  # 操作数派发到**不同**方法时警告 "Incompatible methods" 并落到 default，
  # 于是 As + E 直接报 "non-numeric argument to binary operator"，而两边
  # 各自的代码根本不会被调用。把 Ops.eigenSparse 指向同一个闭包
  # （见 R/sparse.R 末尾）之后，两个操作数派发到同一个函数，R 才肯调。
  if (inherits(e1, "eigenSparse") || (!missing(e2) && inherits(e2, "eigenSparse")))
    return(.ops_sparse(e1, if (missing(e2)) NULL else e2, op, missing(e2)))

  if (missing(e2)) {
    if (op == "+") return(e1)
    if (op == "-") return(reigen_cwise2(eigenMat(0), .as_eigen(e1), 1L))
    # ! 是逻辑取反，对数值矩阵而言等价于 x == 0（TRUE 表示元素为 0）。
    # 注意不要再套一层 !：!E 的值就是 (E == 0)，不是 !(E == 0)。
    if (op == "!") return(.as_r(e1) == 0)
    stop(sprintf("unary '%s' is not defined for eigenMat", op), call. = FALSE)
  }

  # 用 [ + match 而不是 [[：[[ 在名字不存在时抛 "subscript out of bounds"，
  # 而 "==" "<" "%%" 这些本来就不该走 C++，是要落到 base R 的。
  code <- .reigen_cwise2_code[match(op, names(.reigen_cwise2_code))]
  if (!is.na(code))
    return(reigen_cwise2(.as_eigen(e1), .as_eigen(e2), code))

  # 比较、整除、取模：转回 base R。这些不在热路径上，
  # 且 R 的 NA 传播规则比在 C++ 里重写一遍更可信。
  f <- get(op, envir = baseenv())
  f(.as_r(e1), .as_r(e2))
}

# ---------------------------------------------------------------- 矩阵乘法

# %*% 是**组泛型**：它在两个操作数上派发，而且要求两边落到**同一个**方法闭包
# 上。稠密和稀疏各有一个方法时，S %*% E 会警告
#   Incompatible methods ("%*%.eigenSparse", "%*%.eigenMat") for "%*%"
# 然后落到 default 直接报错——两边的函数体一次都不会执行。
# 这和当初 Ops 组泛型上 As + E 失败是**同一个坑**，只是换了个泛型。
#
# 所以 %*%.eigenSparse 必须指向和 %*%.eigenMat **同一个闭包**（见 R/sparse.R
# 末尾），四种组合的逻辑集中在这里。
.matmul <- function(x, y) {
  xs <- inherits(x, "eigenSparse")
  ys <- inherits(y, "eigenSparse")

  if (xs && ys) return(reigen_sparse_matmul_sparse(x, y))

  if (xs) {
    # base R 把无维右操作数当列向量。
    if (is.null(dim(y))) y <- matrix(y, ncol = 1L)
    return(reigen_sparse_matmul(x, y))
  }

  if (ys) {
    # 稠密 %*% 稀疏：Eigen 的 Sparse 模块没有稠密×稀疏的直接乘积，
    # 用 ( (S^T) (D^T) )^T 转一次。维度：D(m x n) S(n x p) -> m x p。
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    return(t(reigen_sparse_matmul(reigen_sparse_transpose(y), t(x))))
  }

  # 稠密 × 稠密。base R 的 %*% 对无维向量是**不对称**的：
  #   b %*% A  把 b 当行向量（等价 t(b) %*% A），结果是 1 x ncol(A)
  #   A %*% b  把 b 当列向量，结果是 nrow(A) x 1
  # 而 .as_eigen() 一律把无维输入当列向量（那是 eigenMat 自己的约定，
  # 对 solve()/算术都合适）。不在这里单独纠正左操作数，b %*% E 就会报
  # "non-conformable arguments: 2x1 %*% 2x2"——base R 却能算出来。
  if (is.null(dim(x)) && !inherits(x, "eigenObject")) x <- matrix(x, nrow = 1L)
  reigen_matmul(.as_eigen(x), .as_eigen(y))
}

#' @export
`%*%.eigenMat` <- function(x, y) .matmul(x, y)

# ------------------------------------------------- crossprod / tcrossprod

# crossprod 和 tcrossprod 是**内部泛型**：body 里既没有 UseMethod 也不在
# Math/Ops 组里，dispatch 发生在 .Internal 内部。靠 grep body 判断"是不是
# 泛型"会得出错误结论——实测（scratch/probe_dispatch2.R）它们确实会派发到
# crossprod.eigenMat，所以这里能挂方法。
#
# 它们同样在**两个**操作数上派发，所以和 %*% 一样，稀疏侧必须共用同一个闭包。
#
# 关键差别：这两个函数把无维输入一律当**列**向量，而 %*% 把左边的无维向量
# 当行。所以不能写成 t(x) %*% y 了事——那对向量操作数会给出转置的结果。
.as_col <- function(z) {
  if (!is.null(dim(z))) return(z)
  matrix(as.vector(z), ncol = 1L)
}

.crossprod <- function(x, y) {
  if (is.null(y)) y <- x
  t(.as_col(x)) %*% .as_col(y)
}

.tcrossprod <- function(x, y) {
  if (is.null(y)) y <- x
  .as_col(x) %*% t(.as_col(y))
}

#' @export
crossprod.eigenMat <- function(x, y = NULL, ...) .crossprod(x, y)

#' @export
tcrossprod.eigenMat <- function(x, y = NULL, ...) .tcrossprod(x, y)

# aperm 是**真**泛型（body 里有 UseMethod），所以能挂。这里只有二维矩阵，
# 因此 perm 只可能是 c(2,1)（转置）或 c(1,2)（恒等）。不做三维推广，是因为
# eigenMat 根本表示不了三维——与其猜，不如在 perm 不合法时报错。
#' @export
aperm.eigenMat <- function(a, perm, ...) {
  d <- dim(a)
  if (missing(perm)) perm <- rev(seq_along(d))
  perm <- as.integer(perm)
  if (length(perm) != length(d) || !identical(sort(perm), seq_along(d)))
    stop("'perm' must be a permutation of seq_along(dim(a))", call. = FALSE)
  if (identical(perm, rev(seq_along(d)))) return(t(a))
  a
}

# is.finite / is.infinite / is.nan 也是**内部泛型**，而且和 Math 组无关：
# 给 Math.eigenMat 加分支是没用的（实测报 "default method not implemented for
# type 'externalptr'"，说明根本没进 Math.eigenMat）。它们各自单独派发。
# 返回逻辑矩阵而不是 eigenMat——eigenMat 只存 double，装不下逻辑值；
# 这一点与 is.na.eigenMat 保持一致。
#' @export
is.finite.eigenMat <- function(x) is.finite(as.matrix(x))

#' @export
is.infinite.eigenMat <- function(x) is.infinite(as.matrix(x))

#' @export
is.nan.eigenMat <- function(x) is.nan(as.matrix(x))

# ---------------------------------------------------------------- Math 组

# Math 组泛型一次覆盖 abs/sqrt/exp/log/log10/sin/cos/tan/.../round/floor/
# ceiling/sign/log1p/expm1。带额外参数（log(x, base)、round(x, 2)）的调用
# 落回 base R。
.reigen_math_code <- c(
  abs = 0L, sqrt = 1L, exp = 2L, log = 3L, log10 = 4L,
  sin = 5L, cos = 6L, tan = 7L, asin = 8L, acos = 9L, atan = 10L,
  sinh = 11L, cosh = 12L, tanh = 13L, round = 14L, floor = 15L,
  ceiling = 16L, sign = 17L, log1p = 22L, expm1 = 23L
)

#' @export
Math.eigenMat <- function(x, ...) {
  op <- .reigen_math_code[match(.Generic, names(.reigen_math_code))]
  if (is.na(op) || length(list(...)) > 0L)
    return(do.call(.Generic, c(list(as.matrix(x)), list(...))))
  reigen_cwise1(x, op)
}

# ---------------------------------------------------------------- Summary 组

.reigen_summary_code <- c(sum = 0L, prod = 1L, min = 3L, max = 4L,
                          all = 6L, any = 7L)

#' @export
Summary.eigenMat <- function(..., na.rm = FALSE) {
  xs <- list(...)
  op <- .reigen_summary_code[match(.Generic, names(.reigen_summary_code))]
  # C++ 的归约是直接对 Eigen 表达式求值的，不处理 NA；na.rm = TRUE 时退回
  # base R（as.matrix 后调基础函数），保证语义正确。代价是一次拷贝，但
  # na.rm 不是热路径。多参数的情形本来就走这条路。
  if (is.na(op) || length(xs) != 1L || isTRUE(na.rm))
    return(do.call(.Generic, c(lapply(xs, .as_r), list(na.rm = na.rm))))
  reigen_reduce(xs[[1]], op)
}

#' @export
mean.eigenMat <- function(x, trim = 0, na.rm = FALSE, ...) {
  if (isTRUE(na.rm) || trim != 0) return(mean(.as_r(x), trim = trim, na.rm = na.rm))
  reigen_reduce(x, 2L)
}

# ---------------------------------------------------------------- 复数

#' @export
Re.eigenMat <- function(z) reigen_cwise1(z, 25L)

#' @export
Im.eigenMat <- function(z) reigen_cwise1(z, 26L)

#' @export
Mod.eigenMat <- function(z) reigen_cwise1(z, 0L)

#' @export
Arg.eigenMat <- function(z) reigen_cwise1(z, 27L)

#' @export
Conj.eigenMat <- function(z) reigen_cwise1(z, 24L)

# ---------------------------------------------------------------- 线性代数

#' @export
solve.eigenMat <- function(a, b, ...) {
  if (missing(b)) {
    d <- dim(a)
    if (d[1] != d[2]) stop("solve() without 'b' needs a square matrix", call. = FALSE)
    return(reigen_solve(a, reigen_mat_identity(d[1])))
  }
  # base R 的 solve(A, b) 在 b 是「无 dim 的向量」时返回向量，不是 n x 1 矩阵。
  # 这里跟着它走，否则 A %*% solve(A, b) 的维度会悄悄变掉。
  drop <- is.null(dim(b)) && !inherits(b, "eigenMat")
  r <- reigen_solve(a, .as_eigen(b))
  if (drop) as.vector(r) else r
}

# determinant() 是标准泛型，可以挂。返回对象必须严格照 base R 的约定：
#   list(modulus, sign)，class = "det"
# 且 logarithm = TRUE（默认）时 modulus 是 log|det|，不是 |det|。
#
# 这一点错了会静默放大成指数级错误：det() 本身不是泛型，它的定义是
#   det <- function(x, ...) c(z$sign * exp(z$modulus))   # z = determinant(x, ...)
# 也就是说 det() 会调用我们这个方法，再把 modulus 取 exp。如果这里返回
# |det| 而声称 logarithm = TRUE，det(E) 就变成 exp(det(A))——2x2 的对角阵
# 6 会返回 403.43，而且不报错。
#
# sign 同样不能恒为 1：行列式为负时 base R 给 -1。复数没有 base R 参考
# （determinant() 对复矩阵直接报错），这里取单位复数 det/|det|，
# 使得 det = sign * exp(modulus) 对复数也成立。
#' @export
determinant.eigenMat <- function(x, logarithm = TRUE, ...) {
  if (!identical(dim(x)[1], dim(x)[2]))
    stop("determinant requires a square matrix", call. = FALSE)
  p <- reigen_det_parts(x)
  sign <- if (p[[3]] == 0) p[[2]] else complex(real = p[[2]], imaginary = p[[3]])
  modulus <- if (isTRUE(logarithm)) p[[1]] else exp(p[[1]])
  # base R 还在 modulus 上挂了一个 logarithm 属性（见 determinant.matrix），
  # 不带上它的话 all.equal(determinant(x), determinant(as.matrix(x))) 会报
  # "Attributes: < target is NULL, current is list >"。
  attr(modulus, "logarithm") <- isTRUE(logarithm)
  structure(list(modulus = modulus, sign = sign), class = "det")
}

#' @export
anyNA.eigenMat <- function(x, recursive = FALSE) reigen_finite(x, 1L)

#' @export
is.na.eigenMat <- function(x) {
  m <- as.matrix(x)
  is.na(m)
}

# ---------------------------------------------------------------- 便捷构造

#' Identity, zero, one and random matrices, and the diagonal
#'
#' Thin wrappers over the Eigen static constructors. `eigen_random()` uses R's
#' RNG, so `set.seed()` makes it reproducible -- Eigen's own `Random()` draws
#' from `std::rand()` and would silently ignore the seed.
#'
#' [eigen_diag()] reads the diagonal of a matrix as a vector, and `eigen_diag<-`
#' writes one. Neither is a base R generic, so `diag()` cannot be hooked; see
#' [eigen-decompositions] for the same caveat about `det()` and `norm()`.
#'
#' @param n Size for [eigen_identity()], and the number of points for
#'   [eigen_linspace()].
#' @param nrow,ncol Dimensions.
#' @param x A matrix or `eigenMat` whose diagonal [eigen_diag()] reads, or the
#'   matrix that `eigen_diag<-` modifies.
#' @param value For [eigen_constant()], the fill value. For `eigen_diag<-`, the
#'   vector written to the diagonal, recycled as in base R.
#' @param from,to Endpoints for [eigen_linspace()]; both are included.
#' @return An `eigenMat`.
#' @name eigen-generators
NULL

#' @rdname eigen-generators
#' @export
eigen_identity <- function(n = 1L) reigen_mat_identity(as.integer(n))

#' @rdname eigen-generators
#' @export
eigen_zeros <- function(nrow = 1L, ncol = 1L) reigen_mat_zeros(as.integer(nrow), as.integer(ncol))

#' @rdname eigen-generators
#' @export
eigen_ones <- function(nrow = 1L, ncol = 1L) reigen_mat_ones(as.integer(nrow), as.integer(ncol))

#' @rdname eigen-generators
#' @export
eigen_constant <- function(value = 0, nrow = 1L, ncol = 1L)
  reigen_mat_constant(as.integer(nrow), as.integer(ncol), as.numeric(value))

#' @rdname eigen-generators
#' @export
eigen_random <- function(nrow = 1L, ncol = 1L) reigen_mat_random(as.integer(nrow), as.integer(ncol))

#' @rdname eigen-generators
#' @export
eigen_linspace <- function(from, to, n = 50L)
  reigen_mat_linspace(as.numeric(from), as.numeric(to), as.integer(n))

#' @rdname eigen-generators
#' @export
eigen_diag <- function(x) {
  # 稀疏侧只读一条对角线，O(k log nnz)，不需要稠密化。返回形状与稠密侧一致
  # （k x 1 的 eigenMat），否则同一个函数在两条路上给出两种类型。
  if (inherits(x, "eigenSparse"))
    return(eigenMat(matrix(reigen_sparse_diag(x), ncol = 1L)))
  reigen_mat_diag(.as_eigen(x))
}

#' @rdname eigen-generators
#' @export
`eigen_diag<-` <- function(x, value) reigen_mat_diag_set(.as_eigen(x), value)
