# Fuzz the sparse layer against base R.
#
# Why a second fuzz file instead of adding sparse cases to fuzz.R: the failure
# modes are different. Dense bugs are numerical; sparse bugs are about *which
# cells exist*. Every sparse bug found by hand this session was of the form
# "the answer is plausible but the implicit zeros were mishandled":
#
#   * quantile(S, 0) returned 0 for a matrix whose minimum was negative,
#     because the code assumed zeros sort to the front;
#   * eigen_sum(S) returned the partial sum 18 where sum(S) was NA, because
#     reigen_sparse_stats() skips missing values for the na.rm = TRUE path and
#     the callers never re-checked;
#   * c(S) returned a length-1 list holding a raw pointer;
#   * eigen_which_min() returned NA@(1,1) when the first cell was NA.
#
# None of those throw. So this file does not test that calls succeed — it tests
# that the *numbers* equal what base R gives on the same dense matrix, over
# random shapes and densities, and it deliberately generates the inputs that
# make implicit zeros matter: mostly-zero matrices, all-negative matrices,
# exact zeros mixed with non-zeros, single non-zeros, and missing values.
#
# A failure prints the seed, so it is reproducible with set.seed(seed).

suppressPackageStartupMessages(library(ROpenEigen))

ok <- 0L; bad <- 0L
chk <- function(label, expr) {
  r <- tryCatch(isTRUE(expr), error = function(e) {
    cat(sprintf("  !! %-46s %s\n", label, conditionMessage(e))); NA
  })
  if (is.na(r)) { bad <<- bad + 1L; return(invisible(FALSE)) }
  if (r) ok <<- ok + 1L else { bad <<- bad + 1L; cat(sprintf("  ?? %s\n", label)) }
  invisible(r)
}

# 相对误差比较，和其它几个测试文件同一口径。
#
# 注意先把 eigenObject 过一遍 as.matrix()：as.numeric()/as.vector() 是
# **primitive，不派发 S3**，直接喂 externalptr 会得到
# "cannot coerce type 'externalptr' to vector of type 'double'"。
# as.matrix() 是泛型，会派发。这个坑第一次写这个文件时就踩了。
rel <- function(a, b, tol = 1e-9) {
  if (inherits(a, "eigenObject")) a <- as.matrix(a)
  if (inherits(b, "eigenObject")) b <- as.matrix(b)
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b)) return(FALSE)
  if (anyNA(a) || anyNA(b)) return(identical(is.na(a), is.na(b)))
  if (!length(a)) return(TRUE)
  d <- max(abs(a - b)); s <- max(1, max(abs(a)), max(abs(b)))
  d <= tol * s
}
wlin <- function(w, nr)
  if (is.na(w$row)) NA_integer_ else as.integer((w$col - 1L) * nr + w$row)

cat("== 随机稀疏矩阵：每个入口对稠密 ==")

SEEDS <- 1:60
for (seed in SEEDS) {
  set.seed(seed)
  nr <- sample(1:12, 1); nc <- sample(1:12, 1)
  cells <- nr * nc
  density <- runif(1, 0.02, 1)
  nz <- max(1L, round(cells * density))
  idx <- sample(cells, min(nz, cells))

  A <- matrix(0, nr, nc)
  vals <- switch(sample(4L, 1),
    rnorm(length(idx)),
    runif(length(idx), -3, 3),
    as.numeric(sample(-5:5, length(idx), replace = TRUE)),
    abs(rnorm(length(idx))))
  A[idx] <- vals
  # 有相当比例的用例带上显式零：存下来的零和隐含零在数值上没区别，
  # 但会走到三元组/统计的不同分支上。
  if (runif(1) < 0.4 && length(idx) > 2)
    A[sample(idx, max(1L, length(idx) %/% 3))] <- 0

  S <- eigen_sparse(A)
  tag <- sprintf("seed %d (%dx%d, %.2f)", seed, nr, nc, density)

  # --- 结构 ---
  chk(sprintf("%s as.matrix", tag), rel(as.matrix(S), A, 0))
  chk(sprintf("%s dim", tag), identical(as.integer(dim(S)), c(as.integer(nr), as.integer(nc))))
  chk(sprintf("%s length", tag), identical(length(S), cells))
  chk(sprintf("%s nnz", tag), identical(as.numeric(eigen_sparse_nnz(S)), as.numeric(sum(A != 0))))

  # --- 转置与逐元素 ---
  chk(sprintf("%s t()", tag), rel(t(S), t(A), 0))
  chk(sprintf("%s S + S", tag), rel(S + S, A + A, 1e-12))
  chk(sprintf("%s S * 2", tag), rel(S * 2, A * 2, 1e-12))
  chk(sprintf("%s S - S", tag), rel(S - S, A - A, 1e-12))
  chk(sprintf("%s -S", tag), rel(-S, -A, 1e-12))
  chk(sprintf("%s abs(S)", tag), rel(abs(S), abs(A), 1e-12))

  # --- 归约 ---
  chk(sprintf("%s sum", tag), rel(eigen_sum(S), sum(A)))
  chk(sprintf("%s prod", tag), rel(eigen_prod(S), prod(A)))
  chk(sprintf("%s mean", tag), rel(eigen_mean(S), mean(A)))
  chk(sprintf("%s sumsq", tag), rel(eigen_sumsq(S), sum(A^2)))
  chk(sprintf("%s trace", tag), rel(eigen_trace(S), sum(diag(A))))
  chk(sprintf("%s norm F", tag), rel(eigen_norm(S, "F"), norm(A, "F")))
  chk(sprintf("%s norm M", tag), rel(eigen_norm(S, "M"), norm(A, "M")))
  chk(sprintf("%s norm O", tag), rel(eigen_norm(S, "O"), norm(A, "O")))
  chk(sprintf("%s norm I", tag), rel(eigen_norm(S, "I"), norm(A, "I")))
  # 谱范数是迭代出来的，容差放宽到 1e-8
  chk(sprintf("%s norm 2", tag), rel(eigen_norm(S, "2"), norm(A, "2"), 1e-8))

  # --- 极值位置 ---
  chk(sprintf("%s which_min 位置", tag),
      identical(wlin(eigen_which_min(S), nr), which.min(A)))
  chk(sprintf("%s which_max 位置", tag),
      identical(wlin(eigen_which_max(S), nr), which.max(A)))

  # --- 对角线 ---
  chk(sprintf("%s diag", tag), rel(as.vector(eigen_diag(S)), diag(A)))

  # --- 分位数 ---
  chk(sprintf("%s median", tag), rel(median(S), median(A)))
  chk(sprintf("%s quantile", tag), rel(as.numeric(quantile(S)), as.numeric(quantile(A))))

  # --- 乘法的四种形状 ---
  B <- matrix(rnorm(nc * 3), nc, 3)
  chk(sprintf("%s S %%*%% B", tag), rel(S %*% B, A %*% B, 1e-10))
  D <- matrix(rnorm(3 * nr), 3, nr)
  chk(sprintf("%s D %%*%% S", tag), rel(D %*% S, D %*% A, 1e-10))
  chk(sprintf("%s crossprod", tag), rel(crossprod(S), crossprod(A), 1e-10))
  chk(sprintf("%s tcrossprod", tag), rel(tcrossprod(S), tcrossprod(A), 1e-10))

  # --- 对称性 ---
  if (nr == nc) chk(sprintf("%s isSymmetric", tag), identical(isSymmetric(S), isSymmetric(A)))

  # --- 缺失值：每个归约都要把 NA 传上去，而不是给一个偏小的数 ---
  if (cells >= 4 && runif(1) < 0.35) {
    A2 <- A; A2[sample(cells, 1)] <- NA
    S2 <- eigen_sparse(A2)
    t2 <- sprintf("%s+NA", tag)
    chk(sprintf("%s sum", t2), identical(eigen_sum(S2), sum(A2)))
    chk(sprintf("%s prod", t2), identical(eigen_prod(S2), prod(A2)))
    chk(sprintf("%s mean", t2), identical(eigen_mean(S2), mean(A2)))
    chk(sprintf("%s sumsq", t2), identical(eigen_sumsq(S2), sum(A2^2)))
    chk(sprintf("%s trace", t2), identical(eigen_trace(S2), sum(diag(A2))))
    chk(sprintf("%s sum() 泛型", t2), identical(sum(S2), sum(A2)))
    chk(sprintf("%s mean() 泛型", t2), identical(mean(S2), mean(A2)))
    chk(sprintf("%s which_min 跳过 NA", t2),
        identical(wlin(eigen_which_min(S2), nr), which.min(A2)))
  }

  # --- NaN 与 NA 不能混为一谈 ---
  if (cells >= 4 && runif(1) < 0.2) {
    A3 <- A; A3[sample(cells, 1)] <- NaN
    S3 <- eigen_sparse(A3)
    t3 <- sprintf("%s+NaN", tag)
    chk(sprintf("%s sum", t3), identical(eigen_sum(S3), sum(A3)))
    chk(sprintf("%s prod", t3), identical(eigen_prod(S3), prod(A3)))
    chk(sprintf("%s mean", t3), identical(eigen_mean(S3), mean(A3)))
  }
}

# 一组刻意构造的极端输入。随机矩阵很少生成这些，而它们正是隐含零最容易
# 出错的地方。
cat("\n== 构造的极端输入 ==")
extreme <- list(
  "全零"          = matrix(0, 5, 4),
  "1x1 零"        = matrix(0, 1, 1),
  "1x1 非零"      = matrix(-7, 1, 1),
  "单列"          = matrix(c(-3, 0, 0, 0), 4, 1),
  "单行"          = matrix(c(-3, 0, 0, 0), 1, 4),
  "单个非零"      = { m <- matrix(0, 7, 7); m[4, 5] <- -2.5; m },
  "全负"          = -matrix(abs(rnorm(20)), 4, 5),
  "全正"          = matrix(abs(rnorm(20)), 4, 5),
  "满阵"          = matrix(rnorm(20), 4, 5),
  "首格即最小"    = { m <- matrix(1, 5, 4); m[1, 1] <- -9; m },
  "最小值在末格"  = { m <- matrix(1, 5, 4); m[5, 4] <- -9; m },
  "首列全非零"    = { m <- matrix(0, 5, 4); m[, 1] <- c(1, 2, 3, 4, 5); m },
  "空列在前"      = { m <- matrix(0, 4, 3); m[, 3] <- -1; m },
  "显式零"        = { m <- matrix(0, 4, 4); m[2, 2] <- 0; m[1, 3] <- -1; m },
  "对角"          = diag(c(3, -7, 1, 9, -2))
)
for (nm in names(extreme)) {
  A <- extreme[[nm]]; S <- eigen_sparse(A)
  chk(sprintf("%-14s as.matrix", nm), rel(as.matrix(S), A, 0))
  chk(sprintf("%-14s sum/prod/mean", nm),
      rel(eigen_sum(S), sum(A)) && rel(eigen_prod(S), prod(A)) && rel(eigen_mean(S), mean(A)))
  chk(sprintf("%-14s which_min", nm), identical(wlin(eigen_which_min(S), nrow(A)), which.min(A)))
  chk(sprintf("%-14s which_max", nm), identical(wlin(eigen_which_max(S), nrow(A)), which.max(A)))
  chk(sprintf("%-14s median/quantile", nm),
      rel(median(S), median(A)) && rel(as.numeric(quantile(S)), as.numeric(quantile(A))))
  chk(sprintf("%-14s diag", nm), rel(as.vector(eigen_diag(S)), diag(A)))
}

# 组泛型要求两个操作数派发到**同一个**闭包。没有这条时 R 会警告
# Incompatible methods 然后落到 default，两边的方法体都不执行——结果不是报错，
# 是一个看起来对的答案。断言"没有 warning"比断言"没有 error"更严：
# warning 出现在 error 之前。
cat("\n== 组泛型不产生 Incompatible methods 警告 ==")
E <- eigenMat(matrix(rnorm(12), 4, 3))
S <- eigen_sparse(matrix(rnorm(12), 4, 3))
for (expr in list(quote(S + E), quote(S %*% t(E)), quote(crossprod(S, E)),
                  quote(tcrossprod(S, E)), quote(S * E), quote(S - E))) {
  lab <- deparse(expr)
  chk(sprintf("无警告: %s", lab), {
    warned <- FALSE
    withCallingHandlers(eval(expr),
      warning = function(w) {
        if (grepl("Incompatible methods", conditionMessage(w))) warned <<- TRUE
        invokeRestart("muffleWarning")
      })
    !warned
  })
}

cat(sprintf("\n通过 %d, 失败 %d\n", ok, bad))
if (bad > 0L) stop(sprintf("%d assertion(s) failed", bad), call. = FALSE)
