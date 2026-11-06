# 对象方法注册表（必须排在所有注册调用之前被求值）
#
# 与 ROpenCV 的 .ropecv_cls 同构：方法表按类名集中放在一个环境里，
# 而不是挂在对象上。理由在这里更直接——分解对象（eigenLLT 等）是
# EXTPTRSXP，R 侧没有地方可以挂函数；而且同一个分解既可能由
# eigen_llt(A) 造出来，也可能由 eigenMat(A)$llt() 造出来，两条路径
# 必须共用同一张方法表。
#
# 表结构：
#   .reigen_cls$eigenLLT$m$solve  ->  function(self, b = NULL) ...
#
# 方法一律以 self 为第一个形参；$.eigenObject 在取用时把接收者绑进去。
#
# 与 ROpenCV 的差别（有意为之）：
#   ROpenCV 里那套 .ropecv_pick 打分分派是为了化解 OpenCV 成百上千的
#   C++ 重载。Eigen 的公开面没有这个问题——每个分解方法的签名是固定的，
#   重载只在 solve()/eigen_svd() 这类少数地方出现，用 inherits() 判断
#   比走一遍打分表更短、也更好读。所以这里保留了注册表和 $ 取用的
#   人机接口，省掉了打分器。

.reigen_cls <- new.env(parent = emptyenv())

.reigen_reg <- function(cls, new = NULL, m = list()) {
  assign(cls, list(new = new, m = m), envir = .reigen_cls)
  invisible(NULL)
}

.reigen_methods <- function(x) {
  e <- .reigen_cls[[class(x)[1]]]
  if (is.null(e)) character(0) else names(e$m)
}

# 把 R 侧的任何输入统一成 eigenMat，供 C++ 入口使用。
# 已经是 eigenMat 就原样返回，不产生拷贝。
.as_eigen <- function(x, what = deparse(substitute(x))) {
  if (inherits(x, "eigenMat")) return(x)
  # 稀疏对象必须在这里拦住。注意 is.atomic(externalptr) 是 TRUE，所以下面那条
  # 通用检查拦不住它，会一路掉进 C++ 报 "expected a numeric matrix or vector,
  # got externalptr"——那句话既没提"稀疏"，也没说该怎么办。
  # 这一条同时覆盖了三十多个只吃稠密矩阵的入口。
  if (inherits(x, "eigenSparse"))
    stop(sprintf(paste0("%s is an eigenSparse, but this function needs a dense matrix. ",
                        "Use as.matrix(%s) to densify, or a sparse-aware entry point ",
                        "(eigen_sum, eigen_prod, eigen_mean, eigen_sumsq, eigen_norm, ",
                        "eigen_trace, eigen_diag, eigen_adjoint, eigen_sparse_solve)."),
                 what, what), call. = FALSE)
  if (is.data.frame(x)) x <- data.matrix(x)
  if (!is.atomic(x) && !is.matrix(x))
    stop(sprintf("%s must be numeric, matrix or eigenMat; got %s", what, class(x)[1]),
         call. = FALSE)
  reigen_mat_new(x)
}

# 把 eigenMat 落回 base R，用于走不通 C++ 的兜底路径（比较、取模等）。
# 判断用 "eigenObject" 而不是 "eigenMat"：稀疏对象的 class 向量里没有
# "eigenMat"，只判 eigenMat 会把它原样（externalptr）交给 base R，于是
# sum(E, S) 报 "invalid 'type' (externalptr) of argument"——而 sum(S, E)
# 却能算，因为 Summary 组泛型**只按第一个参数派发**，两边走的是不同方法。
.as_r <- function(x) {
  if (inherits(x, "eigenObject")) as.matrix(x) else x
}
