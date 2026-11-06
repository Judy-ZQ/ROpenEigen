# 所有 eigen 对象共有的 R 侧行为
#
# 对象分两类：
#   eigenMat   —— 矩阵本身（含 eigenVec）
#   eigenLLT / eigenLDLT / eigenLU / eigenQR / eigenSVD / eigenEigen ...
#              —— 分解结果，算一次、反复用
#
# 分解对象的方法通过 $ 取用：llt$solve(b)、svd$singularValues()。
# 这与 ROpenCV 的 img$apply(...) 是同一套人机接口。

#' Test for an Eigen object
#' @param x Any object.
#' @return A length-1 logical.
#' @export
is_eigen <- function(x) inherits(x, "eigenObject")

#' @export
`$.eigenObject` <- function(x, name) {
  cls <- class(x)[1]
  e <- .reigen_cls[[cls]]
  f <- if (is.null(e)) NULL else e$m[[name]]
  if (is.null(f)) {
    avail <- .reigen_methods(x)
    stop(sprintf("%s has no method '%s'.\n  Available: %s", cls, name,
                 if (length(avail)) paste(avail, collapse = ", ") else "(none)"),
         call. = FALSE)
  }
  function(...) f(x, ...)
}

#' @export
names.eigenObject <- function(x) .reigen_methods(x)

#' @export
print.eigenObject <- function(x, ...) {
  cls <- class(x)[1]
  ms <- .reigen_methods(x)
  cat(sprintf("<%s> Eigen object\n", cls))
  if (length(ms)) cat("  methods: ", paste(ms, collapse = ", "), "\n", sep = "")
  else            cat("  methods: (none)\n")
  invisible(x)
}
