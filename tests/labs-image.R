suppressPackageStartupMessages(library(ROpenEigen))

# 图像实验教程（inst/labs-image）的回归测试。
#
# 与 tests/labs.R 同一思路：那 8 个实验自己会用 imghelpers.R 里的
# expect()/expect_near() 把手算结果和库函数（ROpenCV / ROpenEigen）对账，
# 断言不过就抛错，所以这里只要顺序跑一遍并汇总。
#
# 两点与普通测试不同：
#   1. 实验依赖 ROpenCV。它只是 Suggests，没装就跳过，不能算失败；
#   2. 实验会往磁盘写 PNG。运行前把 LINALG_FIGS 指到临时目录，
#      免得把 imgfigs/ 写进安装目录（R CMD check 会为此报 NOTE）。
#
# 这里**故意不 setwd()**：实验脚本会自己定位 imghelpers.R（脚本路径 -> 当前目录
# 及往上的源码树 -> 已安装的包）。从 tests/ 目录 source 一个绝对路径的脚本，
# 走的是「脚本自身路径」那一条 —— 和 RStudio 里按 Source 完全同一条路。

if (!requireNamespace("ROpenCV", quietly = TRUE)) {
  cat("ROpenCV 未安装，跳过图像实验的回归测试\n")
} else {

  dir <- system.file("labs-image", package = "ROpenEigen")
  if (!nzchar(dir)) stop("inst/labs-image 没有随包安装")
  if (!file.exists(file.path(dir, "imghelpers.R"))) stop("缺少 imghelpers.R")

  files <- sort(list.files(dir, pattern = "^img[0-9][0-9]-.*\\.R$", full.names = TRUE))
  if (length(files) != 8L) {
    stop("应有 8 个图像实验（img00-img07），实际找到 ", length(files), " 个")
  }

  old_figs <- Sys.getenv("LINALG_FIGS", unset = NA_character_)
  Sys.setenv(LINALG_FIGS = file.path(tempdir(), "labs-image-figs"))
  on.exit({
    if (is.na(old_figs)) Sys.unsetenv("LINALG_FIGS")
    else Sys.setenv(LINALG_FIGS = old_figs)
  }, add = TRUE)

  ok <- 0L; bad <- 0L
  for (f in files) {
    nm <- sub("\\.R$", "", basename(f))
    r <- tryCatch({
      utils::capture.output(source(f, echo = FALSE, local = new.env()))
      TRUE
    }, error = function(e) {
      cat("  !! ", nm, ": ", conditionMessage(e), "\n", sep = "")
      FALSE
    })
    if (r) { ok <- ok + 1L; cat("ok   ", nm, "\n", sep = "") }
    else bad <- bad + 1L
  }

  cat(sprintf("\n通过 %d, 失败 %d\n", ok, bad))
  if (bad > 0L) stop(sprintf("%d 个图像实验未通过", bad), call. = FALSE)
}
