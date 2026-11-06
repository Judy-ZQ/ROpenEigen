suppressPackageStartupMessages(library(ROpenEigen))

# 《线性代数》实验教程（inst/labs）的回归测试。
#
# 这份测试和其他测试文件的分工不同：前面那些验的是「封装本身对不对」，
# 这一份验的是「教材配套的 8 个实验还能不能从头跑到尾」。每个实验内部都用
# helpers.R 里的 expect_near() 把手算结果与库函数对账，断言不过就抛错，
# 所以这里只需要顺序跑一遍并汇总。实验本身是给人看的，输出全部捕获掉。

dir <- system.file("labs", package = "ROpenEigen")
if (!nzchar(dir)) stop("inst/labs 没有随包安装")
files <- sort(list.files(dir, pattern = "^lab[0-9][0-9]-.*\\.R$", full.names = TRUE))
if (length(files) != 8L) {
  stop("应有 8 个实验（lab00-lab07），实际找到 ", length(files), " 个")
}
if (!file.exists(file.path(dir, "helpers.R"))) stop("缺少 helpers.R")

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
if (bad > 0L) stop(sprintf("%d 个实验未通过", bad), call. = FALSE)
