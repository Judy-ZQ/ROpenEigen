#!/usr/bin/env Rscript
# ===========================================================================
# run-all.R —— 顺序执行《线性代数》配套的全部实验（lab00 – lab07）
#
# 用法：
#   Rscript run-all.R                 # 跑全部 8 个实验
#   Rscript run-all.R 01 03 05        # 只跑实验一、三、五
#   Rscript run-all.R --quiet         # 只打印每个实验一行摘要
#
# 前置：ROpenEigen 已安装（R CMD INSTALL 或 install.packages(<tarball>)）。
# 包内安装位置：system.file("labs", package = "ROpenEigen")
#
# 说明：每个实验都在独立环境里 source，互不干扰；某个实验失败不会中断
#       其余实验，最后统一汇总，有失败则以状态码 1 退出。
# ===========================================================================

suppressPackageStartupMessages(library(ROpenEigen))

args <- commandArgs(trailingOnly = TRUE)
quiet <- "--quiet" %in% args
want <- setdiff(args, "--quiet")
want <- sub("^lab", "", want)                 # 允许写 lab01 或 01

# 脚本所在目录：Rscript 会把它放在 --file= 里
# 脚本所在目录：Rscript 把路径放进 --file=，source() / RStudio 把它放进 ofile。
# 两者都没有（例如把整段粘进控制台）时退回当前目录。
here <- function() {
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) return(dirname(sub("^--file=", "", f[1])))
  for (fr in rev(sys.frames())) {
    of <- fr$ofile
    if (is.character(of) && length(of) == 1L && file.exists(of)) return(dirname(of))
  }
  getwd()
}
dir <- here()

files <- sort(list.files(dir, pattern = "^lab[0-9][0-9]-.*\\.R$", full.names = TRUE))
if (length(want)) {
  files <- files[substr(basename(files), 4, 5) %in% want]
}
if (!length(files)) stop("在 ", dir, " 下没有找到实验脚本（可用编号 00–07）")

cat("ROpenEigen", as.character(utils::packageVersion("ROpenEigen")), "\n")
cat("R", as.character(getRversion()), "|", .Platform$OS.type, "\n")
cat("实验目录：", dir, "\n", sep = "")
cat("待执行", length(files), "个实验\n\n")

res <- data.frame(lab = character(), secs = numeric(), status = character(),
                  stringsAsFactors = FALSE)

for (f in files) {
  nm <- sub("\\.R$", "", basename(f))
  if (!quiet) cat(strrep("=", 75), "\n")
  cat("[", nm, "] 开始\n", sep = "")
  t0 <- proc.time()[["elapsed"]]
  ok <- tryCatch({
    if (quiet) {
      utils::capture.output(source(f, echo = FALSE, local = new.env()))
    } else {
      source(f, echo = FALSE, local = new.env())
    }
    TRUE
  }, error = function(e) {
    cat("  !! 失败：", conditionMessage(e), "\n", sep = "")
    FALSE
  })
  dt <- round(proc.time()[["elapsed"]] - t0, 2)
  cat("[", nm, "] ", if (ok) "OK" else "FAIL", "  ", dt, "s\n\n", sep = "")
  res <- rbind(res, data.frame(lab = nm, secs = dt,
                               status = if (ok) "OK" else "FAIL",
                               stringsAsFactors = FALSE))
}

cat(strrep("-", 75), "\n")
print(res, row.names = FALSE)
cat("\n合计 ", nrow(res), " 个实验，", round(sum(res$secs), 2), " 秒，失败 ",
    sum(res$status != "OK"), " 个\n", sep = "")

if (any(res$status != "OK")) {
  cat("失败的实验：", paste(res$lab[res$status != "OK"], collapse = ", "), "\n", sep = "")
  quit(status = 1)
}
