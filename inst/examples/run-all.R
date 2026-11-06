#!/usr/bin/env Rscript
# ===========================================================================
# run-all.R —— 顺序执行全部 11 章教程
#
# 用法：
#   Rscript run-all.R                 # 跑全部
#   Rscript run-all.R 03 07           # 只跑第 03、07 章
#   Rscript run-all.R --quiet         # 只打印每章一行摘要
#
# 前置：ROpenEigen 已安装（R CMD INSTALL 或 install.packages(<tarball>)）。
# 包内安装位置：system.file("examples", package = "ROpenEigen")
# ===========================================================================

suppressPackageStartupMessages(library(ROpenEigen))

args <- commandArgs(trailingOnly = TRUE)
quiet <- "--quiet" %in% args
want <- setdiff(args, "--quiet")

# 脚本所在目录：Rscript 会把它放在 --file= 里
self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
dir <- if (length(self)) dirname(sub("^--file=", "", self)) else "."

files <- sort(list.files(dir, pattern = "^[0-9][0-9]-.*\\.R$", full.names = TRUE))
if (length(want)) {
  keep <- substr(basename(files), 1, 2) %in% want
  files <- files[keep]
}
if (!length(files)) stop("在 ", dir, " 下没有找到章节脚本")

cat("ROpenEigen", as.character(utils::packageVersion("ROpenEigen")), "\n")
cat("待执行", length(files), "章\n\n")

res <- data.frame(chapter = character(), secs = numeric(), status = character(),
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
  res <- rbind(res, data.frame(chapter = nm, secs = dt,
                               status = if (ok) "OK" else "FAIL",
                               stringsAsFactors = FALSE))
  if (!ok) break                       # 有一章失败就停，避免后面对不上号
}

cat(strrep("-", 75), "\n")
print(res, row.names = FALSE)
cat("\n合计", nrow(res), "章，", round(sum(res$secs), 2), "秒，失败",
    sum(res$status != "OK"), "章\n", sep = "")

if (any(res$status != "OK")) quit(status = 1)
