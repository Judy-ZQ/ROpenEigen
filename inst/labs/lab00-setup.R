# ===========================================================================
# 实验零 · 熟悉 ROpenEigen
#
# 对应教材：全书通用（第一次上机）
# 建议学时：1 学时
#
# 实验目的
#   1. 能在本机跑通 ROpenEigen，会读它的输出；
#   2. 掌握「Eigen 按行写、R 按列填」这一条最容易出错的地方；
#   3. 知道 ROpenEigen 对象和 base R 矩阵怎么互转，以及为什么不能混着用；
#   4. 建立数值误差的观念：精确为零与 1e-16 的区别。
#
# 运行：Rscript lab00-setup.R   或在 RStudio 里 source()
# ===========================================================================

# 定位 helpers.R：脚本自己所在目录 -> 当前目录（及往上的包源码树）-> 已安装的包。
# 三种入口都能跑：Rscript labXX.R、在 RStudio 里 source()、包内回归测试。
.lab_dir <- local({
  cand <- character(0)
  # (1) 脚本自己在哪里：Rscript 把路径放进 --file=，source() / RStudio 把它放进 ofile
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) cand <- dirname(sub("^--file=", "", f[1]))
  for (fr in rev(sys.frames())) {
    of <- fr$ofile
    if (is.character(of) && length(of) == 1L && file.exists(of)) {
      cand <- c(cand, dirname(of)); break
    }
  }
  # (2) 当前目录，以及从当前目录往上 4 层找包源码树（<root>/inst/labs）
  d <- getwd()
  for (i in 0:4) {
    cand <- c(cand, d, file.path(d, "inst", "labs"), file.path(d, "labs"))
    d <- dirname(d)
  }
  # (3) 已安装的包（安装副本可能比源码树旧，所以放在最后）
  cand <- c(cand, system.file("labs", package = "ROpenEigen"))
  hit <- cand[file.exists(file.path(cand, "helpers.R"))]
  if (!length(hit))
    stop("找不到 helpers.R。三种跑法都支持：在 inst/labs 下 Rscript、",
         "在 RStudio 里 source()、或先安装 ROpenEigen。", call. = FALSE)
  hit[1]
})
source(file.path(.lab_dir, "helpers.R"))
suppressPackageStartupMessages(library(ROpenEigen))

# ===========================================================================
# 一、环境
# ===========================================================================
say("== 0.1 环境 ==")
say("ROpenEigen ", as.character(packageVersion("ROpenEigen")),
    " | ", R.version.string)

# ===========================================================================
# 二、最容易出错的一步：矩阵怎么写进 R
# ===========================================================================
# 教材里的矩阵都是按行写的，例如
#         [ 1  2  3 ]
#    A =  [ 4  5  6 ]
# R 的 matrix() 默认按列填，所以要加 byrow = TRUE。漏了这个参数就得到转置。
A_wrong <- matrix(1:6, 2, 3)                   # 按列：1 3 5 / 2 4 6
A_right <- matrix(1:6, 2, 3, byrow = TRUE)     # 按行：1 2 3 / 4 5 6
say("matrix(1:6, 2, 3) 得到：")
print(A_wrong)
say("matrix(1:6, 2, 3, byrow = TRUE) 得到：")
print(A_right)

# ===========================================================================
# 三、ROpenEigen 的对象是什么
# ===========================================================================
A <- eigenMat(A_right)
say("== 0.2 对象 ==")
say("class(A) = ", paste(class(A), collapse = " / "))
say("typeof(A) = ", typeof(A), "  —— C++ 侧的引用对象，不是 double 向量")

print(A)                                       # print 方法会先报维数与类型

# 取内容必须显式转换：这是 R 与 C++ 之间唯一的一次显式拷贝
say("as.matrix(A) 是普通矩阵：", is.matrix(as.matrix(A)))
expect_near(as.matrix(A), A_right, 0, "as.matrix 往返丢数据")

# ===========================================================================
# 四、base R 的泛型直接派发
# ===========================================================================
# ROpenEigen 的设计是「让 base R 的习惯写法直接用」，所以下面这些不用换函数名：
say("== 0.3 派发给 base R 泛型 ==")
b <- eigenVec(c(1, 1, 1))
say("A %*% b      = ", paste(round(as.vector(A %*% b), 12), collapse = ", "))
say("dim(t(A))    = ", paste(dim(t(A)), collapse = " x "))
say("sum(A)       = ", sum(A))                  # 全部元素之和
say("A[2, 3]      = ", A[2, 3])                 # 1-based 下标
say("A[, 2]       = ", paste(as.vector(A[, 2]), collapse = ", "))

# 方阵才有行列式：非方阵会明确报错，而不是返回一个错误的值
say("对 2x3 矩阵求 det：")
r <- try(det(A), silent = TRUE)
say("  -> ", if (inherits(r, "try-error")) conditionMessage(attr(r, "condition")) else paste(r))

# ===========================================================================
# 五、数值误差：什么时候该写 == ，什么时候不能
# ===========================================================================
say("== 0.4 误差观念 ==")
# 0.1 + 0.2 在浮点下不等于 0.3，这是线性代数实验里天天要面对的事
say("0.1 + 0.2 == 0.3 ？", 0.1 + 0.2 == 0.3)
say("0.1 + 0.2 - 0.3  = ", 0.1 + 0.2 - 0.3)

# 一个看起来「显然成立」的等式，用浮点算出来差多少：
H <- outer(1:8, 1:8, function(i, j) 1 / (i + j - 1))   # Hilbert 矩阵
say("8 阶 Hilbert 矩阵的 det = ", det(H))
say("它的条件数倒数 rcond = ", eigen_rcond(eigenMat(H)))
say("H %*% solve(H) 与单位阵的最大偏差 = ",
    format(max(abs(H %*% solve(H) - diag(8))), scientific = TRUE))

# 结论：线性代数实验里判断「等于」，一律用相对误差，helpers.R 里的
# expect_near() 就是干这件事的。后面每个实验的自动检查都靠它。

# ===========================================================================
# 六、本教程的公共工具
# ===========================================================================
say("== 0.5 helpers.R ==")
say("逆序数 inversions(c(3,4,2,1)) = ", inversions(c(3, 4, 2, 1)))
say("行最简形 rref() 的秩 = ", rref(A_right)$rank, "（2x3 矩阵，秩为 2）")
say("零空间 nullspace() 的维数 = ", ncol(nullspace(A_right)), "（n - r = 3 - 2 = 1）")

# ===========================================================================
# 七、上机小测（做完再看下一节）
# ===========================================================================
# 1. 用一行代码构造 4 阶单位矩阵，要求写出两种写法。
# 2. 构造矩阵 B = [[2, 0], [0, 3]]，验证 B %*% B 与 B * B 是否相同，为什么？
# 3. Hilbert 矩阵的 det 随阶数怎么变？算 4、6、8、10 阶看看，说明这对解方程
#    意味着什么。
#
# 参考写法（取消注释即可运行）：
# print(eigen_identity(4)); print(diag(4))
# B <- eigenMat(matrix(c(2, 0, 0, 3), 2, 2, byrow = TRUE))
# print(as.matrix(B %*% B)); print(as.matrix(B * B))
# for (n in c(4, 6, 8, 10)) {
#   Hn <- outer(1:n, 1:n, function(i, j) 1 / (i + j - 1))
#   cat(n, det(Hn), eigen_rcond(eigenMat(Hn)), "\n")
# }

say("\nlab00 结束：环境、录入约定、对象模型、误差观念，四件事齐了就可以做实验一。")
