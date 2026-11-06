# ===========================================================================
# 综合实验 · 四个把全书串起来的题目
#
# 对应教材：跨章综合（每章都有一节带星号的「应用」，这里把它们做实）
#   一、奇异值分解与图像压缩       （第二、五、六章：秩、特征值、正定性）
#   二、最小二乘曲线拟合           （第四章：不相容方程组）
#   三、主成分分析                 （第五章、第六章：协方差矩阵的特征分解）
#   四、PageRank                   （第五章：特征向量、幂迭代）
# 建议学时：4 学时（可作为课程设计或大作业）
#
# 实验目的
#   1. 把「秩」「特征值」「正定」这些概念真正用在一个完整问题上；
#   2. 学会看误差随参数的变化，而不是只看一个数；
#   3. 体会「同一个问题多种算法」时，选哪个取决于数值稳定性。
#
# 运行：Rscript lab07-project.R   或在 RStudio 里 source()
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
# 一、奇异值分解与图像压缩
# ===========================================================================
say("== 7.1 图像压缩 ==")
# 造一张 240x320 的「图像」：一个亮圆盘 + 条纹 + 噪声。
# 真实照片也是这个结构——低频（大奇异值）承载形状，高频（小奇异值）承载细节。
nr <- 240; nc <- 320
g <- expand.grid(i = seq_len(nr), j = seq_len(nc))
img <- matrix(0, nr, nc)
disk <- (g$i - nr * 0.45)^2 / (nr * 0.28)^2 + (g$j - nc * 0.5)^2 / (nc * 0.35)^2
img[] <- ifelse(disk < 1, 200, 30)
img <- img + 40 * sin(g$j / 12) * (disk < 1.6)     # 条纹
set.seed(1)
img <- img + matrix(rnorm(nr * nc, 0, 8), nr, nc)   # 噪声
say("图像尺寸：", nr, " x ", nc, "，像素值范围 [", round(min(img)), ", ",
    round(max(img)), "]")

S <- eigen_svd(eigenMat(img))
sv <- as.vector(S$singular_values())
say("奇异值个数：", length(sv))
say("前 10 个奇异值：", paste(round(sv[1:10], 3), collapse = ", "))
say("奇异值衰减到 1e-3 以下需要 ", sum(sv > 1e-3 * sv[1]), " 个")

# 低秩逼近：取前 k 个奇异值重构
U <- as.matrix(S$u()); V <- as.matrix(S$v())
approx_k <- function(k) U[, 1:k, drop = FALSE] %*%
  diag(sv[1:k], k) %*% t(V[, 1:k, drop = FALSE])

say("\nk    重构误差(Frobenius)   相对误差   需要存储的数   原图存储")
for (k in c(1, 5, 10, 20, 40, 80, 160)) {
  Ak <- approx_k(k)
  say(sprintf("%-4d %18.4f %11.4f %13d %10d",
              k, eigen_norm(eigenMat(img - Ak), "F"),
              eigen_norm(eigenMat(img - Ak), "F") / eigen_norm(eigenMat(img), "F"),
              k * (nr + nc + 1), nr * nc))
}
# 关键结论：误差 = 被丢掉的奇异值的平方和开根（Eckart-Young 定理）
k <- 20
tail_sv <- sv[(k + 1):length(sv)]
say("\nEckart-Young 定理检验（k = ", k, "）：")
say("  实际重构误差 = ", round(eigen_norm(eigenMat(img - approx_k(k)), "F"), 6))
say("  被丢弃奇异值的平方和开根 = ", round(sqrt(sum(tail_sv^2)), 6))
expect_near(eigen_norm(eigenMat(img - approx_k(k)), "F"),
            sqrt(sum(tail_sv^2)), 1e-8, "Eckart-Young 不成立")

say("压缩比（k=20）：", round(nr * nc / (k * (nr + nc + 1)), 2), " : 1")
say("这说明「图像的大部分信息集中在少数几个奇异值上」。")

# ===========================================================================
# 二、最小二乘曲线拟合
# ===========================================================================
say("\n== 7.2 最小二乘拟合 ==")
# 拟合 y = c0 + c1 x + c2 x^2，数据点故意带噪声
set.seed(2)
x <- seq(0, 4, length.out = 20)
ytrue <- 1 - 2 * x + 0.8 * x^2
y <- ytrue + rnorm(length(x), 0, 0.5)
say("数据：", length(x), " 个点，真实模型 y = 1 - 2x + 0.8x^2（加了噪声）")

# 设计矩阵（Vandermonde 形式）
X <- cbind(1, x, x^2)
say("设计矩阵 X 的维数：", paste(dim(X), collapse = " x "), "，秩 = ", eigen_rank(eigenMat(X)))

# 三种解法
c_ne <- as.vector(solve(eigenMat(t(X) %*% X), eigenVec(t(X) %*% y)))   # 正规方程
c_qr <- as.vector(solve(eigenMat(X), eigenVec(y)))                     # QR
c_svd <- as.vector(eigen_pinv(eigenMat(X)) %*% eigenVec(y))            # SVD

say("正规方程解：", paste(round(c_ne, 6), collapse = ", "))
say("QR 解      ：", paste(round(c_qr, 6), collapse = ", "))
say("SVD 解     ：", paste(round(c_svd, 6), collapse = ", "))
say("三者最大差异：", format(max(abs(c_ne - c_svd)), scientific = TRUE))

resid <- function(cc) sum((X %*% cc - y)^2)
say("残差平方和：", round(resid(c_ne), 8), " / ", round(resid(c_qr), 8),
    " / ", round(resid(c_svd), 8))
say("条件数：cond(X) = ", round(1 / eigen_rcond(eigenMat(X)), 2),
    "，cond(X'X) = ", round(1 / eigen_rcond(eigenMat(t(X) %*% X)), 2),
    "（平方关系，这是正规方程的弱点）")

# 病态情形：在 [0,1] 上拟合 8 次多项式。Vandermonde 型设计矩阵的条件数随
# 次数指数增长，正规方程又把它平方一次，于是「正规方程」开始丢有效数字。
say("\n换成 8 次多项式拟合（15 个点，x 取在 [0,1]）：")
x3 <- seq(0, 1, length.out = 15)
y3 <- cos(3 * x3) + rnorm(15, 0, 1e-3)
X3 <- outer(x3, 0:8, function(a, b) a^b)
c_ne3 <- as.vector(solve(eigenMat(t(X3) %*% X3), eigenVec(t(X3) %*% y3)))
c_svd3 <- as.vector(eigen_pinv(eigenMat(X3)) %*% eigenVec(y3))
say("  cond(X)  = ", format(1 / eigen_rcond(eigenMat(X3)), scientific = TRUE))
say("  cond(X'X)= ", format(1 / eigen_rcond(eigenMat(t(X3) %*% X3)), scientific = TRUE),
    "  ← 平方关系放大到 1e11 量级")
say("  正规方程的前几个系数：", paste(round(c_ne3[1:5], 6), collapse = ", "))
say("  SVD 解的前几个系数    ：", paste(round(c_svd3[1:5], 6), collapse = ", "))
say("  两者最大差 = ", format(max(abs(c_ne3 - c_svd3)), scientific = TRUE))
say("  也就是说：cond ≈ 1e11 时，正规方程只剩 16 - 11 = 5 位十进制有效数字。")

# ===========================================================================
# 三、主成分分析（PCA）
# ===========================================================================
say("\n== 7.3 主成分分析 ==")
# 造二维数据：主要沿 45° 方向展开，另一方向很窄
set.seed(3)
n <- 300
Z <- matrix(rnorm(n * 2), n, 2) %*% matrix(c(2, 1.6, 1.6, 2), 2, 2)
Z <- scale(Z, scale = FALSE)
say("数据：", n, " 个二维点")

Cov <- crossprod(Z) / (n - 1)
say("协方差矩阵：")
print(round(Cov, 4))
say("它是实对称半正定的：特征值 ", paste(round(as.vector(eigen_eigen(eigenMat(Cov))$values()), 6),
                                          collapse = ", "))

ep <- eigen_eigen(eigenMat(Cov))
lam <- as.vector(ep$values()); W <- as.matrix(ep$vectors())
say("特征值（= 各主成分的方差）：", paste(round(lam, 6), collapse = ", "))
say("解释方差比：", paste(round(lam / sum(lam) * 100, 2), "%", collapse = " / "))

# 投影到第一主成分：把二维降到一维，再重构，看损失了多少
score <- Z %*% W[, 2]                      # 最大特征值对应最后一列
recon <- score %*% t(W[, 2])
# 重构误差的平方和 / (n-1) 恰好等于被丢弃的那个特征值
sse <- sum((Z - recon)^2)
say("降到一维并重构：误差平方和 = ", round(sse, 6),
    "，除以 (n-1) 得 ", round(sse / (n - 1), 6))
say("被丢弃的方差（较小的特征值） = ", round(lam[1], 6))
expect_near(sse / (n - 1), lam[1], 1e-8, "PCA 重构误差与特征值不符")

# 验证主成分两两正交、且是协方差矩阵的特征向量
say("W'W = E：", relerr(t(W) %*% W, diag(2)) < 1e-12)
say("Cov %*% w = λ w：", relerr(Cov %*% W[, 2], lam[2] * W[, 2]) < 1e-12)
say("第一主成分方向：", paste(round(W[, 2], 6), collapse = ", "),
    "（数据沿 45° 方向展开，与此相符）")

# ===========================================================================
# 四、PageRank
# ===========================================================================
say("\n== 7.4 PageRank ==")
# 6 个网页的链接关系：M[i,j] = 1 表示 j 链到 i
links <- matrix(c(0, 0, 1, 0, 0, 0,
                  1, 0, 0, 0, 0, 0,
                  0, 1, 0, 0, 0, 0,
                  0, 0, 1, 0, 0, 0,
                  0, 0, 1, 0, 0, 1,
                  0, 0, 0, 0, 1, 0), 6, 6, byrow = TRUE)
say("链接矩阵 M（M[i,j]=1 表示 j 指向 i）：")
print(links)
outdeg <- colSums(links)
say("各页面的出度：", paste(outdeg, collapse = ", "))
say("（有出度为 0 的页面，按惯例把所有页面都当成它的出链，避免「悬挂节点」）")
dangling <- which(outdeg == 0)
for (j in dangling) { links[, j] <- 1; outdeg[j] <- nrow(links) }

# 转移矩阵（列随机）
G <- sweep(links, 2, outdeg, "/")
say("列随机矩阵 G 的列和：", paste(round(colSums(G), 12), collapse = ", "))

# 幂迭代：r_{k+1} = d G r_k + (1-d)/n
damp <- 0.85
n <- nrow(G)
r <- rep(1 / n, n)
say("阻尼系数 d = ", damp, "，迭代过程：")
for (it in 1:200) {
  r_new <- damp * as.vector(G %*% r) + (1 - damp) / n
  if (it %in% c(1, 2, 5, 10, 50, 200)) {
    say(sprintf("  第 %3d 次：max|Δr| = %.3e", it, max(abs(r_new - r))))
  }
  if (max(abs(r_new - r)) < 1e-14) { r <- r_new; break }
  r <- r_new
}
say("收敛后 r = ", paste(round(r, 6), collapse = ", "))
say("总和为 1：", round(sum(r), 12))
say("PageRank 方程 r = dG r + (1-d)/n 的残差：",
    format(max(abs(r - (damp * as.vector(G %*% r) + (1 - damp) / n))), scientific = TRUE))

# 与「解线性方程组」的做法对照：把 r = dGr + (1-d)/n 整理成线性方程组
#   (E - dG) r = (1-d)/n · 1，再补上「和为 1」去掉一个相关方程
LHS <- rbind(diag(n) - damp * G, rep(1, n))
RHS <- c(rep((1 - damp) / n, n), 1)
r_lin <- particular_solution(LHS, RHS)
say("用线性方程组求解得到：", paste(round(r_lin, 6), collapse = ", "))
say("与幂迭代结果的差：", format(max(abs(r_lin - r)), scientific = TRUE))
say("排名（从高到低）：", paste(order(r, decreasing = TRUE), collapse = " > "))

# ===========================================================================
# 五、实验报告要求
# ===========================================================================
# 本实验是综合性的，报告里请回答：
#   1. 图像压缩中，误差随 k 怎么变？为什么不是线性的？
#      压缩比与视觉质量之间的平衡点在哪？（自己定一个「可接受误差」，算 k）
#   2. 三种最小二乘解法在良态与病态两种情形下分别差多少？
#      为什么工程上默认用 QR 或 SVD，而不是正规方程？
#   3. PCA 的「解释方差比」是什么意思？如果只保留一个主成分，
#      信息损失是多少？把数据换成三维，重做一遍。
#   4. PageRank 的收敛速度由什么决定？（提示：第二大特征值的模。）
#      把 d 从 0.85 改成 0.5 和 0.99，观察迭代次数。
#
# 附加题（选做）
#   A. 用真实的灰度图（如 png 读入后转矩阵）重做 7.1，k 取多少时肉眼看不出差别？
#   B. 把 PCA 用到 7.1 的图像矩阵上（对列做中心化），比较与 SVD 压缩的异同。
#   C. 给 PageRank 加上「主题偏向」（personalization vector），观察排名变化。

say("\nlab07 结束。")
