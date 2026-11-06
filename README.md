# ROpenEigen

> ### 平台状态 / Platform status
>
> **目前仅在 Windows 下完成构建与测试**（R 4.6.0 + Rtools，C++17）：
> `R CMD check` 结果为 0 error / 0 warning，包内回归测试与三套配套实验全部通过。
> **Linux 与 macOS 尚未测试**。`src/Makevars` 只用 `-I../inst/include` 与
> `CXX_STD = CXX17`，不含 Windows 专有设置，理论上可以移植，但暂不作保证。
> 若你在其他平台构建成功或遇到问题，欢迎到
> [Issues](https://github.com/Judy-ZQ/ROpenEigen/issues) 反馈。
>
> *Currently tested on Windows only (R 4.6.0 + Rtools, C++17); Linux and macOS builds are untested.*

R bindings for [Eigen 5.0.1](https://eigen.tuxfamily.org/). Eigen itself ships inside the
package (`inst/include`), so nothing external needs to be installed — but the package must
be built from source, because the C++ layer is where all the work happens.

The design follows the same lines as [ROpenCV](https://github.com/Judy-ZQ/ROpenCV): a
matrix lives on the C++ side as a reference object, base R idioms dispatch on it directly,
and the crossing between R and C++ copies once in the direction you asked for.

## Installation

ROpenEigen is a source package. Eigen 5.0.1 ships in `inst/include`, so no system Eigen
is needed — but a **C++17 toolchain** is required (Rtools on Windows, Xcode command line
tools on macOS, a recent g++/clang on Linux).

```r
# 推荐：直接从 GitHub 安装
install.packages("remotes")   # 若尚未安装
remotes::install_github("Judy-ZQ/ROpenEigen")

# 或者下载源码 tarball 后本地构建
# install.packages("ROpenEigen_0.2.0.tar.gz", repos = NULL, type = "source")
```

Dependencies are `Rcpp`, `methods` and `stats`; `MASS`, `Matrix` and `ROpenCV` are listed
under `Suggests` and are only needed by parts of the lab course, not by the package core.

Build time is a few minutes. Eigen is template-heavy, and four translation units each
instantiate the same templates twice (real and complex).

## Tutorial: the official Eigen examples, in R

`inst/examples/` holds an 11-chapter tutorial that rewrites every `.cpp` in Eigen's own
`EigenExample/` tree — 70 examples, 297 snippets, 3 special cases and 7 demo programs —
as runnable R. Each chapter's header comment names the exact files it covers, and every
number printed in the output comes from actually compiling and running the C++ side, not
from the documentation.

```r
# install first, then: eleven chapters in under a second each
source(system.file("examples", "run-all.R", package = "ROpenEigen"))
```

| ch | topic | files |
|----|-------|------:|
| 01 | quick start, comma initialisation, fixed vs dynamic size | 14 |
| 02 | sizes, storage order, `Map`, reshape, cast | 66 |
| 03 | arithmetic and the whole Cwise family | 85 |
| 04 | blocks and slices | 75 |
| 05 | reductions, visitors, broadcasting | 37 |
| 06 | solving linear systems | 23 |
| 07 | LU / QR / Cholesky / SVD / Schur / QZ | 50 |
| 08 | sparse matrices and iterative solvers | 5 |
| 09 | quaternions, angle-axis, transforms, `umeyama` | 6 |
| 10 | `SpecialFunctions` (`erf`, `lgamma`, `igamma`, …) | 3 |
| 11 | what has no R counterpart, and why | 19 |

The point is not transliteration. Chapters 01–10 give the R that works; chapter 11 goes
through expression templates, `noalias`, `Ref<>`, CRTP, compile-time sizes and `Map` and
explains for each why R either does not need it or cannot express it. See
[`inst/examples/README.md`](inst/examples/README.md) for the full index and an
Eigen-to-R cheat sheet.

## A linear algebra lab course

`inst/labs/` is a companion lab course for the textbook *Linear Algebra* (3rd edition,
Duan Zhengmin et al., Higher Education Press, 2019). Eight labs follow its six chapters,
and every starred "applications" section is covered — cross and mixed products,
input–output analysis, orthogonal polynomials, a canteen Markov chain, difference
equations, quadric surfaces — plus a setup lab and a capstone on SVD compression,
least squares, PCA and PageRank.

Each lab does the same thing twice: the hand computation from the book, written out
step by step in R, and then the same quantity from a ROpenEigen routine, with an
assertion that the two agree. In between sits `helpers.R`, which implements the
algorithms the book asks you to carry out by hand and a numerical library is right
not to provide: cofactor expansion, row-reduced echelon form, a nullspace basis,
Faddeev–LeVerrier, Gram–Schmidt (plain and metric-weighted), Lagrange completion of
squares, leading principal minors, inertia.

```bash
Rscript "$(Rscript -e 'cat(system.file("labs","run-all.R",package="ROpenEigen"))')"
```

| lab | topic | textbook |
|-----|-------|----------|
| 00 | getting started: objects, `byrow = TRUE`, numerical zero | — |
| 01 | determinants: permutations, cofactor expansion, geometry | ch. 1, p. 1–24 |
| 02 | matrices: rank by minors, inverses, Cramer, input–output | ch. 2, p. 25–56 |
| 03 | dependence, bases, Gram–Schmidt, change of basis | ch. 3, p. 57–96 |
| 04 | linear systems, least squares, Markov chain, balancing | ch. 4, p. 97–122 |
| 05 | eigenvalues, diagonalisation, difference equations | ch. 5, p. 123–148 |
| 06 | quadratic forms, inertia, definiteness, quadric surfaces | ch. 6, p. 149–166 |
| 07 | capstone: SVD compression, least squares, PCA, PageRank | — |

The labs are also where this package's sharpest edges are documented, because that is
where they hurt: `solve()` on a singular `eigenMat` returns `NaN`/`Inf` instead of
raising an error (hence `safe_solve()` / `safe_inverse()`, which check `rcond` first);
`solve()` on an overdetermined system silently returns the least-squares solution;
`eigen_eigen()` is a `SelfAdjointEigenSolver` that reads only the lower triangle, so a
non-symmetric matrix needs `eigen_geigen()` (a rotation by 90° gives `±i` there, and
nonsense from `eigen_eigen()`). See [`inst/labs/README.md`](inst/labs/README.md) for the
chapter-to-lab map, the per-lab learning outcomes and the suggested 27-hour schedule.

## Seeing the matrices: an image lab course

Numbers prove a theorem; pictures make you believe it. `inst/labs-image/` is a second
pass over the same textbook, one experiment per chapter, where every quantity is a
matrix acting on an image instead of a number on a page. Determinants become areas and
collapses, multiplication order becomes two different warps of the same picture, rank
becomes "how many layers you keep", an ill-conditioned system becomes noise blowing up
in a deblurred image, and eigenvectors become the direction a bright blob is stretched
along.

The image chain is deliberately split between the two packages: ROpenCV turns files
into matrices and matrices back into files (read, resize, colour convert, warp, blur),
ROpenEigen does the linear algebra in between. Warping is also written by hand, once
forward (`x' = Ax`, which still works when `det A = 0`) and once backward (`A⁻¹`, which
does not), because that asymmetry *is* the lesson about invertibility.

```bash
Rscript "$(Rscript -e 'cat(system.file("labs-image","run-all.R",package="ROpenEigen"))')"
```

| experiment | file | textbook | one line |
|------------|------|----------|----------|
| 0 | `img00-matrix.R` | all | an image is three matrices; grey is one linear combination |
| 1 | `img01-determinant.R` | ch. 1 | det is area, the sign is a mirror, zero is a collapse |
| 2 | `img02-operations.R` | ch. 2 | multiplication is composition; elementary matrices are row ops |
| 3 | `img03-vector-space.R` | ch. 3 | an image is a 76800-dimensional vector; projection is the best approximation |
| 4 | `img04-systems.R` | ch. 4 | blurring is `Ax = b`; conditioning amplifies noise |
| 5 | `img05-eigen.R` | ch. 5 | eigen-directions do not turn; `A^k` forgets everything else |
| 6 | `img06-quadratic.R` | ch. 6 | positive definite is an ellipse; the semiaxis is `1/sqrt(λ)` |
| 7 | `img07-svd.R` | ch. 5–6 | SVD splits an image into layers; rank is how many you keep |

The test image is drawn on the fly with R's graphics device and then read back with
`cv_imread()`, so there are no data files to ship and every printed number is
reproducible. `LINALG_IMG=D:/photo.jpg` swaps in your own. ROpenCV is a `Suggests`;
`tests/labs-image.R` skips when it is absent. See
[`inst/labs-image/README.md`](inst/labs-image/README.md) for the per-experiment numbers
and the list of pitfalls this course ran into (CV_64F inputs rejected by `cvtColor`,
`cv_GaussianBlur()` having no destination argument, `eigen_svd()` returning methods
rather than fields, `kappa()` not `cond()`).



## Quick start

```r
library(ROpenEigen)

A <- eigenMat(matrix(c(4, 1, 1, 3), 2))
A
#> <eigenMat> 2 x 2 double
#>      [,1] [,2]
#> [1,]    4    1
#> [2,]    1    3

A %*% A                 # matrix product
t(A)                    # transpose
solve(A, c(1, 2))       # solve, returns a vector like base R
eigen_det(A)            # 11
eigen_norm(A, "2")      # spectral norm
A[1, 2]                 # 1
A + 1                   # scalar broadcast
abs(-A)                 # Math group
sum(A)                  # Summary group
as.matrix(A)            # back to base R
```

The point of `eigenMat` is not that `as.matrix()` is slow — it is that everything between
the two conversions stays inside Eigen.

## Decompositions

A decomposition is computed once and reused. This is the main thing `solve(A, b)` cannot
do: base R refactorises on every call.

```r
S   <- eigenMat(crossprod(matrix(rnorm(300), 100, 3)))   # 3 x 3 SPD
b1  <- rnorm(3); b2 <- rnorm(3)
llt <- eigen_llt(S)

llt$solve(b1)     # reuses the factorisation
llt$solve(b2)     # and again
solve(llt, b1)    # same thing through the base generic
names(llt)        # "solve" "matrixL" "matrixU" "reconstructed" "logdet" "rcond" "is_ok"
```

| constructor | factorisation | notes |
|---|---|---|
| `eigen_llt()` | `A = L Lᵀ` | symmetric positive definite; reads the lower triangle only |
| `eigen_ldlt()` | `A = Pᵀ L D Lᵀ P` | symmetric, not necessarily definite |
| `eigen_lu()` | partial-pivot LU | square |
| `eigen_full_lu()` | complete-pivot LU | gives `rank()`, `rcond()`, `is_invertible()`; also accepts rectangular — but `det()` and `inverse()` stay square-only and raise an error otherwise |
| `eigen_qr()` | Householder QR | least squares |
| `eigen_colpiv_qr()` | column-pivoting QR | rank-revealing; use when the matrix may be rank deficient |
| `eigen_svd()` | SVD (thin U, V) | `singular_values()`, `u()`, `v()`, `pseudo_inverse()`; `method = "auto"/"jacobi"/"bdc"` picks the backend, `$backend()` says which ran |
| `eigen_eigen()` | self-adjoint eigensolver | symmetric/Hermitian; reads the lower triangle only |
| `eigen_geigen()` | complex eigensolver | any square matrix; always returns complex eigenvalues |

One-shot shorthands exist for the common cases: `eigen_values()`, `eigen_vectors()`,
`eigen_singular_values()`, `eigen_pinv()`, `eigen_chol()`, `eigen_det()`, `eigen_inverse()`.

## Sparse matrices

`eigen_sparse()` wraps `Eigen::SparseMatrix<double>`, column-major — the same
layout as a `dgCMatrix`. Indices are **1-based**, as `Matrix::sparseMatrix()`
takes them, and duplicate `(i, j)` pairs are summed the same way.

```r
S <- eigen_sparse(c(1, 2, 3), c(1, 2, 3), c(4, 5, 6), nrow = 3, ncol = 3)
b <- c(1, 2, 3)
dim(S); eigen_sparse_nnz(S)     # 3 3 ; 3
as.matrix(S)                    # back to dense
t(S); S + S; S * 2              # stay sparse
S %*% b                         # sparse x dense is dense, as it must be
solve(S, b, method = "llt")     # factorises sparsely; returns a plain vector
```

`dim()`, `length()`, `t()`, `%*%`, `solve()`, `as.matrix()`, `as.vector()`,
`mean()`, `is.na()`, `anyNA()`, `as.data.frame()`, the arithmetic group, the
`Summary` group and the `Math` group all dispatch on it. (`as.numeric()` does
not — it is a primitive, so it never sees the class.) Conversion works both ways
from dense matrices, `eigenMat`, and (if Matrix is installed) `sparseMatrix`
objects.

Eleven solvers, all reachable through `method` -- four direct, seven iterative:

| method | Eigen class | when |
|---|---|---|
| `"lu"` | `SparseLU` | general square |
| `"llt"` | `SimplicialLLT` | symmetric positive definite |
| `"ldlt"` | `SimplicialLDLT` | symmetric indefinite |
| `"qr"` | `SparseQR` | rectangular / least squares |
| `"cg"` | `ConjugateGradient` | large SPD, iterative |
| `"bicgstab"` | `BiCGSTAB` | **non-symmetric** square, iterative |
| `"lsqr"` | `LeastSquaresConjugateGradient` | rectangular, iterative |
| `"gmres"` | `GMRES` | general non-symmetric, iterative; residual decreases monotonically |
| `"minres"` | `MINRES` | **symmetric** (definite or not), iterative |
| `"dgmres"` | `DGMRES` | `GMRES` plus diagonal scaling, for rows of wildly different scale |
| `"idrs"` | `IDR(s)`, s = 4 | general, iterative; fixed memory, often converges where `"bicgstab"` does not |

`"bicgstab"` is the one that matters most in practice: `ConjugateGradient`
requires symmetry, and most real sparse systems are not symmetric.

`"cg"` and `"bicgstab"` take a `precond` argument, and the defaults differ by
method for a mathematical reason rather than a taste one.

For `"bicgstab"` the default is `"ilut"` (`IncompleteLUT`), and that choice is
measured rather than assumed. Iterations on a 2-d convection–diffusion grid
(`bicgstab`, `tol = 1e-10`):

| grid | convection | `ilut` | `diagonal` |
|---|---|---|---|
| 32×32 | −1.1 | 28 | 464 |
| 32×32 | −1.4 | 85 | does not converge |
| 48×48 | −1.1 | 99 | 1352 |
| 48×48 | −1.2 | 1925 | does not converge |

`IncompleteLUT` performs an incomplete LU **without pivoting**, so it can hit a
zero pivot and produce `NaN` — on a strongly non-symmetric 63×63 grid it does.
When that happens the solve stops with a non-convergence error rather than
handing back the `NaN`. `precond = "diagonal"` is the fallback: weaker, but it
never breaks down.

`"gmres"`, `"dgmres"` and `"idrs"` also default to `"ilut"`.

For `"cg"` and `"minres"` the default is `"diagonal"`, and `precond = "ilut"` is
an **error** for both. ConjugateGradient requires a symmetric positive definite
preconditioner and MINRES requires a symmetric one; `IncompleteLUT` is not
symmetric, so the pairing is not "slower" — it is invalid, and it does not
converge. Rejecting it is better than offering a combination that silently
produces nothing.

One semantic that bites: MINRES reads **only the lower triangle** (it wraps the
matrix in `selfadjointView<Lower>()`). A non-symmetric matrix passed to it is
not rejected; the upper half is simply ignored. `tests/sparse.R` pins that
behaviour so nobody "fixes" it into a disagreement with Eigen.

The iterative methods report what they cost: the result carries `iterations` and
`residual` attributes, and `tol` / `maxiter` are arguments.

```r
# 1-d Laplacian, tridiagonal: the classic CG test problem
n  <- 200
L  <- eigen_sparse(c(1:n, 2:n, 1:(n - 1)), c(1:n, 1:(n - 1), 2:n),
                   c(rep(2, n), rep(-1, 2 * (n - 1))), nrow = n, ncol = n)
b  <- rnorm(n)
x  <- eigen_sparse_solve(L, b, "cg", tol = 1e-10, maxiter = 5000)
attr(x, "iterations"); attr(x, "residual")
```

This is not a marginal speedup. On a 1-d Laplacian with `n = 4000` (11,998
stored entries), against `solve(as.matrix(L), b)`:

| | time | vs dense |
|---|---|---|
| dense `solve` | 6570 ms | — |
| `"lu"` | 0.85 ms | **7729×** |
| `"llt"` | 0.40 ms | **16425×** |
| `"cg"` | 41.25 ms | 159× |

Dense LU is `O(n³)`; banded sparse LU is `O(n)`. Solution error relative to
`‖x‖ ≈ 10⁷` is at machine precision.

**Every solver checks `Eigen::info()`.** Eigen's sparse factorisations do not
throw — they set `info()` to `NumericalIssue` and keep returning values, so an
unchecked `SimplicialLLT` on an indefinite matrix hands back the solution of a
*different* system. That is the same failure that shipped once on the dense
side, so all eleven methods stop instead, and `tests/sparse.R` asserts it.

What no method checks -- because Eigen does not check it either -- is symmetry.
`"llt"`, `"ldlt"` and `"minres"` read only the lower triangle, so a
non-symmetric matrix is silently solved as its symmetric part; `"llt"` and
`"ldlt"` still stop on a bad pivot, which is how an indefinite matrix announces
itself. If you are unsure whether the matrix is symmetric, ask
`isSymmetric()` first, or use `"lu"` / `"gmres"`, which do not care.

Complex sparse matrices are deliberately unsupported: passing complex input
stops rather than silently discarding the imaginary part.

### Reductions do not densify

`sum(S)`, `prod(S)`, `mean(S)`, `min(S)`, `max(S)` work through the `Summary`
group, and the `eigen_*` one-shot functions accept sparse input too:

`eigen_sum()`, `eigen_prod()`, `eigen_mean()`, `eigen_sumsq()`,
`eigen_norm()`, `eigen_trace()`, `eigen_diag()`, `eigen_adjoint()`,
`eigen_is_complex()`, `eigen_all_finite()`, `eigen_has_nan()`, `eigen_has_inf()`,
`eigen_which_min()`, `eigen_which_max()`,
`median()`, `quantile()`, `isSymmetric()`, `c()`, `crossprod()`,
`tcrossprod()`, `aperm()`, `is.finite()`.

All of them come out of one `O(nnz)` pass over the stored entries, and implicit
zeros take part — for an all-positive matrix `prod()` is `0` and `min()` is `0`.

Missing values are a place where that shared pass had the wrong contract for its
callers. `reigen_sparse_stats()` skips `NA` and `NaN` while accumulating, because
that is what `na.rm = TRUE` needs, and it counts them separately so the caller
can decide. The callers did not decide. The result was not an error but a
smaller, entirely plausible number:

```
A = matrix(c(1, 2, NA, 4, 5, 6), 2, 3)

eigen_sum(S)   -> 18      sum(S)  -> NA
eigen_prod(S)  -> 240     prod(S) -> NA
mean(S)        -> 3       mean(A) -> NA
```

Same quantity, two doors, different answers — the exact failure mode
`tests/consistency.R` exists to catch, except nothing there had thought to put an
`NA` in a sparse matrix. The fix is one guard, `.sparse_na()`, applied by every
caller that reads a statistic; `NA` wins over `NaN` when both are present, which
is what `sum(c(NA, NaN))` does.

`eigen_trace()` needed a different guard, and this is the part that is easy to
get wrong: the trace only touches the main diagonal, so an `NA` **off** the
diagonal must not affect it. The matrix-wide counts cannot answer "is there a
missing value on this diagonal", so the statistics pass now counts `trace_na`
and `trace_nan` alongside the trace itself.

`eigen_norm()` covers all five types without densifying: `"F"` and `"M"` are one
scan, `"O"` and `"I"` are `O(nnz)` bin sums over the stored entries, and `"2"` is
a power iteration on `A'A`. Only `"p"` has no sparse meaning and stops.

The Frobenius norm needed a scaling pass, on both sides. `sqrt(sum(x^2))` has two
failure modes at extreme magnitudes, and neither of them raises anything:

```
matrix of 1e-200   -> 0      (the squares underflow; the answer is 1.414e-200)
matrix of 1e300    -> Inf    (the squares overflow;  the answer is 1.414e+300)
```

Returning `0` is the worse of the two, because `0` looks like a perfectly ordinary
answer. `base::norm(x, "F")` does not have this problem — it goes through LAPACK's
`dlange`, which scales by the largest element first. `eigen_norm()` now does the
same: track `max|x|`, accumulate `(x/max)^2`, multiply back. `O(n)`, no temporary
matrix.

`Inf` and `NaN` are handled outside that loop, because the scaling algorithm
computes `Inf/Inf` as soon as the maximum is infinite. Base gives `Inf` for a
matrix containing `Inf` and `NaN` for one containing `NaN`, and so does this.

`eigen_sumsq()` deliberately does **not** get the treatment. It corresponds to
`sum(x^2)`, and `sum(x^2)` really does overflow — "fixing" it would make
`eigen_sumsq()` disagree with base R. Only the norm is scaled.

The default is `"O"`, not `"2"`. `type`'s default is the whole vector
`c("O", "I", "F", "M", "2")` and the function takes its first element, which is
what `base::norm()`'s `match.arg` does too. Worth stating because I got it wrong
once: I wrote in this file that the default was `"2"` and added a test labelled
`eigen_norm(S) 默认 2 报错`. The assertion passed — but only because the old
code rejected *every* type except `"F"`, so the label was never exercised. A
wrong label does not fail; it just misleads the next reader. There is now an
assertion that the default equals `norm(A, "O")`.

The spectral norm is the interesting one, because it is the only type that
genuinely cannot be read off in one pass: it is the square root of the largest
eigenvalue of `A'A`, and Eigen's `SelfAdjointEigenSolver` takes dense input
only. The iteration is `z <- A'(A v)`, `v <- z/|z|`, using `|z|` as the running
estimate of that eigenvalue — that estimate converges at the *square* of the
vector's error rate, which is why the test tolerance can be `1e-8` against
`norm(A, "2")` and the iteration still stops quickly.

Two things about the starting vector. It is a fixed-seed pseudo-random vector,
not all-ones, and it has to be. An all-ones start has a real failure mode: if
the leading singular vector is orthogonal to the ones vector — any matrix whose
columns each sum to zero does it — the first projection is zero and the
iteration sits there forever, returning `0`. A zero spectral norm is a
perfectly plausible-looking answer, which is the worst kind. The fixed seed
also means the same matrix always gives the same digits, so the assertion can
be tight.

If it has not converged in 10000 steps it stops rather than returning a
number that looks converged and is not. A warning would be worse: the contract
of a norm is "a value or an error", and a third outcome invites the value to be
used anyway.

`median()` and `quantile()` are `O(nnz log nnz)`: sort the stored values, then
count how many zeros are implied. One trap is worth writing down, because I fell
into it. **Zeros do not sort to the front.** The ascending order is

```
negatives ... | zeros ... zeros | positives ...
```

so the k-th element has to be looked up in whichever of the three runs it lands
in. A first version assumed the zeros came first and returned `0` for the 0%
quantile of a matrix whose minimum was negative. The unit tests missed it
because their fixture was `abs(A)` — all non-negative, where the wrong
assumption happens to give the right answer. The sparse idiom probe caught it
instead, and the tests now use a fixture with positives, negatives and exact
zeros.

A second bug in the same helper is worth recording because of how small it was.
The interpolation is `kth(q) + r * (kth(q + 1) - kth(q))`. At `p = 1` the
fraction `r` is `0`, so the second term looks like it cannot matter — but the
index `q + 1` is still *evaluated*, and for an all-non-positive matrix it runs
past the empty positives run and yields `NA`. And `0 * NA` is `NA`. So the 100%
quantile came back missing. The fix is to make the lookup total for every `k`:
clamp the index and return `0` when the run is empty. The lesson generalises —
an expression multiplied by zero is not free, and the "obviously harmless" term
is where the off-by-one hides.

`eigen_which_min()` and `eigen_which_max()` used to stop on sparse input, on the
grounds that the extreme value is usually an implicit zero, so its *position*
depends on the storage pattern rather than on the values. That is true, and it
is also the wrong conclusion: the position is well defined — the first
column-major cell whose value is `0`, stored or not — and it is exactly what
`which.min(as.matrix(x))` returns. Refusing just forced the user to densify, at
the cost the sparse representation exists to avoid. They now return it in
`O(nnz + ncol)`.

The interesting part is what looking for that behaviour turned up on the
*dense* side. `eigen_which_min()` disagreed with `which.min()` — but only in one
position out of nine. Eigen's `minCoeff(&r, &c)` takes the first element as its
initial best, and its comparator `numext::mini(b, a)` is `(b < a) ? b : a`; when
that initial `a` is `NaN`, every later `b < NaN` is false, so the `NaN` is kept
for the whole scan. On a 3×3 matrix of ones with a single `NA` moved through all
nine cells, eight positions gave the right answer and only `NA` in the first
column-major cell returned `NA@(1,1)` — while `which.min(c(NA, 1))` is `2`.

A bug that depends on *where* the missing value sits is invisible to a test that
only ever puts it in one place, which is what the old test did. The fix is an
explicit column-major scan that skips `NA`/`NaN` and replaces only on a strict
comparison, and `tests/compat.R` now loops the `NA` through all nine positions.
`Inf` is not skipped — `which.min(c(Inf, 1))` is `2` and `which.min(c(-Inf, 1))`
is `1`, so the test is `ISNAN()`, not `!R_FINITE()`.

Everything else that needs a dense matrix says so, and says what to do:

```r
S <- eigen_sparse(c(1, 2, 3), c(1, 2, 3), c(4, 5, 6), nrow = 3, ncol = 3)
tryCatch(eigen_svd(S), error = function(e) conditionMessage(e))
#> x is an eigenSparse, but this function needs a dense matrix. Use
#> as.matrix(x) to densify, or a sparse-aware entry point (eigen_sum, ...)
```

That message matters because the old one was `x: expected a numeric matrix or
vector, got externalptr` — which neither named the type nor suggested a way
forward. Thirty-five entry points share it, so fixing it in one place fixed all
of them.

## Geometry

Eigen's `Geometry` module is wrapped as reference objects, one class per
concept: `eigenQuat`, `eigenAngleAxis`, `eigenTransform` (3-D affine) and
`eigenRotation2D`. They are not matrices with a `dim` attribute — a quaternion
is not a 4×1 column — so `as.matrix()` is how you get numbers out, and `names()`
lists the `$` methods.

```r
q <- eigen_quat_from_angle_axis(pi / 2, c(0, 0, 1))
q$matrix()                      # 3x3 rotation, z by 90 degrees
eigen_quat_rotate(q, c(1, 0, 0))    # (0, 1, 0)
q$slerp(eigen_quat(), 0.5)          # halfway back to identity
q$toAngleAxis()$angle()             # 1.570796

t <- eigen_translation(c(1, 2, 3))
t$points(rbind(c(0, 0, 0)))         # 1 2 3
eigen_transform_inverse(t)          # -1 -2 -3
eigen_umeyama(src, dst)             # least-squares rigid alignment
```

`q1 %*% q2` composes rotations in the same order as matrix multiplication, and
`t1 %*% t2` does the same for transforms. Every `$` method returns a **new**
object: `t$translate(v)` does not modify `t`, because a shared mutable transform
is a trap in a language where everything else has value semantics.

The one semantic worth stating out loud: `t$rotate(q)$translate(v)` translates in
the **rotated** frame, so the stored translation is `R·v`, not `v`. That is what
`Eigen::Transform::translate()` does, and matching it matters more than
matching intuition.

## Special functions

`unsupported/Eigen/SpecialFunctions` is exposed elementwise, one function per
name, and checked against R's own `stats` implementations in `tests/special.R`:

| `eigen_` | Eigen | R equivalent |
|---|---|---|
| `erf`, `erfc` | `numext::erf` | `2 * pnorm(x * sqrt(2)) - 1` |
| `gammaln` | `numext::lgamma` | `lgamma` |
| `digamma`, `polygamma(n, x)` | `numext::digamma` | `digamma`, `psigamma` |
| `ndtri` | `numext::ndtri` | `qnorm` |
| `igamma(a, x)`, `igammac` | `numext::igamma` | `pgamma(x, a)` and its upper tail |
| `betainc(a, b, x)` | `numext::betainc` | `pbeta(x, a, b)` |
| `zeta(x, q)` | `numext::zeta` | none in base R |
| `bessel(x, type)` | `bessel_i0` … `bessel_y1` | `besselI`, `besselJ`, `besselK`, `besselY` |

`eigen_zeta()` **does not** analytically continue: Eigen implements the Hurwitz
series, so `x <= 1` returns `NaN` where R's `zeta()` would return `-1/2` at 0.
That gap is asserted in the tests rather than papered over.

## Matrix functions

`expm`, `logm`, `sqrtm` and `powm` come from `unsupported/Eigen/MatrixFunctions`
and are **not** elementwise — `eigen_expm(A)` is `sum A^k / k!`. The elementwise
versions are just `exp(A)` on an `eigenMat`, which the `Math` group already
provides.

These four do not trust Eigen's answer; they verify it. Eigen's real matrix
logarithm on a matrix with a negative eigenvalue returned `diag(0, 0.693)` for
`diag(-1, 2)` — every entry finite, every entry wrong. `allFinite()` cannot see
that. So `eigen_logm()` checks `expm(logm(A)) == A`, `eigen_sqrtm()` checks
`sqrtm(A)^2 == A`, integer powers are checked against repeated multiplication,
and non-integer powers are checked against `expm(p * logm(A))` after `logm(A)`
has itself passed its check. If a check fails, the call stops and says which
eigenvalue is the likely cause, instead of returning a plausible-looking matrix.

## What dispatches, and what cannot

`eigenMat` hooks into these base generics:

`Ops`, `Math`, `Summary`, `%*%`, `[`, `[<-`, `c`, `dim`, `dim<-`, `length`, `t`,
`solve`, `determinant`, `mean`, `names`, `print`, `is.na`, `anyNA`, `is.matrix`,
`is.array`, `is.numeric`, `as.matrix`, `as.vector`, `as.data.frame`, `Re`, `Im`,
`Mod`, `Arg`, `Conj`, `sort`, `rev`, `unique`, `duplicated`, `median`, `quantile`,
`diff`, `isSymmetric`, `rbind`, `cbind`.

`eigenSparse` hooks into the same set minus the parts that only make sense once
the matrix has been materialised: `Ops`, `Math`, `Summary`, `%*%`, `[`, `dim`,
`length`, `t`, `solve`, `mean`, `is.na`, `anyNA`, `is.matrix`, `is.array`,
`is.numeric`, `as.matrix`, `as.vector`, `as.data.frame`, `print`. The `Summary`
and `Math` groups matter more than they look: without them `sum(S)` reports
`invalid 'type' (externalptr) of argument`, because the sparse class vector
carries no `eigenMat`, so dispatch never reaches the dense methods. They are
computed from the stored entries in `O(nnz)` and **implicit zeros take part** —
for an all-positive sparse matrix `min()` is `0` and `prod()` is `0`.

Sparse and dense objects can also be mixed: `As + E` returns dense. That needs a
trick, because R's `Ops` group generic warns `Incompatible methods` and falls
back to the default whenever the two operands resolve to *different* methods.
`As + E` therefore errored with `non-numeric argument to binary operator` and
neither method body was ever entered. `Ops.eigenSparse` is now the *same
closure* as `Ops.eigenMat`, with the sparse cases handled by an internal
`.ops_sparse()`; R only calls a method when both sides agree on which one it is.

The same requirement applies to `%*%`, `crossprod()` and `tcrossprod()`, which
all dispatch on **both** operands. Those three were missed for a while, and the
way it surfaced is worth recording: every combination tested was sparse×sparse
or sparse×base-matrix, which works, because only one operand carries a class.
The broken case was `S %*% x` right after `eigen_sparse_solve()` — and that is
the most natural next line a user writes, because the solver returns an
`eigenMat`:

```r
L <- eigen_sparse(c(1, 2, 3), c(1, 2, 3), c(4, 5, 6), nrow = 3, ncol = 3)
b <- c(1, 2, 3)
x <- eigen_sparse_solve(L, b, "lu")
L %*% x          # residual check -- this used to fail
```

`%*%.eigenSparse`, `crossprod.eigenSparse` and `tcrossprod.eigenSparse` are now
the same closures as their dense counterparts, with the four shape combinations
in one internal `.matmul()`.

`Summary` behaves differently again: it dispatches on the **first argument
only**, so `sum(S, E)` and `sum(E, S)` run different methods. Both now work, but
`sum(1, E)` cannot — the first argument is a plain double, so no method is
consulted. Put the Eigen object first.

### Functions that are not generics but work anyway

Several base functions that are not S3 generics still behave correctly, because
their bodies call `as.matrix()` — which *is* a generic and does dispatch. Nobody
designed this; it falls out of the implementation, which is exactly why it needs
to be written down and tested:

`qr()`, `svd()`, `eigen()`, `kappa()`, `rcond()`, `det()`, `diag()`, `max.col()`,
`apply()`, `scale()`, `sweep()`, `prop.table()`, `head()`, `tail()`, `upper.tri()`,
`lower.tri()`, `kronecker()`, `outer()`, `range()`.

### Internal generics: dispatch you cannot see

A second, subtler group dispatches through `.Internal()` with no `UseMethod` in
the body and no membership in `Ops`/`Summary`/`Math`. Reading the source says
"not a generic"; the interpreter disagrees. `crossprod()`, `tcrossprod()`,
`aperm()`, `is.finite()`, `is.infinite()` and `is.nan()` are all in this group,
and all six have methods here.

This was originally a **wrong entry** in the table below — they were listed as
unhookable on the grounds that `body(crossprod)` contains no `UseMethod`. The
tell is that `Math.eigenMat` never received `is.finite(E)` at all: the call died
with `default method not implemented for type 'externalptr'`, which is not
something `Math.eigenMat` can produce. Probing dispatch directly
(`scratch/probe_dispatch2.R`) settles it in seconds; reading source does not.

Two details worth knowing:

* `crossprod()`/`tcrossprod()` treat a dimensionless vector as a **column**,
  while `%*%` treats a left-hand vector as a **row**. So they are not
  `t(x) %*% y` shortcuts — for vector input that spelling gives the transpose.
* They dispatch on **both** operands, so `crossprod(m, E)` works even though `m`
  is a plain matrix. `aperm()` only accepts `c(2,1)` and `c(1,2)`, because a
  two-dimensional object has no other permutation.

`is.finite()`, `is.infinite()` and `is.nan()` return a **logical base matrix**,
not an `eigenMat` — `eigenMat` stores `double` only. `is.na()` behaves the same
way.

`diag()` deserves a note: it takes the diagonal with a **linear** index
(`x[1 + 0:(m-1)*(nrow+1)]`), not `x[cbind(i, i)]`, so it only works if `[`
implements linear subscripting the way base R does. Ours treated a single index
as a *row* index, so `diag(E)` failed with "subscript out of bounds" — and
`E[1]` quietly returned the whole first row instead of the first element. The
hand-written tests only ever used two-index forms, so nothing caught it until
this audit.

### Functions that cannot be hooked

These are plain functions that call `.Internal()` on their argument. There is no
generic to hook, so they error, and `tests/compat.R` asserts that they keep
erroring rather than silently returning something else:

| cannot work | use instead |
|---|---|
| `rowSums()`, `colSums()`, `rowMeans()`, `colMeans()` | `apply(E, 1, sum)`, `eigen_sum()` for the total |
| `norm()` | `eigen_norm()` |
| `which.min()`, `which.max()` | `eigen_which_min()`, `eigen_which_max()` |
| `sd()`, `var()` | `as.matrix()` first |
| `sum(1, E)`, `max(1, E)` | put the Eigen object first: `sum(E, 1)` |
| `det()` | works, but `eigen_det()` avoids the exponentiation |
| `diag()` | works (see above); `eigen_diag()` avoids the copy |

`det()` and `diag()` are listed in both tables on purpose: they work, but for a
reason that is invisible at the call site, and each has a direct equivalent that
is both faster and free of the round-trip through a base matrix.

`as.numeric()` and `is.complex()` are not generics either; use `as.vector()` and
`eigen_is_complex()`. `c()` only dispatches on its **first** argument, so
`c(E, 1)` flattens while `c(1, E)` returns a list; `unlist()` is not a generic at
all.

Three primitives report the *storage* type rather than the logical one, and no
method can change that: `typeof(E)` is `"externalptr"`, `storage.mode(E)` is
`"externalptr"` and `is.double(E)` is `FALSE`. Test with `is.numeric(E)` (which
is `TRUE`) or `is_eigen(E)`. This is inherent to representing the matrix as an
external pointer, not an oversight.

**Vector recycling is not supported.** `E * c(2, 3)` stops with
"non-conformable elementwise operands: 2x2 vs 2x1", where base R would recycle
the vector down the columns. A dimension-less operand is treated as a column
vector, and Eigen does not broadcast beyond scalars. This is deliberate — silent
recycling is a frequent source of shape bugs — but it is a real difference from
base R, so write `E * eigenMat(matrix(c(2, 3)))` (or a scalar) if you meant it.

## Design notes

**Reference objects, not S4.** Matrices are `EXTPTRSXP` with class
`c("eigenMat", "eigenObject")`, released by the `XPtr` finaliser. Decompositions use the
same mechanism with their own class.

**One copy at the boundary.** `eigenMat(x)` copies R's matrix into an `Eigen::MatrixXd`;
`as.matrix()` copies back. Intermediate results never touch R. This keeps R's
copy-on-modify semantics intact — the alternative, mapping R's memory directly, makes
`A[1,1] <- 2` mutate the caller's matrix.

**Real and complex, no integers.** `Eigen::MatrixXi` and `Eigen::MatrixXd` are distinct
types that do not interoperate, so supporting integers would mean writing every operation
twice. Integer input is promoted to double, as `RcppEigen::as<MatrixXd>` does. Complex is
a first-class type throughout.

**No generated overload dispatcher.** ROpenCV needs one because OpenCV has hundreds of C++
overloads per function. Eigen's surface does not: a decomposition method has a fixed
signature, and the few genuine overloads (`solve`, `eigen_svd`) are resolved with
`inherits()`. The class registry and `$` access are kept, because they are what makes
`llt$solve(b)` readable.

**Thresholds are relative.** `JacobiSVD::threshold()` returns `ε · diagSize`, not an
absolute cutoff; `pseudo_inverse()` multiplies it by the largest singular value before
comparing. Getting this wrong turns a zero singular value into `1/0` and the pseudo-inverse
into `1e14`.

**Eigen's `eigen_assert` is compiled out.** R builds packages with `-DNDEBUG`, so every
`eigen_assert` in Eigen disappears — including the ones that guard preconditions.
`FullPivLU::determinant()` and `inverse()` are documented for square input only and assert
`rows() == cols()`; on a rectangular factorisation the release build silently reads out of
bounds and the process dies later at an unrelated allocation. The wrapper therefore checks
the precondition itself and stops with a message. Assume any Eigen method with a stated
precondition needs an explicit check here.

`rcond()` sidesteps the issue entirely: it is `σmin / σmax` from the SVD rather than
Eigen's `FullPivLU::rcond()`. That estimate is square-only *and* returns 1 for an exactly
singular matrix, where the answer is 0 — so it is not used anywhere. Every decomposition
object's `$rcond()` method routes through the same helper, including `LLT` and `LDLT`,
which do ship their own `rcond()` in Eigen but compute a different quantity. One matrix
must not produce two different condition numbers depending on which entry point was
called; `eigen_rcond(x)`, `eigen_llt(x)$rcond()` and `eigen_full_lu(x)$rcond()` now
agree to 1e-8.

**Non-finite input is rejected at the decomposition boundary.** With `eigen_assert`
gone, `Eigen::JacobiSVD` on a matrix containing `NA` does not complain — it returns
all-`NaN` factors. Base R's `svd()`, `qr()` and `eigen()` stop with "infinite or
missing values in 'x'". The nine decomposition constructors therefore call
`from_sexp_finite()` and stop with a message; elementwise operations and
reductions are deliberately *not* gated, because there `NA` should propagate.

**`det()` is not generic, but it dispatches anyway — through `determinant()`.**
Base R defines `det <- function(x, ...) c(z$sign * exp(z$modulus))` where
`z <- determinant(x, ...)`. `determinant()` *is* a generic, so `det(x)` on an
`eigenMat` silently reaches `determinant.eigenMat` and exponentiates whatever
`modulus` it finds. That method must therefore return `log|det|` when
`logarithm = TRUE` — the base R convention — or `det()` comes back as
`exp(det(A))`: 403.43 instead of 6, with no error. `modulus` also carries base
R's `logarithm` attribute, and `sign` has to be `-1` for a negative determinant
rather than a constant `1`.

**Eigen reports failure through `info()`, not by throwing.** `Eigen::LLT` on an
indefinite matrix does not raise anything: it sets `info()` to `NumericalIssue`
and leaves the factors holding a *different* matrix. A 2×2 indefinite input
reconstructs to `[[1,2],[2,5]]` instead of `[[1,2],[2,1]]`, and `$solve()`
happily returns answers for that other matrix. `eigen_llt()` therefore checks
`info()` and stops, as `base::chol()` does. The same pattern applies to every
Eigen class with an `info()` — assume a silent failure unless the wrapper checks.

**Determinants are accumulated in log space.** `det` of a 50×50 random matrix is
around 1e±30; the intermediate product overflows long before it is finished.
`eigen_logdet()` and `determinant()` therefore sum `log|U_ii|` over the LU
factorisation and recover the sign from the diagonal phases times the parity of
the permutation — never forming the product. `eigen_det()` does return the raw
value and can overflow, exactly as `det()` can.

## Performance

The reference-object design exists for one reason: a factorisation should be
computed once. Measured on one machine (Windows 11, R 4.6.0, g++ 14.3), solving
the same symmetric positive-definite system repeatedly:

| n | `solve(A, b)` | `eigen_solve(E, b)` | `llt$solve(b)` | reuse vs base |
| --- | --- | --- | --- | --- |
| 100 | 0.18 ms | 0.09 ms | 0.005 ms | **36×** |
| 400 | 7.55 ms | 2.85 ms | 0.10 ms | **75×** |

So the claim holds where it matters: once the factorisation exists, each solve is
essentially free, and the gap widens with `n` because base R refactorises every
call. Even the one-shot form is 2–2.6× faster than `solve()`, since Eigen's
Cholesky avoids LAPACK's per-call overhead.

Conversion is not free but it is small and linear: at n = 400, `eigenMat(A)` costs
0.45 ms and `as.matrix(E)` 0.20 ms — 16% of a one-shot `eigen_solve`. At n = 100 it
is 39%. This is why the one-shot helpers are worth it for a single call but the
object form is worth it for more than one.

**`eigen_svd()` has two backends.** `JacobiSVD` is the accurate choice for
rank-deficient and ill-conditioned input but has a large constant: a 400×400 SVD
takes about **1.07 s**, against 1.6 ms for `eigen_llt()`. Eigen also ships
`BDCSVD`, a divide-and-conquer variant, so `eigen_svd()` takes
`method = c("auto", "jacobi", "bdc")`:

| n (square) | `jacobi` | `bdc` | speedup | base R `svd()` |
| --- | --- | --- | --- | --- |
| 50 | 1.6 ms | 0.9 ms | 1.9× | <1 ms |
| 100 | 11.3 ms | 2.9 ms | 3.9× | <1 ms |
| 200 | 91.9 ms | 11.9 ms | 7.7× | <1 ms |
| 400 | 1081 ms | 65.6 ms | **16.5×** | 100 ms |

`"auto"` (the default) uses `"bdc"` at `min(dim(x)) >= 16` and `"jacobi"` below
that, where the difference is not worth a second code path. `$backend()` reports
which one ran.

Switching backends is safe here for a specific reason: `threshold()`, `rank()` and
`nonzeroSingularValues()` are **not** implemented per-backend — they live in the
shared `SVDBase` base class (`Eigen/src/SVD/SVDBase.h`). So the rank cut-off and
the pseudo-inverse truncation point are the same code either way, and cannot drift
apart. Measured agreement on the same 400×400 input: singular values to 1e-13
relative, `threshold()` bit-identical, `rank()` identical.

For a symmetric matrix, `eigen_eigen()` is still both faster and more accurate
than an SVD.

## Tests

Eleven files run as part of `R CMD check`, 8,289 assertions in total.  `R CMD check
--no-manual` reports 0 errors and 0 warnings, plus one NOTE from `checking
compiled code`: `reigen_sparse.o` references `rand`, which is `std::rand`
inside `unsupported/Eigen/IterativeSolvers` -- Eigen's IDR(s) builds its
shadow space with it, so the symbol arrives with that solver and cannot be
removed without reimplementing IDR(s) or patching the Eigen header.

| file | assertions | what it covers |
| --- | --- | --- |
| `tests/core.R` | 131 | construction, indexing (linear, row, column, block), arithmetic, reductions, generators, complex, decompositions, printing |
| `tests/edge.R` | 125 | empty and 1×1 matrices, shape mismatches, singular and rank-deficient input, `na.rm`, extreme magnitudes (1e-300 to 1e300, where a naive sum of squares underflows to zero or overflows to `Inf`), the 2-D cross product returning the scalar z component while length 4 keeps erroring, the contract of every decomposition object |
| `tests/compat.R` | 91 | base R functions that must match exactly, and ones that must keep refusing rather than silently returning something else; the `which.min`/`which.max` section moves a single `NA` through all nine cells of a 3×3 matrix, because the bug it guards against only fired when the `NA` sat in the first column-major cell |
| `tests/decomp.R` | 153 | every decomposition's defining identity on random input of assorted shapes and conditioning, plus the two SVD backends agreeing |
| `tests/consistency.R` | 198 | ROpenEigen against *itself*: every entry point that computes the same quantity must return the same answer |
| `tests/sparse.R` | 390 | 1-based/0-based indexing, duplicate triplets, transposes, all four `%*%` combinations, all eleven solvers checked against a dense equivalent, the iterative ones also checked for convergence reporting, reductions and the math group against base R, all five norm types against `norm()`, extreme-value positions checked against `which.min`/`which.max` on fixtures where the extreme is an implicit zero, mixing with a dense `eigenMat`, and a guard that the method names in the signature still match the internal list that maps them to C++ codes |
| `tests/fuzz.R` | 4966 | 250 seeded random matrices × 6 flavours (real, complex, integer, logical with `NA`, ill-conditioned, rank-deficient) × random shapes, each checked against base R |
| `tests/fuzz_sparse.R` | 2032 | 60 seeded random sparse matrices with random shape, density and value distribution — plus 15 hand-built extremes (all zero, single non-zero, all negative, first column dense, empty leading columns, explicit zeros, …) — every entry point compared against base R on the same dense matrix, with a third of the cases carrying an injected `NA` and a fifth a `NaN` |
| `tests/geometry.R` | 72 | quaternions, axis-angle, affine transforms and `Rotation2D` against Eigen's own semantics: `(w, x, y, z)` order, `q1 %*% q2` order (checked by comparing against `R1 %*% R2`), slerp, `umeyama` recovering a known rotation and translation, and the value semantics of `$` methods |
| `tests/special.R` | 57 | every `SpecialFunctions` entry against R's `stats` counterpart, Bessel functions against `besselI`/`besselJ`/`besselK`/`besselY`, and the matrix functions verified by their defining identities rather than by tabulated values; also pins the two known Eigen limits (no analytic continuation for `zeta`, no real log for negative eigenvalues) |
| `tests/eigen_examples.R` | 74 | the official examples in `EigenExample/` redone in R: quick start, arithmetic, arrays, block operations, reductions/broadcasting, the linear-algebra tutorial (QR, LU, LDLT, SVD, rank-revealing, threshold), the sparse Poisson example and matrix-free CG -- the expected numbers are the actual output of compiling and running those `.cpp` files against the bundled headers |

`decomp.R` is the one that checks *contracts* rather than examples. Each
decomposition is defined by an identity — `A = L·Lᴴ`, `P·A = L·U`,
`A = V·Λ·Vᵀ` — and reconstruction is the strongest single check because it
exercises every factor at once. `core.R` proving that `eigen_llt()` works on one
3×3 matrix says nothing about whether it works on an 8×8 one; reconstruction on
random input does.

`consistency.R` covers the one failure mode the others structurally cannot see.
Every other file compares against base R, or a decomposition against its own
identity; both leave each entry point looking fine on its own. But three of the
bugs that actually shipped were **the same quantity computed two ways**:

- `det(E)` returned `exp(det(A))` while `eigen_det(E)` returned `det(A)`
- `llt$rcond()` and `eigen_rcond(E)` returned different numbers
- `eigen_llt()` accepted an indefinite matrix and quietly factorised a
  *different* matrix

So `consistency.R` opens every door and compares them to each other: `det()`,
`eigen_det()`, `$det()`, `exp(eigen_logdet())` and LAPACK; `eigen_logdet()`,
`$logdet()`, `$log_abs_det()` and `determinant()$modulus`; `eigen_rank()` vs
`$rank()` from three decompositions; `rcond()` through four paths; `solve()`
through six; `inverse()` through four; plus the Moore-Penrose identities and
transpose invariance. If two functions claim to compute the same mathematical
object, they must agree.

It also pins the *documented* disagreements, so they cannot drift silently.
`eigen_rank()` uses Eigen's threshold (`eps · min(nrow, ncol)`, ≈1.3e-15) while
`qr(x)$rank` uses `tol = 1e-7`, so a matrix with `σmin/σmax = 1e-12` is full rank
to one and rank-deficient to the other. Both are correct; `eigen_rank(x, 1e-7)`
reproduces `qr(x)$rank`. Likewise `eigen_rcond()` is the 2-norm reciprocal
condition number, matching `1/kappa(x, exact = TRUE)` — *not* base `rcond()`,
which is 1-norm based.

The same section checks every one-shot helper against its base R counterpart
(`eigen_sum` vs `sum`, `eigen_reshape` vs `matrix`, `eigen_rep` vs `cbind`/
`rbind`, `eigen_diag` vs `diag`, `eigen_which_min` vs `which(..., arr.ind)`,
…) and every generic against the one-shot function that bypasses it
(`sum(E)` vs `eigen_sum(E)`, `t(E)` vs `eigen_adjoint(E)`, …).

`sparse.R` is written against dense equivalents rather than against sparse
identities, because indexing is the dominant risk in a sparse layer: R counts
from 1, Eigen from 0, storage is column-major, and a transpose that forgets to
swap `nrow`/`ncol` still produces a plausibly shaped answer. So every sparse
operation is compared with the dense matrix it should equal. The solvers get
their own section, because a sparse factorisation that fails reports it through
`info()` and keeps going — `SimplicialLLT` on an indefinite matrix returns the
solution of a *different* system rather than an error. The suite asserts that
each method refuses instead, and separately that `SparseQR` is the one accepting
rectangular input.

That file found three defects of its own. `Ops.eigenSparse` had an inverted
unary/binary test, so `As + As` returned `As` unchanged and `As - As` became
`-As`. The sparse solvers checked `b` for `NaN`/`Inf` but not `A`, because
sparse input takes a conversion path that bypasses the dense finite check. And
`solve()` on a sparse matrix with a plain vector right-hand side returned a
column `eigenMat` while `solve()` on a dense `eigenMat` returned a vector.

One convention is worth stating plainly: **`eigen_dot()` does not conjugate.**
It is `sum(a * b)`, following R. Eigen's `MatrixBase::dot()` *does* conjugate —
its own documentation says it "returns the hermitian (sesquilinear) dot
product, conjugate-linear in the first variable". For real vectors the two
coincide; for complex ones they do not, and `eigen_dot(x, x)` is therefore not
the squared norm (`sum(x * x)` can be complex). Use `eigen_sumsq(x)` for
`sum(abs(x)^2)`, or `eigen_adjoint(a) %*% b` for the Hermitian inner product.

Everything is compared against base R wherever base R has an equivalent (`%*%`,
`solve`, `determinant`, `svd`, `eigen`, `kappa`, `qr.solve`), which routes to
LAPACK — an implementation independent of Eigen's.

The suite exists because eight bugs only showed up under it: a segfault that
appeared in about one run in four (see the `eigen_assert` note above), `rcond()`
returning 1 for a singular matrix, decompositions silently returning `NaN` for
input containing `NA`, `det()` returning `exp(det(A))`, `E[1]` returning the
first *row* instead of the first element, `!E` being inverted, `eigen_llt()`
returning a factorisation of a different matrix for indefinite input, and
`eigen_solve()` refusing an `eigenMat` because it carries the `eigenObject` class
too. The third
was found only by the randomised file — the hand-written cases used clean input
throughout. The last five were found by `tests/compat.R`, by a benchmark, and by running
the whole suite under GC torture, which exists because a hand-written test can encode the
same misunderstanding as the implementation (`modulus == det(m)` passed for
exactly that reason) while a systematic comparison against base R cannot.

Two of the assertions in `tests/fuzz.R` are deliberately *not* tight, and the
comments say why. `det()` is not compared at all once the condition number passes
`1e8`, because the relative error of a determinant is amplified by `kappa` — a
mismatch there would measure the BLAS, not this package. `pinv()` is checked for
the Penrose conditions only while `kappa < 1e8`, because a pseudo-inverse that
truncates small singular values is not a pseudo-inverse of the original matrix.

## Not covered yet

- **Sparse**: complex-valued sparse matrices. `Eigen::SparseMatrix<complex>`
  itself supports them; the conversion layer, the solvers and the tests would all
  need a second pass, so passing a complex matrix in stops rather than silently
  dropping the imaginary part.
- **Geometry**: `AlignedBox`, `Hyperplane`, `ParametrizedLine` and
  `Umeyama`'s variant that estimates the scale only have no R objects yet —
  `Quaternion`, `AngleAxis`, `Transform`, `Translation`, `Scaling`,
  `Rotation2D`, `umeyama` and `cross` do.
- **Special functions**: `Eigen::numext` also exposes `expint`, `log1p`,
  `expm1` and the `lgamma`-family variants; they are reachable from C++ through
  `#include <ROpenEigen.h>` but not exported to R, because R already has
  `log1p`, `expm1` and `lgamma` and the rest would be name collisions for no
  gain.
- **Eigenvalues of a real non-symmetric matrix** come back complex via
  `eigen_geigen()`. Eigen's `EigenSolver` would keep them real where possible; the
  interface cost was not worth it.

## Using Eigen from another package

`LinkingTo: ROpenEigen` and then `#include <ROpenEigen.h>` gives you Rcpp plus Eigen 5.0.1
(Dense, Sparse, Geometry, and the unsupported modules). Nothing needs to be linked — Eigen is header-only.

## Licence

GPL-2. Eigen itself is MPL-2 and is redistributed unmodified under `inst/include`
(see `inst/COPYRIGHTS`).

The package began as a rename of the RcppEigen glue headers — the pre-rewrite C++ and R
layers are kept in `src/legacy/` and `R/legacy/` for reference, and those files retain
their original Bates / Eddelbuettel / François copyright notices. Nothing from that layer
is compiled or installed any more.
