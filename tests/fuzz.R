# Randomised property tests.
#
# core.R and edge.R cover the cases that came to mind.  This file covers the ones
# that did not: 250 matrices with random shapes (1..8), six flavours of input and
# properties that must hold whatever the input looks like.  Everything is checked
# against base R, which routes to LAPACK, so the two implementations are genuinely
# independent.
#
# Seeded, so a failure is reproducible and the run is deterministic.
#
# Two traps worth remembering when writing assertions here:
#
#   1. A matrix containing NA makes `a - b` NA, so a naive max() of the error
#      reports NA and poisons the assertion.  Compare NA positions explicitly and
#      drop them before measuring the error.
#   2. det() and rcond() of a rank-deficient or ill-conditioned matrix live around
#      1e-17.  Relative comparison there reports a failure for two answers that
#      are both numerically zero.  Use atol + rtol * |expected|, and do not
#      compare det() at all once the condition number is huge -- the relative
#      error of a determinant is amplified by kappa, so the number is not
#      meaningful at 1e-6 and comparing it tests the BLAS, not this package.

suppressPackageStartupMessages(library(ROpenEigen))

set.seed(20261103)
NCASE <- 250
fails <- list()
nchk <- 0L

rec <- function(case, what, detail = "")
  fails[[length(fails) + 1L]] <<- list(case = case, what = what, detail = detail)

ck <- function(cond, case, what, detail = "") {
  nchk <<- nchk + 1L
  if (!isTRUE(cond)) rec(case, what, detail)
}

close_ab <- function(a, b, rtol = 1e-10, atol = 1e-12) {
  a <- as.vector(a); b <- as.vector(b)
  if (length(a) != length(b)) return(FALSE)
  na_a <- is.na(a); na_b <- is.na(b)
  if (!identical(na_a, na_b)) return(FALSE)
  keep <- !na_a
  if (!any(keep)) return(TRUE)
  d <- Mod(a[keep] - b[keep])
  s <- Mod(b[keep])
  all(d <= atol + rtol * s)
}

cat("== randomised properties ==\n")

for (k in seq_len(NCASE)) {
  nr <- sample(1:8, 1); nc <- sample(1:8, 1)
  kind <- sample(c("real", "cplx", "int", "logi", "illcond", "rankdef"), 1)

  m <- switch(kind,
    real    = matrix(rnorm(nr * nc), nr, nc),
    cplx    = matrix(complex(real = rnorm(nr * nc), imaginary = rnorm(nr * nc)), nr, nc),
    int     = matrix(sample(-5:5, nr * nc, TRUE), nr, nc),
    logi    = matrix(sample(c(TRUE, FALSE, NA), nr * nc, TRUE), nr, nc),
    illcond = { s <- 10^seq(0, -12, length.out = min(nr, nc))
                qr.Q(qr(matrix(rnorm(nr * min(nr, nc)), nr)))[, seq_len(min(nr, nc)), drop = FALSE] %*%
                  diag(s, length(s)) %*%
                  t(qr.Q(qr(matrix(rnorm(nc * min(nr, nc)), nc)))) },
    rankdef = { r <- max(1, min(nr, nc) - 1)
                matrix(rnorm(nr * r), nr, r) %*% matrix(rnorm(r * nc), r, nc) }
  )
  if (kind %in% c("illcond", "rankdef")) dim(m) <- c(nr, nc)

  E <- eigenMat(m)
  tag <- sprintf("#%d %s %dx%d", k, kind, nr, nc)
  hasNA <- kind == "logi" && anyNA(m)

  # 1. round trip
  ck(identical(dim(E), c(nr, nc)), tag, "dim")
  ck(close_ab(as.matrix(E), m, 0, 1e-15), tag, "as.matrix round trip")
  ck(length(E) == nr * nc, tag, "length")

  # 2. transpose and conjugate transpose
  ck(close_ab(as.matrix(t(E)), t(m), 0, 1e-15), tag, "t()")
  ck(close_ab(as.matrix(eigen_adjoint(E)), Conj(t(m)), 0, 1e-15), tag, "adjoint")

  # 3. elementwise
  ck(close_ab(as.matrix(E + 2), m + 2, 0, 1e-14), tag, "E + 2")
  ck(close_ab(as.matrix(2 * E), 2 * m, 0, 1e-14), tag, "2 * E")
  ck(close_ab(as.matrix(Conj(E)), Conj(m), 0, 1e-15), tag, "Conj")

  # 4. reductions
  if (!hasNA) {
    ck(close_ab(sum(E), sum(m), 1e-12, 1e-12), tag, "sum")
    ck(close_ab(mean(E), mean(m), 1e-12, 1e-12), tag, "mean")
    # base R's min()/max() reject complex input too, so this is not a gap here
    if (kind != "cplx") {
      ck(close_ab(min(E), min(m), 0, 1e-15), tag, "min")
      ck(close_ab(max(E), max(m), 0, 1e-15), tag, "max")
    }
    ck(close_ab(eigen_trace(E), sum(diag(m)), 1e-12, 1e-12), tag, "trace")
  } else {
    ck(close_ab(sum(E, na.rm = TRUE), sum(m, na.rm = TRUE), 1e-12, 1e-12),
       tag, "sum na.rm")
    ck(close_ab(mean(E, na.rm = TRUE), mean(m, na.rm = TRUE), 1e-12, 1e-12),
       tag, "mean na.rm")
  }

  # 5. rcond is sigma_min / sigma_max, so it equals 1/kappa for any shape.
  #    The two sides come from different SVD implementations (Eigen's Jacobi vs
  #    LAPACK's divide and conquer), so the tolerance cannot be machine epsilon
  #    once kappa is in the 1e12 range.
  kk <- try(kappa(m, exact = TRUE), silent = TRUE)
  if (!inherits(kk, "try-error") && is.finite(kk) && kk > 0) {
    ck(close_ab(eigen_rcond(E), 1 / kk, 1e-4, 1e-14), tag, "rcond = 1/kappa",
       sprintf("got %.6g want %.6g", eigen_rcond(E), 1 / kk))
  }

  # 6. matrix product
  p <- sample(1:6, 1)
  mm <- matrix(rnorm(nc * p), nc, p)
  if (kind == "cplx") mm <- mm + 1i * matrix(rnorm(nc * p), nc, p)
  ck(close_ab(as.matrix(E %*% eigenMat(mm)), m %*% mm, 1e-10, 1e-12), tag, "%*%")

  # 7. square-only operations
  if (hasNA) {
    # Decompositions reject non-finite input, as base R's svd()/qr() do.
    ck(inherits(try(eigen_svd(E), silent = TRUE), "try-error"),
       tag, "decomposition rejects NA")
  } else if (nr == nc) {
    # base R's det() has no complex method, so use the product of the eigenvalues
    d_ref <- if (kind == "cplx") prod(eigen(m, only.values = TRUE)$values) else det(m)
    d <- try(eigen_det(E), silent = TRUE)
    if (!inherits(d, "try-error")) {
      kd <- try(kappa(m, exact = TRUE), silent = TRUE)
      if (inherits(kd, "try-error") || !is.finite(kd) || kd < 1e8) {
        ck(close_ab(d, d_ref, 1e-6, 1e-10), tag, "det",
           sprintf("got %s want %s", format(d), format(d_ref)))
        # det() is not a generic, but it calls determinant(), which dispatches --
        # so it silently becomes sign * exp(modulus).  If modulus were |det|
        # rather than log|det|, this would come back as exp(det) and no other
        # assertion here would notice.
        dm <- determinant(E)
        ck(identical(class(dm), "det"), tag, "determinant class")
        ck(close_ab(det(E), d_ref, 1e-6, 1e-10), tag, "det() via determinant")
        ck(close_ab(dm$sign * exp(dm$modulus), d, 1e-8, 1e-10),
           tag, "det = sign * exp(modulus)")
      } else {
        ck(is.finite(d), tag, "det finite on ill-conditioned input")
        ck(is.finite(determinant(E)$modulus), tag, "determinant modulus finite")
      }
    }

    sv <- try(eigen_svd(E), silent = TRUE)
    if (!inherits(sv, "try-error")) {
      ck(close_ab(as.matrix(sv$reconstructed()), m, 1e-9, 1e-11), tag, "svd rebuild")
      # The pseudo-inverse truncates small singular values, so on an
      # ill-conditioned matrix it does not satisfy the Penrose conditions.  That
      # is the documented behaviour of a truncated inverse, not a defect.
      kk2 <- kappa(m, exact = TRUE)
      if (is.finite(kk2) && kk2 < 1e8) {
        P <- as.matrix(eigen_pinv(E))
        ck(close_ab(m %*% P %*% m, m, 1e-7, 1e-9), tag, "pinv A P A = A")
        ck(close_ab(P %*% m %*% P, P, 1e-7, 1e-9), tag, "pinv P A P = P")
      }
    }
  }

  # 8. least squares.  An overdetermined system has a non-zero residual, so the
  #    comparison has to be between fitted values m %*% sol, not between the
  #    solutions themselves.
  if (!hasNA) {
    b <- matrix(rnorm(nr), nr, 1)
    if (kind == "cplx") b <- b + 1i * matrix(rnorm(nr), nr, 1)
    sol <- try(as.matrix(eigen_svd(E)$solve(b)), silent = TRUE)
    if (!inherits(sol, "try-error")) {
      ck(all(is.finite(sol)), tag, "svd$solve finite")
      if (kind %in% c("real", "cplx", "int")) {
        ref <- try(qr.solve(m, b), silent = TRUE)
        if (!inherits(ref, "try-error"))
          ck(close_ab(m %*% sol, m %*% ref, 1e-6, 1e-8),
             tag, "svd$solve fitted values = qr.solve")
      }
    }
  }

  # 9. subscripting
  if (nr >= 2 && nc >= 2) {
    ck(close_ab(as.matrix(eigen_block(E, 2:nr, 2:nc)), m[2:nr, 2:nc, drop = FALSE], 0, 1e-15),
       tag, "block")
    ck(close_ab(as.matrix(E[-1, ]), m[-1, , drop = FALSE], 0, 1e-15), tag, "negative index")
  }

  # 10. binding
  ck(close_ab(as.matrix(eigen_rbind(E, E)), rbind(m, m), 0, 1e-15), tag, "rbind")
  ck(close_ab(as.matrix(eigen_cbind(E, E)), cbind(m, m), 0, 1e-15), tag, "cbind")
}

cat(sprintf("\nchecked %d, failed %d\n", nchk, length(fails)))
if (length(fails)) {
  cat("first failures:\n")
  for (i in seq_len(min(25L, length(fails)))) {
    f <- fails[[i]]
    cat(sprintf("  %s | %s | %s\n", f$case, f$what, f$detail))
  }
  stop(sprintf("%d randomised assertion(s) failed", length(fails)), call. = FALSE)
}
