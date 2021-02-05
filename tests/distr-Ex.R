
library("mlt")

### 10.1080/15598608.2013.772835
### rho = exp(logrho)
### 1 = rho = exp(0) is identical to .Logistic
lg <- mlt:::.GammaFrailty(logrho = 0)
x <- -30:30 / 20
p <- 1:99 / 100
tol <- sqrt(.Machine$double.eps)
stopifnot(max(abs(lg$p(x) - mlt:::.Logistic()$p(x))) < tol)

stopifnot(max(abs(lg$d(x) - mlt:::.Logistic()$d(x))) < tol)

stopifnot(max(abs(lg$dd(x) - mlt:::.Logistic()$dd(x))) < tol)

stopifnot(max(abs(lg$ddd(x) - mlt:::.Logistic()$ddd(x))) < tol)

stopifnot(max(abs(lg$dd2d(x) - mlt:::.Logistic()$dd2d(x))) < tol)

stopifnot(max(abs(lg$q(p) - mlt:::.Logistic()$q(p))) < tol)

### 0 = rho = exp(-infty) is identical to .MinExtrVal
tol <- (.Machine$double.eps)^(1 / 2.5)
lg <- mlt:::.GammaFrailty(logrho = -18)
stopifnot(max(abs(lg$p(x) - mlt:::.MinExtrVal()$p(x))) < tol)

stopifnot(max(abs(lg$d(x) - mlt:::.MinExtrVal()$d(x))) < tol)

stopifnot(max(abs(lg$dd(x) - mlt:::.MinExtrVal()$dd(x))) < tol)

stopifnot(max(abs(lg$ddd(x) - mlt:::.MinExtrVal()$ddd(x))) < tol)

stopifnot(max(abs(lg$dd2d(x) - mlt:::.MinExtrVal()$dd2d(x))) < tol)

stopifnot(max(abs(lg$q(p) - mlt:::.MinExtrVal()$q(p))) < tol)

### logitalpha = logit(1) is identical to .MinExtrVal
tol <- (.Machine$double.eps)^(1 / 2.5)
lg <- mlt:::.PositiveStableFrailty(logitalpha = 50)
stopifnot(max(abs(lg$p(x) - mlt:::.MinExtrVal()$p(x))) < tol)

stopifnot(max(abs(lg$d(x) - mlt:::.MinExtrVal()$d(x))) < tol)

stopifnot(max(abs(lg$dd(x) - mlt:::.MinExtrVal()$dd(x))) < tol)

stopifnot(max(abs(lg$ddd(x) - mlt:::.MinExtrVal()$ddd(x))) < tol)

stopifnot(max(abs(lg$dd2d(x) - mlt:::.MinExtrVal()$dd2d(x))) < tol)

stopifnot(max(abs(lg$q(p) - mlt:::.MinExtrVal()$q(p))) < tol)

