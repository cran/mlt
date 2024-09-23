
library("mlt")
library("mvtnorm")
set.seed(29)

chk <- function(...) 
    stopifnot(isTRUE(all.equal(..., tol = 1e-6, check.attributes = FALSE)))

thischeck <- expression({
  nm <- LETTERS[1:J]
  ltM <- function(x) ltMatrices(x, diag = FALSE, byrow = TRUE, names = nm)
  ltD <- function(x) ltMatrices(x, diag = TRUE, byrow = TRUE, names = nm)
  prm <- matrix(runif(J * (J - 1) / 2 * N), ncol = N)
  L <- ltM(prm)

  obs <- matrix(rnorm(J * N), ncol = N)
  rownames(obs) <- nm
  lwr <- -2 -abs(obs)
  upr <- 2 + abs(obs)

  fun <- mlt:::.ll(c(J, 0), standardize = FALSE)
  sum(fun$logLik(obs, L))
  s <- fun$score(obs, L)

  f <- function(obs = obs, L = L) 
    sum(fun$logLik(obs, ltM(L)))
  S <- matrix(grad(f, unclass(L), obs = obs), ncol = N)
  chk(S, Lower_tri(s$Lambda))
  S <- matrix(grad(f, obs, L = L), ncol = N)
  chk(S, s$obs)

  fun <- mlt:::.ll(c(J, 0), standardize = TRUE)
  sum(fun$logLik(obs, L))

  LD <- invcholD(L)
  sum(-colSums(Mult(LD, obs)^2 - log(diagonals(LD)^2)))

  C <- solve(L)
  CCt <- Tcrossprod(C, diag_only = TRUE)
  LD <- invcholD(L, D = sqrt(CCt))
  sum(-colSums(Mult(LD, obs)^2 - log(CCt)))

  f <- function(obs, a) {
    LD <- ltMatrices(a, diag = TRUE, byrow = TRUE)
    sum(-colSums(Mult(LD, obs)^2 - log(diagonals(LD)^2)))
  }
  f(obs = obs, LD)

  sLD <- ltMatrices(matrix(grad(f, unclass(LD), obs = obs), ncol = N), 
                    diag = TRUE, byrow = TRUE)
  sobs <- matrix(grad(f, obs, a = unclass(LD)), ncol = N)

  sLDfun <- function(obs, LD) {
    cJ <- dim(LD)[2L]
    Y <- matrix(obs, byrow = TRUE, nrow = cJ, ncol = N * cJ)
    tmp <- -2 * Mult(LD, Mult(LD, obs), transpose = TRUE)
    ret <- - 2 * matrix(Mult(LD, obs)[, rep(1:N, each = cJ)] * Y, ncol = N)
    M <- matrix(1:(cJ^2), nrow = cJ, byrow = FALSE)
    ret <- ltMatrices(ret[M[lower.tri(M, diag = TRUE)],,drop = FALSE], 
                          diag = TRUE, byrow = FALSE)
    ret <- ltMatrices(ret, 
                      diag = TRUE, byrow = TRUE)
    diagonals(ret) <- diagonals(ret) + 2 / diagonals(LD)
    return(list(Lambda = ret, obs = tmp))
  }

  s <- sLDfun(obs, LD)
  chk(s$Lambda, sLD)
  chk(s$obs, sobs)

  s <- fun$score(obs, L)

  f <- function(obs = obs, L = L) 
    sum(fun$logLik(obs, ltM(L)))
  S <- matrix(grad(f, unclass(L), obs = obs), ncol = N)
  chk(S, Lower_tri(s$Lambda))
  S <- matrix(grad(f, obs, L = L), ncol = N)
  chk(S, s$obs)

  w <- matrix(runif((J - 1) * M), ncol = M)
  fun <- mlt:::.ll(c(0, J), standardize = FALSE, list(w = w))
  sum(fun$logLik(lower = lwr, upper = upr, Lambda = L))
  s <- fun$score(lower = lwr, upper = upr, Lambda = L)

  f <- function(lwr = lwr, upr = upr, L = L) 
    sum(fun$logLik(lower = lwr, upper = upr, Lambda = ltM(L)))
  S <- grad(f, unclass(L), lwr = lwr, upr = upr)
  chk(S, c(Lower_tri(s$Lambda)))

  S <- matrix(grad(f, lwr, upr = upr, L = L), ncol = N)
  chk(S, s$lower)
  S <- matrix(grad(f, upr, lwr = lwr, L = L), ncol = N)
  chk(S, s$upper)

  fun <- mlt:::.ll(c(0, J), standardize = TRUE, list(w = w))
  sum(fun$logLik(lower = lwr, upper = upr, Lambda = L))
  s <- fun$score(lower = lwr, upper = upr, Lambda = L)

  f <- function(lwr = lwr, upr = upr, L = L) 
    sum(fun$logLik(lower = lwr, upper = upr, Lambda = ltM(L)))
  S <- matrix(grad(f, unclass(L), lwr = lwr, upr = upr), ncol = N)
  chk(S, Lower_tri(s$Lambda))
  S <- matrix(grad(f, lwr, upr = upr, L = L), ncol = N)
  chk(S, s$lower)
  S <- matrix(grad(f, upr, lwr = lwr, L = L), ncol = N)
  chk(S, s$upper)

  w <- matrix(runif((dJ - 1) * M), ncol = M)
  fun <- mlt:::.ll(c(cJ, dJ), standardize = FALSE, list(w = w))
  sum(fun$logLik(obs[1:cJ,,drop = FALSE], lwr[-(1:cJ),,drop = FALSE], upr[-(1:cJ),,drop = FALSE], L))
  s <- fun$score(obs[1:cJ,,drop = FALSE], lwr[-(1:cJ),,drop = FALSE], upr[-(1:cJ),,drop = FALSE], L)
  f <- function(obs = obs[1:cJ,,drop = FALSE], lwr = lwr[-(1:cJ),,drop = FALSE], upr = upr[-(1:cJ),,drop = FALSE], L = L) 
    sum(fun$logLik(obs, lwr, upr, ltM(L)))
  S <- matrix(grad(f, unclass(L), obs = obs[1:cJ,,drop = FALSE], lwr = lwr[-(1:cJ),,drop = FALSE], 
                 upr = upr[-(1:cJ),,drop = FALSE]), ncol = N)
  chk(S, Lower_tri(s$Lambda))
  S <- matrix(grad(f, obs[1:cJ,,drop = FALSE], lwr = lwr[-(1:cJ),,drop = FALSE], 
                    upr = upr[-(1:cJ),,drop = FALSE], L = L), ncol = N)
  chk(S, s$obs)
  S <- matrix(grad(f, lwr[-(1:cJ),,drop = FALSE], obs = obs[1:cJ,,drop = FALSE], 
                    upr = upr[-(1:cJ),,drop = FALSE], L = L), ncol = N)
  chk(S, s$lower)
  S <- matrix(grad(f, upr[-(1:cJ),,drop = FALSE], obs = obs[1:cJ,,drop = FALSE], 
                    lwr = lwr[-(1:cJ),,drop = FALSE], L = L), ncol = N)
  chk(S, s$upper)

  fun <- mlt:::.ll(c(cJ, dJ), standardize = TRUE, list(w = w))
  sum(fun$logLik(obs[1:cJ,,drop = FALSE], lwr[-(1:cJ),,drop = FALSE], upr[-(1:cJ),,drop = FALSE], L))
  s <- fun$score(obs[1:cJ,,drop = FALSE], lwr[-(1:cJ),,drop = FALSE], 
                 upr[-(1:cJ),,drop = FALSE], L)

  f <- function(obs = obs[1:cJ,,drop = FALSE], lwr = lwr[-(1:cJ),,drop = FALSE], upr = upr[-(1:cJ),,drop = FALSE], L = L) 
    sum(fun$logLik(obs, lwr, upr, ltM(L)))
  S <- matrix(grad(f, unclass(L), obs = obs[1:cJ,,drop = FALSE], lwr = lwr[-(1:cJ),,drop = FALSE], 
                 upr = upr[-(1:cJ),,drop = FALSE]), ncol = N)
  chk(S, Lower_tri(s$Lambda))

  S <- matrix(grad(f, obs[1:cJ,,drop = FALSE], lwr = lwr[-(1:cJ),,drop = FALSE], 
                    upr = upr[-(1:cJ),,drop = FALSE], L = L), ncol = N)
  chk(S, s$obs)
  S <- matrix(grad(f, lwr[-(1:cJ),,drop = FALSE], obs = obs[1:cJ,,drop = FALSE], 
                    upr = upr[-(1:cJ),,drop = FALSE], L = L), ncol = N)
  chk(S, s$lower)
  S <- matrix(grad(f, upr[-(1:cJ),,drop = FALSE], obs = obs[1:cJ,,drop = FALSE], 
                    lwr = lwr[-(1:cJ),,drop = FALSE], L = L), ncol = N)
  chk(S, s$upper)

})


if (require("numDeriv", quietly = TRUE)) {

J <- (cJ <- 5) + (dJ <- 6)
N <- 3
M <- 10

eval(thischeck)

J <- (cJ <- 1) + (dJ <- 1)

eval(thischeck)

J <- (cJ <- 1) + (dJ <- 4)

eval(thischeck)

J <- (cJ <- 4) + (dJ <- 1)

eval(thischeck)

}
