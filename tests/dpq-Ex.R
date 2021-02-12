
library("mlt")
library("numDeriv")
set.seed(29)
options(digits = 4)

n <- 20
### just for interface checking
### we need something better!
d <- data.frame(x1 = 1:n, x2 = sample(1:n) + 1, y = rnorm(n))
m <- ctm(polynomial_basis(numeric_var("y", support = range(d$y)),
                            coef = c(TRUE, TRUE), ci = c(-Inf, 0)),
           shift = ~ x1 + x2, data = d)
mod <- mlt(m, data = d)
coef(mod)

p <- predict(mod, newdata = d)

K <- 15
q <- mkgrid(m, n = K)[["y"]]
p1 <- predict(mod, newdata = d[, c("x1", "x2")], q = q)
p2 <- predict(mod, newdata = d[, c("x1", "x2")], K = K)
stopifnot(all.equal(p1, p2))

p0 <- predict(mod$model$model, 
    newdata = expand.grid(d), coef = coef(mod))
p1 <- predict(mod, newdata = as.list(d))
p2 <- predict(mod, newdata = d, q = d$y[1])

max(abs(p0 - as.vector(p1)))

all.equal(p1[cbind(1:n, 1:n, 1), drop = TRUE],
          drop(p2))

all.equal(p1[cbind(1:n, 1:n, 1:n), drop = TRUE],
          drop(p), check.attributes = FALSE)

predict(mod, newdata = list(x1 = 1:3, x2 = 2:3), prob = c(.25, .5), type = "quantile")

simulate(mod, nsim = 1, seed = 291, interpolate = FALSE)

d$y <- gl(3, 1, ordered = TRUE)[rep(1:3, length = n)]

r <- as.basis(d$y) #as.basis(~ y, data = d, remove_intercept = TRUE,
#              contrasts.arg = list(y = function(n)
#                  contr.treatment(n, base = 3)),
#              ui = diff(diag(2)), ci = 0)

mod2 <- mlt(ctm(r, shift = ~ x1 + x2, data = d), data = d)

predict(mod2, q = unique(d$y))

predict(mod2, prob = 1:9 / 10, type = "quantile")

simulate(mod2, nsim = 3, seed = 29)

predict(mod2, q = unique(d$y), type = "density")

predict(mod2, list(y = unique(d$y), x1 = 1:3, x2 = 2:3), type = "density")

### some basis checks: continuous
d <- data.frame(x1 = 1:n, x2 = sample(1:n) + 1, y = rnorm(n))
m <- ctm(polynomial_basis(numeric_var("y", support = range(d$y)),
                            coef = c(TRUE, TRUE), ci = c(-Inf, 0)),
           shift = ~ x1 + x2, data = d)
mod <- mlt(m, data = d)

.chk <- function(x)
    stopifnot(isTRUE(max(abs(x), na.rm = TRUE) < sqrt(.Machine$double.eps)))

cont <- quote({
nd <- d
nd$y <- NULL
q <- mkgrid(mod, 10)[[1]]
p <- predict(mod, newdata = nd, q = q, type = "distribution")
s <- predict(mod, newdata = nd, q = q, type = "survivor")
.chk(predict(mod, newdata = nd, q = q, type = "distribution", log = TRUE) - log(p))
.chk(predict(mod, newdata = nd, q = q, type = "distribution", lower.tail = FALSE) - s)
.chk(predict(mod, newdata = nd, q = q, type = "distribution", lower.tail =
FALSE, log = TRUE) - log(s))

o <- predict(mod, newdata = nd, q = q, type = "odds")
.chk(o - p / s)

df <- function(q)
    predict(mod, newdata = nd[1,], q = q, type = "distribution")
da <- sapply(q, function(q) grad(df, q))
dd <- predict(mod, newdata = nd[1,], q = q, type = "density")
.chk(da - dd)

h <- predict(mod, newdata = nd, q = q, type = "hazard")
.chk(predict(mod, newdata = nd, q = q, type = "loghazard") - log(h))

H <- predict(mod, newdata = nd, q = q, type = "cumhazard")
.chk(H + log(s))

.chk(predict(mod, newdata = nd, q = q, type = "logcumhazard") - log(H))

dh <- function(q)
    predict(mod, newdata = nd[1,], q = q, type = "cumhazard")

da <- sapply(q, function(q) grad(dh, q))
dd <- predict(mod, newdata = nd[1,], q = q, type = "hazard")
.chk(da - dd)
})

m <- ctm(polynomial_basis(numeric_var("y", support = range(d$y)),
                            coef = c(TRUE, TRUE), ci = c(-Inf, 0)),
           shift = ~ x1 + x2, data = d, todistr = "Normal")
mod <- mlt(m, data = d)

eval(cont)

m <- ctm(polynomial_basis(numeric_var("y", support = range(d$y)),
                            coef = c(TRUE, TRUE), ci = c(-Inf, 0)),
           shift = ~ x1 + x2, data = d, todistr = "Logistic")
mod <- mlt(m, data = d)

eval(cont)

m <- ctm(polynomial_basis(numeric_var("y", support = range(d$y)),
                            coef = c(TRUE, TRUE), ci = c(-Inf, 0)),
           shift = ~ x1 + x2, data = d, todistr = "MinExtrVal")
mod <- mlt(m, data = d)

eval(cont)

m <- ctm(polynomial_basis(numeric_var("y", support = range(d$y)),
                            coef = c(TRUE, TRUE), ci = c(-Inf, 0)),
           shift = ~ x1 + x2, data = d, todistr = "MaxExtrVal")
mod <- mlt(m, data = d)

eval(cont)

### some basis checks: discrete

disc <- quote({
nd <- d
nd$y <- NULL
q <- mkgrid(mod, 10)[[1]]
p <- predict(mod, newdata = nd, q = q, type = "distribution")
s <- predict(mod, newdata = nd, q = q, type = "survivor")
.chk(predict(mod, newdata = nd, q = q, type = "distribution", log = TRUE) - log(p))
.chk(predict(mod, newdata = nd, q = q, type = "distribution", lower.tail = FALSE) - s)
.chk(predict(mod, newdata = nd, q = q, type = "distribution", lower.tail =
FALSE, log = TRUE) - log(s))

o <- predict(mod, newdata = nd, q = q, type = "odds")
.chk(o - p / s)

dd <- predict(mod, newdata = nd, q = q, type = "density")

.chk(apply(p, 2, function(x) diff(c(0, x))) - dd)

h <- predict(mod, newdata = nd, q = q, type = "hazard")

.chk(dd / (1 - (p - dd)) - h)

.chk(apply(h, 2, function(x) cumprod(1 - x)) - s)

H <- predict(mod, newdata = nd, q = q, type = "cumhazard")

.chk(H + log(s))
})

d <- data.frame(x1 = 1:n, x2 = sample(1:n) + 1, y = sample(gl(4, 5, ordered = TRUE)))
m <- ctm(as.basis(d$y), shift = ~ x1 + x2, data = d, todistr = "Normal")
mod <- mlt(m, data = d)

eval(disc)

d <- data.frame(x1 = 1:n, x2 = sample(1:n) + 1, y = sample(gl(4, 5, ordered = TRUE)))
m <- ctm(as.basis(d$y), shift = ~ x1 + x2, data = d, todistr = "Logistic")
mod <- mlt(m, data = d)

eval(disc)

d <- data.frame(x1 = 1:n, x2 = sample(1:n) + 1, y = sample(gl(4, 5, ordered = TRUE)))
m <- ctm(as.basis(d$y), shift = ~ x1 + x2, data = d, todistr = "MinExtrVal")
mod <- mlt(m, data = d)

eval(disc)

d <- data.frame(x1 = 1:n, x2 = sample(1:n) + 1, y = sample(gl(4, 5, ordered = TRUE)))
m <- ctm(as.basis(d$y), shift = ~ x1 + x2, data = d, todistr = "MaxExtrVal")
mod <- mlt(m, data = d)

eval(disc)

