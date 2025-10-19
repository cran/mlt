
library("mlt")
library("survival")

set.seed(290875)

chk <- function(x, y) stopifnot(all.equal(x, y, tol = 1e-3, check.attributes = FALSE))

### check -Inf, Inf interval censoring for handling missing
### response values
### numeric response
N <- 100
x <- rnorm(N)
ina <- sample(1:length(x))[1:floor(N / 10)]
x[ina] <- NA

d <- data.frame(x = x)
d$Rx <- d$Rxi <- R(x)
d$Rxi$cleft <- -Inf
d$Rxi$cright <- Inf
d$Rx2 <- R(x, as.R.interval = TRUE)

m <- mlt(ctm(response = Bernstein_basis(numeric_var("x"), order = 1)),
         data = d[-ina,,drop = FALSE])
mx <- mlt(ctm(response = Bernstein_basis(numeric_var("x"), order = 1)),
          data = d)
mRx <- mlt(ctm(response = Bernstein_basis(numeric_var("Rx"), order = 1)),
           data = d)
mRxi <- mlt(ctm(response = Bernstein_basis(numeric_var("Rxi"), order = 1)),
            data = d)
mRx2 <- mlt(ctm(response = Bernstein_basis(numeric_var("Rx2"), order = 1)),
            data = d)
mx2 <- mlt(ctm(response = Bernstein_basis(numeric_var("Rx2"), order = 1)),
           data = d[-ina,,drop = FALSE])

chk(coef(mx), coef(m))
chk(coef(mRx), coef(m))
chk(coef(mRxi), coef(m))

chk(logLik(mx), logLik(m))
chk(logLik(mRx), logLik(m))
chk(logLik(mRxi), logLik(m))

chk(coef(mRx2), coef(mx2))

chk(logLik(mRx2), logLik(mx2))

### ordered factors
x <- gl(4, N, ordered = TRUE)
x[ina] <- NA
d <- data.frame(x = x)

m <- mlt(ctm(response = as.basis(d$x)),
         data = d[-ina,,drop = FALSE])

mx <- mlt(ctm(response = as.basis(d$x)),
          data = d)

chk(coef(m), coef(mx))

chk(logLik(m), logLik(mx))

### integers
x <- sample(1:N)
x[ina] <- NA
d <- data.frame(x = x)
d$Rx <- R(x)

mx <- mlt(ctm(response = Bernstein_basis(numeric_var("Rx"), order = 1)),
          data = d)

mx2 <- mlt(ctm(response = Bernstein_basis(numeric_var("Rx"), order = 1)),
           data = d[-ina,,drop = FALSE])

chk(coef(mx), coef(mx2))

chk(logLik(mx), logLik(mx2))

### survival
x <- exp(rnorm(N))
x[ina] <- NA

### right censoring
d <- data.frame(y = Surv(x, event = x > 1))
tol <- sqrt(.Machine$double.eps)
m <- mlt(ctm(response = Bernstein_basis(
    numeric_var("y", bounds = c(tol, Inf), support = c(tol, 10)), 
    order = 1, log_first = TRUE)),
    data = d[-ina,,drop = FALSE])
mx <- mlt(ctm(response = Bernstein_basis(
    numeric_var("y", bounds = c(tol, Inf), support = c(tol, 10)), 
    order = 1, log_first = TRUE)),
    data = d)

chk(coef(m), coef(mx))

chk(logLik(m), logLik(mx))

### left censoring
d <- data.frame(y = Surv(x, event = x > 1, type = "left"))
tol <- sqrt(.Machine$double.eps)
m <- mlt(ctm(response = Bernstein_basis(
    numeric_var("y", bounds = c(tol, Inf), support = c(tol, 10)), 
    order = 1, log_first = TRUE)),
    data = d[-ina,,drop = FALSE])
mx <- mlt(ctm(response = Bernstein_basis(
    numeric_var("y", bounds = c(tol, Inf), support = c(tol, 10)), 
    order = 1, log_first = TRUE)),
    data = d)

chk(coef(m), coef(mx))

chk(logLik(m), logLik(mx))

### empirical likelihood
d$y <- R(Surv(x, event = x > 1), as.R.interval = TRUE)
m <- mlt(ctm(response = Bernstein_basis(
    numeric_var("y", bounds = c(tol, Inf), support = c(tol, 10)), 
    order = 1, log_first = TRUE)),
    data = d[-ina,,drop = FALSE])
mx <- mlt(ctm(response = Bernstein_basis(
    numeric_var("y", bounds = c(tol, Inf), support = c(tol, 10)), 
    order = 1, log_first = TRUE)),
    data = d)

chk(coef(m), coef(mx))

chk(logLik(m), logLik(mx))
