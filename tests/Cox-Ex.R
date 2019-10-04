
library("mlt")
library("survival")
set.seed(29)

### true dgp
rY <- function(n, ...) rexp(n, ...)
pY <- function(x, ...) pexp(x, ...)
dY <- function(x, ...) dexp(x, ...)

### tree groups
gf <- gl(3, 1)
g <- rep(gf, 100)
y <- rY(length(g), rate = (1:nlevels(g))[g])
mydata <- data.frame(y = y, g = g)

boxplot(y ~ g, data = mydata)

### uncensored, Cox model, h = bernstein
Bb <- Bernstein_basis(numeric_var("y", support = c(0, max(y) + .1), bounds = c(0, Inf)), 
                      order = 5, ui = "increasing")
s <- as.basis(~ g, data = data.frame(g = gf), remove_intercept = TRUE)
m <- ctm(response = Bb, shifting = s, todist = "MinExtrVal")
(cf1 <- coef(opt <- mlt(m, data = mydata)))
coef(cph <- coxph(Surv(y, rep(TRUE, nrow(mydata))) ~ g, data = mydata))
yn <- mkgrid(Bb, 50)$y
yn <- yn[yn > 0]
a <- predict(opt, newdata = data.frame(g = gf[1]), q = yn)
layout(matrix(1:4, ncol = 2))
plot(yn, a, type = "l", col = "red")
lines(yn, log(yn))
a <- predict(opt, newdata = data.frame(g = gf), q = yn, type = "survivor")
plot(yn, a[,1], type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[1])))
plot(yn, a[,2], type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[2])))
plot(yn, a[,3], type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[3])))

### h = bernstein(log())
logBb <- Bernstein_basis(numeric_var("y", support = c(1, max(y) + .1), bounds = c(min(y) / 2, Inf)), 
                         order = 5, ui = "increasing", log_first = TRUE)
m <- ctm(response = logBb, shifting = s, todist = "MinExtrVal")
(cf1 <- coef(opt <- mlt(m, data = mydata)))
## sample from this model
sam <- simulate(opt, newdata = data.frame(g = gf), nsim = 100)
nd <- data.frame(y = unlist(sam), g = rep(gf, length(sam)))
opt2 <- mlt(m, data = nd)
## visualise
yn <- mkgrid(Bb, 50)$y
yn <- yn[yn > 0]
a <- predict(opt, newdata = data.frame(g = gf[1]), q = yn)
layout(matrix(1:4, ncol = 2))
plot(yn, a, type = "l", col = "red")
lines(yn, log(yn))
a <- predict(opt, newdata = data.frame(g = gf), q = yn, type = "survivor")
plot(yn, a[,1], type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[1])))
plot(yn, a[,2], type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[2])))
plot(yn, a[,3], type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[3])))

### right censoring
mydata <- data.frame(y = Surv(y, sample(0:1, length(y), replace = TRUE)), g = g)
coef(opt <- mlt(m, data = mydata))
coef(cph <- coxph(y ~ g, data = mydata))

### left censoring
mydata <- data.frame(y = Surv(y, sample(0:1, length(y), replace = TRUE), type = "left"), g = g)
coef(opt <- mlt(m, data = mydata))

### interval censoring
mydata <- data.frame(y = Surv(y, y + 1, sample(0:3, length(y), replace = TRUE), type = "interval"), 
                     g = g)
coef(opt <- mlt(m, data = mydata))

### uncensored, time-varying coefficients in both groups
mydata <- data.frame(y = y, g = g)
m <- ctm(response = logBb, 
           interacting = as.basis(~ g, data = mydata),
           todist = "MinExtrVal")
coef(opt <- mlt(m, data = mydata, maxit = 5000))
coef(cph <- coxph(Surv(y, rep(TRUE, nrow(mydata))) ~ g, data = mydata))
## visualize
a <- predict(opt, newdata = data.frame(g = gf[1]), q = yn)
layout(matrix(1:4, ncol = 2))
plot(yn, a, type = "l", col = "red")
lines(yn, log(yn))
a <- predict(opt, newdata = data.frame(g = gf), q = yn, type = "survivor")
plot(yn, a[,1], type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[1])))
plot(yn, a[,2], type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[2])))
plot(yn, a[,3], type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[3])))

