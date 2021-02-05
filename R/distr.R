
.Normal <- function()
    list(parm = function(x) NULL,
         p = pnorm, d = dnorm, q = qnorm, 
         ### see also MiscTools::ddnorm
         dd = function(x) -dnorm(x = x) * x,
         ddd = function(x) dnorm(x = x) * (x^2 - 1), 
         dd2d = function(x) -x,
         call = ".Normal",
         name = "normal")

.Exponential <- function()
    list(parm = function(x) NULL,
         p = pexp, d = dexp, q = qexp, 
         dd = function(x) -dexp(x = x),
         ddd = function(x) dexp(x = x),
         dd2d = function(x) -1,
         call = ".Exponential",
         name = "exponential")

.Logistic <- function()
    list(parm = function(x) NULL,
         p = plogis, d = dlogis, q = qlogis,
         dd = function(x) {
             ex <- exp(x)
             (ex - ex^2) / (1 + ex)^3
         },
         ddd = function(x) {
             ex <- exp(x)
             (ex - 4 * ex^2 + ex^3) / (1 + ex)^4
         },
         dd2d = function(x) {
             ex <- exp(x)
             (1 - ex) / (1 + ex)
         },
         call = ".Logistic",
         name = "logistic")

### Gompertz distribution
.MinExtrVal <- function()
    list(parm = function(x) NULL,
         p = function(x, lower.tail = TRUE, log.p = FALSE) {
             ### p = 1 - exp(-exp(x))
             ret <- exp(-exp(x))
             if (log.p) {
                 if (lower.tail)
                     return(log1p(-ret))
                 return(-exp(x))
             }
             if (lower.tail)
                 return(1 - exp(-exp(x)))
             return(ret)             
         },
         q = function(p) log(-log1p(- p)),
         d = function(x, log = FALSE) {
             ret <- x - exp(x)
             if (!log) return(exp(ret))
             ret
         },
         dd = function(x) {
             ex <- exp(x)
             (ex - ex^2) / exp(ex)
         },
         ddd = function(x) {
             ex <- exp(x)
             (ex - 3*ex^2 + ex^3) / exp(ex)
         },
         dd2d = function(x)
             1 - exp(x),
         call = ".MinExtrVal",
         name = "minimum extreme value")

### Gumbel distribution
.MaxExtrVal <- function()
    list(parm = function(x) NULL,
         p = function(x, lower.tail = TRUE, log.p = FALSE) {
             ### p = exp(-exp(-x))
             if (log.p) {
                 if (lower.tail)
                     return(-exp(-x))
                 return(log1p(-exp(-exp(-x))))
             }
             if (lower.tail)
                 return(exp(-exp(-x)))
             1 - exp(-exp(-x))
         },
         q = function(p) -log(-log(p)),
         d = function(x, log = FALSE) {
             ret <- - x - exp(-x)
             if (!log) return(exp(ret))
             ret
         },
         dd = function(x) {
             ex <- exp(-x)
             exp(-ex - x) * (ex - 1)
         },
         ddd = function(x) {
             ex <- exp(-x)
             exp(-x - ex) * (ex - 1)^2 - exp(-ex - 2 * x)
         },
         dd2d = function(x)
             exp(-x) - 1,
         call = ".MaxExtrVal",
         name = "maximum extreme value")

### see 10.1080/15598608.2013.772835
.GammaFrailty <- function(logrho = 0) {
    logrho <- pmax(logrho, log(sqrt(.Machine$double.eps)))
    list(parm = function() c("logrho" = logrho),
         ### note: p(x) is 1 - LaplaceTransform(exp(x))
         p = function(x, lower.tail = TRUE, log.p = FALSE) {
             ### p = 1 - (1 + exp(x + logrho))^(-exp(-logrho))
             ret <- (1 + exp(x + logrho))^(-exp(-logrho))
             if (log.p) {
                 if (lower.tail)
                     return(log1p(-ret))
                 return(log(ret))
             }
             if (lower.tail)
                 return(1 - ret)
             return(ret)
         },
         q = function(p)
             log((1 - p)^(-exp(logrho)) - 1) - logrho,
         d = .d <- function(x, log = FALSE) {
             ret <- x + (-exp(-logrho) - 1) * log(exp(x + logrho) + 1)
             if (!log) return(exp(ret))
             ret
         },
         dd = .dd <- function(x) {
             exlr <- exp(x + logrho)
             memlr <- -exp(-logrho)
             (memlr - 1) * (exlr + 1)^(memlr - 2) * exp(2 * x + logrho) + 
               exp(x) * (exlr + 1)^(memlr - 1)
         },
         ddd = function(x) {
             exlr <- exp(x + logrho)
             memlr <- -exp(-logrho)
             (memlr - 2) * (memlr - 1) * (exlr + 1)^(memlr - 3) * exp(3 * x + 2 * logrho) +
               3 * (memlr - 1) * (exlr + 1)^(memlr - 2) * exp(2 * x + logrho) +
               exp(x) * (exlr + 1)^(memlr - 1)
         },
         dd2d = function(x) {
             .dd(x) / .d(x)
         },
         support = c(-25, log(10)),
         call = ".GammaFrailty",
         name = paste0("GammaFrailty(rho = ", 
                       round(exp(logrho), 
                       options("digits")$digits), ")"))
}

### see 10.1002/sim.687
.InvGaussFrailty <- function(logtheta = 0) {
    logtheta <- pmax(logtheta, log(sqrt(.Machine$double.eps)))
    list(parm = function() c("logtheta" = logtheta),
         ### note: p(x) is 1 - LaplaceTransform(exp(x))
         p = function(x, lower.tail = TRUE, log.p = FALSE) {
             ### p = 1 - exp(- sqrt(4 * exp(logtheta) * (exp(logtheta) + exp(x))) + 2 * exp(logtheta))
             ret <- exp(- sqrt(4 * exp(logtheta) * (exp(logtheta) + exp(x))) + 2 * exp(logtheta))
             if (log.p) {
                 if (lower.tail)
                     return(log1p(-ret))
                 return(log(ret))
             }
             if (lower.tail)
                 return(1 - ret)
             return(ret)
         },
         q = function(p) {
             theta <- exp(logtheta)
             log((-log1p(- p) + 2 * theta)^2 / (4 * theta) - theta)
         },
         d = .d <- function(x, log = FALSE) {
             ret <- (logtheta - 2 * sqrt(exp(logtheta) * (exp(x) + exp(logtheta))) + x + 2 * exp(logtheta)) - 
                     .5 * (logtheta + log(exp(x) + exp(logtheta)))
             if (!log) return(exp(ret))
             ret
         },
         dd = .dd <- function(x) {
             theta <- exp(logtheta)
             stx <- sqrt(theta * (exp(x) + theta))
             ret <- theta * (1 - theta * exp(x) / stx) * exp(-2 * stx + x + 2 * theta)
             ret <- ret / stx
             ret <- ret - theta^2 * exp(-2 * stx + 2 * x + 2 * theta) / (2 * stx * stx^2)
             ret
         },
         ddd = function(x) {
             theta <- exp(logtheta)
             txt <- theta * (exp(x) + theta)
             stx <- sqrt(txt)
             ret <- 3 * theta^3 * exp(-2 * stx + 3 * x + 2 * theta) / (4 * txt^(5/2))
             ret <- ret - theta^2 * (2 - theta * exp(x) / stx) * exp(-2 * stx + 2 * x + 2 * theta) / (2 * txt^(3 / 2))
             ret <- ret - theta^2 * (1 - theta * exp(x) / stx) * exp(-2 * stx + 2 * x + 2 * theta) / (2 * txt^(3 / 2))
             ret <- ret + theta * (theta^2 * exp(2 * x) / (2 * txt^(3/2)) - theta * exp(x) / stx) * exp(-2 * stx + x + 2 * theta) / stx
             ret <- ret + theta * (1 - theta * exp(x) / stx)^2 * exp(-2 * stx + x + 2 * theta) / stx
             ret
         },
         dd2d = function(x) {
             .dd(x) / .d(x)
         },
         call = ".InvGaussFrailty",
         support = c(-25, log(10)),
         name = paste0("InvGaussFrailty(theta = ", 
                       round(exp(logtheta), 
                       options("digits")$digits), ")"))
}

### see 10.18637/jss.v051.i11
.PositiveStableFrailty <- function(logitalpha = 0) {
    list(parm = function() c("logitalpha" = logitalpha),
         ### note: p(x) is 1 - LaplaceTransform(exp(x))
         p = function(x, lower.tail = TRUE, log.p = FALSE) {
             ### p = 1 - exp(-exp(plogis(logitalpha) * x))
             ret <- exp(-exp(plogis(logitalpha) * x))
             if (log.p) {
                 if (lower.tail)
                     return(log1p(-ret))
                 return(log(ret))
             }
             if (lower.tail)
                 return(1 - ret)
             return(ret)
         },
         q = function(p)
             log(-log1p(-p)) / plogis(logitalpha),
         d = .d <- function(x, log = FALSE) {
             alpha <- plogis(logitalpha)
             ret <- plogis(logitalpha, log.p = TRUE) + 
                    alpha * x - exp(alpha * x) 
             if (!log) return(exp(ret))
             ret
         },
         dd = .dd <- function(x) {
             alpha <- plogis(logitalpha)
             eax <- exp(alpha * x)
             alpha * (alpha - alpha * eax) * exp(alpha * x - eax)
         },
         ddd = function(x) {
             alpha <- plogis(logitalpha)
             eax <- exp(alpha * x)
             ret <- alpha * (alpha - alpha * eax)^2 * exp(alpha * x - eax)
             ret <- ret - alpha^3 * exp(2 * alpha * x - eax)
             ret
         },
         dd2d = function(x) {
             .dd(x) / .d(x)
         },
         support = qlogis(sqrt(.Machine$double.eps), 
                          lower.tail = FALSE) * c(-1, 1),
         call = ".PositiveStableFrailty",
         name = paste0("PositiveStableFrailty(alpha = ", 
                       round(plogis(logitalpha), 
                       options("digits")$digits), ")"))
}


### Cure rate mixture models: 10.1002/sim.687
.CureRate <- function(logitrho = 0, ..., distr = .MinExtrVal) {
    d <- distr(...)
    list(parm = function() c("logitrho" = logitrho, d$parm()),
         p = function(x, lower.tail = TRUE, log.p = FALSE) {
             if (log.p) {
                 if (lower.tail)
                     return(plogis(logitrho, log.p = TRUE) + d$p(x, log.p = TRUE))
                 return(log1p(-plogis(logitrho) * d$p(x)))
             }
             if (lower.tail)
                 return(plogis(logitrho) * d$p(x))
             return(1 - plogis(logitrho) * d$p(x))
         },
         q = function(p)
             d$q(p / plogis(logitrho)),
         d = .d <- function(x, log = FALSE) {
             if (log)
                 return(plogis(logitrho, log.p = TRUE) + d$d(x, log = TRUE))
             plogis(logitrho) * d$d(x)
         },
         dd = .dd <- function(x)
             plogis(logitrho) * d$dd(x),
         ddd = function(x)
             plogis(logitrho) * d$ddd(x),
         dd2d = function(x)
             d$dd2d(x),
         call = ".CureRate",
         support = qlogis(sqrt(.Machine$double.eps), 
                          lower.tail = FALSE) * c(-1, 1),
         name = paste0("CureRate(rho = ", 
                       round(plogis(logitrho), 
                       options("digits")$digits), ")"))
}

.distr <- function(which = c("Normal", "Logistic", 
                             "MinExtrVal", "MaxExtrVal", "Exponential")) {
    which <- match.arg(which)
    do.call(paste(".", which, sep = ""), list())
}

